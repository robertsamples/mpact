"""
MPACT
Copyright 2022, Robert M. Samples, Sara P. Puckett, and Marcy J. Balunas

Qt-free extraction of the data-quality metrics shown on the CV (coefficient of
variation) rarefaction plot tab. These numbers -- Reproducibility, Skewness,
and an Overall quality score -- already existed, but the math was buried inside
``plotting.prev_cv.plot()`` (entangled with the matplotlib drawing and the
``lbl_spllist_3`` label update) and had no test coverage.

Moving the computation here (the same pattern as ``ordination.py`` /
``biogroups.py`` / ``dbsearch.py`` / ``clusterpurity.py``) makes it unit-
testable and keeps ``prev_cv`` as a thin draw-the-result wrapper. The logic is
a faithful port of the original inline code -- ``tests/test_qualityscore.py``
pins it against a verbatim copy of that original to guarantee the extraction
didn't change any displayed number.

Definitions (as originally implemented):
- The CV rarefaction curve plots, for the mean-CV and median-CV orderings, the
  cumulative percentage of features (0-100) against their CV. ``aucav``/
  ``aucmed`` are the areas under those two curves (percentage x CV).
- ``modelstdev`` is the maximum theoretical CV expected from pure
  count-statistics noise given the average number of injections per sample
  (``[1] + [0]*(n-1)`` treated as a sample -> its CV); it sets the CV axis
  scale and normalises the AUC.
- ``rep`` (Reproducibility) = mean AUC / (modelstdev x 100): how far left
  (low-CV) the curve sits relative to the noise-model scale.
- ``sumskew`` (Skewness) = normalised integrated gap between the mean-CV and
  median-CV curves: how asymmetric the CV distribution is.
- ``qualscore`` (Overall) = (1 - skew) x rep x 100.
"""

import math

import numpy as np
import pandas as pd


class CVQualityResult:
    """Bundle of everything ``prev_cv`` needs to label and draw the CV plot.

    Attributes:
        iondictmean: features sorted by 'average CV', with column 0 replaced by
            the cumulative percentage (0-100) -- the mean-CV rarefaction curve.
        iondictmed: same, sorted/ranked by 'median CV'.
        modelstdev: the count-statistics noise-model CV (axis scale).
        rep: reproducibility fraction (0-1); ``100*rep`` is the displayed %.
        sumskew: skewness fraction (0-1); ``100*sumskew`` is the displayed %.
        qualscore: overall quality score (already on a 0-100 scale).
    """
    __slots__ = ('iondictmean', 'iondictmed', 'modelstdev', 'rep', 'sumskew', 'qualscore')

    def __init__(self, iondictmean, iondictmed, modelstdev, rep, sumskew, qualscore):
        self.iondictmean = iondictmean
        self.iondictmed = iondictmed
        self.modelstdev = modelstdev
        self.rep = rep
        self.sumskew = sumskew
        self.qualscore = qualscore


def noise_model_cv(average_n):
    """Maximum theoretical CV from pure presence/absence count noise given
    ``average_n`` injections per sample: the CV of ``[1] + [0]*(n-1)``.

    Falls back to 1.7 when undefined (e.g. a single injection per sample makes
    the model series a single value with no spread) -- matching the original.
    """
    modelstdevlist = [1] + [0] * (int(average_n) - 1)
    series = pd.Series(modelstdevlist)
    modelstdev = series.std() / series.mean()
    if math.isnan(modelstdev):
        modelstdev = 1.7
    return modelstdev


def compute_cv_quality(iondict, average_n):
    """Compute the CV-plot quality metrics from an ion dictionary.

    Args:
        iondict: DataFrame with 'average CV' and 'median CV' columns (the
            ``iondict.csv`` produced by the CV filter). Rows with NaN
            'average CV' are dropped, matching the plot.
        average_n: average number of injections per sample (used for the
            count-statistics noise model that scales the CV axis).

    Returns:
        CVQualityResult.
    """
    iondict = iondict[~np.isnan(iondict['average CV'])]

    # Cumulative-percentage rarefaction curves for the mean- and median-CV
    # orderings. The double reset_index reproduces the original exactly:
    # the first moves the Compound index to a column, the second materialises
    # the 0..n-1 rank as column 0, which is then rescaled to a 0-100 percentage.
    iondictmean = iondict.sort_values(['average CV']).reset_index()
    iondictmed = iondict.sort_values(['median CV']).reset_index()
    iondictmean = iondictmean.reset_index()
    iondictmed = iondictmed.reset_index()
    # Replace column 0 (the integer rank) with the 0-100 cumulative percentage.
    # Assign by column LABEL (not in-place via .iloc) so the column's dtype is
    # replaced wholesale rather than an incompatible float cast into an int64
    # column -- the latter raises a FutureWarning on pandas 2.x and is slated
    # to become a hard error (same class as the LossySetitemError fixed
    # elsewhere). Numerically identical to the original.
    mean_col0, med_col0 = iondictmean.columns[0], iondictmed.columns[0]
    iondictmean[mean_col0] = 100 * iondictmean.iloc[:, 0] / len(iondictmean['average CV'])
    iondictmed[med_col0] = 100 * iondictmed.iloc[:, 0] / len(iondictmed['median CV'])

    modelstdev = noise_model_cv(average_n)

    # Area under each rarefaction curve (percentage integrated over CV).
    prevav = 0
    aucav = 0
    prevmed = 0
    aucmed = 0
    for pos in range(0, len(iondictmean.iloc[:, 0])):
        dist = iondictmean.iloc[pos, :]['average CV'] - prevav
        aucav += dist * iondictmean.iloc[pos, 0]
        prevav = iondictmean.iloc[pos, :]['average CV']

        dist = iondictmed.iloc[pos, :]['median CV'] - prevmed
        aucmed += dist * iondictmed.iloc[pos, 0]
        prevmed = iondictmed.iloc[pos, :]['median CV']

    # Integrated gap between the mean and median curves (distribution skew).
    sumskew = 0
    for val in range(1, int((modelstdev * 100))):
        pos = val / 100
        meanav = iondictmean[abs(iondictmean['average CV'] - pos - modelstdev / 200) < modelstdev / 200].iloc[:, 0].mean()
        meanmed = iondictmed[abs(iondictmed['average CV'] - pos - modelstdev / 200) < modelstdev / 200].iloc[:, 0].mean()
        skew = abs(meanmed - meanav)
        if not np.isnan(skew):
            sumskew += skew * modelstdev / 100

    sumskew = sumskew / ((aucmed + aucav) / 2)
    rep = ((aucmed + aucav) / 2) / (modelstdev * 100)
    qualscore = (1 - sumskew) * rep * 100

    return CVQualityResult(iondictmean, iondictmed, modelstdev, rep, sumskew, qualscore)
