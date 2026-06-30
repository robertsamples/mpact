"""Unit tests for the CV-plot quality metrics (``qualityscore.py``).

The key test is a *faithfulness* check: ``_reference_inline`` is a verbatim
copy of the original computation that used to live in ``plotting.prev_cv.plot``
(before it was extracted into ``qualityscore``), and the extracted function is
asserted to reproduce its numbers exactly -- so the refactor provably did not
change any displayed Reproducibility / Skewness / Overall value. Plus a couple
of sanity checks and a run against the real bundled example dataset's iondict.
"""

import math
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

import qualityscore

REPO_ROOT = Path(__file__).resolve().parents[2]
EXAMPLE_DIR = REPO_ROOT / 'rawdata' / 'PTY087I2'


def _reference_inline(iondict, average_n):
    """Verbatim copy of the original prev_cv.plot() computation (pre-extraction),
    returning (rep, sumskew, qualscore). Do not 'clean up' the ALGORITHM -- it
    exists to pin the extracted module to the exact historical behaviour.

    One mechanical exception: the original assigned the rescaled percentage
    back via ``.iloc[:, 0] = <float Series>`` directly into an int64-dtype
    rank column. That was already a FutureWarning on the pandas pinned in
    qualityscore.py's own fix (see there); on a newer pandas resolved by CI
    (this repo's tests.yml installs an unpinned ``pandas`` for Python 3.11)
    the same line raises a hard TypeError instead, which would make this
    reference function -- not the code under test -- unable to run at all.
    Assigning by column LABEL instead of positional .iloc (replacing the
    column's dtype wholesale rather than casting a float into the existing
    int64 column) is dtype-mechanics only and was already proven
    value-identical by qualityscore.py's own equivalent fix.
    """
    iondict = iondict[~np.isnan(iondict['average CV'])]
    iondictmean = iondict.sort_values(['average CV']).reset_index()
    iondictmed = iondict.sort_values(['median CV']).reset_index()
    iondictmean = iondictmean.reset_index()
    iondictmed = iondictmed.reset_index()
    mean_col0, med_col0 = iondictmean.columns[0], iondictmed.columns[0]
    iondictmean[mean_col0] = 100 * iondictmean.iloc[:, 0] / len(iondictmean['average CV'])
    iondictmed[med_col0] = 100 * iondictmed.iloc[:, 0] / len(iondictmed['median CV'])
    modelstdevlist = [1] + [0] * (int(average_n) - 1)
    modelstdev = pd.Series(modelstdevlist).std() / pd.Series(modelstdevlist).mean()
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
    sumskew = 0
    if math.isnan(modelstdev):
        modelstdev = 1.7
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
    return rep, sumskew, qualscore


def _synthetic_iondict(n=200, seed=0):
    rng = np.random.RandomState(seed)
    # Plausible CVs: mostly low, a tail of high ones; median a touch below mean.
    avg = np.abs(rng.normal(0.15, 0.08, size=n))
    med = np.clip(avg - np.abs(rng.normal(0.01, 0.01, size=n)), 0, None)
    return pd.DataFrame({'average CV': avg, 'median CV': med},
                        index=[f'f{i}' for i in range(n)])


# --------------------------------------------------------------------------- #
# noise model
# --------------------------------------------------------------------------- #

def test_noise_model_cv_matches_count_statistics():
    # [1, 0, 0] -> std/mean of that series.
    s = pd.Series([1, 0, 0])
    assert qualityscore.noise_model_cv(3) == pytest.approx(s.std() / s.mean())


def test_noise_model_cv_falls_back_when_undefined():
    # average_n = 1 -> series is just [1] -> std is NaN -> fallback 1.7.
    assert qualityscore.noise_model_cv(1) == 1.7


# --------------------------------------------------------------------------- #
# faithful extraction
# --------------------------------------------------------------------------- #

@pytest.mark.parametrize('average_n', [3, 4, 6])
def test_matches_original_inline_on_synthetic(average_n):
    iondict = _synthetic_iondict()
    result = qualityscore.compute_cv_quality(iondict, average_n)
    ref_rep, ref_skew, ref_qual = _reference_inline(iondict, average_n)
    assert result.rep == pytest.approx(ref_rep)
    assert result.sumskew == pytest.approx(ref_skew)
    assert result.qualscore == pytest.approx(ref_qual)


def test_result_fields_are_finite_and_sensible():
    result = qualityscore.compute_cv_quality(_synthetic_iondict(), average_n=3)
    assert np.isfinite(result.rep) and result.rep > 0
    assert np.isfinite(result.sumskew)
    assert np.isfinite(result.qualscore)
    # The rarefaction curves run from ~0 to 100% of features.
    assert result.iondictmean.iloc[:, 0].max() == pytest.approx(100, abs=1)


def test_matches_original_inline_on_real_iondict(tmp_path):
    """If the example dataset is present, run the real pipeline once and check
    the extracted metrics match the original inline computation on the real
    iondict.csv (the strongest faithfulness guarantee)."""
    csv = EXAMPLE_DIR / '200826_PTY087I2codingdataset.csv'
    if not csv.exists():
        pytest.skip('example dataset not present')

    from MSFaST import analysis_parameters, run_MSFaST
    from groupsets import GroupSetModel, build_query_dict

    params = analysis_parameters()
    params.filename = csv
    params.samplelistfilename = EXAMPLE_DIR / 'samplelist.csv'
    params.extractmetadatafilename = EXAMPLE_DIR / 'extractmetadata.csv'
    params.outputdir = tmp_path / params.filename.stem
    params.outputdir.mkdir(parents=True)
    params.relfil = True; params.merge = True
    params.ringingwin = 0.5; params.isopeakwin = 0.01; params.dimerpeakwin = 0.01
    params.RTwin = 0.005; params.maxisowin = 3
    params.grpave = True; params.prperr = True
    params.blnkfltr = True; params.blnkgrp = 'Blanks'; params.blankfilthresh = 0.01
    params.CVfil = True; params.cvthresh = 0.2; params.cvparam = 'median CV'
    params.decon = True; params.deconthresh = 0.95
    params.FC = False; params.Ttest = False
    params.statstgrps = ['250um_Ce', '0um_Ce']
    params.graphfilters = ['cv', 'rel', 'insource']
    params.MZRTplt = params.FC3Dplt = params.KMD = False
    params.PCA = params.Dendrogram = params.Volcanoplt = False
    model = GroupSetModel()
    model.add('Features not in blanks', all_groups=['Blanks', 'Media', '0um_Ce', '250um_Ce'])
    model.update(0, src=['Media', '0um_Ce', '250um_Ce'], excl=['Blanks'])
    params.querydict = build_query_dict(model, params.graphfilters)
    params.querylist = list(params.querydict.keys())
    run_MSFaST(params)

    iondict = pd.read_csv(params.outputdir / 'iondict.csv', header=0, index_col=0)
    header = pd.read_csv(params.outputdir / (params.filename.stem + '_filtered.csv'),
                         sep=',', header=None, index_col=[0, 1, 2]).iloc[:3, :].transpose()
    header.columns = ['Biolgroup', 'Sample', 'Injection']
    average_n = header['Injection'].nunique() / header['Sample'].nunique()

    result = qualityscore.compute_cv_quality(iondict, average_n)
    ref_rep, ref_skew, ref_qual = _reference_inline(iondict, average_n)
    assert result.qualscore == pytest.approx(ref_qual)
    assert result.rep == pytest.approx(ref_rep)
    assert result.sumskew == pytest.approx(ref_skew)
