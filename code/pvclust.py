"""
MPACT
Copyright 2026, Robert M. Samples

Support values for hierarchical clustering by multiscale bootstrap resampling.

For every cluster in the dendrogram of the columns of a data matrix, two
p-values are computed:

- BP, the bootstrap probability: how often the cluster appears when the rows
  are resampled with replacement;
- AU, the approximately unbiased p-value, from the multiscale bootstrap.

Method: H. Shimodaira (2004), Approximately unbiased tests of regions using
multistep-multiscale bootstrap resampling, Annals of Statistics 32, 2616-2641;
R. Suzuki and H. Shimodaira (2006), Pvclust: an R package for assessing the
uncertainty in hierarchical clustering, Bioinformatics 22, 1540-1542.

The rows are resampled at several sample sizes n' = r * n. At each scale r
the share of replicates containing a cluster, BP(r), gives

    z(r) = qnorm(1 - BP(r)),

and the model z(r) = v * sqrt(r) + c / sqrt(r) is fitted to the scales where
0 < BP(r) < 1 by weighted least squares, each z(r) weighted by the inverse of
its binomial variance. Then

    AU = 1 - pnorm(v - c),    BP = 1 - pnorm(v + c),

with standard errors by the delta method and a chi-square test of the fit.
The values agree with the R package pvclust to within bootstrap error.
"""

from multiprocessing import Pool

import numpy as np
import pandas as pd
from scipy.cluster import hierarchy
from scipy.spatial.distance import pdist
from scipy.stats import chi2, norm

try:
    import fastcluster
except ImportError:  # optional: scipy is used without it
    fastcluster = None

# Scales r = n'/n, as in pvclust.
DEFAULT_SCALES = np.round(np.arange(0.5, 1.45, 0.1), 1)
COLUMNS = ['AU', 'BP', 'SE.AU', 'SE.BP', 'pchi', 'v', 'c']


def linkage(data, method='ward', metric='euclidean'):
    """Linkage matrix of the columns of ``data`` (variables x objects)."""
    distances = pdist(np.asarray(data, dtype=float).T, metric=metric)
    if fastcluster is not None:
        return fastcluster.linkage(distances, method=method)
    return hierarchy.linkage(distances, method=method)


def cluster_masks(linkage_matrix, n_objects):
    """The objects in each cluster of a linkage matrix, one bit per object,
    in linkage row order."""
    members = [1 << i for i in range(n_objects)]
    for left, right in np.asarray(linkage_matrix)[:, :2].astype(int):
        members.append(members[left] | members[right])
    return members[n_objects:]


def _count_scale(job):
    """For one scale: how many of ``nboot`` replicates contain each target
    cluster. Run in a worker when ``parallel`` is set."""
    data, n_rows, targets, nboot, method, metric, seed = job
    rng = np.random.default_rng(seed)
    n_objects = data.shape[1]
    counts = np.zeros(len(targets), dtype=np.int64)
    for _ in range(nboot):
        rows = rng.integers(0, data.shape[0], size=n_rows)
        found = cluster_masks(linkage(data[rows], method, metric), n_objects)
        for mask in found:
            position = targets.get(mask)
            if position is not None:
                counts[position] += 1
    return counts


def fit(scales, bp, nboot):
    """Fit z(r) = v sqrt(r) + c / sqrt(r) to one cluster's bootstrap
    probabilities. Returns the values for ``COLUMNS``."""
    scales = np.asarray(scales, dtype=float)
    bp = np.asarray(bp, dtype=float)
    usable = (bp > 0) & (bp < 1)
    if usable.sum() < 2:
        # Too few scales to fit: always (or never) found at almost every
        # scale. Report what was seen at the scale nearest r = 1.
        seen = float(bp[np.argmin(np.abs(scales - 1))])
        return [seen, seen, 0.0, 0.0, np.nan, np.nan, np.nan]

    r, p = scales[usable], bp[usable]
    z = norm.ppf(1 - p)
    variance = p * (1 - p) / (norm.pdf(z) ** 2 * nboot)
    weight = 1 / variance
    design = np.column_stack([np.sqrt(r), 1 / np.sqrt(r)])
    information = design.T @ (design * weight[:, None])
    covariance = np.linalg.inv(information)
    v, c = covariance @ (design.T @ (weight * z))

    residual = z - design @ np.array([v, c])
    statistic = float(np.sum(weight * residual ** 2))
    df = int(usable.sum()) - 2
    pchi = float(chi2.sf(statistic, df)) if df > 0 else 1.0

    au_z, bp_z = v - c, v + c
    se_au = norm.pdf(au_z) * np.sqrt(max(covariance[0, 0] + covariance[1, 1]
                                         - 2 * covariance[0, 1], 0.0))
    se_bp = norm.pdf(bp_z) * np.sqrt(max(covariance[0, 0] + covariance[1, 1]
                                         + 2 * covariance[0, 1], 0.0))
    return [float(norm.sf(au_z)), float(norm.sf(bp_z)), float(se_au), float(se_bp),
            pchi, float(v), float(c)]


class PvClust:
    """AU and BP support for every cluster of a hierarchical clustering.

    ``data`` is variables x objects, as in pvclust: the columns are clustered
    and the rows are resampled. ``r`` is the list of scales, as fractions of
    the number of rows (whole numbers are read as tenths, so the older
    ``range(5, 15)`` still means 0.5 to 1.4). ``seed`` makes a run
    repeatable; results do not depend on ``parallel``.

    After construction:
      ``linkage_matrix``  the clustering of the full data
      ``clusters``        each cluster's objects, by column label
      ``bp_by_scale``     clusters x scales bootstrap probabilities
      ``result``          one row per cluster, in linkage row order, with
                          the columns AU, BP, SE.AU, SE.BP, pchi, v, c
    """

    def __init__(self, data, method='ward', metric='euclidean', nboot=1000,
                 r=None, parallel=False, seed=None):
        frame = pd.DataFrame(data)
        self.data = frame
        self.method = method
        self.metric = metric
        self.nboot = int(nboot)
        self.parallel = bool(parallel)

        scales = DEFAULT_SCALES if r is None else np.asarray(r, dtype=float)
        if np.all(scales >= 2) and np.allclose(scales, np.round(scales)):
            scales = scales / 10
        values = frame.to_numpy(dtype=float)
        self.n = values.shape[0]
        self.n_scaled = np.maximum(np.round(scales * self.n).astype(int), 1)
        # The scales actually used: n' is a whole number of rows.
        self.r = self.n_scaled / self.n

        self.linkage_matrix = linkage(values, method, metric)
        masks = cluster_masks(self.linkage_matrix, values.shape[1])
        labels = list(frame.columns)
        self.clusters = [[labels[i] for i in range(len(labels)) if mask >> i & 1]
                         for mask in masks]
        targets = {mask: position for position, mask in enumerate(masks)}

        seeds = np.random.SeedSequence(seed).spawn(len(self.n_scaled))
        jobs = [(values, int(n_rows), targets, self.nboot, method, metric, child)
                for n_rows, child in zip(self.n_scaled, seeds)]
        if self.parallel and len(jobs) > 1:
            with Pool() as pool:
                counts = pool.map(_count_scale, jobs)
        else:
            counts = [_count_scale(job) for job in jobs]

        self.bp_by_scale = pd.DataFrame(np.column_stack(counts) / self.nboot,
                                        columns=np.round(self.r, 3))
        self.result = pd.DataFrame(
            [fit(self.r, row, self.nboot) for row in self.bp_by_scale.to_numpy()],
            columns=COLUMNS)

    def node_positions(self, leaves):
        """(x, y) of each cluster's join in a scipy dendrogram whose leaf
        order is ``leaves``, in linkage row order. scipy places the i-th leaf
        at x = 5 + 10 i and each join midway between its two branches."""
        n_objects = len(leaves)
        x = np.empty(2 * n_objects - 1)
        x[np.asarray(leaves)] = 5 + 10 * np.arange(n_objects)
        for row, (left, right) in enumerate(self.linkage_matrix[:, :2].astype(int)):
            x[n_objects + row] = (x[left] + x[right]) / 2
        return x[n_objects:], self.linkage_matrix[:, 2]

    def plot(self, ax, labels=None, link_color_func=None, fontsize=7):
        """Draw the dendrogram on a matplotlib ``ax`` with AU (red, left) and
        BP (green, right) in percent at each join, as pvclust does. Returns
        scipy's dendrogram dictionary."""
        dend = hierarchy.dendrogram(
            self.linkage_matrix, ax=ax, labels=labels, leaf_rotation=90,
            color_threshold=0, above_threshold_color='black',
            link_color_func=link_color_func)
        x, y = self.node_positions(dend['leaves'])
        lift = 0.01 * (y.max() if len(y) else 1)
        for (au, bp), xi, yi in zip(self.result[['AU', 'BP']].to_numpy(), x, y):
            ax.text(xi - 1, yi + lift, f'{round(100 * au)}', color='red',
                    fontsize=fontsize, ha='right', va='bottom')
            ax.text(xi + 1, yi + lift, f'{round(100 * bp)}', color='green',
                    fontsize=fontsize, ha='left', va='bottom')
        return dend
