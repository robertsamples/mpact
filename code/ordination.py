"""
MPACT
Copyright 2022, Robert M. Samples, Sara P. Puckett, and Marcy J. Balunas

Qt-free multivariate ordination backend: PCA, NMDS, and PLS-DA on the
samples x features intensity matrix, plus the data prep (technical-replicate
collapsing) and loadings-selection logic the "multivariate" plot tab needs.

OPLS-DA is intentionally not implemented here: scikit-learn has no native
support, and the only alternatives (the unmaintained ``pyopls`` package, or a
from-scratch orthogonal-signal-correction implementation) are both
meaningfully riskier than PCA/NMDS/PLS-DA without a reference dataset to
validate against. Logged as future work, not started.

This module is Qt-free and unit-tested (see ``tests/test_ordination.py``).
"""

import numpy as np
import pandas as pd
from sklearn.cross_decomposition import PLSRegression
from sklearn.decomposition import PCA
from sklearn.metrics import pairwise_distances
from sklearn import manifold


def load_ordination_matrix(file, raw_msdata_header, collapse_replicates):
    """Load the samples x features intensity matrix used for ordination.

    This is a near-verbatim port of the data-loading half of the original
    (dead-checkbox) ``plot_PCA.plot()`` -- deliberately not redesigned, since
    the original's row-grouping math is correct (verified empirically in
    ``test_ordination.py``/the scratch script, not re-derived by inspection
    here -- this header-juggling is genuinely easy to get subtly wrong by
    reasoning about it instead of testing it). Only the hardcoded
    ``collapsereps = False`` is now a real parameter.

    Args:
        file: path to the canonical ``_filtered.csv`` peak table (3-row
            header: Biolgroup, Sample, Injection; see devnotes.md).
        raw_msdata_header: the same peak table's 3 header rows, read
            *raw* (``header=None, index_col=[0,1,2]).iloc[:3,:].transpose()``,
            exactly as the original code reads it -- NOT yet renamed or
            re-indexed; that happens inside this function (for the
            ``collapse_replicates=True`` case, a *different* header --
            read from the freshly-collapsed intermediate file -- is used
            instead, matching the original's control flow exactly).
        collapse_replicates: if True, average technical replicates (multiple
            Injections under the same Sample) together, keeping biological
            replicates (distinct Samples) separate. If False, every
            Injection is its own row, as-is.

    Returns:
        (X, biolgroup): ``X`` is a DataFrame indexed by sample identifier
        ('File'; an Injection name, or a Sample name when collapsed),
        columns = features. ``biolgroup`` is a Series, same index as ``X``,
        mapping each sample to its biological group.
    """
    if collapse_replicates:
        # Average technical replicates (Injection) while keeping biological
        # replicates (Sample) distinct -- level order is Compound, m/z, RT,
        # Biolgroup, Sample, Injection (MSFaST.py's msdata_header.columns
        # assignment). Round-trips through a CSV (matching the original)
        # so the relabeled 3-row header can be read back the same way the
        # uncollapsed path reads the real file, rather than hand-deriving
        # unstack()'s resulting column-level order.
        msdata = pd.read_csv(file, sep=',', header=[0, 1, 2], index_col=[0, 1, 2])
        try:
            msdata = msdata.stack([0, 1, 2], future_stack=True)
        except TypeError:
            msdata = msdata.stack([0, 1, 2])
        msdata = msdata.groupby(level=[0, 1, 2, 3, 4]).mean().unstack(level=[-1, -2])
        collapsed_columns = msdata.columns.to_list()
        msdata = msdata.reset_index()
        header = [('', '', 'Compound'), ('', '', 'm/z'), ('', '', 'Retention time (min)')]
        for elem in collapsed_columns:
            header.append((elem[1], '', elem[0]))
        msdata.columns = pd.MultiIndex.from_tuples(header)
        msdata.to_csv('averagepca.csv', header=True, index=False)

        msdata_header = pd.read_csv('averagepca.csv', sep=',', header=None,
                                     index_col=[0, 1, 2]).iloc[:3, :].transpose()
        pcadf = (pd.read_csv('averagepca.csv', sep=',', header=[2], index_col=[0])
                 .drop(['m/z', 'Retention time (min)'], axis=1)
                 .transpose().astype(float).reset_index().rename(columns={'index': 'File'}))
    else:
        msdata_header = raw_msdata_header
        pcadf = (pd.read_csv(file, sep=',', header=[2], index_col=[0])
                 .drop(['m/z', 'Retention time (min)'], axis=1)
                 .transpose().astype(float).reset_index().rename(columns={'index': 'File'}))

    msdata_header.columns = ['Biolgroup', 'Sample', 'Injection']
    msdata_header = msdata_header.set_index('Injection')

    x = pcadf.set_index('File')
    biolgroup = pd.Series(
        {file_id: msdata_header.loc[file_id, 'Biolgroup'] for file_id in pcadf['File']},
        name='Biolgroup',
    )
    biolgroup.index.name = 'File'
    return x, biolgroup


def autoscale(x):
    """Mean-center and scale each feature to unit variance ("UV-scaling" /
    "autoscaling" in chemometrics terminology -- the standard pre-treatment
    for PCA/PLS-DA on mass-spec intensity data).

    Without this, raw intensities (which can span several orders of
    magnitude between features -- confirmed on real example data: feature
    standard deviations ranged from ~1.8 to ~10,000, a ~5800x spread) let a
    handful of high-abundance features dominate both the apparent
    explained-variance and the loadings, regardless of which features
    actually separate the biological groups. NMDS is deliberately NOT put
    through this -- its Bray-Curtis dissimilarity is conventionally computed
    on raw (or relative) abundances, not standardized ones.
    """
    std = x.std(axis=0)
    std = std.replace(0, 1)  # a zero-variance feature would divide by zero; leave it at 0 instead
    return (x - x.mean(axis=0)) / std


def run_pca(x, n_components):
    """PCA on the samples x features matrix, after autoscaling (see
    ``autoscale()``) so the result isn't dominated by whichever features
    happen to have the largest raw intensity.

    Returns:
        (scores, loadings, explained_variance_ratio): ``scores`` is a
        DataFrame (index=samples, columns=PC1..PCn); ``loadings`` is a
        DataFrame (index=features, columns=PC1..PCn) of each feature's
        contribution to each component; ``explained_variance_ratio`` is an
        ndarray of length ``n_components``.
    """
    x_scaled = autoscale(x)
    pca = PCA(n_components=n_components)
    scores = pca.fit_transform(x_scaled.values)
    columns = [f'PC{i + 1}' for i in range(n_components)]
    scores = pd.DataFrame(scores, index=x.index, columns=columns)
    loadings = pd.DataFrame(pca.components_.T, index=x.columns, columns=columns)
    return scores, loadings, pca.explained_variance_ratio_


def run_nmds(x, n_components):
    """Non-metric MDS on Bray-Curtis sample dissimilarities, warm-started
    from a metric MDS solution, then rotated onto principal axes purely for
    a stable/sensible orientation (this rotation is NOT a second ordination
    of the original features -- it doesn't change the NMDS embedding's
    inter-point distances, only its axis orientation).

    Returns:
        (scores, explained_variance_ratio, stress): ``explained_variance_ratio``
        here is the variance of the *embedded* (already-reduced) NMDS
        coordinates explained by each rotated axis -- not, unlike PCA's, a
        measure of how much of the original feature-space variance is
        captured. Callers should label this distinctly (e.g. "% of
        embedding variance") rather than implying it's the same quantity as
        PCA's explained variance.
    """
    similarities = pairwise_distances(x.values - x.values.mean(), metric='braycurtis')

    mds = manifold.MDS(n_components=n_components, max_iter=3000, eps=1e-9,
                        random_state=1, dissimilarity="precomputed", n_jobs=1)
    pos = mds.fit(similarities).embedding_

    nmds = manifold.MDS(n_components=n_components, metric=False, max_iter=3000,
                         eps=1e-12, dissimilarity="precomputed", random_state=1,
                         n_jobs=1, n_init=1)
    npos = nmds.fit_transform(similarities, init=pos)
    stress = nmds.stress_

    pca = PCA(n_components=n_components)
    rotated = pca.fit_transform(npos)
    columns = [f'NMDS{i + 1}' for i in range(n_components)]
    scores = pd.DataFrame(rotated, index=x.index, columns=columns)
    return scores, pca.explained_variance_ratio_, stress


def nmds_loading_proxy(x, scores):
    """Per-feature correlation with each NMDS axis, as a loadings-equivalent.

    NMDS has no linear feature loadings (it's a rank-based embedding, not a
    linear projection of the original features) -- this is the standard
    ecology "vector fitting" approach (cf. vegan::envfit): correlate each
    original feature with each ordination axis and use that as the
    loadings-plot substitute.

    Returns:
        DataFrame (index=features, columns=same as ``scores``) of Pearson
        correlation coefficients.
    """
    return pd.DataFrame(
        {col: x.corrwith(scores[col]) for col in scores.columns},
        index=x.columns,
    )


def run_plsda(x, y, n_components):
    """PLS-DA: PLS regression of the samples x features matrix against
    one-hot-encoded group labels.

    Args:
        x: samples x features DataFrame.
        y: Series of group labels, indexed the same as ``x``.
        n_components: number of PLS components.

    Returns:
        (scores, loadings, explained_variance_ratio): shapes match
        ``run_pca``'s. scikit-learn doesn't expose an explained-variance
        ratio for PLS directly, so it's computed manually here as each
        component's X-score variance divided by the total variance of
        (autoscaled) ``x`` -- the standard approach for reporting
        %-explained on a PLS biplot.
    """
    x_scaled = autoscale(x)
    y_dummies = pd.get_dummies(y)
    # scale=False: we already autoscaled x ourselves (above), consistent
    # with run_pca -- letting PLSRegression's default scale=True scale it
    # AGAIN (and scale the 0/1 dummy y columns, which doesn't make sense for
    # group-membership indicators) would both double-scale x and distort y.
    # (scale=False also previously fixed a bug where x_scores_ lived on a
    # different scale than the unscaled total-variance denominator below,
    # before autoscaling was added -- see devnotes.md.)
    pls = PLSRegression(n_components=n_components, scale=False)
    pls.fit(x_scaled.values, y_dummies.values)
    x_scores = pls.x_scores_
    columns = [f'PLS{i + 1}' for i in range(n_components)]
    scores = pd.DataFrame(x_scores, index=x.index, columns=columns)
    loadings = pd.DataFrame(pls.x_loadings_, index=x.columns, columns=columns)

    total_variance = np.sum(x_scaled.values ** 2)
    component_variance = np.sum(x_scores ** 2, axis=0)
    explained_variance_ratio = component_variance / total_variance
    return scores, loadings, explained_variance_ratio


def top_loadings(loadings, n=25, always_include=()):
    """Subset of ``loadings`` to actually draw on a loadings plot.

    High-dimensional data (thousands of features) can't all be plotted
    legibly, so this returns only the top ``n`` features by vector magnitude
    (Euclidean norm across all of ``loadings``'s columns) -- plus any
    feature named in ``always_include``, even if its magnitude wouldn't
    otherwise make the cut. That's what lets the app highlight a specific
    (possibly tiny) feature on demand without changing the default view.

    Args:
        loadings: DataFrame (index=features, columns=components).
        n: how many top-magnitude features to include by default.
        always_include: iterable of feature names (must be a subset of
            ``loadings.index``) to include regardless of magnitude.

    Returns:
        DataFrame: subset of ``loadings`` (same columns), index order
        preserved from ``loadings``, with at most ``n + len(always_include)``
        rows (fewer if there's overlap or ``loadings`` itself is smaller).
    """
    magnitude = np.sqrt((loadings ** 2).sum(axis=1))
    top_n_index = magnitude.nlargest(min(n, len(loadings))).index
    forced = [feat for feat in always_include if feat in loadings.index]
    keep = top_n_index.union(forced, sort=False)
    # Preserve loadings' original row order rather than magnitude-sorted order.
    keep = [feat for feat in loadings.index if feat in keep]
    return loadings.loc[keep]
