"""Tests for pvclust: AU/BP support by multiscale bootstrap."""

import numpy as np
import pandas as pd
import pytest
from scipy.stats import norm

import pvclust


def separated_groups(seed=3):
    rng = np.random.default_rng(seed)
    centres = rng.normal(scale=3, size=(3, 200))
    objects = np.vstack([centres[g] + rng.normal(size=(4, 200)) for g in range(3)])
    names = [f'g{g}_{i}' for g in range(3) for i in range(4)]
    return pd.DataFrame(objects.T, columns=names)


def test_the_fit_recovers_v_and_c_from_exact_probabilities():
    scales = pvclust.DEFAULT_SCALES
    v, c = 0.4, 0.9
    bp = norm.sf(v * np.sqrt(scales) + c / np.sqrt(scales))
    au, bp1, _, _, pchi, fitted_v, fitted_c = pvclust.fit(scales, bp, 1000)
    assert (fitted_v, fitted_c) == pytest.approx((v, c), abs=1e-9)
    assert au == pytest.approx(norm.sf(v - c))
    assert bp1 == pytest.approx(norm.sf(v + c))
    assert pchi == pytest.approx(1.0)


def test_a_cluster_found_every_time_has_full_support():
    assert pvclust.fit(pvclust.DEFAULT_SCALES, np.ones(10), 100)[:2] == [1.0, 1.0]
    assert pvclust.fit(pvclust.DEFAULT_SCALES, np.zeros(10), 100)[:2] == [0.0, 0.0]


def test_cluster_masks_follow_the_linkage():
    # objects 0 and 2 join first, then 1 joins them
    linkage = np.array([[0, 2, 1.0, 2], [1, 3, 2.0, 3]])
    assert pvclust.cluster_masks(linkage, 3) == [0b101, 0b111]


def test_true_groups_are_supported():
    data = separated_groups()
    pv = pvclust.PvClust(data, nboot=200, seed=0)
    assert list(pv.result.columns) == pvclust.COLUMNS
    assert len(pv.result) == len(pv.clusters) == data.shape[1] - 1
    groups = [{f'g{g}_{i}' for i in range(4)} for g in range(3)]
    support = {frozenset(c): au for c, au in zip(pv.clusters, pv.result['AU'])}
    assert all(support[frozenset(group)] > 0.95 for group in groups)


def test_the_tree_is_the_ward_tree_of_the_columns():
    from scipy.cluster.hierarchy import linkage
    data = separated_groups()
    pv = pvclust.PvClust(data, nboot=10, seed=0)
    assert np.allclose(pv.linkage_matrix, linkage(data.T.to_numpy(), method='ward'))


def test_a_seed_makes_a_run_repeatable():
    data = separated_groups()
    first = pvclust.PvClust(data, nboot=50, seed=5).result
    second = pvclust.PvClust(data, nboot=50, seed=5).result
    pd.testing.assert_frame_equal(first, second)


def test_scales_given_in_tenths_are_read_as_fractions():
    data = separated_groups()
    tenths = pvclust.PvClust(data, nboot=5, r=np.array(range(5, 15)), seed=0)
    assert np.allclose(tenths.r, pvclust.PvClust(data, nboot=5, seed=0).r)
    assert tenths.n_scaled[0] == 100 and tenths.n_scaled[-1] == 280


def test_plot_labels_every_join_at_its_place_in_the_drawing():
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    data = separated_groups()
    pv = pvclust.PvClust(data, nboot=20, seed=0)
    fig, ax = plt.subplots()
    dend = pv.plot(ax, labels=list(data.columns))
    x, y = pv.node_positions(dend['leaves'])
    # Each join drawn by scipy is a segment whose top midpoint is (x, y).
    drawn = sorted((round((seg_x[1] + seg_x[2]) / 2, 6), round(seg_y[1], 6))
                   for seg_x, seg_y in zip(dend['icoord'], dend['dcoord']))
    assert drawn == sorted(zip(np.round(x, 6), np.round(y, 6)))
    texts = [t.get_text() for t in ax.texts]
    assert len(texts) == 2 * len(pv.result)
    plt.close(fig)
