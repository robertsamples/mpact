"""Unit tests for the multivariate ordination backend (``ordination.py``).

Covers data loading/replicate-collapsing (verified against a synthetic
3-header-row peak table with a known technical/biological-replicate
structure -- see the plan for why this is empirically checked rather than
trusted by inspection) and the PCA/NMDS/PLS-DA/top_loadings math against
small synthetic matrices.
"""

import numpy as np
import pandas as pd
import pytest

from ordination import (
    load_ordination_matrix, nmds_loading_proxy, replicate_label_components,
    run_nmds, run_pca, run_plsda, similarity_matrix, top_loadings,
)


# --------------------------------------------------------------------------- #
# load_ordination_matrix / collapse_replicates
# --------------------------------------------------------------------------- #

def _write_synthetic_filtered_csv(path):
    """3 samples (S1, S1b in groupA; S2 in groupB), 3 technical-replicate
    injections each (9 injections total) -- enough to tell "collapsed to
    one row per Sample" apart from both "per Injection" (9) and "per
    Biolgroup" (2, since there are only 2 biolgroups but 3 samples).
    """
    with open(path, 'w') as f:
        f.write(',,,groupA,groupA,groupA,groupA,groupA,groupA,groupB,groupB,groupB\n')
        f.write(',,,S1,S1,S1,S1b,S1b,S1b,S2,S2,S2\n')
        f.write('Compound,m/z,Retention time (min),inj1,inj2,inj3,inj4,inj5,inj6,inj7,inj8,inj9\n')
        f.write('feat1,100.0,1.0,10,12,11,30,32,31,50,52,51\n')
        f.write('feat2,200.0,2.0,5,6,4,15,16,14,20,19,21\n')


def _raw_header(path):
    return pd.read_csv(path, sep=',', header=None, index_col=[0, 1, 2]).iloc[:3, :].transpose()


def test_uncollapsed_keeps_one_row_per_injection(tmp_path):
    path = tmp_path / 'example_filtered.csv'
    _write_synthetic_filtered_csv(path)
    x, biolgroup = load_ordination_matrix(path, _raw_header(path), collapse_replicates=False)
    assert x.shape == (9, 2)
    assert len(biolgroup) == 9


def test_collapsed_averages_technical_not_biological_replicates(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)  # 'averagepca.csv' lands here, not the repo
    path = tmp_path / 'example_filtered.csv'
    _write_synthetic_filtered_csv(path)
    x, biolgroup = load_ordination_matrix(path, _raw_header(path), collapse_replicates=True)

    # 3 distinct Samples (S1, S1b, S2) -- not 9 (uncollapsed) and not 2
    # (the number of Biolgroups, which would mean biological replicates got
    # wrongly merged too).
    assert x.shape[0] == 3
    assert x.shape[1] == 2
    assert biolgroup.nunique() == 2
    assert sorted(biolgroup.unique()) == ['groupA', 'groupB']
    # Two of the three collapsed rows belong to groupA (S1, S1b).
    assert (biolgroup == 'groupA').sum() == 2
    assert (biolgroup == 'groupB').sum() == 1


def test_collapsed_values_are_the_mean_of_their_technical_replicates(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    path = tmp_path / 'example_filtered.csv'
    _write_synthetic_filtered_csv(path)
    x, _ = load_ordination_matrix(path, _raw_header(path), collapse_replicates=True)

    # S1's feat1 replicates are 10, 12, 11 -> mean 11.
    s1_row = x.loc[x.index.str.contains('S1') & ~x.index.str.contains('S1b')]
    assert s1_row['feat1'].iloc[0] == pytest.approx(11.0)


# --------------------------------------------------------------------------- #
# replicate_label_components
# --------------------------------------------------------------------------- #

def test_replicate_label_components_numbers_bio_and_tech_reps(tmp_path):
    # Reuses the same fixture as the collapse tests: groupA has 2 Samples
    # (S1, S1b -- BioRep 1 and 2), groupB has 1 Sample (S2 -- BioRep 1, the
    # "only one biological replicate" edge case), every Sample has 3
    # Injections (TechRep 1-3).
    path = tmp_path / 'example_filtered.csv'
    _write_synthetic_filtered_csv(path)
    components = replicate_label_components(_raw_header(path))

    assert components.loc['inj1', ['Biolgroup', 'Sample', 'BioRep', 'TechRep']].tolist() == ['groupA', 'S1', 1, 1]
    assert components.loc['inj3', ['Biolgroup', 'Sample', 'BioRep', 'TechRep']].tolist() == ['groupA', 'S1', 1, 3]
    assert components.loc['inj4', ['Biolgroup', 'Sample', 'BioRep', 'TechRep']].tolist() == ['groupA', 'S1b', 2, 1]
    # groupB has only one Sample -- still BioRep 1, not skipped/blank.
    assert components.loc['inj7', ['Biolgroup', 'Sample', 'BioRep', 'TechRep']].tolist() == ['groupB', 'S2', 1, 1]


def test_replicate_label_components_single_technical_replicate(tmp_path):
    # Edge case: a Sample with only one Injection should still get
    # TechRep=1, not be skipped or raise.
    path = tmp_path / 'single_techrep_filtered.csv'
    with open(path, 'w') as f:
        f.write(',,,groupA,groupA,groupB\n')
        f.write(',,,S1,S1b,S2\n')
        f.write('Compound,m/z,Retention time (min),inj1,inj2,inj3\n')
        f.write('feat1,100.0,1.0,10,30,50\n')
    components = replicate_label_components(_raw_header(path))

    assert components.loc['inj1', ['Sample', 'BioRep', 'TechRep']].tolist() == ['S1', 1, 1]
    assert components.loc['inj2', ['Sample', 'BioRep', 'TechRep']].tolist() == ['S1b', 2, 1]
    assert components.loc['inj3', ['Sample', 'BioRep', 'TechRep']].tolist() == ['S2', 1, 1]


# --------------------------------------------------------------------------- #
# run_pca / run_nmds / run_plsda / top_loadings
# --------------------------------------------------------------------------- #

def _make_low_rank_matrix():
    # 12 samples, 5 features, but only 2 real underlying dimensions of
    # variation -- PCA on this should recover ~100% explained variance in
    # the first 2 components.
    rng = np.random.RandomState(0)
    latent = rng.normal(size=(12, 2))
    loading_matrix = rng.normal(size=(2, 5))
    x = pd.DataFrame(
        latent @ loading_matrix + rng.normal(scale=0.01, size=(12, 5)),
        index=[f's{i}' for i in range(12)],
        columns=[f'f{i}' for i in range(5)],
    )
    return x


def test_pca_recovers_known_low_rank_structure():
    x = _make_low_rank_matrix()
    scores, loadings, expvar = run_pca(x, n_components=3)
    assert scores.shape == (12, 3)
    assert loadings.shape == (5, 3)
    # Two real dimensions of variation + tiny noise -> first two components
    # should capture almost all the variance.
    assert expvar[:2].sum() > 0.99


def test_plsda_separates_two_groups_along_first_component():
    rng = np.random.RandomState(1)
    group_a = rng.normal(loc=0, scale=0.5, size=(10, 6))
    group_b = rng.normal(loc=5, scale=0.5, size=(10, 6))
    x = pd.DataFrame(
        np.vstack([group_a, group_b]),
        index=[f's{i}' for i in range(20)],
        columns=[f'f{i}' for i in range(6)],
    )
    y = pd.Series(['A'] * 10 + ['B'] * 10, index=x.index)

    scores, loadings, expvar = run_plsda(x, y, n_components=2)
    assert scores.shape == (20, 2)
    assert loadings.shape == (6, 2)
    # The groups are cleanly separated along PLS1: every A score should be
    # on one side of 0 and every B score on the other (sign is arbitrary).
    pls1 = scores['PLS1']
    assert (pls1[:10] > 0).all() != (pls1[10:] > 0).all()
    # A real, well-separated signal should explain a meaningful share of
    # variance -- catches the scale=True/scale=False bug (manually
    # confirmed against real data: that bug produced ratios on the order of
    # 1e-6 instead of comparable-to-PCA's ~0.7).
    assert expvar[0] > 0.1


def test_nmds_smoke_test_on_clustered_data():
    rng = np.random.RandomState(2)
    cluster_a = rng.normal(loc=0, scale=0.2, size=(6, 8))
    cluster_b = rng.normal(loc=10, scale=0.2, size=(6, 8))
    x = pd.DataFrame(
        np.vstack([cluster_a, cluster_b]),
        index=[f's{i}' for i in range(12)],
        columns=[f'f{i}' for i in range(8)],
    )
    scores, expvar, stress = run_nmds(x, n_components=2)
    assert scores.shape == (12, 2)
    assert len(expvar) == 2
    assert np.isfinite(stress)
    assert stress >= 0

    proxy = nmds_loading_proxy(x, scores)
    assert proxy.shape == (8, 2)
    assert proxy.values.min() >= -1.0001 and proxy.values.max() <= 1.0001


# --------------------------------------------------------------------------- #
# similarity_matrix
# --------------------------------------------------------------------------- #

def test_similarity_matrix_spearman_self_correlation_is_one():
    x = pd.DataFrame(
        [[1.0, 2.0, 3.0], [3.0, 2.0, 1.0], [1.0, 5.0, 2.0]],
        index=['s1', 's2', 's3'], columns=['f1', 'f2', 'f3'],
    )
    sim = similarity_matrix(x, 'Spearman')
    assert sim.shape == (3, 3)
    assert np.allclose(np.diag(sim.values), 1.0)
    # s1 and s2 are perfectly rank-anticorrelated.
    assert sim.loc['s1', 's2'] == pytest.approx(-1.0)


def test_similarity_matrix_jaccard_identical_presence_is_one():
    # s1/s2 detect exactly the same features (different abundances);
    # s3 detects none of them.
    x = pd.DataFrame(
        [[5.0, 0.0, 2.0], [50.0, 0.0, 20.0], [0.0, 0.0, 0.0]],
        index=['s1', 's2', 's3'], columns=['f1', 'f2', 'f3'],
    )
    sim = similarity_matrix(x, 'Jaccard')
    assert sim.loc['s1', 's2'] == pytest.approx(1.0)
    assert np.allclose(np.diag(sim.values)[:2], 1.0)


def test_similarity_matrix_braycurtis_identical_profiles_is_one():
    x = pd.DataFrame(
        [[1.0, 2.0, 3.0], [1.0, 2.0, 3.0], [10.0, 0.0, 0.0]],
        index=['s1', 's2', 's3'], columns=['f1', 'f2', 'f3'],
    )
    sim = similarity_matrix(x, 'Bray-Curtis')
    assert sim.loc['s1', 's2'] == pytest.approx(1.0)
    assert sim.loc['s1', 's3'] < sim.loc['s1', 's2']


def test_similarity_matrix_unknown_method_raises():
    x = pd.DataFrame([[1.0, 2.0]], index=['s1'], columns=['f1', 'f2'])
    with pytest.raises(ValueError):
        similarity_matrix(x, 'Pearson')


# --------------------------------------------------------------------------- #
# top_loadings
# --------------------------------------------------------------------------- #

def _make_loadings(n=30):
    rng = np.random.RandomState(3)
    return pd.DataFrame(
        rng.normal(size=(n, 2)),
        index=[f'feat{i}' for i in range(n)],
        columns=['PC1', 'PC2'],
    )


def test_top_loadings_returns_n_rows_by_default():
    loadings = _make_loadings(30)
    top = top_loadings(loadings, n=10)
    assert len(top) == 10


def test_top_loadings_includes_forced_feature_outside_top_n():
    loadings = _make_loadings(30)
    top = top_loadings(loadings, n=5)
    magnitude = np.sqrt((loadings ** 2).sum(axis=1))
    smallest_feature = magnitude.idxmin()
    assert smallest_feature not in top.index

    top_forced = top_loadings(loadings, n=5, always_include=[smallest_feature])
    assert len(top_forced) == 6
    assert smallest_feature in top_forced.index


def test_top_loadings_forced_feature_already_in_top_n_is_not_duplicated():
    loadings = _make_loadings(30)
    magnitude = np.sqrt((loadings ** 2).sum(axis=1))
    largest_feature = magnitude.idxmax()
    top = top_loadings(loadings, n=10, always_include=[largest_feature])
    assert len(top) == 10  # already in the top 10, no duplicate row added


def test_top_loadings_n_larger_than_available_returns_everything():
    loadings = _make_loadings(5)
    top = top_loadings(loadings, n=100)
    assert len(top) == 5
