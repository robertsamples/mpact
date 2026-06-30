import pandas as pd

from csvcache import cached_read_csv
from dbsearch import search_npatlas


def make_atlas():
    # Two NPAtlas-shaped rows: one whose [M+H] mass is close to our test
    # feature, one whose [M+Na] mass is close, one that's nowhere near.
    return pd.DataFrame({
        'compound_name': ['near_h', 'near_na', 'far_off'],
        'compound_m_plus_h': [101.0, 500.0, 900.0],
        'compound_m_plus_na': [600.0, 123.0, 950.0],
    })


def setup_files(tmp_path):
    outputdir = tmp_path
    stem = 'example'

    iondict = pd.DataFrame({'Compound': ['c1', 'c2'], 'other': [1, 2]})
    iondict.to_csv(outputdir / 'iondict.csv', index=False)

    # _filtered.csv read with header=[2] -- mirror the canonical 3-header-row
    # shape, with the third row as the usable single-level header.
    filtered_path = outputdir / (stem + '_filtered.csv')
    with open(filtered_path, 'w') as f:
        f.write(',,\n')
        f.write(',,\n')
        f.write('Compound,m/z,Retention time (min)\n')
        f.write('c1,101.0005,1.0\n')   # matches near_h within ppm
        f.write('c2,123.0002,2.0\n')   # matches near_na within ppm

    return outputdir, stem


def test_matches_found_within_ppm_window(tmp_path):
    outputdir, stem = setup_files(tmp_path)
    atlas = make_atlas()

    hitdb, iondict = search_npatlas(outputdir, stem, atlas, ppm_threshold=50)

    assert 'near_h' in hitdb['c1']['compound_name'].values
    assert 'far_off' not in hitdb['c1']['compound_name'].values
    assert 'near_na' in hitdb['c2']['compound_name'].values


def test_hit_counts_written_back_to_iondict_csv(tmp_path):
    outputdir, stem = setup_files(tmp_path)
    atlas = make_atlas()

    search_npatlas(outputdir, stem, atlas, ppm_threshold=50)

    on_disk = pd.read_csv(outputdir / 'iondict.csv', index_col=0)
    assert on_disk.loc['c1', 'hits'] == 1
    assert on_disk.loc['c2', 'hits'] == 1


def test_no_match_outside_ppm_window(tmp_path):
    outputdir, stem = setup_files(tmp_path)
    atlas = make_atlas()

    hitdb, _ = search_npatlas(outputdir, stem, atlas, ppm_threshold=0.001)

    assert hitdb['c1'].empty
    assert hitdb['c2'].empty


def _write_single_feature(tmp_path, mz):
    """One-feature filtered table + matching iondict, for the ordering tests."""
    stem = 'example'
    pd.DataFrame({'Compound': ['feat'], 'other': [1]}).to_csv(tmp_path / 'iondict.csv', index=False)
    with open(tmp_path / (stem + '_filtered.csv'), 'w') as f:
        f.write(',,\n,,\nCompound,m/z,Retention time (min)\n')
        f.write('feat,%s,1.0\n' % mz)
    return tmp_path, stem


def test_hits_sorted_by_ppm_across_both_adducts(tmp_path):
    # m+h matches A (~20 ppm) and B (~30 ppm); m+na matches C (~5 ppm).
    # The combined result must be ascending by ppm: C, A, B.
    atlas = pd.DataFrame({
        'compound_name': ['A', 'B', 'C'],
        'compound_m_plus_h': [200.000, 200.010, 999.0],
        'compound_m_plus_na': [999.0, 999.0, 200.005],
    })
    outputdir, stem = _write_single_feature(tmp_path, 200.004)
    hitdb, _ = search_npatlas(outputdir, stem, atlas, ppm_threshold=100)
    hits = hitdb['feat']
    assert list(hits['compound_name']) == ['C', 'A', 'B']
    assert list(hits['ppm']) == sorted(hits['ppm'])  # ascending


def test_single_atlas_row_matching_both_adducts_appears_twice(tmp_path):
    # One atlas row whose [M+H] and [M+Na] are both at the feature mass must
    # appear once per adduct (two rows), matching the old concat behaviour.
    atlas = pd.DataFrame({
        'compound_name': ['D'],
        'compound_m_plus_h': [300.000],
        'compound_m_plus_na': [300.000],
    })
    outputdir, stem = _write_single_feature(tmp_path, 300.000)
    hitdb, _ = search_npatlas(outputdir, stem, atlas, ppm_threshold=10)
    assert len(hitdb['feat']) == 2
    assert list(hitdb['feat']['compound_name']) == ['D', 'D']
    on_disk = pd.read_csv(outputdir / 'iondict.csv', index_col=0)
    assert on_disk.loc['feat', 'hits'] == 2


def test_invalidates_stale_cached_reads_under_other_shapes(tmp_path):
    """Regression guard: fillfttree() (main.py) reads iondict.csv with
    header=[0], index_col=None -- a different cache key than
    search_npatlas's own internal read (index_col=[0]). Before
    search_npatlas() ran, something else (e.g. _finish_analysis) may have
    already cached that other shape *without* the 'hits' column. If
    search_npatlas doesn't invalidate the whole cache after writing, that
    other shape keeps serving the pre-search copy forever -- exactly the
    KeyError: 'hits' crash this guards against.
    """
    outputdir, stem = setup_files(tmp_path)
    atlas = make_atlas()

    # Simulate _finish_analysis's earlier read, before any hits exist,
    # under the *other* shape fillfttree() will later use.
    stale = cached_read_csv(outputdir / 'iondict.csv', sep=',', header=[0], index_col=None)
    assert 'hits' not in stale.columns

    search_npatlas(outputdir, stem, atlas, ppm_threshold=50)

    fresh = cached_read_csv(outputdir / 'iondict.csv', sep=',', header=[0], index_col=None)
    assert 'hits' in fresh.columns
