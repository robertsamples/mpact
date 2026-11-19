"""Tests for ``mzmineimport.rename_duplicates``."""

import pandas as pd

import mzmineimport


def _write_table(path, ids):
    columns = pd.MultiIndex.from_tuples(
        [(f'Unnamed: {i}_level_0', f'Unnamed: {i}_level_1', name)
         for i, name in enumerate(['Compound', 'm/z', 'Retention time (min)'])]
        + [('A', 's1', 'i1')])
    rows = [[cid, 100.0 + n, 1.0, 10.0] for n, cid in enumerate(ids)]
    pd.DataFrame(rows, columns=columns).to_csv(path, index=False)


def test_rename_duplicates_with_text_ids(tmp_path):
    path = tmp_path / 'x.csv'
    _write_table(path, ['a', 'b', 'a'])
    mzmineimport.rename_duplicates(path)
    ids = pd.read_csv(path, header=[0, 1, 2], index_col=0).index
    assert list(ids) == ['a', 'b', 'a1']


def test_rename_duplicates_with_numeric_ids(tmp_path):
    # MS-DIAL alignment IDs read back as integers; adding the text suffix
    # used to raise TypeError even when nothing was duplicated.
    path = tmp_path / 'x.csv'
    _write_table(path, [7, 8, 7])
    mzmineimport.rename_duplicates(path)
    ids = pd.read_csv(path, header=[0, 1, 2], index_col=0).index
    assert [str(i) for i in ids] == ['7', '8', '71']
