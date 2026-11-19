"""Unit tests for the statistical helpers in ``stats.py``.

``tstat`` and ``ws`` operate elementwise on pandas Series (they call
``.replace`` to scrub inf/nan), so the tests feed single-element Series and
check against hand-computed values.
"""

from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pandas as pd
import pytest

import stats


def _s(value):
    return pd.Series([float(value)])


def test_tstat_matches_welch_numerator_over_pooled_se():
    # m1=10, m2=8, sd=2, n=4 for both groups.
    # den1 = sd^2/n = 1, den2 = 1, denom = sqrt(2); t = |10-8| / sqrt(2).
    out = stats.tstat(_s(10), _s(2), _s(4), _s(8), _s(2), _s(4))
    np.testing.assert_allclose(out.values, [2.0 / np.sqrt(2)])


def test_tstat_is_symmetric_in_group_order():
    a = stats.tstat(_s(10), _s(2), _s(4), _s(8), _s(2), _s(4))
    b = stats.tstat(_s(8), _s(2), _s(4), _s(10), _s(2), _s(4))
    np.testing.assert_allclose(a.values, b.values)


def test_ws_effective_sample_size_for_equal_groups():
    # sd=2, n=4 for both groups -> Welch-Satterthwaite effective n = 6.
    out = stats.ws(_s(2), _s(4), _s(2), _s(4))
    np.testing.assert_allclose(out.values, [6.0])


def _groupave_by_stacking(formatted_path):
    """The previous groupave, unchunked: one row per value, then groupby."""
    data = pd.read_csv(formatted_path, header=[0, 1, 2], index_col=[0, 1, 2])
    levels = ['Compound', 'm/z', 'Retention time', 'Group', 'Sample', 'Injection']
    long = data.stack(level=[0, 1, 2], future_stack=True)
    long.index.names = levels
    mean_values = long.groupby(level=levels).sum() / long.groupby(level=levels).count()
    tech = levels[:5]
    tech_mean = mean_values.groupby(level=tech).mean()
    tech_rsd = (mean_values.groupby(level=tech).std().fillna(0) / tech_mean).fillna(0)
    techn = mean_values.groupby(level=tech).size()
    biol = levels[:4]
    biol_mean = tech_mean.groupby(level=biol).mean()
    return pd.DataFrame({
        'average': biol_mean,
        'biolRSD': (tech_mean.groupby(level=biol).std().fillna(0) / biol_mean).fillna(0),
        'bioln': tech_mean.groupby(level=biol).size(),
        'techRSD': tech_rsd.groupby(level=biol).mean(),
        'techn': techn.groupby(level=biol).mean(),
    })


def test_groupave_matches_the_stacked_calculation(tmp_path):
    # Group A: two samples with 3 and 2 injections. Group B: one sample.
    # Includes a missing value and an all-zero feature.
    injections = [('A', 's1', 'i1'), ('A', 's1', 'i2'), ('A', 's1', 'i3'),
                  ('A', 's2', 'i4'), ('A', 's2', 'i5'), ('B', 's3', 'i6')]
    columns = pd.MultiIndex.from_tuples(
        [(f'Unnamed: {i}_level_0', f'Unnamed: {i}_level_1', name)
         for i, name in enumerate(['Compound', 'm/z', 'Retention time (min)'])]
        + injections)
    rows = [
        ['f1', 100.1, 1.0, 10.0, 12.0, 11.0, 20.0, 22.0, 5.0],
        ['f2', 200.2, 2.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
        ['f3', 300.3, 3.0, 7.0, np.nan, 9.0, 1.0, 3.0, 4.0],
    ]
    pd.DataFrame(rows, columns=columns).to_csv(tmp_path / 'x_formatted.csv', index=False)
    stats.groupave(SimpleNamespace(outputdir=tmp_path, filename=Path('x.csv')))

    written = pd.read_csv(tmp_path / 'x_summarydata.csv', index_col=[0, 1, 2, 3])
    expected = _groupave_by_stacking(tmp_path / 'x_formatted.csv')
    pd.testing.assert_frame_equal(written, expected, check_dtype=False, check_names=False)
    assert written.loc[('f1', 100.1, 1.0, 'A'), 'average'] == pytest.approx((11.0 + 21.0) / 2)
    assert written.loc[('f1', 100.1, 1.0, 'A'), 'techn'] == 2.5
    assert written.loc[('f1', 100.1, 1.0, 'B'), 'bioln'] == 1
