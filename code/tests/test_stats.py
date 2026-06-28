"""Unit tests for the statistical helpers in ``stats.py``.

``tstat`` and ``ws`` operate elementwise on pandas Series (they call
``.replace`` to scrub inf/nan), so the tests feed single-element Series and
check against hand-computed values.
"""

import numpy as np
import pandas as pd

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
