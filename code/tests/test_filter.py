"""Unit tests for the pure-logic helpers in ``filter.py``.

These exercise the data-wrangling functions that don't need the Qt GUI or any
on-disk analysis outputs.
"""

import numpy as np
import pandas as pd

import filter


class _Params:
    """Minimal stand-in for the analysis_parameters object."""
    def __init__(self, blankfilthresh=0.5):
        self.blankfilthresh = blankfilthresh


def test_listfilter_include_keeps_listed_rows():
    df = pd.DataFrame({'Compound': ['a', 'b', 'c'], 'val': [1, 2, 3]})
    out = filter.listfilter(df, ['a', 'c'], True)
    assert out['Compound'].tolist() == ['a', 'c']


def test_listfilter_exclude_drops_listed_rows():
    df = pd.DataFrame({'Compound': ['a', 'b', 'c'], 'val': [1, 2, 3]})
    out = filter.listfilter(df, ['a', 'c'], False)
    assert out['Compound'].tolist() == ['b']


def test_listfilter_include_empty_list_yields_nothing():
    df = pd.DataFrame({'Compound': ['a', 'b'], 'val': [1, 2]})
    out = filter.listfilter(df, [], True)
    assert out.empty


def test_listfilter_exclude_empty_list_keeps_everything():
    df = pd.DataFrame({'Compound': ['a', 'b'], 'val': [1, 2]})
    out = filter.listfilter(df, [], False)
    assert out['Compound'].tolist() == ['a', 'b']


def test_reformatcorr_returns_lower_triangle_distances():
    # Symmetric correlation matrix; reformatcorr walks the sub-diagonal of each
    # column, clips negatives to zero, and returns 1 - value (a distance).
    corr = pd.DataFrame([[1.0, 0.5, -0.2],
                         [0.5, 1.0, 0.8],
                         [-0.2, 0.8, 1.0]])
    out = filter.reformatcorr(corr)
    # sub-diagonal column-major order: (1,0)=0.5, (2,0)=-0.2, (2,1)=0.8
    # clip negative -> 0, then 1 - x  =>  [0.5, 1.0, 0.2]
    np.testing.assert_allclose(out, np.array([0.5, 1.0, 0.2]))


def test_groupfilter_returns_ions_above_threshold():
    df = pd.DataFrame({'A': [0.2, 0.6, 0.9]}, index=['x', 'y', 'z'])
    out = filter.groupfilter('A', df, _Params(blankfilthresh=0.5))
    assert out == ['y', 'z']


def test_groupfilter_threshold_is_strict():
    df = pd.DataFrame({'A': [0.5, 0.5001]}, index=['x', 'y'])
    out = filter.groupfilter('A', df, _Params(blankfilthresh=0.5))
    assert out == ['y']
