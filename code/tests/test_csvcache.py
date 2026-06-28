from unittest.mock import patch

import pandas as pd

import csvcache
from csvcache import cached_read_csv, invalidate


def setup_function(_):
    # Each test gets a clean cache regardless of execution order.
    invalidate()


def write_csv(tmp_path, name='data.csv'):
    path = tmp_path / name
    pd.DataFrame({'a': [1, 2], 'b': [3, 4]}).to_csv(path, index=False)
    return path


def test_second_call_does_not_hit_disk(tmp_path):
    path = write_csv(tmp_path)
    with patch.object(csvcache.pd, 'read_csv', wraps=pd.read_csv) as spy:
        first = cached_read_csv(path, sep=',')
        second = cached_read_csv(path, sep=',')
        assert spy.call_count == 1
    pd.testing.assert_frame_equal(first, second)


def test_different_kwargs_are_cached_independently(tmp_path):
    path = write_csv(tmp_path)
    with patch.object(csvcache.pd, 'read_csv', wraps=pd.read_csv) as spy:
        cached_read_csv(path, sep=',', header=0)
        cached_read_csv(path, sep=',', header=None)
        assert spy.call_count == 2
        # repeating either exact signature again should not hit disk
        cached_read_csv(path, sep=',', header=0)
        cached_read_csv(path, sep=',', header=None)
        assert spy.call_count == 2


def test_list_valued_kwargs_are_hashable_in_the_cache_key(tmp_path):
    path = write_csv(tmp_path)
    # header=[0] previously would have failed: lists aren't hashable, so a
    # naive tuple(sorted(kwargs.items())) cache key would raise TypeError.
    df = cached_read_csv(path, sep=',', header=[0], index_col=[0])
    assert list(df.columns) == ['b']


def test_returned_frame_is_a_copy_not_shared_state(tmp_path):
    path = write_csv(tmp_path)
    first = cached_read_csv(path, sep=',')
    first.iloc[0, 0] = 999
    second = cached_read_csv(path, sep=',')
    assert second.iloc[0, 0] != 999


def test_invalidate_forces_a_fresh_read(tmp_path):
    path = write_csv(tmp_path)
    with patch.object(csvcache.pd, 'read_csv', wraps=pd.read_csv) as spy:
        cached_read_csv(path, sep=',')
        invalidate()
        cached_read_csv(path, sep=',')
        assert spy.call_count == 2
