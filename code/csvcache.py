"""
MPACT
Copyright 2022, Robert M. Samples, Sara P. Puckett, and Marcy J. Balunas

Memoizing cache for repeated ``pd.read_csv`` calls against the static,
per-run output files (``iondict.csv``, ``<name>_filtered.csv``,
``<name>_summarydata.csv``) that several interactive code paths -- feature
selection/highlighting, heatmap W/S navigation, the abundance plot -- were
re-reading from disk on every single click/keypress, even though the
underlying file doesn't change between analysis runs.

Caches the *exact* ``pd.read_csv`` call (path + kwargs), not a derived
"view" of one universal parse: several call sites read the same file with
genuinely different ``header``/``index_col`` shapes (e.g. ``_filtered.csv``
is read with ``header=[2]`` in one place and ``header=[0, 1, 2]`` in
another) -- those are different parses, not different slices of the same
dataframe, and deriving one from the other by hand would risk subtly wrong
dtypes/column structure with no way to catch it without running the live
app against real data. Keying on the exact call signature instead means
every call site still gets back something equivalent to what
``pd.read_csv`` would have returned, just served from memory after the
first call instead of re-read from disk every time.
"""

import pandas as pd

_cache = {}


def _hashable(value):
    if isinstance(value, list):
        return tuple(_hashable(v) for v in value)
    return value


def cached_read_csv(path, **kwargs):
    """``pd.read_csv(path, **kwargs)``, memoized by ``(path, kwargs)``.

    Returns a fresh ``.copy()`` each call so callers can slice/mutate their
    own result without affecting other callers or the cached copy --
    matching the previous behaviour where every call site got its own
    independent DataFrame straight from disk.
    """
    key = (str(path), tuple(sorted((k, _hashable(v)) for k, v in kwargs.items())))
    if key not in _cache:
        _cache[key] = pd.read_csv(path, **kwargs)
    return _cache[key].copy()


def invalidate():
    """Drop every cached read.

    Call this whenever a new analysis run completes or a different
    session/dataset is loaded -- the underlying files may now be different
    even though their paths could coincidentally repeat (e.g. ``iondict.csv``
    is always named the same within its own output directory).
    """
    _cache.clear()
