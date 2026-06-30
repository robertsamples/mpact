"""
MPACT
Copyright 2022, Robert M. Samples, Sara P. Puckett, and Marcy J. Balunas
Modifications Copyright 2026, Robert M. Samples

Qt-free Natural Products Atlas (NPAtlas) mass-matching search, extracted
from MainWindow.fulldbsearch so the ppm-window adduct matching logic can be
unit-tested without a GUI. fulldbsearch() was already entirely Qt-free
internally (no self.ui.* calls) -- this only changes its signature from
"take the whole window" to "take the specific inputs it needs," same as the
run_MSFaST(self) -> run_MSFaST(params) change.
"""

import numpy as np

from csvcache import cached_read_csv, invalidate


def search_npatlas(outputdir, filename_stem, atlas, ppm_threshold):
    """Search every filtered feature's mass against NPAtlas for sodiated/
    protonated adduct matches within ``ppm_threshold``.

    Returns ``(hitdb, iondict)``: ``hitdb`` maps Compound -> a DataFrame of
    matching NPAtlas rows (sorted by ppm error); ``iondict`` is the updated
    ion dictionary with a ``hits`` count column, already written back to
    ``outputdir / 'iondict.csv'`` (matching the previous behaviour).
    """
    hitdb = {}
    iondict = cached_read_csv(outputdir / 'iondict.csv', sep=',', header=[0], index_col=[0])
    # NaN, not '' -- an empty-string default makes pandas infer a strict
    # string dtype for this column on newer pandas, which then rejects the
    # int assignment below (hits.shape[0]) outright. NaN keeps the column
    # numeric from the start (compatible with int assignment) and, for any
    # iondict row never touched by the loop below (a feature outside the
    # filtered/searched set), the comparison fillfttree() does later
    # (iondict['hits'] >= 0) correctly evaluates to False rather than
    # erroring on a string/int comparison.
    iondict['hits'] = np.nan
    msdata = cached_read_csv(outputdir / (filename_stem + '_filtered.csv'),
                              sep=',', header=[2], index_col=None).iloc[:, :3]

    # Pre-sort the two adduct-mass columns once so each feature only tests a
    # tiny m/z window (via searchsorted) instead of scanning all ~36k atlas
    # rows twice -- the old per-feature ``atlas[boolean_mask]`` over the whole
    # table was O(features x atlas_rows). The exact original ppm test
    # (``abs(1e6*(atlas_mz - mass)/atlas_mz) < ppm_threshold``) is re-applied to
    # the windowed candidates, so the matched set is bit-for-bit identical; the
    # window (mass*(1 +/- 2*t)) is a safe superset of the true ppm window for
    # the small tolerances used here. Verified output-identical (hitdb frames,
    # incl. row order + ppm, and the iondict 'hits' column) against the old
    # implementation on the real example dataset (~5x faster there).
    mph = atlas['compound_m_plus_h'].to_numpy(dtype=float)
    mna = atlas['compound_m_plus_na'].to_numpy(dtype=float)
    order_h = np.argsort(mph, kind='stable'); sorted_h = mph[order_h]
    order_na = np.argsort(mna, kind='stable'); sorted_na = mna[order_na]
    t = ppm_threshold / 1e6

    def _match(mass, sorted_vals, order, col_vals):
        # Atlas positions whose ppm error vs `mass` is below threshold, in
        # ascending atlas-position (i.e. original boolean-mask) order, plus
        # their ppm values.
        lo = np.searchsorted(sorted_vals, mass * (1 - 2 * t), side='left')
        hi = np.searchsorted(sorted_vals, mass * (1 + 2 * t), side='right')
        cand = order[lo:hi]
        if cand.size == 0:
            return cand, cand.astype(float)
        cv = col_vals[cand]
        sel = np.sort(cand[np.abs(1e6 * (cv - mass) / cv) < ppm_threshold])
        sel_cv = col_vals[sel]
        return sel, np.abs(1e6 * (sel_cv - mass) / sel_cv)

    masses = msdata['m/z'].to_numpy(dtype=float)
    compounds = msdata.iloc[:, 0].to_numpy()
    counts = np.empty(len(masses), dtype=float)
    for i in range(len(masses)):
        mass = masses[i]
        # m+h matches then m+na matches, concatenated in that order (matching
        # the old ``pd.concat([hits_h, hits_na])``) and slicing the atlas once.
        pos_h, ppm_h = _match(mass, sorted_h, order_h, mph)
        pos_na, ppm_na = _match(mass, sorted_na, order_na, mna)
        positions = np.concatenate([pos_h, pos_na])
        hits = atlas.iloc[positions].copy()
        hits['ppm'] = np.concatenate([ppm_h, ppm_na])
        hits = hits.sort_values(by=['ppm'])
        hitdb[compounds[i]] = hits
        counts[i] = positions.size
    # One vectorised column assignment instead of a per-feature ``.loc`` scalar
    # set (msdata's Compound ids are unique, a subset of iondict's index).
    iondict.loc[compounds, 'hits'] = counts

    iondict.to_csv(outputdir / 'iondict.csv', header=True, index=True)
    # iondict.csv just changed on disk (gained/updated the 'hits' column) --
    # drop every cached read of it, regardless of which header/index_col
    # shape it was cached under, so the next read (by anyone, any shape)
    # sees this update instead of a stale pre-search copy.
    invalidate()
    return hitdb, iondict
