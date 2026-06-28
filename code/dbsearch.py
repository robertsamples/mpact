"""
MPACT
Copyright 2022, Robert M. Samples, Sara P. Puckett, and Marcy J. Balunas

Qt-free Natural Products Atlas (NPAtlas) mass-matching search, extracted
from MainWindow.fulldbsearch so the ppm-window adduct matching logic can be
unit-tested without a GUI. fulldbsearch() was already entirely Qt-free
internally (no self.ui.* calls) -- this only changes its signature from
"take the whole window" to "take the specific inputs it needs," same as the
run_MSFaST(self) -> run_MSFaST(params) change.
"""

import pandas as pd

from csvcache import cached_read_csv


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
    iondict['hits'] = ''
    msdata = cached_read_csv(outputdir / (filename_stem + '_filtered.csv'),
                              sep=',', header=[2], index_col=None).iloc[:, :3]

    for _, mrow in msdata.iterrows():
        # Iterates over iondict, filters DB matches within window.
        # Repeats for adducts, uses length of concat DF for feature hits
        mass = mrow['m/z']
        hits_h = atlas[abs(1000000 * (atlas['compound_m_plus_h'] - mass) / atlas['compound_m_plus_h']) < ppm_threshold].copy()
        hits_h['ppm'] = abs(1000000 * (hits_h['compound_m_plus_h'] - mass) / hits_h['compound_m_plus_h'])
        hits_na = atlas[abs(1000000 * (atlas['compound_m_plus_na'] - mass) / atlas['compound_m_plus_na']) < ppm_threshold].copy()
        hits_na['ppm'] = abs(1000000 * (hits_na['compound_m_plus_na'] - mass) / hits_na['compound_m_plus_na'])
        hits = pd.concat([hits_h, hits_na])
        hits = hits.sort_values(by=['ppm'])
        hitdb[mrow['Compound']] = hits
        iondict.loc[mrow['Compound'], 'hits'] = hits.shape[0]

    iondict.to_csv(outputdir / 'iondict.csv', header=True, index=True)
    return hitdb, iondict
