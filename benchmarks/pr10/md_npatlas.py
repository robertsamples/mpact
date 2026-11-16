"""False-rejection rate of the PR's mass-defect filter on known natural products.

Every NPAtlas compound is a real, published microbial natural product, so any
[M+H]+ (or [M+Na]+, [M-H]-) m/z that md_plausible() rejects is a real compound
the filter would delete.

    python md_npatlas.py <repo root with the PR> <npatlas.tsv>
"""
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd

repo, npatlas = Path(sys.argv[1]).resolve(), Path(sys.argv[2])
sys.path.insert(0, str(repo / 'code'))
import filter as mpfilter  # noqa: E402

DEFAULT = dict(hi_slope=0.00112, hi_int=0.03, lo_slope=-0.00025, lo_int=-0.03)
HALOGEN = dict(hi_slope=0.00112, hi_int=0.03, lo_slope=-0.0011, lo_int=-0.05)
PROTON, NA = 1.007276, 22.989218

df = pd.read_csv(npatlas, sep='\t', low_memory=False,
                 usecols=['compound_molecular_formula', 'compound_accurate_mass'])
df = df.dropna()
mass = df['compound_accurate_mass'].to_numpy(float)
formula = df['compound_molecular_formula'].astype(str)
halogen = formula.str.contains(r'Cl|Br|I(?![a-z])', regex=True).to_numpy()
elements = formula.map(lambda f: frozenset(re.findall(r'[A-Z][a-z]?', f)))
chnops = elements.map(lambda e: e <= {'C', 'H', 'N', 'O', 'P', 'S'}).to_numpy()

ions = {'[M+H]+': mass + PROTON, '[M+Na]+': mass + NA, '[M-H]-': mass - PROTON,
        '[M+2H]2+': (mass + 2 * PROTON) / 2}

print(f'NPAtlas compounds: {len(mass)}  (halogenated: {halogen.sum()}, CHNOPS-only: {chnops.sum()})')
print(f'{"ion":10} {"bounds":8} {"rejected":>9} {"%all":>7} {"%CHNOPS-only":>13} {"%halogenated":>13}')
for ion, mz in ions.items():
    for label, b in (('default', DEFAULT), ('halogen', HALOGEN)):
        ok = mpfilter.md_plausible(mz, **b)
        rej = ~ok
        print(f'{ion:10} {label:8} {rej.sum():9d} {100 * rej.mean():6.2f}% '
              f'{100 * rej[chnops].mean():12.2f}% {100 * rej[halogen].mean():12.2f}%')

print('\n[M+H]+ rejection by mass range, default bounds:')
mz = ions['[M+H]+']
rej = ~mpfilter.md_plausible(mz, **DEFAULT)
for lo, hi in ((0, 200), (200, 300), (300, 400), (400, 500), (500, 600), (600, 700), (700, 5000)):
    m = (mz >= lo) & (mz < hi)
    print(f'  {lo:4d}-{hi:<5d} n={m.sum():6d}  rejected {100 * rej[m].mean():5.2f}%')

width_one = (1 - (DEFAULT['hi_int'] - DEFAULT['lo_int'])) / (DEFAULT['hi_slope'] - DEFAULT['lo_slope'])
print(f'\nFilter can reject nothing above m/z {width_one:.0f} (allowed window >= 1 Da wide).')
