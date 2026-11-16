"""What does the peak-correction filter actually remove on PTY087I2?

Runs relationalfilter from one code tree on the example dataset, then describes
every (parent, removed feature) pair: mass difference, intensity ratio, and
whether the mass-defect filter would have removed the feature anyway.

    python relfil_pairs.py <repo to take relationalfilter from> <repo with md_plausible> <work dir> [mode]

mode is passed as analysis_params.relfil_mode (legacy | beta) when given.
"""
import importlib.util
import sys
from pathlib import Path

import numpy as np
import pandas as pd

repo, pr_repo, work = (Path(a).resolve() for a in sys.argv[1:4])
mode = sys.argv[4] if len(sys.argv) > 4 else None

sys.path.insert(0, str(repo / 'code'))
spec = importlib.util.spec_from_file_location(
    'pipeline_fixture', repo / 'code' / 'tests' / 'test_msfast_pipeline.py')
fixture = importlib.util.module_from_spec(spec)
spec.loader.exec_module(fixture)
import filter as mpfilter  # noqa: E402
from MSFaST import run_MSFaST  # noqa: E402

spec = importlib.util.spec_from_file_location('filter_pr', pr_repo / 'code' / 'filter.py')
prfilter = importlib.util.module_from_spec(spec)
spec.loader.exec_module(prfilter)

import shutil  # noqa: E402
if work.exists():
    shutil.rmtree(work)
work.mkdir(parents=True)
params = fixture.make_params(work)
if mode:
    params.relfil_mode = mode

captured = {}
original = mpfilter.relationalfilter


def capture(analysis_params, ionfilterlist):
    f = analysis_params.outputdir / (analysis_params.filename.stem + '_formatted.csv')
    captured['table'] = pd.read_csv(f, header=[0, 1, 2], index_col=0, low_memory=False)
    out = original(analysis_params, ionfilterlist)
    captured['merge'] = {k: list(np.atleast_1d(v.sources)) for k, v in out['relfil'].merge.items()}
    return out


mpfilter.relationalfilter = capture
run_MSFaST(params)

t = captured['table']
t.index = t.index.astype(str)
mz = pd.to_numeric(t.iloc[:, 0]).to_dict()
keep = [c for c in t.columns[2:] if 'BLANK' not in str(c[0]).upper()]
mx = t[keep].apply(pd.to_numeric, errors='coerce').fillna(0).max(axis=1).to_dict()

rows = []
for parent, kids in captured['merge'].items():
    for kid in kids:
        p, k = str(parent), str(kid)
        rows.append((k, mz[k] - mz[p], mx[k] / mx[p] if mx[p] > 0 else np.inf, mz[k]))
d = pd.DataFrame(rows, columns=['removed', 'dm', 'ratio', 'mz'])
frac = d['dm'] - np.floor(d['dm'])
n13 = np.round(d['dm'] / 1.00335)
iso = (n13 >= 1) & ((d['dm'] - n13 * 1.00335).abs() < 0.01)
half = (frac - 0.5017).abs() < 0.01
near = d['dm'] < 0.5
d['kind'] = np.select([iso, half, near], ['13C isotope spacing (n x 1.0034 +/- 0.01)',
                                          'half-integer spacing (2+ ion)',
                                          'within 0.5 Da of parent (ringing band)'],
                      'other spacing in 1.0-1.5 / 2.0-2.5 Da band')
d['md_fail'] = ~prfilter.md_plausible(d['mz'].to_numpy(), 0.00112, 0.03, -0.00025, -0.03)

print(f'\n{len(d)} features removed')
print(f'{"kind":46} {"n":>4} {"weaker than parent":>19} {"<1% of parent":>14} {"fails mass defect":>18}')
for kind, g in d.groupby('kind'):
    print(f'{kind:46} {len(g):4d} {int((g.ratio < 1).sum()):19d} {int((g.ratio < 0.01).sum()):14d} {int(g.md_fail.sum()):18d}')
print(f'{"all":46} {len(d):4d} {int((d.ratio < 1).sum()):19d} {int((d.ratio < 0.01).sum()):14d} {int(d.md_fail.sum()):18d}')
print('removed feature STRONGER than the parent it was merged into:', int((d.ratio > 1).sum()))
print('intensity ratio removed/parent, quartiles:', np.percentile(d.ratio.replace(np.inf, 1e9), [10, 25, 50, 75, 90]).round(3))
