"""Run the MPACT pipeline on the bundled PTY087I2 example with one code tree
and dump which features each filter removed, plus timings, as JSON.

Run once per code tree / variant, then compare the JSON files:

    python ab_pipeline.py <repo root> <out dir> <variant>

variant: default | mdfil | adduct | both | relfil_beta | all
(the non-default variants only mean something on the PR tree)
"""
import importlib.util
import json
import sys
import time
from pathlib import Path

repo = Path(sys.argv[1]).resolve()
out = Path(sys.argv[2]).resolve()
variant = sys.argv[3] if len(sys.argv) > 3 else 'default'

sys.path.insert(0, str(repo / 'code'))
spec = importlib.util.spec_from_file_location(
    'pipeline_fixture', repo / 'code' / 'tests' / 'test_msfast_pipeline.py')
fixture = importlib.util.module_from_spec(spec)
spec.loader.exec_module(fixture)

import filter as mpfilter  # noqa: E402
import pandas as pd  # noqa: E402
from MSFaST import run_MSFaST  # noqa: E402

out.mkdir(parents=True, exist_ok=True)
params = fixture.make_params(out)
if variant in ('mdfil', 'both', 'all'):
    params.mdfil = True
if variant in ('adduct', 'both', 'all'):
    params.adductgrp = True
    params.adduct_mode = 'pos'
if variant in ('relfil_beta', 'all'):
    params.relfil_mode = 'beta'

timings = {}


def timed(name):
    original = getattr(mpfilter, name, None)
    if original is None:
        return

    def wrapper(*a, **k):
        t0 = time.perf_counter()
        try:
            return original(*a, **k)
        finally:
            timings[name] = timings.get(name, 0.0) + time.perf_counter() - t0
    setattr(mpfilter, name, wrapper)


for fn in ('relationalfilter', 'mdfilter', 'adductgroup', 'decon', 'cvfilter', 'applyfilters'):
    timed(fn)

t0 = time.perf_counter()
result = run_MSFaST(params)
timings['total'] = time.perf_counter() - t0

stem = params.filename.stem
unfiltered = pd.read_csv(params.filename, sep=',', header=[0, 1, 2], index_col=[0, 1, 2])
filtered = pd.read_csv(params.outputdir / (stem + '_filtered.csv'),
                       sep=',', header=[0, 1, 2], index_col=[0, 1, 2])
iondict = pd.read_csv(params.outputdir / 'iondict.csv', index_col=0)

summary = {
    'variant': variant,
    'repo': str(repo),
    'n_input': int(len(unfiltered.index)),
    'n_filtered': int(len(filtered.index)),
    'kept': [str(i[0]) for i in filtered.index],
    'filters': {k: [str(x) for x in v.ions] for k, v in result.ionfilters.items()},
    'pass_columns': {c: int((iondict[c].astype(str).str.upper() == 'TRUE').sum())
                     for c in iondict.columns if str(c).startswith('pass')},
    'iondict_rows': int(len(iondict.index)),
    'timings_s': {k: round(v, 3) for k, v in timings.items()},
    'outputdir': str(params.outputdir),
}
(out / 'summary.json').write_text(json.dumps(summary, indent=1))
print(json.dumps({k: v for k, v in summary.items() if k not in ('kept', 'filters')}, indent=1))
print({k: len(v) for k, v in summary['filters'].items()})
