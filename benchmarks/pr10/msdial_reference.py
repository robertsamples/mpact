"""Benchmark the peak-correction and adduct code against MS-DIAL's own calls.

Uses the CHM350 Area_1 alignment (31,652 features, 48 injections), which exists
both as the raw MS-DIAL export (with MS-DIAL's isotope tracking and charge
assignment) and as the MPACT-layout table derived from it, row for row.

  1. Isotopes: MS-DIAL marks 13C isotopologue rows ("Isotope tracking weight
     number" > 0). Old and new relationalfilter are scored against that.
  2. Adducts: target/decoy. adductgroup() is run with the real adduct mass
     differences and then with shifted, chemically meaningless ones; decoy hits
     estimate how many "adducts" are chance matches.

    python msdial_reference.py <main repo> <PR repo> <CHM350_composite dir> <work dir>
"""
import importlib.util
import shutil
import sys
import time
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pandas as pd

main_repo, pr_repo, data_dir, work = (Path(a).resolve() for a in sys.argv[1:5])
STEM = 'Area_1_2024_11_19_19_38_54'


def load_filter(repo, name):
    sys.path.insert(0, str(repo / 'code'))
    spec = importlib.util.spec_from_file_location(name, repo / 'code' / 'filter.py')
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    sys.path.pop(0)
    return mod


old = load_filter(main_repo, 'filter_main')
new = load_filter(pr_repo, 'filter_pr')

# ---- reference: MS-DIAL ---------------------------------------------------
ref = pd.read_csv(data_dir / (STEM + '.txt'), sep='\t', skiprows=4, low_memory=False)
table = pd.read_csv(data_dir / (STEM + '.csv'), header=[0, 1, 2], index_col=None, low_memory=False)
table.columns = pd.MultiIndex.from_tuples(
    [('', '', c[2]) if i < 3 else c for i, c in enumerate(table.columns)])
ids = table.iloc[:, 0].astype(str).to_numpy()
mz = table.iloc[:, 1].to_numpy(float)
assert len(ref) == len(table) and np.allclose(mz, ref['Average Mz'].to_numpy(float), atol=1e-4), \
    'MPACT table and MS-DIAL export are not row-aligned'
assert len(set(ids)) == len(ids), 'duplicate compound IDs'

weight = ref['Isotope tracking weight number'].to_numpy(int)
parent = ref['Isotope tracking parent ID'].to_numpy(int)
iso = weight > 0
dm = mz[iso] - ref['Average Mz'].to_numpy(float)[parent[iso]]
per_step = dm / weight[iso]
print(f'MS-DIAL: {len(ref)} features, {iso.sum()} marked as isotopologues')
print(f'  sanity: mass step per isotope number, median {np.median(per_step):.4f} '
      f'(1.0034 = 13C, 0.5017 = doubly charged); '
      f'{np.mean(np.abs(per_step - 1.00336) < 0.01) * 100:.1f}% within 0.01 of 13C, '
      f'{np.mean(np.abs(per_step - 0.50168) < 0.01) * 100:.1f}% within 0.01 of 13C/2')
iso_ids = set(ids[iso])
doubly = set(ids[(ref['Adduct type'] == '[M+2H]2+').to_numpy()])


# ---- build the pipeline's intermediate files ------------------------------
def group_of(injection):
    name = injection.split('_', 1)[1]
    if name.lower().startswith('blank'):
        return 'Blank'
    return name.rsplit('_r', 1)[0]


def make_outdir(tag):
    out = work / tag
    if out.exists():
        shutil.rmtree(out)
    out.mkdir(parents=True)
    t = table.copy()
    t.columns = pd.MultiIndex.from_tuples(
        [(f'Unnamed: {i}_level_0', f'Unnamed: {i}_level_1', c[2]) if i < 3
         else (group_of(c[2]), group_of(c[2]), c[2]) for i, c in enumerate(t.columns)])
    t.to_csv(out / (STEM + '_formatted.csv'), index=False)
    iondict = pd.DataFrame({'Compound': ids, 'm/z': mz, 'Retention time (min)': table.iloc[:, 2].to_numpy(float)})
    iondict.to_csv(out / 'iondict.csv', index=False)
    return out


def params_for(out, rtwin, mode='legacy'):
    return SimpleNamespace(relfil_mode=mode, outputdir=out, filename=Path(STEM + '.csv'), RTwin=rtwin,
                           maxisowin=3, ringingwin=0.5, isopeakwin=0.01, dimerpeakwin=0.01,
                           adduct_mode='pos')


def score(label, removed, seconds):
    removed = set(map(str, removed))
    tp = len(removed & iso_ids)
    prec = tp / len(removed) if removed else float('nan')
    print(f'  {label:26} removed {len(removed):5d} | isotopologues caught {tp:4d}/{len(iso_ids)} '
          f'(recall {100 * tp / len(iso_ids):5.1f}%) | share of removals that are MS-DIAL isotopes '
          f'{100 * prec:5.1f}% | {seconds:6.1f}s')
    return removed


print('\n== 1. relationalfilter vs MS-DIAL isotope tracking ==')
for rtwin in (0.005, 0.02):
    print(f' RT window {rtwin} min')
    for label, mod, mode in (('main (published)', old, 'legacy'), ('PR, legacy mode', new, 'legacy'), ('PR, beta mode', new, 'beta')):
        out = make_outdir(f'relfil_{label.split()[0]}_{mode}_{rtwin}')
        t0 = time.perf_counter()
        res = mod.relationalfilter(params_for(out, rtwin, mode), {})
        score(label, res['relfil'].ions, time.perf_counter() - t0)

# ---- 2. adduct grouping: target vs decoy ----------------------------------
print('\n== 2. adductgroup target/decoy (RT window 0.02 min, 3 ppm) ==')
H, NA, K, NH4, H2O = 1.007276, 22.989770, 38.963707, 18.033823, 18.010565
TARGET = [('M+Na', NA - H), ('M+K', K - H), ('M+NH4', NH4 - H), ('M-H2O+H', -H2O), ('M-H4O2+H', -2 * H2O)]


def run_adduct(tag, adduct_list):
    out = make_outdir(tag)
    p = params_for(out, 0.02)
    t0 = time.perf_counter()
    new.adductgroup(p, adduct_list)
    seconds = time.perf_counter() - t0
    iondict = pd.read_csv(out / 'iondict.csv', index_col=0)
    forms = iondict['adduct_form'].fillna('').astype(str)
    counts = forms[forms.isin([n for n, _ in adduct_list])].value_counts()
    return counts.reindex([n for n, _ in adduct_list]).fillna(0).astype(int), seconds, iondict


target_counts, seconds, target_iondict = run_adduct('adduct_target', TARGET)
print(f'  target: {int(target_counts.sum())} ions assigned an adduct form in {seconds:.1f}s')

rng = np.random.default_rng(0)
decoys = []
for k in range(10):
    shift = rng.uniform(0.15, 0.45, size=len(TARGET)) * rng.choice([-1, 1], size=len(TARGET))
    counts, _, _ = run_adduct('adduct_decoy', [(n, d + s) for (n, d), s in zip(TARGET, shift)])
    decoys.append(counts.to_numpy())
decoys = np.array(decoys)
print(f'  {"form":10} {"target":>7} {"decoy mean":>11} {"decoy sd":>9} {"est. FDR":>9}')
for i, (name, _) in enumerate(TARGET):
    t, d = target_counts.iloc[i], decoys[:, i].mean()
    print(f'  {name:10} {t:7d} {d:11.1f} {decoys[:, i].std():9.1f} {100 * min(1, d / t) if t else float("nan"):8.1f}%')
tot_t, tot_d = target_counts.sum(), decoys.sum(axis=1).mean()
print(f'  {"all":10} {tot_t:7d} {tot_d:11.1f} {decoys.sum(axis=1).std():9.1f} {100 * tot_d / tot_t:8.1f}%')

# how often does an "adduct" call land on something MS-DIAL says is an isotope or 2+?
forms = target_iondict['adduct_form'].fillna('').astype(str)
called = set(forms[forms.isin([n for n, _ in TARGET])].index.astype(str))
print(f'  adduct calls that MS-DIAL marks as isotopologues: {len(called & iso_ids)}, as [M+2H]2+: {len(called & doubly)}')
