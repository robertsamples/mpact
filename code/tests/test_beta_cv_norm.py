"""Tests for the beta CV sample cutoff (``filter.cvfilter`` with ``cvgate``)
and the beta sample normalization (``stats.samplenormalize``)."""

from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pandas as pd
import pytest

import filter
import stats
from MSFaST import run_MSFaST
from test_msfast_pipeline import make_params


def write_formatted(tmp_path, features, injections):
    """features: [(id, m/z, RT, [one value per injection]), ...];
    injections: [(group, sample, injection), ...]."""
    columns = pd.MultiIndex.from_tuples(
        [(f'Unnamed: {i}_level_0', f'Unnamed: {i}_level_1', name)
         for i, name in enumerate(['Compound', 'm/z', 'Retention time (min)'])]
        + injections)
    rows = [[fid, mz, rt, *values] for fid, mz, rt, values in features]
    pd.DataFrame(rows, columns=columns).to_csv(tmp_path / 'x_formatted.csv', index=False)
    pd.DataFrame({'Compound': [f[0] for f in features], 'm/z': [f[1] for f in features],
                  'Retention time (min)': [f[2] for f in features]}).to_csv(
        tmp_path / 'iondict.csv', index=False)
    return SimpleNamespace(outputdir=tmp_path, filename=Path('x.csv'), cvthresh=0.2,
                           cvparam='median CV')


# Six samples with two injections each.
INJECTIONS = [(g, f's{n}', f's{n}_{k}') for n, g in enumerate('AABBCC') for k in (1, 2)]
# Real in two samples (CV 0.05); near noise level in the other four (CV ~0.9).
SPARSE = ('sparse', 300.1, 1.0,
          [1e7, 1.07e7, 1e7, 1.07e7, 50, 500, 40, 600, 100, 900, 30, 400])
STEADY = ('steady', 200.1, 2.0, [1e6, 1.05e6] * 6)


def test_cv_cutoff_off_keeps_the_published_result(tmp_path):
    params = write_formatted(tmp_path, [SPARSE, STEADY], INJECTIONS)
    failed = filter.cvfilter(params, {}, 0.2)['cv'].ions
    assert set(failed) == {'sparse'}
    assert 'CV samples used' not in pd.read_csv(tmp_path / 'iondict.csv').columns


def test_cv_cutoff_leaves_out_noise_level_samples(tmp_path):
    params = write_formatted(tmp_path, [SPARSE, STEADY], INJECTIONS)
    params.cvgate, params.cvgate_frac = True, 0.01
    failed = filter.cvfilter(params, {}, 0.2)['cv'].ions
    assert failed == []
    iondict = pd.read_csv(tmp_path / 'iondict.csv', index_col=0)
    assert iondict.loc['sparse', 'CV samples used'] == 2
    assert iondict.loc['steady', 'CV samples used'] == 6
    assert iondict.loc['sparse', 'median CV'] == pytest.approx(0.07 / 1.035 / np.sqrt(2), rel=1e-3)


def wide(values, injections):
    return pd.DataFrame(np.array(values, dtype=float),
                        columns=pd.MultiIndex.from_tuples(injections))


def test_factors_recover_a_dilution():
    rng = np.random.default_rng(1)
    base = rng.lognormal(10, 2, size=(200, 1)) * np.ones((1, 12))
    dilution = np.repeat([1, 1, 0.25, 0.25, 1, 1], 2)
    data = wide(base * dilution, INJECTIONS)
    for method in ('sum', 'pqn'):
        factors = stats.normalization_factors(data, method)
        assert factors.loc[('B', 's2'), 'relative signal'] == pytest.approx(0.25)
        assert factors.loc[('B', 's2'), 'factor'] == pytest.approx(4.0)
        assert factors.loc[('A', 's0'), 'factor'] == pytest.approx(1.0)
        assert (factors['injections'] == 2).all()


def test_pqn_is_not_fooled_by_a_few_large_differences():
    rng = np.random.default_rng(2)
    base = rng.lognormal(10, 1, size=(200, 1)) * np.ones((1, 12))
    base[:10, 8:] *= 50    # group C has ten compounds the others barely have
    data = wide(base, INJECTIONS)
    pqn = stats.normalization_factors(data, 'pqn')
    total = stats.normalization_factors(data, 'sum')
    assert pqn.loc[('C', 's4'), 'relative signal'] == pytest.approx(1.0)
    assert total.loc[('C', 's4'), 'relative signal'] > 2


def test_skipped_groups_keep_factor_one():
    data = wide(np.ones((5, 12)) * np.repeat([1, 1, 2, 2, 0.001, 0.001], 2), INJECTIONS)
    factors = stats.normalization_factors(data, 'sum', skip_groups=['C'])
    assert factors.loc[('C', 's4'), 'factor'] == 1.0
    assert np.isnan(factors.loc[('C', 's4'), 'relative signal'])
    # Totals 5, 5, 10, 10 (C left out): the median sample is 7.5.
    assert factors.loc[('B', 's2'), 'factor'] == pytest.approx(0.75)


def test_samplenormalize_scales_tables_and_keeps_the_raw_copy(tmp_path):
    rows = [(f'f{i}', 100.0 + i, 1.0, list(np.repeat([1, 1, 4, 4, 1, 1], 2) * (i + 1) * 1000.0))
            for i in range(20)]
    params = write_formatted(tmp_path, rows, INJECTIONS)
    (tmp_path / 'x_filtered.csv').write_bytes((tmp_path / 'x_formatted.csv').read_bytes())
    params.samplenorm, params.blnkgrp = 'pqn', ''

    factors = stats.samplenormalize(params)

    assert factors.loc[('B', 's2'), 'factor'] == pytest.approx(0.25)
    raw = pd.read_csv(tmp_path / 'x_filtered_raw.csv', header=[0, 1, 2], index_col=[0, 1, 2])
    scaled = pd.read_csv(tmp_path / 'x_filtered.csv', header=[0, 1, 2], index_col=[0, 1, 2])
    assert raw.iloc[0, 4] == 4000 and scaled.iloc[0, 4] == pytest.approx(1000)
    assert np.allclose(scaled.values, np.repeat(scaled.values[:, :1], 12, axis=1))
    formatted = pd.read_csv(tmp_path / 'x_formatted.csv', header=[0, 1, 2], index_col=[0, 1, 2])
    assert np.allclose(formatted.values, scaled.values)
    saved = pd.read_csv(tmp_path / 'x_normfactors.csv', index_col=[0, 1])
    assert list(saved.columns) == ['injections', 'relative signal', 'factor']


def test_samplenormalize_leaves_chosen_groups_out_and_warns(tmp_path, capsys):
    # C is a media blank with a trace of everything: far below the rest.
    rows = [(f'f{i}', 100.0 + i, 1.0, list(np.repeat([1, 1, 1.2, 1.2, 0.01, 0.01], 2) * (i + 1) * 1000.0))
            for i in range(20)]
    params = write_formatted(tmp_path, rows, INJECTIONS)
    (tmp_path / 'x_filtered.csv').write_bytes((tmp_path / 'x_formatted.csv').read_bytes())
    params.samplenorm, params.blnkgrp = 'pqn', ''

    factors = stats.samplenormalize(params)
    assert factors.loc[('C', 's5'), 'factor'] > 5
    out = capsys.readouterr().out
    assert 'more than 5x from the median sample in s4 (0.01x), s5 (0.01x)' in out

    # Again from the unscaled table, with C left out.
    (tmp_path / 'x_filtered.csv').write_bytes((tmp_path / 'x_filtered_raw.csv').read_bytes())
    params.samplenorm_skip = ['C']
    factors = stats.samplenormalize(params)
    out = capsys.readouterr().out
    assert factors.loc[('C', 's5'), 'factor'] == 1.0
    assert np.isnan(factors.loc[('C', 's5'), 'relative signal'])
    assert 'left unscaled: C' in out and 'warning' not in out
    # A and B are the only samples used: the median of 1 and 1.2.
    assert factors.loc[('A', 's1'), 'relative signal'] == pytest.approx(1 / 1.1)


def technical_cv(table):
    by_sample = table.T.groupby(level=[0, 1])
    return (by_sample.std() / by_sample.mean()).T


@pytest.mark.parametrize('method', ['sum', 'pqn'])
def test_pipeline_with_cv_cutoff_and_normalization(tmp_path, method):
    params = make_params(tmp_path)
    params.cvgate, params.cvgate_frac = True, 0.01
    params.samplenorm = method

    run_MSFaST(params)

    out = params.outputdir
    stem = params.filename.stem
    factors = pd.read_csv(out / f'{stem}_normfactors.csv', index_col=[0, 1])
    assert (factors['factor'] > 0).all()
    raw = pd.read_csv(out / f'{stem}_filtered_raw.csv', header=[0, 1, 2], index_col=[0, 1, 2])
    scaled = pd.read_csv(out / f'{stem}_filtered.csv', header=[0, 1, 2], index_col=[0, 1, 2])
    assert raw.shape == scaled.shape and raw.index.equals(scaled.index)
    # Technical CVs are the same before and after.
    a, b = technical_cv(raw).values, technical_cv(scaled).values
    assert np.allclose(a[~np.isnan(a)], b[~np.isnan(a)])
    info = (out / 'analysisinfo.txt').read_text()
    assert 'CV sample cutoff (beta): 0.01' in info
    assert 'Sample normalization (beta): ' in info
