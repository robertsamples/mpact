"""Tests for the opt-in beta filters in ``filter.py``: the mass defect filter,
adduct / neutral-loss grouping, and the beta peak-correction mode.

The small cases build the two intermediate files those functions read
(``<stem>_formatted.csv`` and ``iondict.csv``) by hand, so each rule can be
checked on a table where the right answer is known. The last test switches all
three on for a real run of the bundled PTY087I2 example.
"""

from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pandas as pd
import pytest

import filter
from MSFaST import run_MSFaST
from test_msfast_pipeline import make_params

RISING = [100.0, 200.0, 300.0, 400.0, 500.0, 600.0]
SHUFFLED = [600.0, 100.0, 500.0, 200.0, 400.0, 300.0]   # uncorrelated with RISING
NA_DELTA = 22.989218 - 1.007276     # [M+Na]+ - [M+H]+
K_DELTA = 38.963158 - 1.007276      # [M+K]+ - [M+H]+
C13 = 1.0033548


def write_tables(tmp_path, features):
    """features: [(id, m/z, RT, [six intensities]), ...] -> params for the filters."""
    groups = ['A', 'A', 'A', 'B', 'B', 'B']
    columns = pd.MultiIndex.from_tuples(
        [(f'Unnamed: {i}_level_0', f'Unnamed: {i}_level_1', name)
         for i, name in enumerate(['Compound', 'm/z', 'Retention time (min)'])]
        + [(group, f'sample{i}', f'injection{i}') for i, group in enumerate(groups)])
    rows = [[fid, mz, rt, *intensities] for fid, mz, rt, intensities in features]
    pd.DataFrame(rows, columns=columns).to_csv(tmp_path / 'x_formatted.csv', index=False)
    pd.DataFrame({'Compound': [f[0] for f in features],
                  'm/z': [f[1] for f in features],
                  'Retention time (min)': [f[2] for f in features]}).to_csv(tmp_path / 'iondict.csv', index=False)
    return SimpleNamespace(outputdir=tmp_path, filename=Path('x.csv'), RTwin=0.02,
                           maxisowin=3, ringingwin=0.5, isopeakwin=0.01, dimerpeakwin=0.01)


def scaled(values, factor):
    return [v * factor for v in values]


# ---- mass defect ---------------------------------------------------------

def test_md_plausible_accepts_ordinary_ions_and_rejects_half_integer_ones():
    ok = filter.md_plausible([419.1521, 210.0803, 165.5441, 281.5609],
                             0.00112, 0.03, -0.00025, -0.03)
    assert ok.tolist() == [True, True, False, False]


def test_md_plausible_cannot_reject_anything_at_high_mass():
    mz = np.linspace(700, 2000, 500)
    assert filter.md_plausible(mz, 0.00112, 0.03, -0.00025, -0.03).all()


def test_mdfilter_flags_and_drops_implausible_ions(tmp_path):
    params = write_tables(tmp_path, [('a', 419.1521, 1.0, RISING),
                                     ('b', 165.5441, 2.0, RISING),
                                     ('c', 210.0803, 3.0, RISING)])
    filters = filter.mdfilter(params, {})

    assert filters['massdefect'].ions == ['b']
    iondict = pd.read_csv(tmp_path / 'iondict.csv', index_col=0)
    assert iondict['pass_mdfil'].to_dict() == {'a': True, 'b': False, 'c': True}
    # dropped from the working table, which must still have its three header rows
    table = pd.read_csv(tmp_path / 'x_formatted.csv', header=[0, 1, 2], index_col=0)
    assert table.index.tolist() == ['a', 'c']
    assert table.columns[0][2] == 'm/z'


# ---- adduct grouping -----------------------------------------------------

@pytest.mark.parametrize('ids', [('base', 'sodium', 'elsewhere', 'unrelated'), (1, 2, 3, 4)],
                         ids=['text ids', 'integer ids'])
def test_adductgroup_needs_mass_rt_and_correlation(tmp_path, ids):
    base, sodium, elsewhere, unrelated = ids
    params = write_tables(tmp_path, [
        (base, 300.1000, 5.00, RISING),
        (sodium, 300.1000 + NA_DELTA, 5.00, scaled(RISING, 0.3)),      # real [M+Na]+
        (elsewhere, 300.1000 + NA_DELTA, 7.00, scaled(RISING, 0.3)),   # right mass, wrong RT
        (unrelated, 300.1000 + K_DELTA, 5.00, SHUFFLED),               # right mass and RT, uncorrelated
    ])
    params.adduct_mode = 'pos'

    removed = filter.adductgroup(params)

    assert removed == {sodium}
    iondict = pd.read_csv(tmp_path / 'iondict.csv', index_col=0)
    assert iondict['pass_adduct'].tolist() == [True, False, True, True]
    assert iondict.loc[base, 'adduct_form'] == 'M+H'
    assert iondict.loc[sodium, 'adduct_form'] == 'M+Na'
    attribution = pd.read_csv(tmp_path / 'adduct_attribution.csv').dropna(how='all')
    assert attribution['Compound'].tolist() == [base, sodium]
    assert attribution['is_main'].fillna('').tolist() == ['main', '']


def test_default_sodium_and_potassium_deltas_are_ion_masses(tmp_path):
    # At 1 ppm the neutral-atom masses the defaults once used (0.00055 Da high)
    # miss both of these.
    params = write_tables(tmp_path, [('base', 300.1000, 5.0, RISING),
                                     ('sodium', 300.1000 + NA_DELTA, 5.0, scaled(RISING, 0.3)),
                                     ('potassium', 300.1000 + K_DELTA, 5.0, scaled(RISING, 0.2))])
    params.adduct_mode = 'pos'
    params.adduct_ppm = 1.0
    assert filter.adductgroup(params) == {'sodium', 'potassium'}


def test_adductgroup_keeps_the_most_abundant_member(tmp_path):
    params = write_tables(tmp_path, [('weak_base', 300.1000, 5.0, RISING),
                                     ('strong_sodium', 300.1000 + NA_DELTA, 5.0, scaled(RISING, 4))])
    params.adduct_mode = 'pos'
    assert filter.adductgroup(params) == {'weak_base'}


def test_adductgroup_negative_mode_uses_the_given_list_and_base_label(tmp_path):
    params = write_tables(tmp_path, [('base', 300.1000, 5.0, RISING),
                                     ('chloride', 300.1000 + 35.9767, 5.0, scaled(RISING, 0.2))])
    params.adduct_mode = 'neg'
    assert filter.adductgroup(params, [('M+Cl', 35.9767)]) == {'chloride'}
    iondict = pd.read_csv(tmp_path / 'iondict.csv', index_col=0)
    assert iondict['adduct_form'].to_dict() == {'base': 'M-H', 'chloride': 'M+Cl'}


# ---- peak correction: legacy stays the default ----------------------------

def peak_correction_tables(tmp_path):
    return write_tables(tmp_path, [
        ('parent', 300.1000, 5.0, RISING),
        ('isotope', 300.1000 + C13, 5.0, scaled(RISING, 0.2)),
        # 0.2 Da above the parent and five times stronger: inside the legacy
        # ringing band, but not something an intensity check calls an artefact
        ('stronger_neighbour', 300.3000, 5.0, scaled(RISING, 5)),
    ])


def test_peak_correction_defaults_to_the_legacy_filter(tmp_path):
    params = peak_correction_tables(tmp_path)
    filters = filter.relationalfilter(params, {})
    assert sorted(filters['relfil'].ions) == ['isotope', 'stronger_neighbour']


def test_beta_peak_correction_requires_the_removed_peak_to_be_weaker(tmp_path):
    params = peak_correction_tables(tmp_path)
    params.relfil_mode = 'beta'
    filters = filter.relationalfilter(params, {})
    assert filters['relfil'].ions == ['isotope']
    assert filters['relfil'].merge['parent'].sources == ['isotope']
    iondict = pd.read_csv(tmp_path / 'iondict.csv', index_col=0)
    assert iondict['pass_relfil'].to_dict() == {'parent': True, 'isotope': False, 'stronger_neighbour': True}


# ---- all three on, real data ----------------------------------------------

def test_pipeline_runs_with_every_beta_feature_on(tmp_path):
    params = make_params(tmp_path)
    params.mdfil = True
    params.adductgrp = True
    params.adduct_mode = 'pos'
    params.relfil_mode = 'beta'

    result = run_MSFaST(params)

    assert {'massdefect', 'relfil', 'cv', 'insource', 'adduct'} <= set(result.ionfilters)
    assert len(result.ionfilters['massdefect'].ions) > 0
    iondict = pd.read_csv(params.outputdir / 'iondict.csv', index_col=0)
    assert {'pass_mdfil', 'pass_relfil', 'pass_adduct', 'adduct_form'} <= set(iondict.columns)
    filtered = pd.read_csv(params.outputdir / (params.filename.stem + '_filtered.csv'),
                           header=[0, 1, 2], index_col=[0, 1, 2])
    survivors = set(filtered.index.get_level_values(0).astype(str))
    assert survivors and not survivors & set(map(str, result.ionfilters['massdefect'].ions))


def test_mass_defect_filter_runs_without_peak_correction(tmp_path):
    params = make_params(tmp_path)
    params.relfil = False
    params.merge = False
    params.mdfil = True

    result = run_MSFaST(params)

    assert 'relfil' not in result.ionfilters
    assert len(result.ionfilters['massdefect'].ions) > 0
    assert result.groupionlists['massdefect'] == result.ionfilters['massdefect'].ions


def test_beta_removed_features_are_in_no_plot_feature_set(tmp_path):
    params = make_params(tmp_path)
    params.mdfil = True
    params.adductgrp = True
    params.adduct_mode = 'pos'

    result = run_MSFaST(params)

    removed = set(map(str, result.ionfilters['massdefect'].ions)) | set(map(str, result.ionfilters['adduct'].ions))
    plotted = {str(ion) for groupset in result.groupsets.values() for ion in groupset.ionlist}
    assert removed and plotted and not removed & plotted
    info = (params.outputdir / 'analysisinfo.txt').read_text()
    assert 'Features failing mass defect filtering: ' in info
    assert 'Features removed by adduct grouping: ' in info
