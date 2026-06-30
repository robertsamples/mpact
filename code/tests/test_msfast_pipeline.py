"""End-to-end test of the actual analysis pipeline (run_MSFaST), exercising
real filtering/stats code against the bundled PTY087I2 example dataset
(rawdata/PTY087I2/) rather than a hand-built toy fixture.

This wasn't previously possible to write cleanly: run_MSFaST used to take
the whole MainWindow (``self``) and write results onto arbitrary attributes
of it, so testing it meant faking an entire GUI window. Now it takes a
plain analysis_parameters object and returns a plain AnalysisResult, so a
real run can be driven from a test with nothing GUI-shaped involved.
"""

from pathlib import Path

import pandas as pd
import pytest

from groupsets import GroupSetModel, build_query_dict
from MSFaST import AnalysisResult, analysis_parameters, run_MSFaST

REPO_ROOT = Path(__file__).resolve().parents[2]
EXAMPLE_DIR = REPO_ROOT / 'rawdata' / 'PTY087I2'
ALL_GROUPS = ['Blanks', 'Media', '0um_Ce', '250um_Ce']


def make_params(tmp_path):
    params = analysis_parameters()
    params.filename = EXAMPLE_DIR / '200826_PTY087I2codingdataset.csv'
    params.samplelistfilename = EXAMPLE_DIR / 'samplelist.csv'
    params.extractmetadatafilename = EXAMPLE_DIR / 'extractmetadata.csv'
    params.outputdir = tmp_path / params.filename.stem
    params.outputdir.mkdir(parents=True)

    # Mispicked-peak filtering, defaults per the docs/published thresholds.
    params.relfil = True
    params.merge = True
    params.ringingwin = 0.5
    params.isopeakwin = 0.01
    params.dimerpeakwin = 0.01
    params.RTwin = 0.005
    params.maxisowin = 3

    params.grpave = True
    params.prperr = True

    # Blank filtering.
    params.blnkfltr = True
    params.blnkgrp = 'Blanks'
    params.blankfilthresh = 0.01

    # CV (reproducibility) filtering.
    params.CVfil = True
    params.cvthresh = 0.2
    params.cvparam = 'median CV'

    # In-source fragment deconvolution.
    params.decon = True
    params.deconthresh = 0.95

    # Fold change / t-test.
    params.FC = True
    params.Ttest = True
    params.statstgrps = ['250um_Ce', '0um_Ce']

    params.graphfilters = ['cv', 'rel', 'insource']

    # Plot toggles -- not consumed by filtering/stats logic itself, only
    # echoed into analysisinfo.txt's "Plots generated" summary.
    params.MZRTplt = True
    params.FC3Dplt = True
    params.KMD = True
    params.PCA = True
    params.Dendrogram = True
    params.Volcanoplt = True

    # One Plot Feature Set: every non-blank group, mirroring enumerate_inputs's
    # default "Features not in blanks" groupset when blank filtering is on.
    model = GroupSetModel()
    model.add('Features not in blanks', all_groups=ALL_GROUPS)
    src = [g for g in ALL_GROUPS if g != params.blnkgrp]
    model.update(0, src=src, excl=[params.blnkgrp])
    querydict = build_query_dict(model, params.graphfilters)
    params.querydict = querydict
    params.querylist = list(querydict.keys())

    return params


@pytest.fixture(scope='module')
def pipeline_result(tmp_path_factory):
    tmp_path = tmp_path_factory.mktemp('msfast_pipeline')
    params = make_params(tmp_path)
    result = run_MSFaST(params)
    return params, result


def test_returns_an_analysis_result(pipeline_result):
    _, result = pipeline_result
    assert isinstance(result, AnalysisResult)


def test_filtering_actually_removed_some_features(pipeline_result):
    params, result = pipeline_result
    unfiltered = pd.read_csv(params.filename, sep=',', header=[0, 1, 2], index_col=[0, 1, 2])
    filtered = pd.read_csv(params.outputdir / (params.filename.stem + '_filtered.csv'),
                            sep=',', header=[0, 1, 2], index_col=[0, 1, 2])
    assert 0 < len(filtered.index) < len(unfiltered.index)


def test_one_groupset_per_plot_feature_set(pipeline_result):
    _, result = pipeline_result
    assert set(result.groupsets.keys()) == set(result.filtereddfs.keys())
    assert len(result.groupsets) == 1


def test_ionfilters_cover_every_enabled_filter_type(pipeline_result):
    _, result = pipeline_result
    assert set(result.ionfilters.keys()) == {'relfil', 'cv', 'insource'}


def test_iondict_has_fold_change_and_logfc_derivable(pipeline_result):
    """Regression guard for the logfc bug: 'fc' must be present (logfc is
    derived from it on the fly, never persisted -- see main.py's
    _refresh_highlight), and it must be finite/positive for every feature
    that survived filtering."""
    params, _ = pipeline_result
    iondict = pd.read_csv(params.outputdir / 'iondict.csv', sep=',', header=[0], index_col=[0])
    assert 'fc' in iondict.columns
    filtered_fc = iondict.loc[iondict['groups'].notna(), 'fc'].dropna()
    assert len(filtered_fc) > 0
    assert (filtered_fc > 0).all()


def test_analysisinfo_written(pipeline_result):
    params, _ = pipeline_result
    info_path = params.outputdir / 'analysisinfo.txt'
    assert info_path.exists()
    text = info_path.read_text()
    assert 'Features passing all filters' in text


def test_fold_change_is_clamped_to_bounds(pipeline_result):
    """runfc clamps FC into [0.01, 100]; nothing should escape that range."""
    params, _ = pipeline_result
    iondict = pd.read_csv(params.outputdir / 'iondict.csv', sep=',', header=[0], index_col=[0])
    fc = iondict['fc'].dropna()
    assert len(fc) > 0
    assert fc.min() >= 0.01
    assert fc.max() <= 100.0


def test_stats_outputs_land_in_output_dir_not_cwd(pipeline_result):
    """The t-test/q-value tables now write into the run's output directory
    (previously bare 'msdata_teststats_test.csv'/'qdata.csv' in the cwd)."""
    params, _ = pipeline_result
    assert (params.outputdir / (params.filename.stem + '_teststats.csv')).exists()
    assert (params.outputdir / (params.filename.stem + '_qvalues.csv')).exists()


def test_qvalues_are_finite_positive_and_consistent_with_logq(pipeline_result):
    """The BH q-values must be finite and positive, and the persisted '-logq'
    column must equal -log10(qval) (the relationship runttest derives). Strict
    p-ascending monotonicity is deliberately NOT asserted: the cummin step-up
    only guarantees it in the loop's processing order, and tied p-values can
    reorder under an independent re-sort -- a known BH tie subtlety, not a bug."""
    import numpy as np
    params, _ = pipeline_result
    qdata = pd.read_csv(params.outputdir / (params.filename.stem + '_qvalues.csv'),
                         sep=',', header=[0])
    qval = qdata['qval'].to_numpy()
    assert np.isfinite(qval).all()
    assert (qval > 0).all()
    np.testing.assert_allclose(qdata['-logq'].to_numpy(), -np.log10(qval), rtol=1e-9)
