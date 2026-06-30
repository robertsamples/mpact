"""Regression test for the ``groupionlists`` defensive-init fix in
``run_MSFaST`` (MSFaST.py).

``groupionlists`` is only *populated* inside ``if analysis_params.grpave:``,
but is referenced unconditionally afterwards (the ``groupionlists['cv'/...]``
assignments and the iondict ``groups``-column loop). The GUI hardcodes
``grpave=True`` so this never fired in production, but a minimal run with
``grpave=False`` would have raised ``NameError`` before the fix. This drives
exactly that path against the bundled example dataset, with every optional
filter/stat turned off, and asserts the run completes and returns a result.
"""

from pathlib import Path

import pytest

from MSFaST import AnalysisResult, analysis_parameters, run_MSFaST

REPO_ROOT = Path(__file__).resolve().parents[2]
EXAMPLE_DIR = REPO_ROOT / 'rawdata' / 'PTY087I2'
ALL_GROUPS = ['Blanks', 'Media', '0um_Ce', '250um_Ce']


def _minimal_params(tmp_path):
    """Everything that gates a filtering/stats stage turned OFF, so the run
    exercises the no-grpave branch. Threshold/echo-only fields still have to
    be present because analysisinfo.txt prints them verbatim."""
    params = analysis_parameters()
    params.filename = EXAMPLE_DIR / '200826_PTY087I2codingdataset.csv'
    params.samplelistfilename = EXAMPLE_DIR / 'samplelist.csv'
    params.extractmetadatafilename = EXAMPLE_DIR / 'extractmetadata.csv'
    params.outputdir = tmp_path / params.filename.stem
    params.outputdir.mkdir(parents=True)

    # All optional stages OFF -- this is the configuration that used to crash.
    params.relfil = False
    params.merge = False
    params.grpave = False
    params.prperr = False
    params.blnkfltr = False
    params.CVfil = False
    params.decon = False
    params.FC = False
    params.Ttest = False

    # Thresholds / echo-only fields (printed into analysisinfo.txt).
    params.ringingwin = 0.5
    params.isopeakwin = 0.01
    params.dimerpeakwin = 0.01
    params.RTwin = 0.005
    params.maxisowin = 3
    params.blnkgrp = ''
    params.cvthresh = 0.2
    params.statstgrps = ['250um_Ce', '0um_Ce']
    params.graphfilters = []
    params.MZRTplt = params.FC3Dplt = params.KMD = False
    params.PCA = params.Dendrogram = params.Volcanoplt = False

    # No Plot Feature Sets configured -> empty querylist/querydict.
    params.querydict = {}
    params.querylist = []
    return params


def test_run_msfast_with_grpave_off_does_not_raise(tmp_path):
    params = _minimal_params(tmp_path)
    result = run_MSFaST(params)  # used to raise NameError on groupionlists
    assert isinstance(result, AnalysisResult)
    assert isinstance(result.groupionlists, dict)
    # The three filter keys are always added, each empty since every filter is off.
    assert result.groupionlists == {'cv': [], 'relfil': [], 'insource': []}
    # With no filters applied, the filtered table is written and analysisinfo exists.
    assert (params.outputdir / (params.filename.stem + '_filtered.csv')).exists()
    assert (params.outputdir / 'analysisinfo.txt').exists()
