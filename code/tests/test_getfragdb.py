"""Unit tests for the MSP fragmentation-database importer (``getfragdb.py``).

Covers the two parsers (``importfrag_v1`` Progenesis-style, ``importfrag_v2``
MS-DIAL-style) and the format auto-detection wrapper (``importfrag``), which
previously had no coverage at all. Synthetic MSP fixtures are written to
``tmp_path``; a couple of smoke checks run against the real example files at
the repo root (skipped when absent).
"""

from pathlib import Path

import pytest

import getfragdb

REPO_ROOT = Path(__file__).resolve().parents[2]


# Progenesis-style: "Name: Unknown (<id>)" with parenthetical id, no pipes.
PROGENESIS_MSP = (
    "Name: Unknown (0.80_627.2171n)\n"
    "PrecursorMZ: 627.2171\n"
    "Num Peaks: 2\n"
    "418.1451 257254\n"
    "200.1000 5000\n"
    "\n"
    "Name: Unknown (1.20_300.1000n)\n"
    "PrecursorMZ: 300.1\n"
    "Num Peaks: 1\n"
    "150.0500 100\n"
)

# MS-DIAL-style: "NAME: ...|ID=|MZ=|RT=" with pipes, PRECURSORMZ + RETENTIONTIME.
MSDIAL_MSP = (
    "NAME: Unknown|ID=0|MZ=150.0267|RT=9.09\n"
    "PRECURSORMZ: 150.0267\n"
    "RETENTIONTIME: 9.0898957\n"
    "Num Peaks: 2\n"
    "56.0500 334\n"
    "70.0600 120\n"
    "\n"
    "NAME: Unknown|ID=1|MZ=200.1000|RT=3.50\n"
    "PRECURSORMZ: 200.1\n"
    "RETENTIONTIME: 3.5\n"
    "Num Peaks: 1\n"
    "99.0000 50\n"
)


def _write(tmp_path, name, text):
    p = tmp_path / name
    p.write_text(text)
    return p


def test_importfrag_v1_parses_progenesis_ids_and_peaks(tmp_path):
    db = getfragdb.importfrag_v1(_write(tmp_path, 'p.msp', PROGENESIS_MSP))
    assert set(db.ions.keys()) == {'0.80_627.2171n', '1.20_300.1000n'}
    # First entry has 2 peaks parsed into an (n, 2) array.
    first = db.ions['0.80_627.2171n']
    assert first.pattern.shape == (2, 2)
    assert first.pattern[0][0] == pytest.approx(418.1451)
    assert db.ions['1.20_300.1000n'].pattern.shape == (1, 2)


def test_importfrag_v2_parses_msdial_rt_mz_keyed_ids(tmp_path):
    db = getfragdb.importfrag_v2(_write(tmp_path, 'd.msp', MSDIAL_MSP))
    # name = f"{round(rt,3)}_{precursormz}" using the raw PRECURSORMZ string.
    assert '9.09_150.0267' in db.ions
    assert '3.5_200.1' in db.ions
    assert db.ions['9.09_150.0267'].pattern.shape == (2, 2)
    assert db.ions['9.09_150.0267'].fragparams['PRECURSORMZ'] == '150.0267'


def test_importfrag_autodetects_progenesis(tmp_path):
    db = getfragdb.importfrag(_write(tmp_path, 'p.msp', PROGENESIS_MSP))
    # Progenesis ids are the parenthetical comment, not RT_mz keys.
    assert '0.80_627.2171n' in db.ions


def test_importfrag_autodetects_msdial(tmp_path):
    db = getfragdb.importfrag(_write(tmp_path, 'd.msp', MSDIAL_MSP))
    assert '9.09_150.0267' in db.ions


def test_importfrag_v2_skips_entries_without_rt_or_precursor(tmp_path):
    # An entry missing RETENTIONTIME/PRECURSORMZ must be dropped, not crash.
    msp = (
        "NAME: Unknown|ID=0|MZ=150\n"
        "Num Peaks: 1\n"
        "56.05 334\n"
        "\n"
        "NAME: Unknown|ID=1|MZ=200.1|RT=3.5\n"
        "PRECURSORMZ: 200.1\n"
        "RETENTIONTIME: 3.5\n"
        "Num Peaks: 1\n"
        "99.0 50\n"
    )
    db = getfragdb.importfrag_v2(_write(tmp_path, 'd.msp', msp))
    assert list(db.ions.keys()) == ['3.5_200.1']


@pytest.mark.parametrize('name', ['progenesis.msp', 'msdial.msp'])
def test_importfrag_on_real_example_files(name):
    path = REPO_ROOT / name
    if not path.exists():
        pytest.skip(name + ' not present')
    db = getfragdb.importfrag(path)
    assert len(db.ions) > 0
    # Every parsed ion's peak array is either empty (an entry with 0 peaks --
    # which does occur in the real MS-DIAL export) or a proper (n, 2) array.
    for entry in db.ions.values():
        assert entry.pattern.size == 0 or (
            entry.pattern.ndim == 2 and entry.pattern.shape[1] == 2)
