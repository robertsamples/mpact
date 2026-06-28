"""Unit tests for the import/export translator framework (``translators.py``).

Uses small synthetic fixtures written to ``tmp_path`` for deterministic checks,
plus a few smoke checks against the real example files at the repo root (skipped
automatically when those files aren't present).
"""

from pathlib import Path

import pytest

import translators as t

REPO_ROOT = Path(__file__).resolve().parents[2]


# --------------------------------------------------------------------------- #
# fixtures
# --------------------------------------------------------------------------- #

CANONICAL_PEAKTABLE = (
    ",,,gA,gA\n"
    ",,,s1,s2\n"
    "Compound,m/z,Retention time (min),inj1,inj2\n"
    "1.0_200.1,200.1,1.0,10,20\n"
    "2.0_300.2,300.2,2.0,30,40\n"
    "3.0_400.3,400.3,3.0,50,60\n"
)


def _write(tmp_path, name, text):
    p = tmp_path / name
    p.write_text(text)
    return p


# --------------------------------------------------------------------------- #
# detection
# --------------------------------------------------------------------------- #

def test_detect_progenesis(tmp_path):
    p = _write(tmp_path, 'pk.csv', CANONICAL_PEAKTABLE)
    assert t.detect_peaktable_format(p) == t.PROGENESIS


def test_detect_mzmine(tmp_path):
    p = _write(tmp_path, 'mz.csv', "row ID,row m/z,row retention time,A.mzML Peak area\n1,200.1,1.0,10\n")
    assert t.detect_peaktable_format(p) == t.MZMINE


def test_detect_metaboscape(tmp_path):
    p = _write(tmp_path, 'ms.csv', "Bucket label,foo,bar\nx,1,2\n")
    assert t.detect_peaktable_format(p) == t.METABOSCAPE


def test_detect_msdial_txt(tmp_path):
    text = "\t\t\n\t\t\n\t\t\n\t\t\nAlignment ID\tAverage Rt(min)\tAverage Mz\n0\t1.0\t200.1\n"
    p = _write(tmp_path, 'd.txt', text)
    assert t.detect_peaktable_format(p) == t.MSDIAL


def test_detect_unknown(tmp_path):
    p = _write(tmp_path, 'u.csv', "a,b,c\n1,2,3\n")
    assert t.detect_peaktable_format(p) == t.UNKNOWN


# --------------------------------------------------------------------------- #
# fragmentation parsing
# --------------------------------------------------------------------------- #

def test_parse_mgf_mzmine_style_has_mz_no_rt(tmp_path):
    mgf = "BEGIN IONS\nSCANS=1\nPEPMASS=335.12648\nCHARGE=1\n139.04 8898\nEND IONS\n"
    entries = t.parse_mgf(_write(tmp_path, 'a.mgf', mgf))
    assert len(entries) == 1
    assert entries[0].mz == pytest.approx(335.12648)
    assert entries[0].rt is None


def test_parse_mgf_msdial_style_has_rt(tmp_path):
    mgf = ("BEGIN IONS\nTITLE=Methionine||PEAKID=18\nPEPMASS=150.05829\n"
           "RTINMINUTES=3.826\nNum Peaks: 1\n56.05 334\nEND IONS\n")
    entries = t.parse_mgf(_write(tmp_path, 'b.mgf', mgf))
    assert entries[0].mz == pytest.approx(150.05829)
    assert entries[0].rt == pytest.approx(3.826)
    assert 'Methionine' in entries[0].comment


def test_parse_msp_progenesis_precursor_and_comment(tmp_path):
    msp = ("Name: Unknown (0.80_627.2171n)\nCharge:\nPrecursor: 627.2171\n"
           "Comment: 0.80_627.2171n\nNum peaks: 1\n418.1451 257254\n")
    entries = t.parse_msp(_write(tmp_path, 'a.msp', msp))
    assert len(entries) == 1
    assert entries[0].mz == pytest.approx(627.2171)
    assert entries[0].comment == '0.80_627.2171n'


def test_parse_msp_msdial_precursormz_and_rt(tmp_path):
    msp = ("NAME: Unknown|ID=0|MZ=150.0267|RT=9.09\nPRECURSORMZ: 150.0267\n"
           "RETENTIONTIME: 9.0898957\nNum Peaks: 0\n")
    entries = t.parse_msp(_write(tmp_path, 'b.msp', msp))
    assert entries[0].mz == pytest.approx(150.0267)
    assert entries[0].rt == pytest.approx(9.0898957)


# --------------------------------------------------------------------------- #
# GNPS2 re-indexing
# --------------------------------------------------------------------------- #

def test_reindex_msp_matches_by_compound_id_and_renumbers(tmp_path):
    pk = _write(tmp_path, 'pk.csv', CANONICAL_PEAKTABLE)
    # Entries out of order; precursor m/z deliberately NOT matching the table to
    # prove the compound-id (Comment) path is used.
    msp = (
        "Name: Unknown (3.0_400.3)\nPrecursor: 999.9\nComment: 3.0_400.3\nNum peaks: 1\n100.0 5\n\n"
        "Name: Unknown (1.0_200.1)\nPrecursor: 888.8\nComment: 1.0_200.1\nNum peaks: 1\n50.0 9\n"
    )
    src = _write(tmp_path, 'frags.msp', msp)
    out = tmp_path / 'out.msp'
    n = t.reindex_fragments(pk, src, out)
    assert n == 2
    entries = t.parse_msp(out)
    # Sorted by assigned row id: 1.0_200.1 (row 1) then 3.0_400.3 (row 3).
    assert [e.comment for e in entries] == ['1.0_200.1', '3.0_400.3']
    text = out.read_text()
    assert 'SCANNUMBER: 1' in text and 'SCANNUMBER: 3' in text


def test_reindex_mgf_matches_by_mz_and_rewrites_scans(tmp_path):
    pk = _write(tmp_path, 'pk.csv', CANONICAL_PEAKTABLE)
    mgf = (
        "BEGIN IONS\nSCANS=50\nPEPMASS=400.3\nCHARGE=1\n10 100\nEND IONS\n"
        "BEGIN IONS\nSCANS=51\nPEPMASS=200.1\nCHARGE=1\n11 110\nEND IONS\n"
    )
    src = _write(tmp_path, 'frags.mgf', mgf)
    out = tmp_path / 'out.mgf'
    n = t.reindex_fragments(pk, src, out, mz_tol=0.01)
    assert n == 2
    entries = t.parse_mgf(out)
    # First written entry is row 1 (m/z 200.1), second is row 3 (m/z 400.3).
    assert [e.mz for e in entries] == [pytest.approx(200.1), pytest.approx(400.3)]
    text = out.read_text()
    assert 'SCANS=1' in text and 'SCANS=3' in text
    assert 'SCANS=50' not in text and 'SCANS=51' not in text


# --------------------------------------------------------------------------- #
# filtered source-format export
# --------------------------------------------------------------------------- #

def test_filter_source_peaktable_keeps_only_surviving_rows(tmp_path):
    source = _write(tmp_path, 'src.csv', CANONICAL_PEAKTABLE)
    # Internal filtered table dropping the middle feature (300.2).
    filtered = _write(tmp_path, 'flt.csv',
                      ",,,gA,gA\n,,,s1,s2\n"
                      "Compound,m/z,Retention time (min),inj1,inj2\n"
                      "1.0_200.1,200.1,1.0,10,20\n"
                      "3.0_400.3,400.3,3.0,50,60\n")
    out = tmp_path / 'src_filtered.csv'
    kept = t.filter_source_peaktable(source, filtered, out)
    assert kept == 2
    lines = [ln for ln in out.read_text().splitlines() if ln.strip()]
    # 3 header rows + 2 data rows
    assert len(lines) == 5
    assert '200.1' in out.read_text() and '400.3' in out.read_text()
    assert '300.2' not in out.read_text()


# --------------------------------------------------------------------------- #
# matching helper
# --------------------------------------------------------------------------- #

def test_match_index_respects_tolerance():
    features = [('a', 200.1, 1.0), ('b', 300.2, 2.0)]
    assert t._match_index(200.105, None, features, 0.01, 0.05) == 0
    assert t._match_index(250.0, None, features, 0.01, 0.05) is None


# --------------------------------------------------------------------------- #
# real example files (skipped when absent)
# --------------------------------------------------------------------------- #

@pytest.mark.parametrize('name,expected', [
    ('nativehead.csv', t.PROGENESIS),
    ('progenesis.csv', t.PROGENESIS),
    ('mzminehead.csv', t.MZMINE),
    ('msdial_unformatted.txt', t.MSDIAL),
])
def test_detect_real_example_files(name, expected):
    path = REPO_ROOT / name
    if not path.exists():
        pytest.skip(name + ' not present')
    assert t.detect_peaktable_format(path) == expected


@pytest.mark.parametrize('name', ['progenesis.msp', 'msdial.msp', 'mzmine.mgf', 'msdial.mgf'])
def test_parse_real_fragment_files(name):
    path = REPO_ROOT / name
    if not path.exists():
        pytest.skip(name + ' not present')
    entries = t.parse_fragments(path)
    assert len(entries) > 0
    assert all(e.mz is not None for e in entries)
