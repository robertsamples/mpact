"""Unit tests for the NPAtlas updater (``npatlasupdate.py``).

The network is never touched: ``download_atlas`` takes an injectable
``opener``, so tests feed canned bytes through a fake response object and
assert the staleness logic, header validation, and atomic-replace behaviour.
"""

import io
import os
import time

import pytest

import npatlasupdate as nu


# --------------------------------------------------------------------------- #
# fixtures / helpers
# --------------------------------------------------------------------------- #

VALID_HEADER = (
    'npaid\tcompound_id\tcompound_name\tcompound_m_plus_h\tcompound_m_plus_na\t'
    'compound_smiles\torigin_type\tgenus\n'
)
VALID_TSV = VALID_HEADER + '1\t0.80_418n\tFoo\t419.1\t441.1\tCCO\tBacterium\tStreptomyces\n'


class _FakeResponse(io.BytesIO):
    """A BytesIO that also works as a context manager (like urlopen's result)."""
    def __enter__(self):
        return self

    def __exit__(self, *exc):
        self.close()
        return False


def _opener_returning(content_bytes, record=None):
    def opener(url, timeout=None):
        if record is not None:
            record.append(url)
        return _FakeResponse(content_bytes)
    return opener


# --------------------------------------------------------------------------- #
# staleness
# --------------------------------------------------------------------------- #

def test_age_is_none_when_missing(tmp_path):
    assert nu.atlas_age_days(tmp_path / 'nope.tsv') is None


def test_update_due_when_missing(tmp_path):
    assert nu.is_update_due(tmp_path / 'nope.tsv') is True


def test_update_not_due_for_fresh_file(tmp_path):
    p = tmp_path / 'npatlas.tsv'
    p.write_text(VALID_TSV)
    # Just created -> age ~0 days -> not due.
    assert nu.is_update_due(p, max_age_days=30) is False


def test_update_due_for_old_file(tmp_path):
    p = tmp_path / 'npatlas.tsv'
    p.write_text(VALID_TSV)
    # Backdate mtime to 45 days ago.
    old = time.time() - 45 * 86400
    os.utime(p, (old, old))
    assert nu.is_update_due(p, max_age_days=30) is True
    assert nu.atlas_age_days(p) == pytest.approx(45, abs=0.1)


# --------------------------------------------------------------------------- #
# header validation
# --------------------------------------------------------------------------- #

def test_validate_header_accepts_full_header():
    assert nu.validate_tsv_header(VALID_HEADER) is True


def test_validate_header_rejects_missing_columns():
    assert nu.validate_tsv_header('npaid\tcompound_id\tgenus\n') is False


def test_validate_header_rejects_html_error_page():
    assert nu.validate_tsv_header('<!DOCTYPE html>') is False


# --------------------------------------------------------------------------- #
# download (atomic + validated)
# --------------------------------------------------------------------------- #

def test_download_writes_validated_tsv(tmp_path):
    dest = tmp_path / 'npatlas.tsv'
    seen = []
    n = nu.download_atlas(dest, url='http://example/atlas.tsv',
                          opener=_opener_returning(VALID_TSV.encode(), record=seen))
    assert dest.exists()
    assert n == len(VALID_TSV.encode())
    assert nu.validate_tsv_header(dest.read_text().splitlines(keepends=True)[0])
    assert seen == ['http://example/atlas.tsv']


def test_download_overwrites_existing_atlas_atomically(tmp_path):
    dest = tmp_path / 'npatlas.tsv'
    dest.write_text('OLD CONTENT')
    nu.download_atlas(dest, opener=_opener_returning(VALID_TSV.encode()))
    assert 'OLD CONTENT' not in dest.read_text()
    assert dest.read_text() == VALID_TSV


def test_invalid_download_is_rejected_and_existing_atlas_preserved(tmp_path):
    dest = tmp_path / 'npatlas.tsv'
    dest.write_text(VALID_TSV)  # a good existing atlas
    bad = b'<html>503 Service Unavailable</html>'
    with pytest.raises(ValueError):
        nu.download_atlas(dest, opener=_opener_returning(bad))
    # The good atlas must be untouched, and no temp files left behind.
    assert dest.read_text() == VALID_TSV
    leftovers = [f for f in os.listdir(tmp_path) if f.startswith('.npatlas_')]
    assert leftovers == []


def test_empty_download_is_rejected(tmp_path):
    dest = tmp_path / 'npatlas.tsv'
    dest.write_text(VALID_TSV)
    with pytest.raises(ValueError):
        nu.download_atlas(dest, opener=_opener_returning(b''))
    assert dest.read_text() == VALID_TSV


def test_network_error_leaves_existing_atlas(tmp_path):
    dest = tmp_path / 'npatlas.tsv'
    dest.write_text(VALID_TSV)

    def failing_opener(url, timeout=None):
        raise OSError('connection refused')

    with pytest.raises(OSError):
        nu.download_atlas(dest, opener=failing_opener)
    assert dest.read_text() == VALID_TSV
    assert [f for f in os.listdir(tmp_path) if f.startswith('.npatlas_')] == []
