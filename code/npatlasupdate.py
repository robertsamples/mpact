"""
MPACT
Copyright 2026, Robert M. Samples

Qt-free updater for the bundled Natural Products Atlas database
(``npatlas.tsv``). The Natural Products Atlas (https://www.npatlas.org)
publishes periodic full-database downloads; the copy MPACT ships with goes
stale over time as new compounds are deposited.

This module provides the staleness check and the download/validation logic,
with no Qt dependency, so it can be unit-tested headlessly (see
``tests/test_npatlasupdate.py``). The GUI side (``main.py``) only has to ask
the user a yes/no question and call :func:`download_atlas`.

Format decision (do not "upgrade" without reason): MPACT reads the atlas as a
tab-separated table via ``pd.read_csv('npatlas.tsv', sep='\\t', ...)`` and
``dbsearch.search_npatlas`` accesses specific columns
(``compound_m_plus_h``/``compound_m_plus_na``/``compound_smiles``/
``origin_type``/``genus``). The published ``NPAtlas_download.tsv`` already has
exactly those columns, so the TSV is a drop-in replacement and is what we
fetch. The ``NPAtlas_download.json`` is the same data in a nested JSON shape
that would need flattening before pandas/dbsearch could use it -- there is no
benefit to switching formats and a real cost (rewriting the read + column
access), so we deliberately stay on the TSV.
"""

import os
import shutil
import tempfile
import time
import urllib.request

DEFAULT_TSV_URL = 'https://www.npatlas.org/static/downloads/NPAtlas_download.tsv'
DEFAULT_JSON_URL = 'https://www.npatlas.org/static/downloads/NPAtlas_download.json'

# Columns the app actually consumes (main.py's atlas read + dbsearch). A
# download missing any of these is rejected rather than allowed to clobber a
# working atlas -- guards against the server returning an HTML error page or a
# truncated/renamed export.
REQUIRED_COLUMNS = frozenset({
    'compound_id',
    'compound_m_plus_h',
    'compound_m_plus_na',
    'compound_smiles',
    'origin_type',
    'genus',
})

DEFAULT_MAX_AGE_DAYS = 30


def atlas_age_days(path, now=None):
    """Age of the atlas file in days based on its mtime, or ``None`` if the
    file doesn't exist."""
    if not os.path.exists(path):
        return None
    now = time.time() if now is None else now
    return (now - os.path.getmtime(path)) / 86400.0


def is_update_due(path, max_age_days=DEFAULT_MAX_AGE_DAYS, now=None):
    """True if the atlas is missing or older than ``max_age_days``.

    This is the cheap check the app runs at startup before deciding whether to
    prompt the user -- it only stats the file, it never touches the network.

    Caveat: staleness is judged by file *mtime*, which is the requested
    "last modified over a month ago" behaviour but reflects when the file was
    last written locally, not the vintage of the data inside it. A fresh
    ``git clone`` stamps the checked-out file with the clone time, so a
    just-cloned-but-data-old atlas will read as "fresh" until 30 days pass.
    Embedding the NPAtlas release date in a sidecar and comparing that would
    be more accurate if this ever matters.
    """
    age = atlas_age_days(path, now=now)
    return age is None or age > max_age_days


def validate_tsv_header(first_line):
    """True if a TSV header line contains every :data:`REQUIRED_COLUMNS`."""
    columns = {c.strip() for c in first_line.rstrip('\n').split('\t')}
    return REQUIRED_COLUMNS.issubset(columns)


def download_atlas(dest, url=DEFAULT_TSV_URL, opener=None, validate=True,
                   timeout=60):
    """Download the NPAtlas TSV to ``dest`` atomically.

    The download streams to a temporary file in ``dest``'s directory, is
    validated (the header must contain :data:`REQUIRED_COLUMNS`), and only then
    ``os.replace``-d over ``dest`` -- so a network error, an HTML error page,
    or a partial transfer can never leave a corrupt or truncated atlas in
    place; the previous file is untouched on any failure.

    Args:
        dest: path the atlas should end up at (e.g. ``code/npatlas.tsv``).
        url: download URL (defaults to the published full TSV).
        opener: callable ``opener(url, timeout=...)`` returning a readable,
            context-managed response (defaults to ``urllib.request.urlopen``).
            Injectable so tests can supply canned content without a network.
        validate: if True, reject a download whose header is missing required
            columns (raises ``ValueError``, leaving ``dest`` unchanged).
        timeout: per-request timeout in seconds (urllib only).

    Returns:
        Number of bytes written to ``dest``.
    """
    dest = os.fspath(dest)
    dest_dir = os.path.dirname(os.path.abspath(dest)) or '.'
    opener = opener if opener is not None else urllib.request.urlopen

    fd, tmp_path = tempfile.mkstemp(prefix='.npatlas_', suffix='.tsv', dir=dest_dir)
    try:
        with os.fdopen(fd, 'wb') as tmp_file:
            try:
                response = opener(url, timeout=timeout)
            except TypeError:
                # Injected openers in tests may not accept a timeout kwarg.
                response = opener(url)
            with response as resp:
                shutil.copyfileobj(resp, tmp_file)
        bytes_written = os.path.getsize(tmp_path)

        if bytes_written == 0:
            raise ValueError('Downloaded atlas is empty')
        if validate:
            with open(tmp_path, 'r', encoding='utf-8', errors='replace') as check:
                first_line = check.readline()
            if not validate_tsv_header(first_line):
                raise ValueError(
                    'Downloaded file is not a valid NPAtlas TSV '
                    '(missing required columns); keeping existing atlas')

        os.replace(tmp_path, dest)
        return bytes_written
    except BaseException:
        # Clean up the temp file on ANY failure (including the validation
        # errors above); never disturb the existing dest.
        try:
            if os.path.exists(tmp_path):
                os.remove(tmp_path)
        except OSError:
            pass
        raise
