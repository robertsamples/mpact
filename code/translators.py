"""
MPACT
Copyright 2022, Robert M. Samples, Sara P. Puckett, and Marcy J. Balunas

Flexible import/export translator framework.

MPACT operates on one *canonical* peak-table layout (see ``mpact-file-formats``
notes): a CSV with three header rows where row 3 (index 2) is
``Compound, m/z, Retention time (min), <injection names...>`` and each data row
is ``<id>, <m/z>, <rt>, <abundances...>``. Different tools (Progenesis -- the
native format -- MZmine, MS-DIAL, Metaboscape) export different layouts; this
module detects the source format and converts it to the canonical layout.

It also provides the inverse direction that downstream tools need:
- writing the *filtered* peak table back out in the original source layout
  (a row subset of the source file), and
- re-indexing an MSP/MGF fragmentation file so its entry IDs line up with the
  peak-table row order/IDs, which GNPS2 requires.

These functions are deliberately free of any Qt/GUI dependency so they can be
unit-tested headlessly.
"""

import re
from pathlib import Path

import pandas as pd


# --------------------------------------------------------------------------- #
# Peak-table format detection
# --------------------------------------------------------------------------- #

PROGENESIS = 'progenesis'
MZMINE = 'mzmine'
MSDIAL = 'msdial'
METABOSCAPE = 'metaboscape'
UNKNOWN = 'unknown'


def _first_cell(value):
    """Normalise a cell to a stripped string ('' for NaN/None)."""
    if value is None:
        return ''
    try:
        if pd.isna(value):
            return ''
    except (TypeError, ValueError):
        pass
    return str(value).strip()


def detect_peaktable_format(path):
    """Return the source format of a peak-table file by inspecting its content.

    Detection is content-based (not extension- or position-locked) so it keeps
    working when tools shuffle metadata columns between versions.

    Returns one of ``PROGENESIS``, ``MZMINE``, ``MSDIAL``, ``METABOSCAPE`` or
    ``UNKNOWN``.
    """
    path = Path(path)

    # MS-DIAL alignment exports are tab-separated .txt files.
    if path.suffix.lower() == '.txt':
        head = pd.read_csv(path, sep='\t', header=None, dtype=str, nrows=8, engine='python')
        if _scan_for(head, 'Alignment ID') or _scan_for(head, 'Average Mz'):
            return MSDIAL
        return UNKNOWN

    head = pd.read_csv(path, sep=',', header=None, dtype=str, nrows=8)
    top_left = _first_cell(head.iloc[0, 0]) if head.shape[0] else ''

    if top_left == 'row ID':
        return MZMINE
    if top_left == 'Bucket label':
        return METABOSCAPE
    if _scan_for(head, 'Alignment ID'):
        return MSDIAL
    # Progenesis / canonical: 'Compound' + 'Retention time (min)' in the header block.
    if _scan_for(head, 'Compound') and _scan_for(head, 'Retention time (min)'):
        return PROGENESIS
    return UNKNOWN


def _scan_for(df, target):
    """True if any cell within the first few header rows equals ``target``."""
    target = target.strip()
    for _, row in df.head(6).iterrows():
        for cell in row:
            if _first_cell(cell) == target:
                return True
    return False


# --------------------------------------------------------------------------- #
# MSP / MGF fragmentation-file parsing
# --------------------------------------------------------------------------- #

class FragmentEntry:
    """One parsed fragmentation entry.

    Attributes:
        mz: precursor m/z (float) or None.
        rt: retention time in minutes (float) or None when the format omits it.
        comment: the entry's compound-id/comment string (e.g. ``0.80_627.2171n``)
            when present -- this matches the peak table ``Compound`` id exactly
            and is the most reliable way to line entries up with features.
        lines: the raw text lines of the entry (kept verbatim for rewriting).
    """
    __slots__ = ('mz', 'rt', 'comment', 'lines')

    def __init__(self, mz=None, rt=None, comment=None, lines=None):
        self.mz = mz
        self.rt = rt
        self.comment = comment
        self.lines = lines if lines is not None else []


def _to_float(text):
    try:
        return float(str(text).strip())
    except (TypeError, ValueError):
        return None


def parse_mgf(path):
    """Parse an MGF file into a list of :class:`FragmentEntry`.

    Handles both MZmine-style (``SCANS``/``PEPMASS``, no RT) and MS-DIAL-style
    (``RTINMINUTES``/``RTINSECONDS``) entries.
    """
    entries = []
    current = None
    with open(path, 'r') as handle:
        for line in handle:
            stripped = line.strip()
            if stripped == 'BEGIN IONS':
                current = FragmentEntry(lines=[line])
            elif stripped == 'END IONS':
                if current is not None:
                    current.lines.append(line)
                    entries.append(current)
                    current = None
            elif current is not None:
                current.lines.append(line)
                upper = stripped.upper()
                if upper.startswith('PEPMASS='):
                    current.mz = _to_float(stripped.split('=', 1)[1].split()[0])
                elif upper.startswith('RTINMINUTES='):
                    current.rt = _to_float(stripped.split('=', 1)[1])
                elif upper.startswith('RTINSECONDS='):
                    seconds = _to_float(stripped.split('=', 1)[1])
                    current.rt = seconds / 60.0 if seconds is not None else None
                elif upper.startswith('TITLE='):
                    current.comment = stripped.split('=', 1)[1].strip()
    return entries


def parse_msp(path):
    """Parse an MSP file into a list of :class:`FragmentEntry`.

    Handles both Progenesis-style (``PrecursorMZ``, RT encoded in the name) and
    MS-DIAL-style (``PRECURSORMZ``/``RETENTIONTIME``) entries.
    """
    entries = []
    current = None
    with open(path, 'r') as handle:
        for line in handle:
            if line.split(':', 1)[0].strip().upper() == 'NAME':
                if current is not None:
                    entries.append(current)
                # Fall back to the parenthetical in the Name line as the id.
                paren = re.search(r'\(([^)]*)\)', line)
                current = FragmentEntry(lines=[line],
                                        comment=paren.group(1) if paren else None)
                continue
            if current is None:
                continue
            current.lines.append(line)
            key = line.split(':', 1)[0].strip().upper()
            if key in ('PRECURSORMZ', 'PRECURSOR'):
                current.mz = _to_float(line.split(':', 1)[1])
            elif key == 'RETENTIONTIME':
                current.rt = _to_float(line.split(':', 1)[1])
            elif key == 'COMMENT':
                value = line.split(':', 1)[1].strip()
                if value:
                    current.comment = value
    if current is not None:
        entries.append(current)
    return entries


def parse_fragments(path):
    """Parse a fragmentation file (``.mgf`` or ``.msp``) into entries."""
    suffix = Path(path).suffix.lower()
    if suffix == '.mgf':
        return parse_mgf(path)
    if suffix == '.msp':
        return parse_msp(path)
    raise ValueError('Unsupported fragmentation file type: ' + suffix)


# --------------------------------------------------------------------------- #
# Feature matching helpers
# --------------------------------------------------------------------------- #

def _read_feature_table(filtered_path):
    """Read (id, m/z, rt) for every feature in a canonical/internal peak table.

    Returns a list of ``(compound_id, mz, rt)`` tuples in file order.
    """
    df = pd.read_csv(filtered_path, sep=',', header=[2], index_col=None)
    features = []
    for _, row in df.iterrows():
        mz = _to_float(row.get('m/z'))
        rt = _to_float(row.get('Retention time (min)'))
        compound = row.iloc[0]
        if mz is not None:
            features.append((compound, mz, rt))
    return features


def _match_index(mz, rt, features, mz_tol, rt_tol):
    """Index of the closest feature within tolerance, or None.

    RT is only used when both the entry and the feature provide it.
    """
    best_idx = None
    best_score = None
    for idx, (_compound, fmz, frt) in enumerate(features):
        if abs(fmz - mz) > mz_tol:
            continue
        if rt is not None and frt is not None and abs(frt - rt) > rt_tol:
            continue
        score = abs(fmz - mz) + (abs(frt - rt) if (rt is not None and frt is not None) else 0.0)
        if best_score is None or score < best_score:
            best_score = score
            best_idx = idx
    return best_idx


# --------------------------------------------------------------------------- #
# GNPS2 re-indexing of fragmentation files
# --------------------------------------------------------------------------- #

_FEATURE_ID_RE = re.compile(r'^(FEATURE_ID|SCANS)=', re.IGNORECASE)


def reindex_fragments(peaktable_path, frag_path, out_path,
                      mz_tol=0.01, rt_tol=0.05):
    """Re-index an MSP/MGF file so entry IDs match peak-table row numbers.

    GNPS2 requires the fragmentation file's per-entry scan/feature id to equal
    the (1-based) row number of the matching feature in the peak table. Each
    fragmentation entry is matched to a feature by m/z (and RT when available),
    assigned that feature's row number, and the entries are written out sorted
    by row number.

    Args:
        peaktable_path: canonical/internal peak table whose row order defines IDs.
        frag_path: source ``.mgf`` or ``.msp`` file.
        out_path: destination path (same type as ``frag_path``).
        mz_tol: m/z match tolerance (Da).
        rt_tol: RT match tolerance (min).

    Returns:
        The number of fragmentation entries successfully matched and written.
    """
    features = _read_feature_table(peaktable_path)
    id_map = {str(compound): idx for idx, (compound, _mz, _rt) in enumerate(features)}
    entries = parse_fragments(frag_path)
    suffix = Path(frag_path).suffix.lower()

    matched = []
    for entry in entries:
        # Prefer an exact compound-id/comment match (robust to neutral-mass vs
        # adduct-m/z differences); fall back to m/z (+RT) tolerance matching.
        idx = id_map.get(str(entry.comment)) if entry.comment is not None else None
        if idx is None and entry.mz is not None:
            idx = _match_index(entry.mz, entry.rt, features, mz_tol, rt_tol)
        if idx is not None:
            matched.append((idx + 1, entry))  # 1-based row id

    matched.sort(key=lambda pair: pair[0])

    if suffix == '.mgf':
        _write_mgf(matched, out_path)
    else:
        _write_msp(matched, out_path)
    return len(matched)


def _write_mgf(matched, out_path):
    with open(out_path, 'w') as handle:
        for new_id, entry in matched:
            wrote_scans = False
            for line in entry.lines:
                if _FEATURE_ID_RE.match(line.strip()):
                    handle.write('SCANS=' + str(new_id) + '\n')
                    wrote_scans = True
                elif line.strip() == 'END IONS' and not wrote_scans:
                    handle.write('SCANS=' + str(new_id) + '\n')
                    handle.write(line if line.endswith('\n') else line + '\n')
                    wrote_scans = True
                else:
                    handle.write(line if line.endswith('\n') else line + '\n')
            handle.write('\n')


def _write_msp(matched, out_path):
    with open(out_path, 'w') as handle:
        for new_id, entry in matched:
            wrote_id = False
            for line in entry.lines:
                key = line.split(':', 1)[0].strip().upper()
                if key in ('NAME',):
                    handle.write(line if line.endswith('\n') else line + '\n')
                    handle.write('SCANNUMBER: ' + str(new_id) + '\n')
                    wrote_id = True
                elif key == 'SCANNUMBER':
                    continue  # replaced above
                else:
                    handle.write(line if line.endswith('\n') else line + '\n')
            if not wrote_id:
                handle.write('SCANNUMBER: ' + str(new_id) + '\n')
            handle.write('\n')


# --------------------------------------------------------------------------- #
# Writing the filtered peak table back in the source layout
# --------------------------------------------------------------------------- #

def filter_source_peaktable(source_path, filtered_internal_path, out_path,
                            mz_tol=0.01, rt_tol=0.01):
    """Write the *source* peak table filtered to surviving features.

    Rather than reconstructing each tool's layout, the surviving features (read
    from the internal ``_filtered.csv``) are matched back to rows of the
    original source file by m/z + RT, and only those rows are written out -- so
    the result keeps the exact original format, just filtered.

    Args:
        source_path: the original source peak table (any supported format).
        filtered_internal_path: MPACT's internal ``_filtered.csv``.
        out_path: destination CSV/TXT path.
        mz_tol / rt_tol: matching tolerances.

    Returns:
        Number of rows written.
    """
    fmt = detect_peaktable_format(source_path)
    sep = '\t' if Path(source_path).suffix.lower() == '.txt' else ','
    raw = pd.read_csv(source_path, sep=sep, header=None, dtype=str, engine='python')

    mz_col, rt_col, data_start = _source_mz_rt_columns(raw, fmt)
    if mz_col is None or rt_col is None:
        raise ValueError('Could not locate m/z and RT columns in source file (format=' + fmt + ')')

    surviving = _read_feature_table(filtered_internal_path)
    kept_ids = {str(cid) for cid, _mz, _rt in surviving}
    kept = [(mz, rt) for _id, mz, rt in surviving]

    keep_rows = list(range(data_start))  # header rows always kept
    for row_idx in range(data_start, len(raw)):
        # Fast path: exact compound-id match on column 0 (shared between the
        # source and internal tables); fall back to m/z (+RT) tolerance.
        if _first_cell(raw.iloc[row_idx, 0]) in kept_ids:
            keep_rows.append(row_idx)
            continue
        mz = _to_float(raw.iloc[row_idx, mz_col])
        rt = _to_float(raw.iloc[row_idx, rt_col])
        if mz is None:
            continue
        for fmz, frt in kept:
            if abs(fmz - mz) <= mz_tol and (frt is None or rt is None or abs(frt - rt) <= rt_tol):
                keep_rows.append(row_idx)
                break

    out = raw.iloc[keep_rows, :]
    out.to_csv(out_path, sep=sep, header=False, index=False)
    return len(keep_rows) - data_start


def _source_mz_rt_columns(raw, fmt):
    """Locate (mz_col, rt_col, data_start_row) in a raw source table."""
    if fmt == MZMINE:
        header = raw.iloc[0]
        mz_col = rt_col = None
        for idx, val in enumerate(header):
            low = _first_cell(val).lower()
            if low == 'row m/z':
                mz_col = idx
            elif low == 'row retention time':
                rt_col = idx
        return mz_col, rt_col, 1

    if fmt in (PROGENESIS, METABOSCAPE):
        # Canonical: header row holding 'm/z' and 'Retention time (min)'.
        for row_idx in range(min(4, len(raw))):
            row = raw.iloc[row_idx]
            mz_col = rt_col = None
            for idx, val in enumerate(row):
                cell = _first_cell(val)
                if cell == 'm/z':
                    mz_col = idx
                elif cell == 'Retention time (min)':
                    rt_col = idx
            if mz_col is not None and rt_col is not None:
                return mz_col, rt_col, row_idx + 1
        return None, None, 0

    if fmt == MSDIAL:
        for row_idx in range(min(8, len(raw))):
            row = raw.iloc[row_idx]
            mz_col = rt_col = None
            for idx, val in enumerate(row):
                cell = _first_cell(val)
                if cell == 'Average Mz':
                    mz_col = idx
                elif cell == 'Average Rt(min)':
                    rt_col = idx
            if mz_col is not None and rt_col is not None:
                return mz_col, rt_col, row_idx + 1
        return None, None, 0

    return None, None, 0
