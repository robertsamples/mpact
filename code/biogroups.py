"""
MPACT
Copyright 2022, Robert M. Samples, Sara P. Puckett, and Marcy J. Balunas

Qt-free derivation of the biological-group list from the metadata/sample
list/peak-table files, extracted from MainWindow.getgroups() so the
metadata-join logic (and its duplicate-injection/missing-metadata edge
cases) can be unit-tested without a GUI.
"""

import pandas as pd


def compute_biological_groups(extractmetadatafilename, samplelistfilename, filename):
    """Derive the distinct biological groups present in a peak table.

    Returns ``(groups, unresolved)``: ``groups`` is the list of distinct
    ``Biological_Group`` values (first-seen order) found across every
    injection column in the peak table; ``unresolved`` is the list of
    injection names that had no matching row in the joined metadata (e.g. a
    sample list/metadata mismatch).

    Raises ``ValueError`` (chained from the original exception) if the
    metadata/sample-list join itself fails -- callers report that to the
    user the same way as before, just without this function needing to know
    about Qt.
    """
    extractmetadata = pd.read_csv(extractmetadatafilename, sep=',', header=[0], index_col=None)
    samplelist = pd.read_csv(samplelistfilename, sep=',', header=[0], index_col=None)

    try:
        combinedmetadata = extractmetadata.set_index('Sample_Code').join(samplelist.set_index('Sample_Code')) \
                                                .reset_index().set_index('Injection')
    except Exception as exc:
        raise ValueError('Data read failure: Check input files') from exc

    msdata = pd.read_csv(filename, sep=',', header=None, index_col=[0, 1, 2], low_memory=False)
    groups = []
    unresolved = []

    for elem in msdata.iloc[2]:
        try:
            biolgroup = combinedmetadata.loc[elem, 'Biological_Group']
            if isinstance(biolgroup, pd.Series):
                # A duplicate Injection entry in the joined metadata makes
                # .loc[] return a Series instead of a scalar -- take the
                # first rather than letting the `in`/append below raise on
                # an ambiguous Series truth value.
                biolgroup = biolgroup.iloc[0]
            if biolgroup not in groups:
                groups.append(biolgroup)
        except Exception:
            unresolved.append(elem)

    return groups, unresolved
