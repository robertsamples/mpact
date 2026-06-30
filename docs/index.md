# MPACT

MPACT is a desktop tool for processing and visually exploring untargeted
mass-spectrometry / natural-products metabolomics data. It takes a peak
table (Progenesis QI, MZmine, MS-DIAL, or Bruker Metaboscape), a sample
list, and a metadata file, and turns them into a filtered, statistically
annotated dataset with a full suite of interactive plots: data-quality
review, group-level set/correlation analysis, hierarchical clustering,
multivariate analysis (NMDS), m/z-vs-RT and mass-defect views, volcano
plots, heatmaps, and per-feature spectral/database-match lookup.

This site covers installing and running MPACT, the file formats it
expects, the analysis and filtering options, and what each plot/tab
shows. If you're contributing to MPACT itself, see
[`devnotes.md`](https://github.com/robertsamples/mpact/blob/main/devnotes.md)
in the repo root.

## Where to start

- New to MPACT? Start with [Installation](installation.md) and
  [Getting Started](getting-started.md).
- Preparing input files? See [File Selection](user-guide/file-selection.md).
- Looking for what a specific plot means? See [Plots & Results](plots/data-review.md).
- Something broken? See [Troubleshooting](troubleshooting.md).

!!! note
    This documentation is adapted from the original MPACT user guide
    (2022) and updated to reflect the current codebase (mid-2026),
    including the import/export translator framework for MZmine/MS-DIAL/
    Metaboscape peak tables, the background-threaded analysis run, the
    groupset (Plot Feature Sets) editor, the multivariate ordination
    rework (PCA/NMDS/PLS-DA with scores and loadings views), and the
    dendrogram rework (purity coloring, view/bootstrap/label switchers).
    Some screenshots referenced in the original guide have not been
    re-captured yet — see
    [`devnotes.md`](https://github.com/robertsamples/mpact/blob/main/devnotes.md)
    if you'd like to contribute updated images.
