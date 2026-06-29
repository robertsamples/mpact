# Changelog

## 2026 refactor

- Background-threaded analysis run (`QThread` worker) so the GUI no longer
  freezes during the heavy compute step.
- New import/export translator framework (`translators.py`) with
  content-based peak-table format detection (Progenesis, MZmine, MS-DIAL,
  Metaboscape) and automatic GNPS2-compatible re-indexing of MSP/MGF
  fragment databases on export.
- MVC refactor of "Plot Feature Sets" (`groupsets.py`).
- Multivariate ordination rework: the former NMDS-only tab now supports
  **PCA**, **NMDS**, and **PLS-DA** via a method switcher bar, plus a
  scores/loadings view toggle and a loadings highlight synced to feature
  selection elsewhere in the app. PCA/PLS-DA features are autoscaled
  (mean-center + unit variance) before fitting; NMDS stays on raw
  abundances (conventional for Bray-Curtis). A new Qt-free backend
  (`ordination.py`) handles all ordination math and is covered by
  headless unit tests.
- Dendrogram rework:
  - Per-plot switcher bar for **View** (Technical vs Biological
    Replicates), **Color** (Purity / None), **Bootstrap**, and
    **Use Sample/Group Names**.
  - Purity coloring: green branches are cleanly within-sample/group;
    red branches mark the exact merge point where two groups' leaves
    overlap (proven non-monophyly). A title reports how many
    samples/groups are fully correctly clustered.
  - "Use Sample/Group Names" checkbox: replaces raw injection/file names
    with `<group>_b<#>_s<#>` (Technical Replicates view) or
    `<group>_b<#>` (Biological Replicates view).
  - Bootstrap can be toggled off for a faster undecorated dendrogram.
  - Fixed AU/BP annotation alignment: labels now stay a constant pixel
    gap apart regardless of leaf count, with leaf-count-scaled font size.
  - `fastcluster` optional acceleration for bootstrap linkage.
- UpSet and treemap plots are now rendered directly on a Qt canvas
  (replacing a PNG round-trip).
- Headless unit test suite (`code/tests/`).
- Hardened dependency installer to prevent NumPy 2.x environment breaks.
- Numerous bugfixes: Spearman double-colorbar on re-run, highlight not
  clearing on feature-info tab switches, plot objects not (re)created when
  an optional output newly turns on mid-session, overplotted features not
  selecting reliably, an abundance-plot crash under newer pandas/seaborn,
  and a stale launcher script/shortcut.

## Rev 23.05.15

- Added data export tab.
- Added GNPS peak table filtering functionality (experimental, only tested
  with FBMN export in MS-DIAL).
- Removed analysis info tab (later reinstated as part of Feature Info).

## Rev 23.02.26

- Added filtering of database hits by kingdom and genus.
- Dependency check now shows a window if any dependencies were installed,
  informing the user the kernel will need to restart.
- Bugfixes in group-set interactions and Bruker data import.

## Rev 23.02.19

- Added support for Bruker Metaboscape peak lists.
- Added support for MS-DIAL MSP files.
- Added MSP file writer to export in-source fragmentation patterns.
- Reduced lag when selecting features.
- Added multithreading when generating feature abundance plots.
- Fixed a bug where feature highlights weren't visible after data
  re-analysis.
- Fixed an issue preventing samples/injections from being selected and
  highlighted in NMDS plots.
- Database hit images are now saved to a folder to avoid re-rendering.
