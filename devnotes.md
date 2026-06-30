# devnotes.md — MPACT

MPACT is a PyQt5 desktop tool for mass-spectrometry / natural-products data analysis. 
Code runs in ./code.

## Running the app

- **Windows:** double-click the `MPACT` shortcut, or `run.bat`. Non-fatal errors
  terminate the program when launched this way. Both resolve paths relative to
  the repo's actual location (`%~dp0` in `run.bat`) — if you move/re-clone the
  repo, nothing needs editing. `run.bat` activates Anaconda then `cd`s into
  `code\` (not the repo root) before launching, since several relative paths
  in the app (`npatlas.tsv`, `compoundimages\`, etc.) assume that cwd.
- **Mac / debugging:** launch via Spyder (Anaconda) → open `code/main.py` (NOT
  `__main__.py`) → Run.
- Entry point is `code/main.py`. It must be run as `__main__`; `import main`
  fails standalone because `main` ⇄ `ui_functions` is a circular import that
  only resolves when run as a script. So the GUI stack cannot be imported
  headlessly.

## Portable build (no Anaconda required)

`build/mpact.spec` + `build/build_windows.bat` produce a standalone Windows
folder build via PyInstaller that end users can unzip and run without
installing Anaconda, Python, or any packages. See `build/README.md` for the
build steps and the reasoning behind the two non-obvious spec choices
(`--onedir` not `--onefile`; `contents_directory='.'`) and the dependency
gotchas (`epam.indigo`'s native DLL needs `collect_all`, not just
`hiddenimports`; `requirements.txt` intentionally allows NumPy 2.x, unlike the
Anaconda-env cap below). Built and smoke-tested 2026-06-29 (offscreen launch,
all imports including indigo's native lib resolved cleanly); not yet run
against real data end-to-end by a human.

### Portable builds for Linux / Mac (not yet implemented)

Not attempted yet — logged here for whoever picks this up next, since the
Windows build surfaced platform-specific issues that will very likely recur:

- **PyInstaller itself is cross-platform but builds aren't** — you cannot
  cross-compile a Linux or Mac build from this Windows dev machine; each
  needs to be built (and ideally smoke-tested per `build/README.md`) on that
  OS. The existing `build/mpact.spec` should mostly carry over unchanged
  (`pathex`, `datas`, `hiddenimports` are all OS-agnostic Python), but the
  `epam.indigo` fix is the one item proven to be OS-specific: its DLL lookup
  path is `indigo/lib/<platform>/...` (`windows-x86_64` here) — confirm
  `epam.indigo` actually ships Linux (`linux-x86_64`) and Mac
  (`mac-x86_64`/`mac-arm64`) native libraries in its PyPI wheel for the
  target platform/architecture before assuming `collect_all('indigo')` just
  works there too; if not, that dependency may need a different
  installation path (e.g. a system package) on that OS.
- **macOS app bundling**: a plain PyInstaller onedir folder runs on Mac, but
  for something distributable as a normal double-clickable app, the spec
  would need a `BUNDLE(...)` step (PyInstaller supports this directly in the
  same spec file) to produce a `MPACT.app`. Unsigned/unnotarized `.app`
  bundles get Gatekeeper warnings on other people's Macs — fine for the
  team's own use, a blocker for wider distribution without an Apple
  Developer signing certificate.
- **Linux**: PyInstaller onedir output runs fine as a plain executable +
  folder (no installer step strictly required), but Qt's platform plugin
  selection is more failure-prone than on Windows/Mac — expect to need
  `QT_QPA_PLATFORM`/`xcb` library checks on whatever distro it's tested on
  first, the same class of issue the offscreen-platform smoke test exists to
  catch early (`build/README.md`).
- **`run.bat`'s cwd-relative-path assumption** (`code/main.py` reading
  `npatlas.tsv`/`compoundimages/` relative to cwd, not `__file__`) is
  platform-agnostic in principle — a Linux/Mac launcher script just needs
  the equivalent `cd` into the build's own folder before exec'ing the binary,
  same as `--contents-directory .` solves it for the Windows PyInstaller
  build.
- Given all of the above is unverified, treat any of it as a starting
  hypothesis to confirm on real Linux/Mac hardware, not a finished plan.

##  NumPy 2.x / Anaconda dependency hazard (read before touching deps)

This runs in a conda env whose pandas/matplotlib/scipy are compiled against
**NumPy 1.x**. A bare `pip install <pkg>` can upgrade NumPy to 2.x, which breaks
*every* conda-compiled module ("A module compiled using NumPy 1.x cannot be run
in NumPy 2.5.0"). If that happens, recover with:

```
pip install "numpy<2"
```

`code/importdependencies.py` `install()` now pins `numpy<2` and uses
`--upgrade-strategy only-if-needed` so auto-install can't break the env. Keep it
that way. Required deps (gate startup): `epam.indigo`→`indigo`, `UpSetPlot`→
`upsetplot`, `squarify`. Optional/perf (never gate): `fastcluster`.
`ensure_dependencies()` is called at the top of `main.py` before those imports.

## Architecture / file map

- **Generated — DO NOT EDIT** (regenerated from Qt Designer; overwritten on
  update): `ui_main.py`, `ui_main1.py`, `ui_featureinfo.py`, `ui_plotparam.py`,
  and the resource blobs `files.py`, `files_rc.py`. Despite its name,
  **`ui_functions.py` is hand-written** (the `UIFunctions` controller) and IS
  editable.
- **Hand-written app code:** `main.py` (MainWindow, run/save/load, db search),
  `plotting.py` (plot classes; base `ui_plot`), `filter.py`, `stats.py`,
  `MSFaST.py` (analysis driver `run_MSFaST` + `analysis_parameters`),
  `pvclust.py` (bootstrap dendrogram), `translators.py` (import/export),
  `mzmineimport.py` (format conversion), `getfragdb.py`, `mspwriter.py`.
- **Qt-free helper modules extracted from `MainWindow` methods** (so the
  logic is unit-testable without a GUI — `main.py` can't be imported
  standalone, see above): `csvcache.py` (memoizing `pd.read_csv` for static
  per-run output files), `groupsets.py` (Plot Feature Set MVC),
  `plotslots.py` (per-plot widget-state registry), `paramfields.py`
  (shared save/restore schema for simple `analysis_parameters` checkbox
  fields), `biogroups.py` (`getgroups()`'s metadata-join/group-derivation
  core), `dbsearch.py` (`fulldbsearch()`'s NPAtlas ppm-window matching
  core), `ordination.py` (PCA/NMDS/PLS-DA + technical-replicate collapsing
  + top-N loadings selection for the multivariate plot tab), `clusterpurity.py`
  (dendrogram branch-purity coloring for the dendrogram tab). Each
  corresponding `MainWindow` method is now a thin wrapper: call the module
  function, then apply the result to widgets/`self`.
- **Runtime widget substitution into a Designer placeholder** is an
  established pattern here, not a one-off — `plotting.py` does it for every
  matplotlib canvas (inserted into a Designer-created `QFrame`), and
  `searchtree.py`'s `SearchTreePanel` does it for the feature-search tab's
  tree: it removes `ui_main.py`'s Designer-created `QTreeWidget` from its
  layout at runtime and inserts a real `QTreeView` + per-column filter bar
  in the same slot. No generated file is edited; re-opening the `.ui` file
  in Qt Designer would still show the old `QTreeWidget`, only the running
  app uses the substitute. `SearchTreePanel`'s filter-matching rules
  (per-column text/numeric-range/multi-select-category) live in the
  Qt-free, unit-tested `treefilters.py`; `searchtree.py` is just the Qt
  model/proxy/widget plumbing wired to it. `MainWindow.fillfttree()`/
  `on_tree_item_selection_changed()` talk to `self.searchtree` (`set_rows()`/
  `selected_compound()`) instead of `QTreeWidgetItem`s now.
- **Canonical peak table** (what MPACT consumes; Progenesis = native): CSV, 3
  header rows; row 2 = `Compound,m/z,Retention time (min),<injections...>`;
  col0 = `RT_mz` id, col1 = m/z, col2 = RT. Rows 0–1 are overwritten by
  `MSFaST.importdata` from the metadata CSVs.

## Threading model (do not break)

`run_MSFaST` is Qt-free and runs on a `QThread` worker (`AnalysisWorker` in
`main.py`) so the GUI stays responsive during the heavy compute.
`MainWindow.run_analysis` reads widgets on the main thread → starts the worker →
`_on_compute_finished` → `_finish_analysis` does ALL matplotlib/Qt plotting on
the **main thread** (matplotlib is not thread-safe). Never create Qt/matplotlib
objects on the worker thread.

## Importer / translator framework (`translators.py`)

Qt-free, unit-tested. `detect_peaktable_format`, `parse_msp`/`parse_mgf`
(→`FragmentEntry`), `reindex_fragments` (GNPS2: matches fragments to peak-table
rows by **compound-id first**, then m/z+RT — Progenesis MSP stores neutral mass,
not adduct m/z), `filter_source_peaktable` (row-subset the source to survivors).
`mzmineimport.format_check` delegates detection here and reports errors (no more
silent `except: pass`). `MainWindow.export_filtered_outputs()` auto-runs at the
end of `_finish_analysis`, writing `_filtered_source.<ext>` + `<frag>_reindexed`.
Source-format export is Progenesis-only (other formats are converted in place at
import).

## Testing

Headless unit tests in `code/tests/` (mostly Qt-free logic, plus standalone
Qt widget/model tests that don't import `main.py` — see below). Run:

```
python -m pytest code/tests -q
```

Covers `filter`, `stats`, `importdependencies`, `translators`, `groupsets`,
`searchtree`, `ordination`, `clusterpurity`. Add tests here for any new
Qt-free logic.

`conftest.py` sets `QT_QPA_PLATFORM=offscreen` and provides a session-scoped
`qapp` fixture: PyQt5 widgets/models/signals *can* be exercised headlessly via
Qt's offscreen platform plugin, as long as the module under test doesn't
import `main.py` (the documented `main`⇄`ui_functions` circular import is a
separate, unrelated constraint and still applies — `main.py` itself still
can't be imported standalone). `searchtree.py`'s tests are a real example:
they construct an actual `QTreeWidget` in an actual layout and exercise the
full runtime widget swap, not just the Qt-free filter logic underneath it.
Full end-to-end GUI behaviour (anything that requires `main.py`/`MainWindow`)
still must be verified by running the app.

## Groupsets (`groupsets.py`)

The user-defined "groupsets" (colour/presence-absence rules for plotting
features) are an MVC model: `GroupSet` (data) + `GroupSetModel` (collection,
bounds-safe selection, CRUD) + `build_query_dict()` (the {descriptive_name:
GroupSet} mapping MSFaST/plotting consume). `MainWindow.groupsetmodel` replaces
the old bare `self.querys` list + `self.selset` index. `ui_functions.py`'s
`addgroup`/`removegroup`/`updatesets`/`updategroups`/`writegroups`/
`colour_picker1` are thin view-sync controllers over the model — same method
names, so call sites elsewhere didn't change. `main.py`'s `query` class is kept
ONLY as the unpickle target for old `.mpct` files (pickled by qualified name
`main.query`); `GroupSet.from_legacy`/`GroupSetModel.from_legacy_list` convert
on load. UI uses `QListWidget` (generated, off-limits), so this is not a true
`QAbstractListModel`/`QListView` setup — the "view" side is the existing
hand-written widget-sync code in `ui_functions.py`, kept thin.

## Multivariate ordination plot (`plotting.plot_ordination`, `ordination.py`)

What used to be called "PCA" (`plot_PCA`, `checkBox_pca`/`btn_pca`'s old
tooltip) actually only ran NMDS, with a PCA rotation applied to the NMDS
coordinates purely to orient the axes — not a second ordination of the
original features. `plot_ordination` now genuinely supports PCA, NMDS, and
PLS-DA, switchable via a combo-box bar inserted above the plot canvas (same
runtime widget-substitution pattern as `searchtree.py`'s filter bar — see
above), plus a Scores/Loadings view toggle. The math lives in the Qt-free
`ordination.py` (unit-tested in `tests/test_ordination.py`); `plotting.py`
only handles the combo boxes, axes, and pick events.

- **Save-file compatibility preserved on purpose**: `analysis_params.PCA`
  and `checkBox_pca`'s objectName are unchanged (still pickled into `.mpct`
  saves) — only the visible checkbox text/tooltip changed (set at runtime in
  `MainWindow.__init__`, same mechanism as `label_credits.setText`). Only
  the hand-written class name (`plot_PCA` → `plot_ordination`) changed,
  since that's never pickled.
- **"Collapse Technical Replicates" used to be dead** (`plotting.py` had
  `parent.collapsereps = False#parent.dialog.ui.checkBox_collapsereps.isChecked()`
  — hardcoded off, the real read commented out). Now wired for real via
  `ordination.load_ordination_matrix(..., collapse_replicates=...)`. The
  collapse logic itself (average technical replicates/Injections, keep
  biological replicates/Samples distinct) was ported verbatim from the
  original rather than rewritten — its header-relabeling-via-CSV-round-trip
  is easy to get subtly wrong by inspection, so it's verified empirically
  instead (`test_ordination.py`'s synthetic-replicate-structure test, cross-
  checked against real example data with a scratch script during
  development).
- **The checkbox itself later moved off the global plot-config dialog**:
  `checkBox_collapsereps` only ever affected this one plot, so it's now
  `plot_ordination`'s own "Collapse Replicates" checkbox (`self.collapse_replicates`,
  default `True`) in its switcher bar, same move as `plot_dendrogram`'s
  "Bootstrap" checkbox (see the dendrogram section). The dialog checkbox's
  containing frame (`frame_2`) is hidden at runtime rather than edited out
  of generated code. This field was never in `paramfields.CHECKBOX_FIELDS`
  (wasn't pickled before either), so no save/load behavior changed.
- **Loadings view and high-dimensional data**: thousands of features can't
  all be drawn legibly, so only the top-25 by loading-vector magnitude are
  shown by default (`ordination.top_loadings()`). Whichever feature is
  currently highlighted elsewhere in the app (`MainWindow.pickedfeature`) is
  always included regardless of magnitude — `plot_ordination.highlight_loading()`,
  called from `MainWindow._refresh_highlight()`, follows the same
  pre-create-an-empty-artist/update-via-`set_data()` convention every other
  plot's highlight marker already uses. This is a *feature* highlight
  (Loadings view), a different concept from the Scores view's existing
  *sample* highlight (`parent.pickedsample`, set by clicking a sample point)
  — the two views show different kinds of points and were never the same
  selection concept.
- **NMDS has no linear feature loadings** (it's a rank-based embedding, not
  a linear projection) — its Loadings view uses `ordination.nmds_loading_proxy()`,
  per-feature correlation with each NMDS axis (the standard ecology "vector
  fitting"/`envfit` approach), not real loadings. **NMDS's axis labels don't
  show percent-explained at all** (just "NMDS1"/"NMDS2") — NMDS is a
  rank-based embedding, not a linear decomposition of the feature space, so
  it doesn't canonically have a %-variance-explained quantity the way
  PCA/PLS-DA do; the plot title shows stress instead, the conventional NMDS
  fit-quality metric.
- **PCA/PLS-DA are autoscaled** (`ordination.autoscale()`: mean-center +
  scale each feature to unit variance) before fitting — without this, raw
  intensities (confirmed on real example data: feature standard deviations
  ranged from ~1.8 to ~10,000, a ~5800x spread) let a handful of
  high-abundance features dominate both the apparent explained-variance and
  the loadings, drowning out features that actually separate the
  biological groups but happen to have lower raw intensity. This is the
  standard chemometrics pretreatment for PCA/PLS-DA on this kind of data;
  NMDS is deliberately NOT autoscaled (its Bray-Curtis dissimilarity is
  conventionally computed on raw/relative abundances).
- **PLS-DA's explained-variance gotcha**: `sklearn.cross_decomposition.PLSRegression`
  defaults to `scale=True` (standardizes X internally), which -- before
  autoscaling was added -- silently produced explained-variance ratios off
  by ~6 orders of magnitude when compared against unscaled total variance,
  caught only by running against real data, not by inspection. Fixed with
  `scale=False` and autoscaling `x` ourselves first instead (matching PCA's
  treatment, and avoiding PLSRegression's `scale=True` also incorrectly
  scaling the 0/1 group-membership dummy columns).
- **Loadings-view rendering gotchas** (both only surfaced by checking
  against real data, not by inspection): (1) `ax.annotate()`-drawn arrows
  don't reliably participate in matplotlib's autoscale the way
  `ax.scatter()`/`ax.plot()` do — confirmed empirically that plotted arrow
  tips could fall outside the axis' auto-picked view limits — so
  `plot_ordination._plot_loadings()` now sets `ax.set_xlim`/`set_ylim`
  explicitly from the actually-plotted subset's coordinates, symmetric
  around 0. (2) `ordination.top_loadings()` must be called with only the 2
  displayed components (`loadings.iloc[:, :2]`), not the full (up to
  10-component) loadings — ranking by overall magnitude across all
  components could let a feature into the "top 25" purely from a large
  contribution to some unplotted component while barely showing up in the
  actual PC1-vs-PC2 view, displacing a feature that's genuinely prominent
  there.
- **OPLS-DA intentionally not implemented**: no native scikit-learn support;
  the alternatives (the unmaintained `pyopls` package, or a from-scratch
  orthogonal-signal-correction implementation) are both riskier than
  shipping PCA/NMDS/PLS-DA without a reference dataset to validate against.
  Logged here as the next ordination method to add if ever revisited, not
  started.

## Dendrogram purity coloring (`plotting.plot_dendrogram`, `clusterpurity.py`)

The dendrogram tab has a switcher bar (same runtime-widget-substitution
pattern as `plot_ordination`'s method/view bar) with two combo boxes (View,
Color) and two checkboxes (Bootstrap, Use Sample/Group Names -- both
documented further down, formerly/newly local to this tab respectively):

- **View** — which leaves to cluster:
  - **Technical Replicates** (default — matches the tab's original
    behaviour): every Injection is its own leaf, purity judged against
    Sample membership — a tight monophyletic clump means that sample's
    replicates agree.
  - **Biological Replicates**: technical replicates are averaged first
    (same `ordination.load_ordination_matrix(..., collapse_replicates=True)`
    used by the multivariate tab's checkbox), so leaves are Samples, purity
    judged against Biolgroup instead — a monophyletic clade means a whole
    biological group's samples cluster together before meeting another
    group, i.e. the groups are separable.
- **Color** — how to render purity:
  - **Purity** (default): green wherever a branch's leaves are entirely one
    group (correctly clustered), magenta wherever a branch mixes more than
    one group (polyphyletic) — a QC judgment visible at a glance rather than
    read off leaf labels one at a time. The plot title reports
    `n_pure/n_total` (e.g. "7/9 samples' replicates clustered together",
    "3/3 biological groups separable") via `clusterpurity.purity_summary()`.
  - **None**: plain black, no title — deliberately reproduces the tab's
    appearance from *before* purity coloring existed (there was no title at
    all previously), for anyone who just wants the dendrogram shape without
    the QC overlay. Implemented as `link_color_func=None` with
    `color_threshold=0` still set (dropping `color_threshold=0` here was a
    real regression caught while testing: without it, scipy falls back to
    its own default 0.7-of-max-height threshold and a multi-color palette
    instead of plain black).

Both views' purity math is the same Qt-free linkage-traversal logic in
`clusterpurity.py`, unit-tested in `tests/test_clusterpurity.py`.

- **`false_color` marks proven non-monophyly (overlap), not "any impure
  merge"**: two earlier attempts both got this wrong in opposite directions.
  First, every impure merge was colored `false_color`, including every
  ancestor above a single mixing event all the way to the root -- since
  almost any real dataset has *some* mixing somewhere, this painted most of
  the tree's upper structure regardless of how localized the problem was.
  The second attempt ("impure but at least one child was pure = bridge =
  `false_color`, both children already impure = neutral") fixed the worst
  of the cascading but still mis-colored real data: it could still mark a
  high-level merge `false_color` merely because one side happened to be a
  single freshly-introduced pure clade, *and* it could miss real tangles
  where two already-impure children share a label without one side being
  trivially pure.

  `purity_link_color_func()` now compares the two children's label sets
  directly at each merge:
  - identical and a single label -> monophyletic (`true_color`/green).
  - **disjoint** (no label in common) -> neutral (`neutral_color`/black) --
    a clean join of two regions that don't contradict each other, *even if
    one or both children are themselves impure from a different label's
    tangle further down*. This is what stops a low-level tangle from
    cascading: once a tangled label has nothing more of itself left to fold
    in, every merge above it only ever joins disjoint regions, so it goes
    back to black.
  - **overlap** (share >=1 label, without being identical-and-singleton) ->
    polyphyletic (`false_color`/magenta) -- definitive proof that some
    label's leaves are split across this exact merge (present on both
    sides), not just "still mixed from an earlier merge".

  Verified against the real example dataset's bootstrap dendrogram (the
  case that exposed both earlier bugs): only the two merges that actually
  re-unite a scattered sample's replicates (e.g. one sample's reps split
  into two non-sister sub-clades that only meet again higher up) render
  `false_color`; the higher-level merges joining that region with
  cleanly-resolved, unrelated samples stay black, same as a hand-built
  synthetic linkage (`tests/test_clusterpurity.py`'s
  `_scattered_pair_linkage`) reproducing the same pattern deterministically.

  `true_color`/`false_color` default to green/magenta, not the more
  conventional green/red: red-green colorblindness (the most common form)
  can't distinguish red from green, while magenta stays distinguishable
  from green under all common forms of color vision deficiency. (Changed
  from an original green/red default after user feedback; see
  `clusterpurity.py`'s `purity_link_color_func()` default args and
  `plotting.py`'s `plot_dendrogram.plot()` call site.)
- **Bootstrap is now a per-tab checkbox, not a global one**: the
  plot-config dialog's "Bootstrap Analysis" checkbox (`checkBox_bootstrap`)
  only ever affected this one plot, so it moved into `plot_dendrogram`'s own
  switcher bar (`self.bootstrap`, default `True` -- matching the effective
  startup default the old checkbox was forced to in `UIFunctions`, which
  differed from its own Designer-set default of `False`). The dialog
  checkbox's containing frame (`frame_bootstrap`) is hidden at runtime in
  `UIFunctions` rather than edited out of the generated `ui_plotparam.py`.
  `('bootstrap', ...)` was also dropped from `paramfields.CHECKBOX_FIELDS`,
  so it's no longer saved into `.mpct` files -- consistent with the
  dendrogram's other per-tab state (View, Color), none of which persist
  across save/load either.

- **Purity is a strict, whole-group check, not "any uniform subset"**: a
  label only counts as pure if *every* leaf carrying it ends up in one clade
  before that clade touches a different label — 2 of a Sample's 3 replicates
  merging together does NOT make that Sample pure if the third replicate
  clusters elsewhere. An earlier version of `purity_summary()` got this
  wrong (counted a label pure as soon as ANY uniform-label merge occurred,
  which is right for `purity_link_color_func`'s per-branch coloring but wrong
  for the whole-group summary count) — caught by a test built from a
  deliberately "rogue" planted point, not by inspection.
- **PvClust orientation gotcha**: `pvclust.PvClust` expects "variables x
  objects" (rows = the things bootstrapped over, i.e. features; it
  transposes internally before clustering the columns) — the *opposite*
  orientation from `scipy.cluster.hierarchy.linkage`, which expects "objects
  x variables". `plot_dendrogram.plot()` builds both orientations
  (`data_for_linkage`, `data_for_pvclust`) from the same scaled data rather
  than reusing one array, since which one is "transposed" flips between the
  Technical (features x injections is the natural read) and Biological
  (samples x features, from `load_ordination_matrix`) views.
- **`link_color_func` threaded through the bootstrap path too**:
  `PvClust.plot()` and the free function `pvclust.plot_dendrogram()` both
  gained a `link_color_func=None` passthrough parameter into their inner
  `scipy.cluster.hierarchy.dendrogram()` call, so the AU/BP bootstrap
  dendrogram gets the same purity coloring as the regular one (`scipy`'s own
  precedence rule: `link_color_func`, when given, overrides
  `color_threshold`/`above_threshold_color`).
- **Multiprocessing safety note (re-learned, not new)**: validating the
  bootstrap path's wiring during development used `parallel=False` and a
  tiny `nboot`, never `PvClust(..., parallel=True)` in an ad hoc script —
  `multiprocessing.Pool()` re-executes a script's top-level code in each
  spawned child on Windows unless the call site is guarded by
  `if __name__ == '__main__':` (the same class of hazard as the frozen-exe
  fork-bomb bug elsewhere in this file, just without needing
  `freeze_support()` specifically). The real app is fine — `main.py` already
  guards its entry point — but throwaway test scripts need the same
  discipline.
- **"Use Sample/Group Names" leaf labels** (`ordination.replicate_label_components()`):
  swaps the raw file/injection names for `<Biolgroup>_b<BioRep#>_s<TechRep#>`
  (Technical Replicates view) or `<Biolgroup>_b<BioRep#>` (Biological
  Replicates view -- no TechRep#, since that view already collapsed
  technical replicates), for when the real file names are long or
  uninformative. BioRep# is the 1-based rank of a Sample within its
  Biolgroup (first-seen order); TechRep# is the 1-based rank of an
  Injection within its Sample -- both numbers are assigned unconditionally,
  so a Biolgroup with only one Sample still shows `_b1` and a Sample with
  only one Injection still shows `_s1` (no special-casing needed for either
  edge case, verified in `test_ordination.py`). This only changes the
  `labels=` argument passed to `dendrogram()`/`PvClust.plot()` -- the
  underlying data orientation, clustering, and purity-coloring lookups all
  still key off the raw names internally.
- **AU/BP label scaling, regardless of leaf count**: the bootstrap
  dendrogram's per-node AU/BP annotations used to be positioned with a
  fixed x-shift in *icoord* units (e.g. `x-7`, with a separate `x-10` for
  3-digit "100" values). icoord spacing is always 10 units per leaf no
  matter how many leaves there are, but the axes' actual pixel width isn't
  -- with more leaves squeezed into the same plot width, each icoord unit
  maps to fewer and fewer pixels, so that fixed icoord offset shrinks to an
  ever-smaller *pixel* gap, eventually merging "AU"/"BP" into "AUBP" (and
  every node's AU/BP pair into illegible overlapping text) once there are
  enough leaves. Fixed by switching to `ax.annotate(..., xytext=(±2, 0),
  textcoords='offset points', ha='right'/'left')`: a constant gap in
  *points* (real pixels-at-a-given-DPI) stays a constant gap regardless of
  icoord scale or leaf count, and `ha='right'`/`ha='left'` anchoring makes
  the old digit-width-dependent branching (-7 vs -10) unnecessary entirely.
  Per-node fontsize is also now scaled down as leaf count grows
  (`max(5, min(8, 140 / n_leaves))`) so neighbouring *different* nodes'
  labels -- which do have a fixed minimum icoord (and therefore pixel)
  separation -- don't run into each other either. Verified by rendering
  both a 6-leaf and a 27-leaf synthetic tree (matching the real 9-sample
  x3-techrep dataset's leaf count) and visually confirming no overlap in
  either. Also removed a `plt.figure(figsize=(12,8))`/`plt.tight_layout()`
  pair that created and immediately abandoned an unused Figure on every
  call -- it never affected the actual target `axis` and was a real (if
  small) per-redraw resource leak.

## Treemap / upset plot canvases (`plotting.plot_treemap`, `plotting.plot_upset`)

These two tabs used to be the only plots in the app that weren't real
matplotlib canvases: `gen_treemap`/`gen_upsetplt` (free functions, not
`ui_plot` subclasses) drew with `squarify`/`upsetplot`, `savefig()`'d a PNG
to the repo root (`treemap.png`/`test_upsetplt.png`), then loaded that PNG
into a `QPixmap` on the Designer-placed `label_treemap`/`label_upset`. That
meant no zoom/pan/save-at-resolution toolbar, a flat raster rewritten from
scratch on every run, and files left sitting at the repo root.

Both are now `ui_plot`-style classes (`plot_treemap`/`plot_upset`) drawing
directly onto a persistent `FigureCanvas`, wired into `MainWindow._generate_plots()`
via `_create_or_reset()` exactly like every other plot — so they're created
once and `.reset()` afterward, regenerating on both a fresh analysis run and
the dialog's "Apply" button (`regenerateplts()`), the same as every other
plot. They previously were NOT regenerated by Apply at all (`gen_treemap`/
`gen_upsetplt` were only ever called once, directly from `_finish_analysis`)
— a small behavior change, but one that brings them in line with how every
other plot already worked, not a new inconsistency.

- **`frame_treemap`/`frame_upset` needed a different substitution trick**:
  unlike most plot frames (empty in Designer, so `ui_plot.__init__` can just
  call `frame.setLayout(...)`), these two already have a Designer-built
  layout holding the old placeholder `QLabel` — Qt refuses `setLayout()` on
  a frame that already has one. `plotting._detach_placeholder_widget()`
  removes the old label and reparents the old layout onto a throwaway
  widget (the standard Qt "delete this layout" trick) before the normal
  `ui_plot.__init__`/manual canvas setup runs — same runtime
  widget-substitution pattern as `searchtree.py`'s filter-bar swap, just
  with an extra detach step first. Verified headlessly (offscreen Qt) that
  this doesn't raise and the frame ends up with exactly the new layout.
- **`plot_upset` doesn't subclass `ui_plot`**, same as `plot_heatmap` and for
  the same reason: `upsetplot.plot()` lays out several axes (matrix,
  totals, intersections, shading) via its own gridspec on whatever figure
  it's given — there's no single "ax" to hand callers the way every
  scatter/line plot here has. Unlike `plot_heatmap` (which has to transplant
  axes from a brand-new seaborn figure onto the persistent one, since
  `sns.clustermap()` doesn't accept an existing figure), `upsetplot.plot()`
  takes a `fig=` kwarg directly — `reset()` just `fig.clf()`s and re-plots
  onto the same figure, no axes-transplant needed.
- **Verified against real data, not just import-checked**: a throwaway
  headless-Qt script built fake Designer-style frames (QFrame + layout +
  placeholder QLabel, matching `frame_treemap`/`frame_upset`'s actual
  structure), ran both classes against the real example dataset, and
  asserted (1) the new canvas/toolbar actually replaced the old layout, (2)
  axes count is stable across `.reset()` calls (not growing — would mean
  old axes/figures were leaking), and (3) no PNG got written to disk by
  either plot anymore.

## Sample correlation matrix (`plotting.plot_samplecorr`, `ordination.similarity_matrix`)

Used to be a hardcoded Spearman-only heatmap with technical replicates
always pre-averaged and no way to relabel the raw injection/sample names.
Now has a Method (Spearman/Jaccard/Bray-Curtis) switcher, a View
(Biological Replicates/Individual Injections/Biological Groups) switcher,
and a "Use Sample/Group Names" checkbox — same nomenclature and
`ordination.replicate_label_components()` reuse as `plot_dendrogram`'s.

- **`ordination.similarity_matrix(x, method)`** is the new Qt-free backend
  (covered by `test_ordination.py`): `x` is samples x features, same
  convention as `run_pca`/`run_nmds`/`run_plsda`. Spearman is
  `x.transpose().corr(method='spearman')`; Jaccard/Bray-Curtis go through
  `sklearn.metrics.pairwise_distances` (`metric='jaccard'`/`'braycurtis'`)
  and return `1 - distance`. Jaccard is computed on `x > 0` (presence/
  absence of detection, ignoring abundance) — deliberately *not* derived
  from the groupset query-dict machinery the user floated as a possible
  source, since that's per-feature-list bookkeeping for the UpSet/treemap
  tabs, a different concept from per-sample/group detection.
  Pearson/Kendall were considered and rejected: Pearson assumes
  normally-distributed abundances (the wrong fit for heavy-tailed LC-MS
  intensities, same reasoning that makes Spearman the established choice
  here), Kendall is a slower, largely redundant rank-correlation
  alternative to Spearman.
- **Controls live in the *shared* `frame_12`/`horizontalLayout_25` nav
  bar** (the one holding the pre-existing `btn_upsetplt`/`btn_samplecorr`
  buttons that switch `stackedWidget_grpanalysis`), not in this plot's own
  canvas frame — unlike `plot_dendrogram`/`plot_ordination`'s per-canvas
  switcher bars, these controls are specific to the Sample Correlations
  page but the nav bar is shared with the unrelated UpSet Plot page.
  `plot_samplecorr._build_grpanalysis_controls()` appends a stretch then
  its own control widget onto the existing layout (no `ui_main.py` edit)
  so the new controls sit to the right of the two buttons and the bar
  stays a single row regardless of window width.
- **Greying out on the UpSet Plot tab**: `ui_functions.py`'s
  `btn_upsetplt`/`btn_samplecorr` click handlers used to call
  `stackedWidget_grpanalysis.setCurrentIndex()` directly; they now route
  through `UIFunctions.switch_grpanalysis_tab(self, idx)`, which also calls
  `self.samplecorr.set_controls_enabled(idx == 1)` (guarded by
  `getattr(self, 'samplecorr', None) is not None`, since this can fire
  before any analysis has run and created the plot object). The Designer
  default for `stackedWidget_grpanalysis` is already index 1
  (Sample Correlations), so the controls start enabled, matching the
  default active page.
- **View → row/column construction**: "Biological Replicates" and
  "Individual Injections" reuse `ordination.load_ordination_matrix()`
  exactly like `plot_dendrogram` does (`collapse_replicates=True`/`False`).
  "Biological Groups" takes the collapsed (Biological Replicates) matrix
  and does one more `x.groupby(biolgroup).mean()` to average across
  biological replicates too — deliberately not a third mode inside
  `load_ordination_matrix` itself, since it's a trivial one-line reduction
  of an already-correct, already-tested intermediate result.
- **Heatmap `vmin` is `0` for all three methods**, including Spearman.
  Spearman is mathematically capable of going negative, but real
  sample-vs-sample correlations in this kind of data cluster tightly
  positive (e.g. 0.7-1.0) — a `-1..1` scale (tried first) compressed all of
  that meaningful variation into a sliver of the colour range, making the
  heatmap look uniformly dark/uninformative. `0..1` keeps the full colour
  range usable for the variation that actually occurs.
- Dropped the dead `iondict = cached_read_csv(...)` read that was never
  actually used by the old `plot()` body.

## Conventions

- Don't edit the generated UI files (above). Put behaviour in `main.py` /
  `ui_functions.py` / `plotting.py`.
- Plot generation is wrapped by `MainWindow.safe_generate` so one failing plot
  doesn't abort the rest. `.mpct` saves are atomic (temp file + `os.replace`)
  with per-component guards (`write_save`).
- `loadsession` restores each parameter independently (a bad/missing field
  can't abort the rest); add new analysis params to BOTH `enumerate_inputs`
  (save) and `loadsession` (restore).
- Plot objects (`self.ftplt`, `self.kmd`, `self.spec`, ...) are created the
  first time they're needed and `.reset()` afterwards, via
  `MainWindow._create_or_reset()` / `_generate_plots()` — never gate
  create-vs-reset on `self.analysisrun`, since an optional output (frag file,
  KMD, PCA...) can newly turn on in the same session after being off for a
  prior dataset, and the object would never have been created.
- Use `MainWindow._refresh_highlight()` (not `highlight_feature()`) to redraw
  the currently-selected feature's displays without changing selection state
  (e.g. on a feature-info dialog tab switch). `highlight_feature(newfeature)`
  is for real selection events and toggles the highlight off if the same
  feature is clicked twice — calling it with the already-selected feature
  re-triggers that toggle, which is a bug, not a refresh.
- Every matplotlib `pick_event` handler (`ui_plot.onpick`, heatmap's
  `onpick8`, PCA's `picksample`) must call `plotting._is_duplicate_pick(parent,
  event)` first and bail if it returns True. Matplotlib fires one pick_event
  per artist that registers a hit, not one per click — a feature plotted in
  more than one groupset/colour layer (separate `ax.scatter()` calls at the
  same coordinates) otherwise fires the handler twice for one click, and the
  second call's toggle-off logic undoes the first call's selection.
- `importdependencies.checkdep()` is silent when nothing needs installing —
  it runs on every launch (including every Spyder "Run File", which
  re-executes `main.py`'s top level), so don't reintroduce a routine
  "checking/already installed" print; only report actual installs/failures.
- Per-plot widget state (`self.fig`/`self.canvas`/`self.pltlayout`/
  `self.toolbar`/`self.ax`/`self.highlight`) is backed by one
  `PlotSlotRegistry` (`plotslots.py`), not six independent dicts — but the
  six attributes are still dict-like views for backward compatibility, so
  existing `self.canvas['heatmap']`-style call sites don't need to change.
  When adding a 7th piece of per-plot state, add it to `PlotSlot`/`FIELDS`
  in `plotslots.py` rather than a new bare dict on `MainWindow`.

## Known refactor backlog

Lower-priority siblings of the groupset/plot-slot refactors above:

- ~~`self.groups` can silently end up incomplete~~ — **done.** `getgroups()`
  now catches the lookup itself (not just the append), handles a
  duplicate-`Injection`-row edge case that previously raised inside the
  `try`, and reports once via `self.error()` plus a single consolidated
  console line instead of either crashing or printing silently per row.
- ~~`enumerate_inputs`/`loadsession` have no shared schema~~ — **partially
  done.** `paramfields.py` now holds a shared table (`CHECKBOX_FIELDS`) for
  the plain 1:1 boolean toggle fields (`PCA`, `Dendrogram`, `bootstrap`,
  `MZRTplt`, `KMD`, `mdguide`, `FC`, `FC3Dplt`, `Ttest`, `Volcanoplt`,
  `FDR`) — adding one of *this kind* of field means one entry in
  `paramfields.py`, not a hand-mirrored line in both `main.py` and
  `ui_functions.py`. Deliberately left out of the table (and still
  hand-written): fields with backward-compat `getattr(..., default)`
  fallbacks, compound/derived fields (`statstgrps` from two combo boxes),
  and branching fields (blank-filter abs/rel, max-iso-shift combo-by-text) —
  abstracting those into a generic table would obscure logic that's
  clearer inline, for marginal duplication savings.
- ~~Redundant disk I/O on every click~~ — **done.** See `csvcache.py`.
- ~~`graphfilters` stringly-typed flag~~ — **done.** See
  `groupsets.normalize_graphfilters()`.
- **`exportgnps()` (`main.py:309-522`, ~210 lines) duplicates, less
  robustly, logic `translators.py` already has tested.** It hand-rolls
  GNPS/MGF compound matching via an O(n·m) pairwise RT/m·z tolerance loop
  with no compound-ID fast path — exactly what
  `translators.reindex_fragments`/`filter_source_peaktable` already do, by
  matching on compound ID first and falling back to tolerance, with test
  coverage. Not dead code: it's a separately-triggered button
  (`btn_exportgnps`) for a user-supplied GNPS-format file, distinct from
  the automatic `export_filtered_outputs()` export. Worth migrating to
  reuse the tested matching logic instead of maintaining two
  implementations — logged here, not yet started.
- **`run_MSFaST`'s blank-filter block re-reads `_formatted.csv` right
  after `importdata()` already had the identical data in memory**
  (`MSFaST.py`, the `if analysis_params.blnkfltr:` block) — same shape as
  the `importdata()`/`iondict` fix already done, but **tried this one and
  found it's the risky kind, not the safe kind**: `importdata()`'s
  in-memory frame is indexed by Compound (`index_col=[0]`); the blank-filter
  block needs Compound as a plain column instead (`index_col=None`).
  Deriving that via `.reset_index()` produces a *different* MultiIndex
  label for that column than `pd.read_csv` does (pandas' own
  `Unnamed: 0_level_0`/`Unnamed: 0_level_1` placeholder convention for
  blank multi-row header cells vs. `reset_index()`'s `'index'`/`''`/`''`).
  That wouldn't matter for a leaf computation, but this block's result
  gets written *back* into the canonical `_formatted.csv` — and 12
  separate places elsewhere read that file with `header=[2]` (flattening
  to the bottom header level), almost certainly expecting that column to
  come out named literally `'Compound'`. Getting the label wrong here
  would silently corrupt a shared, widely-read file rather than just one
  derived value. Not done; would need either a verified-safe way to
  reconstruct that exact column label, or confirmation that none of those
  12 read sites actually depend on it (not yet checked across all 12).
- **`iondict.csv` is read-modify-written as shared, accumulating state
  across `filter.py`'s `relationalfilter`/`cvfilter`/`decon` and
  `stats.py`'s `properr`/`runfc`/`runttest`** (7-8 read+write round-trips
  for one analysis run) — each function adds its own column(s), and
  genuinely needs the *previous* function's writes already applied, so
  this isn't a "redundant read of unchanged data" the way the importdata
  cases were; it's the disk literally being used as the shared mutable
  state connecting a chain of functions across two modules. Fixing this
  for real means threading an `iondict` DataFrame through every one of
  those function signatures (`filter.py`/`stats.py`/`MSFaST.py` call
  sites all change together) rather than a one-function patch — this is
  the "bigger, multi-session" item flagged in earlier architecture
  discussion, not a quick incremental slice. Logged for visibility, not
  started.

## Refactor status (Jun 2026)

All 5 stages done & validated: (1) save/load param round-trip fix, (2) QThread
worker so the UI doesn't freeze, (3) `fastcluster` optional accel +
`low_memory` warning cleanup, (4) translator framework + GNPS2 reindex +
auto-export, (5) MVC refactor for groupsets. Plus several user-reported
bugfixes: overplotted-feature pick-event duplication, a pandas-version-
dependent `plot_abund` crash (`sns.barplot(x="index", ...)`), dependency-check
console spam, and a stale `run.bat`/`MPACT.lnk` pointing at an old checkout
location. Earlier bugfixes: highlight toggling on feature-info tab switches,
and plot objects never (re)created when an optional output turns on
mid-session. Most recent follow-ups: heatmap W/S selection had no bounds
clamping (`mv_heatmap`, could crash or silently wrap past either end of the
feature list); the six per-plot dicts were consolidated into
`PlotSlotRegistry` (`plotslots.py`). 65 passing tests.
