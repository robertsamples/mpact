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

## Code review pass (dev branch, 2026-06-30)

Full read-through of every hand-written, Qt-free module plus the docs and
TODO block. The codebase is in good shape; findings were modest. Test count
is now **159 passing** (the count above is stale).

### Fixes applied on this branch (low-risk, test-validated)

- **`MSFaST.py` `analysisinfo.txt` decon label was a copy-paste bug**: the
  `if analysis_params.decon:` branch wrote "Features failing **blank**
  filtering" (a verbatim copy of the blank-filter line above it). Corrected
  to "Features failing in-source/deconvolution filtering". Confirmed against
  `main.py`'s parallel data-review summary writer (`_finish_analysis`,
  ~line 1208), which already labels the same quantity correctly as
  "in-source ion filtering" — so the two writers now agree. Also fixed two
  user-facing typos in the same file: "Runetime" -> "Runtime" and
  "PCA unfitlered" -> "PCA unfiltered". These are pure string-label changes,
  not the risky re-read logic the analysisinfo backlog item warns about.
- **`MSFaST.run_MSFaST` latent `NameError` on `groupionlists`**: it was
  only initialised inside `if analysis_params.grpave:`, but referenced
  unconditionally further down (the `groupionlists['cv'/'relfil'/'insource']`
  writes and the groups-column loop) and inside the blank-filter block.
  The GUI hardcodes `grpave = True` (`main.py:~1335`), so this never fired
  in practice, but a loaded session or a test with `grpave=False` would
  crash. Added a defensive `groupionlists = {}` next to `ionfilters = {}`.
  Behaviour unchanged when `grpave=True` (`parsionlists` reassigns it).
- **Stray debug CSVs written to the current working directory** (which is
  `code/` in the deployed app, per `run.bat`): `stats.py` wrote
  `msdata_teststats_test.csv` (a debug-named file, never read back) and
  `qdata.csv` (never read back — the canonical `-logq` goes into
  `iondict.csv`), and `ordination.py` wrote `averagepca.csv` (an internal
  collapse round-trip scratch file). All three now write into the run's
  output directory: `<stem>_teststats.csv`, `<stem>_qvalues.csv`, and
  `averagepca.csv` next to the input peak table respectively. The
  pre-existing leftover copies sitting untracked in `code/`
  (`qdata.csv`, `msdata_teststats_test.csv`, `averagepca.csv`) are now
  obsolete and safe to delete — they will no longer be regenerated there.
- **Dead code removal** (`stats.py` `groupave`): a per-injection
  `variance_values`/`stddev_values` was computed but never used (the
  technical/biological RSDs are derived from grouped means, not from it).
  Removing it made the entire sum-of-squares accumulation chain dead too
  (`sum_squares_list`/`sum_squares_chunk`/`all_sum_squares`/`sum_squares_df`),
  so that's gone as well — a small but real per-chunk optimization (drops a
  `(chunk ** 2).groupby(...).sum()` on every chunk of the formatted table).
  Validated by `tests/test_msfast_pipeline.py`, which runs the real
  `groupave` against the bundled example dataset. Also dropped an unused
  `from pathlib import Path` in `MSFaST.py`.

### Findings NOT changed (need a decision or live-GUI validation)

- **The "Or Groups" (`src`) control not being applied is intended, NOT a
  bug** (confirmed by the developer, 2026-06-30). The groupset editor has
  three lists — And (`listWidget_andgrps` -> `incl`, feature must be in all),
  Exclude (`listWidget_allgrps` -> `excl`, feature must not be in any), and
  Or (`listWidget_orgrps` -> `src`, the groups a feature is *allowed* to
  appear in). `MSFaST.groupset.__init__` deliberately filters only on `incl`
  and `excl`: a feature that already satisfies And/Exclude is a member of the
  groupset, and `src` ("allowed in") by design doesn't further remove it, so
  there's nothing for `src` to do at filter time. This matches the observed
  behaviour. (Earlier in this review pass it was mis-flagged as an inert
  control — that was wrong; leaving the note here so it isn't re-flagged.)
- **`mspwriter.convert_to_msp` num-peaks loop is fragile.** `for frags in
  sources: numpeaks = len(frags)` overwrites rather than accumulates, and
  assumes `sources` is a list-of-one-list. It happens to be correct for the
  only live caller (the decon path, where `ionmerge.sources == [[frag,...]]`),
  but would silently miscount if ever called on a `relationalfilter`-shaped
  merge (flat list of id strings) — `len(frags)` would then be a string
  length and the inner `for fragment in frags:` would iterate characters.
  Left as-is (single caller, wrapped in try/except), but worth hardening if
  the MSP writer is ever reused.
- ~~Docs repo-URL inconsistency~~ — **resolved.** Canonical repo is
  `github.com/robertsamples/mpact` (confirmed by the developer);
  `docs/installation.md`'s `git clone` line was corrected from `BalunasLab` to
  `robertsamples` to match `mkdocs.yml`/`docs/index.md`. `docs/index.md`'s
  stale "multivariate analysis (NMDS)" blurb was also updated to
  "(PCA/NMDS/PLS-DA)".

- **Two orphaned/broken scratch scripts in `code/`.**
  `npatlassearch.py` reads `npatlas.csv` (the real file is `npatlas.tsv`) at
  module top level and references an undefined `indigo`/`renderer`, so it
  would crash if ever imported/run — but nothing imports it.
  `masstdriver.py` is referenced only by a commented-out import in
  `ui_functions.py`. Both are dead leftover dev scratch, not part of the
  running app. Flagged rather than deleted (pre-existing files; the auto-mode
  classifier has blocked deleting UI-adjacent files before). Safe to remove
  once you confirm you don't want them as references.

### Already-logged items re-confirmed still open (see backlog above)

- `exportgnps()` duplicating `translators.reindex_fragments` matching logic.
- `iondict.csv` read-modify-write chain across `filter.py`/`stats.py`.
- `run_MSFaST`'s blank-filter `_formatted.csv` re-read (the "risky kind").
- Lazy per-tab plot updates; generalizing the `ui_plot` subclasses.

### Test-suite assessment

The suite is well-targeted and not redundant — each test guards a specific
behaviour or a previously-fixed bug (the PLS-DA `scale=` regression, the
replicate-collapse structure, the dendrogram purity edge cases, the
end-to-end pipeline). No tests are recommended for removal. Gaps worth
filling when convenient (all Qt-free, so headless-testable):
- `translators.reindex_fragments` / `filter_source_peaktable` end-to-end on
  the bundled MSP/MGF + peak tables (currently only smaller-unit coverage).
- `getfragdb.importfrag` format auto-detection (Progenesis vs MS-DIAL MSP).
- A `run_MSFaST` variant with `grpave=False`/minimal filters to lock in the
  `groupionlists` defensive-init fix above.
- `stats.runfc`/`runttest` numeric outputs (FC clamping, q-value monotonicity)
  against a tiny synthetic `iondict.csv`.

**Update (2026-06-30, second pass):** all four gaps above are now filled —
`tests/test_translators_e2e.py`, `tests/test_getfragdb.py`,
`tests/test_msfast_grpave_off.py`, `tests/test_stats_numeric.py`. Plus the
three new subsystems below ship with their own Qt-free tests
(`test_npatlasupdate.py`, `test_mpactupdate.py`, `test_crashreport.py`).

## Future feature dev plan (post-review, 2026-06-30)

Candidate features, ordered roughly by value-to-effort. None started; all
need the GUI runnable against real data to validate. Several already appear
in `main.py`'s TODO block — this is the triaged version.

1. ~~Wire up the "Or Groups" groupset constraint~~ — **withdrawn**: not a
   bug, the `src` "allowed in" semantics are intended (see finding above).
2. ~~Data-quality score / summary~~ — **partially done.** The score the TODO
   asked for already existed (Reproducibility / Skewness / Overall, from the
   AUC of the CV rarefaction curve), but the math was buried untested inside
   `plotting.prev_cv.plot()`. Extracted verbatim into the Qt-free, unit-tested
   `qualityscore.py` (`compute_cv_quality`), pinned against a copy of the
   original by `tests/test_qualityscore.py` so no displayed number changed;
   `prev_cv` is now a thin draw-the-result wrapper. (Also fixed a latent pandas
   FutureWarning in the extracted percentage assignment.) Remaining/optional:
   surface the score outside the CV tab (e.g. Data Review summary), and fold in
   the other available signals (per-group RSDs from `_summarydata.csv`, the
   dendrogram purity `n_pure/n_total`) if a richer composite is wanted -- both
   are scientific-design calls to make with the lab, not coded blind.
3. **OPLS-DA ordination method** (next item after the PCA/NMDS/PLS-DA rework,
   already deferred — see "Multivariate ordination plot"). Needs either the
   unmaintained `pyopls` or a from-scratch OSC implementation plus a
   reference dataset to validate against.
4. **Status-bar terminal/log viewer** (TODO). Replace the static status
   strings with a live log line + an expandable full-output pane. Mostly a
   Qt plumbing task (route the existing `print()` progress through a
   `QPlainTextEdit`/signal); no scientific risk.
5. **Additional databases beyond NPAtlas** (TODO: HMDB etc.). `dbsearch.py`
   is already a clean Qt-free ppm-window matcher taking an `atlas` DataFrame
   — adding a second source is mostly a loader + a column-name adapter, and
   the matching core is reusable as-is.
6. **`exportgnps()` migration onto `translators`** (backlog). Correctness +
   maintenance win, not a new feature: replace the ~210-line hand-rolled
   O(n·m) MGF matcher with the tested `reindex_fragments`/
   `filter_source_peaktable` path.
7. **Specificity/sensitivity & comparison-mode plots** (TODO, "likely items
   that need more thought"). Larger scientific-design questions; needs spec
   work with the lab before implementation.

## New subsystems (2026-06-30, second pass)

Three new Qt-free, unit-tested modules plus thin GUI wiring in `main.py`. The
cores are fully testable headlessly (network/git/dialog all injected); the
GUI wiring (`MainWindow._run_startup_checks`/`_check_atlas_freshness`/
`_check_app_update` and the `__main__` crash-dialog) is the only part that
needs a live launch to verify. **No new hard dependencies** — all three use
only the stdlib (`urllib`, `json`, `subprocess`, `webbrowser`, `platform`)
plus `packaging` (already present, with a tuple-comparison fallback), so
`requirements.txt`/the portable build are unaffected.

### NPAtlas auto-updater (`npatlasupdate.py`, `tests/test_npatlasupdate.py`)

On startup (deferred via `QTimer.singleShot` so the window paints first), if
`npatlas.tsv` is missing or its mtime is > 30 days old, the user is asked
whether to re-download it from
`https://www.npatlas.org/static/downloads/NPAtlas_download.tsv`. The download
streams to a temp file, is **validated** (header must contain the columns the
app uses — `compound_id`/`compound_m_plus_h`/`compound_m_plus_na`/
`compound_smiles`/`origin_type`/`genus`) and only then `os.replace`-d over the
existing file, so a server error page / partial transfer / network drop can
never clobber a working atlas.

- **Format decision (asked: would changing format help?): no — stay on TSV.**
  `main.py` reads the atlas with `pd.read_csv(sep='\t')` and `dbsearch` keys
  off the specific columns above; the published `NPAtlas_download.tsv` already
  has exactly those, so it's a drop-in. The `NPAtlas_download.json` is the
  same data in a nested shape that would need flattening before pandas/dbsearch
  could touch it — pure cost, no benefit. The `.json` URL is recorded in the
  module (`DEFAULT_JSON_URL`) only for completeness.
- **Refactor evaluation (asked): minimal and not needed now.** `dbsearch.py`
  is already the clean Qt-free matcher; the only related cleanup is that the
  atlas read in `main.py:enumerate_inputs` (`pd.read_csv('npatlas.tsv', ...)`)
  is hardcoded to that filename/cwd — the updater writes to the same path, so
  no change required. If a second database is added later (HMDB etc., dev-plan
  item 5), factor the atlas load + column-name mapping into a small loader then.
- **Threading caveat:** the 33 MB download currently runs on the main thread
  behind a wait cursor. It's a user-confirmed, infrequent (>30-day-gated)
  action so blocking briefly is acceptable, but moving it onto a `QThread`
  worker (like `AnalysisWorker`) is the obvious future improvement — left out
  here because GUI threading can't be verified headlessly.

### MPACT self-update checker (`mpactupdate.py`, `tests/test_mpactupdate.py`)

On startup, queries the GitHub Releases API for the configured repo
(`robertsamples/mpact` by default — Robert's fork), compares the latest
published release tag against the running version (`__version__`, kept in
`mpactupdate.py`; **keep it in sync with `main.py`'s `label_credits`** string,
currently `v1.00.01` -> `__version__ = '1.0.01'`), and if newer offers a
`git pull --ff-only` update (with a "please restart" prompt on success, or
opens the release page on failure). Version compare uses `packaging.version`
(PEP 440, numeric — so 2.10 > 2.9) with a dotted-int fallback; an unparseable
tag is treated as "not newer" (never nags). Every failure mode (offline, no
releases yet/404, malformed JSON, no git) is non-fatal and silent.

- **Updater-framework evaluation (asked): no off-the-shelf framework.** The
  standard option, `pyupdater`, targets *frozen* PyInstaller/cx_Freeze apps
  and needs its own patch-server + signing setup — heavyweight for a tool run
  from a git clone. For a source checkout the meaningful update is `git pull`,
  and "is there a newer release" is one API call + a version compare, which is
  all this module is. **Action needed from you:** tag releases on the fork
  (e.g. `v1.0.1`) and bump `__version__` per release, or this finds nothing.
- For the *portable PyInstaller build* (no git), `apply_git_update` will fail
  gracefully and the user is sent to the release page to download manually —
  a real auto-updater for the frozen build is the `pyupdater`-shaped project
  to consider only if/when that distribution channel matters.

### Crash / error reporter (`crashreport.py`, `tests/test_crashreport.py`)

Installs a `sys.excepthook` (after `QApplication` exists) that, on any
unhandled exception: chains to the default hook (traceback still hits the
console), formats a full report (traceback + MPACT/Python/platform versions +
timestamp + optional context), writes it to a timestamped file under
`~/.mpact/crashlogs/`, and shows a dialog offering to open a **prefilled
GitHub issue** (title + fenced traceback body) in the browser. Nothing is sent
without the user clicking through. The excepthook is hardened to never raise.

- **Crash-logger-framework evaluation (asked): Sentry is the off-the-shelf
  option, deliberately not used.** `sentry-sdk` is built for hosted/web
  services: it sends events to a Sentry project by default (silent cloud
  egress — wrong default for a desktop research tool), needs a DSN/account
  provisioned, and *still* needs a custom `before_send` hook + dialog to honour
  "ask the user first." The local-log + prefilled-GitHub-issue flow gives the
  maintainer the same thing (a complete traceback) with zero infrastructure
  and no privacy surprise. If MPACT later ships to many non-technical users and
  a central error feed becomes worthwhile, Sentry with `before_send` gating is
  the documented upgrade path (noted in `crashreport.py`).
- **PyQt5 note to verify live:** PyQt5 routes unhandled exceptions raised
  inside Qt slots through `sys.excepthook` (then may abort), so this should
  catch most in-GUI crashes — but the exact abort-after-hook behaviour is
  PyQt5-version-dependent and is the one thing to confirm by actually
  triggering an error in the running app.

### Dialog styling (`dialogs.py`, `tests/test_dialogs.py`)

The three subsystems above all pop `QMessageBox` dialogs. On the live app these
first rendered **black-on-black** (an unstyled dark background with invisible
black text — confirmed from a user screenshot): a `QMessageBox` inherits the
app's dark look but ships no text/background colours of its own.
`dialogs.styled_message_box()` applies a stylesheet matching the GUI palette
(background `rgb(40,40,40)`, text `rgb(212,212,212)`, detailed-text `QTextEdit`
darkened too) so every app dialog is legible and on-theme. Kept in its own
module (not `main.py`) so the box construction is headless-testable via
offscreen Qt (`build_message_box` returns the box without the blocking
`exec_`); `main.py`'s atlas/update/crash prompts all route through it.

Two follow-ups after the first attempt (from a second user screenshot):
- **Buttons stayed black-on-black/borderless.** The `QMessageBox QPushButton`
  *descendant* selector did not take effect on the standard buttons even
  though the box/label rules did. Fixed by styling each button object
  directly (`for b in box.buttons(): b.setStyleSheet(...)`) with a clearly
  visible border (`rgb(120,120,120)`) — selector-independent and reliable.
- **Native title bar was light + rounded** (Win11). `apply_dark_titlebar()`
  sets the DWM window attributes (immersive dark mode `20`/`19`, corner
  preference `33` = do-not-round) via `ctypes`/`dwmapi`, best-effort and
  Windows-only (no-op elsewhere, all failures swallowed). Called from
  `styled_message_box` after `winId()` realises the handle but before
  `exec_()` (dark mode must be set pre-show). **Verify live on Win11** — this
  is the part that can't be checked headlessly.

## Performance pass (2026-06-30, measurement-driven)

Profiled `run_MSFaST` on the bundled example dataset (cProfile + wall timing;
scratch scripts not committed) and benchmarked the algorithmic sections that
scale with feature/DB size. **Every change below was verified output-identical
against the original on real data, not just "looks equivalent"** — the bar the
user set ("functionally identical in terms of I/O").

Finding: on the small example the *pipeline* is dominated by pandas CSV
I/O (the inter-stage `iondict.csv`/`_formatted.csv` round-trips, ~0.6s of
to_csv + ~0.5s of read_csv out of ~2.3s), not by Python loops. That I/O chain
is the already-logged "bigger, multi-session" refactor (threading an `iondict`
DataFrame through `filter`/`stats`); left alone here as too invasive/risky for
this pass. The wins below are in the per-feature/per-DB-row algorithmic code,
which is what actually scales badly on large real datasets.

- **`dbsearch.search_npatlas`: ~5x faster, output identical.** Was
  O(features x atlas_rows): per feature it scanned all ~36k atlas rows twice
  (once per adduct) with a full-DataFrame boolean mask, then `.copy()` +
  `pd.concat` + `sort_values` + a scalar `.loc` write. Now pre-sorts the two
  adduct-mass columns once and uses `np.searchsorted` to test only a tiny m/z
  window per feature; the **exact original ppm test is re-applied to the
  windowed candidates** so the matched set is bit-identical (the window
  `mass*(1 ± 2·ppm/1e6)` is a proven superset of the true ppm window). Also:
  build one DataFrame per feature from concatenated m+h/m+na positions (no
  per-feature `pd.concat`), iterate numpy arrays instead of `iterrows`, and
  assign the `hits` column once instead of 979 scalar `.loc` sets. Verified on
  the real example (979 feats × 36,454 atlas rows): 1.41s → 0.28s, **0 hitdb
  DataFrame mismatches** (incl. row order + `ppm` values) and an identical
  `iondict['hits']` column. New edge-case tests in `test_dbsearch.py`
  (ppm-sort across both adducts; a single atlas row matching both adducts
  appearing twice).
- **`qualityscore.compute_cv_quality`: ~6.5x faster, output identical.** The
  AUC-under-the-CV-curve step was a per-feature Python loop doing
  `iondict.iloc[pos, :]['col']` scalar lookups (the classic slow pandas
  pattern) over thousands of rows. Replaced with the vectorised equivalent
  `np.sum(np.diff(cv, prepend=0) * pct)`. ~0.4s → ~0.06s per call (n≈5000).
  The faithfulness test (`test_qualityscore.py`, which pins against a verbatim
  copy of the original loop) confirms identical values; np.sum's pairwise
  summation can differ from the sequential loop by <1 ULP, far below the
  0.1%-rounded display precision.
- **`stats.groupave`: dead sum-of-squares chain removed** (pass 1) — dropped a
  `(chunk**2).groupby().sum()` per CSV chunk that only fed an unused variance.
- **`filter.relationalfilter`: measured, left alone.** Looks O(n²) but the
  early `break` once past the max isotope window makes it O(n·k) with small k:
  benchmarked at 0.017 / 0.077 / 0.371 s for 2k / 8k / 20k synthetic features
  (near-linear). Not a bottleneck; its intricate ringing/dimer-band logic
  isn't worth the regression risk to micro-optimize.
- **`filter.decon` / `stats.groupave` remaining cost is the per-stage CSV
  round-trips**, i.e. the same I/O-chain refactor noted above — not addressed
  here.
