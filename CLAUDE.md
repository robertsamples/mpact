# CLAUDE.md — MPACT

Guidance for Claude Code (and humans) working in this repo. MPACT is a PyQt5
desktop tool for mass-spectrometry / natural-products data analysis. Code lives
in `code/`.

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

## ⚠️ NumPy 2.x / Anaconda dependency hazard (read before touching deps)

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
  core). Each corresponding `MainWindow` method is now a thin wrapper:
  call the module function, then apply the result to widgets/`self`.
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
`searchtree`. Add tests here for any new Qt-free logic.

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
