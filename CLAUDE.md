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

Headless unit tests in `code/tests/` (pure-logic only — no GUI). Run:

```
python -m pytest code/tests -q
```

Covers `filter`, `stats`, `importdependencies`, `translators`, `groupsets`. Add
tests here for any new Qt-free logic. GUI behaviour must be verified by running
the app (can't be tested headlessly).

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
mid-session. 57 passing tests.
