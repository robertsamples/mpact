# Development

This page is for people contributing to MPACT itself, not end users. For
the authoritative, most up-to-date version of these notes (kept alongside
the code), see [`devnotes.md`](https://github.com/BalunasLab/mpact/blob/main/devnotes.md)
in the repo root.

## Architecture

- **Generated — do not edit:** `ui_main.py`, `ui_main1.py`,
  `ui_featureinfo.py`, `ui_plotparam.py`, `files.py`, `files_rc.py`. These
  are Qt Designer output and get overwritten on regeneration. (Despite the
  name, `ui_functions.py` is hand-written and fully editable — it's the
  `UIFunctions` controller class.)
- **Hand-written app code:** `main.py` (`MainWindow`, run/save/load,
  database search), `plotting.py` (plot classes), `filter.py`, `stats.py`,
  `MSFaST.py` (analysis driver), `pvclust.py` (bootstrap dendrogram),
  `ordination.py` (Qt-free PCA/NMDS/PLS-DA backend),
  `clusterpurity.py` (dendrogram branch-purity logic),
  `csvcache.py` (cached CSV reads for the ordination data path),
  `translators.py` (import/export framework), `mzmineimport.py` (format
  conversion), `getfragdb.py`, `mspwriter.py`.
- **Canonical peak table** format (what `MSFaST` consumes internally;
  Progenesis is the native/baseline format): CSV with 3 header rows, row 2
  = `Compound,m/z,Retention time (min),<injections...>`, col0 = `RT_mz` id,
  col1 = m/z, col2 = RT.

## Threading model

`run_MSFaST` is Qt-free and runs on a `QThread` worker (`AnalysisWorker` in
`main.py`), so the GUI stays responsive during the heavy compute.
`MainWindow.run_analysis` reads widgets on the main thread, starts the
worker, and `_finish_analysis` does all matplotlib/Qt plotting back on the
**main thread** (matplotlib is not thread-safe). Never create Qt/matplotlib
objects on the worker thread.

## Importer/translator framework (`translators.py`)

Qt-free and unit-tested: `detect_peaktable_format`, `parse_msp`/
`parse_mgf` (→ `FragmentEntry`), `reindex_fragments` (matches fragments to
peak-table rows by compound ID first, then m/z+RT — Progenesis MSP stores
neutral mass, not adduct m/z), `filter_source_peaktable` (row-subsets the
source peak table to surviving features). `mzmineimport.format_check`
delegates detection to this module.

## Groupsets MVC (`groupsets.py`)

`GroupSet` (data) + `GroupSetModel` (collection, bounds-safe selection,
CRUD) + `build_query_dict()` replace what used to be a bare list +
selected-index pair. `MainWindow.groupsetmodel` is the live state;
`ui_functions.py`'s `addgroup`/`removegroup`/`updatesets`/`updategroups`/
`writegroups`/`colour_picker1` are thin view-sync controllers over it.
`main.py`'s `query` class still exists, but **only** as the unpickle target
for old `.mpct` files — `GroupSet.from_legacy`/`GroupSetModel.from_legacy_list`
convert on load.

## Testing

Headless unit tests live in `code/tests/` (pure-logic only — no Qt):

```
python -m pytest code/tests -q
```

Covers `filter`, `stats`, `importdependencies`, `translators`,
`groupsets`, and `ordination`. Add tests here for any new Qt-free logic.
GUI behaviour can't be tested headlessly — verify it by running the app.

## Conventions

- Never edit the generated UI files listed above.
- Plot generation goes through `MainWindow.safe_generate`, so one failing
  plot doesn't abort the rest.
- `.mpct` saves are atomic (temp file + `os.replace`), with per-component
  guards (`write_save`).
- `loadsession` restores each saved parameter independently — a bad/missing
  field can't cascade and abort restoration of the rest. Add new analysis
  parameters to **both** `enumerate_inputs` (save) and `loadsession`
  (restore).
- Plot objects (`self.ftplt`, `self.kmd`, `self.spec`, ...) are created the
  first time they're needed and `.reset()` afterward, via
  `MainWindow._create_or_reset()` / `_generate_plots()` — never gate
  create-vs-reset on a whole-session flag like `self.analysisrun`, since an
  optional output can newly turn on mid-session for a dataset that didn't
  have it before, and the object would never get created.
- Use `MainWindow._refresh_highlight()` (not `highlight_feature()`) to
  redraw the current selection without changing it (e.g. on a tab switch).
  `highlight_feature(newfeature)` is for real selection events and toggles
  the highlight off if the same feature is clicked twice — calling it with
  the already-selected feature re-triggers that toggle, which is a bug, not
  a refresh.
- Every matplotlib `pick_event` handler must call
  `plotting._is_duplicate_pick(parent, event)` first and bail if it
  returns `True`. Matplotlib fires one `pick_event` per artist that
  registers a hit, not one per click — a feature plotted in more than one
  groupset/colour layer otherwise fires the handler twice per click.
- `importdependencies.checkdep()` should stay silent when nothing needs
  installing — it runs on every launch, including every Spyder "Run File"
  (which re-executes `main.py`'s top level). Only report actual
  installs/failures.

## Building the docs site locally

```
pip install mkdocs mkdocs-material
mkdocs serve
```

Then open `http://127.0.0.1:8000`. See [Hosting](#hosting-this-site) below
for deployment.

## Hosting this site

This site is plain static HTML generated by MkDocs — there's no backend,
database, or server-side logic, so it needs essentially no infrastructure.
**GitHub Pages is the right fit here**: it's free, MkDocs has built-in
support for deploying to it (`mkdocs gh-deploy`), and a low-traffic docs
site for a research tool doesn't need a dedicated host. A `gh-pages`
deploy workflow is included in this repo (`.github/workflows/docs.yml`) —
once GitHub Pages is enabled for this repo (Settings → Pages → Source:
`gh-pages` branch), the site builds and publishes automatically on every
push to `main` that touches `docs/` or `mkdocs.yml`.
