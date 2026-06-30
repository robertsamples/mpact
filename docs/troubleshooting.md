# Troubleshooting

## NumPy 2.x environment break

**Symptom:** MPACT (or anything else in the same Anaconda environment that
imports pandas/matplotlib/scipy) fails to import with an error like:

```
ImportError: A module that was compiled using NumPy 1.x cannot be run in
NumPy 2.5.0
```

**Cause:** Anaconda's base environment ships pandas/matplotlib/scipy
compiled against NumPy 1.x. A plain `pip install <package>` can silently
upgrade NumPy to 2.x as a side effect, which breaks every one of those
conda-compiled modules at once.

**Fix:**

```
pip install "numpy<2"
```

MPACT's own dependency installer (`importdependencies.py`) always pins
`numpy<2` and uses `--upgrade-strategy only-if-needed` specifically so it
can't cause this — but a manual `pip install` of something unrelated, in
the same environment, still can. If you need to install something into
the MPACT Anaconda environment, pin `numpy<2` alongside it.

## "Checking dependencies" printed repeatedly in Spyder

This is expected and harmless. Spyder's User Module Reloader re-executes a
script's full top level — including MPACT's startup dependency check — on
every "Run File," not just the first time in a session. The dependency
check itself only prints output when it actually needs to install
something or a dependency install fails; if you're seeing the same
"already installed" message repeated, you're likely on an older version —
update to the latest `code/importdependencies.py`.

## A dependency was installed but MPACT still won't import it

This is the typical first-run case under Spyder/Anaconda: a freshly
`pip install`-ed package usually can't be imported into an already-running
interpreter/kernel. MPACT detects this and shows a "please restart" dialog
— restart MPACT (or restart the Spyder kernel) once and it will resolve.

## `run.bat` / desktop shortcut doesn't launch MPACT

Both `run.bat` and the `MPACT` shortcut resolve paths relative to the
repo's actual location, so moving or re-cloning the repo shouldn't require
editing either. If you previously edited either by hand (or copied the
repo from a different machine/account), check that:

- `run.bat` is being run from inside the repo (it uses `%~dp0`, its own
  location, to find `code/`).
- The shortcut's target and working directory both still point at this
  repo's `run.bat` / repo root — re-create the shortcut if in doubt.

`run.bat` also assumes Anaconda is installed at
`%USERPROFILE%\anaconda3`; edit the `activate.bat` line near the top if
yours is installed elsewhere.

## A plot tab won't generate / `AttributeError` on `.reset()`

If you're running an older build, this could be the "plot object never
created" bug: switching to a dataset that newly enables an optional output
(e.g. a fragment file, after a prior dataset in the same session had none)
could crash because the plot object was only ever created on the *first*
analysis run of the session. This is fixed in the current version — plot
objects are created lazily the first time they're actually needed,
independent of whether this is the first analysis run in the session.

## A feature won't select / immediately deselects when clicked

If a feature is plotted in more than one colour layer (e.g. it belongs to
overlapping [Plot Feature Sets](user-guide/analysis-settings.md#plot-feature-sets-groupsets)),
older builds could fire the click handler twice for a single click,
which combined with "click an already-selected feature to deselect it"
made the feature appear to never select. This is fixed in the current
version (`plotting._is_duplicate_pick`).

## Still stuck?

Check `code/tests/` — the pure-logic modules (`filter`, `stats`,
`translators`, `groupsets`, `importdependencies`) have headless unit
tests you can run to rule out a logic bug:

```
python -m pytest code/tests -q
```

GUI behaviour itself can't be tested headlessly and needs to be checked by
running the app — see
[`devnotes.md`](https://github.com/robertsamples/mpact/blob/main/devnotes.md)
in the repo root.
