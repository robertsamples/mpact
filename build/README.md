# Portable Windows build (no Anaconda required)

Produces a folder you can zip and hand to someone with no Python, no
Anaconda, and no `pip install` of their own -- they just unzip and double-click
`MPACT.exe`.

## Build it

```
build\build_windows.bat
```

This creates a throwaway venv at `build\.buildvenv\` (never touches your dev
Anaconda env), installs `requirements.txt` + PyInstaller into it, and runs
`pyinstaller build\mpact.spec`. Output: `dist\MPACT\` -- a folder, **not** a
single exe. Ship the whole folder; `dist\MPACT\MPACT.exe` is the launcher.

Re-running the script reuses the existing venv. Delete `build\.buildvenv\` to
force a clean dependency reinstall; delete `dist\` and `build\mpact\` to force
a clean PyInstaller rebuild (the batch script does the build step
unconditionally, but PyInstaller itself caches via its own `build\mpact\*.toc`
files).

## Why onedir, not onefile, and why `--contents-directory .`

- **`--onedir` (the PyInstaller default this spec uses), not `--onefile`.**
  `code/main.py` reads/writes `npatlas.tsv` and `compoundimages/*.png` via
  bare relative paths (e.g. `'compoundimages/' + cmpd + '.png'`), resolved
  against the process's current working directory -- not against the
  installed package location. A onefile build extracts everything into a
  fresh temp directory on *every launch* and deletes it on exit, so anything
  main.py wrote into `compoundimages/` (newly-rendered structure PNGs) would
  vanish the moment the app closed, and onefile startup is also slower
  (extract-then-run vs. run-in-place). Onedir keeps those files persistently
  next to the exe, matching how `run.bat` already works (it `cd`s into
  `code\` before launching, for the same cwd-relative-paths reason -- see
  CLAUDE.md).
- **`contents_directory='.'` on the `EXE(...)` call in `mpact.spec`.**
  PyInstaller 6.0 changed the onedir default to nest everything except the
  exe itself under `dist/MPACT/_internal/`, to keep the top-level folder
  clean. That breaks the same cwd-relative-path assumption above (cwd is
  `dist/MPACT/`, but `npatlas.tsv` would land in `dist/MPACT/_internal/`).
  Setting `contents_directory='.'` restores the pre-6.0 flat layout instead
  of patching main.py's path handling to be build-system-aware.

## Dependency notes specific to this build

- `requirements.txt` deliberately does **not** carry the `numpy<2` cap that
  `importdependencies.py` / CLAUDE.md document for the Anaconda dev env. That
  cap exists only to protect conda-compiled pandas/matplotlib/scipy binaries
  built against NumPy 1.x's ABI -- it doesn't apply to a clean pip venv, where
  every wheel is built against NumPy 2.x from the start. Capping it here
  would actually make pip's resolve fail outright: the PyPI `fastcluster`
  wheel for Python 3.12 requires `numpy>=2`. (The codebase has no NumPy
  1.x-only API usage -- `np.float`/`np.int`/etc. were checked and aren't
  used -- so NumPy 2.x is safe here.)
- `epam.indigo` ships a compiled cheminformatics engine
  (`indigo/lib/<platform>/*.dll`) that it loads via its own hand-rolled path
  lookup in `indigo/_common/lib.py`, completely invisible to PyInstaller's
  static import analysis -- a plain `hiddenimports=['indigo']` is not
  enough; the DLL itself never gets copied and the app fails at
  `Indigo()` construction with `RuntimeError: Could not find native
  libraries`. Fixed in `mpact.spec` via
  `PyInstaller.utils.hooks.collect_all('indigo')`, which pulls in the
  package's data files (the DLLs included) and binaries, not just its
  Python modules.
- scikit-learn and scipy ship compiled extension submodules that their
  PyInstaller hooks don't always discover via static analysis (import-time
  C-level plugin lookups). A few are listed explicitly in
  `hiddenimports` in `mpact.spec`; if a future scikit-learn/scipy upgrade
  introduces an `ImportError`/`ModuleNotFoundError` only in the *frozen* exe
  (not when running `python main.py` normally), that's almost always this
  category of problem -- check `build\mpact\warn-mpact.txt` after a build
  for "not found" hidden-import warnings as a starting point.

## Smoke-testing a build

The frozen exe can be sanity-checked headlessly (no display needed) the same
way the test suite's `conftest.py` does it, via Qt's offscreen platform:

```
set QT_QPA_PLATFORM=offscreen
dist\MPACT\MPACT.exe
```

It won't get far without real peak-table data, but a clean launch with no
`Traceback`/`ImportError` in the console confirms every dependency (PyQt5,
pandas, indigo's native DLL, etc.) is actually bundled and importable in the
frozen environment -- which is the failure mode this kind of build most
commonly hits, not application logic bugs.

## Linux / Mac builds

Not yet built or tested on those platforms. See the "Portable builds for
Linux / Mac" section in the repo-root `CLAUDE.md` for what's expected to
carry over from this Windows build vs. what would need platform-specific
work, for whoever picks this up next.
