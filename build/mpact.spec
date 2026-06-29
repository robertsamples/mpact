# -*- mode: python ; coding: utf-8 -*-
#
# PyInstaller spec for a portable MPACT build (no Anaconda / conda env needed
# on the end user's machine). See build/README.md for how to run this and
# devnotes.md for the background on why --onedir (not --onefile) was chosen.
#
# Build from the repo root with:
#   pyinstaller build/mpact.spec --noconfirm
#
# Output lands in dist/MPACT/ (a folder, not a single exe) -- ship the whole
# folder. dist/MPACT/MPACT.exe is the launcher.

import os

from PyInstaller.utils.hooks import collect_all

block_cipher = None

REPO_ROOT = os.path.abspath(os.path.join(os.path.dirname(SPEC), os.pardir))
CODE_DIR = os.path.join(REPO_ROOT, "code")

# Data files main.py/plotting.py read via cwd-relative paths (run.bat / this
# spec both assume the process cwd is the folder main.py itself lives in --
# here, the dist/MPACT/ folder next to the exe). npatlas.tsv is read-only
# reference data; compoundimages/ is a read+write cache of rendered structure
# PNGs that main.py creates new entries in at runtime (code/main.py ~line
# 664), so it must land in a writable location next to the exe, not inside
# the PyInstaller onefile temp-extraction dir -- this is the main reason this
# build uses --onedir instead of --onefile. cog.ico is the plot-toolbar
# "Configure" button's icon (plotting.py's NavigationToolbar -- opens the
# bootstrap/collapse-technical-replicates dialog); it's a bare cwd-relative
# QIcon("cog.ico") load there too, same as the other two, and was missed in
# the original spec -- PyInstaller's static analysis can't see a file path
# inside a string literal, so anything loaded that way must be added here
# explicitly or it silently fails to load (a QToolButton with no icon and no
# text collapses to looking like it's missing entirely, rather than erroring).
datas = [
    (os.path.join(CODE_DIR, "npatlas.tsv"), "."),
    (os.path.join(CODE_DIR, "compoundimages"), "compoundimages"),
    (os.path.join(CODE_DIR, "cog.ico"), "."),
]
binaries = []

# Packages whose PyInstaller hooks don't always pick up every native/data
# file automatically. epam.indigo ships its own compiled cheminformatics
# engine (indigo/lib/<platform>/*.dll, loaded by a hand-rolled path lookup in
# indigo's own _common/lib.py -- not detected by static import analysis at
# all, so it must be pulled in explicitly via collect_all); scikit-learn/scipy
# ship compiled extension modules with import-time plugin discovery that
# static analysis can miss.
hiddenimports = [
    "sklearn.utils._cython_blas",
    "sklearn.neighbors._quad_tree",
    "sklearn.tree._utils",
    "scipy.special.cython_special",
]
for pkg in ("indigo",):
    pkg_datas, pkg_binaries, pkg_hiddenimports = collect_all(pkg)
    datas += pkg_datas
    binaries += pkg_binaries
    hiddenimports += pkg_hiddenimports

a = Analysis(
    [os.path.join(CODE_DIR, "main.py")],
    pathex=[CODE_DIR],
    binaries=binaries,
    datas=datas,
    hiddenimports=hiddenimports,
    hookspath=[],
    hooksconfig={},
    runtime_hooks=[],
    excludes=[],
    noarchive=False,
    cipher=block_cipher,
)

pyz = PYZ(a.pure, a.zipped_data, cipher=block_cipher)

exe = EXE(
    pyz,
    a.scripts,
    [],
    exclude_binaries=True,
    name="MPACT",
    debug=False,
    bootloader_ignore_signals=False,
    strip=False,
    upx=False,
    # Keep a console window, matching run.bat's existing behavior: MPACT
    # prints non-fatal errors/progress to stdout and there's no in-app log
    # viewer yet (see the main.py TODO about a status-bar terminal view), so
    # console=False would silently swallow diagnostics a user might need to
    # report a bug.
    console=True,
    icon=os.path.join(CODE_DIR, "cog.ico") if os.path.isfile(os.path.join(CODE_DIR, "cog.ico")) else None,
    # PyInstaller >=6.0 defaults onedir builds to nesting everything except
    # the exe under an _internal/ subfolder. main.py's data paths
    # ('npatlas.tsv', 'compoundimages/...') are bare and resolved against the
    # process cwd, which is dist/MPACT/ (the exe's own folder) when launched
    # normally -- so keep the flat pre-6.0 layout, putting datas/binaries
    # directly alongside MPACT.exe instead of in _internal/.
    contents_directory=".",
)

coll = COLLECT(
    exe,
    a.binaries,
    a.zipfiles,
    a.datas,
    strip=False,
    upx=False,
    name="MPACT",
)
