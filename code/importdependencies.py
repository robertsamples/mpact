"""
MPACT
Copyright 2022, Robert M. Samples, Sara P. Puckett, and Marcy J. Balunas

Dependency bootstrap.

A few of MPACT's dependencies are not part of a stock Anaconda install
(``epam.indigo``, ``UpSetPlot`` and ``squarify``).  These are imported at the
top level of ``main.py`` / ``plotting.py``, so if any of them is missing the
program crashes on import before it can do anything useful.  This module
detects the missing packages, installs them with pip, and -- because a freshly
pip-installed package usually cannot be imported into an already-running
interpreter/kernel (e.g. when launched through Spyder on Anaconda) -- tells the
user to restart when that happens.

``ensure_dependencies()`` is called once, as early as possible, from
``main.py`` (before the packages it installs are imported).
"""

import importlib
import subprocess
import sys

# pip package name -> importable module name (required: gate startup if missing)
DEPENDENCIES = {
    "epam.indigo": "indigo",
    "UpSetPlot": "upsetplot",
    "squarify": "squarify",
}

# Performance-only packages: installed best-effort, but a failure here must
# never block startup (the code falls back to slower implementations).
OPTIONAL_DEPENDENCIES = {
    "fastcluster": "fastcluster",  # faster hierarchical clustering (heatmap, bootstrap dendrogram)
}


def _is_available(import_name):
    """Return True if ``import_name`` can be imported in this interpreter."""
    try:
        return importlib.util.find_spec(import_name) is not None
    except (ImportError, ValueError):
        # find_spec can raise for half-installed / namespace edge cases.
        return False


def install(package):
    """Install a single package with the current interpreter's pip.

    IMPORTANT: this runs inside a conda (Anaconda) environment whose pandas,
    matplotlib, scipy, etc. are compiled against NumPy 1.x. A bare
    ``pip install`` can upgrade NumPy to 2.x as a side effect, which makes every
    one of those conda-compiled modules fail to import ("module compiled using
    NumPy 1.x cannot be run in NumPy 2.5.0"). We therefore constrain NumPy to
    <2 and only upgrade dependencies when strictly necessary, so installing an
    optional package can never break the environment.
    """
    subprocess.check_call([
        sys.executable, "-m", "pip", "install",
        "--upgrade-strategy", "only-if-needed",
        package, "numpy<2",
    ])


def checkdep():
    """Install any missing dependencies.

    Returns a list of the pip package names that were freshly installed (these
    may require an interpreter/kernel restart before they can be imported).
    """
    print("Checking dependencies")
    installed = []
    for pip_name, import_name in DEPENDENCIES.items():
        if _is_available(import_name):
            print(import_name + " already installed")
            continue
        print("Installing " + pip_name)
        try:
            install(pip_name)
            installed.append(pip_name)
        except Exception as exc:  # pragma: no cover - network/pip dependent
            print("Failed to install " + pip_name + ": " + str(exc))

    # Optional/perf packages: try to install, but never treat a failure as fatal
    # and never report them as gating a restart.
    for pip_name, import_name in OPTIONAL_DEPENDENCIES.items():
        if _is_available(import_name):
            continue
        print("Installing optional " + pip_name)
        try:
            install(pip_name)
        except Exception as exc:  # pragma: no cover - network/pip dependent
            print("Optional dependency " + pip_name + " not installed: " + str(exc))

    # Allow newly installed packages to be discovered without a restart where
    # the interpreter supports it.
    importlib.invalidate_caches()
    return installed


def _notify_restart(message):
    """Best-effort, GUI-agnostic "please restart" notice."""
    print(message)
    try:
        from tkinter import Tk, messagebox

        root = Tk()
        root.withdraw()  # don't show an empty root window
        messagebox.showinfo("Restart required", message)
        root.destroy()
    except Exception:
        # No display / tkinter unavailable -- the printed message is enough.
        pass


def ensure_dependencies():
    """Make sure every dependency is importable, installing as needed.

    If a package had to be installed but still isn't importable in this
    process (the typical Spyder/Anaconda first-run case), the user is asked to
    restart and the process exits cleanly instead of crashing with an
    ImportError.
    """
    installed = checkdep()

    still_missing = [name for name in DEPENDENCIES.values() if not _is_available(name)]
    if still_missing:
        _notify_restart(
            "Dependencies were installed but require a restart.\n"
            "Please restart MPACT (or restart the kernel in Spyder)."
        )
        sys.exit(0)
