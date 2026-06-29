"""Pytest configuration for MPACT unit tests.

Adds the ``code/`` directory (the parent of this ``tests/`` folder) to
``sys.path`` so the application modules can be imported as top-level modules
(e.g. ``import filter``) exactly as they are at runtime.
"""

import os
import sys

CODE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if CODE_DIR not in sys.path:
    sys.path.insert(0, CODE_DIR)

# PyQt5's "offscreen" platform plugin lets QApplication/QWidget subclasses be
# constructed and exercised with no real display -- this does NOT make
# main.py importable headlessly (the main<->ui_functions circular import
# documented in devnotes.md is unrelated and still applies), but it does mean
# standalone Qt-coupled modules that don't import main.py (e.g. searchtree.py)
# can have their actual widget/model/signal behaviour covered by real tests
# instead of relying solely on manual testing. Set before any test imports
# PyQt5, and only if the environment hasn't already chosen a platform.
os.environ.setdefault('QT_QPA_PLATFORM', 'offscreen')

import pytest


@pytest.fixture(scope='session')
def qapp():
    """One shared QApplication for the whole test session -- constructing
    more than one in the same process crashes, so any test touching real
    Qt widgets/models should depend on this fixture rather than creating
    its own."""
    from PyQt5.QtWidgets import QApplication
    app = QApplication.instance()
    if app is None:
        app = QApplication([])
    return app
