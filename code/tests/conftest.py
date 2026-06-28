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
