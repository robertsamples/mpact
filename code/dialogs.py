"""
MPACT
Copyright 2026, Robert M. Samples

Dark-themed QMessageBox helpers, matching the main GUI palette so the app's
dialogs don't render with invisible black-on-black text (the default when a
QMessageBox inherits the app's dark styling but has no colours of its own).

Kept out of ``main.py`` (which can't be imported headlessly -- the documented
main<->ui_functions circular import) so the box construction/styling can be
unit-tested via Qt's offscreen platform (see ``tests/test_dialogs.py``), the
same approach used for ``searchtree.py``.

Palette mirrors ui_main.py: background rgb(40,40,40), text rgb(212,212,212),
buttons rgb(62,62,62) / hover rgb(75,75,75).

The push-buttons are styled *per widget* (``box.buttons()``) rather than via a
``QMessageBox QPushButton`` descendant selector: in practice that descendant
rule did not take effect on the standard buttons (they rendered borderless with
black text), while the box/label rules did -- styling each button object
directly is selector-independent and reliable. On Windows the native title bar
is also switched to dark + square corners via the DWM API so the dialog frame
matches the app's dark theme instead of the default light, rounded Win11 bar.
"""

import sys

from PyQt5 import QtWidgets

_BG = 'rgb(40,40,40)'
_TEXT = 'rgb(212,212,212)'

DIALOG_STYLE = """
QMessageBox { background-color: %s; }
QMessageBox QLabel { color: %s; }
QMessageBox QTextEdit { background-color: rgb(35,35,35); color: %s; }
""" % (_BG, _TEXT, _TEXT)

_BUTTON_STYLE = """
QPushButton {
    background-color: rgb(62,62,62);
    color: rgb(212,212,212);
    border: 1px solid rgb(120,120,120);
    border-radius: 3px;
    padding: 4px 16px;
    min-width: 64px;
}
QPushButton:hover { background-color: rgb(75,75,75); }
QPushButton:pressed { background-color: rgb(55,55,55); }
QPushButton:default { border: 1px solid rgb(160,160,160); }
"""


def apply_dark_titlebar(widget):
    """Best-effort: make a top-level window's title bar dark with square
    corners on Windows 11 (no-op / silently ignored everywhere else).

    Uses the DWM window attributes (immersive dark mode + corner preference).
    Must run before the window is first shown to take effect cleanly, so call
    it after ``winId()`` realises the native handle but before ``exec_()``.
    """
    if sys.platform != 'win32':
        return
    try:
        import ctypes
        hwnd = int(widget.winId())
        dwm = ctypes.windll.dwmapi
        flag = ctypes.c_int(1)
        # DWMWA_USE_IMMERSIVE_DARK_MODE: 20 on current Win10/11, 19 on early
        # 20H1 builds -- set both; the wrong one just returns a failure code.
        for attr in (20, 19):
            dwm.DwmSetWindowAttribute(hwnd, attr, ctypes.byref(flag), ctypes.sizeof(flag))
        # DWMWA_WINDOW_CORNER_PREFERENCE = 33, DWMWCP_DONOTROUND = 1 (Win11).
        corner = ctypes.c_int(1)
        dwm.DwmSetWindowAttribute(hwnd, 33, ctypes.byref(corner), ctypes.sizeof(corner))
    except Exception:
        pass


def build_message_box(parent, icon, title, text, buttons=None, default=None,
                      detailed=None):
    """Construct a QMessageBox styled to match the MPACT GUI (does NOT exec).

    Separated from :func:`styled_message_box` so tests can inspect the
    configured box without blocking on a modal ``exec_()``.
    """
    box = QtWidgets.QMessageBox(parent)
    box.setIcon(icon)
    box.setWindowTitle(title)
    box.setText(text)
    if buttons is not None:
        box.setStandardButtons(buttons)
    if default is not None:
        box.setDefaultButton(default)
    if detailed is not None:
        box.setDetailedText(detailed)
    box.setStyleSheet(DIALOG_STYLE)
    # Style each button object directly -- reliable where the descendant
    # selector wasn't applied to the standard buttons.
    for button in box.buttons():
        button.setStyleSheet(_BUTTON_STYLE)
    return box


def styled_message_box(parent, icon, title, text, buttons=None, default=None,
                       detailed=None):
    """Build the styled box, show it modally, and return the clicked button."""
    box = build_message_box(parent, icon, title, text, buttons=buttons,
                            default=default, detailed=detailed)
    apply_dark_titlebar(box)
    return box.exec_()
