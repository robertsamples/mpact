"""Headless tests for the dark-themed dialog helpers (``dialogs.py``).

Uses the offscreen Qt platform (the ``qapp`` fixture in conftest.py) to build
the message boxes without a display and without blocking on a modal exec.
Guards the regression that prompted this module: a QMessageBox with no
explicit colours rendered black-on-black under the app's dark styling.
"""

import sys

from PyQt5 import QtWidgets

import dialogs


def test_style_sets_visible_label_and_dark_background(qapp):
    # The style must define a light label colour and a non-default background
    # (the fix for the invisible black-on-black text).
    assert 'color: rgb(212,212,212)' in dialogs.DIALOG_STYLE
    assert 'background-color: rgb(40,40,40)' in dialogs.DIALOG_STYLE


def test_build_message_box_applies_style_and_content(qapp):
    box = dialogs.build_message_box(
        None, QtWidgets.QMessageBox.Question, 'Title here', 'Body text',
        buttons=QtWidgets.QMessageBox.Yes | QtWidgets.QMessageBox.No,
        default=QtWidgets.QMessageBox.No)
    assert isinstance(box, QtWidgets.QMessageBox)
    assert box.text() == 'Body text'
    if sys.platform != 'darwin':
        # Qt's Cocoa (macOS) integration treats QMessageBox as a native
        # alert panel -- per Apple HIG, alerts have no title bar -- and
        # doesn't retain the windowTitle property for it specifically (other
        # widget types aren't affected). build_message_box still calls
        # setWindowTitle unconditionally since it's meaningful on every other
        # platform (and harmless here); only the readback assertion is
        # platform-gated. Reproduced on stock PyQt5 with no styling applied
        # at all, so this is not something dialogs.py's theming can fix.
        assert box.windowTitle() == 'Title here'
    # The dark theme is actually applied to this box.
    assert 'rgb(212,212,212)' in box.styleSheet()
    assert box.standardButtons() == (QtWidgets.QMessageBox.Yes | QtWidgets.QMessageBox.No)
    assert box.defaultButton() == box.button(QtWidgets.QMessageBox.No)
    box.deleteLater()


def test_each_button_is_styled_directly(qapp):
    # Regression: the descendant selector didn't reach the standard buttons,
    # so each button must carry the button stylesheet itself (visible text +
    # border).
    box = dialogs.build_message_box(
        None, QtWidgets.QMessageBox.Question, 't', 'b',
        buttons=QtWidgets.QMessageBox.Yes | QtWidgets.QMessageBox.No)
    buttons = box.buttons()
    assert len(buttons) == 2
    for button in buttons:
        sheet = button.styleSheet()
        assert 'color: rgb(212,212,212)' in sheet
        assert 'border: 1px solid' in sheet
    box.deleteLater()


def test_apply_dark_titlebar_never_raises(qapp):
    # No-op off Windows; on Windows it best-effort sets DWM attributes and must
    # swallow any failure (e.g. an offscreen/invalid HWND).
    box = dialogs.build_message_box(None, QtWidgets.QMessageBox.Information, 't', 'b')
    dialogs.apply_dark_titlebar(box)  # must not raise
    box.deleteLater()


def test_build_message_box_supports_detailed_text(qapp):
    box = dialogs.build_message_box(
        None, QtWidgets.QMessageBox.Critical, 'Crash', 'Something failed',
        buttons=QtWidgets.QMessageBox.Yes | QtWidgets.QMessageBox.No,
        detailed='full traceback here')
    assert box.detailedText() == 'full traceback here'
    # The detailed-text pane is a QTextEdit, also styled for visibility.
    assert 'QTextEdit' in box.styleSheet()
    box.deleteLater()
