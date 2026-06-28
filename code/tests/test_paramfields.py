import pytest

from paramfields import CHECKBOX_FIELDS, restore_checkbox_fields, save_checkbox_fields


class FakeCheckbox:
    def __init__(self, checked=False):
        self._checked = checked

    def isChecked(self):
        return self._checked

    def setChecked(self, value):
        self._checked = value


class FakeUi:
    pass


class FakeDialog:
    def __init__(self):
        self.ui = FakeUi()


class FakeWindow:
    def __init__(self):
        self.ui = FakeUi()
        self.dialog = FakeDialog()
        for attr, (container, name) in CHECKBOX_FIELDS:
            obj = self.ui if container == 'ui' else self.dialog.ui
            setattr(obj, name, FakeCheckbox())


class FakeParams:
    pass


def test_no_duplicate_attrs_or_widget_paths():
    attrs = [attr for attr, _ in CHECKBOX_FIELDS]
    paths = [path for _, path in CHECKBOX_FIELDS]
    assert len(attrs) == len(set(attrs))
    assert len(paths) == len(set(paths))


def test_save_then_restore_roundtrip():
    window = FakeWindow()
    for attr, (container, name) in CHECKBOX_FIELDS:
        obj = window.ui if container == 'ui' else window.dialog.ui
        getattr(obj, name).setChecked(True)

    params = FakeParams()
    save_checkbox_fields(window, params)
    for attr, _ in CHECKBOX_FIELDS:
        assert getattr(params, attr) is True

    # flip every widget, then restore from params and confirm they flip back
    for attr, (container, name) in CHECKBOX_FIELDS:
        obj = window.ui if container == 'ui' else window.dialog.ui
        getattr(obj, name).setChecked(False)

    def restore(what, action):
        action()

    restore_checkbox_fields(window, params, restore)
    for attr, (container, name) in CHECKBOX_FIELDS:
        obj = window.ui if container == 'ui' else window.dialog.ui
        assert getattr(obj, name).isChecked() is True


def test_restore_isolates_failures_per_field():
    window = FakeWindow()
    params = FakeParams()
    # Only set the first field's attribute -- the rest are "missing", as if
    # restoring from an older save file.
    setattr(params, CHECKBOX_FIELDS[0][0], True)

    calls = []

    def restore(what, action):
        try:
            action()
            calls.append((what, True))
        except Exception:
            calls.append((what, False))

    restore_checkbox_fields(window, params, restore)

    assert calls[0] == (CHECKBOX_FIELDS[0][0], True)
    assert all(ok is False for _, ok in calls[1:])
