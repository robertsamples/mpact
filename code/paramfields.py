"""
MPACT
Copyright 2022, Robert M. Samples, Sara P. Puckett, and Marcy J. Balunas

Shared schema for the ``analysis_parameters`` fields that are plain 1:1
checkbox reads/writes -- no branching, no derived values, no side effects.
``MainWindow.enumerate_inputs`` (save, ``main.py``) and ``loadsession``
(restore, ``ui_functions.py``) both iterate ``CHECKBOX_FIELDS`` via the
helpers below instead of each hand-writing a mirrored line per field, so
adding one of this kind of field means adding one entry here instead of
editing both functions.

Deliberately narrow in scope: fields with conditional logic, derived
values, backward-compatibility defaults, or side effects (graphfilters
construction, blank-filter abs/rel branching, max-iso-shift combo-by-text
mapping, groupset construction, NPAtlas filtering, output-directory
creation, ...) are intentionally NOT here and stay hand-written in
enumerate_inputs/loadsession, where that logic is easier to follow inline
than abstracted into a generic table.
"""

# (analysis_parameters attribute name, widget path as ('container', 'name'))
CHECKBOX_FIELDS = (
    ('PCA', ('ui', 'checkBox_pca')),
    ('Dendrogram', ('ui', 'checkBox_dend')),
    ('MZRTplt', ('ui', 'checkBox_mzrt')),
    ('KMD', ('ui', 'checkBox_kmd')),
    # 'mdguide' moved to the mass-defect plot's own live control bar
    # (plotting.kendrick); no longer a dialog checkbox or a persisted field.
    ('FC', ('ui', 'checkBox_fc')),
    ('FC3Dplt', ('ui', 'checkBox_3dfc')),
    ('Ttest', ('ui', 'checkBox_ttest')),
    ('Volcanoplt', ('ui', 'checkBox_volcano')),
    ('FDR', ('ui', 'checkBox_FDR')),
)


def _resolve_widget(window, path):
    """Resolve a dotted container path (e.g. ``'dialog.ui'``) plus widget
    name against ``window`` (the ``MainWindow`` instance)."""
    container, name = path
    obj = window
    for part in container.split('.'):
        obj = getattr(obj, part)
    return getattr(obj, name)


def save_checkbox_fields(window, params):
    """Copy every ``CHECKBOX_FIELDS`` widget's checked state into ``params``."""
    for attr, path in CHECKBOX_FIELDS:
        setattr(params, attr, _resolve_widget(window, path).isChecked())


def restore_checkbox_fields(window, params, restore):
    """Restore every ``CHECKBOX_FIELDS`` widget's checked state from ``params``.

    ``restore(what, action)`` is the caller's per-field failure-isolation
    helper (see ``loadsession``) -- each field restores independently so a
    missing attribute on an older save can't abort the rest.
    """
    for attr, path in CHECKBOX_FIELDS:
        restore(attr, lambda attr=attr, path=path:
                _resolve_widget(window, path).setChecked(getattr(params, attr)))
