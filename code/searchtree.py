"""
MPACT
Copyright 2022, Robert M. Samples, Sara P. Puckett, and Marcy J. Balunas

Replaces the feature-search tab's Designer-created QTreeWidget
(``ui_main.py``'s ``self.treeWidget``) with a real QTreeView backed by a
QAbstractTableModel + QSortFilterProxyModel, plus a small per-column filter
bar -- so each column gets its own search/range/checkbox filter instead of
one global search box.

This does NOT edit any generated UI file. ``SearchTreePanel`` is
constructed once, at runtime, against the already-built Designer widget
(see ``MainWindow.__init__``): it removes the Designer QTreeWidget from its
layout and inserts a container (filter bar + QTreeView) in the same slot --
the same "swap a Designer placeholder for hand-built content at runtime"
pattern ``plotting.py`` already uses for matplotlib canvases inserted into
Designer-created QFrame placeholders. Re-opening the .ui file in Qt
Designer would still show the old QTreeWidget; only the running app uses
this.

The actual filter-matching rules live in ``treefilters.py`` (Qt-free,
unit-tested) -- this module is just the Qt plumbing that calls into them.
"""

from PyQt5 import QtCore, QtGui, QtWidgets

from treefilters import CategoryFilter, RangeFilter, TextFilter, distinct_category_tokens, row_passes

COLUMNS = ('Compound', 'm/z', 'TR', 'Max', 'Sets', 'Groups', 'FC', 'Hits')
TEXT_COLUMNS = {0}
NUMERIC_COLUMNS = {1, 2, 3, 6, 7}
CATEGORY_COLUMNS = {4, 5}

_NUMERIC_FORMATS = {
    1: '{:.4f}',  # m/z
    2: '{:.3f}',  # TR (retention time)
    3: '{:.0f}',  # Max
    6: '{:.2f}',  # FC
    7: '{:.0f}',  # Hits
}

# Matches the dark palette ui_main.py already sets on the Designer treeWidget
# (rgb(70,70,70)-family backgrounds, rgb(200,200,200) light-grey text) --
# these widgets are newly created at runtime, so they don't pick that up
# automatically the way the treeWidget itself did from setupUi().
_FILTER_BAR_STYLE = """
QWidget {
    background-color: rgba(70,70,70,25);
}
QLineEdit, QToolButton {
    background-color: rgb(50,50,50);
    color: rgb(200,200,200);
    border: 1px solid rgb(70,70,70);
    border-radius: 2px;
    padding: 2px;
}
QLabel {
    color: rgb(200,200,200);
    background: transparent;
}
QToolButton::menu-indicator {
    image: none;
}
QMenu {
    background-color: rgb(50,50,50);
    color: rgb(200,200,200);
    border: 1px solid rgb(70,70,70);
}
QMenu::item:selected {
    background-color: rgb(70,70,70);
}
"""


class IonTableModel(QtCore.QAbstractTableModel):
    """One row per feature; columns match COLUMNS. Stores raw (unformatted,
    correctly typed) values so sorting and filtering compare real numbers,
    not the display strings -- ``RAW_ROLE`` exposes them for that purpose.
    """

    RAW_ROLE = QtCore.Qt.UserRole + 1

    def __init__(self, parent=None):
        super().__init__(parent)
        self._rows = []

    def set_rows(self, rows):
        """``rows``: list of tuples, one per feature, in COLUMNS order."""
        self.beginResetModel()
        self._rows = list(rows)
        self.endResetModel()

    def raw_row(self, row):
        return self._rows[row]

    def rowCount(self, parent=QtCore.QModelIndex()):
        return 0 if parent.isValid() else len(self._rows)

    def columnCount(self, parent=QtCore.QModelIndex()):
        return 0 if parent.isValid() else len(COLUMNS)

    def data(self, index, role=QtCore.Qt.DisplayRole):
        if not index.isValid():
            return None
        value = self._rows[index.row()][index.column()]
        if role == self.RAW_ROLE:
            return value
        if role == QtCore.Qt.DisplayRole:
            fmt = _NUMERIC_FORMATS.get(index.column())
            if fmt is not None:
                try:
                    return fmt.format(float(value))
                except (TypeError, ValueError):
                    return str(value)
            return str(value)
        return None

    def headerData(self, section, orientation, role=QtCore.Qt.DisplayRole):
        if role == QtCore.Qt.DisplayRole and orientation == QtCore.Qt.Horizontal:
            return COLUMNS[section]
        return None


class IonFilterProxyModel(QtCore.QSortFilterProxyModel):
    """Per-column filters (column index -> a TextFilter/RangeFilter/
    CategoryFilter from treefilters.py) combined with AND, evaluated via
    treefilters.row_passes(). Sorting compares IonTableModel.RAW_ROLE
    values instead of the default display-text comparison, so numeric
    columns sort numerically without any string-parsing tricks.
    """

    def __init__(self, parent=None):
        super().__init__(parent)
        self._filters = {}

    def set_filter(self, column, criterion):
        if criterion is None:
            self._filters.pop(column, None)
        else:
            self._filters[column] = criterion
        self.invalidateFilter()

    def filterAcceptsRow(self, source_row, source_parent):
        model = self.sourceModel()
        if model is None:
            return True
        return row_passes(model.raw_row(source_row), self._filters)

    def lessThan(self, left, right):
        model = self.sourceModel()
        left_value = model.data(left, IonTableModel.RAW_ROLE)
        right_value = model.data(right, IonTableModel.RAW_ROLE)
        try:
            return left_value < right_value
        except TypeError:
            return str(left_value) < str(right_value)


class SearchTreePanel:
    """Owns the model/proxy/view/filter-bar trio that replaces a Designer
    QTreeWidget. Construct once with the Designer widget to replace; use
    ``set_rows()``/``selected_compound()`` afterward in place of the old
    ``treeWidget.clear()``/``addTopLevelItem()``/``selectedItems()`` calls.
    """

    def __init__(self, old_tree_widget):
        layout = old_tree_widget.parentWidget().layout()
        index = layout.indexOf(old_tree_widget)
        parent = old_tree_widget.parentWidget()

        # Reuse the Designer-set dark-theme stylesheet rather than guessing
        # colours -- the new QTreeView is a genuinely different widget
        # instance, so it doesn't inherit what setupUi() applied to the one
        # being replaced.
        old_stylesheet = old_tree_widget.styleSheet()

        layout.removeWidget(old_tree_widget)
        old_tree_widget.setParent(None)
        old_tree_widget.deleteLater()

        self.model = IonTableModel(parent)
        self.proxy = IonFilterProxyModel(parent)
        self.proxy.setSourceModel(self.model)

        self.view = QtWidgets.QTreeView(parent)
        self.view.setModel(self.proxy)
        self.view.setSortingEnabled(True)
        self.view.setRootIsDecorated(False)
        # Not the original QTreeWidget's behaviour -- it never set this, and
        # the app's underlying QPalette was never themed dark for the
        # AlternateBase role (only the QSS stylesheets are dark), so Qt's
        # default (light grey) alternate-row colour clashes badly here.
        self.view.setAlternatingRowColors(False)
        self.view.setSelectionBehavior(QtWidgets.QAbstractItemView.SelectRows)
        self.view.setUniformRowHeights(True)
        if old_stylesheet:
            self.view.setStyleSheet(old_stylesheet)

        container = QtWidgets.QWidget(parent)
        outer = QtWidgets.QVBoxLayout(container)
        outer.setContentsMargins(0, 0, 0, 0)
        outer.setSpacing(2)

        # A single-line toggle rather than always showing all 8 filters at
        # once -- one text box, five min/max range pairs, and two category
        # dropdowns doesn't fit comfortably side by side in a narrow side
        # panel. Collapsed by default so it doesn't change anything for
        # someone who never needs to filter.
        toggle_row = QtWidgets.QWidget(container)
        toggle_row.setStyleSheet(_FILTER_BAR_STYLE)
        toggle_layout = QtWidgets.QHBoxLayout(toggle_row)
        toggle_layout.setContentsMargins(4, 2, 4, 2)
        self.toggle_button = QtWidgets.QToolButton(toggle_row)
        self.toggle_button.setCheckable(True)
        self.toggle_button.setText('Filters ▸')
        self.toggle_button.toggled.connect(self._on_toggle_filters)
        toggle_layout.addWidget(self.toggle_button)
        toggle_layout.addStretch()

        self.filter_panel = QtWidgets.QWidget(container)
        self.filter_panel.setStyleSheet(_FILTER_BAR_STYLE)
        self._filter_layout = QtWidgets.QFormLayout(self.filter_panel)
        self._filter_layout.setContentsMargins(6, 4, 6, 4)
        self._filter_layout.setSpacing(4)
        self.filter_panel.setVisible(False)

        self._category_buttons = {}
        for column in range(len(COLUMNS)):
            self._filter_layout.addRow(COLUMNS[column], self._build_filter_widget(column))

        outer.addWidget(toggle_row)
        outer.addWidget(self.filter_panel)
        outer.addWidget(self.view)

        layout.insertWidget(index, container)

    def _on_toggle_filters(self, checked):
        self.filter_panel.setVisible(checked)
        self.toggle_button.setText('Filters ▾' if checked else 'Filters ▸')

    def _build_filter_widget(self, column):
        if column in TEXT_COLUMNS:
            return self._build_text_filter(column)
        if column in NUMERIC_COLUMNS:
            return self._build_range_filter(column)
        return self._build_category_filter(column)

    def _build_text_filter(self, column):
        edit = QtWidgets.QLineEdit()
        edit.setPlaceholderText('contains…')
        edit.setClearButtonEnabled(True)
        edit.textChanged.connect(
            lambda text, c=column: self.proxy.set_filter(c, TextFilter(text) if text else None))
        return edit

    def _build_range_filter(self, column):
        box = QtWidgets.QWidget()
        layout = QtWidgets.QHBoxLayout(box)
        layout.setContentsMargins(0, 0, 0, 0)
        layout.setSpacing(4)

        min_edit = QtWidgets.QLineEdit()
        min_edit.setPlaceholderText('min')
        min_edit.setValidator(QtGui.QDoubleValidator())
        max_edit = QtWidgets.QLineEdit()
        max_edit.setPlaceholderText('max')
        max_edit.setValidator(QtGui.QDoubleValidator())

        def update_filter(*_args):
            minimum = float(min_edit.text()) if min_edit.text() else None
            maximum = float(max_edit.text()) if max_edit.text() else None
            criterion = RangeFilter(minimum, maximum) if (minimum is not None or maximum is not None) else None
            self.proxy.set_filter(column, criterion)

        min_edit.textChanged.connect(update_filter)
        max_edit.textChanged.connect(update_filter)

        layout.addWidget(min_edit)
        layout.addWidget(max_edit)
        return box

    def _build_category_filter(self, column):
        button = QtWidgets.QToolButton()
        button.setText('Select…')
        button.setPopupMode(QtWidgets.QToolButton.InstantPopup)
        menu = QtWidgets.QMenu(button)
        button.setMenu(menu)
        self._category_buttons[column] = (button, menu, {})
        return button

    def _refresh_category_options(self, column):
        button, menu, _old_actions = self._category_buttons[column]
        menu.clear()
        actions = {}
        values = [self.model.raw_row(r)[column] for r in range(self.model.rowCount())]
        for token in distinct_category_tokens(values):
            action = menu.addAction(token)
            action.setCheckable(True)
            action.toggled.connect(lambda _checked, c=column: self._update_category_filter(c))
            actions[token] = action
        self._category_buttons[column] = (button, menu, actions)

    def _update_category_filter(self, column):
        _button, _menu, actions = self._category_buttons[column]
        selected = frozenset(token for token, action in actions.items() if action.isChecked())
        self.proxy.set_filter(column, CategoryFilter(selected) if selected else None)

    def set_rows(self, rows):
        """``rows``: list of tuples, one per feature, in COLUMNS order
        (Compound, m/z, TR, Max, Sets, Groups, FC, Hits)."""
        self.model.set_rows(rows)
        for column in CATEGORY_COLUMNS:
            self._refresh_category_options(column)

    def selected_compound(self):
        """Return the Compound id of the currently selected row, or None."""
        indexes = self.view.selectionModel().selectedRows()
        if not indexes:
            return None
        source_index = self.proxy.mapToSource(indexes[0])
        return self.model.raw_row(source_index.row())[0]
