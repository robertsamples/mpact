from PyQt5 import QtCore, QtWidgets

from searchtree import COLUMNS, IonFilterProxyModel, IonTableModel, SearchTreePanel
from treefilters import CategoryFilter, RangeFilter, TextFilter

SAMPLE_ROWS = [
    ('c1', 101.5, 1.0, 1000, ' Blanks 0um_Ce', ' Blanks 0um_Ce', 1.5, 2),
    ('c2', 250.25, 2.0, 2000, ' 250um_Ce', ' 250um_Ce', 3.0, 0),
    ('c3', 99.9, 0.5, 500, ' Media', ' Media', 0.5, 5),
]


def make_model_and_proxy():
    model = IonTableModel()
    model.set_rows(SAMPLE_ROWS)
    proxy = IonFilterProxyModel()
    proxy.setSourceModel(model)
    return model, proxy


# --------------------------------------------------------------------------- #
# IonTableModel
# --------------------------------------------------------------------------- #

def test_model_row_and_column_counts(qapp):
    model, _ = make_model_and_proxy()
    assert model.rowCount() == 3
    assert model.columnCount() == len(COLUMNS)


def test_model_header_labels(qapp):
    model, _ = make_model_and_proxy()
    for i, name in enumerate(COLUMNS):
        assert model.headerData(i, 1) == name  # Qt.Horizontal == 1


def test_model_raw_role_returns_unformatted_value(qapp):
    from PyQt5.QtCore import QModelIndex
    model, _ = make_model_and_proxy()
    index = model.index(0, 1)  # c1's m/z
    assert model.data(index, IonTableModel.RAW_ROLE) == 101.5


def test_model_display_role_formats_numeric_columns(qapp):
    model, _ = make_model_and_proxy()
    mz_index = model.index(0, 1)
    assert model.data(mz_index) == '101.5000'
    hits_index = model.index(0, 7)
    assert model.data(hits_index) == '2'


# --------------------------------------------------------------------------- #
# IonFilterProxyModel
# --------------------------------------------------------------------------- #

def visible_compounds(proxy):
    return [proxy.data(proxy.index(r, 0)) for r in range(proxy.rowCount())]


def test_no_filters_shows_every_row(qapp):
    _, proxy = make_model_and_proxy()
    assert set(visible_compounds(proxy)) == {'c1', 'c2', 'c3'}


def test_text_filter_on_compound_column(qapp):
    _, proxy = make_model_and_proxy()
    proxy.set_filter(0, TextFilter('c2'))
    assert visible_compounds(proxy) == ['c2']


def test_range_filter_on_numeric_column(qapp):
    _, proxy = make_model_and_proxy()
    proxy.set_filter(7, RangeFilter(minimum=1))  # Hits >= 1
    assert set(visible_compounds(proxy)) == {'c1', 'c3'}


def test_category_filter_on_groups_column(qapp):
    _, proxy = make_model_and_proxy()
    proxy.set_filter(4, CategoryFilter(frozenset({'Media'})))
    assert visible_compounds(proxy) == ['c3']


def test_combined_filters_across_columns(qapp):
    _, proxy = make_model_and_proxy()
    proxy.set_filter(6, RangeFilter(minimum=1, maximum=2))   # FC in [1,2]
    proxy.set_filter(4, CategoryFilter(frozenset({'0um_Ce'})))
    assert visible_compounds(proxy) == ['c1']


def test_clearing_a_filter_restores_rows(qapp):
    _, proxy = make_model_and_proxy()
    proxy.set_filter(0, TextFilter('c2'))
    assert len(visible_compounds(proxy)) == 1
    proxy.set_filter(0, None)
    assert len(visible_compounds(proxy)) == 3


def test_sorting_numeric_column_is_numeric_not_lexicographic(qapp):
    _, proxy = make_model_and_proxy()
    proxy.sort(1, 0)  # sort by m/z ascending (Qt.AscendingOrder == 0)
    mz_values = [proxy.data(proxy.index(r, 1), IonTableModel.RAW_ROLE) for r in range(proxy.rowCount())]
    assert mz_values == sorted(mz_values)


# --------------------------------------------------------------------------- #
# SearchTreePanel -- the actual Designer-widget-swap, end to end
# --------------------------------------------------------------------------- #

def make_host_with_tree_widget():
    """Mimic the relevant slice of ui_main.py's setupUi(): a frame containing
    a layout with a QTreeWidget in it -- the exact shape SearchTreePanel
    expects to replace."""
    host = QtWidgets.QWidget()
    layout = QtWidgets.QVBoxLayout(host)
    tree_widget = QtWidgets.QTreeWidget(host)
    layout.addWidget(tree_widget)
    return host, layout, tree_widget


def test_panel_removes_old_widget_and_inserts_view_in_same_slot(qapp):
    host, layout, tree_widget = make_host_with_tree_widget()
    assert layout.indexOf(tree_widget) == 0

    panel = SearchTreePanel(tree_widget)

    assert layout.indexOf(tree_widget) == -1  # removed
    assert layout.count() == 1  # the new container took its place
    assert panel.view.parentWidget() is not None


def test_panel_copies_the_old_widgets_stylesheet_onto_the_new_view(qapp):
    """The new QTreeView is a different widget instance, so it doesn't
    automatically inherit whatever Designer's setupUi() applied to the
    widget it's replacing -- confirm it's carried across explicitly rather
    than defaulting to unstyled (light-themed) widgets in a dark-themed app."""
    host, layout, tree_widget = make_host_with_tree_widget()
    tree_widget.setStyleSheet('QTreeView { background-color: rgb(50,50,50); }')

    panel = SearchTreePanel(tree_widget)

    assert panel.view.styleSheet() == 'QTreeView { background-color: rgb(50,50,50); }'


def test_panel_set_rows_populates_the_view(qapp):
    host, layout, tree_widget = make_host_with_tree_widget()
    panel = SearchTreePanel(tree_widget)
    panel.set_rows(SAMPLE_ROWS)
    assert panel.proxy.rowCount() == 3


def test_panel_selected_compound_maps_through_proxy_to_source(qapp):
    host, layout, tree_widget = make_host_with_tree_widget()
    panel = SearchTreePanel(tree_widget)
    panel.set_rows(SAMPLE_ROWS)

    # Filter down to just c3, select the only visible row, and confirm
    # selected_compound() correctly maps the (filtered, possibly reordered)
    # proxy row back to the right source row rather than assuming row 0
    # of the proxy lines up with row 0 of the unfiltered data.
    panel.proxy.set_filter(0, TextFilter('c3'))
    proxy_index = panel.proxy.index(0, 0)
    panel.view.selectionModel().select(
        proxy_index,
        QtCore.QItemSelectionModel.Select | QtCore.QItemSelectionModel.Rows,
    )
    assert panel.selected_compound() == 'c3'


def test_panel_selected_compound_none_when_nothing_selected(qapp):
    host, layout, tree_widget = make_host_with_tree_widget()
    panel = SearchTreePanel(tree_widget)
    panel.set_rows(SAMPLE_ROWS)
    assert panel.selected_compound() is None


def test_panel_category_filter_options_repopulate_on_set_rows(qapp):
    host, layout, tree_widget = make_host_with_tree_widget()
    panel = SearchTreePanel(tree_widget)
    panel.set_rows(SAMPLE_ROWS)

    _button, _menu, actions = panel._category_buttons[4]
    assert set(actions.keys()) == {'Blanks', '0um_Ce', '250um_Ce', 'Media'}
