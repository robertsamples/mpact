import pytest

from plotslots import PlotSlotRegistry


def test_set_and_get_roundtrip():
    reg = PlotSlotRegistry()
    reg.fig['heatmap'] = 'FIG'
    reg.canvas['heatmap'] = 'CANVAS'
    assert reg.fig['heatmap'] == 'FIG'
    assert reg.canvas['heatmap'] == 'CANVAS'


def test_independent_fields_share_one_slot():
    reg = PlotSlotRegistry()
    reg.fig['heatmap'] = 'FIG'
    reg.canvas['heatmap'] = 'CANVAS'
    slot = reg.get_slot('heatmap')
    assert slot.fig == 'FIG'
    assert slot.canvas == 'CANVAS'
    assert 'heatmap' not in reg.ax


def test_missing_key_raises_keyerror_like_a_dict():
    reg = PlotSlotRegistry()
    with pytest.raises(KeyError):
        reg.fig['missing']


def test_field_view_only_yields_keys_that_field_was_set_for():
    """Regression guard: a plot with fig/canvas set but no highlight marker
    must not show up when iterating reg.highlight, even though the same
    plot name is present in reg.fig/reg.canvas. This previously broke
    _refresh_highlight()'s `for plot in self.highlight: ...` loop -- it
    iterated every plot name that had ANY field set, not just ones that
    actually had a highlight artist, and crashed on the resulting None.
    """
    reg = PlotSlotRegistry()
    reg.fig['heatmap'] = 'FIG'
    reg.canvas['heatmap'] = 'CANVAS'
    reg.highlight['volcano'] = 'MARKER'

    assert 'heatmap' not in reg.highlight
    assert set(reg.highlight) == {'volcano'}
    with pytest.raises(KeyError):
        reg.highlight['heatmap']

    # but the registry as a whole knows about both plot names
    assert set(reg) == {'heatmap', 'volcano'}


def test_iteration_yields_only_keys_set_on_that_field():
    reg = PlotSlotRegistry()
    reg.fig['heatmap'] = 'FIG'
    reg.canvas['volcano'] = 'CANVAS'
    assert set(reg.fig) == {'heatmap'}
    assert set(reg.canvas) == {'volcano'}
    assert set(reg) == {'heatmap', 'volcano'}


def test_contains():
    reg = PlotSlotRegistry()
    reg.fig['heatmap'] = 'FIG'
    assert 'heatmap' in reg
    assert 'heatmap' in reg.fig
    assert 'volcano' not in reg


def test_delitem_on_field_view_clears_only_that_field():
    reg = PlotSlotRegistry()
    reg.fig['heatmap'] = 'FIG'
    reg.canvas['heatmap'] = 'CANVAS'
    del reg.fig['heatmap']
    assert 'heatmap' not in reg.fig
    with pytest.raises(KeyError):
        reg.fig['heatmap']
    assert reg.canvas['heatmap'] == 'CANVAS'


def test_delitem_on_registry_removes_whole_slot():
    reg = PlotSlotRegistry()
    reg.fig['heatmap'] = 'FIG'
    del reg['heatmap']
    assert 'heatmap' not in reg
    with pytest.raises(KeyError):
        reg.fig['heatmap']


def test_len():
    reg = PlotSlotRegistry()
    reg.fig['heatmap'] = 'FIG'
    reg.canvas['volcano'] = 'CANVAS'
    assert len(reg) == 2
    assert len(reg.fig) == 1
    assert len(reg.canvas) == 1
