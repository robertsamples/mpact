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
    assert slot.ax is None


def test_missing_key_raises_keyerror_like_a_dict():
    reg = PlotSlotRegistry()
    with pytest.raises(KeyError):
        reg.fig['missing']


def test_iteration_yields_plot_names():
    reg = PlotSlotRegistry()
    reg.fig['heatmap'] = 'FIG'
    reg.canvas['volcano'] = 'CANVAS'
    assert set(reg.fig) == {'heatmap', 'volcano'}
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
    assert reg.fig['heatmap'] is None
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
    assert len(reg.fig) == 2
