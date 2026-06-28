"""
MPACT
Copyright 2022, Robert M. Samples, Sara P. Puckett, and Marcy J. Balunas

Per-plot widget/render state, consolidated.

MainWindow used to track one piece of state per plot name across six
separate dicts (fig/canvas/pltlayout/toolbar/ax/highlight), kept in sync
purely by convention -- nothing enforced that all six existed for a given
plot name. PlotSlotRegistry owns a single dict of PlotSlot objects (one per
plot name) instead. Existing call sites that index by name and field (e.g.
``self.canvas['heatmap']`` or ``parent.fig[currplt] = Figure()``) keep
working completely unchanged, since each field is exposed as a dict-like
view over the same underlying slots.
"""

from collections.abc import MutableMapping
from dataclasses import dataclass

FIELDS = ('fig', 'canvas', 'pltlayout', 'toolbar', 'ax', 'highlight')


@dataclass
class PlotSlot:
    fig: object = None
    canvas: object = None
    pltlayout: object = None
    toolbar: object = None
    ax: object = None
    highlight: object = None


class _PlotFieldView(MutableMapping):
    """Dict-like view of one PlotSlot field across every plot name."""

    def __init__(self, slots, field_name):
        self._slots = slots
        self._field = field_name

    def __getitem__(self, key):
        return getattr(self._slots[key], self._field)

    def __setitem__(self, key, value):
        if key not in self._slots:
            self._slots[key] = PlotSlot()
        setattr(self._slots[key], self._field, value)

    def __delitem__(self, key):
        setattr(self._slots[key], self._field, None)

    def __iter__(self):
        return iter(self._slots)

    def __len__(self):
        return len(self._slots)


class PlotSlotRegistry:
    """Owns one PlotSlot per plot name; exposes a dict-like view per field."""

    def __init__(self):
        self._slots = {}
        for field_name in FIELDS:
            setattr(self, field_name, _PlotFieldView(self._slots, field_name))

    def __contains__(self, key):
        return key in self._slots

    def __iter__(self):
        return iter(self._slots)

    def __len__(self):
        return len(self._slots)

    def __delitem__(self, key):
        del self._slots[key]

    def get_slot(self, key):
        """Return the PlotSlot for ``key``, or None if it doesn't exist."""
        return self._slots.get(key)
