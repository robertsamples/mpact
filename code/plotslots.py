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

Each field view's keys are exactly the plot names that field has actually
been assigned for -- matching how the six dicts behaved independently
before (a plot could have a canvas but no highlight marker yet, and
``'name' in self.highlight`` would be False even though
``'name' in self.canvas`` was True). A bare ``None`` default can't express
that distinction, since a field that genuinely was never set and one that's
deliberately holding ``None`` would be indistinguishable -- so unset fields
use the ``_UNSET`` sentinel instead, and every view treats a field holding
it as "this key doesn't exist," consistent with plain dict semantics.
"""

from collections.abc import MutableMapping
from dataclasses import dataclass

FIELDS = ('fig', 'canvas', 'pltlayout', 'toolbar', 'ax', 'highlight')

_UNSET = object()


@dataclass
class PlotSlot:
    fig: object = _UNSET
    canvas: object = _UNSET
    pltlayout: object = _UNSET
    toolbar: object = _UNSET
    ax: object = _UNSET
    highlight: object = _UNSET


class _PlotFieldView(MutableMapping):
    """Dict-like view of one PlotSlot field across every plot name."""

    def __init__(self, slots, field_name):
        self._slots = slots
        self._field = field_name

    def __getitem__(self, key):
        value = getattr(self._slots[key], self._field)
        if value is _UNSET:
            raise KeyError(key)
        return value

    def __setitem__(self, key, value):
        if key not in self._slots:
            self._slots[key] = PlotSlot()
        setattr(self._slots[key], self._field, value)

    def __delitem__(self, key):
        if key not in self._slots or getattr(self._slots[key], self._field) is _UNSET:
            raise KeyError(key)
        setattr(self._slots[key], self._field, _UNSET)

    def __iter__(self):
        return (key for key, slot in self._slots.items()
                if getattr(slot, self._field) is not _UNSET)

    def __len__(self):
        return sum(1 for _ in self)


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
