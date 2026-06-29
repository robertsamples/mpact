"""
MPACT
Copyright 2022, Robert M. Samples, Sara P. Puckett, and Marcy J. Balunas

Model for the user-defined "groupsets" used to colour/group features in plots
by their presence/absence across biological groups.

This used to be a bare list of ``query`` objects (``main.py``) plus a
``selset`` integer index mutated directly by UI callbacks (``ui_functions.py``)
with no bounds checking. ``GroupSetModel`` owns that data and every mutation,
so the UI layer only has to keep widgets in sync with it (a thin
view/controller), rather than implement the bookkeeping itself.

This module is Qt-free and unit-tested (see ``tests/test_groupsets.py``).
"""


class GroupSet:
    """One groupset: a colour-coded presence/absence rule for plotting features.

    Attributes:
        name: display name shown in the groupset list.
        src: biological groups a feature must be present in at least one of
            ("or" groups). Empty means no constraint from this rule.
        incl: biological groups a feature must be present in ALL of
            ("and" groups).
        excl: biological groups a feature must NOT be present in.
        colour: hex colour string used to plot features matching this groupset.
    """

    __slots__ = ('name', 'src', 'incl', 'excl', 'colour')

    def __init__(self, name='New Feature Set', src=None, incl=None, excl=None, colour='#000000'):
        self.name = name
        self.src = list(src) if src else []
        self.incl = list(incl) if incl else []
        self.excl = list(excl) if excl else []
        self.colour = colour

    @classmethod
    def from_legacy(cls, obj):
        """Build from an older pickled ``query``-shaped object (or anything
        exposing the same four attributes), for ``.mpct`` backward compatibility.
        """
        return cls(
            name=getattr(obj, 'name', 'New Feature Set'),
            src=getattr(obj, 'src', []),
            incl=getattr(obj, 'incl', []),
            excl=getattr(obj, 'excl', []),
            colour=getattr(obj, 'colour', '#000000'),
        )

    def __repr__(self):
        return ('GroupSet(name=%r, src=%r, incl=%r, excl=%r, colour=%r)'
                % (self.name, self.src, self.incl, self.excl, self.colour))

    def __eq__(self, other):
        if not isinstance(other, GroupSet):
            return NotImplemented
        return (self.name, self.src, self.incl, self.excl, self.colour) == \
               (other.name, other.src, other.incl, other.excl, other.colour)


class GroupSetModel:
    """Owns the ordered collection of GroupSets and which one is selected.

    All mutation goes through this class (add/remove/update/select) so the
    selected index can never point out of bounds -- the bug class that
    afflicted the old bare ``self.selset`` integer (e.g. after removing the
    last item, or before any groupset existed).
    """

    def __init__(self):
        self._items = []
        self._selected = -1

    def __len__(self):
        return len(self._items)

    def __iter__(self):
        return iter(self._items)

    def __getitem__(self, index):
        return self._items[index]

    @property
    def selected_index(self):
        return self._selected

    @property
    def selected(self):
        if 0 <= self._selected < len(self._items):
            return self._items[self._selected]
        return None

    def select(self, index):
        """Set the selected index, clamped to a valid range (-1 if empty).

        Returns the resulting (clamped) index.
        """
        if not self._items:
            self._selected = -1
        else:
            self._selected = max(0, min(index, len(self._items) - 1))
        return self._selected

    def add(self, name='New Feature Set', all_groups=None):
        """Create a new groupset that excludes every known biological group by
        default (matching the previous ``addgroup`` behaviour), select it, and
        return it.
        """
        groupset = GroupSet(name=name, excl=list(all_groups or []))
        self._items.append(groupset)
        self._selected = len(self._items) - 1
        return groupset

    def remove(self, index=None):
        """Remove a groupset (the selected one by default)."""
        if index is None:
            index = self._selected
        if 0 <= index < len(self._items):
            del self._items[index]
        self.select(self._selected)

    def move(self, from_index, to_index):
        """Reorder the groupset at ``from_index`` to ``to_index``.

        Both indices are clamped to the valid range; out-of-range or equal
        indices are a no-op. Selection follows the moved item, so the
        groupset that was selected before the move is still selected after
        (by identity, not by index) -- a drag-and-drop reorder shouldn't
        change which groupset is being edited.
        """
        if not self._items:
            return
        from_index = max(0, min(from_index, len(self._items) - 1))
        to_index = max(0, min(to_index, len(self._items) - 1))
        if from_index == to_index:
            return
        selected_item = self.selected
        groupset = self._items.pop(from_index)
        self._items.insert(to_index, groupset)
        if selected_item is not None:
            # Identity, not '==' -- GroupSet.__eq__ is value-based, and two
            # distinct groupsets can compare equal (e.g. freshly added ones
            # before either is edited), so list.index() could pick the wrong
            # one.
            for i, item in enumerate(self._items):
                if item is selected_item:
                    self._selected = i
                    break

    def update(self, index, *, name=None, src=None, incl=None, excl=None, colour=None):
        """Overwrite the given fields of the groupset at ``index``."""
        groupset = self._items[index]
        if name is not None:
            groupset.name = name
        if src is not None:
            groupset.src = list(src)
        if incl is not None:
            groupset.incl = list(incl)
        if excl is not None:
            groupset.excl = list(excl)
        if colour is not None:
            groupset.colour = colour

    def to_legacy_list(self):
        """Plain list of ``GroupSet`` objects, for ``.mpct`` pickling."""
        return list(self._items)

    @classmethod
    def from_legacy_list(cls, items):
        """Build a model from a list of ``query``-shaped objects (older
        ``.mpct`` saves) or ``GroupSet`` objects (current saves). The first
        item becomes selected, matching the previous load behaviour.
        """
        model = cls()
        for item in items or []:
            model._items.append(item if isinstance(item, GroupSet) else GroupSet.from_legacy(item))
        if model._items:
            model._selected = 0
        return model


def normalize_graphfilters(graphfilters):
    """Coerce ``graphfilters`` to a list of filter-name tokens.

    Current code populates ``analysis_parameters.graphfilters`` as a list
    (``['cv', 'rel', 'insource']``); older ``.mpct`` saves pickled it as a
    space-joined string instead, so a loaded session can hand back either
    shape here.
    """
    if isinstance(graphfilters, str):
        return graphfilters.split()
    return list(graphfilters)


def build_query_dict(model, graphfilters=''):
    """Build the ``{descriptive_name: GroupSet}`` mapping MSFaST/plotting consume.

    Each groupset's exclusion list additionally gets the currently active
    graph-level filters (cv/rel/insource, from ``graphfilters``) appended, so
    features removed by those filters stay excluded from every groupset.

    The descriptive name (used as the dict key, and written verbatim into
    ``analysisinfo.txt`` as a human-readable summary of the groupset's rule)
    has the form ``"<src> +<incl> -<excl> c=<colour> n=<name>"``.
    """
    extra_excl = normalize_graphfilters(graphfilters)
    result = {}
    for groupset in model:
        merged_excl = groupset.excl + extra_excl
        descriptive_name = (
            ' '.join(groupset.src) + ' +' + ' '.join(groupset.incl) +
            ' -' + ' '.join(merged_excl) +
            ' c=' + str(groupset.colour) + ' n=' + groupset.name
        )
        result[descriptive_name] = GroupSet(
            name=groupset.name, src=groupset.src, incl=groupset.incl,
            excl=merged_excl, colour=groupset.colour,
        )
    return result
