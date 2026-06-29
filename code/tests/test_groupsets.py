"""Unit tests for the groupset model (``groupsets.py``).

These cover the MVC refactor of the previous ``self.querys`` list +
``self.selset`` index pair: selection bounds-safety, CRUD operations, legacy
``.mpct`` (pickled ``query`` object) compatibility, and the query-dict builder
consumed by ``MSFaST``/``plotting``.
"""

import pytest

from groupsets import GroupSet, GroupSetModel, build_query_dict, normalize_graphfilters


class LegacyQuery:
    """Stand-in for the old ``main.query`` class, as found in old .mpct files."""
    def __init__(self):
        self.name = ''
        self.src = ''
        self.incl = ''
        self.excl = ''
        self.colour = '#000000'


# --------------------------------------------------------------------------- #
# GroupSet
# --------------------------------------------------------------------------- #

def test_groupset_defaults():
    gs = GroupSet()
    assert gs.name == 'New Feature Set'
    assert gs.src == [] and gs.incl == [] and gs.excl == []
    assert gs.colour == '#000000'


def test_groupset_copies_input_lists_not_aliases():
    src = ['A', 'B']
    gs = GroupSet(src=src)
    gs.src.append('C')
    assert src == ['A', 'B']  # original untouched


def test_groupset_from_legacy_query_object():
    legacy = LegacyQuery()
    legacy.name = 'Old set'
    legacy.src = ['A']
    legacy.incl = ['B']
    legacy.excl = ['C']
    legacy.colour = '#ff0000'
    gs = GroupSet.from_legacy(legacy)
    assert gs == GroupSet(name='Old set', src=['A'], incl=['B'], excl=['C'], colour='#ff0000')


def test_groupset_from_legacy_missing_attrs_uses_defaults():
    class Bare:
        pass
    gs = GroupSet.from_legacy(Bare())
    assert gs == GroupSet()


# --------------------------------------------------------------------------- #
# GroupSetModel: selection safety
# --------------------------------------------------------------------------- #

def test_empty_model_selection_is_negative_one():
    model = GroupSetModel()
    assert model.selected_index == -1
    assert model.selected is None
    assert model.select(5) == -1  # still empty


def test_select_clamps_to_valid_range():
    model = GroupSetModel()
    model.add('a')
    model.add('b')
    model.add('c')
    assert model.select(99) == 2
    assert model.select(-5) == 0
    assert model.select(1) == 1
    assert model.selected.name == 'b'


def test_remove_last_item_reclamps_selection():
    model = GroupSetModel()
    model.add('a')
    model.add('b')
    model.select(1)
    model.remove()  # removes 'b' (selected)
    assert len(model) == 1
    assert model.selected_index == 0
    assert model.selected.name == 'a'


def test_remove_all_items_leaves_selection_at_negative_one():
    model = GroupSetModel()
    model.add('a')
    model.remove(0)
    assert len(model) == 0
    assert model.selected_index == -1
    assert model.selected is None


# --------------------------------------------------------------------------- #
# GroupSetModel: move (reordering)
# --------------------------------------------------------------------------- #

def test_move_reorders_items():
    model = GroupSetModel()
    model.add('a')
    model.add('b')
    model.add('c')
    model.move(0, 2)
    assert [g.name for g in model] == ['b', 'c', 'a']


def test_move_keeps_selection_on_the_moved_item():
    model = GroupSetModel()
    model.add('a')
    model.add('b')
    model.add('c')
    model.select(0)  # 'a' selected
    model.move(0, 2)
    assert model.selected.name == 'a'
    assert model.selected_index == 2


def test_move_keeps_selection_on_a_different_item_that_shifted_position():
    model = GroupSetModel()
    model.add('a')
    model.add('b')
    model.add('c')
    model.select(1)  # 'b' selected
    model.move(0, 2)  # moves 'a' past 'b', so 'b' shifts from index 1 to 0
    assert model.selected.name == 'b'
    assert model.selected_index == 0


def test_move_with_equal_indices_is_a_noop():
    model = GroupSetModel()
    model.add('a')
    model.add('b')
    model.move(1, 1)
    assert [g.name for g in model] == ['a', 'b']


def test_move_clamps_out_of_range_indices():
    model = GroupSetModel()
    model.add('a')
    model.add('b')
    model.add('c')
    model.move(-5, 99)
    assert [g.name for g in model] == ['b', 'c', 'a']


def test_move_on_empty_model_is_a_noop():
    model = GroupSetModel()
    model.move(0, 1)  # must not raise
    assert len(model) == 0


def test_move_disambiguates_value_equal_groupsets_by_identity():
    # Two freshly-added default groupsets compare equal (GroupSet.__eq__ is
    # value-based, and both start with identical fields), so move() must
    # track the selected item by identity, not by list.index()'s '=='.
    model = GroupSetModel()
    model.add('dup')
    model.add('dup')
    model.select(0)
    first = model.selected
    model.move(0, 1)
    assert model.selected is first
    assert model.selected_index == 1


# --------------------------------------------------------------------------- #
# GroupSetModel: CRUD
# --------------------------------------------------------------------------- #

def test_add_excludes_all_known_groups_by_default():
    model = GroupSetModel()
    gs = model.add('Features not in blanks', all_groups=['Ctrl', 'Treated', 'Blanks'])
    assert gs.excl == ['Ctrl', 'Treated', 'Blanks']
    assert gs.src == [] and gs.incl == []
    assert model.selected is gs


def test_add_then_update_models_blank_filter_default_groupset():
    # Mirrors enumerate_inputs' "no groupsets yet, blank filter on" branch:
    # src = all groups except blanks, excl = [blank group].
    model = GroupSetModel()
    model.add('Features not in blanks', all_groups=['Ctrl', 'Treated', 'Blanks'])
    gs = model.selected
    src = list(gs.excl)
    src.remove('Blanks')
    model.update(model.selected_index, src=src, excl=['Blanks'])
    assert model.selected.src == ['Ctrl', 'Treated']
    assert model.selected.excl == ['Blanks']


def test_update_only_overwrites_given_fields():
    model = GroupSetModel()
    model.add('a', all_groups=['G1'])
    model.update(0, colour='#abcdef')
    assert model[0].colour == '#abcdef'
    assert model[0].excl == ['G1']  # untouched


def test_len_iter_getitem():
    model = GroupSetModel()
    model.add('a')
    model.add('b')
    assert len(model) == 2
    assert [g.name for g in model] == ['a', 'b']
    assert model[1].name == 'b'


# --------------------------------------------------------------------------- #
# legacy round-trip (.mpct backward compatibility)
# --------------------------------------------------------------------------- #

def test_from_legacy_list_of_query_objects_selects_first():
    legacy_items = []
    for nm in ('one', 'two'):
        q = LegacyQuery()
        q.name = nm
        legacy_items.append(q)
    model = GroupSetModel.from_legacy_list(legacy_items)
    assert [g.name for g in model] == ['one', 'two']
    assert model.selected_index == 0


def test_from_legacy_list_of_groupsets_passthrough():
    items = [GroupSet(name='x'), GroupSet(name='y')]
    model = GroupSetModel.from_legacy_list(items)
    assert [g.name for g in model] == ['x', 'y']


def test_from_legacy_list_empty():
    model = GroupSetModel.from_legacy_list([])
    assert len(model) == 0
    assert model.selected_index == -1


def test_to_legacy_list_roundtrip():
    model = GroupSetModel()
    model.add('a', all_groups=['G1'])
    model.add('b')
    items = model.to_legacy_list()
    assert [i.name for i in items] == ['a', 'b']
    restored = GroupSetModel.from_legacy_list(items)
    assert [g.name for g in restored] == ['a', 'b']


# --------------------------------------------------------------------------- #
# build_query_dict
# --------------------------------------------------------------------------- #

def test_build_query_dict_basic_shape():
    model = GroupSetModel()
    model.add('Set1')
    model.update(0, src=['Ctrl'], incl=['Treated'], excl=['Blanks'], colour='#112233')
    qd = build_query_dict(model, graphfilters='')
    assert len(qd) == 1
    name, gs = next(iter(qd.items()))
    assert name == 'Ctrl +Treated -Blanks c=#112233 n=Set1'
    assert gs.src == ['Ctrl'] and gs.incl == ['Treated'] and gs.excl == ['Blanks']


def test_build_query_dict_appends_graphfilters_to_excl():
    model = GroupSetModel()
    model.add('Set1')
    model.update(0, excl=['Blanks'])
    qd = build_query_dict(model, graphfilters='cv rel')
    gs = next(iter(qd.values()))
    assert gs.excl == ['Blanks', 'cv', 'rel']


def test_build_query_dict_does_not_mutate_source_model():
    model = GroupSetModel()
    model.add('Set1')
    model.update(0, excl=['Blanks'])
    build_query_dict(model, graphfilters='cv')
    assert model[0].excl == ['Blanks']  # unaffected by graphfilters merge


def test_build_query_dict_multiple_groupsets_distinct_keys():
    model = GroupSetModel()
    model.add('Set1')
    model.update(0, src=['A'], colour='#000000')
    model.add('Set2')
    model.update(1, src=['B'], colour='#ffffff')
    qd = build_query_dict(model)
    assert len(qd) == 2


def test_build_query_dict_accepts_list_graphfilters():
    """Current code (main.py's enumerate_inputs) builds graphfilters as a
    list (['cv', 'rel']), not the older space-joined string -- both shapes
    must merge into excl identically."""
    model = GroupSetModel()
    model.add('Set1')
    model.update(0, excl=['Blanks'])
    qd = build_query_dict(model, graphfilters=['cv', 'rel'])
    gs = next(iter(qd.values()))
    assert gs.excl == ['Blanks', 'cv', 'rel']


# --------------------------------------------------------------------------- #
# normalize_graphfilters
# --------------------------------------------------------------------------- #

def test_normalize_graphfilters_splits_legacy_string():
    """Older .mpct saves pickled graphfilters as a space-joined string."""
    assert normalize_graphfilters('cv rel') == ['cv', 'rel']


def test_normalize_graphfilters_handles_trailing_space():
    # main.py's old string-building always left a trailing space unless
    # 'insource' (added last, with no trailing space) was checked --
    # str.split() (no args) must not produce a spurious empty token.
    assert normalize_graphfilters('cv ') == ['cv']
    assert normalize_graphfilters('cv rel ') == ['cv', 'rel']


def test_normalize_graphfilters_passes_through_list():
    assert normalize_graphfilters(['cv', 'rel']) == ['cv', 'rel']


def test_normalize_graphfilters_empty_string():
    assert normalize_graphfilters('') == []
