import pytest

from treefilters import (
    CategoryFilter,
    RangeFilter,
    TextFilter,
    distinct_category_tokens,
    row_passes,
)


# --------------------------------------------------------------------------- #
# TextFilter
# --------------------------------------------------------------------------- #

def test_text_filter_blank_matches_everything():
    assert TextFilter('').matches('anything')
    assert TextFilter('').matches('')


def test_text_filter_case_insensitive_substring():
    f = TextFilter('gran')
    assert f.matches('Granaticin C')
    assert f.matches('GRANATICIN')
    assert not f.matches('Coelichelin')


# --------------------------------------------------------------------------- #
# RangeFilter
# --------------------------------------------------------------------------- #

def test_range_filter_no_bounds_matches_everything():
    f = RangeFilter()
    assert f.matches(0)
    assert f.matches(-999.5)
    assert f.matches(1e6)


def test_range_filter_inclusive_bounds():
    f = RangeFilter(minimum=1, maximum=5)
    assert f.matches(1)
    assert f.matches(5)
    assert f.matches(3)
    assert not f.matches(0.999)
    assert not f.matches(5.001)


def test_range_filter_one_sided_bounds():
    assert RangeFilter(minimum=10).matches(1000)
    assert not RangeFilter(minimum=10).matches(9)
    assert RangeFilter(maximum=10).matches(-1000)
    assert not RangeFilter(maximum=10).matches(11)


def test_range_filter_non_numeric_value_never_matches():
    f = RangeFilter(minimum=0, maximum=10)
    assert not f.matches('not a number')
    assert not f.matches(None)


# --------------------------------------------------------------------------- #
# CategoryFilter
# --------------------------------------------------------------------------- #

def test_category_filter_no_selection_matches_everything():
    assert CategoryFilter(None).matches(' Blanks 0um_Ce')
    assert CategoryFilter(frozenset()).matches(' Blanks 0um_Ce')


def test_category_filter_matches_any_overlapping_token():
    f = CategoryFilter(frozenset({'0um_Ce'}))
    assert f.matches(' Blanks 0um_Ce')
    assert not f.matches(' Blanks Media')


def test_category_filter_multiple_selected_tokens_is_an_or():
    f = CategoryFilter(frozenset({'0um_Ce', '250um_Ce'}))
    assert f.matches(' 0um_Ce')
    assert f.matches(' 250um_Ce')
    assert not f.matches(' Media')


# --------------------------------------------------------------------------- #
# row_passes
# --------------------------------------------------------------------------- #

def test_row_passes_with_no_filters():
    assert row_passes(['anything', 1, 2], {})


def test_row_passes_requires_all_active_filters_to_match():
    row = ['Granaticin C', 101.5, ' 0um_Ce']
    filters = {
        0: TextFilter('gran'),
        1: RangeFilter(minimum=100, maximum=200),
        2: CategoryFilter(frozenset({'0um_Ce'})),
    }
    assert row_passes(row, filters)


def test_row_fails_if_any_single_filter_fails():
    row = ['Granaticin C', 101.5, ' 0um_Ce']
    filters = {
        0: TextFilter('gran'),
        1: RangeFilter(minimum=200, maximum=300),  # out of range
        2: CategoryFilter(frozenset({'0um_Ce'})),
    }
    assert not row_passes(row, filters)


# --------------------------------------------------------------------------- #
# distinct_category_tokens
# --------------------------------------------------------------------------- #

def test_distinct_category_tokens_sorted_and_deduplicated():
    values = [' Blanks 0um_Ce', ' 0um_Ce Media', ' Media', '']
    assert distinct_category_tokens(values) == ['0um_Ce', 'Blanks', 'Media']


def test_distinct_category_tokens_empty_input():
    assert distinct_category_tokens([]) == []
