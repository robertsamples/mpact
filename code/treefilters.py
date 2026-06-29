"""
MPACT
Copyright 2022, Robert M. Samples, Sara P. Puckett, and Marcy J. Balunas

Qt-free per-column filter criteria and row-matching logic for the feature
search tree (main.py's fillfttree()/treeWidget). Kept separate from the Qt
model/proxy/widget layer (searchtree.py) so the actual matching rules --
substring search, numeric range, multi-select category -- are testable
without a GUI.
"""

from dataclasses import dataclass
from typing import Optional


@dataclass(frozen=True)
class TextFilter:
    """Case-insensitive substring match. Blank text matches every row."""
    text: str = ''

    def matches(self, value):
        if not self.text:
            return True
        return self.text.lower() in str(value).lower()


@dataclass(frozen=True)
class RangeFilter:
    """Inclusive numeric range. Either bound may be ``None`` (unbounded).
    A value that can't be parsed as a number never matches a range filter
    (it's not "in range" if it isn't a number at all).
    """
    minimum: Optional[float] = None
    maximum: Optional[float] = None

    def matches(self, value):
        try:
            number = float(value)
        except (TypeError, ValueError):
            return False
        if self.minimum is not None and number < self.minimum:
            return False
        if self.maximum is not None and number > self.maximum:
            return False
        return True


@dataclass(frozen=True)
class CategoryFilter:
    """For columns whose value is a space-joined set of category tokens
    (e.g. the ``groups`` column: ``' Blanks 0um_Ce'``). Matches if the row
    has at least one token in common with ``allowed``. An empty/``None``
    ``allowed`` matches every row -- "nothing selected" means "no filter
    applied," matching how a checkbox filter list is normally expected to
    behave (vs. "nothing selected" hiding everything).
    """
    allowed: Optional[frozenset] = None

    def matches(self, value):
        if not self.allowed:
            return True
        row_tokens = set(str(value).split())
        return bool(row_tokens & self.allowed)


def row_passes(row_values, filters):
    """``row_values``: a sequence of cell values for one row, indexed by
    column position. ``filters``: ``{column_index: criterion}`` where each
    criterion has a ``.matches(value)`` method (``TextFilter``/
    ``RangeFilter``/``CategoryFilter``, or any object with that interface).
    Columns with no entry in ``filters`` are unfiltered.

    Returns ``True`` iff every active filter's column passes.
    """
    for column, criterion in filters.items():
        if not criterion.matches(row_values[column]):
            return False
    return True


def distinct_category_tokens(values):
    """Given a column of space-joined category strings (e.g. every row's
    ``groups`` value), return the sorted set of distinct individual tokens
    -- used to populate a per-column checkbox filter's option list."""
    tokens = set()
    for value in values:
        tokens.update(str(value).split())
    return sorted(tokens)
