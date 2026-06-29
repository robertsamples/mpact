"""Unit tests for dendrogram purity coloring (``clusterpurity.py``)."""

import numpy as np
from scipy.cluster.hierarchy import linkage

from clusterpurity import purity_link_color_func, purity_summary


def _two_clean_groups():
    """6 points: 3 tightly clustered near (0, 0) labeled 'A', 3 tightly
    clustered near (10, 10) labeled 'B' -- each group should merge with
    itself long before the two groups merge with each other.
    """
    data = np.array([
        [0.0, 0.0], [0.1, 0.0], [0.0, 0.1],
        [10.0, 10.0], [10.1, 10.0], [10.0, 10.1],
    ])
    labels = ['A', 'A', 'A', 'B', 'B', 'B']
    return data, labels


def test_purity_summary_both_groups_pure():
    data, labels = _two_clean_groups()
    Z = linkage(data, method='ward')
    n_pure, n_total = purity_summary(Z, labels)
    assert (n_pure, n_total) == (2, 2)


def test_purity_link_color_func_roots_to_false_color_leaves_to_true_color():
    data, labels = _two_clean_groups()
    Z = linkage(data, method='ward')
    n_leaves = len(labels)
    color_func = purity_link_color_func(Z, labels)

    # The final merge (root) joins group A's clade with group B's clade --
    # that link must NOT be the "pure" color.
    root_node_id = n_leaves + len(Z) - 1
    assert color_func(root_node_id) == 'black'

    # Every internal node strictly below the root is a within-group merge
    # for this dataset (each group's 3 points cluster before the cross-group
    # merge) -- those links must be the "pure" color.
    for i in range(len(Z) - 1):
        node_id = n_leaves + i
        assert color_func(node_id) == 'green'


def test_purity_link_color_func_custom_colors():
    data, labels = _two_clean_groups()
    Z = linkage(data, method='ward')
    color_func = purity_link_color_func(Z, labels, true_color='cyan', false_color='grey')
    n_leaves = len(labels)
    root_node_id = n_leaves + len(Z) - 1
    assert color_func(root_node_id) == 'grey'
    assert color_func(n_leaves) == 'cyan'


def test_purity_summary_one_mismatched_leaf_breaks_purity_for_its_group():
    # Same as the clean two-group case, but one of group A's points is
    # actually closest to group B -- A should no longer be reported pure
    # (its leaves don't all merge together before meeting a 'B' leaf),
    # while B (unaffected) should still be pure.
    data = np.array([
        [0.0, 0.0], [0.1, 0.0], [9.9, 9.9],  # last "A" point sits with B
        [10.0, 10.0], [10.1, 10.0], [10.0, 10.1],
    ])
    labels = ['A', 'A', 'A', 'B', 'B', 'B']
    Z = linkage(data, method='ward')
    n_pure, n_total = purity_summary(Z, labels)
    assert n_pure == 1
    assert n_total == 2


def test_purity_link_color_func_unknown_link_id_falls_back_to_false_color():
    data, labels = _two_clean_groups()
    Z = linkage(data, method='ward')
    color_func = purity_link_color_func(Z, labels)
    assert color_func(99999) == 'black'
