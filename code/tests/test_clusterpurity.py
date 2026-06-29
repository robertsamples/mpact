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


def test_purity_link_color_func_root_bridges_two_pure_clades():
    data, labels = _two_clean_groups()
    Z = linkage(data, method='ward')
    n_leaves = len(labels)
    color_func = purity_link_color_func(Z, labels)

    # The final merge (root) joins group A's whole clade with group B's
    # whole clade -- both children are themselves pure, so this is exactly
    # the "bridge" merge (the one and only point the two groups meet) and
    # must be the false/bridge color, not the neutral one.
    root_node_id = n_leaves + len(Z) - 1
    assert color_func(root_node_id) == 'red'

    # Every internal node strictly below the root is a within-group merge
    # for this dataset (each group's 3 points cluster before the cross-group
    # merge) -- those links must be the "pure" color.
    for i in range(len(Z) - 1):
        node_id = n_leaves + i
        assert color_func(node_id) == 'green'


def test_purity_link_color_func_does_not_cascade_red_up_the_whole_tree():
    # Hand-built linkage (not derived from real coordinates, so the merge
    # order is exact and unambiguous): 4 groups of 2 leaves each --
    # A=(0,1), B=(2,3), C=(4,5), D=(6,7). Merge order: each group merges
    # with itself first (pure), then A+B bridge, then C+D bridge, then the
    # root merges the two already-impure (A+B) and (C+D) clades together.
    # Z columns: [child1, child2, distance (unused), count (unused)].
    Z = np.array([
        [0, 1, 0.1, 2],    # node 8:  A+A  (pure)
        [2, 3, 0.1, 2],    # node 9:  B+B  (pure)
        [8, 9, 5.0, 4],    # node 10: A+B  (bridge -- both children pure)
        [4, 5, 0.1, 2],    # node 11: C+C  (pure)
        [6, 7, 0.1, 2],    # node 12: D+D  (pure)
        [11, 12, 5.0, 4],  # node 13: C+D  (bridge -- both children pure)
        [10, 13, 50.0, 8],  # node 14: (A+B)+(C+D) -- root
    ])
    labels = ['A', 'A', 'B', 'B', 'C', 'C', 'D', 'D']
    color_func = purity_link_color_func(Z, labels)

    assert color_func(8) == 'green'   # A+A
    assert color_func(9) == 'green'   # B+B
    assert color_func(10) == 'red'    # A+B bridge
    assert color_func(11) == 'green'  # C+C
    assert color_func(12) == 'green'  # D+D
    assert color_func(13) == 'red'    # C+D bridge
    # The root combines two ALREADY-impure clades -- no new bridge event,
    # so it must NOT also render red (that's the "entire tree turns red"
    # behaviour this function is specifically built to avoid).
    assert color_func(14) == 'black'


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


def test_purity_link_color_func_unknown_link_id_falls_back_to_neutral_color():
    data, labels = _two_clean_groups()
    Z = linkage(data, method='ward')
    color_func = purity_link_color_func(Z, labels)
    assert color_func(99999) == 'black'
