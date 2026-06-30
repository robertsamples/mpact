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


def _scattered_pair_linkage():
    """Hand-built linkage (not derived from real coordinates, so the merge
    order is exact and unambiguous) reproducing the real-data pattern that
    motivated the overlap-based coloring rule: labels P and Q are each
    split across two separate leaves that DON'T merge with each other
    first (P0, Q0, Q1, P1 -- interleaved, not P0+P1 then Q0+Q1), so neither
    P nor Q is monophyletic -- plus an unrelated label R that cleanly joins
    in afterward and should NOT show as part of the tangle.

    Leaves: 0=P, 1=Q, 2=Q, 3=P, 4=R.
    Merge order: (1,2)=Q+Q pure; (0, that)=P+{Q} disjoint bridge;
    (3, that)=P+{P,Q} overlap -- the actual tangle; (4, that)=R+{P,Q}
    disjoint again (R was never part of the P/Q mixing).
    """
    Z = np.array([
        [1, 2, 0.1, 2],   # node 5: Q+Q            -> {Q}        (pure)
        [0, 5, 1.0, 3],   # node 6: P + {Q}         -> {P,Q}      (disjoint)
        [3, 6, 2.0, 4],   # node 7: P + {P,Q}       -> {P,Q}      (overlap!)
        [4, 7, 3.0, 5],   # node 8: R + {P,Q}       -> {P,Q,R}    (disjoint)
    ])
    labels = ['P', 'Q', 'Q', 'P', 'R']
    return Z, labels


def test_purity_summary_both_groups_pure():
    data, labels = _two_clean_groups()
    Z = linkage(data, method='ward')
    n_pure, n_total = purity_summary(Z, labels)
    assert (n_pure, n_total) == (2, 2)


def test_purity_link_color_func_clean_disjoint_groups_stay_neutral_even_at_root():
    data, labels = _two_clean_groups()
    Z = linkage(data, method='ward')
    n_leaves = len(labels)
    color_func = purity_link_color_func(Z, labels)

    # The final merge (root) joins group A's whole clade with group B's
    # whole clade -- their label sets are disjoint ({A} vs {B}, no overlap),
    # so this is a clean join, not evidence either group is non-
    # monophyletic. It must be the neutral color, NOT the polyphyletic one
    # -- two cleanly-resolved groups simply existing in the same tree isn't
    # itself a problem.
    root_node_id = n_leaves + len(Z) - 1
    assert color_func(root_node_id) == 'black'

    # Every internal node strictly below the root is a within-group merge
    # for this dataset -- those links must be the monophyletic color.
    for i in range(len(Z) - 1):
        node_id = n_leaves + i
        assert color_func(node_id) == 'green'


def test_purity_link_color_func_overlap_is_the_only_false_color_and_it_does_not_cascade():
    Z, labels = _scattered_pair_linkage()
    color_func = purity_link_color_func(Z, labels)

    assert color_func(5) == 'green'    # Q+Q, monophyletic
    assert color_func(6) == 'black'    # P + {Q}: disjoint, clean bridge
    assert color_func(7) == 'magenta'  # P + {P,Q}: OVERLAP -- the actual tangle
    # R joining afterward is disjoint from {P,Q} -- R was never part of the
    # P/Q mixing, so this must NOT also render false_color just because it's
    # above (contains) the node-7 tangle. This is the specific behaviour
    # this rule exists for: a real, low-level tangle must not paint every
    # ancestor false_color all the way to the root.
    assert color_func(8) == 'black'


def test_purity_link_color_func_custom_colors():
    Z, labels = _scattered_pair_linkage()
    color_func = purity_link_color_func(Z, labels, true_color='cyan', false_color='magenta', neutral_color='grey')
    assert color_func(5) == 'cyan'
    assert color_func(6) == 'grey'
    assert color_func(7) == 'magenta'
    assert color_func(8) == 'grey'


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
