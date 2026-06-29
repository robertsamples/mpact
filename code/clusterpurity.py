"""
MPACT
Copyright 2022, Robert M. Samples, Sara P. Puckett, and Marcy J. Balunas

Qt-free dendrogram "purity" coloring: a branch is colored green if every
leaf beneath it shares the same group label -- i.e. that group is a
monophyletic clade, it clustered together before merging with anything
else -- and left at the default color otherwise. Used by the dendrogram tab
to make it visually obvious whether technical replicates of one Sample
cluster tightly together, and separately whether biological replicates of
one Biolgroup are well separated from other groups.

This module is Qt-free and unit-tested (see ``tests/test_clusterpurity.py``).
"""


def purity_link_color_func(Z, leaf_labels, true_color='green', false_color='black'):
    """Build a ``link_color_func`` for ``scipy.cluster.hierarchy.dendrogram``.

    Args:
        Z: linkage matrix (``scipy.cluster.hierarchy.linkage`` or
            fastcluster's drop-in) built on observations in the same order
            as ``leaf_labels``.
        leaf_labels: sequence, length == number of observations clustered by
            ``Z``, giving each leaf's group label (e.g. its Sample or
            Biolgroup), in the same order as the data passed to ``linkage``.

    Returns:
        callable: ``link_color_func(k)`` as expected by ``dendrogram``'s
        ``link_color_func`` argument -- for link index ``k``
        (``len(leaf_labels) <= k``), returns ``true_color`` if every leaf
        descending from that link shares one label, else ``false_color``.
    """
    n_leaves = len(leaf_labels)
    leaf_label_sets = {i: {leaf_labels[i]} for i in range(n_leaves)}
    colors = {}
    for i, row in enumerate(Z):
        a, b = int(row[0]), int(row[1])
        node_id = n_leaves + i
        merged = leaf_label_sets[a] | leaf_label_sets[b]
        leaf_label_sets[node_id] = merged
        colors[node_id] = true_color if len(merged) == 1 else false_color
    return lambda k: colors.get(k, false_color)


def purity_summary(Z, leaf_labels):
    """Count how many distinct group labels form one pure clade each.

    A label is "pure" only if *every* leaf carrying that label ends up
    together in one clade before that clade merges with any other leaf --
    i.e. the group is exactly monophyletic in the dendrogram. (A node whose
    descendants are a uniform-but-incomplete subset of a label -- e.g. 2 of
    a Sample's 3 technical replicates -- does NOT count: the third
    replicate clustering elsewhere means that Sample isn't really pure.)

    Returns:
        (n_pure, n_total): number of distinct labels that are fully pure
        clades, out of the total number of distinct labels in
        ``leaf_labels``.
    """
    n_leaves = len(leaf_labels)
    leaf_index_sets = {i: frozenset((i,)) for i in range(n_leaves)}
    target_sets = {
        label: frozenset(i for i in range(n_leaves) if leaf_labels[i] == label)
        for label in set(leaf_labels)
    }
    pure_labels = set()
    for i, row in enumerate(Z):
        a, b = int(row[0]), int(row[1])
        node_id = n_leaves + i
        merged = leaf_index_sets[a] | leaf_index_sets[b]
        leaf_index_sets[node_id] = merged
        for label, target in target_sets.items():
            if merged == target:
                pure_labels.add(label)
    return len(pure_labels), len(target_sets)
