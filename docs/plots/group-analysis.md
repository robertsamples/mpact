# Group Analysis

Visualizes differences between treatment groups via set analysis,
correlation analysis, and hierarchical clustering.

## UpSet Plot

Shows feature-set overlap across treatment groups: a bar chart of feature
counts per group (left), and the number of features present in each
combination of groups (top bar chart + dot matrix).

![UpSet plot](../images/upset-plot.png)
*MPACT UpSet plot showing the distribution of features across sample sets.*

## Spearman Correlation Matrix

Pairwise Spearman correlation between every group, useful for evaluating
overall metabolomic similarity at a glance. Colour scheme is configurable
in the plot options dialog.

![Spearman correlation matrix](../images/spearman-matrix.png)
*MPACT Spearman correlation matrix.*

## Dendrogram

Hierarchical clustering analysis (Ward's method) of samples and treatment
groups — useful both for assessing metabolomic similarity and for
sanity-checking your CV filter threshold (see
[Filtering Settings](../user-guide/filtering-settings.md#cv-reproducibility-filtering)).
For datasets where biological/treatment differences should exceed
technical noise, technical replicates should cluster together after
filtering.

A settings bar above the plot controls how it's drawn — all four options
are local to this tab (not the general plot options dialog) and redraw
immediately when changed:

**View** — which leaves to cluster:

- **Technical Replicates** (default): every injection is its own leaf,
  letting you see whether each sample's individual injections agree with
  each other.
- **Biological Replicates**: technical replicates are averaged together
  first (one leaf per sample), so the plot reflects clustering of
  biological/treatment groups without technical noise.

**Color** — how branches are colored:

- **Purity** (default): a branch is colored **green** if every leaf beneath
  it belongs to the same sample (Technical Replicates view) or the same
  treatment group (Biological Replicates view) — i.e. it's correctly,
  unambiguously clustered. A branch is colored **red** if it's the specific
  point where two different samples/groups' leaves are proven to overlap
  (some of that sample's/group's replicates are on each side of the
  split) — a real sign of poor clustering, not just "still mixed from
  somewhere lower in the tree." Every other branch (a clean join of two
  unrelated, already-resolved regions) stays black, even if it sits above a
  red branch elsewhere in the tree — so a single tangled sample doesn't
  paint the whole tree red. The plot title reports how many
  samples/groups are *fully* correctly clustered (e.g. "7/9 samples'
  replicates clustered together").
- **None**: a plain, uncolored dendrogram with no title — useful if you
  just want the clustering shape without the QC overlay.

**Bootstrap** — when checked (default on), runs bootstrap resampling
(1000 iterations) and annotates the dendrogram with approximately-unbiased
(AU) p-values and bootstrap probabilities (BP) at each branch point. AU
values above 95 are generally considered statistically significant.
Bootstrap computation uses `fastcluster` if it's installed (falling back to
SciPy's hierarchical clustering otherwise) for substantially faster linkage
on large datasets. Uncheck for a faster, unannotated dendrogram.

**Use Sample/Group Names** — when checked, leaf labels switch from the raw
injection/file names (which can be long or uninformative) to
`<group>_b<#>_s<#>` (Technical Replicates view) or `<group>_b<#>`
(Biological Replicates view), where `b<#>` numbers each biological
replicate (sample) within its group and `s<#>` numbers each technical
replicate (injection) within its sample.

![Dendrogram](../images/dendrogram.png)
*MPACT dendrogram after filtering, showing correct clustering of most
technical replicates and all biological groups. Approximately-unbiased
(AU) p-values greater than 95 are considered statistically significant.*
