# Data Review

Open the **Plots** pane (left sidebar) to access all data visualization
tools. Plots can be zoomed, panned, and saved using the toolbar beneath
each one. Clicking a feature in most plots highlights it and updates the
other plots/feature-info display to match — see
[Feature Info](../feature-info.md) for selection details.

The Data Review tab summarizes overall data quality.

## Summary

Reports an estimate of data quality based on where the dataset's CV50 (the
CV of the median feature) falls relative to the theoretical maximum
reproducibility, plus the number of features removed by each filtering
step.

## CV Plot

Shows the distribution of mean and median CV values across the dataset —
useful for assessing reproducibility and spotting peak-picking/alignment
issues.

- A large gap between the mean and median CV curves indicates a highly
  skewed distribution (some sets of technical replicates are much less
  reproducible than others).
- A very shallow curve, with a low percentage of features meeting common
  CV thresholds (0.2–0.3), can indicate alignment problems.
- Marked multimodality (a sharp rise partway through the distribution) can
  indicate inconsistent detection across technical replicate sets.

## Feature Plot

Plots every feature (m/z vs. retention time), coloured by which filtering
stage removed it, if any. Features that passed all filters are shown in
black as high-quality features.

## Treemap

Shows the number and percentage of features removed by each filtering
step, as a treemap.
