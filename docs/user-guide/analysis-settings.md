# Analysis Settings & Plot Feature Sets

## Group presence/absence threshold

A relative or absolute group-parsing threshold determines whether a
feature counts as "present" in a given biological group, and is also used
to filter features that are really just present in a blank/control group.
A relative threshold is a fraction of the feature's abundance in whichever
group it's most abundant in; an absolute threshold is a raw or normalized
abundance value.

## Selecting which plots/calculations to generate

You can select which plots and calculations to generate based on
preference — deselecting most optional outputs does **not** meaningfully
speed up analysis, since most of the runtime is spent in filtering, not
plot generation.

For **Volcano plots**, both fold-change and t-test options must be
selected, with an experimental and a control/treatment group specified
(available groups are detected automatically from your metadata file).
Benjamini-Hochberg false-discovery-rate correction (based on the number of
plotted features) can optionally be enabled to reduce type II error from
multiple-hypothesis testing.

## Plot Feature Sets (groupsets)

"Plot Feature Sets" control which colour a feature is plotted in across
MPACT's plots, based on which biological groups it was detected in.

- Click **+** to create a new feature set; double-click its default name
  to rename it. Click **−** to remove the selected set.
- Click between sets in the list to switch which one you're editing.
- When a set is created, the biological groups from your metadata file
  populate an "Available groups" list. Click to select a group
  (Shift-click / Ctrl-click for multiple, or a range), then drag into
  **"May be in"** or **"Must be in"**.
- Click the colour button to open the colour picker and assign that set's
  plotting colour.

A feature is plotted in a set's colour if it is detected in **all** of that
set's "Must be in" groups, and **none** of the groups still in "Available
groups" (i.e. not dragged into either list). Groups in "May be in" are
optional — a feature can be present or absent in them without affecting
whether it matches the set.

To plot every feature in the dataset regardless of group, drag every group
(except the blank/control group, if blank filtering is enabled) into "May
be in".

!!! tip "Behind the scenes"
    Internally, each Plot Feature Set is a `GroupSet` object managed by a
    small model/collection class (`GroupSetModel`) rather than a bare list
    + selected-index pair. This is purely an implementation detail (see
    [Development](../development.md)) — old `.mpct` save files still load
    correctly, with their saved feature sets converted into the current
    representation automatically.
