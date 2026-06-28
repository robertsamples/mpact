# Running an Analysis

Once File Selection, Filtering Settings, and Analysis Settings are
configured, click **Run** (bottom-left of the main window). The status bar
displays "Analysis Complete" when finished.

The heavy data-processing/statistics step runs on a background thread, so
the rest of the UI stays responsive while it works — for a large dataset
(thousands of features, multiple groups, bootstrap dendrogram support
enabled) this can take on the order of a few minutes; smaller datasets
typically finish in seconds. Plot rendering itself happens on the main
thread once the background computation finishes (matplotlib isn't
thread-safe), so you'll briefly see the UI busy while plots draw, even
though the long compute step itself doesn't block you.

The first time you open the [Feature Info](../feature-info.md) tab after a
run, MPACT also searches the bundled Natural Products Atlas database for
mass matches, which can take a moment.

At the end of a successful run, MPACT automatically:

- Writes the `.mpct` save file and processed-data outputs to your chosen
  output directory (see [Outputs](outputs.md)).
- If a fragment database was provided, writes a re-indexed copy of it and
  a filtered copy of the source peak table for GNPS2 submission (see
  [File Selection](file-selection.md#gnps2-export-re-indexing)).

If an individual plot fails to generate (e.g. due to a missing optional
input), MPACT logs the error and continues generating the rest — a single
broken plot won't abort the whole run or your saved session.
