# Outputs

MPACT saves data into the selected output directory, in a subfolder named
after the peak table file. Most output files share the peak table's
filename with a suffix indicating their content, and index on feature ID,
m/z, and retention time.

| File | Contents |
|---|---|
| `<name>_formatted` | The original peak table, reformatted for MPACT processing. |
| `<name>_merged` | The formatted data with mispicked peaks merged (if that filter was enabled). |
| `<name>_filtered` | The formatted data after all filtering steps. |
| `<name>_groupaverages` | Group-average abundance per feature, in long format. |
| `<name>_summarydata` | All group-level statistics in wide format: average abundance per group, biological/technical relative standard deviations and n, combined relative/absolute standard deviation (error propagation), and effective n (`neff`, Welch–Satterthwaite). |
| `<name>_filtered_source` | The original source-format peak table, row-subset to surviving features (Progenesis-format sources only). See [File Selection](file-selection.md#gnps2-export-re-indexing). |
| `<frag>_reindexed` | The fragment database (MSP/MGF), re-numbered to match the filtered peak table's row order. |
| `iondict` | Everything in the feature-info table plus extra per-feature data for **all** features, including ones that didn't pass filtering: compound ID, m/z, RT, mass defect (`kmd`), which filters it passed, median CV, group presence/absence, max group abundance (and its log, used for plot point opacity), fold change, -log p, -log q, and database-match count. |
| `analysisinfo.txt` | A text summary of the run: basic statistics, filtering results, and all parameters needed to regenerate the analysis later. Also viewable in the [Analysis Info tab](../feature-info.md#analysis-info). |
| `<name>.mpct` | The save file used by **Import from .MPCT session file** to reload this exact analysis (raw data + all parameters). Written atomically (via a temp file + rename), so a failure mid-write can't corrupt a previously good save. |
