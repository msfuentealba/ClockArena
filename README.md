# ClockArena
The two R scripts reproduce the main figures and associated supplementary figures and tables from saved analysis inputs. These inputs include the curated clock catalogue, clock–CpG membership matrix, integrated sample metadata, previously calculated clock estimates and cell fractions, and benchmark results. The scripts summarize these inputs and perform additional analyses needed for individual panels; they do not process raw methylation files or calculate the original epigenetic clock estimates.

**Figure 1: overview of the benchmark and its data**

`build_figure1.R` generates the eight panels describing the study workflow, clock characteristics and sample composition:

- **Panel A** illustrates the workflow from clock and metadata curation through methylation processing, clock calculation, cell-composition estimation, integration and benchmarking.
- **Panels B–D** summarize publication year and training tissue, examine the relationship between publication year and the number of required CpGs, and identify the most frequently shared CpGs with their nearest-gene annotations.
- **Panels E–H** summarize sample age and sex distributions, geographic coverage, and the most frequent disease and intervention/exposure labels with their comparison-group counts.

The script also examines the association between publication year and median age-prediction error in Supplementary Figure 1. It exports Figure 1 as PDF and PNG, Supplementary Figure 1 as PDF, and Supplementary Tables 1–2 as CSV files.

**Figure 2: comparisons across benchmark domains**

`build_figure2.R` combines saved benchmark results with calculations from the integrated data to generate eleven panels:

- **Panel A** calculates within-study Pearson correlations between adjusted age-acceleration estimates, averages valid correlations equally across studies, and clusters clocks by their agreement.
- **Panel B** displays the mean absolute change in age-acceleration residuals associated with cell-composition adjustment.
- **Panel C** displays age-group-specific mean absolute errors for clocks ranking among the three most accurate in at least one age bin.
- **Panel D** presents disease-discrimination AUCs for selected study–disease comparisons.
- **Panels E–F** compare prediction errors across study–country combinations and show country-specific results for the clock achieving the lowest error most frequently.
- **Panels G–H** summarize within-person deviations from the expected change over the sampling interval, including a detailed one-week comparison in GSE227809.
- **Panel I** compares case-minus-control age-acceleration differences for clocks with the ten highest and ten lowest responsiveness scores.
- **Panel J** calculates Spearman correlations between rankings across the seven benchmark domains.
- **Panel K** compares domain rankings between clocks trained to predict chronological age and those trained on mortality-related outcomes, displaying individual clocks, medians and interquartile ranges.

Additional outputs examine cell-composition influence, reliability, responsiveness and analysis coverage. The script exports Figure 2 as PDF and PNG, Supplementary Figures 2–6 as PDFs, and Supplementary Tables 3–15 as CSV files. The coverage tables additionally use saved integration and cell-composition audit records.

Both scripts run from the project root and save their outputs in the paper output directory. Figure 1 should be generated first because Figure 2 uses its Supplementary Table 1 to classify clocks by training outcome. Combining the supplementary CSV tables into the final Excel workbook is a separate step; workbook-specific column selections do not alter the source tables or figure calculations.
