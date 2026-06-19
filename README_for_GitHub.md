# M. abscessus MOM/CW proteomics analysis

This repository contains the R script used for the R-based proteomics analysis associated with the manuscript:

**Integrated subcellular fractionation proteomics and structural prediction analysis identifies putative mycolic acid outer membrane proteins and ESX-associated toxin systems in *Mycobacterium abscessus*.**

## Files

- `Mabscessus_MOM_proteomics_analysis_clean.R`  
  Cleaned R script for FragPipe MaxLFQ processing, replicate-level filtering, QC/PCA outputs, limma differential abundance analysis, and volcano plot generation.

## What this script does

The script performs the following R-based steps:

1. Imports FragPipe `combined_protein` output.
2. Extracts protein-level MaxLFQ intensity columns.
3. Parses sample metadata from sample names.
4. Treats zero MaxLFQ intensity values as missing values.
5. Filters proteins based on reproducible detection in at least 2 of 3 biological replicates in at least one strain-fraction condition.
6. Generates QC outputs, including quantified protein counts per sample, PCA plots, scree values, and a Pearson sample correlation matrix.
7. Performs limma differential abundance analysis for manuscript-relevant contrasts.
8. Exports limma contrast output tables.
9. Generates volcano plots for:
   - P8 versus P100 membrane fraction comparisons.
   - Fig. 2-style P8/P100 marker volcano plots highlighting Antigen 85, MspA-like proteins, ATP synthase, and NDH-1 markers.
   - Mutant culture filtrate versus wild-type culture filtrate comparisons.
   - Culture filtrate versus total lysate comparisons.

## What this script does not do

Presence/absence heatmaps, hybrid heatmap assembly, and ESX-dependency classifications were performed in Microsoft Excel using FragPipe detection outputs and limma contrast tables. Those Excel-based steps are intentionally not reproduced in this R script.

## Input data

Place the FragPipe `combined_protein` file in a folder called `data/` and name it either:

- `combined_protein.tsv`, or
- update `combined_protein_file` in the script to point to your file.

Optional: to highlight all curated pMOMPs in the candidate P8/P100 volcano plot, add a file:

```text
data/pMOMP_ids.txt
```

with one UniProt ID per line.

## Output

Running the script creates an `outputs/` folder containing:

- QC tables and plots.
- Filtered log2 MaxLFQ matrix.
- limma contrast tables.
- Volcano plot PDFs/PNGs and source data tables.
- `sessionInfo.txt`.

## Notes

The script is intended to document the analysis workflow used for the manuscript. It may require minor edits to file names depending on how the FragPipe output files are named locally.
