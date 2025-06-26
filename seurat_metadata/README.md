# Seurat Metadata Utilities

This folder contains scripts for exporting and summarizing metadata from the main Seurat object used in the project.

## Scripts

- **export_full_metadata.R**
  - Exports the full metadata (`@meta.data`) from the Seurat object to a CSV file (`seurat_metadata_full.csv`).
  - Output: All cell-level metadata for downstream analysis or inspection.

- **seurat_metadata_summary.R**
  - Summarizes the Seurat object and writes a summary to `seurat_summary.txt`.
    - Lists assays, number of cells, number of features (genes), metadata columns, identity levels, and dimensional reductions.
  - Also saves a sample (first 5 rows) of the metadata as `seurat_metadata_sample.csv` for quick inspection.

## Usage

Run these scripts in R or submit them as part of a batch job to generate metadata summaries and exports for quality control, reporting, or downstream analysis.

---

These scripts help ensure transparency and reproducibility by documenting the structure and content of the Seurat object used in the analysis.
