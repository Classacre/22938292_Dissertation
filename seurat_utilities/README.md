# Seurat Utilities

This folder contains scripts for visualization, filtering, and cell counting using the main Seurat object in the project.

## Scripts

- **Seurat_Visualization.R**
  - Generates a variety of plots from the Seurat object, including UMAPs, violin plots, dot plots, heatmaps, and ridge plots.
  - Plots are grouped by cell type, genotype, cluster, and other metadata where available.
  - All plots are saved to the `seurat_plots` directory.

- **cell_counts_by_type_and_lineage.R**
  - Counts the number of cells per cell type (`identity`) and per lineage (`lineage`) in the Seurat object.
  - Outputs: `cell_counts_by_celltype.csv` and `cell_counts_by_lineage.csv`.

- **remove_na_genotype.R**
  - Filters the Seurat object to only include cells with genotype == 'col' (removes NA and other genotypes).
  - Overwrites the original Seurat object with the filtered version.

## Usage

Run these scripts in R or submit them as batch jobs to generate plots, filter cells, or count cells by type/lineage for downstream analysis and reporting.

---

These utilities help with quality control, visualization, and data exploration for single-cell RNA-seq analysis in Arabidopsis root.
