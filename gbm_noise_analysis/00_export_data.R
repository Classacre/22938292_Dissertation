# ==============================================================================
# SCRIPT: 00_export_data.R
#
# PURPOSE:
#   This script serves as the first step in the analysis pipeline. Its sole
#   purpose is to extract the necessary data from the primary Seurat object
#   and save it in a simple, portable format (CSV). This decouples the
#   computationally heavy Seurat object from the downstream analysis scripts,
#   making them faster and more modular.
#
# OUTPUTS:
#   1. expression_matrix.csv: A matrix of log-normalized gene expression values
#      with genes as rows and cells as columns.
#   2. cell_metadata.csv: A table containing all metadata for each cell,
#      including cell type annotations, UMAP coordinates, and any batch info.
# ==============================================================================

# --- 1. SETUP: Load libraries and define paths ---

# Ensure required packages are installed and loaded
if (!requireNamespace("Seurat", quietly = TRUE)) {
  install.packages("Seurat", repos = "https://cloud.r-project.org/")
}
library(Seurat)

# --- Define I/O Paths ---
# Using variables for paths makes the script more portable and easier to manage.
# INPUT
seurat_path <- "/group/sms029/Oliva_dataset/integrated_col_trajectories_colonly.rds"

# OUTPUTS
output_dir <- "/group/sms029/mnieuwenh/gbm_noise_analysis"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

expr_matrix_outfile <- file.path(output_dir, "expression_matrix.csv")
cell_meta_outfile <- file.path(output_dir, "seurat_metadata_full.csv")


# --- 2. DATA LOADING ---

message(paste("Loading Seurat object from:", seurat_path))
if (!file.exists(seurat_path)) {
  stop("Error: Seurat object file not found. Please check the path.")
}
seurat_obj <- readRDS(seurat_path)
message("Seurat object loaded successfully.")


# --- 3. EXPORT EXPRESSION MATRIX ---

message("Extracting the expression matrix from the 'RNA' assay...")

# Use the recommended GetAssayData function for Seurat v5.
# This extracts the log-normalized counts from the 'data' layer.
expr_mat <- GetAssayData(seurat_obj, assay = "RNA", layer = "data")

# Validation: Ensure the matrix is not empty and has correct dimensions
if (nrow(expr_mat) == 0 || ncol(expr_mat) == 0) {
  stop("Error: The extracted expression matrix is empty.")
}
if (!all(colnames(expr_mat) == colnames(seurat_obj))) {
    stop("Error: Mismatch between matrix column names and Seurat object cell names.")
}

message(paste("Writing expression matrix with", nrow(expr_mat), "genes and", ncol(expr_mat), "cells."))
# write.csv is a standard, robust way to save this data.
# row.names=TRUE ensures gene IDs are saved.
write.csv(as.matrix(expr_mat), file = expr_matrix_outfile, row.names = TRUE)
message(paste("Expression matrix saved to:", expr_matrix_outfile))


# --- 4. EXPORT CELL METADATA ---

message("Extracting cell metadata...")
# The metadata is stored in the meta.data slot of the Seurat object.
cell_metadata <- seurat_obj@meta.data

# Validation: Ensure the metadata corresponds to the expression matrix
if (nrow(cell_metadata) != ncol(expr_mat)) {
  stop("Error: Number of cells in metadata does not match the expression matrix.")
}

message(paste("Writing metadata for", nrow(cell_metadata), "cells."))
# write.csv is also used here for consistency.
# row.names=TRUE is important as the rownames are the cell barcodes.
write.csv(cell_metadata, file = cell_meta_outfile, row.names = TRUE)
message(paste("Cell metadata saved to:", cell_meta_outfile))


# --- 5. FINAL CONFIRMATION ---
message("Script finished. All necessary data has been exported.")