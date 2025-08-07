# ==============================================================================
# SCRIPT: 02_prepare_expression_data.R
#
# PURPOSE:
#   This script is the second foundational step in the analysis pipeline. Its
#   sole purpose is to extract the necessary single-cell expression data from
#   the primary Seurat object and save it in a simple, portable format (CSV).
#
#   This decouples the computationally heavy Seurat object from the downstream
#   analysis scripts, improving performance and modularity. It provides a
#   standardized and portable data format for all subsequent noise and
#   expression analyses.
#
# OUTPUTS:
#   1. 02_expression_matrix.csv: A matrix of log-normalized gene expression
#      values with genes as rows and cells as columns.
#   2. 02_cell_metadata.csv: A table containing all metadata for each cell,
#      including cell type annotations, UMAP coordinates, and batch information.
# ==============================================================================

# --- 1. SETUP: Load libraries and define paths ---

# Ensure the Seurat package is installed and loaded
if (!requireNamespace("Seurat", quietly = TRUE)) {
  # Install Seurat from the official CRAN repository
  install.packages("Seurat", repos = "https://cloud.r-project.org/")
}
library(Seurat)

# --- Define I/O Paths ---
# Using variables for paths makes the script more portable and easier to manage.

# INPUT: The path to the integrated Seurat object
seurat_path <- "/group/sms029/Oliva_dataset/integrated_col_trajectories_colonly.rds"

# OUTPUT: Define the directory for this script's outputs
output_dir <- "/group/sms029/mnieuwenh/gbm_noise_analysis/results/02_prepare_expression_data"
# Create the directory if it doesn't exist, without showing a warning if it does
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# Define the full paths for the output files, named according to the script number
expr_matrix_outfile <- file.path(output_dir, "02_expression_matrix.csv")
cell_meta_outfile <- file.path(output_dir, "02_cell_metadata.csv")


# --- 2. DATA LOADING ---

message(paste("Loading Seurat object from:", seurat_path))
# Check if the input file exists before attempting to load it
if (!file.exists(seurat_path)) {
  stop("Error: Seurat object file not found. Please check the path.")
}
seurat_obj <- readRDS(seurat_path)
message("Seurat object loaded successfully.")


# --- 3. EXPORT EXPRESSION MATRIX ---

message("Extracting the expression matrix from the 'RNA' assay...")

# Use the recommended GetAssayData function.
# This extracts the log-normalized counts from the 'data' layer of the 'RNA' assay.
expr_mat <- GetAssayData(seurat_obj, assay = "RNA", layer = "data")

# Validation: Ensure the matrix is not empty and has correct dimensions
if (nrow(expr_mat) == 0 || ncol(expr_mat) == 0) {
  stop("Error: The extracted expression matrix is empty.")
}
# Validation: Confirm that the cell barcodes in the matrix match the Seurat object
if (!all(colnames(expr_mat) == colnames(seurat_obj))) {
    stop("Error: Mismatch between matrix column names and Seurat object cell names.")
}

message(paste("Writing expression matrix with", nrow(expr_mat), "genes and", ncol(expr_mat), "cells."))
# write.csv is a standard, robust way to save this data.
# We convert to a standard matrix to ensure compatibility.
# row.names=TRUE is essential to save the gene IDs.
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
# row.names=TRUE is important as the rownames are the cell barcodes, which link to the expression matrix.
write.csv(cell_metadata, file = cell_meta_outfile, row.names = TRUE)
message(paste("Cell metadata saved to:", cell_meta_outfile))


# --- 5. FINAL CONFIRMATION ---
message("Script finished. All necessary expression data has been exported.")