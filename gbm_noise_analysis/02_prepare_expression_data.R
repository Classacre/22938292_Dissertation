# ==============================================================================
# SCRIPT: 02_prepare_expression_data.R
#
# PURPOSE:
#   This script extracts single-cell expression data and metadata from the
#   primary Seurat object and saves them in portable CSV format. This decouples
#   heavy Seurat objects from downstream analysis scripts.
#
#   **MODIFICATION**: This version also parses the 'identity' column to create
#   separate columns for the base cell type and the ERS state, which is
#   critical for the new ERS-focused analysis.
# ==============================================================================

# --- 1. SETUP: Load libraries and define paths ---
if (!requireNamespace("Seurat", quietly = TRUE)) {
  install.packages("Seurat", repos = "https://cloud.r-project.org/")
}
if (!requireNamespace("dplyr", quietly = TRUE)) {
  install.packages("dplyr", repos = "https://cloud.r-project.org/")
}
library(Seurat)
library(dplyr)

# --- Define I/O Paths ---
seurat_path <- "/group/sms029/Oliva_dataset/integrated_col_trajectories_colonly.rds"
output_dir <- "/group/sms029/mnieuwenh/gbm_noise_analysis/results/02_prepare_expression_data"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
expr_matrix_outfile <- file.path(output_dir, "02_expression_matrix.csv")
cell_meta_outfile <- file.path(output_dir, "02_cell_metadata_annotated.csv")
summary_file_path <- file.path(output_dir, "02_expression_data_summary.txt")

# --- 2. DATA LOADING ---
message(paste("Loading Seurat object from:", seurat_path))
if (!file.exists(seurat_path)) {
  stop("Error: Seurat object file not found.")
}
seurat_obj <- readRDS(seurat_path)
message("Seurat object loaded successfully.")

# --- 3. EXPORT EXPRESSION MATRIX ---
message("Extracting the log-normalized expression matrix from the 'RNA' assay...")
expr_mat <- GetAssayData(seurat_obj, assay = "RNA", layer = "data")
write.csv(as.matrix(expr_mat), file = expr_matrix_outfile, row.names = TRUE)
message(paste("Expression matrix saved to:", expr_matrix_outfile))

# --- 4. PARSE METADATA AND EXPORT ---
message("Extracting and parsing cell metadata...")
cell_metadata <- seurat_obj@meta.data

# Parse identity to get base cell type and ERS state
cell_metadata <- cell_metadata %>%
  mutate(
    identity = as.character(identity), # Ensure it's a character vector
    ers_state = case_when(
      grepl("ERS\\+$", identity) ~ "ERS+",
      grepl("ERS-$", identity) ~ "ERS-",
      TRUE ~ "NA" # For cell types without an ERS state
    ),
    base_cell_type = trimws(gsub("ERS\\+?/?-?$", "", identity))
  )

message("Writing annotated metadata for ", nrow(cell_metadata), " cells.")
write.csv(cell_metadata, file = cell_meta_outfile, row.names = TRUE)
message(paste("Annotated cell metadata saved to:", cell_meta_outfile))

# --- 5. GENERATE STATISTICAL SUMMARY ---
message("Generating statistical summary file...")

# Use sink to capture print output of tables
sink(summary_file_path)

cat("=================================================================\n")
cat("       STATISTICAL SUMMARY: 02_prepare_expression_data.R\n")
cat("=================================================================\n\n")

cat(paste("Summary generated on:", Sys.Date()), "\n\n")

cat("--- Exported Data Dimensions ---\n")
cat(paste("Expression Matrix: ", nrow(expr_mat), " genes x ", ncol(expr_mat), " cells\n", sep=""))
cat(paste("Cell Metadata Table: ", nrow(cell_metadata), " cells x ", ncol(cell_metadata), " attributes\n\n", sep=""))

cat("--- Metadata Annotation Summary ---\n")
cat("A summary of the 'ers_state' column created from 'identity':\n")
print(table(cell_metadata$ers_state, useNA = "ifany"))
cat("\n")

cat("A summary of the 'base_cell_type' column created from 'identity':\n")
print(table(cell_metadata$base_cell_type, useNA = "ifany"))
cat("\n")

sink() # Stop capturing output

message(paste("Statistical summary saved to:", summary_file_path))
message("Script 02 finished successfully.")