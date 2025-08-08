# ==============================================================================
# SCRIPT: diagnose_seurat_object.R
#
# PURPOSE:
#   This is a standalone diagnostic script to inspect the contents of a Seurat
#   object. It is designed to be run non-interactively to help debug issues,
#   such as missing data slots or incorrect naming.
#
# IT PERFORMS THE FOLLOWING CHECKS:
#   1. Loads the Seurat object.
#   2. Prints the default Seurat object summary.
#   3. Explicitly lists all available dimensional reductions (e.g., pca, tsne, umap).
#   4. Lists all available assays (e.g., RNA, integrated).
#   5. Shows the first few rows of the metadata to check for column names.
#   6. Provides a summary of the 'identity' column.
# ==============================================================================

# --- 1. SETUP: Load libraries and define paths ---
suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(dplyr))

# INPUT: The path to the integrated Seurat object
seurat_path <- "/group/sms029/Oliva_dataset/integrated_col_trajectories_colonly.rds"

# --- 2. SCRIPT EXECUTION ---
message("--- Starting Seurat Object Diagnostic ---")

# Check if the file exists before attempting to load
if (!file.exists(seurat_path)) {
  stop(paste("FATAL ERROR: Seurat object file not found at:", seurat_path))
}

message(paste("\nLoading Seurat object from:", seurat_path))
seurat_obj <- readRDS(seurat_path)
message("Object loaded successfully.")

# --- 3. PERFORM DIAGNOSTIC CHECKS ---

message("\n\n--- 1. Default Seurat Object Printout ---")
# This printout often shows the active assay and available reductions.
print(seurat_obj)

message("\n\n--- 2. Available Dimensional Reductions ---")
# This is the most important check for the DimPlot error.
# It explicitly lists the names of all reduction objects.
reductions_available <- names(seurat_obj@reductions)
if (length(reductions_available) > 0) {
  message("The following dimensional reductions were found in the object:")
  print(reductions_available)
  message("\nNOTE: The 'reduction' parameter in DimPlot MUST be one of these names.")
} else {
  message("WARNING: No dimensional reductions were found in this Seurat object.")
}


message("\n\n--- 3. Available Assays ---")
assays_available <- names(seurat_obj@assays)
message("The following assays were found:")
print(assays_available)


message("\n\n--- 4. Metadata Head ---")
message("Showing the first 6 rows and first 10 columns of the metadata:")
# This helps verify that columns like 'identity' exist.
print(head(seurat_obj@meta.data[, 1:min(10, ncol(seurat_obj@meta.data))]))


message("\n\n--- 5. Summary of Cell Identities ---")
# This confirms the 'identity' column used for coloring the plot is present and valid.
if ("identity" %in% colnames(seurat_obj@meta.data)) {
  message("Counts for each level in the 'identity' column:")
  print(table(seurat_obj$identity))
} else {
  message("WARNING: The 'identity' column was not found in the metadata.")
}

message("\n\n--- Diagnostic Complete ---")