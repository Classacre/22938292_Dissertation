# Load required libraries
library(Seurat)
library(dplyr)
if (!requireNamespace("knitr", quietly = TRUE)) install.packages("knitr")
library(knitr)

# Path to Seurat object
seurat_path <- "/group/sms029/Oliva_dataset/integrated_col_trajectories.rds"

# Load the Seurat object
seurat_obj <- readRDS(seurat_path)

# Save summary to a text file
sink("/group/sms029/mnieuwenh/seurat_metadata/seurat_summary.txt")

cat("=== Seurat Object Summary ===\n")
cat("Assays:\n")
print(Assays(seurat_obj))

cat("\nNumber of cells:", ncol(seurat_obj), "\n")
cat("Number of features (genes):", nrow(seurat_obj), "\n")

cat("\nMetadata columns:\n")
print(colnames(seurat_obj@meta.data))

cat("\nIdent levels:\n")
print(levels(seurat_obj))

cat("\nDimensional reductions:\n")
print(names(seurat_obj@reductions))

sink() # End writing to text file

# Save sample of metadata as CSV
write.csv(head(seurat_obj@meta.data, 5), "/group/sms029/mnieuwenh/seurat_metadata/seurat_metadata_sample.csv", row.names = TRUE)