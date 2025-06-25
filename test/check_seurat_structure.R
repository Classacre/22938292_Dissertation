# Diagnostic script for Seurat object structure
library(Seurat)

seurat_path <- "/group/sms029/Oliva_dataset/integrated_col_trajectories_colonly.rds"
seurat_obj <- readRDS(seurat_path)

cat("\n--- Structure of seurat_obj@assays$RNA@features ---\n")
str(seurat_obj@assays$RNA@features)
cat("Length:", length(seurat_obj@assays$RNA@features), "\n")

cat("\n--- Structure of seurat_obj@assays$RNA@layers ---\n")
if ("layers" %in% slotNames(seurat_obj@assays$RNA)) {
  print(names(seurat_obj@assays$RNA@layers))
  expr_mat <- as.matrix(seurat_obj@assays$RNA@layers$data)
  cat("Expression matrix (layers$data):\n")
} else {
  expr_mat <- as.matrix(seurat_obj@assays$RNA@data)
  cat("Expression matrix (@data):\n")
}
cat("Dimensions:", paste(dim(expr_mat), collapse = ", "), "\n")

cat("\nFirst 5 feature names:\n")
print(head(seurat_obj@assays$RNA@features, 5))

cat("\nFirst 5 rownames of expression matrix:\n")
print(head(rownames(expr_mat), 5))

cat("\nNumber of features:", length(seurat_obj@assays$RNA@features), "\n")
cat("Number of rows in expression matrix:", nrow(expr_mat), "\n")
