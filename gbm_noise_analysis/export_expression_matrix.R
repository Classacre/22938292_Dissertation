# Export gene expression matrix from Seurat object
library(Seurat)

seurat_path <- "/group/sms029/Oliva_dataset/integrated_col_trajectories_colonly.rds"
outfile <- "/group/sms029/mnieuwenh/gbm_noise_analysis/expression_matrix.csv"

seurat_obj <- readRDS(seurat_path)
# Use the new Seurat v5 Assay5 structure: data is in layers
expr_mat <- as.matrix(seurat_obj@assays$RNA@layers$data) # log-normalized counts

# Get gene names from the rownames of the LogMap features object
feature_names <- rownames(seurat_obj@assays$RNA@features)
if (length(feature_names) == nrow(expr_mat)) {
  rownames(expr_mat) <- feature_names
} else {
  stop("Number of features does not match number of rows in expression matrix. Cannot assign gene names.")
}

# Get cell barcodes from the Seurat object and assign as column names
cell_names <- colnames(seurat_obj)
if (length(cell_names) == ncol(expr_mat)) {
  colnames(expr_mat) <- cell_names
} else {
  stop("Number of cells does not match number of columns in expression matrix. Cannot assign cell barcodes.")
}

# Write as CSV (genes as rows, cells as columns)
write.csv(expr_mat, file=outfile, row.names=TRUE)
cat("Expression matrix saved to:", outfile, "\n")
