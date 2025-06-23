# Export gene expression matrix from Seurat object
library(Seurat)

seurat_path <- "/group/sms029/Oliva_dataset/integrated_col_trajectories.rds"
outfile <- "/group/sms029/mnieuwenh/gbm_noise_analysis/expression_matrix.csv"

seurat_obj <- readRDS(seurat_path)
expr_mat <- as.matrix(seurat_obj@assays$RNA@data) # log-normalized counts
gene_names <- rownames(expr_mat)
cell_names <- colnames(expr_mat)

# Write as CSV (genes as rows, cells as columns)
write.csv(expr_mat, file=outfile, row.names=TRUE)
cat("Expression matrix saved to:", outfile, "\n")
