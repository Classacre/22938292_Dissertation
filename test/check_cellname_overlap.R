# Check cell name matching between expression matrix and metadata
expr_mat <- as.matrix(read.csv("/group/sms029/mnieuwenh/gbm_noise_analysis/expression_matrix.csv", row.names=1, check.names=FALSE))
cell_metadata <- read.csv("/group/sms029/mnieuwenh/seurat_metadata/seurat_metadata_sample.csv", row.names=1, check.names=FALSE)

cat("First 5 colnames in expression matrix:\n")
print(head(colnames(expr_mat), 5))
cat("\nFirst 5 rownames in metadata:\n")
print(head(rownames(cell_metadata), 5))

cat("\nNumber of overlapping cell names:\n")
print(length(intersect(colnames(expr_mat), rownames(cell_metadata))))

cat("\nExample of non-overlapping cell names in expr_mat:\n")
print(setdiff(colnames(expr_mat), rownames(cell_metadata))[1:5])
cat("\nExample of non-overlapping cell names in metadata:\n")
print(setdiff(rownames(cell_metadata), colnames(expr_mat))[1:5])
