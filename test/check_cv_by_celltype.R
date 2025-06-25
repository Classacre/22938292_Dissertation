# Diagnostic script for cv_by_celltype structure
library(Seurat)
library(dplyr)
library(tidyr)

expr_mat <- as.matrix(read.csv("/group/sms029/mnieuwenh/gbm_noise_analysis/expression_matrix.csv", row.names=1, check.names=FALSE))
gbm_genes <- read.csv("/group/sms029/mnieuwenh/gbM_data/unique_gbm_genes.csv", stringsAsFactors=FALSE)$Gene
cell_metadata <- read.csv("/group/sms029/mnieuwenh/seurat_metadata/seurat_metadata_sample.csv", row.names=1, check.names=FALSE)

common_cells <- intersect(colnames(expr_mat), rownames(cell_metadata))
expr_mat <- expr_mat[, common_cells]
cell_metadata <- cell_metadata[common_cells, ]

cv_by_celltype <- lapply(unique(cell_metadata$identity), function(celltype) {
  cells_in_type <- rownames(cell_metadata)[cell_metadata$identity == celltype]
  if (length(cells_in_type) < 2) return(NULL)
  sub_mat <- expr_mat[, cells_in_type, drop=FALSE]
  mean_expr <- rowMeans(sub_mat)
  var_expr <- apply(sub_mat, 1, var)
  cv_expr <- sqrt(var_expr) / mean_expr
  data.frame(
    gene = rownames(sub_mat),
    celltype = celltype,
    gbm = rownames(sub_mat) %in% gbm_genes,
    cv_expr = cv_expr
  )
})

# Check structure
cat("cv_by_celltype is a list of length:", length(cv_by_celltype), "\n")
cat("First element structure:\n")
str(cv_by_celltype[[which(!sapply(cv_by_celltype, is.null))[1]]])

# Combine and check for non-finite values
cv_by_celltype_df <- bind_rows(cv_by_celltype)
cat("\nSummary of cv_expr column:\n")
print(summary(cv_by_celltype_df$cv_expr))
cat("\nNumber of non-finite values in cv_expr:", sum(!is.finite(cv_by_celltype_df$cv_expr)), "\n")
cat("\nFirst 10 non-finite values (if any):\n")
print(head(cv_by_celltype_df$cv_expr[!is.finite(cv_by_celltype_df$cv_expr)], 10))
