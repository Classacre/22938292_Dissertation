# Inspect the expression matrix and GBM gene overlap
expr_mat <- as.matrix(read.csv("/group/sms029/mnieuwenh/gbm_noise_analysis/expression_matrix.csv", row.names=1, check.names=FALSE))
gbm_genes <- read.csv("/group/sms029/mnieuwenh/gbM_data/unique_gbm_genes.csv", stringsAsFactors=FALSE)$Gene

cat("Expression matrix dimensions: ", dim(expr_mat)[1], "genes x", dim(expr_mat)[2], "cells\n")
cat("First 5 gene names in matrix:\n")
print(head(rownames(expr_mat), 5))
cat("Number of unique GBM genes in list: ", length(unique(gbm_genes)), "\n")
cat("Number of genes in matrix that are GBM: ", sum(rownames(expr_mat) %in% gbm_genes), "\n")
cat("Number of genes in matrix that are non-GBM: ", sum(!(rownames(expr_mat) %in% gbm_genes)), "\n")
cat("Are there both GBM and non-GBM genes in the matrix? ", 
    (sum(rownames(expr_mat) %in% gbm_genes) > 0) & (sum(!(rownames(expr_mat) %in% gbm_genes)) > 0), "\n")
