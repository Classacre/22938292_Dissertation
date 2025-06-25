# Identify responsive genes and characterize their noise
# Responsive genes: highly expressed in one cell type (identity), low in others
# Output: responsive gene list and their noise characteristics

library(dplyr)
library(tidyr)
library(ggplot2)

# Load expression matrix and metadata
expr_mat <- as.matrix(read.csv("/group/sms029/mnieuwenh/gbm_noise_analysis/expression_matrix.csv", row.names=1, check.names=FALSE))
cell_metadata <- read.csv("/group/sms029/mnieuwenh/seurat_metadata/seurat_metadata_full.csv", row.names=1, check.names=FALSE)

# Ensure columns match between expr_mat and metadata
common_cells <- intersect(colnames(expr_mat), rownames(cell_metadata))
expr_mat <- expr_mat[, common_cells]
cell_metadata <- cell_metadata[common_cells, ]

# Calculate average expression per gene per cell type
cell_types <- unique(cell_metadata$identity)
gene_means <- sapply(cell_types, function(ct) {
  cells <- rownames(cell_metadata)[cell_metadata$identity == ct]
  if (length(cells) < 1) return(rep(NA, nrow(expr_mat)))
  rowMeans(expr_mat[, cells, drop=FALSE])
})
rownames(gene_means) <- rownames(expr_mat)

# For each gene, find the cell type with max mean expression and the second highest
max_ct <- apply(gene_means, 1, function(x) names(which.max(x)))
max_val <- apply(gene_means, 1, max, na.rm=TRUE)
second_val <- apply(gene_means, 1, function(x) sort(x, decreasing=TRUE)[2])

# Define responsive genes: max mean > 2x second highest mean and max mean > 1 (adjust threshold as needed)
responsive <- (max_val > 2 * second_val) & (max_val > 1)
responsive_genes <- data.frame(
  gene = rownames(expr_mat),
  is_responsive = responsive,
  max_celltype = max_ct,
  max_mean = max_val,
  second_mean = second_val
)

# Save responsive gene list
write.csv(responsive_genes, "/group/sms029/mnieuwenh/gbm_noise_analysis/results/responsive_genes/responsive_genes.csv", row.names=FALSE)

# Characterize noise for responsive vs non-responsive genes
# Calculate CV for each gene (across all cells)
mean_expr <- rowMeans(expr_mat)
var_expr <- apply(expr_mat, 1, var)
cv_expr <- sqrt(var_expr) / mean_expr
responsive_genes$cv_expr <- cv_expr

# Boxplot: CV for responsive vs non-responsive genes
p_resp <- ggplot(responsive_genes, aes(x=is_responsive, y=cv_expr, fill=is_responsive)) +
  geom_boxplot(outlier.size=0.5) +
  scale_x_discrete(labels=c("Non-Responsive", "Responsive")) +
  labs(title="Expression Noise (CV) for Responsive vs Non-Responsive Genes", x="Gene Type", y="Coefficient of Variation (CV)") +
  theme_bw()
ggsave("/group/sms029/mnieuwenh/gbm_noise_analysis/results/responsive_genes/responsive_vs_nonresponsive_noise.png", plot=p_resp, width=8, height=6)

# Wilcoxon test
stat_resp <- wilcox.test(cv_expr ~ is_responsive, data=responsive_genes)
cat("Wilcoxon test p-value:", stat_resp$p.value, "\n", file="/group/sms029/mnieuwenh/gbm_noise_analysis/results/responsive_genes/responsive_noise_stats.txt")
