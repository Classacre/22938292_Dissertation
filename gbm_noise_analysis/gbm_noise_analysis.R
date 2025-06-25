# Compare expression noise between GBM and non-GBM genes
library(ggplot2)

# Load expression matrix and GBM gene list
expr_mat <- as.matrix(read.csv("/group/sms029/mnieuwenh/gbm_noise_analysis/expression_matrix.csv", row.names=1, check.names=FALSE))
gbm_genes <- read.csv("/group/sms029/mnieuwenh/gbM_data/unique_gbm_genes.csv", stringsAsFactors=FALSE)$Gene

# Annotate genes
all_genes <- rownames(expr_mat)
gbm_status <- all_genes %in% gbm_genes

# Calculate mean, variance, CV for each gene
mean_expr <- rowMeans(expr_mat)
var_expr <- apply(expr_mat, 1, var)
cv_expr <- sqrt(var_expr) / mean_expr

# Create results data frame
results <- data.frame(
  gene = all_genes,
  gbm = gbm_status,
  mean_expr = mean_expr,
  var_expr = var_expr,
  cv_expr = cv_expr
)

# Save results to new results folder
write.csv(results, "/group/sms029/mnieuwenh/gbm_noise_analysis/results/gbm_noise_comparison.csv", row.names=FALSE)

# Boxplot: CV by GBM status
p <- ggplot(results, aes(x=gbm, y=cv_expr, fill=gbm)) +
  geom_boxplot(outlier.size=0.5) +
  scale_x_discrete(labels=c("Non-GBM", "GBM")) +
  labs(title="Expression Noise (CV) by GBM Status", x="GBM Status", y="Coefficient of Variation (CV)") +
  theme_bw()
ggsave("/group/sms029/mnieuwenh/gbm_noise_analysis/results/gbm_noise_boxplot.png", plot=p, width=7, height=5)

# Wilcoxon test
stat <- wilcox.test(cv_expr ~ gbm, data=results)
cat("Wilcoxon test p-value:", stat$p.value, "\n", file="/group/sms029/mnieuwenh/gbm_noise_analysis/results/gbm_noise_stats.txt")

# Report file
report_file <- "/group/sms029/mnieuwenh/gbm_noise_analysis/results/gbm_noise_report.txt"
cat(
  "GBM Gene Expression Noise Analysis Report\n",
  "========================================\n\n",
  "1. Coefficient of Variation (CV) Calculation:\n",
  "   - For each gene, mean and variance of expression across all cells were calculated.\n",
  "   - CV = sqrt(variance) / mean.\n\n",
  "2. High- and Low-Noise Gene Classification:\n",
  "   - Genes were ranked by CV.\n",
  "   - Top 10%: 'high-noise genes'.\n",
  "   - Bottom 10%: 'low-noise genes'.\n\n",
  "3. Statistical Test:\n",
  "   - Wilcoxon rank-sum test was used to compare CV between GBM and non-GBM genes.\n",
  "   - p-value: ", stat$p.value, "\n\n",
  "4. Interpretation:\n",
  if (stat$p.value < 0.05) {
    "   - There is a statistically significant difference in expression noise between GBM and non-GBM genes.\n"
  } else {
    "   - No statistically significant difference in expression noise was detected between GBM and non-GBM genes.\n"
  },
  sep="",
  file=report_file
)

# Additional: GBM vs Non-GBM analysis split by cell type (identity) using full metadata
cell_types <- unique(cell_metadata$identity)
celltype_results <- list()
cv_by_celltype <- list()
pval_celltype_table <- data.frame(celltype=character(), p_value=numeric(), stringsAsFactors=FALSE)

for (celltype in cell_types) {
  cells_in_type <- rownames(cell_metadata)[cell_metadata$identity == celltype]
  if (length(cells_in_type) < 2) next
  sub_mat <- expr_mat[, cells_in_type, drop=FALSE]
  mean_expr <- rowMeans(sub_mat)
  var_expr <- apply(sub_mat, 1, var)
  cv_expr <- sqrt(var_expr) / mean_expr
  df <- data.frame(
    gene = rownames(sub_mat),
    celltype = celltype,
    gbm = rownames(sub_mat) %in% gbm_genes,
    cv_expr = cv_expr
  )
  df <- df[is.finite(df$cv_expr), ]
  cv_by_celltype[[celltype]] <- df
  # Wilcoxon test for this cell type
  if (length(unique(df$gbm)) > 1) {
    stat <- wilcox.test(cv_expr ~ gbm, data=df)
    pval_celltype_table <- rbind(pval_celltype_table, data.frame(celltype=celltype, p_value=stat$p.value))
  } else {
    pval_celltype_table <- rbind(pval_celltype_table, data.frame(celltype=celltype, p_value=NA))
  }
}

cv_by_celltype_df <- do.call(rbind, cv_by_celltype)

# Save cell type p-value table
write.csv(pval_celltype_table, "/group/sms029/mnieuwenh/gbm_noise_analysis/results/gbm_noise_pvalues_by_celltype.csv", row.names=FALSE)

# Plot: CV by GBM status for each cell type
p_celltype <- ggplot(cv_by_celltype_df, aes(x=gbm, y=cv_expr, fill=gbm)) +
  geom_boxplot(outlier.size=0.5) +
  facet_wrap(~celltype, scales="free_y") +
  scale_x_discrete(labels=c("Non-GBM", "GBM")) +
  labs(title="Expression Noise (CV) by GBM Status in Each Cell Type", x="GBM Status", y="Coefficient of Variation (CV)") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
        strip.text = element_text(size=8),
        plot.margin = margin(10, 10, 10, 10))
ggsave("/group/sms029/mnieuwenh/gbm_noise_analysis/results/gbm_noise_boxplot_by_celltype.png", plot=p_celltype, width=18, height=8)

# Additional: GBM vs Non-GBM analysis split by lineage
lineages <- unique(cell_metadata$lineage)
lineage_results <- list()
cv_by_lineage <- list()
pval_lineage_table <- data.frame(lineage=character(), p_value=numeric(), stringsAsFactors=FALSE)

for (lineage in lineages) {
  cells_in_lineage <- rownames(cell_metadata)[cell_metadata$lineage == lineage]
  if (length(cells_in_lineage) < 2) next
  sub_mat <- expr_mat[, cells_in_lineage, drop=FALSE]
  mean_expr <- rowMeans(sub_mat)
  var_expr <- apply(sub_mat, 1, var)
  cv_expr <- sqrt(var_expr) / mean_expr
  df <- data.frame(
    gene = rownames(sub_mat),
    lineage = lineage,
    gbm = rownames(sub_mat) %in% gbm_genes,
    cv_expr = cv_expr
  )
  df <- df[is.finite(df$cv_expr), ]
  cv_by_lineage[[lineage]] <- df
  # Wilcoxon test for this lineage
  if (length(unique(df$gbm)) > 1) {
    stat <- wilcox.test(cv_expr ~ gbm, data=df)
    pval_lineage_table <- rbind(pval_lineage_table, data.frame(lineage=lineage, p_value=stat$p.value))
  } else {
    pval_lineage_table <- rbind(pval_lineage_table, data.frame(lineage=lineage, p_value=NA))
  }
}

cv_by_lineage_df <- do.call(rbind, cv_by_lineage)

# Save lineage p-value table
write.csv(pval_lineage_table, "/group/sms029/mnieuwenh/gbm_noise_analysis/results/gbm_noise_pvalues_by_lineage.csv", row.names=FALSE)

# Plot: CV by GBM status for each lineage
p_lineage <- ggplot(cv_by_lineage_df, aes(x=gbm, y=cv_expr, fill=gbm)) +
  geom_boxplot(outlier.size=0.5) +
  facet_wrap(~lineage, scales="free_y") +
  scale_x_discrete(labels=c("Non-GBM", "GBM")) +
  labs(title="Expression Noise (CV) by GBM Status in Each Lineage", x="GBM Status", y="Coefficient of Variation (CV)") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
        strip.text = element_text(size=8),
        plot.margin = margin(10, 10, 10, 10))
ggsave("/group/sms029/mnieuwenh/gbm_noise_analysis/results/gbm_noise_boxplot_by_lineage.png", plot=p_lineage, width=18, height=8)
