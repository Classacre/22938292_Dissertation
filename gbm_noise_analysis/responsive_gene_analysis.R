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

# Load methylation annotation and merge
meth_anno <- read.csv("/group/sms029/mnieuwenh/Methylation_Data/combined_methylation_data.csv", stringsAsFactors=FALSE)
responsive_genes <- merge(responsive_genes, meth_anno, by.x="gene", by.y="Gene_ID", all.x=TRUE)

# Summarize: proportion of responsive genes by methylation group (Cahn)
table_cahn <- table(responsive_genes$is_responsive, responsive_genes$Cahn_Methylation_status)
prop_cahn <- prop.table(table_cahn, margin=2)
write.csv(prop_cahn, "/group/sms029/mnieuwenh/gbm_noise_analysis/results/responsive_genes/prop_responsive_by_cahn.csv")

# Summarize: proportion of responsive genes by methylation group (Bewick)
table_bewick <- table(responsive_genes$is_responsive, responsive_genes$Bewick_Classification)
prop_bewick <- prop.table(table_bewick, margin=2)
write.csv(prop_bewick, "/group/sms029/mnieuwenh/gbm_noise_analysis/results/responsive_genes/prop_responsive_by_bewick.csv")

# Boxplot: CV for responsive/non-responsive genes within each methylation group (Cahn)
p_cahn <- ggplot(responsive_genes, aes(x=Cahn_Methylation_status, y=cv_expr, fill=is_responsive)) +
  geom_boxplot(outlier.size=0.5, outlier.shape=16, outlier.alpha=0.3, position=position_dodge(width=0.8)) +
  scale_fill_manual(values=c("#999999", "#E69F00"), labels=c("Non-Responsive", "Responsive")) +
  labs(title="CV by Cahn Methylation and Responsiveness", x="Cahn Methylation", y="CV") +
  theme_bw()
ggsave("/group/sms029/mnieuwenh/gbm_noise_analysis/results/responsive_genes/cv_by_cahn_and_responsive.png", plot=p_cahn, width=8, height=6)

# Boxplot: CV for responsive/non-responsive genes within each methylation group (Bewick)
p_bewick <- ggplot(responsive_genes, aes(x=Bewick_Classification, y=cv_expr, fill=is_responsive)) +
  geom_boxplot(outlier.size=0.5, outlier.shape=16, outlier.alpha=0.3, position=position_dodge(width=0.8)) +
  scale_fill_manual(values=c("#999999", "#E69F00"), labels=c("Non-Responsive", "Responsive")) +
  labs(title="CV by Bewick Methylation and Responsiveness", x="Bewick Methylation", y="CV") +
  theme_bw()
ggsave("/group/sms029/mnieuwenh/gbm_noise_analysis/results/responsive_genes/cv_by_bewick_and_responsive.png", plot=p_bewick, width=8, height=6)

# Fisher's exact test: responsive vs non-responsive by methylation group (Cahn)
fisher_cahn <- fisher.test(table_cahn)
cat("Fisher's exact test p-value (Cahn):", fisher_cahn$p.value, "\n", file="/group/sms029/mnieuwenh/gbm_noise_analysis/results/responsive_genes/fisher_cahn.txt")

# Fisher's exact test: responsive vs non-responsive by methylation group (Bewick)
fisher_bewick <- fisher.test(table_bewick)
cat("Fisher's exact test p-value (Bewick):", fisher_bewick$p.value, "\n", file="/group/sms029/mnieuwenh/gbm_noise_analysis/results/responsive_genes/fisher_bewick.txt")

# Add H2A.Z annotation columns
responsive_genes$H2AZ_Depleted <- responsive_genes$H2AZ_Depleted == TRUE
responsive_genes$H2AZ_Enriched <- responsive_genes$H2AZ_Enriched == TRUE
responsive_genes$H2AZ_Group <- NA
responsive_genes$H2AZ_Group[responsive_genes$H2AZ_Depleted] <- "H2A.Z-Depleted"
responsive_genes$H2AZ_Group[responsive_genes$H2AZ_Enriched] <- "H2A.Z-Enriched"
responsive_genes$H2AZ_Group[is.na(responsive_genes$H2AZ_Group)] <- "Other/NA"
responsive_genes$H2AZ_Group <- factor(responsive_genes$H2AZ_Group, levels=c("H2A.Z-Depleted", "H2A.Z-Enriched", "Other/NA"))

# Proportion of responsive genes by H2A.Z group
prop_h2az <- prop.table(table(responsive_genes$is_responsive, responsive_genes$H2AZ_Group), margin=2)
write.csv(prop_h2az, "/group/sms029/mnieuwenh/gbm_noise_analysis/results/responsive_genes/prop_responsive_by_h2az.csv")

# Boxplot: CV for responsive/non-responsive genes within each H2A.Z group
p_h2az <- ggplot(responsive_genes, aes(x=H2AZ_Group, y=cv_expr, fill=is_responsive)) +
  geom_boxplot(outlier.size=0.5, outlier.shape=16, outlier.alpha=0.3, position=position_dodge(width=0.8)) +
  scale_fill_manual(values=c("#999999", "#E69F00"), labels=c("Non-Responsive", "Responsive")) +
  labs(title="CV by H2A.Z Group and Responsiveness", x="H2A.Z Group", y="CV") +
  theme_bw()
ggsave("/group/sms029/mnieuwenh/gbm_noise_analysis/results/responsive_genes/cv_by_h2az_and_responsive.png", plot=p_h2az, width=8, height=6)

# Fisher's exact test: responsive vs non-responsive by H2A.Z group
fisher_h2az <- fisher.test(table(responsive_genes$is_responsive, responsive_genes$H2AZ_Group))
cat("Fisher's exact test p-value (H2A.Z):", fisher_h2az$p.value, "\n", file="/group/sms029/mnieuwenh/gbm_noise_analysis/results/responsive_genes/fisher_h2az.txt")

# Combined methylation + H2A.Z group (optional, for interaction)
responsive_genes$meth_h2az_group <- paste0(responsive_genes$Cahn_Methylation_status, "+", responsive_genes$H2AZ_Group)

# Boxplot: CV by methylation + H2A.Z group
p_meth_h2az <- ggplot(responsive_genes, aes(x=meth_h2az_group, y=cv_expr, fill=is_responsive)) +
  geom_boxplot(outlier.size=0.5, outlier.shape=16, outlier.alpha=0.3, position=position_dodge(width=0.8)) +
  scale_fill_manual(values=c("#999999", "#E69F00"), labels=c("Non-Responsive", "Responsive")) +
  labs(title="CV by Cahn Methylation + H2A.Z Group and Responsiveness", x="Cahn Methylation + H2A.Z Group", y="CV") +
  theme_bw() +
  theme(axis.text.x = element_text(angle=45, hjust=1))
ggsave("/group/sms029/mnieuwenh/gbm_noise_analysis/results/responsive_genes/cv_by_meth_h2az_and_responsive.png", plot=p_meth_h2az, width=12, height=6)
