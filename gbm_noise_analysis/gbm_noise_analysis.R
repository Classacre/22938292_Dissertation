# Compare expression noise between gbM and TE-like methylation genes (Cahn and Bewick)
library(ggplot2)

# Load expression matrix
expr_mat <- as.matrix(read.csv("/group/sms029/mnieuwenh/gbm_noise_analysis/expression_matrix.csv", row.names=1, check.names=FALSE))

# Load combined methylation annotation
gene_anno <- read.csv("/group/sms029/mnieuwenh/Methylation_Data/combined_methylation_data.csv", stringsAsFactors=FALSE)
rownames(gene_anno) <- gene_anno$Gene_ID

# Annotate genes for each methylation status
genes <- rownames(expr_mat)
anno <- gene_anno[genes, ]

cahn_gbm <- !is.na(anno$Cahn_Methylation_status) & anno$Cahn_Methylation_status == "gbM"
bewick_gbm <- !is.na(anno$Bewick_Classification) & anno$Bewick_Classification == "gbM"
bewick_te_like <- !is.na(anno$Bewick_Classification) & (anno$Bewick_Classification == "mCHG" | anno$Bewick_Classification == "mCHH")

# Calculate mean, variance, CV for each gene
mean_expr <- rowMeans(expr_mat)
var_expr <- apply(expr_mat, 1, var)
cv_expr <- sqrt(var_expr) / mean_expr

# Create results data frame
results <- data.frame(
  gene = genes,
  cahn_gbm = cahn_gbm,
  bewick_gbm = bewick_gbm,
  bewick_te_like = bewick_te_like,
  mean_expr = mean_expr,
  var_expr = var_expr,
  cv_expr = cv_expr
)

write.csv(results, "/group/sms029/mnieuwenh/gbm_noise_analysis/results/high_low_noise/gbm_noise_comparison.csv", row.names=FALSE)

# Boxplots: CV by each methylation status
gbm_labels <- c("FALSE"="Non-gbM", "TRUE"="gbM")
te_labels <- c("FALSE"="Other", "TRUE"="TE-like")

p1 <- ggplot(results, aes(x=factor(cahn_gbm), y=cv_expr, fill=factor(cahn_gbm))) +
  geom_boxplot(outlier.size=0.5) +
  scale_x_discrete(labels=gbm_labels) +
  labs(title="Expression Noise (CV) by Cahn gbM Status", x="Cahn gbM", y="CV") +
  theme_bw()
ggsave("/group/sms029/mnieuwenh/gbm_noise_analysis/results/high_low_noise/gbm_noise_boxplot_cahn.png", plot=p1, width=7, height=5)

p2 <- ggplot(results, aes(x=factor(bewick_gbm), y=cv_expr, fill=factor(bewick_gbm))) +
  geom_boxplot(outlier.size=0.5) +
  scale_x_discrete(labels=gbm_labels) +
  labs(title="Expression Noise (CV) by Bewick gbM Status", x="Bewick gbM", y="CV") +
  theme_bw()
ggsave("/group/sms029/mnieuwenh/gbm_noise_analysis/results/high_low_noise/gbm_noise_boxplot_bewick.png", plot=p2, width=7, height=5)

p3 <- ggplot(results, aes(x=factor(bewick_te_like), y=cv_expr, fill=factor(bewick_te_like))) +
  geom_boxplot(outlier.size=0.5) +
  scale_x_discrete(labels=te_labels) +
  labs(title="Expression Noise (CV) by Bewick TE-like Status", x="Bewick TE-like", y="CV") +
  theme_bw()
ggsave("/group/sms029/mnieuwenh/gbm_noise_analysis/results/high_low_noise/gbm_noise_boxplot_te_like.png", plot=p3, width=7, height=5)

# Assign each gene to a methylation group for combined boxplot
results$methylation_group <- NA
results$methylation_group[!is.na(anno$Cahn_Methylation_status) & anno$Cahn_Methylation_status == "gbM"] <- "Cahn gbM"
results$methylation_group[!is.na(anno$Bewick_Classification) & anno$Bewick_Classification == "gbM"] <- "Bewick gbM"
results$methylation_group[!is.na(anno$Bewick_Classification) & anno$Bewick_Classification == "mCHH"] <- "CHH"
results$methylation_group[!is.na(anno$Bewick_Classification) & anno$Bewick_Classification == "mCHG"] <- "CHG"
results$methylation_group[!is.na(anno$Cahn_Methylation_status) & anno$Cahn_Methylation_status == "TE-like methylation"] <- "TE-Like Methylation"
results$methylation_group[!is.na(anno$Cahn_Methylation_status) & anno$Cahn_Methylation_status == "Unmethylated" & is.na(results$methylation_group)] <- "Unmethylated"
results$methylation_group[is.na(results$methylation_group)] <- "Non-gbM"

# Combined boxplot
results$methylation_group <- factor(results$methylation_group, levels=c("Unmethylated", "Non-gbM", "Bewick gbM", "Cahn gbM", "CHH", "CHG", "TE-Like Methylation"))
p_combined <- ggplot(results, aes(x=methylation_group, y=cv_expr, color=methylation_group)) +
  geom_jitter(width=0.25, alpha=0.5, size=0.7) +
  labs(title="Univariant Scatter Plot: Expression Noise (CV) by Methylation Group", x="Methylation Group", y="CV") +
  theme_bw() +
  theme(axis.text.x = element_text(angle=45, hjust=1))
ggsave("/group/sms029/mnieuwenh/gbm_noise_analysis/results/high_low_noise/gbm_noise_scatter_combined.png", plot=p_combined, width=10, height=6)

# Assign methylation group for Cahn
results$cahn_group <- NA
results$cahn_group[!is.na(anno$Cahn_Methylation_status) & anno$Cahn_Methylation_status == "gbM"] <- "gbM"
results$cahn_group[!is.na(anno$Cahn_Methylation_status) & anno$Cahn_Methylation_status == "TE-like methylation"] <- "TE-like"
results$cahn_group[!is.na(anno$Cahn_Methylation_status) & anno$Cahn_Methylation_status == "Unmethylated"] <- "Unmethylated"
results$cahn_group <- factor(results$cahn_group, levels=c("Unmethylated", "gbM", "TE-like"))

# Assign methylation group for Bewick
results$bewick_group <- NA
results$bewick_group[!is.na(anno$Bewick_Classification) & anno$Bewick_Classification == "gbM"] <- "gbM"
results$bewick_group[!is.na(anno$Bewick_Classification) & anno$Bewick_Classification == "mCHH"] <- "CHH"
results$bewick_group[!is.na(anno$Bewick_Classification) & anno$Bewick_Classification == "mCHG"] <- "CHG"
results$bewick_group[!is.na(anno$Bewick_Classification) & anno$Bewick_Classification == "Unmethylated"] <- "Unmethylated"
results$bewick_group <- factor(results$bewick_group, levels=c("Unmethylated", "gbM", "CHH", "CHG"))

# Add sd_expr column before filtering for plotting
results$sd_expr <- sqrt(results$var_expr)

# Filter out non-finite CV values for plotting
plot_results <- results[is.finite(results$cv_expr), ]

# Remove extreme outliers (e.g., above 99th percentile) for plotting
cv_cap <- quantile(plot_results$cv_expr, 0.99, na.rm=TRUE)
plot_results <- plot_results[plot_results$cv_expr <= cv_cap, ]

# Cahn boxplot
p_cahn <- ggplot(plot_results[!is.na(plot_results$cahn_group),], aes(x=cahn_group, y=cv_expr, fill=cahn_group)) +
  geom_boxplot(outlier.size=0.5, outlier.shape=16, outlier.alpha=0.3) +
  labs(title="Expression Noise (CV) by Cahn Methylation Group", x="Cahn Methylation Group", y="CV") +
  theme_bw() +
  theme(axis.text.x = element_text(angle=45, hjust=1))
ggsave("/group/sms029/mnieuwenh/gbm_noise_analysis/results/high_low_noise/gbm_noise_boxplot_cahn.png", plot=p_cahn, width=8, height=6)

# Bewick boxplot
p_bewick <- ggplot(plot_results[!is.na(plot_results$bewick_group),], aes(x=bewick_group, y=cv_expr, fill=bewick_group)) +
  geom_boxplot(outlier.size=0.5, outlier.shape=16, outlier.alpha=0.3) +
  labs(title="Expression Noise (CV) by Bewick Methylation Group", x="Bewick Methylation Group", y="CV") +
  theme_bw() +
  theme(axis.text.x = element_text(angle=45, hjust=1))
ggsave("/group/sms029/mnieuwenh/gbm_noise_analysis/results/high_low_noise/gbm_noise_boxplot_bewick.png", plot=p_bewick, width=8, height=6)

# Additional scatter plots: mean vs sd and mean vs CV, split by Cahn and Bewick methylation groups

# Cahn: Mean vs SD
p_mean_sd_cahn <- ggplot(plot_results[!is.na(plot_results$cahn_group),], aes(x=mean_expr, y=sd_expr, color=cahn_group)) +
  geom_point(alpha=0.5, size=0.7) +
  labs(title="Mean Expression vs SD by Cahn Methylation Group", x="Mean Expression", y="Standard Deviation") +
  theme_bw()
ggsave("/group/sms029/mnieuwenh/gbm_noise_analysis/results/high_low_noise/mean_vs_sd_by_cahn_group.png", plot=p_mean_sd_cahn, width=8, height=6)

# Cahn: Mean vs CV
p_mean_cv_cahn <- ggplot(plot_results[!is.na(plot_results$cahn_group),], aes(x=mean_expr, y=cv_expr, color=cahn_group)) +
  geom_point(alpha=0.5, size=0.7) +
  labs(title="Mean Expression vs CV by Cahn Methylation Group", x="Mean Expression", y="Coefficient of Variation (CV)") +
  theme_bw()
ggsave("/group/sms029/mnieuwenh/gbm_noise_analysis/results/high_low_noise/mean_vs_cv_by_cahn_group.png", plot=p_mean_cv_cahn, width=8, height=6)

# Bewick: Mean vs SD
p_mean_sd_bewick <- ggplot(plot_results[!is.na(plot_results$bewick_group),], aes(x=mean_expr, y=sd_expr, color=bewick_group)) +
  geom_point(alpha=0.5, size=0.7) +
  labs(title="Mean Expression vs SD by Bewick Methylation Group", x="Mean Expression", y="Standard Deviation") +
  theme_bw()
ggsave("/group/sms029/mnieuwenh/gbm_noise_analysis/results/high_low_noise/mean_vs_sd_by_bewick_group.png", plot=p_mean_sd_bewick, width=8, height=6)

# Bewick: Mean vs CV
p_mean_cv_bewick <- ggplot(plot_results[!is.na(plot_results$bewick_group),], aes(x=mean_expr, y=cv_expr, color=bewick_group)) +
  geom_point(alpha=0.5, size=0.7) +
  labs(title="Mean Expression vs CV by Bewick Methylation Group", x="Mean Expression", y="Coefficient of Variation (CV)") +
  theme_bw()
ggsave("/group/sms029/mnieuwenh/gbm_noise_analysis/results/high_low_noise/mean_vs_cv_by_bewick_group.png", plot=p_mean_cv_bewick, width=8, height=6)

# Annotate genes for each methylation status and H2A.Z status
h2az_depleted <- !is.na(anno$H2AZ_Depleted) & anno$H2AZ_Depleted == TRUE
h2az_enriched <- !is.na(anno$H2AZ_Enriched) & anno$H2AZ_Enriched == TRUE
h2az_group <- rep(NA, length(genes))
h2az_group[h2az_depleted] <- "H2A.Z-Depleted"
h2az_group[h2az_enriched] <- "H2A.Z-Enriched"
h2az_group[is.na(h2az_group)] <- "Other/NA"

# Add H2A.Z columns to results
gbm_labels <- c("FALSE"="Non-gbM", "TRUE"="gbM")
te_labels <- c("FALSE"="Other", "TRUE"="TE-like")

results$h2az_depleted <- h2az_depleted
results$h2az_enriched <- h2az_enriched
results$h2az_group <- factor(h2az_group, levels=c("H2A.Z-Depleted", "H2A.Z-Enriched", "Other/NA"))

# Boxplot: CV by H2A.Z group
p_h2az <- ggplot(results, aes(x=h2az_group, y=cv_expr, fill=h2az_group)) +
  geom_boxplot(outlier.size=0.5) +
  labs(title="Expression Noise (CV) by H2A.Z Group", x="H2A.Z Group", y="CV") +
  theme_bw()
ggsave("/group/sms029/mnieuwenh/gbm_noise_analysis/results/high_low_noise/gbm_noise_boxplot_h2az.png", plot=p_h2az, width=8, height=6)

# Scatter plots: mean vs SD and mean vs CV by H2A.Z group
p_mean_sd_h2az <- ggplot(results, aes(x=mean_expr, y=sd_expr, color=h2az_group)) +
  geom_point(alpha=0.5, size=0.7) +
  labs(title="Mean Expression vs SD by H2A.Z Group", x="Mean Expression", y="Standard Deviation") +
  theme_bw()
ggsave("/group/sms029/mnieuwenh/gbm_noise_analysis/results/high_low_noise/mean_vs_sd_by_h2az_group.png", plot=p_mean_sd_h2az, width=8, height=6)

p_mean_cv_h2az <- ggplot(results, aes(x=mean_expr, y=cv_expr, color=h2az_group)) +
  geom_point(alpha=0.5, size=0.7) +
  labs(title="Mean Expression vs CV by H2A.Z Group", x="Mean Expression", y="Coefficient of Variation (CV)") +
  theme_bw()
ggsave("/group/sms029/mnieuwenh/gbm_noise_analysis/results/high_low_noise/mean_vs_cv_by_h2az_group.png", plot=p_mean_cv_h2az, width=8, height=6)

# Combined methylation + H2A.Z group (optional, for interaction)
results$meth_h2az_group <- paste0(results$methylation_group, "+", results$h2az_group)

# Boxplot: CV by methylation + H2A.Z group
p_meth_h2az <- ggplot(results, aes(x=meth_h2az_group, y=cv_expr, fill=meth_h2az_group)) +
  geom_boxplot(outlier.size=0.5) +
  labs(title="Expression Noise (CV) by Methylation + H2A.Z Group", x="Methylation + H2A.Z Group", y="CV") +
  theme_bw() +
  theme(axis.text.x = element_text(angle=45, hjust=1))
ggsave("/group/sms029/mnieuwenh/gbm_noise_analysis/results/high_low_noise/gbm_noise_boxplot_meth_h2az.png", plot=p_meth_h2az, width=12, height=6)

# Wilcoxon tests
stat_cahn <- wilcox.test(cv_expr ~ cahn_gbm, data=results)
stat_bewick <- wilcox.test(cv_expr ~ bewick_gbm, data=results)
stat_te_like <- wilcox.test(cv_expr ~ bewick_te_like, data=results)
stat_h2az <- wilcox.test(cv_expr ~ h2az_group, data=results, subset=h2az_group %in% c("H2A.Z-Depleted", "H2A.Z-Enriched"))

cat("Wilcoxon test p-value (Cahn gbM):", stat_cahn$p.value, "\n", file="/group/sms029/mnieuwenh/gbm_noise_analysis/results/high_low_noise/gbm_noise_stats.txt")
cat("Wilcoxon test p-value (Bewick gbM):", stat_bewick$p.value, "\n", file="/group/sms029/mnieuwenh/gbm_noise_analysis/results/high_low_noise/gbm_noise_stats.txt", append=TRUE)
cat("Wilcoxon test p-value (Bewick TE-like):", stat_te_like$p.value, "\n", file="/group/sms029/mnieuwenh/gbm_noise_analysis/results/high_low_noise/gbm_noise_stats.txt", append=TRUE)
cat("Wilcoxon test p-value (H2A.Z-Depleted vs H2A.Z-Enriched):", stat_h2az$p.value, "\n", file="/group/sms029/mnieuwenh/gbm_noise_analysis/results/high_low_noise/gbm_noise_stats.txt", append=TRUE)
