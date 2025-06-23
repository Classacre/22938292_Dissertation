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

# Save results
write.csv(results, "/group/sms029/mnieuwenh/gbm_noise_analysis/gbm_noise_comparison.csv", row.names=FALSE)

# Boxplot: CV by GBM status
p <- ggplot(results, aes(x=gbm, y=cv_expr, fill=gbm)) +
  geom_boxplot(outlier.size=0.5) +
  scale_x_discrete(labels=c("Non-GBM", "GBM")) +
  labs(title="Expression Noise (CV) by GBM Status", x="GBM Status", y="Coefficient of Variation (CV)") +
  theme_bw()
ggsave("/group/sms029/mnieuwenh/gbm_noise_analysis/gbm_noise_boxplot.png", plot=p, width=7, height=5)

# Wilcoxon test
stat <- wilcox.test(cv_expr ~ gbm, data=results)
cat("Wilcoxon test p-value:", stat$p.value, "\n", file="/group/sms029/mnieuwenh/gbm_noise_analysis/gbm_noise_stats.txt")
