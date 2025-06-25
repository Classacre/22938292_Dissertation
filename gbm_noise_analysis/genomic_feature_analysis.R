# Genomic feature analysis for high/low noise and responsive genes
# - Identify features (e.g., gbM status) of high/low noise and high/low responsiveness genes
# - Determine shared features between noisy and responsive genes

library(dplyr)

# Load results from previous steps
noise_results <- read.csv("/group/sms029/mnieuwenh/gbm_noise_analysis/results/high_low_noise/gbm_noise_comparison.csv")
responsive_results <- read.csv("/group/sms029/mnieuwenh/gbm_noise_analysis/results/responsive_genes/responsive_genes.csv")

# Classify high- and low-noise genes (top/bottom 10% by CV)
noise_results <- noise_results %>% arrange(desc(cv_expr))
n <- nrow(noise_results)
high_noise_genes <- noise_results$gene[1:ceiling(0.1*n)]
low_noise_genes <- noise_results$gene[(n-ceiling(0.1*n)+1):n]

# Classify high- and low-responsiveness genes (top/bottom 10% by max_mean - second_mean)
responsive_results <- responsive_results %>% mutate(resp_score = max_mean - second_mean)
responsive_results <- responsive_results %>% arrange(desc(resp_score))
nr <- nrow(responsive_results)
high_resp_genes <- responsive_results$gene[1:ceiling(0.1*nr)]
low_resp_genes <- responsive_results$gene[(nr-ceiling(0.1*nr)+1):nr]

# Annotate genes with noise and responsiveness class
all_genes <- unique(c(noise_results$gene, responsive_results$gene))
genomic_features <- data.frame(
  gene = all_genes,
  high_noise = all_genes %in% high_noise_genes,
  low_noise = all_genes %in% low_noise_genes,
  high_responsive = all_genes %in% high_resp_genes,
  low_responsive = all_genes %in% low_resp_genes
)

# Add gbM status
genomic_features$gbm <- all_genes %in% noise_results$gene[noise_results$gbm]

# Save annotated table
write.csv(genomic_features, "/group/sms029/mnieuwenh/gbm_noise_analysis/results/genomic_features/genomic_features_annotated.csv", row.names=FALSE)

# Overlap analysis: shared features between noisy and responsive genes
shared_high <- sum(genomic_features$high_noise & genomic_features$high_responsive)
shared_low <- sum(genomic_features$low_noise & genomic_features$low_responsive)

# Fisher's exact test: gbM frequency in high/low noise and high/low responsive genes
fisher_noise <- fisher.test(table(genomic_features$gbm, genomic_features$high_noise))
fisher_resp <- fisher.test(table(genomic_features$gbm, genomic_features$high_responsive))

# Save summary report
report_file <- "/group/sms029/mnieuwenh/gbm_noise_analysis/results/genomic_features/genomic_features_report.txt"
cat(
  "Genomic Feature Analysis Report\n",
  "==============================\n\n",
  "1. High/Low Noise Genes:\n",
  "   - Top 10% CV: high-noise (n=", sum(genomic_features$high_noise), ")\n",
  "   - Bottom 10% CV: low-noise (n=", sum(genomic_features$low_noise), ")\n\n",
  "2. High/Low Responsive Genes:\n",
  "   - Top 10% responsiveness: high-responsive (n=", sum(genomic_features$high_responsive), ")\n",
  "   - Bottom 10% responsiveness: low-responsive (n=", sum(genomic_features$low_responsive), ")\n\n",
  "3. Shared Features:\n",
  "   - Genes that are both high-noise and high-responsive: ", shared_high, "\n",
  "   - Genes that are both low-noise and low-responsive: ", shared_low, "\n\n",
  "4. Fisher's Exact Test Results:\n",
  "   - gbM frequency in high- vs low-noise genes: p-value = ", fisher_noise$p.value, "\n",
  "   - gbM frequency in high- vs low-responsive genes: p-value = ", fisher_resp$p.value, "\n",
  sep="",
  file=report_file
)
