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

# Add methylation and H2A.Z annotation for each gene
meth_anno <- read.csv("/group/sms029/mnieuwenh/Methylation_Data/combined_methylation_data.csv", stringsAsFactors=FALSE)
genomic_features <- merge(genomic_features, meth_anno, by.x="gene", by.y="Gene_ID", all.x=TRUE)

# Add H2A.Z group
genomic_features$H2AZ_Group <- NA
genomic_features$H2AZ_Group[genomic_features$H2AZ_Depleted == TRUE] <- "H2A.Z-Depleted"
genomic_features$H2AZ_Group[genomic_features$H2AZ_Enriched == TRUE] <- "H2A.Z-Enriched"
genomic_features$H2AZ_Group[is.na(genomic_features$H2AZ_Group)] <- "Other/NA"
genomic_features$H2AZ_Group <- factor(genomic_features$H2AZ_Group, levels=c("H2A.Z-Depleted", "H2A.Z-Enriched", "Other/NA"))

# Add methylation group (Cahn and Bewick)
genomic_features$Cahn_Group <- genomic_features$Cahn_Methylation_status
genomic_features$Bewick_Group <- genomic_features$Bewick_Classification

# Overlap analysis: shared features between noisy and responsive genes by methylation and H2A.Z group
shared_table <- genomic_features %>%
  group_by(Cahn_Group, Bewick_Group, H2AZ_Group) %>%
  summarise(
    n = n(),
    high_noise = sum(high_noise),
    low_noise = sum(low_noise),
    high_responsive = sum(high_responsive),
    low_responsive = sum(low_responsive),
    shared_high = sum(high_noise & high_responsive),
    shared_low = sum(low_noise & low_responsive)
  )
write.csv(shared_table, "/group/sms029/mnieuwenh/gbm_noise_analysis/results/genomic_features/genomic_features_overlap_summary.csv", row.names=FALSE)

# Fisher's exact tests for enrichment in high/low noise and responsiveness by methylation and H2A.Z group
fisher_results <- list()
for (grp in c("Cahn_Group", "Bewick_Group", "H2AZ_Group")) {
  for (val in unique(genomic_features[[grp]])) {
    if (is.na(val)) next
    tab_noise <- table(genomic_features[[grp]] == val, genomic_features$high_noise)
    tab_resp <- table(genomic_features[[grp]] == val, genomic_features$high_responsive)
    if (all(dim(tab_noise) == c(2,2))) {
      fisher_results[[paste0(grp, "_", val, "_noise")]] <- fisher.test(tab_noise)
    }
    if (all(dim(tab_resp) == c(2,2))) {
      fisher_results[[paste0(grp, "_", val, "_resp")]] <- fisher.test(tab_resp)
    }
  }
}

# Save all Fisher's test results in a table
fisher_table <- data.frame(
  Group = character(),
  Value = character(),
  Test = character(),
  P_value = numeric(),
  Odds_Ratio = numeric(),
  stringsAsFactors=FALSE
)
for (nm in names(fisher_results)) {
  parts <- strsplit(nm, "_")[[1]]
  grp <- parts[1]
  val <- paste(parts[2:(length(parts)-1)], collapse="_")
  test <- parts[length(parts)]
  res <- fisher_results[[nm]]
  fisher_table <- rbind(fisher_table, data.frame(
    Group=grp, Value=val, Test=test, P_value=res$p.value, Odds_Ratio=res$estimate))
}
write.csv(fisher_table, "/group/sms029/mnieuwenh/gbm_noise_analysis/results/genomic_features/genomic_features_fisher_summary.csv", row.names=FALSE)

# Save annotated table
write.csv(genomic_features, "/group/sms029/mnieuwenh/gbm_noise_analysis/results/genomic_features/genomic_features_annotated.csv", row.names=FALSE)

# Save a summary text file
overlap_high <- sum(genomic_features$high_noise & genomic_features$high_responsive)
overlap_low <- sum(genomic_features$low_noise & genomic_features$low_responsive)
summary_text <- c(
  "Genomic Feature Analysis Summary:",
  sprintf("Genes both high-noise and high-responsive: %d", overlap_high),
  sprintf("Genes both low-noise and low-responsive: %d", overlap_low),
  "See overlap and Fisher's test summary tables for details."
)
writeLines(summary_text, "/group/sms029/mnieuwenh/gbm_noise_analysis/results/genomic_features/genomic_features_summary.txt")
