# ==============================================================================
# SCRIPT: 04_subsampling_validation.R
#
# PURPOSE:
#   This script performs a mean-expression-matched subsampling analysis to
#   validate that differences in the corrected noise metric are not an
#   artifact of residual correlation with mean expression.
# ==============================================================================

# --- 1. SETUP: Load libraries and define paths ---

packages_to_load <- c("dplyr", "tidyr", "ggplot2")
for (pkg in packages_to_load) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, repos = "https://cloud.r-project.org/")
  }
  library(pkg, character.only = TRUE)
}

# --- Define I/O Paths ---
master_summary_path <- "/group/sms029/mnieuwenh/gbm_noise_analysis/results/03_master_summary/master_gene_summary_table.csv"
complete_noise_data_path <- "/group/sms029/mnieuwenh/gbm_noise_analysis/results/01_noise_analysis/noise_analysis_complete_data.csv"

# OUTPUT
output_dir <- "/group/sms029/mnieuwenh/gbm_noise_analysis/results/04_subsampling_validation"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

set.seed(42)

# --- 2. DATA PREPARATION ---

message("Loading and preparing data for subsampling...")
master_data <- read.csv(master_summary_path)
complete_noise_data <- read.csv(complete_noise_data_path)

mean_expr_for_matching <- complete_noise_data %>%
  filter(mean_expr > 0.01, pct_expressing > 0.1) %>%
  group_by(gene) %>%
  summarise(median_mean_expr = median(mean_expr, na.rm = TRUE))

# CORRECTED: Filter using 'variance.standardized'
analysis_data <- master_data %>%
  left_join(mean_expr_for_matching, by = "gene") %>%
  filter(!is.na(median_mean_expr) & !is.na(variance.standardized))


# --- 3. CORE FUNCTION: Perform Binned Sub-Sampling ---

perform_matched_subsampling <- function(data, classification_col, group1_name, group2_name, match_col, num_bins = 20) {
  message(paste("--- Performing matching for", classification_col, "on", match_col, "---"))

  group1_pool <- data %>% filter(.data[[classification_col]] == group1_name)
  group2_pool <- data %>% filter(.data[[classification_col]] == group2_name)
  message(paste("Initial pool sizes:", group1_name, "=", nrow(group1_pool), ";", group2_name, "=", nrow(group2_pool)))

  combined_range <- range(data[[match_col]], na.rm = TRUE)
  breaks <- seq(combined_range[1], combined_range[2], length.out = num_bins + 1)
  breaks[length(breaks)] <- breaks[length(breaks)] + 1e-9

  group1_pool$bin <- cut(group1_pool[[match_col]], breaks = breaks, include.lowest = TRUE, labels = FALSE)
  group2_pool$bin <- cut(group2_pool[[match_col]], breaks = breaks, include.lowest = TRUE, labels = FALSE)

  matched_genes <- lapply(1:num_bins, function(i) {
    genes1_in_bin <- group1_pool %>% filter(bin == i)
    genes2_in_bin <- group2_pool %>% filter(bin == i)
    n_to_sample <- min(nrow(genes1_in_bin), nrow(genes2_in_bin))
    if (n_to_sample > 0) {
      bind_rows(
        sample_n(genes1_in_bin, n_to_sample),
        sample_n(genes2_in_bin, n_to_sample)
      )
    } else {
      NULL
    }
  })

  final_matched_df <- bind_rows(matched_genes)
  message(paste("Final matched set size:", nrow(final_matched_df), "genes total."))
  return(final_matched_df)
}


# --- 4. EXECUTE ANALYSIS & VISUALIZE RESULTS ---

run_and_visualize_analysis <- function(data, classification_col, group1, group2) {
  match_col <- "median_mean_expr"
  # CORRECTED: Use 'variance.standardized' as the noise column
  noise_col <- "variance.standardized"

  matched_df <- perform_matched_subsampling(data, classification_col, group1, group2, match_col)

  p_verify <- ggplot(matched_df, aes(x = .data[[match_col]], fill = .data[[classification_col]])) +
    geom_density(alpha = 0.7) +
    scale_x_log10(labels = scales::label_number(accuracy = 0.01)) +
    labs(
      title = paste("Verification of Expression Matching (", classification_col, ")"),
      subtitle = "Distributions of matched sets should be nearly identical.",
      x = "Median Mean Expression (log10 scale)", y = "Density", fill = "Gene Group"
    ) +
    theme_bw(base_size = 14)
  ggsave(file.path(output_dir, paste0("verification_plot_", classification_col, ".png")), plot = p_verify, width = 8, height = 6)

  p_compare <- ggplot(matched_df, aes(x = .data[[classification_col]], y = .data[[noise_col]], fill = .data[[classification_col]])) +
    geom_violin(trim = FALSE, alpha = 0.7) +
    geom_boxplot(width = 0.1, outlier.size = 0.5, fill = "white") +
    labs(
      title = paste("Corrected Noise of Matched Sets (", classification_col, ")"),
      subtitle = "Comparison after controlling for mean expression.",
      x = "Gene Group", y = "Corrected Noise (Standardized Variance)"
    ) +
    theme_bw(base_size = 14) + theme(legend.position = "none")
  ggsave(file.path(output_dir, paste0("noise_comparison_plot_", classification_col, ".png")), plot = p_compare, width = 8, height = 7)

  group1_data <- matched_df %>% filter(.data[[classification_col]] == group1)
  group2_data <- matched_df %>% filter(.data[[classification_col]] == group2)
  stat_test <- wilcox.test(group1_data[[noise_col]], group2_data[[noise_col]], alternative = "two.sided")
  return(stat_test)
}

cahn_test_result <- run_and_visualize_analysis(analysis_data, "cahn_group", "gbM", "Unmethylated")
bewick_test_result <- run_and_visualize_analysis(analysis_data, "bewick_group", "gbM", "Unmethylated")


# --- 5. SAVE SUMMARY OF FINDINGS ---

sink(file.path(output_dir, "subsampling_validation_summary.txt"))
cat("====================================================================\n")
cat(" Expression-Matched Subsampling Validation - Summary\n")
cat("====================================================================\n\n")
cat("CONTEXT: The primary analysis uses the corrected noise score (standardized variance),\n")
cat("which already accounts for the mean-variance relationship. This subsampling\n")
cat("analysis serves as an independent, rigorous validation of those findings.\n\n")
cat("METHOD: We compare the corrected noise scores between gene groups that have\n")
cat("been explicitly matched for their median mean expression level.\n\n")
cat("--------------------------------------------------------------------\n")
cat(" Cahn (2016) Classification: gbM vs. Unmethylated\n")
cat("--------------------------------------------------------------------\n")
cat("Wilcoxon rank-sum test on corrected noise of matched sets:\n")
print(cahn_test_result)
cat("\nCONCLUSION:\n")
if (cahn_test_result$p.value < 0.05) {
  cat("A significant difference in corrected noise remains even after matching for\n")
  cat("expression level. This strongly supports the conclusion that the association\n")
  cat("is independent of mean expression.\n")
} else {
  cat("No significant difference in corrected noise was found after matching.\n")
}
cat("\n\n--------------------------------------------------------------------\n")
cat(" Bewick (2017) Classification: gbM vs. Unmethylated\n")
cat("--------------------------------------------------------------------\n")
cat("Wilcoxon rank-sum test on corrected noise of matched sets:\n")
print(bewick_test_result)
cat("\nCONCLUSION:\n")
if (bewick_test_result$p.value < 0.05) {
  cat("A significant difference in corrected noise remains even after matching for\n")
  cat("expression level. This strongly supports the conclusion that the association\n")
  cat("is independent of mean expression.\n")
} else {
  cat("No significant difference in corrected noise was found after matching.\n")
}
sink()

message("Script finished successfully. All outputs are in the '04_subsampling_validation' directory.")