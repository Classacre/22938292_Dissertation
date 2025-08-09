# ==============================================================================
# SCRIPT: 08_bootstrap_analysis.R
#
# PURPOSE:
#   This script performs a bootstrap analysis to validate the stability and
#   robustness of the final multivariate linear model.
#
#   FIX:
#   - Corrects the 97.5% quantile typo in CI computation.
# ==============================================================================

# --- 1. SETUP: Load libraries and define paths ---
packages_to_load <- c("dplyr", "ggplot2", "patchwork", "broom", "purrr")
for (pkg in packages_to_load) {
  if (!requireNamespace(pkg, quietly = TRUE)) install.packages(pkg, repos = "https://cloud.r-project.org/")
  library(pkg, character.only = TRUE)
}

# --- Define I/O Paths ---
# INPUT
annotated_noise_data_path <- "/group/sms029/mnieuwenh/gbm_noise_analysis/results/05_analyze_responsive_genes/05_noise_data_with_responsiveness.csv"
# OUTPUT
output_dir <- "/group/sms029/mnieuwenh/gbm_noise_analysis/results/08_bootstrap_analysis"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
bootstrap_plot_path <- file.path(output_dir, "08_bootstrap_coefficient_distributions.png")
bootstrap_summary_path <- file.path(output_dir, "08_bootstrap_summary_statistics.csv")
summary_file_path <- file.path(output_dir, "08_bootstrap_summary.txt")

# --- 2. DATA PREPARATION FOR MODELING ---
message("Loading and preparing data for bootstrapping...")
if (!file.exists(annotated_noise_data_path)) stop("Annotated noise data file not found!")
all_data <- read.csv(annotated_noise_data_path)

# Perform the same aggregation logic as in script 06 for consistency
model_data <- all_data %>%
  filter(
    mean_expr > 0.01, pct_expressing > 0.1,
    is.finite(variance.standardized), !is.na(cahn_group)
  ) %>%
  mutate(avg_log2FC = ifelse(is.na(avg_log2FC), 0, avg_log2FC)) %>%
  group_by(gene, cahn_group, has_TATA_box, avg_log2FC) %>%
  summarise(
    variance.standardized = mean(variance.standardized, na.rm = TRUE),
    log10_mean_expr = log10(mean(mean_expr, na.rm = TRUE)),
    .groups = "drop"
  ) %>%
  filter(is.finite(log10_mean_expr))

# Set reference levels for factors
model_data$cahn_group <- factor(model_data$cahn_group, levels = c("Unmethylated", "gbM", "TE-like"))

# --- 3. RUN BOOTSTRAP ANALYSIS ---
N_BOOTSTRAPS <- 1000
message(paste("Running bootstrap analysis with", N_BOOTSTRAPS, "iterations..."))

safe_lm <- purrr::safely(lm)

bootstrap_results <- purrr::map_dfr(seq_len(N_BOOTSTRAPS), function(i) {
  boot_sample <- dplyr::sample_n(model_data, size = nrow(model_data), replace = TRUE)
  fit <- safe_lm(
    variance.standardized ~ log10_mean_expr + has_TATA_box + cahn_group * avg_log2FC,
    data = boot_sample
  )
  if (is.null(fit$error)) {
    broom::tidy(fit$result) %>% dplyr::mutate(bootstrap_id = i, .before = 1)
  } else {
    NULL
  }
})

message("Bootstrap analysis complete.")

# --- 4. SUMMARIZE BOOTSTRAP RESULTS ---
message("Summarizing bootstrap results...")
bootstrap_summary <- bootstrap_results %>%
  group_by(term) %>%
  summarise(
    n_success = n(),
    mean_estimate = mean(estimate, na.rm = TRUE),
    std_error = sd(estimate, na.rm = TRUE),
    ci_2.5 = quantile(estimate, 0.025, na.rm = TRUE),
    ci_97.5 = quantile(estimate, 0.975, na.rm = TRUE),
    .groups = "drop"
  )

# --- 5. VISUALIZE BOOTSTRAP DISTRIBUTIONS ---
message("Visualizing coefficient distributions...")
p_bootstrap <- ggplot(bootstrap_results, aes(x = estimate)) +
  geom_density(fill = "skyblue", alpha = 0.7) +
  geom_vline(data = bootstrap_summary, aes(xintercept = mean_estimate), color = "red", linetype = "dashed") +
  geom_vline(data = bootstrap_summary, aes(xintercept = ci_2.5), color = "black", linetype = "dotted") +
  geom_vline(data = bootstrap_summary, aes(xintercept = ci_97.5), color = "black", linetype = "dotted") +
  facet_wrap(~ term, scales = "free", ncol = 2) +
  labs(
    title = "Bootstrap Distributions of Model Coefficients",
    subtitle = "Confidence intervals from 1000 bootstrap samples",
    x = "Coefficient Estimate",
    y = "Density"
  ) +
  theme_bw()

ggsave(bootstrap_plot_path, plot = p_bootstrap, width = 10, height = 12, dpi = 300)

# --- 6. SAVE RESULTS AND SUMMARY ---
message("Saving final bootstrap summary...")
write.csv(bootstrap_summary, bootstrap_summary_path, row.names = FALSE)

sink(summary_file_path)
cat("=================================================================\n")
cat("         STATISTICAL SUMMARY: 08_bootstrap_analysis.R\n")
cat("=================================================================\n\n")
cat(paste("Summary generated on:", Sys.Date()), "\n\n")
cat("--- Methodology ---\n")
cat(paste("A bootstrap analysis with", N_BOOTSTRAPS, "iterations was performed on the final interaction model.\n"))
cat("This analysis repeatedly samples the data and refits the model to assess the stability\n")
cat("and confidence in each coefficient's estimate.\n\n")
cat("Final Model Formula:\n")
cat("lm(Noise ~ log10_mean_expr + Has_TATA_Box + Methylation_Group * Responsiveness_Score)\n\n")
cat("--- Bootstrap Summary Statistics ---\n")
cat("Mean estimate and 95% CI for each term across successful bootstrap fits.\n\n")
print(bootstrap_summary, n = 20)
sink()

message(paste("Script 08 finished. Outputs saved in:", output_dir))