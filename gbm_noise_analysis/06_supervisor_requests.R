# ==============================================================================
# SCRIPT: 06_supervisor_requests_diag.R
#
# PURPOSE:
#   This is a DIAGNOSTIC version of the script to identify why execution might
#   be halting. It adds extensive logging to track progress and memory usage.
#
# TASKS:
#   1. Correlate the raw noise metric (CV²) with the corrected VST-based
#      noise metric (variance.standardized).
#   2. Create diagnostic plots showing data distributions before and after
#      expression-matched subsampling.
# ==============================================================================

# --- 1. SETUP: Load libraries and define paths ---
message("Step 1: SETUP - Loading required libraries...")

# Added 'readr' for faster file reading
packages_to_load <- c("dplyr", "tidyr", "ggplot2", "patchwork", "ggpubr", "readr")
for (pkg in packages_to_load) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    message(paste("Installing package:", pkg))
    install.packages(pkg, repos = "https://cloud.r-project.org/")
  }
  library(pkg, character.only = TRUE)
}

message("Libraries loaded successfully.")

# --- Define I/O Paths ---
master_summary_path <- "/group/sms029/mnieuwenh/gbm_noise_analysis/results/03_master_summary/master_gene_summary_table.csv"
complete_noise_data_path <- "/group/sms029/mnieuwenh/gbm_noise_analysis/results/01_noise_analysis/noise_analysis_complete_data.csv"

# OUTPUTS
output_dir <- "/group/sms029/mnieuwenh/gbm_noise_analysis/results/06_supervisor_requests"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

set.seed(42) # for reproducibility
message("Setup complete. Starting Task 1.")


# --- 2. TASK 1: Correlate CV vs. VST Corrected Noise ---

message("\n--- TASK 1: Correlate CV vs. VST score ---")
message(paste("Loading master summary data from:", master_summary_path))

if (!file.exists(master_summary_path)) {
  stop("FATAL ERROR: Master summary file not found. Halting execution.")
}
# Using readr::read_csv for better performance
master_data <- readr::read_csv(master_summary_path, show_col_types = FALSE)
message("Master summary data loaded.")
message(paste("Dimensions of master_data:", paste(dim(master_data), collapse = " x ")))
message(paste("Memory usage of master_data:", format(object.size(master_data), units = "auto")))


message("Preparing data for CV vs. VST plot...")
cv_vst_data <- master_data %>%
  filter(!is.na(median_within_celltype_cv2) & !is.na(variance.standardized)) %>%
  filter(median_within_celltype_cv2 > 0)

message("Data prepared. Generating plot...")
p_cv_vst <- ggplot(cv_vst_data, aes(x = median_within_celltype_cv2, y = variance.standardized)) +
  geom_point(alpha = 0.2, size = 1) +
  scale_x_log10(
    name = "Median Within-Cell-Type CV² (log10 scale)",
    labels = scales::trans_format("log10", scales::math_format(10^.x))
  ) +
  scale_y_continuous(name = "Corrected Noise (VST Standardized Variance)") +
  geom_smooth(method = "gam", formula = y ~ s(x, bs = "cs"), color = "red") +
  stat_cor(method = "spearman", label.x.npc = 0.05, label.y.npc = 0.95, cor.coef.name = "rho") +
  labs(
    title = "Correlation between Raw and Corrected Noise Metrics",
    subtitle = "Comparing median CV² to the VST-based standardized variance for each gene."
  ) +
  theme_bw(base_size = 14)

message("Plot generated. Saving CV vs. VST plot...")
plot_path_1 <- file.path(output_dir, "correlation_cv2_vs_vst_variance.png")
ggsave(plot_path_1, plot = p_cv_vst, width = 8, height = 7)
message(paste("Plot saved to:", plot_path_1))
message("--- TASK 1 COMPLETE ---")


# --- 3. TASK 2: Subsampling Before-and-After Visualization ---

message("\n--- TASK 2: Subsampling Before-and-After Visualization ---")

# --- 3a. Data Preparation ---
message("Preparing data for subsampling analysis...")
message(paste("Loading complete noise data from:", complete_noise_data_path))

if (!file.exists(complete_noise_data_path)) {
  stop("FATAL ERROR: Complete noise data file not found. Halting execution.")
}
complete_noise_data <- readr::read_csv(complete_noise_data_path, show_col_types = FALSE)
message("Complete noise data loaded.")
message(paste("Dimensions of complete_noise_data:", paste(dim(complete_noise_data), collapse = " x ")))
message(paste("Memory usage of complete_noise_data:", format(object.size(complete_noise_data), units = "auto")))


message("Calculating median mean expression per gene...")
mean_expr_for_matching <- complete_noise_data %>%
  filter(mean_expr > 0.01, pct_expressing > 0.1) %>%
  group_by(gene) %>%
  summarise(median_mean_expr = median(mean_expr, na.rm = TRUE))
message("Median mean expression calculated.")

message("Joining master data with expression data to create final analysis table...")
analysis_data <- master_data %>%
  left_join(mean_expr_for_matching, by = "gene") %>%
  filter(!is.na(median_mean_expr) & !is.na(variance.standardized))
message("Final analysis table created.")
message(paste("Dimensions of analysis_data:", paste(dim(analysis_data), collapse = " x ")))
message(paste("Memory usage of analysis_data:", format(object.size(analysis_data), units = "auto")))

# Clean up large objects to free memory
rm(complete_noise_data, master_data)
gc() # Force garbage collection
message("Cleaned up intermediate large objects to free memory.")


# --- 3b. Re-usable Subsampling Function ---
perform_matched_subsampling <- function(data, classification_col, group1_name, group2_name, match_col, num_bins = 20) {
  # ... (function is unchanged)
  group1_pool <- data %>% filter(.data[[classification_col]] == group1_name)
  group2_pool <- data %>% filter(.data[[classification_col]] == group2_name)
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
  return(bind_rows(matched_genes))
}


# --- 3c. Plotting Function for Before-and-After ---
generate_subsampling_plots <- function(data, classification_col, group1, group2) {
  # ... (function is unchanged)
  match_col <- "median_mean_expr"
  noise_col <- "variance.standardized"
  before_df <- data %>% filter(.data[[classification_col]] %in% c(group1, group2))
  message(paste("Performing subsampling for", classification_col, "..."))
  after_df <- perform_matched_subsampling(data, classification_col, group1, group2, match_col)
  message("Subsampling complete. Generating plots...")
  p_before_match <- ggplot(before_df, aes(x = .data[[match_col]], fill = .data[[classification_col]])) +
    geom_density(alpha = 0.7) + scale_x_log10(labels = scales::label_number(accuracy = 0.01)) +
    labs(title = "Before Matching", x = "Median Mean Expression", y = "Density") +
    theme_bw() + theme(legend.position = "bottom")
  p_before_noise <- ggplot(before_df, aes(x = .data[[classification_col]], y = .data[[noise_col]], fill = .data[[classification_col]])) +
    geom_violin(trim = FALSE, alpha = 0.7, show.legend = FALSE) +
    geom_boxplot(width = 0.1, fill = "white", show.legend = FALSE) +
    labs(title = "Before Matching", x = "Gene Group", y = "Corrected Noise") + theme_bw()
  p_after_match <- ggplot(after_df, aes(x = .data[[match_col]], fill = .data[[classification_col]])) +
    geom_density(alpha = 0.7) + scale_x_log10(labels = scales::label_number(accuracy = 0.01)) +
    labs(title = "After Matching", x = "Median Mean Expression", y = "Density") +
    theme_bw() + theme(legend.position = "bottom")
  p_after_noise <- ggplot(after_df, aes(x = .data[[classification_col]], y = .data[[noise_col]], fill = .data[[classification_col]])) +
    geom_violin(trim = FALSE, alpha = 0.7, show.legend = FALSE) +
    geom_boxplot(width = 0.1, fill = "white", show.legend = FALSE) +
    labs(title = "After Matching", x = "Gene Group", y = "Corrected Noise") + theme_bw()
  combined_plot <- (p_before_match + p_before_noise) / (p_after_match + p_after_noise) +
    plot_layout(guides = "collect") & theme(legend.position = 'bottom')
  combined_plot <- combined_plot + plot_annotation(
    title = paste("Subsampling Validation:", group1, "vs.", group2),
    subtitle = "Top row shows all data; bottom row shows data after matching for mean expression.",
    theme = theme(plot.title = element_text(size = 18, face = "bold"))
  )
  return(combined_plot)
}

# --- 3d. Execute and Save ---
message("\nStarting plot generation for Cahn classification...")
cahn_plots <- generate_subsampling_plots(analysis_data, "cahn_group", "gbM", "Unmethylated")
message("Cahn plots generated. Saving file...")
plot_path_2 <- file.path(output_dir, "subsampling_before_after_cahn.png")
ggsave(plot_path_2, plot = cahn_plots, width = 10, height = 10)
message(paste("Cahn subsampling plot saved to:", plot_path_2))

message("\nStarting plot generation for Bewick classification...")
bewick_plots <- generate_subsampling_plots(analysis_data, "bewick_group", "gbM", "Unmethylated")
message("Bewick plots generated. Saving file...")
plot_path_3 <- file.path(output_dir, "subsampling_before_after_bewick.png")
ggsave(plot_path_3, plot = bewick_plots, width = 10, height = 10)
message(paste("Bewick subsampling plot saved to:", plot_path_3))

message("\n--- TASK 2 COMPLETE ---")
message("Script finished successfully.")