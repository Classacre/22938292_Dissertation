# ==============================================================================
# SCRIPT: 05_visualize_filtering.R
#
# PURPOSE:
#   This script creates a diagnostic plot to visualize the effect of the data
#   filtering steps applied in script 01.
# ==============================================================================

# --- 1. SETUP: Load libraries and define paths ---

packages_to_load <- c("dplyr", "ggplot2", "scales")
for (pkg in packages_to_load) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, repos = "https://cloud.r-project.org/")
  }
  library(pkg, character.only = TRUE)
}

# --- Define I/O Paths ---
complete_noise_data_path <- "/group/sms029/mnieuwenh/gbm_noise_analysis/results/01_noise_analysis/noise_analysis_complete_data.csv"

# OUTPUT
output_dir <- "/group/sms029/mnieuwenh/gbm_noise_analysis/results/05_diagnostic_plots"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

filter_plot_path <- file.path(output_dir, "mean_vs_cv_filtering_visualization.png")


# --- 2. DATA LOADING AND PREPARATION ---

message("Loading complete noise data...")
if (!file.exists(complete_noise_data_path)) {
  stop(paste("Error: Input data file not found at", complete_noise_data_path))
}
noise_data <- read.csv(complete_noise_data_path)

message("Preparing data for plotting...")
plot_data <- noise_data %>%
  filter(is.finite(mean_expr) & mean_expr > 0 & is.finite(cv2) & cv2 > 0) %>%
  mutate(
    cv = sqrt(cv2),
    filter_status = case_when(
      mean_expr > 0.01 & pct_expressing > 0.1 ~ "Passed",
      mean_expr <= 0.01 & pct_expressing > 0.1 ~ "Filtered (Low Mean Expression)",
      mean_expr > 0.01 & pct_expressing <= 0.1 ~ "Filtered (Low Detection Rate)",
      TRUE ~ "Filtered (Both)"
    )
  ) %>%
  mutate(filter_status = factor(
    filter_status,
    levels = c("Passed", "Filtered (Low Detection Rate)", "Filtered (Low Mean Expression)", "Filtered (Both)")
  ))

passed_data <- plot_data %>% filter(filter_status == "Passed")
filtered_data <- plot_data %>% filter(filter_status != "Passed")

set.seed(123)
if (nrow(filtered_data) > 200000) {
  sampled_filtered_data <- filtered_data %>% sample_n(200000)
} else {
  sampled_filtered_data <- filtered_data
}

plot_data_sampled <- bind_rows(passed_data, sampled_filtered_data)

message(paste(
  "Plotting", nrow(plot_data_sampled), "data points (all",
  nrow(passed_data), "passing points, and a sample of filtered points)."
))


# --- 3. CREATE AND SAVE THE VISUALIZATION ---

message("Generating the mean-CV plot...")

filter_colors <- c(
  "Passed" = "black",
  "Filtered (Low Detection Rate)" = "lightskyblue",
  "Filtered (Low Mean Expression)" = "salmon",
  "Filtered (Both)" = "plum"
)

p <- ggplot(plot_data_sampled, aes(x = mean_expr, y = cv, color = filter_status)) +
  geom_point(size = 0.5, alpha = 0.4) +
  geom_function(
    fun = function(x) 1/sqrt(x),
    color = "darkred",
    linetype = "dashed",
    linewidth = 1
  ) +
  geom_vline(xintercept = 0.01, linetype = "dotted", color = "gray40", linewidth = 0.8) +
  scale_x_log10(
    labels = trans_format("log10", math_format(10^.x)),
    breaks = 10^seq(-4, 4, by = 2)
  ) +
  scale_y_log10(
    labels = trans_format("log10", math_format(10^.x)),
    breaks = 10^seq(-2, 4, by = 2)
  ) +
  scale_color_manual(values = filter_colors) +
  guides(color = guide_legend(override.aes = list(size = 3, alpha = 1))) +
  annotate("text", x = 10, y = 0.05, label = "Shot Noise\n(CV = 1/sqrt(Mean))", color = "darkred", hjust = 0) +
  annotate("text", x = 0.008, y = 1000, label = "Mean Expr = 0.01", color = "gray40", angle = 90, hjust = 1, vjust = -0.5) +
  labs(
    title = "Visualization of Gene Filtering based on Expression Level and Detection Rate",
    subtitle = paste("Showing", format(nrow(passed_data), big.mark=","), "gene-celltype pairs that passed filters"),
    x = "Mean Gene Expression (Linear Scale)",
    y = "Coefficient of Variation (CV = SD / Mean)",
    color = "Filter Status"
  ) +
  theme_bw(base_size = 14) +
  theme(
    legend.position = "bottom",
    legend.title = element_text(face = "bold"),
    plot.title = element_text(face = "bold"),
    panel.grid.minor = element_blank()
  )

message(paste("Saving plot to:", filter_plot_path))
ggsave(filter_plot_path, plot = p, width = 12, height = 9, dpi = 300)

# --- 4. FINAL CONFIRMATION ---
summary_text <- c(
  "Filter Visualization Script Complete.",
  "",
  paste("A diagnostic plot showing the effect of data filtering has been generated."),
  "The plot visualizes all gene-celltype combinations, coloring them based on whether they passed or failed the filtering criteria:",
  "  - Mean Expression > 0.01",
  "  - Percentage of Expressing Cells > 0.1",
  "",
  paste("The plot has been saved to:", filter_plot_path)
)
writeLines(summary_text, file.path(output_dir, "filtering_visualization_summary.txt"))

message("Script finished successfully.")