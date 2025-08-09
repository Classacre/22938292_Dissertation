# ==============================================================================
# SCRIPT: 04_validate_gene_filtering.R
#
# PURPOSE:
#   Validate gene filtering and visualize mean–variance behavior and corrected
#   noise metrics. This version is robust to cases where strict filters leave
#   zero rows, and when variance.standardized and/or dm contain many NAs.
#
# CHANGES:
#   - Decouples base filtering (mean_expr, pct_expressing) from noise-metric
#     validity to avoid empty datasets.
#   - Falls back to relaxed thresholds if the base filter yields zero rows.
#   - Builds each panel from the subset that has the required columns available.
#   - If a panel has no data, it renders a placeholder with an explanatory note.
#   - Summary file now reports counts for VST-only, DM-only, and overlap.
# ==============================================================================

# --- 1. SETUP: Load libraries and define paths ---
packages_to_load <- c("dplyr", "ggplot2", "patchwork", "viridis")
for (pkg in packages_to_load) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, repos = "https://cloud.r-project.org/")
  }
  library(pkg, character.only = TRUE)
}

# --- Define I/O Paths ---
# INPUT
noise_data_path <- "/group/sms029/mnieuwenh/gbm_noise_analysis/results/03_calculate_expression_noise/03_gene_noise_metrics_complete.csv"
# OUTPUT
output_dir <- "/group/sms029/mnieuwenh/gbm_noise_analysis/results/04_validate_gene_filtering"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
output_plot_path <- file.path(output_dir, "04_filtering_and_validation_figure.png")
summary_file_path <- file.path(output_dir, "04_filtering_summary.txt")

# --- 2. DATA LOADING AND CLEANING ---
message("Loading complete noise data...")
if (!file.exists(noise_data_path)) stop("Noise data file not found!")
noise_data <- read.csv(noise_data_path, check.names = FALSE)

# Coerce critical columns to numeric if needed
num_cols <- c("mean_expr", "cv2", "variance.standardized", "dm", "pct_expressing")
for (nm in num_cols) {
  if (!nm %in% names(noise_data)) next
  if (!is.numeric(noise_data[[nm]])) {
    noise_data[[nm]] <- suppressWarnings(as.numeric(noise_data[[nm]]))
  }
}

# Replace Inf/-Inf with NA across numeric columns
noise_data[] <- lapply(noise_data, function(x) {
  if (is.numeric(x)) {
    x[!is.finite(x)] <- NA_real_
  }
  x
})

# Define a consistent plotting theme
publication_theme <- theme_bw() + theme(
  plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
  axis.title = element_text(size = 12),
  axis.text = element_text(size = 10),
  legend.title = element_text(size = 12),
  legend.text = element_text(size = 10),
  panel.grid.minor = element_blank()
)

empty_panel <- function(title, subtitle = "No data available after filtering", xlab = "", ylab = "") {
  ggplot() +
    theme_void() +
    labs(title = title, subtitle = subtitle, x = xlab, y = ylab) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"))
}

# --- 3. APPLY BASE FILTERING (robust) ---
message("Applying filtering criteria...")

# Primary thresholds
MIN_MEAN_EXPR <- 0.01
MIN_PCT_EXPRESSING <- 0.1
used_relaxed_thresholds <- FALSE

has_base_cols <- all(c("mean_expr", "pct_expressing") %in% names(noise_data))
if (!has_base_cols) stop("Required columns 'mean_expr' and/or 'pct_expressing' not found in noise_data.")

filtered_base <- noise_data %>%
  filter(
    is.finite(mean_expr), is.finite(pct_expressing),
    mean_expr > MIN_MEAN_EXPR,
    pct_expressing > MIN_PCT_EXPRESSING
  )

# If too few rows remain, relax thresholds so the plots won’t be empty.
if (nrow(filtered_base) == 0) {
  message("No rows after base filtering; relaxing thresholds (mean_expr > 0 & pct_expressing > 0)...")
  filtered_base <- noise_data %>%
    filter(
      is.finite(mean_expr), is.finite(pct_expressing),
      mean_expr > 0,
      pct_expressing > 0
    )
  used_relaxed_thresholds <- TRUE
}

# Counts for availability of metrics
n_total_before <- nrow(noise_data)
n_after_base <- nrow(filtered_base)
n_vst_avail <- sum(is.finite(filtered_base$`variance.standardized`), na.rm = TRUE)
n_dm_avail  <- sum(is.finite(filtered_base$dm), na.rm = TRUE)
n_both_avail <- sum(is.finite(filtered_base$`variance.standardized`) & is.finite(filtered_base$dm), na.rm = TRUE)

# --- 4. GENERATE VISUALIZATIONS (panel-wise robustness) ---
message("Generating visualization panels...")

# Helper to choose hex bin count based on sample size
choose_bins <- function(n) {
  if (is.na(n) || n <= 0) return(30L)
  b <- max(10L, min(100L, round(sqrt(n))))
  as.integer(b)
}

# Panel A: Mean vs. Variance (CV^2) – raw noise trend
panel_a_data <- filtered_base %>%
  filter(is.finite(mean_expr), is.finite(cv2), mean_expr > 0, cv2 > 0)

if (nrow(panel_a_data) > 0) {
  bins_a <- choose_bins(nrow(panel_a_data))
  panel_a <- ggplot(panel_a_data, aes(x = log10(mean_expr), y = log10(cv2))) +
    geom_hex(bins = bins_a) +
    scale_fill_viridis(option = "C", name = "Density") +
    geom_smooth(method = "gam", color = "red3", se = FALSE) +
    labs(
      title = "A. Raw Noise (CV²) vs. Mean Expression",
      subtitle = "Demonstrates mean–variance dependency",
      x = "Log10(Mean Expression)",
      y = "Log10(CV²)"
    ) +
    publication_theme
} else {
  panel_a <- empty_panel(
    title = "A. Raw Noise (CV²) vs. Mean Expression",
    subtitle = "Insufficient non-zero mean or CV² after filtering",
    xlab = "Log10(Mean Expression)", ylab = "Log10(CV²)"
  )
}

# Panel B: Mean vs. Corrected Noise (VST) – requires variance.standardized
panel_b_data <- filtered_base %>%
  filter(is.finite(`variance.standardized`), is.finite(mean_expr), mean_expr > 0)

if (nrow(panel_b_data) > 0) {
  bins_b <- choose_bins(nrow(panel_b_data))
  panel_b <- ggplot(panel_b_data, aes(x = log10(mean_expr), y = `variance.standardized`)) +
    geom_hex(bins = bins_b) +
    scale_fill_viridis(option = "C", name = "Density") +
    geom_smooth(method = "gam", color = "red3", se = FALSE) +
    labs(
      title = "B. Corrected Noise vs. Mean Expression (VST)",
      subtitle = "VST should reduce mean–variance trend",
      x = "Log10(Mean Expression)",
      y = "Corrected Noise (VST)"
    ) +
    publication_theme
} else {
  panel_b <- empty_panel(
    title = "B. Corrected Noise vs. Mean Expression (VST)",
    subtitle = "No finite VST-corrected noise available after filtering",
    xlab = "Log10(Mean Expression)", ylab = "Corrected Noise (VST)"
  )
}

# Panel C: Distribution of Corrected Noise
# Prefer VST if available; otherwise, fall back to DM and indicate it
use_vst_for_hist <- n_vst_avail > 0
if (use_vst_for_hist) {
  panel_c_data <- filtered_base %>% filter(is.finite(`variance.standardized`))
  metric_label <- "Corrected Noise (VST)"
  x_col <- "variance.standardized"
} else {
  panel_c_data <- filtered_base %>% filter(is.finite(dm))
  metric_label <- "Corrected Noise (DM)"
  x_col <- "dm"
}

if (nrow(panel_c_data) > 0) {
  panel_c <- ggplot(panel_c_data, aes(x = .data[[x_col]])) +
    geom_histogram(aes(y = after_stat(density)), bins = 100, fill = "skyblue", color = "black") +
    geom_density(color = "navy", linewidth = 1) +
    labs(
      title = "C. Distribution of Corrected Noise",
      subtitle = paste0("Metric: ", metric_label),
      x = metric_label,
      y = "Density"
    ) +
    publication_theme
} else {
  panel_c <- empty_panel(
    title = "C. Distribution of Corrected Noise",
    subtitle = "Neither VST nor DM is available for histogram"
  )
}

# Panel D: Metric Agreement – VST vs. DM correlation (requires both)
panel_d_data <- filtered_base %>%
  filter(is.finite(`variance.standardized`) & is.finite(dm))

if (nrow(panel_d_data) > 2) {
  correlation <- suppressWarnings(cor(panel_d_data$`variance.standardized`, panel_d_data$dm, method = "pearson"))
  bins_d <- choose_bins(nrow(panel_d_data))
  panel_d <- ggplot(panel_d_data, aes(x = `variance.standardized`, y = dm)) +
    geom_hex(bins = bins_d) +
    scale_fill_viridis(option = "C", name = "Density") +
    geom_smooth(method = "lm", color = "red3", se = FALSE) +
    labs(
      title = "D. Noise Metric Agreement",
      subtitle = paste0("Pearson r = ", ifelse(is.finite(correlation), round(correlation, 3), "NA")),
      x = "Corrected Noise (VST)",
      y = "Corrected Noise (DM)"
    ) +
    publication_theme
} else {
  correlation <- NA_real_
  panel_d <- empty_panel(
    title = "D. Noise Metric Agreement",
    subtitle = "Insufficient overlap (VST ∩ DM) to compute correlation",
    xlab = "Corrected Noise (VST)",
    ylab = "Corrected Noise (DM)"
  )
}

# --- 5. ASSEMBLE AND SAVE THE FINAL FIGURE ---
message("Assembling and saving the final figure...")
final_figure <- (panel_a | panel_b) / (panel_c | panel_d) +
  plot_annotation(
    title = "Validation of Gene Filtering and Noise Correction Pipeline",
    caption = "Panels adaptively use available metrics and state when data is unavailable."
  )

ggsave(output_plot_path, plot = final_figure, width = 14, height = 10, dpi = 300)

# --- 6. GENERATE SUMMARY FILE ---
message("Generating summary file...")
summary_lines <- c(
  "=================================================================",
  "         STATISTICAL SUMMARY: 04_validate_gene_filtering.R",
  "=================================================================",
  "",
  paste("Summary generated on:", Sys.Date()),
  "",
  "--- Filtering Summary ---",
  paste("Total gene-identity pairs before filtering:", n_total_before),
  if (!used_relaxed_thresholds) {
    paste("Filtering criteria: mean_expr >", MIN_MEAN_EXPR, "& pct_expressing >", MIN_PCT_EXPRESSING)
  } else {
    paste0("Filtering criteria: mean_expr > 0 & pct_expressing > 0 (relaxed; original ",
           "mean_expr >", MIN_MEAN_EXPR, " & pct_expressing >", MIN_PCT_EXPRESSING, " yielded 0 rows)")
  },
  paste("Total gene-identity pairs after base filtering:", n_after_base),
  paste("Rows with finite VST (variance.standardized):", n_vst_avail),
  paste("Rows with finite DM:", n_dm_avail),
  paste("Rows with both VST and DM:", n_both_avail),
  "",
  "--- Noise Metric Validation ---",
  paste("Pearson correlation (VST vs DM) on overlapping rows:", ifelse(is.finite(correlation), round(correlation, 4), "NA")),
  "Panels automatically fall back or display informative placeholders when specific",
  "metrics are not available after filtering."
)

writeLines(summary_lines, con = summary_file_path)
message(paste("Script finished. Output saved in:", output_dir))