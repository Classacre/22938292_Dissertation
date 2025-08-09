# ==============================================================================
# SCRIPT: 06_generate_summary_figures.R
#
# PURPOSE:
#   This script synthesizes all previous analyses into a final summary figure
#   and a comprehensive multivariate linear model. It aims to explain the
#   drivers of gene expression noise, both globally and within each cell type.
#
#   ENHANCEMENTS/FIXES:
#   - Adds Panel E to visualize a key interaction effect.
#   - Implements a linear mixed-effects model (LMM) with gene random intercepts.
#   - Uses I(log10(mean_expr)) inline in models to avoid missing column errors.
#   - Adds purrr for safely() and fixes a typo when printing model errors.
# ==============================================================================

# --- 1. SETUP: Load libraries and define paths ---
packages_to_load <- c("dplyr", "ggplot2", "patchwork", "ggpubr", "lme4", "purrr")
for (pkg in packages_to_load) {
  if (!requireNamespace(pkg, quietly = TRUE)) install.packages(pkg, repos = "https://cloud.r-project.org/")
  library(pkg, character.only = TRUE)
}

# --- Define I/O Paths ---
# INPUT
annotated_noise_data_path <- "/group/sms029/mnieuwenh/gbm_noise_analysis/results/05_analyze_responsive_genes/05_noise_data_with_responsiveness.csv"
# OUTPUT
output_dir <- "/group/sms029/mnieuwenh/gbm_noise_analysis/results/06_generate_summary_figures"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# Create the new subdirectory for cell type specific results
celltype_output_dir <- file.path(output_dir, "celltype_analyses")
dir.create(celltype_output_dir, showWarnings = FALSE, recursive = TRUE)

final_figure_path <- file.path(output_dir, "06_summary_figure_noise_drivers_GLOBAL.png")
celltype_figure_path <- file.path(celltype_output_dir, "06_summary_figure_noise_drivers_BY_CELL_TYPE.png")
model_summary_path <- file.path(output_dir, "06_multivariate_model_summary_GLOBAL.txt")
celltype_model_summary_path <- file.path(celltype_output_dir, "06_multivariate_model_summary_BY_CELL_TYPE.txt")

# --- 2. DATA LOADING AND PREPARATION ---
message("Loading annotated noise data...")
if (!file.exists(annotated_noise_data_path)) stop("Annotated noise data file not found!")
all_data <- read.csv(annotated_noise_data_path)

filtered_data <- all_data %>%
  filter(
    mean_expr > 0.01, pct_expressing > 0.1,
    is.finite(variance.standardized), !is.na(cahn_group), !is.na(base_cell_type)
  ) %>%
  mutate(avg_log2FC = ifelse(is.na(avg_log2FC), 0, avg_log2FC))

publication_theme <- theme_bw() + theme(
  plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
  axis.title = element_text(size = 12), axis.text = element_text(size = 10),
  legend.position = "bottom", panel.grid.minor = element_blank()
)

# --- 3. GENERATE GLOBAL SUMMARY FIGURE PANELS ---
message("Generating panels for the final global summary figure...")
y_zoom <- quantile(filtered_data$variance.standardized, c(0.01, 0.99), na.rm = TRUE)

panel_a <- ggplot(filtered_data, aes(x = cahn_group, y = variance.standardized, fill = cahn_group)) +
  geom_violin(trim = FALSE) +
  geom_boxplot(width = 0.1, outlier.shape = NA, fill = "white") +
  scale_fill_manual(values = c("gbM" = "#0072B2", "TE-like" = "#D55E00", "Unmethylated" = "grey")) +
  coord_cartesian(ylim = y_zoom) +
  labs(
    title = "A. gbM is associated with low noise",
    x = "Methylation Status", y = "Corrected Noise (VST)"
  ) +
  ggpubr::stat_compare_means(
    comparisons = list(c("gbM", "Unmethylated"), c("gbM", "TE-like")),
    method = "wilcox.test", label = "p.signif"
  ) +
  theme(legend.position = "none")

panel_b <- ggplot(filtered_data, aes(x = h2az_group, y = variance.standardized, fill = h2az_group)) +
  geom_violin(trim = FALSE) +
  geom_boxplot(width = 0.1, outlier.shape = NA, fill = "white") +
  scale_fill_manual(values = c("H2A.Z-Enriched" = "#009E73", "H2A.Z-Depleted" = "#CC79A7", "Other" = "grey")) +
  coord_cartesian(ylim = y_zoom) +
  labs(
    title = "B. H2A.Z status modulates noise",
    x = "H2A.Z Occupancy", y = "Corrected Noise (VST)"
  ) +
  theme(legend.position = "none")

panel_c <- ggplot(filtered_data, aes(x = is_responsive, y = variance.standardized, fill = is_responsive)) +
  geom_violin(trim = FALSE) +
  geom_boxplot(width = 0.1, outlier.shape = NA, fill = "white") +
  scale_x_discrete(labels = c("FALSE" = "Non-responsive", "TRUE" = "Responsive")) +
  scale_fill_manual(values = c("FALSE" = "cornflowerblue", "TRUE" = "coral1")) +
  coord_cartesian(ylim = y_zoom) +
  labs(
    title = "C. Responsive genes are noisier",
    x = "Functional State", y = "Corrected Noise (VST)"
  ) +
  ggpubr::stat_compare_means(method = "wilcox.test", label = "p.signif") +
  theme(legend.position = "none")

panel_d <- ggplot(filtered_data, aes(x = has_TATA_box, y = variance.standardized, fill = has_TATA_box)) +
  geom_violin(trim = FALSE) +
  geom_boxplot(width = 0.1, outlier.shape = NA, fill = "white") +
  scale_x_discrete(labels = c("FALSE" = "TATA-less", "TRUE" = "TATA-containing")) +
  scale_fill_manual(values = c("FALSE" = "grey70", "TRUE" = "#F0E442")) +
  coord_cartesian(ylim = y_zoom) +
  labs(
    title = "D. TATA-containing genes are noisier",
    x = "Promoter Architecture", y = "Corrected Noise (VST)"
  ) +
  ggpubr::stat_compare_means(method = "wilcox.test", label = "p.signif") +
  theme(legend.position = "none")

# Aggregate for model visualization (Panel E)
model_data <- filtered_data %>%
  group_by(gene, cahn_group, h2az_group, has_TATA_box, is_responsive, avg_log2FC) %>%
  summarise(
    variance.standardized = mean(variance.standardized, na.rm = TRUE),
    log10_mean_expr = log10(mean(mean_expr, na.rm = TRUE)),
    .groups = "drop"
  ) %>%
  filter(is.finite(log10_mean_expr))

panel_e <- ggplot(model_data, aes(x = avg_log2FC, y = variance.standardized, color = cahn_group)) +
  geom_smooth(method = "lm", se = TRUE, aes(fill = cahn_group), alpha = 0.1) +
  scale_color_manual(values = c("Unmethylated" = "grey50", "gbM" = "#0072B2", "TE-like" = "#D55E00"), name = "Methylation:") +
  scale_fill_manual(values = c("Unmethylated" = "grey50", "gbM" = "#0072B2", "TE-like" = "#D55E00"), name = "Methylation:") +
  labs(
    title = "E. gbM Dampens Noise Primarily in Non-Responsive Genes",
    subtitle = "Visualizing the interaction between methylation and responsiveness",
    x = "Gene Responsiveness (avg_log2FC)",
    y = "Average Corrected Noise (VST)"
  ) +
  publication_theme

final_figure <- (panel_a | panel_d) / (panel_c | panel_b) / (panel_e) +
  plot_annotation(title = "Global Drivers of Gene Expression Noise", tag_levels = "A")
ggsave(final_figure_path, plot = final_figure, width = 12, height = 15, dpi = 300)
message(paste("Global summary figure saved to:", final_figure_path))

# --- 4. GENERATE CELL-TYPE-SPECIFIC FACETED FIGURE ---
message("Generating cell-type-specific faceted summary figure...")

panel_a_facet <- ggplot(filtered_data, aes(x = cahn_group, y = variance.standardized, fill = cahn_group)) +
  geom_violin(trim = FALSE) +
  scale_fill_manual(values = c("gbM" = "#0072B2", "TE-like" = "#D55E00", "Unmethylated" = "grey")) +
  coord_cartesian(ylim = y_zoom) +
  labs(x = "Methylation Status", y = "Corrected Noise (VST)") +
  facet_wrap(~base_cell_type) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

panel_d_facet <- ggplot(filtered_data, aes(x = has_TATA_box, y = variance.standardized, fill = has_TATA_box)) +
  geom_violin(trim = FALSE) +
  scale_fill_manual(values = c("FALSE" = "grey70", "TRUE" = "#F0E442")) +
  coord_cartesian(ylim = y_zoom) +
  labs(x = "Promoter Architecture", y = "Corrected Noise (VST)") +
  facet_wrap(~base_cell_type) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

celltype_figure <- (panel_a_facet / panel_d_facet) +
  plot_annotation(
    title = "Noise Drivers by Cell Type",
    caption = "Distribution of noise across epigenetic and promoter groups for each cell type."
  )
ggsave(celltype_figure_path, plot = celltype_figure, width = 16, height = 12, dpi = 300)
message(paste("Cell-type faceted figure saved to:", celltype_figure_path))

# --- 5. GLOBAL MULTIVARIATE MODELING ---
message("Building and summarizing global multivariate models...")

filtered_data$cahn_group <- factor(filtered_data$cahn_group, levels = c("Unmethylated", "gbM", "TE-like"))

mixed_effects_model <- lmer(
  variance.standardized ~ I(log10(mean_expr)) + has_TATA_box + cahn_group * avg_log2FC + (1 | gene),
  data = filtered_data
)

sink(model_summary_path)
cat("=================================================================\n"); cat("       GLOBAL STATISTICAL SUMMARY: 06_multivariate_modeling\n"); cat("=================================================================\n\n")
cat(paste("Summary generated on:", Sys.Date()), "\n\n")
cat("--- LINEAR MIXED-EFFECTS MODEL (LMM) - PLATINUM STANDARD ---\n")
cat("Formula: lmer(Noise ~ log10(MeanExpr) + TATA + Methylation * Responsiveness + (1 | gene))\n\n")
print(summary(mixed_effects_model))
sink()
message(paste("Global model summary saved to:", model_summary_path))

# --- 6. CELL-TYPE-SPECIFIC MULTIVARIATE MODELING ---
message("Building and summarizing cell-type-specific multivariate models...")
unique_cell_types <- unique(filtered_data$base_cell_type)

safe_lmer <- purrr::safely(lmer)

sink(celltype_model_summary_path)
cat("=================================================================\n"); cat("   CELL-TYPE-SPECIFIC MODEL SUMMARY: 06_multivariate_modeling\n"); cat("=================================================================\n\n")
cat(paste("Summary generated on:", Sys.Date()), "\n\n")
cat("This file contains a separate Linear Mixed-Effects Model for each cell type.\n")

for (cell_type in unique_cell_types) {
  message(paste("... Modeling cell type:", cell_type))
  ct_data <- filtered_data %>% filter(base_cell_type == cell_type)
  
  cat("\n\n-----------------------------------------------------------------\n")
  cat(paste("                MODEL FOR CELL TYPE:", toupper(cell_type)), "\n")
  cat("-----------------------------------------------------------------\n")
  
  if (nrow(ct_data) < 1000 || n_distinct(ct_data$gene) < 50) {
    cat("Skipping model due to insufficient data points or unique genes.\n")
    next
  }
  
  ct_model_fit <- safe_lmer(
    variance.standardized ~ I(log10(mean_expr)) + has_TATA_box + cahn_group * avg_log2FC + (1 | gene),
    data = ct_data
  )
  
  if (is.null(ct_model_fit$error)) {
    print(summary(ct_model_fit$result))
  } else {
    cat("Model failed to converge for this cell type.\n")
    cat("Error: ", ct_model_fit$error$message, "\n")
  }
}
sink()
message(paste("Cell-type-specific model summaries saved to:", celltype_model_summary_path))
message("Script 06 finished successfully.")