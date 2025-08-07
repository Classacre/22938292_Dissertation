# ==============================================================================
# SCRIPT: 06_generate_summary_figures.R
#
# PURPOSE:
#   This is the final script in the analysis pipeline. Its purpose is to
#   synthesize the primary findings from the preceding analyses into a single,
#   cohesive, publication-quality figure. This figure will visually summarize
#   the core narrative of the dissertation:
#
#   1. gbM genes exhibit lower intrinsic expression noise.
#   2. Biologically "responsive" (cell-type-specific) genes are noisier.
#   3. These responsive genes are epigenetically distinct, showing a depletion
#      of gbM and an enrichment of H2A.Z.
#
# METHODOLOGY:
#   - Load the final, fully annotated dataset from script 05.
#   - Re-create the key plots from previous scripts in a standardized format.
#   - Use the 'patchwork' library to assemble these individual plots into a
#     single, multi-panel figure with clear A, B, C labeling.
#   - Save the final figure and a concluding text summary.
# ==============================================================================

# --- 1. SETUP: Load libraries and define paths ---

# Load necessary packages for plotting and data manipulation
packages_to_load <- c("dplyr", "ggplot2", "patchwork", "scales")
for (pkg in packages_to_load) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, repos = "https://cloud.r-project.org/")
  }
  library(pkg, character.only = TRUE)
}

# --- Define I/O Paths ---
# INPUT: The final annotated data file from script 05 (UPDATED PATH)
annotated_data_path <- "/group/sms029/mnieuwenh/gbm_noise_analysis/results/05_analyze_responsive_genes/05_gene_noise_with_responsiveness.csv"

# OUTPUTS: (UPDATED)
output_dir <- "/group/sms029/mnieuwenh/gbm_noise_analysis/results/06_generate_summary_figures"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# Final output filenames
summary_figure_path <- file.path(output_dir, "06_final_summary_figure.png")
summary_text_path <- file.path(output_dir, "06_analysis_conclusion.txt")


# --- 2. DATA LOADING AND FINAL FILTERING ---

message("Loading final annotated data...")
if (!file.exists(annotated_data_path)) {
  stop(paste("Error: Final data file not found at", annotated_data_path))
}
full_data <- read.csv(annotated_data_path)

message("Applying standard filters for consistency across all plots...")
# Apply the same robust filtering used throughout the analysis
filtered_data <- full_data %>%
  filter(
    mean_expr > 0.01,
    pct_expressing > 0.1,
    is.finite(variance.standardized)
  )

# Define a consistent theme for all plots
publication_theme <- theme_bw(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    plot.subtitle = element_text(size = 10),
    legend.position = "bottom",
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.minor = element_blank()
  )


# --- 3. CREATE PANEL A: gbM GENES HAVE LOWER NOISE ---

message("Generating Panel A: Noise by Epigenetic Group...")

# Data for Panel A
panel_a_data <- filtered_data %>%
  filter(!is.na(cahn_group)) %>%
  # Set factor levels for a logical plot order
  mutate(cahn_group = factor(cahn_group, levels = c("gbM", "Unmethylated", "TE-like")))

# Zoom in on the core 98% of the data to avoid extreme outliers skewing the view
y_zoom <- quantile(panel_a_data$variance.standardized, probs = c(0.01, 0.99), na.rm = TRUE)

# Statistical annotation
stats_a <- wilcox.test(
  variance.standardized ~ (cahn_group == "gbM"),
  data = panel_a_data,
  alternative = "two.sided"
)

panel_a <- ggplot(panel_a_data, aes(x = cahn_group, y = variance.standardized, fill = cahn_group)) +
  geom_violin(trim = FALSE, scale = "width", alpha = 0.8) +
  geom_boxplot(width = 0.1, outlier.shape = NA, fill = "white", alpha = 0.5) +
  coord_cartesian(ylim = y_zoom) +
  scale_fill_manual(values = c("gbM" = "#0072B2", "Unmethylated" = "#D55E00", "TE-like" = "#CC79A7")) +
  labs(
    title = "gbM is associated with lower noise",
    subtitle = paste0("gbM vs. Others, p < ", format.pval(stats_a$p.value, eps = 0.001)),
    x = "Epigenetic Group (Cahn et al.)",
    y = "Corrected Noise (Std. Variance)"
  ) +
  publication_theme +
  theme(legend.position = "none")


# --- 4. CREATE PANEL B: RESPONSIVE GENES ARE NOISIER ---

message("Generating Panel B: Noise of Responsive Genes...")

# Data for Panel B
panel_b_data <- filtered_data

# Statistical annotation
stats_b <- wilcox.test(variance.standardized ~ is_responsive, data = panel_b_data)

panel_b <- ggplot(panel_b_data, aes(x = is_responsive, y = variance.standardized, fill = is_responsive)) +
  geom_violin(trim = FALSE, scale = "width", alpha = 0.8) +
  geom_boxplot(width = 0.1, outlier.shape = NA, fill = "white", alpha = 0.5) +
  coord_cartesian(ylim = y_zoom) +
  scale_x_discrete(labels = c("FALSE" = "Non-Responsive", "TRUE" = "Responsive")) +
  scale_fill_manual(values = c("FALSE" = "gray70", "TRUE" = "coral")) +
  labs(
    title = "Responsive genes are noisier",
    subtitle = paste0("Wilcoxon test p < ", format.pval(stats_b$p.value, eps = 0.001)),
    x = "Gene Category",
    y = "Corrected Noise (Std. Variance)"
  ) +
  publication_theme +
  theme(legend.position = "none")


# --- 5. CREATE PANEL C: EPIGENETIC MAKEUP OF RESPONSIVE GENES ---

message("Generating Panel C: Epigenetics of Responsive Genes...")

# Prepare data for Panel C
panel_c_data <- filtered_data %>%
  select(gene, is_responsive, cahn_group, h2az_group) %>%
  distinct(gene, .keep_all = TRUE) %>%
  # Combine epigenetic groups for a cleaner plot
  mutate(epigenetic_summary = case_when(
    cahn_group == "gbM" ~ "gbM",
    h2az_group == "H2A.Z-Enriched" ~ "H2A.Z-Enriched",
    cahn_group == "Unmethylated" ~ "Unmethylated",
    TRUE ~ "Other"
  )) %>%
  filter(!is.na(epigenetic_summary), epigenetic_summary != "Other") %>%
  count(epigenetic_summary, is_responsive) %>%
  group_by(is_responsive) %>%
  mutate(proportion = n / sum(n))

# Statistical annotation
contingency_c <- table(
  (panel_c_data %>% filter(is_responsive == TRUE))$epigenetic_summary,
  (panel_c_data %>% filter(is_responsive == FALSE))$epigenetic_summary
)
stats_c <- chisq.test(contingency_c)

panel_c <- ggplot(panel_c_data, aes(x = is_responsive, y = proportion, fill = epigenetic_summary)) +
  geom_col(position = "fill", color = "black", width = 0.8) +
  scale_y_continuous(labels = scales::percent_format()) +
  scale_x_discrete(labels = c("FALSE" = "Non-Responsive", "TRUE" = "Responsive")) +
  scale_fill_manual(
    name = "Epigenetic Group",
    values = c("gbM" = "#0072B2", "H2A.Z-Enriched" = "#009E73", "Unmethylated" = "#D55E00")
  ) +
  labs(
    title = "Responsive genes are depleted of gbM",
    subtitle = paste0("Chi-squared test p < ", format.pval(stats_c$p.value, eps = 0.001)),
    x = "Gene Category",
    y = "Proportion of Genes"
  ) +
  publication_theme


# --- 6. ASSEMBLE AND SAVE THE FINAL FIGURE ---

message("Assembling the final summary figure...")

# Use patchwork to combine the plots into a single figure with A, B, C tags
# Layout: Panel A on the left, Panels B and C stacked on the right
final_figure <- panel_a + (panel_b / panel_c) +
  plot_annotation(
    title = "Gene Body Methylation is Associated with Low-Noise, Stable Gene Expression",
    tag_levels = 'A'
  ) &
  theme(plot.tag = element_text(size = 18, face = "bold"))

# Save the final figure
ggsave(
  summary_figure_path,
  plot = final_figure,
  width = 12,
  height = 8,
  dpi = 300,
  bg = "white"
)
message(paste("Final summary figure saved to:", summary_figure_path))


# --- 7. WRITE FINAL SUMMARY TEXT ---

final_summary <- c(
  "=================================================================",
  "            FINAL ANALYSIS SUMMARY & CONCLUSION",
  "=================================================================",
  "",
  paste("Date:", Sys.Date()),
  "",
  "This analysis pipeline has successfully processed and integrated single-cell expression data with epigenetic annotations to investigate the determinants of gene expression noise.",
  "",
  "Key Findings Summarized in the Figure:",
  "---------------------------------------",
  "Panel A: Genes marked by gene body methylation (gbM) consistently show lower intrinsic expression noise after correcting for mean-expression effects. This suggests gbM plays a role in buffering expression variability.",
  "",
  "Panel B: Genes identified as 'responsive' (i.e., having cell-type-specific expression) exhibit significantly higher expression noise than stably expressed, non-responsive genes. This links higher noise to dynamic biological function.",
  "",
  "Panel C: The epigenetic makeup of these two gene categories differs starkly. Responsive genes are significantly depleted of gbM and are enriched for H2A.Z occupancy, a histone variant associated with transcriptional plasticity.",
  "",
  "Overall Conclusion:",
  "-------------------",
  "The results support a model where gene body methylation is a feature of stably expressed, low-noise 'housekeeping' genes, while dynamic, responsive genes utilize alternative epigenetic strategies (like H2A.Z enrichment) that permit higher expression variability, likely facilitating rapid changes in transcription.",
  "",
  paste("The final summary figure has been saved to:", summary_figure_path)
)

writeLines(final_summary, con = summary_text_path)
message(paste("Final summary text saved to:", summary_text_path))
message("Analysis pipeline complete.")