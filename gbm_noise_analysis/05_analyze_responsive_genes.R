# ==============================================================================
# SCRIPT: 05_analyze_responsive_genes.R
#
# PURPOSE:
#   To investigate the relationship between gene expression noise and biological
#   function, this script defines and analyzes "responsive" genes. We define
#   responsive genes as those showing statistically significant cell-type-specific
#   expression.
#
# METHODOLOGY:
#   1. Load the full Seurat object to perform differential expression analysis.
#   2. Use Seurat's `FindAllMarkers` function to identify genes significantly
#      upregulated in each cell type compared to all others. These are our
#      "responsive" genes for that cell type.
#   3. Load the complete noise data from script 03.
#   4. Integrate the responsiveness status into the main noise data table.
#   5. Compare the corrected noise profiles and epigenetic enrichments between
#      responsive genes and their non-responsive counterparts.
#   6. Perform statistical tests to quantify these differences.
#
# OUTPUTS:
#   - A table of all identified responsive genes.
#   - An updated master data file annotated with responsiveness status.
#   - Plots comparing noise and epigenetic features of responsive vs. non-responsive genes.
#   - A summary of all statistical tests.
# ==============================================================================

# --- 1. SETUP: Load libraries and define paths ---

# Load necessary packages
packages_to_load <- c("Seurat", "dplyr", "tidyr", "ggplot2", "patchwork")
for (pkg in packages_to_load) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, repos = "https://cloud.r-project.org/")
  }
  library(pkg, character.only = TRUE)
}

# --- Define I/O Paths ---
# INPUTS
# Full Seurat object is required for FindAllMarkers
seurat_object_path <- "/group/sms029/Oliva_dataset/integrated_col_trajectories_colonly.rds"
# Complete noise data from script 03 (UPDATED PATH)
complete_noise_data_path <- "/group/sms029/mnieuwenh/gbm_noise_analysis/results/03_calculate_expression_noise/03_gene_noise_metrics_complete.csv"

# OUTPUTS (UPDATED)
output_dir <- "/group/sms029/mnieuwenh/gbm_noise_analysis/results/05_analyze_responsive_genes"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# Consistent output filenames
responsive_markers_path <- file.path(output_dir, "05_responsive_gene_markers.csv")
data_with_responsiveness_path <- file.path(output_dir, "05_gene_noise_with_responsiveness.csv")
noise_plot_path <- file.path(output_dir, "05_noise_responsive_vs_nonresponsive.png")
epigenetics_plot_path <- file.path(output_dir, "05_epigenetics_of_responsive_genes.png")
stats_summary_path <- file.path(output_dir, "05_responsive_gene_statistics.txt")


# --- 2. IDENTIFY RESPONSIVE GENES ---

message("Loading Seurat object to identify responsive genes...")
seurat_obj <- readRDS(seurat_object_path)
DefaultAssay(seurat_obj) <- "RNA"
# Ensure cell identities are set correctly for marker finding
Idents(seurat_obj) <- "identity"

message("Running FindAllMarkers to find cell-type-specific (responsive) genes...")
# We use stringent criteria to define responsive genes:
# - logfc.threshold: Only test genes with at least a modest fold-change.
# - min.pct: Only test genes expressed in a minimum fraction of cells.
# - only.pos: We are interested in genes that are *upregulated* in a cell type.
all_markers <- FindAllMarkers(
  seurat_obj,
  only.pos = TRUE,
  min.pct = 0.1,
  logfc.threshold = 0.25
)
rm(seurat_obj); gc() # Free up memory

# Filter for significance and save the list of responsive genes
responsive_genes <- all_markers %>%
  filter(p_val_adj < 0.05)

message(paste("Identified", nrow(responsive_genes), "significant responsive gene markers across all cell types."))
write.csv(responsive_genes, responsive_markers_path, row.names = FALSE)


# --- 3. INTEGRATE RESPONSIVENESS STATUS WITH NOISE DATA ---

message("Loading complete noise data from script 03...")
noise_data <- read.csv(complete_noise_data_path)

message("Annotating main data table with responsiveness status...")
# Create a unique key for matching (gene + cell_type)
responsive_genes_key <- responsive_genes %>%
  mutate(key = paste(gene, cluster, sep = "_")) %>%
  pull(key)

# Add a new boolean column 'is_responsive' to the main noise data table.
# A gene is responsive if it was identified as a significant positive marker
# for that specific cell type.
noise_data_annotated <- noise_data %>%
  mutate(
    key = paste(gene, cell_type, sep = "_"),
    is_responsive = key %in% responsive_genes_key
  ) %>%
  select(-key) # Remove the temporary key

# Apply the same robust filtering as in previous scripts for analysis
filtered_data <- noise_data_annotated %>%
  filter(
    mean_expr > 0.01,
    pct_expressing > 0.1,
    is.finite(variance.standardized)
  )

# Save this newly annotated table, as it will be crucial for the final summary
write.csv(noise_data_annotated, data_with_responsiveness_path, row.names = FALSE)
message(paste("Annotated data saved to:", data_with_responsiveness_path))


# --- 4. ANALYSIS: COMPARE NOISE PROFILES ---

message("Analyzing noise differences between responsive and non-responsive genes...")

# Statistical test: Wilcoxon test to compare the two groups
noise_stat_test <- wilcox.test(
  variance.standardized ~ is_responsive,
  data = filtered_data,
  alternative = "two.sided"
)

# Visualization
y_zoom_noise <- quantile(filtered_data$variance.standardized, probs = c(0.01, 0.99), na.rm = TRUE)

p_noise <- ggplot(filtered_data, aes(x = is_responsive, y = variance.standardized, fill = is_responsive)) +
  geom_violin(trim = FALSE, scale = "width", alpha = 0.8) +
  geom_boxplot(width = 0.1, outlier.size = 0.5, fill = "white") +
  coord_cartesian(ylim = y_zoom_noise) +
  scale_x_discrete(labels = c("FALSE" = "Non-Responsive", "TRUE" = "Responsive")) +
  scale_fill_manual(values = c("FALSE" = "gray70", "TRUE" = "coral")) +
  labs(
    title = "Responsive Genes Exhibit Higher Expression Noise",
    subtitle = paste0("Wilcoxon test p-value: ", format.pval(noise_stat_test$p.value, digits = 2)),
    x = "Gene Category",
    y = "Corrected Noise (Standardized Variance)"
  ) +
  theme_bw(base_size = 14) +
  theme(legend.position = "none")

ggsave(noise_plot_path, plot = p_noise, width = 8, height = 7)
message(paste("Noise comparison plot saved to:", noise_plot_path))


# --- 5. ANALYSIS: CHARACTERIZE EPIGENETIC FEATURES ---

message("Analyzing epigenetic differences between responsive and non-responsive genes...")

# Prepare data for epigenetic analysis
# We need one entry per gene, so we'll focus on genes that are responsive in *at least one* cell type
epigenetic_data <- filtered_data %>%
  select(gene, is_responsive, cahn_group, bewick_group, h2az_group) %>%
  group_by(gene) %>%
  # A gene is considered 'responsive' overall if it's responsive in any cell type
  summarise(
    is_responsive = any(is_responsive),
    cahn_group = first(cahn_group),
    bewick_group = first(bewick_group),
    h2az_group = first(h2az_group)
  ) %>%
  ungroup()

# Function to create proportion plots and run Chi-squared tests
analyze_epigenetic_proportions <- function(data, group_col) {
  # Prepare data for plotting and testing
  plot_data <- data %>%
    filter(!is.na(.data[[group_col]])) %>%
    count(.data[[group_col]], is_responsive) %>%
    group_by(is_responsive) %>%
    mutate(proportion = n / sum(n)) %>%
    ungroup()

  # Chi-squared test for independence
  contingency_table <- table(data$is_responsive, data[[group_col]])
  chi_sq_test <- chisq.test(contingency_table)

  # Plot
  p <- ggplot(plot_data, aes(x = is_responsive, y = proportion, fill = .data[[group_col]])) +
    geom_col(position = "fill", color = "black") +
    scale_y_continuous(labels = scales::percent) +
    scale_x_discrete(labels = c("FALSE" = "Non-Responsive", "TRUE" = "Responsive")) +
    labs(
      subtitle = paste0("Chi-squared p: ", format.pval(chi_sq_test$p.value, digits = 2)),
      x = "Gene Category",
      y = "Proportion of Genes",
      fill = "Epigenetic Group"
    ) +
    theme_bw(base_size = 12) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

  return(list(plot = p, test = chi_sq_test, table = contingency_table))
}

# Run analysis for each epigenetic classification
cahn_analysis <- analyze_epigenetic_proportions(epigenetic_data, "cahn_group")
bewick_analysis <- analyze_epigenetic_proportions(epigenetic_data, "bewick_group")
h2az_analysis <- analyze_epigenetic_proportions(epigenetic_data, "h2az_group")

# Combine plots
p_epigenetics <- (cahn_analysis$plot + labs(title="Cahn")) +
                 (bewick_analysis$plot + labs(title="Bewick")) +
                 (h2az_analysis$plot + labs(title="H2A.Z")) +
                 plot_layout(guides = "collect") & theme(legend.position = "bottom")

p_epigenetics <- p_epigenetics + plot_annotation(title = "Epigenetic Composition of Responsive vs. Non-Responsive Genes")

ggsave(epigenetics_plot_path, plot = p_epigenetics, width = 12, height = 8)
message(paste("Epigenetic composition plot saved to:", epigenetics_plot_path))


# --- 6. SAVE STATISTICAL SUMMARY ---

sink(stats_summary_path)
cat("=================================================================\n")
cat(" STATISTICAL SUMMARY FOR RESPONSIVE GENE ANALYSIS\n")
cat("=================================================================\n\n")

cat("--- Noise Comparison (Corrected Noise) ---\n")
cat("Comparison of 'variance.standardized' between responsive and non-responsive genes.\n")
print(noise_stat_test)

cat("\n\n--- Epigenetic Enrichment (Cahn Classification) ---\n")
cat("Contingency Table:\n"); print(cahn_analysis$table)
cat("\nChi-squared Test:\n"); print(cahn_analysis$test)

cat("\n\n--- Epigenetic Enrichment (Bewick Classification) ---\n")
cat("Contingency Table:\n"); print(bewick_analysis$table)
cat("\nChi-squared Test:\n"); print(bewick_analysis$test)

cat("\n\n--- Epigenetic Enrichment (H2A.Z Occupancy) ---\n")
cat("Contingency Table:\n"); print(h2az_analysis$table)
cat("\nChi-squared Test:\n"); print(h2az_analysis$test)

sink()
message(paste("Statistical summary saved to:", stats_summary_path))
message("Script finished successfully.")