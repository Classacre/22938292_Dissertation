# ==============================================================================
# SCRIPT: 01_calculate_and_analyze_noise.R
#
# PURPOSE:
#   This script performs the foundational analysis of gene expression noise.
#   It corrects the critical methodological flaw of using a global coefficient
#   of variation (CV) by instead calculating noise (mean, variance, CV) for
#   each gene *within each homogenous cell type*. This approach disentangles
#   true stochastic noise from structured biological variation between cell types.
#
# MAJOR CORRECTIONS IMPLEMENTED:
#   1.  **Noise Definition:** All global CV calculations have been removed. Noise is
#       now exclusively calculated on a per-cell-type basis.
#   2.  **Data Structure:** The primary output is a long-format data frame where
#       each row is a gene-in-a-cell-type, facilitating robust statistical testing.
#   3.  **Statistical Analysis:** Tests are now performed on the distributions of
#       these valid within-cell-type CV values.
#   4.  **Modularity & Reproducibility:** The script is now structured to first
#       calculate and save the clean, processed noise data, which subsequent
#       scripts in the pipeline should use as their input.
#   5.  **Code Clarity:** The script is streamlined, with clear sections for
#       setup, calculation, analysis, and visualization.
#
# ==============================================================================

# --- 1. SETUP: Load libraries and define paths ---

# Ensure required packages are installed and loaded
packages_to_load <- c("ggplot2", "dplyr", "tidyr", "ggridges", "gridExtra", "grid")
for (pkg in packages_to_load) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, repos = "https://cloud.r-project.org/")
  }
  library(pkg, character.only = TRUE)
}

# --- Define I/O Paths ---
# Using variables for paths makes the script more portable and easier to manage.
# INPUTS
expr_matrix_path <- "/group/sms029/mnieuwenh/gbm_noise_analysis/expression_matrix.csv"
gene_anno_path <- "/group/sms029/mnieuwenh/Methylation_Data/combined_methylation_data.csv"
cell_meta_path <- "/group/sms029/mnieuwenh/seurat_metadata/seurat_metadata_full.csv"

# OUTPUTS
output_dir <- "/group/sms029/mnieuwenh/gbm_noise_analysis/results/01_noise_analysis"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# This is the most important output file for downstream scripts
processed_noise_data_path <- file.path(output_dir, "within_celltype_noise_data.csv")


# --- 2. DATA PREPARATION: Load and preprocess all required data ---

# Load expression matrix from Seurat's log-normalized data slot
message("Loading and transforming expression matrix...")
expr_mat_log <- as.matrix(read.csv(expr_matrix_path, row.names = 1, check.names = FALSE))

# CRITICAL: Reverse the log1p transformation to get linear-scale data for CV calculation
expr_mat <- expm1(expr_mat_log)
rm(expr_mat_log) # Free up memory

# Load and prepare cell metadata
message("Loading cell metadata...")
cell_metadata <- read.csv(cell_meta_path, row.names = 1, check.names = FALSE)

# Load and prepare master gene annotation table
message("Loading and preparing gene annotations...")
gene_anno_raw <- read.csv(gene_anno_path, stringsAsFactors = FALSE)

# Create a clean, unified gene annotation data frame
gene_annotations <- gene_anno_raw %>%
  mutate(
    # Cahn Group
    cahn_group = case_when(
      Cahn_Methylation_status == "gbM" ~ "gbM",
      Cahn_Methylation_status == "TE-like methylation" ~ "TE-like",
      Cahn_Methylation_status == "Unmethylated" ~ "Unmethylated",
      TRUE ~ NA_character_
    ),
    # Bewick Group
    bewick_group = case_when(
      Bewick_Classification == "gbM" ~ "gbM",
      Bewick_Classification %in% c("mCHH", "mCHG") ~ "TE-like", # Combine for easier comparison
      Bewick_Classification == "Unmethylated" ~ "Unmethylated",
      TRUE ~ NA_character_
    ),
    # H2A.Z Group
    h2az_group = case_when(
      H2AZ_Depleted == TRUE ~ "H2A.Z-Depleted",
      H2AZ_Enriched == TRUE ~ "H2A.Z-Enriched",
      TRUE ~ "Other"
    )
  ) %>%
  # Select only the necessary columns for merging later
  select(gene = Gene_ID, cahn_group, bewick_group, h2az_group) %>%
  # Ensure gene IDs are unique
  distinct(gene, .keep_all = TRUE)

# Ensure consistency between expression matrix and metadata
common_cells <- intersect(colnames(expr_mat), rownames(cell_metadata))
expr_mat <- expr_mat[, common_cells]
cell_metadata <- cell_metadata[common_cells, ]

message("Data preparation complete.")


# --- 3. CORE CALCULATION: Calculate noise within each cell type ---

message("Calculating noise within each cell type. This may take a while...")
cell_types <- unique(cell_metadata$identity)
MIN_CELLS_FOR_NOISE_CALC <- 20 # A reasonable threshold for robust variance calculation

# Use lapply to iterate over cell types and collect results in a list
noise_data_list <- lapply(cell_types, function(ct) {
  # Skip if cell type is NA or empty string
  if (is.na(ct) || ct == "") return(NULL)

  # Get cells belonging to the current cell type
  cells_in_ct <- rownames(cell_metadata)[cell_metadata$identity == ct]

  # Proceed only if there are enough cells
  if (length(cells_in_ct) < MIN_CELLS_FOR_NOISE_CALC) {
    message(paste("Skipping cell type '", ct, "' - only", length(cells_in_ct), "cells."))
    return(NULL)
  }

  # Subset the expression matrix
  expr_mat_ct <- expr_mat[, cells_in_ct, drop = FALSE]

  # Calculate mean, variance, and CV for each gene within this cell type
  mean_expr_ct <- rowMeans(expr_mat_ct)
  var_expr_ct <- apply(expr_mat_ct, 1, var)
  cv_expr_ct <- sqrt(var_expr_ct) / mean_expr_ct

  # Create a data frame for this cell type's results
  data.frame(
    gene = rownames(expr_mat_ct),
    cell_type = ct,
    n_cells = length(cells_in_ct),
    mean_expr = mean_expr_ct,
    cv_expr = cv_expr_ct,
    stringsAsFactors = FALSE
  )
})

# Combine the list of data frames into a single long-format data frame
noise_data_long <- bind_rows(noise_data_list)

# Merge with gene annotations
noise_data_long <- noise_data_long %>%
  left_join(gene_annotations, by = "gene") %>%
  # Filter out rows with infinite or NA CV values, which are uninformative
  filter(is.finite(cv_expr))

message(paste("Noise calculation complete. Processed", nrow(noise_data_long), "gene-in-cell-type observations."))


# --- 4. SAVE PROCESSED DATA: Export the clean noise data for other scripts ---

message(paste("Saving processed noise data to:", processed_noise_data_path))
write.csv(noise_data_long, processed_noise_data_path, row.names = FALSE)
message("This file should be the starting point for all subsequent noise-related analyses.")


# --- 5. OVERALL NOISE ANALYSIS: Compare noise distributions across all cell types ---

message("Performing overall noise analysis...")

# Create a filtered dataset for robust plotting and stats
# Here, we filter based on the mean expression *within each cell type*
noise_data_filtered <- noise_data_long %>%
  filter(mean_expr > 0.1) # A sensible baseline filter to remove unexpressed genes

# --- Visualization ---
plot_overall_noise <- function(data, group_col, title) {
  # Remove NA values for the specified group column for cleaner plots
  plot_data <- data %>% filter(!is.na(.data[[group_col]]))

  ggplot(plot_data, aes(x = .data[[group_col]], y = cv_expr, fill = .data[[group_col]])) +
    geom_violin(trim = FALSE, scale = "width", alpha = 0.8) +
    geom_boxplot(width = 0.1, outlier.size = 0.5, fill = "white") +
    scale_y_log10() + # Log scale is essential for visualizing CV
    labs(
      title = title,
      subtitle = "Comparing within-cell-type CV across all valid cell types",
      x = "Gene Group",
      y = "Coefficient of Variation (CV) (log10 scale)"
    ) +
    theme_bw(base_size = 14) +
    theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust = 1))
}

p_cahn_overall <- plot_overall_noise(noise_data_filtered, "cahn_group", "Overall Noise Comparison by Cahn Classification")
p_bewick_overall <- plot_overall_noise(noise_data_filtered, "bewick_group", "Overall Noise Comparison by Bewick Classification")
p_h2az_overall <- plot_overall_noise(noise_data_filtered, "h2az_group", "Overall Noise Comparison by H2A.Z Occupancy")

# Save the primary summary plots
ggsave(file.path(output_dir, "overall_noise_cahn.png"), plot = p_cahn_overall, width = 8, height = 7)
ggsave(file.path(output_dir, "overall_noise_bewick.png"), plot = p_bewick_overall, width = 8, height = 7)
ggsave(file.path(output_dir, "overall_noise_h2az.png"), plot = p_h2az_overall, width = 8, height = 7)


# --- Statistical Testing ---
run_overall_stats <- function(data, group_col) {
  # Filter data to ensure we have valid groups to compare
  stat_data <- data %>% filter(!is.na(.data[[group_col]]))
  
  # Ensure there are at least two groups to compare
  if (n_distinct(stat_data[[group_col]]) < 2) return(NULL)

  # Kruskal-Wallis test for an overall difference
  kruskal_test <- kruskal.test(cv_expr ~ get(group_col), data = stat_data)

  # Pairwise Wilcoxon tests for specific comparisons
  pairwise_tests <- pairwise.wilcox.test(stat_data$cv_expr, stat_data[[group_col]], p.adjust.method = "BH")

  return(list(kruskal = kruskal_test, pairwise = pairwise_tests))
}

stats_cahn <- run_overall_stats(noise_data_filtered, "cahn_group")
stats_bewick <- run_overall_stats(noise_data_filtered, "bewick_group")
stats_h2az <- run_overall_stats(noise_data_filtered, "h2az_group")

# Save statistical results to a text file
sink(file.path(output_dir, "overall_noise_statistics.txt"))
cat("====================================================\n")
cat(" OVERALL NOISE ANALYSIS - STATISTICAL SUMMARY\n")
cat("====================================================\n\n")
cat("NOTE: All tests performed on within-cell-type CV values for genes with mean_expr > 0.1 in that cell type.\n\n")

cat("--- Cahn Classification ---\n")
print(stats_cahn$kruskal)
cat("\nPairwise comparisons (BH-adjusted p-values):\n")
print(stats_cahn$pairwise)

cat("\n\n--- Bewick Classification ---\n")
print(stats_bewick$kruskal)
cat("\nPairwise comparisons (BH-adjusted p-values):\n")
print(stats_bewick$pairwise)

cat("\n\n--- H2A.Z Occupancy ---\n")
print(stats_h2az$kruskal)
cat("\nPairwise comparisons (BH-adjusted p-values):\n")
print(stats_h2az$pairwise)
sink()

message("Overall noise analysis complete.")


# --- 6. CELL-TYPE-SPECIFIC ANALYSIS: Visualize noise patterns in each cell type ---

message("Generating multi-page PDF with cell-type-specific noise plots...")

# Function to create a list of plots, one for each cell type
generate_celltype_plots <- function(data, group_col) {
  plot_list <- lapply(unique(data$cell_type), function(ct) {
    plot_data_ct <- data %>% filter(cell_type == ct, !is.na(.data[[group_col]]))

    # Skip if not enough data to plot
    if (nrow(plot_data_ct) < 10 || n_distinct(plot_data_ct[[group_col]]) < 2) return(NULL)

    ggplot(plot_data_ct, aes(x = .data[[group_col]], y = cv_expr, fill = .data[[group_col]])) +
      geom_boxplot(outlier.size = 0.5) +
      scale_y_log10() +
      labs(
        title = paste("Noise in:", ct),
        x = NULL, y = "CV (log10 scale)"
      ) +
      theme_bw(base_size = 10) +
      theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust = 1))
  })
  # Remove NULL entries from the list
  plot_list[!sapply(plot_list, is.null)]
}

# Generate plots for Cahn and Bewick classifications
plots_cahn_by_celltype <- generate_celltype_plots(noise_data_filtered, "cahn_group")
plots_bewick_by_celltype <- generate_celltype_plots(noise_data_filtered, "bewick_group")

# Save plots to multi-page PDFs
pdf(file.path(output_dir, "celltype_specific_noise_cahn.pdf"), width = 11, height = 8.5)
marrangeGrob(plots_cahn_by_celltype, nrow = 2, ncol = 3, top = "Within-Celltype Noise by Cahn Classification")
dev.off()

pdf(file.path(output_dir, "celltype_specific_noise_bewick.pdf"), width = 11, height = 8.5)
marrangeGrob(plots_bewick_by_celltype, nrow = 2, ncol = 3, top = "Within-Celltype Noise by Bewick Classification")
dev.off()

message("Cell-type-specific analysis complete.")
message("Script finished successfully.")