# ==============================================================================
# SCRIPT: 02_analyze_responsive_genes.R
#
# PURPOSE:
#   This script identifies cell-type-specific ("responsive") genes using a
#   statistically robust differential expression (DE) approach. It then
#   characterizes these gene sets by analyzing their epigenetic features and
#   their intrinsic expression noise, using the corrected noise data from script 01.
#
# METHODOLOGY:
#   1. Load the full Seurat object to run DE analysis.
#   2. Use FindAllMarkers to statistically define responsive genes.
#   3. Load the comprehensive, corrected noise data from script 01.
#   4. Compare the CORRECTED noise levels of responsive vs. non-responsive genes.
#   5. Perform Fisher's exact tests to check for epigenetic enrichment.
# ==============================================================================

# --- 1. SETUP: Load libraries and define paths ---

# Ensure required packages are installed and loaded
packages_to_load <- c("Seurat", "dplyr", "tidyr", "ggplot2", "ComplexUpset", "pheatmap")
for (pkg in packages_to_load) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, repos = "https://cloud.r-project.org/")
  }
  library(pkg, character.only = TRUE)
}

# --- Define I/O Paths ---
# INPUTS
seurat_object_path <- "/group/sms029/Oliva_dataset/integrated_col_trajectories_colonly.rds"
# CRITICAL: Use the output from the new script 01
processed_noise_data_path <- "/group/sms029/mnieuwenh/gbm_noise_analysis/results/01_noise_analysis/noise_analysis_complete_data.csv"

# OUTPUTS
output_dir <- "/group/sms029/mnieuwenh/gbm_noise_analysis/results/02_responsive_gene_analysis"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# --- 2. IDENTIFY RESPONSIVE GENES (CELL-TYPE MARKERS) ---

message("Loading Seurat object to identify cell-type markers...")
seurat_obj <- readRDS(seurat_object_path)
DefaultAssay(seurat_obj) <- "RNA"
Idents(seurat_obj) <- "identity" # Set the main identity for marker finding

message("Running FindAllMarkers to define responsive genes. This may take some time...")
# This is the standard, statistically robust method for finding cell-type-specific genes.
# `only.pos = TRUE` finds genes that are upregulated in a cell type compared to all others.
all_markers <- FindAllMarkers(seurat_obj, only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0.25)
rm(seurat_obj); gc() # Free up memory

# Save the full list of markers for reference
write.csv(all_markers, file.path(output_dir, "all_celltype_markers_full_list.csv"), row.names = FALSE)

# Define "responsive" genes as those that are statistically significant (adjusted p-value < 0.05)
responsive_genes_df <- all_markers %>%
  filter(p_val_adj < 0.05)

# Create a master list of unique responsive gene IDs
responsive_gene_list <- unique(responsive_genes_df$gene)
message(paste("Identified", length(responsive_gene_list), "unique responsive genes across all cell types."))
write.csv(data.frame(gene = responsive_gene_list), file.path(output_dir, "responsive_genes_list.csv"), row.names = FALSE)


# --- 3. CHARACTERIZE NOISE OF RESPONSIVE VS. NON-RESPONSIVE GENES ---

message("Loading pre-calculated and corrected noise data...")
# This file contains the 'variance.standardized' score
complete_noise_data <- read.csv(processed_noise_data_path)

# Add a column indicating if a gene is responsive or not
complete_noise_data <- complete_noise_data %>%
  mutate(gene_type = ifelse(gene %in% responsive_gene_list, "Responsive", "Non-Responsive"))

# Apply the same robust filtering as in script 01
# CORRECTED: Filter using 'variance.standardized'
noise_data_filtered <- complete_noise_data %>%
  filter(
    mean_expr > 0.01,
    pct_expressing > 0.1,
    is.finite(variance.standardized)
  )

# --- Visualization using the CORRECTED noise metric ---
# CORRECTED: Use 'variance.standardized'
p_noise_comparison_corrected <- ggplot(noise_data_filtered, aes(x = gene_type, y = variance.standardized, fill = gene_type)) +
  geom_violin(trim = FALSE, scale = "width", alpha = 0.8) +
  geom_boxplot(width = 0.1, outlier.size = 0.5, fill = "white") +
  labs(
    title = "Corrected Noise of Responsive vs. Non-Responsive Genes",
    subtitle = "Noise metric is standardized variance, calculated within each cell type",
    x = "Gene Type",
    y = "Corrected Noise (Standardized Variance)"
  ) +
  theme_bw(base_size = 14) +
  theme(legend.position = "none")

ggsave(file.path(output_dir, "noise_responsive_vs_nonresponsive_corrected.png"), plot = p_noise_comparison_corrected, width = 8, height = 7)

# --- Statistical Test using the CORRECTED noise metric ---
# CORRECTED: Use 'variance.standardized'
stat_noise_resp_corrected <- wilcox.test(variance.standardized ~ gene_type, data = noise_data_filtered)


# --- 4. EPIGENETIC ENRICHMENT ANALYSIS ---

message("Performing epigenetic enrichment analysis for responsive genes...")

# The "universe" or "background" is all genes that were successfully processed in script 01
background_genes <- unique(complete_noise_data$gene)

# Function to perform and store Fisher's exact test results
perform_enrichment_test <- function(feature_group) {
  # Get the set of genes belonging to the feature group (e.g., all 'gbM' genes)
  feature_genes <- complete_noise_data %>%
    filter(!is.na(.data[[feature_group]]) & .data[[feature_group]] != "Other") %>%
    distinct(gene, .data[[feature_group]])

  # Iterate over each category within the feature (e.g., 'gbM', 'Unmethylated')
  lapply(unique(feature_genes[[feature_group]]), function(category) {
    genes_in_category <- feature_genes$gene[feature_genes[[feature_group]] == category]

    # Create the 2x2 contingency table for Fisher's Exact Test
    a <- length(intersect(responsive_gene_list, genes_in_category))
    b <- length(genes_in_category) - a
    c <- length(responsive_gene_list) - a
    d <- length(background_genes) - a - b - c

    fisher_test <- fisher.test(matrix(c(a, c, b, d), nrow = 2), alternative = "greater")

    data.frame(
      feature_group = feature_group,
      category = category,
      odds_ratio = fisher_test$estimate,
      p_value = fisher_test$p.value
    )
  }) %>% bind_rows()
}

# Run tests for each epigenetic feature
enrichment_cahn <- perform_enrichment_test("cahn_group")
enrichment_bewick <- perform_enrichment_test("bewick_group")
enrichment_h2az <- perform_enrichment_test("h2az_group")

# Combine all results and adjust p-values for multiple testing
final_enrichment_df <- bind_rows(enrichment_cahn, enrichment_bewick, enrichment_h2az)
final_enrichment_df$p_adj <- p.adjust(final_enrichment_df$p_value, method = "BH")


# --- 5. SAVE ALL RESULTS ---

message("Saving all analysis outputs...")
write.csv(final_enrichment_df, file.path(output_dir, "responsive_gene_epigenetic_enrichment.csv"), row.names = FALSE)

sink(file.path(output_dir, "responsive_analysis_summary.txt"))
cat("====================================================\n")
cat(" RESPONSIVE GENE ANALYSIS - SUMMARY\n")
cat("====================================================\n\n")
cat("1. Definition of Responsive Genes:\n")
cat("   - Based on Seurat::FindAllMarkers (only.pos=T, min.pct=0.1, logfc.threshold=0.25)\n")
cat(paste("   - Genes with p_val_adj < 0.05 were considered responsive.\n"))
cat(paste("   - Found", length(responsive_gene_list), "unique responsive genes in total.\n\n"))
cat("----------------------------------------------------\n\n")
cat("2. Noise Comparison (Responsive vs. Non-Responsive):\n")
cat("   - CRITICAL: Test performed on the CORRECTED noise metric (standardized variance).\n")
cat("   - Data was filtered for mean_expr > 0.01 and pct_expressing > 0.1.\n")
print(stat_noise_resp_corrected)
cat("\n----------------------------------------------------\n\n")
cat("3. Epigenetic Enrichment in Responsive Genes (Fisher's Exact Test):\n")
cat("   - 'p_adj' is the Benjamini-Hochberg corrected p-value across all tests.\n\n")
print(final_enrichment_df)
sink()


# --- 6. VISUALIZATION OF OVERLAPS ---

message("Generating UpSet plot for overlaps...")

# Create a boolean data frame for the UpSet plot
upset_data <- complete_noise_data %>%
  distinct(gene, .keep_all = TRUE) %>%
  mutate(
    Is_Responsive = gene %in% responsive_gene_list,
    Is_gbM_Cahn = cahn_group == 'gbM',
    Is_H2AZ_Depleted = h2az_group == 'H2A.Z-Depleted'
  ) %>%
  select(Is_Responsive, Is_gbM_Cahn, Is_H2AZ_Depleted) %>%
  mutate_all(~ifelse(is.na(.), FALSE, .)) # Replace NA with FALSE for plotting

p_upset <- ComplexUpset::upset(
  upset_data,
  colnames(upset_data),
  name = "Gene Sets",
  width_ratio = 0.1,
  min_size = 10 # Only show intersections with at least 10 genes
)

ggsave(file.path(output_dir, "upset_responsive_gbm_h2az.png"), plot = p_upset, width = 10, height = 7)

message("Script finished successfully.")