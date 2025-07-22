# ==============================================================================
# SCRIPT: 02_analyze_responsive_genes.R
#
# PURPOSE:
#   This script identifies cell-type-specific ("responsive") genes using a
#   statistically robust differential expression (DE) approach. It then
#   characterizes these gene sets by analyzing their epigenetic features and
#   their intrinsic expression noise, using the pre-calculated, valid noise
#   data from the previous script (01_calculate_and_analyze_noise.R).
#
# MAJOR CORRECTIONS IMPLEMENTED:
#   1.  **Statistical Definition of Responsive Genes:** Replaced the arbitrary
#       fold-change heuristic with a formal DE analysis (`Seurat::FindAllMarkers`)
#       to identify genes significantly upregulated in each cell type.
#   2.  **Correct Noise Data:** The script NO LONGER calculates noise. It loads the
#       validated `within_celltype_noise_data.csv` file, ensuring all noise
#       comparisons are statistically sound.
#   3.  **Modular & Efficient:** The script acts as a self-contained analysis
#       module, relying on pre-processed inputs for both noise and expression data.
#   4.  **Rigorous Enrichment Analysis:** Implements a systematic framework for
#       testing the enrichment of epigenetic features (gbM, H2A.Z) within the
#       newly defined responsive gene sets, including p-value adjustments.
#
# ==============================================================================

# --- 1. SETUP: Load libraries and define paths ---

# Ensure required packages are installed and loaded
packages_to_load <- c("Seurat", "ggplot2", "dplyr", "tidyr", "ComplexUpset", "pheatmap")
for (pkg in packages_to_load) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, repos = "https://cloud.r-project.org/")
  }
  library(pkg, character.only = TRUE)
}

# --- Define I/O Paths ---
# INPUTS
# The original Seurat object is needed for the DE analysis.
seurat_object_path <- "/group/sms029/Oliva_dataset/integrated_col_trajectories_colonly.rds"
# This is the CRITICAL input from the previous script.
processed_noise_data_path <- "/group/sms029/mnieuwenh/gbm_noise_analysis/results/01_noise_analysis/within_celltype_noise_data.csv"

# OUTPUTS
output_dir <- "/group/sms029/mnieuwenh/gbm_noise_analysis/results/02_responsive_gene_analysis"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# --- 2. IDENTIFY RESPONSIVE GENES (CELL-TYPE MARKERS) ---

message("Loading Seurat object to identify cell-type markers...")
# Note: This assumes the Seurat object has an 'identity' column in its metadata
# that matches the cell types used in the noise analysis.
seurat_obj <- readRDS(seurat_object_path)

# Set the main identity for marker finding
Idents(seurat_obj) <- "identity"

# Find markers for every cell type compared to all other cells.
# This is the standard, statistically robust method for finding cell-type-specific genes.
# We use a positive log-fold change threshold to find upregulated genes.
message("Running FindAllMarkers to define responsive genes. This may take some time...")
all_markers <- FindAllMarkers(seurat_obj, only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0.25)

# Save the full list of markers for reference
write.csv(all_markers, file.path(output_dir, "all_celltype_markers_full_list.csv"), row.names = FALSE)

# Define "responsive" genes as those that are statistically significant markers
# A common threshold is an adjusted p-value < 0.05.
responsive_genes_df <- all_markers %>%
  filter(p_val_adj < 0.05) %>%
  select(gene, responsive_cell_type = cluster, avg_log2FC, p_val_adj)

# Create a master list of unique responsive genes
responsive_gene_list <- unique(responsive_genes_df$gene)

message(paste("Identified", length(responsive_gene_list), "unique responsive genes across all cell types."))
write.csv(data.frame(gene = responsive_gene_list), file.path(output_dir, "responsive_genes_list.csv"), row.names = FALSE)


# --- 3. CHARACTERIZE NOISE OF RESPONSIVE VS. NON-RESPONSIVE GENES ---

message("Loading pre-calculated noise data...")
noise_data <- read.csv(processed_noise_data_path)

# Add a column to the noise data indicating if a gene is responsive or not
noise_data <- noise_data %>%
  mutate(gene_type = ifelse(gene %in% responsive_gene_list, "Responsive", "Non-Responsive"))

# Filter for robust comparison (e.g., mean expression > 0.1)
noise_data_filtered <- noise_data %>%
  filter(mean_expr > 0.1)

# --- Visualization ---
p_noise_comparison <- ggplot(noise_data_filtered, aes(x = gene_type, y = cv_expr, fill = gene_type)) +
  geom_violin(trim = FALSE, scale = "width", alpha = 0.8) +
  geom_boxplot(width = 0.1, outlier.size = 0.5, fill = "white") +
  scale_y_log10() +
  labs(
    title = "Expression Noise of Responsive vs. Non-Responsive Genes",
    subtitle = "Noise (CV) calculated within each cell type",
    x = "Gene Type",
    y = "Within-Celltype Coefficient of Variation (CV) (log10 scale)"
  ) +
  theme_bw(base_size = 14) +
  theme(legend.position = "none")

ggsave(file.path(output_dir, "noise_responsive_vs_nonresponsive.png"), plot = p_noise_comparison, width = 8, height = 7)

# --- Statistical Test ---
stat_noise_resp <- wilcox.test(cv_expr ~ gene_type, data = noise_data_filtered)


# --- 4. EPIGENETIC ENRICHMENT ANALYSIS ---

message("Performing epigenetic enrichment analysis for responsive genes...")

# We need a data frame of all genes included in the noise analysis as the "background" set
background_genes <- unique(noise_data$gene)

# Create a summary data frame for Fisher's exact tests
enrichment_results <- list()

# Function to perform and store Fisher's exact test results
perform_enrichment_test <- function(feature_group) {
  # Get the set of genes belonging to the feature group (e.g., all 'gbM' genes)
  feature_genes <- noise_data %>%
    filter(!is.na(.data[[feature_group]]) & .data[[feature_group]] != "Other") %>%
    distinct(gene, .data[[feature_group]])

  # Iterate over each category within the feature (e.g., 'gbM', 'Unmethylated')
  lapply(unique(feature_genes[[feature_group]]), function(category) {
    genes_in_category <- feature_genes$gene[feature_genes[[feature_group]] == category]

    # Create the 2x2 contingency table
    #               Is Responsive | Is Not Responsive
    # Is Feature    a             | b
    # Is Not Feature  c             | d
    a <- length(intersect(responsive_gene_list, genes_in_category))
    b <- length(genes_in_category) - a
    c <- length(responsive_gene_list) - a
    d <- length(background_genes) - a - b - c

    contingency_table <- matrix(c(a, c, b, d), nrow = 2)

    # Perform the test
    fisher_test <- fisher.test(contingency_table, alternative = "greater")

    data.frame(
      feature_group = feature_group,
      category = category,
      odds_ratio = fisher_test$estimate,
      p_value = fisher_test$p.value
    )
  })
}

# Run tests for each epigenetic feature
enrichment_results[["cahn"]] <- perform_enrichment_test("cahn_group")
enrichment_results[["bewick"]] <- perform_enrichment_test("bewick_group")
enrichment_results[["h2az"]] <- perform_enrichment_test("h2az_group")

# Combine all results into a single data frame
final_enrichment_df <- bind_rows(unlist(enrichment_results, recursive = FALSE))

# Adjust p-values for multiple testing across all tests performed
final_enrichment_df$p_adj <- p.adjust(final_enrichment_df$p_value, method = "BH")


# --- 5. SAVE ALL RESULTS ---

message("Saving all analysis outputs...")

# Save the enrichment table
write.csv(final_enrichment_df, file.path(output_dir, "responsive_gene_epigenetic_enrichment.csv"), row.names = FALSE)

# Save a summary text file
sink(file.path(output_dir, "responsive_analysis_summary.txt"))
cat("====================================================\n")
cat(" RESPONSIVE GENE ANALYSIS - SUMMARY\n")
cat("====================================================\n\n")
cat("1. Definition of Responsive Genes:\n")
cat("   - Based on FindAllMarkers (Seurat) with parameters: only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0.25\n")
cat(paste("   - Genes with p_val_adj < 0.05 were considered responsive.\n"))
cat(paste("   - Found", length(responsive_gene_list), "unique responsive genes in total.\n\n"))
cat("----------------------------------------------------\n\n")
cat("2. Noise Comparison (Responsive vs. Non-Responsive):\n")
cat("   - Test: Wilcoxon rank-sum test on within-cell-type CVs (mean_expr > 0.1).\n")
print(stat_noise_resp)
cat("\n----------------------------------------------------\n\n")
cat("3. Epigenetic Enrichment in Responsive Genes (Fisher's Exact Test):\n")
cat("   - 'p_adj' is the Benjamini-Hochberg corrected p-value across all tests.\n\n")
print(final_enrichment_df)
sink()

# --- 6. VISUALIZATION OF OVERLAPS ---

message("Generating UpSet plot for overlaps...")

# Create a boolean data frame for the UpSet plot
upset_data <- noise_data %>%
  distinct(gene, .keep_all = TRUE) %>%
  mutate(
    Is_Responsive = gene %in% responsive_gene_list,
    Is_gbM_Cahn = cahn_group == "gbM",
    Is_H2AZ_Depleted = h2az_group == "H2A.Z-Depleted"
  ) %>%
  select(Is_Responsive, Is_gbM_Cahn, Is_H2AZ_Depleted) %>%
  # Replace NA with FALSE for plotting
  mutate_all(~ifelse(is.na(.), FALSE, .))

p_upset <- ComplexUpset::upset(
  upset_data,
  colnames(upset_data),
  name = "Gene Sets",
  width_ratio = 0.1
)

ggsave(file.path(output_dir, "upset_responsive_gbm_h2az.png"), plot = p_upset, width = 10, height = 7)

message("Script finished successfully.")