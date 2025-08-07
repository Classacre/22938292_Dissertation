# ==============================================================================
# SCRIPT: 07_perform_go_enrichment.R
#
# PURPOSE:
#   This script serves as a supplementary analysis to provide biological context
#   to the findings from the core pipeline. Having identified distinct gene
#   sets based on expression noise and responsiveness, we now investigate their
#   functional roles using Gene Ontology (GO) enrichment analysis.
#
# METHODOLOGY:
#   1. Load the final annotated data from script 05.
#   2. Define distinct gene sets of interest:
#      - Responsive vs. Non-Responsive genes.
#      - High-Noise vs. Low-Noise genes (top and bottom quartiles).
#   3. Use the `clusterProfiler` package to perform GO enrichment analysis for
#      "Biological Process" (BP) terms for each gene set.
#   4. The background for the enrichment test will be all genes that passed
#      our expression filters, ensuring a fair comparison.
#   5. Visualize the most significantly enriched GO terms using dot plots.
#   6. Save the full enrichment results as CSV files for further exploration.
# ==============================================================================

# --- 1. SETUP: Load libraries and define paths ---

# Load packages for analysis and plotting.
# Note: clusterProfiler and its annotation database are from Bioconductor.
packages_to_load <- c("dplyr", "ggplot2", "clusterProfiler", "org.Hs.eg.db")
for (pkg in packages_to_load) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    if (pkg %in% c("clusterProfiler", "org.Hs.eg.db")) {
      # Install Bioconductor packages
      if (!requireNamespace("BiocManager", quietly = TRUE)) {
        install.packages("BiocManager")
      }
      BiocManager::install(pkg)
    } else {
      # Install CRAN packages
      install.packages(pkg, repos = "https://cloud.r-project.org/")
    }
  }
  library(pkg, character.only = TRUE)
}

# --- Define I/O Paths ---
# INPUT: The final annotated data file from script 05 (UPDATED PATH)
annotated_data_path <- "/group/sms029/mnieuwenh/gbm_noise_analysis/results/05_analyze_responsive_genes/05_gene_noise_with_responsiveness.csv"

# OUTPUTS: (UPDATED)
output_dir <- "/group/sms029/mnieuwenh/gbm_noise_analysis/results/07_perform_go_enrichment"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# Define output file paths
go_responsive_plot_path <- file.path(output_dir, "07_go_enrichment_responsive.png")
go_noise_plot_path <- file.path(output_dir, "07_go_enrichment_noise.png")
go_responsive_table_path <- file.path(output_dir, "07_go_results_responsive.csv")
go_nonresponsive_table_path <- file.path(output_dir, "07_go_results_nonresponsive.csv")
go_high_noise_table_path <- file.path(output_dir, "07_go_results_high_noise.csv")
go_low_noise_table_path <- file.path(output_dir, "07_go_results_low_noise.csv")


# --- 2. DATA LOADING AND GENE SET DEFINITION ---

message("Loading final annotated data...")
if (!file.exists(annotated_data_path)) {
  stop(paste("Error: Final data file not found at", annotated_data_path))
}
full_data <- read.csv(annotated_data_path)

message("Applying filters and defining gene sets for analysis...")
# Apply the standard robust filtering
filtered_data <- full_data %>%
  filter(
    mean_expr > 0.01,
    pct_expressing > 0.1,
    is.finite(variance.standardized)
  )

# Define the "universe" of genes for the background set.
# This should be all genes that passed our filters. We use ENTREZ IDs as required by clusterProfiler.
gene_universe_entrez <- bitr(
  unique(filtered_data$gene),
  fromType = "SYMBOL",
  toType = "ENTREZID",
  OrgDb = org.Hs.eg.db
)$ENTREZID

# --- Define Gene Sets of Interest ---
# We create one entry per gene, summarizing its properties across cell types.
gene_summary <- filtered_data %>%
  group_by(gene) %>%
  summarise(
    is_responsive = any(is_responsive),
    mean_noise = mean(variance.standardized, na.rm = TRUE)
  ) %>%
  ungroup()

# Convert gene symbols to ENTREZ IDs for each set
convert_symbols_to_entrez <- function(symbols) {
  bitr(symbols, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)$ENTREZID
}

# Set 1: Responsive Genes
responsive_genes_entrez <- convert_symbols_to_entrez(
  (gene_summary %>% filter(is_responsive == TRUE))$gene
)
# Set 2: Non-Responsive Genes
nonresponsive_genes_entrez <- convert_symbols_to_entrez(
  (gene_summary %>% filter(is_responsive == FALSE))$gene
)

# Set 3 & 4: High- and Low-Noise Genes
noise_quartiles <- quantile(gene_summary$mean_noise, probs = c(0.25, 0.75), na.rm = TRUE)
low_noise_genes_entrez <- convert_symbols_to_entrez(
  (gene_summary %>% filter(mean_noise <= noise_quartiles[1]))$gene
)
high_noise_genes_entrez <- convert_symbols_to_entrez(
  (gene_summary %>% filter(mean_noise >= noise_quartiles[2]))$gene
)

message(paste("Defined Gene Sets: Universe =", length(gene_universe_entrez), "genes,",
              "Responsive =", length(responsive_genes_entrez), ",",
              "High-Noise =", length(high_noise_genes_entrez)))


# --- 3. RUN GO ENRICHMENT ANALYSIS ---

message("Running GO enrichment for 'Biological Process'...")

# Reusable function for running the analysis
run_go_enrichment <- function(gene_set, universe) {
  enrichGO(
    gene = gene_set,
    universe = universe,
    OrgDb = org.Hs.eg.db,
    keyType = "ENTREZID",
    ont = "BP", # Biological Process
    pAdjustMethod = "BH",
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.1,
    readable = TRUE # Converts ENTREZ IDs back to symbols in the output
  )
}

# Run for each gene set
go_responsive <- run_go_enrichment(responsive_genes_entrez, gene_universe_entrez)
go_nonresponsive <- run_go_enrichment(nonresponsive_genes_entrez, gene_universe_entrez)
go_high_noise <- run_go_enrichment(high_noise_genes_entrez, gene_universe_entrez)
go_low_noise <- run_go_enrichment(low_noise_genes_entrez, gene_universe_entrez)

# Save the full results tables
write.csv(as.data.frame(go_responsive), go_responsive_table_path, row.names = FALSE)
write.csv(as.data.frame(go_nonresponsive), go_nonresponsive_table_path, row.names = FALSE)
write.csv(as.data.frame(go_high_noise), go_high_noise_table_path, row.names = FALSE)
write.csv(as.data.frame(go_low_noise), go_low_noise_table_path, row.names = FALSE)

message("GO analysis complete. Full results saved to CSV files.")


# --- 4. VISUALIZE ENRICHMENT RESULTS ---

message("Generating dot plots for visualization...")

# Reusable function to create dot plots
create_go_dotplot <- function(go_results, title, n_terms = 15) {
  if (is.null(go_results) || nrow(go_results) == 0) {
    return(ggplot() + labs(title = title, subtitle = "No significant terms found") + theme_void())
  }
  # Simplify the results for plotting
  go_simplified <- simplify(go_results, cutoff = 0.7, by = "p.adjust", select_fun = min)
  dotplot(go_simplified, showCategory = n_terms, font.size = 10) +
    labs(title = title) +
    theme(plot.title = element_text(face = "bold"))
}

# Generate plots for responsive vs. non-responsive
p_responsive <- create_go_dotplot(go_responsive, "Functional Profile of Responsive Genes")
p_nonresponsive <- create_go_dotplot(go_nonresponsive, "Functional Profile of Non-Responsive Genes")

# Combine and save the first figure
responsive_figure <- p_responsive / p_nonresponsive +
  plot_annotation(tag_levels = 'A')

ggsave(
  go_responsive_plot_path,
  plot = responsive_figure,
  width = 10, height = 12, dpi = 300, bg = "white"
)

# Generate plots for high-noise vs. low-noise
p_high_noise <- create_go_dotplot(go_high_noise, "Functional Profile of High-Noise Genes")
p_low_noise <- create_go_dotplot(go_low_noise, "Functional Profile of Low-Noise Genes")

# Combine and save the second figure
noise_figure <- p_high_noise / p_low_noise +
  plot_annotation(tag_levels = 'A')

ggsave(
  go_noise_plot_path,
  plot = noise_figure,
  width = 10, height = 12, dpi = 300, bg = "white"
)

message(paste("GO enrichment plots saved to:", output_dir))
message("Script finished successfully.")