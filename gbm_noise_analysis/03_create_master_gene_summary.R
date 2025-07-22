# ==============================================================================
# SCRIPT: 03_create_master_gene_summary.R
#
# PURPOSE:
#   This script consolidates all previously generated results into a single,
#   comprehensive summary table. It loads the validated noise data and the
#   statistically defined responsive gene list, summarizes the noise for each
#   gene, and merges all annotations into one master file.
#
#   This final table is a key deliverable, perfect for supplementary data in a
#   publication and for downstream exploratory analysis.
#
# MAJOR CORRECTIONS IMPLEMENTED:
#   1.  **Removed Arbitrary Cutoffs:** The flawed "top/bottom 10%" analysis for
#       high/low noise and responsiveness has been completely removed.
#   2.  **Focus on Consolidation:** The script's new, clear purpose is to
#       integrate results, not perform new statistical tests (which are already
#       done correctly in scripts 01 and 02).
#   3.  **Correct Inputs:** The script now loads the validated outputs from the
#       previous two scripts, ensuring the entire pipeline is connected and sound.
#   4.  **Summarized Noise Metric:** It calculates a robust summary of each gene's
#       intrinsic noise (e.g., median within-cell-type CV) to provide a single,
#       useful noise value per gene.
#
# ==============================================================================

# --- 1. SETUP: Load libraries and define paths ---

# Ensure required packages are installed and loaded
packages_to_load <- c("dplyr", "tidyr")
for (pkg in packages_to_load) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, repos = "https://cloud.r-project.org/")
  }
  library(pkg, character.only = TRUE)
}

# --- Define I/O Paths ---
# INPUTS (from the previous two corrected scripts)
processed_noise_data_path <- "/group/sms029/mnieuwenh/gbm_noise_analysis/results/01_noise_analysis/within_celltype_noise_data.csv"
responsive_genes_path <- "/group/sms029/mnieuwenh/gbm_noise_analysis/results/02_responsive_gene_analysis/responsive_genes_list.csv"
full_marker_list_path <- "/group/sms029/mnieuwenh/gbm_noise_analysis/results/02_responsive_gene_analysis/all_celltype_markers_full_list.csv"

# OUTPUTS
output_dir <- "/group/sms029/mnieuwenh/gbm_noise_analysis/results/03_master_summary"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

master_summary_path <- file.path(output_dir, "master_gene_summary_table.csv")


# --- 2. LOAD AND PROCESS DATA ---

message("Loading processed data from previous steps...")

# Load the validated, within-cell-type noise data
noise_data <- read.csv(processed_noise_data_path)

# Load the list of statistically-defined responsive genes
responsive_gene_list <- read.csv(responsive_genes_path)$gene

# Load the full marker list to get details about responsiveness (e.g., which cell type)
marker_details <- read.csv(full_marker_list_path) %>%
  # Filter for only the significant markers to match our responsive definition
  filter(p_val_adj < 0.05) %>%
  # In case a gene is a marker for multiple cell types, we'll keep the one where it has the highest log2FC
  group_by(gene) %>%
  arrange(desc(avg_log2FC)) %>%
  slice(1) %>%
  ungroup() %>%
  select(gene, responsive_cell_type = cluster, responsive_avg_log2FC = avg_log2FC)


# --- 3. SUMMARIZE NOISE PER GENE ---

message("Summarizing noise metrics for each gene...")

# To create a single noise value per gene, we can calculate the median or mean
# of its within-cell-type CVs. The median is generally more robust to outliers.
# We also count how many cell types each gene was expressed in.
gene_noise_summary <- noise_data %>%
  # We should only consider noise values where the gene is reasonably expressed
  filter(mean_expr > 0.1) %>%
  group_by(gene) %>%
  summarise(
    # Median CV is a robust measure of a gene's typical intrinsic noise
    median_within_celltype_cv = median(cv_expr, na.rm = TRUE),
    # Mean CV can also be useful
    mean_within_celltype_cv = mean(cv_expr, na.rm = TRUE),
    # Number of cell types where the gene is expressed and noise could be calculated
    num_celltypes_expressed = n()
  )


# --- 4. CREATE AND ANNOTATE THE MASTER TABLE ---

message("Creating the master gene summary table...")

# Start with a unique list of all genes from the noise analysis
all_genes <- unique(noise_data$gene)

# Create the base data frame
master_table <- data.frame(gene = all_genes)

# --- Merge all annotations ---

# 1. Add the summarized noise metrics
master_table <- master_table %>%
  left_join(gene_noise_summary, by = "gene")

# 2. Add responsive gene status and details
master_table <- master_table %>%
  mutate(is_responsive = gene %in% responsive_gene_list) %>%
  left_join(marker_details, by = "gene")

# 3. Add epigenetic annotations (from the `noise_data` table, which already has them)
epigenetic_annotations <- noise_data %>%
  distinct(gene, cahn_group, bewick_group, h2az_group)

master_table <- master_table %>%
  left_join(epigenetic_annotations, by = "gene")


# --- Reorder and clean the final table ---
master_table <- master_table %>%
  select(
    gene,
    is_responsive,
    responsive_cell_type,
    responsive_avg_log2FC,
    median_within_celltype_cv,
    mean_within_celltype_cv,
    num_celltypes_expressed,
    cahn_group,
    bewick_group,
    h2az_group
  ) %>%
  # Arrange by a meaningful metric, e.g., noise or responsiveness
  arrange(desc(median_within_celltype_cv))

message("Master table created successfully.")


# --- 5. SAVE THE FINAL OUTPUT ---

message(paste("Saving master gene summary table to:", master_summary_path))
write.csv(master_table, master_summary_path, row.names = FALSE)

# Optionally, save a summary text file describing the final output
summary_text <- c(
  "Master Gene Summary Table Generation Complete.",
  "",
  paste("A master table with", nrow(master_table), "genes and", ncol(master_table), "columns has been generated."),
  "This table integrates information on gene identity, responsiveness (statistically defined),",
  "summarized intrinsic noise (median within-cell-type CV), and epigenetic annotations.",
  "",
  paste("The final file is located at:", master_summary_path)
)
writeLines(summary_text, file.path(output_dir, "master_table_generation_summary.txt"))

message("Script finished successfully.")