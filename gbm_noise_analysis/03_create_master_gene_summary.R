# ==============================================================================
# SCRIPT: 03_create_master_gene_summary.R
#
# PURPOSE:
#   This script consolidates all previously generated results into a single,
#   comprehensive summary table. It integrates the corrected noise metrics,
#   HVG status, responsive gene status, and epigenetic annotations for every
#   gene in the analysis.
#
#   This final table is a key deliverable, perfect for supplementary data and
#   for downstream exploratory analysis.
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
# INPUTS (from the previous two rewritten scripts)
complete_noise_data_path <- "/group/sms029/mnieuwenh/gbm_noise_analysis/results/01_noise_analysis/noise_analysis_complete_data.csv"
full_marker_list_path <- "/group/sms029/mnieuwenh/gbm_noise_analysis/results/02_responsive_gene_analysis/all_celltype_markers_full_list.csv"

# OUTPUTS
output_dir <- "/group/sms029/mnieuwenh/gbm_noise_analysis/results/03_master_summary"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

master_summary_path <- file.path(output_dir, "master_gene_summary_table.csv")


# --- 2. LOAD AND PROCESS DATA ---

message("Loading processed data from previous steps...")

# Load the comprehensive noise data from script 01
complete_noise_data <- read.csv(complete_noise_data_path)

# Load the full marker list from script 02 to get details about responsiveness
responsive_gene_details <- read.csv(full_marker_list_path) %>%
  filter(p_val_adj < 0.05) %>%
  group_by(gene) %>%
  arrange(desc(avg_log2FC)) %>%
  slice(1) %>%
  ungroup() %>%
  select(gene, responsive_cell_type = cluster, responsive_avg_log2FC = avg_log2FC)

responsive_gene_list <- unique(responsive_gene_details$gene)


# --- 3. CREATE THE MASTER GENE TABLE ---

message("Creating the master gene summary table...")

# CORRECTED: Select 'variance.standardized'
gene_summary_base <- complete_noise_data %>%
  distinct(gene, .keep_all = TRUE) %>%
  select(
    gene,
    variance.standardized,
    is_hvg,
    cahn_group,
    bewick_group,
    h2az_group
  )

raw_noise_summary <- complete_noise_data %>%
  filter(mean_expr > 0.01, pct_expressing > 0.1) %>%
  group_by(gene) %>%
  summarise(
    median_within_celltype_cv2 = median(cv2, na.rm = TRUE),
    num_celltypes_expressed = n()
  )

master_table <- gene_summary_base %>%
  left_join(raw_noise_summary, by = "gene") %>%
  left_join(responsive_gene_details, by = "gene") %>%
  mutate(is_responsive = gene %in% responsive_gene_list)


# --- 4. FINALIZE AND REORDER THE TABLE FOR CLARITY ---

# CORRECTED: Reorder and arrange using 'variance.standardized'
master_table <- master_table %>%
  select(
    gene,
    # Core Noise Metrics
    variance.standardized, # This is the primary, corrected noise metric
    is_hvg,
    # Responsiveness Info
    is_responsive,
    responsive_cell_type,
    responsive_avg_log2FC,
    # Raw Noise Summary (for reference)
    median_within_celltype_cv2,
    num_celltypes_expressed,
    # Epigenetic Annotations
    cahn_group,
    bewick_group,
    h2az_group
  ) %>%
  arrange(desc(variance.standardized))

message("Master table created successfully.")


# --- 5. SAVE THE FINAL OUTPUT ---

message(paste("Saving master gene summary table to:", master_summary_path))
write.csv(master_table, master_summary_path, row.names = FALSE)

summary_text <- c(
  "Master Gene Summary Table Generation Complete.",
  "",
  paste("A master table with", nrow(master_table), "genes and", ncol(master_table), "columns has been generated."),
  "",
  "This table integrates information on:",
  "  - Gene identity (gene)",
  "  - Corrected noise (variance.standardized) <- PRIMARY NOISE METRIC",
  "  - High-Variability status (is_hvg)",
  "  - Responsiveness status and details (is_responsive, responsive_cell_type, etc.)",
  "  - Summarized raw noise for reference (median_within_celltype_cv2)",
  "  - Epigenetic annotations (cahn_group, bewick_group, h2az_group)",
  "",
  paste("The final file is located at:", master_summary_path)
)
writeLines(summary_text, file.path(output_dir, "master_table_generation_summary.txt"))

message("Script finished successfully.")