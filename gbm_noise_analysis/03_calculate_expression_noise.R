# ==============================================================================
# SCRIPT: 03_calculate_expression_noise.R
#
# PURPOSE:
#   This script performs the central analysis: quantifying gene expression noise.
#   It uses Seurat's VST to calculate a corrected noise score and includes
#   the DM (Distance to Median) metric for validation. It now also performs
#   statistical tests on noise drivers within each individual cell type.
#
# METHODOLOGY:
#   1. Load expression data and the enhanced epigenetic/promoter annotations.
#   2. Calculate raw noise metrics (mean, CV²) and DM for each gene within each cell identity.
#   3. Use VST to extract the corrected noise metric ('variance.standardized').
#   4. Integrate all metrics and annotations into a single data frame.
#   5. Perform statistical tests on the full dataset and on cell-type subsets.
# ==============================================================================

# --- 1. SETUP: Load libraries and define paths ---
packages_to_load <- c("Seurat", "dplyr", "tidyr", "ggplot2", "patchwork", "SingleCellExperiment")
for (pkg in packages_to_load) {
  if (!requireNamespace(pkg, quietly = TRUE)) install.packages(pkg, repos = "https://cloud.r-project.org/")
  library(pkg, character.only = TRUE)
}
# Add scran for alternative noise metric
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager", repos = "https://cloud.r-project.org/")
if (!requireNamespace("scran", quietly = TRUE)) BiocManager::install("scran")
library(scran)

# --- Define I/O Paths ---
expr_matrix_path <- "/group/sms029/mnieuwenh/gbm_noise_analysis/results/02_prepare_expression_data/02_expression_matrix.csv"
cell_meta_path <- "/group/sms029/mnieuwenh/gbm_noise_analysis/results/02_prepare_expression_data/02_cell_metadata_annotated.csv"
gene_anno_path <- "/group/sms029/mnieuwenh/gbm_noise_analysis/results/01_prepare_epigenetic_data/01_epigenetic_annotations.csv"
seurat_object_path <- "/group/sms029/Oliva_dataset/integrated_col_trajectories_colonly.rds"
output_dir <- "/group/sms029/mnieuwenh/gbm_noise_analysis/results/03_calculate_expression_noise"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# Create the new subdirectory for cell type specific results
celltype_output_dir <- file.path(output_dir, "celltype_analyses")
dir.create(celltype_output_dir, showWarnings = FALSE, recursive = TRUE)

processed_noise_data_path <- file.path(output_dir, "03_gene_noise_metrics_complete.csv")
hvg_plot_path <- file.path(output_dir, "03_hvg_identification_plot.png")
stats_outfile <- file.path(output_dir, "03_noise_statistics_summary.txt")
celltype_stats_outfile <- file.path(celltype_output_dir, "03_celltype_noise_statistics_summary.txt")


# --- 2. DATA LOADING AND PREPARATION ---
message("Loading and preparing data...")
expr_mat_log <- as.matrix(read.csv(expr_matrix_path, row.names = 1, check.names = FALSE))
expr_mat <- expm1(expr_mat_log)
rm(expr_mat_log); gc()
cell_metadata <- read.csv(cell_meta_path, row.names = 1, check.names = FALSE)
gene_annotations <- read.csv(gene_anno_path, stringsAsFactors = FALSE) %>%
  mutate(
    cahn_group = case_when(
      Cahn_Methylation_status == "gbM" ~ "gbM",
      Cahn_Methylation_status == "TE-like methylation" ~ "TE-like",
      Cahn_Methylation_status == "Unmethylated" ~ "Unmethylated",
      TRUE ~ NA_character_
    ),
    h2az_group = case_when(
      H2AZ_Depleted == TRUE ~ "H2A.Z-Depleted",
      H2AZ_Enriched == TRUE ~ "H2A.Z-Enriched",
      TRUE ~ "Other"
    )
  ) %>%
  select(gene = Gene_ID, cahn_group, h2az_group, has_TATA_box) %>%
  distinct(gene, .keep_all = TRUE)

# --- 3. CALCULATE RAW & ALTERNATIVE WITHIN-IDENTITY NOISE METRICS ---
message("Calculating raw (mean, CV^2) and alternative (DM) noise metrics within each cell identity...")
cell_identities <- unique(cell_metadata$identity)
MIN_CELLS_FOR_NOISE_CALC <- 20
noise_data_list <- lapply(cell_identities, function(id) {
  if (is.na(id) || id == "") return(NULL)
  cells_in_id <- rownames(cell_metadata)[cell_metadata$identity == id]
  if (length(cells_in_id) < MIN_CELLS_FOR_NOISE_CALC) return(NULL)
  expr_mat_id <- expr_mat[, cells_in_id, drop = FALSE]
  mean_expr_id <- rowMeans(expr_mat_id)
  var_expr_id <- apply(expr_mat_id, 1, var)
  sce_id <- SingleCellExperiment(assays = list(counts = expr_mat_id))
  dm_id <- tryCatch({DM(sce_id)}, error = function(e) {rep(NA_real_, nrow(sce_id))})
  data.frame(
    gene = rownames(expr_mat_id), identity = id, n_cells = length(cells_in_id),
    mean_expr = mean_expr_id, variance = var_expr_id, cv2 = var_expr_id / (mean_expr_id^2),
    dm = dm_id, pct_expressing = rowSums(expr_mat_id > 0) / length(cells_in_id),
    stringsAsFactors = FALSE
  )
})
raw_noise_data <- bind_rows(noise_data_list)

# --- 4. STATISTICAL HVG IDENTIFICATION & CORRECTED NOISE METRIC ---
message("Loading full Seurat object for VST modeling...")
seurat_obj <- readRDS(seurat_object_path)
DefaultAssay(seurat_obj) <- "RNA"
message("Calculating corrected noise scores using VST...")
seurat_obj <- FindVariableFeatures(seurat_obj, selection.method = "vst", nfeatures = 2000)
hvg_info <- HVFInfo(seurat_obj, assay = "RNA")
hvg_info$gene <- rownames(hvg_info)
p_hvg <- VariableFeaturePlot(seurat_obj)
p_hvg_labeled <- LabelPoints(plot = p_hvg, points = head(VariableFeatures(seurat_obj), 10), repel = TRUE) +
  ggtitle("Identification of Highly Variable Genes (VST Method)")
ggsave(hvg_plot_path, plot = p_hvg_labeled, width = 10, height = 6)
n_variable_features <- length(VariableFeatures(seurat_obj))
rm(seurat_obj); gc()

# --- 5. INTEGRATE DATA AND SAVE ---
message("Integrating all metrics into a final data table...")

identity_metadata_map <- cell_metadata %>%
  select(identity, base_cell_type, ers_state) %>%
  distinct(identity, .keep_all = TRUE)

complete_noise_data <- raw_noise_data %>%
  left_join(hvg_info %>% select(gene, variance.standardized), by = "gene") %>%
  left_join(gene_annotations, by = "gene") %>%
  left_join(identity_metadata_map, by = "identity")

message(paste("Saving final complete data to:", processed_noise_data_path))
write.csv(complete_noise_data, processed_noise_data_path, row.names = FALSE)

# --- 6. PERFORM OVERALL STATISTICAL ANALYSIS ---
message("Performing overall statistical analysis...")
noise_data_filtered <- complete_noise_data %>%
  filter(mean_expr > 0.01, pct_expressing > 0.1, is.finite(variance.standardized))
run_stats <- function(data, group_col) {
  stat_data <- data %>% filter(!is.na(.data[[group_col]]))
  if (n_distinct(stat_data[[group_col]]) < 2) return(NULL)
  kruskal_test <- kruskal.test(variance.standardized ~ get(group_col), data = stat_data)
  pairwise_tests <- pairwise.wilcox.test(stat_data$variance.standardized, stat_data[[group_col]], p.adjust.method = "BH")
  return(list(kruskal = kruskal_test, pairwise = pairwise_tests))
}
stats_cahn <- run_stats(noise_data_filtered, "cahn_group")
stats_h2az <- run_stats(noise_data_filtered, "h2az_group")
stats_tata <- run_stats(noise_data_filtered, "has_TATA_box")

sink(stats_outfile)
cat("=================================================================\n")
cat("       STATISTICAL SUMMARY: 03_calculate_expression_noise.R\n")
cat("=================================================================\n\n")
cat(paste("Summary generated on:", Sys.Date()), "\n\n")
cat("--- Data Processing and Filtering Summary ---\n")
cat("NOTE: The dataset now includes an alternative noise metric 'dm' for validation.\n")
cat(paste("Gene-identity pairs retained after filtering for min cells (>", MIN_CELLS_FOR_NOISE_CALC, "):", format(nrow(raw_noise_data), big.mark=",")), "\n")
cat(paste("Total gene-identity pairs in the final complete dataset:", format(nrow(complete_noise_data), big.mark=",")), "\n")
cat(paste("Number of Highly Variable Features identified by VST:", n_variable_features), "\n\n")
cat("--- Overall Noise Analysis (Corrected Metric on Filtered Data) ---\n")
cat(paste("Number of pairs used for statistical tests (mean_expr > 0.01 & pct > 0.1):", format(nrow(noise_data_filtered), big.mark=",")), "\n\n")
cat("--- Cahn Classification ---\n"); print(stats_cahn$kruskal); cat("\nPairwise (BH-adjusted p-values):\n"); print(stats_cahn$pairwise); cat("\n\n")
cat("--- H2A.Z Occupancy ---\n"); print(stats_h2az$kruskal); cat("\nPairwise (BH-adjusted p-values):\n"); print(stats_h2az$pairwise); cat("\n\n")
cat("--- TATA Box Presence ---\n"); print(stats_tata$kruskal); cat("\nPairwise (BH-adjusted p-values):\n"); print(stats_tata$pairwise); cat("\n\n")
sink()

message(paste("Statistical summary saved to:", stats_outfile))

# ==============================================================================
# --- 7. CELL-TYPE-SPECIFIC STATISTICAL ANALYSIS ---
# ==============================================================================
message("Performing cell-type-specific statistical analysis...")

# Get the list of unique base cell types to loop through
unique_cell_types <- unique(na.omit(noise_data_filtered$base_cell_type))

# Open the output file once
sink(celltype_stats_outfile)
cat("=================================================================\n")
cat("   CELL-TYPE-SPECIFIC STATISTICAL SUMMARY: 03_calculate_expression_noise.R\n")
cat("=================================================================\n\n")
cat(paste("Summary generated on:", Sys.Date()), "\n\n")
cat("This file contains statistical tests on the drivers of expression noise,\n")
cat("run independently for each base cell type.\n\n")

# Loop through each cell type
for(cell_type in unique_cell_types) {
  
  message(paste("... Analyzing cell type:", cell_type))
  
  # Filter data for the current cell type
  cell_type_data <- noise_data_filtered %>%
    filter(base_cell_type == cell_type)
    
  cat("-----------------------------------------------------------------\n")
  cat(paste("                ANALYSIS FOR CELL TYPE:", toupper(cell_type)), "\n")
  cat("-----------------------------------------------------------------\n")
  cat(paste("Number of gene-identity pairs used for tests:", format(nrow(cell_type_data), big.mark=",")), "\n\n")

  # Run stats for Cahn classification
  stats_cahn_ct <- run_stats(cell_type_data, "cahn_group")
  cat("--- Cahn Classification ---\n")
  if (!is.null(stats_cahn_ct)) {
    print(stats_cahn_ct$kruskal)
    cat("\nPairwise (BH-adjusted p-values):\n"); print(stats_cahn_ct$pairwise)
  } else {
    cat("Not enough data to perform test.\n")
  }
  cat("\n\n")

  # Run stats for H2A.Z occupancy
  stats_h2az_ct <- run_stats(cell_type_data, "h2az_group")
  cat("--- H2A.Z Occupancy ---\n")
  if (!is.null(stats_h2az_ct)) {
    print(stats_h2az_ct$kruskal)
    cat("\nPairwise (BH-adjusted p-values):\n"); print(stats_h2az_ct$pairwise)
  } else {
    cat("Not enough data to perform test.\n")
  }
  cat("\n\n")

  # Run stats for TATA box presence
  stats_tata_ct <- run_stats(cell_type_data, "has_TATA_box")
  cat("--- TATA Box Presence ---\n")
  if (!is.null(stats_tata_ct)) {
    print(stats_tata_ct$kruskal)
    cat("\nPairwise (BH-adjusted p-values):\n"); print(stats_tata_ct$pairwise)
  } else {
    cat("Not enough data to perform test.\n")
  }
  cat("\n\n")
}

# Close the file connection
sink()
message(paste("Cell-type-specific statistical summary saved to:", celltype_stats_outfile))
message("Script 03 finished successfully.")