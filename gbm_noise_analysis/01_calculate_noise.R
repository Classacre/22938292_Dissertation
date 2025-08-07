# ==============================================================================
# SCRIPT: 03_calculate_noise.R
#
# PURPOSE:
#   This script performs the core noise analysis based on the Cortijo et al.
#   (2019) methodology. It moves beyond raw CV to a statistically robust,
#   corrected measure of gene expression noise.
#
# METHODOLOGY:
#   1. Load data exported from script 00.
#   2. Calculate raw noise metrics (mean, CV²) within each cell type.
#   3. Statistically identify Highly Variable Genes (HVGs) using Seurat's
#      Variance-Stabilizing Transformation (VST), which models the mean-variance
#      relationship.
#   4. Generate the final, corrected noise metric ('variance.standardized')
#      for each gene. This metric is decoupled from mean expression.
#   5. Integrate all metrics and apply robust filtering.
#   6. Perform initial analysis and visualization using the CORRECTED noise metric.
# ==============================================================================

# --- 1. SETUP: Load libraries and define paths ---

# Ensure required packages are installed and loaded
packages_to_load <- c("Seurat", "dplyr", "tidyr", "ggplot2", "patchwork", "gridExtra")
for (pkg in packages_to_load) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, repos = "https://cloud.r-project.org/")
  }
  library(pkg, character.only = TRUE)
}

# --- Define I/O Paths ---
# INPUTS

# Data exported from script 02
expr_matrix_path <- "/group/sms029/mnieuwenh/gbm_noise_analysis/02_expression_matrix.csv"
# CORRECTED PATH: Pointing to the correct location in seurat_metadata (update if needed)
cell_meta_path <- "/group/sms029/mnieuwenh/gbm_noise_analysis/02_seurat_metadata_full.csv"
# Other inputs
gene_anno_path <- "/group/sms029/mnieuwenh/Methylation_Data/01_combined_methylation_data.csv"
# We need the full Seurat object for VST modeling.
seurat_object_path <- "/group/sms029/Oliva_dataset/integrated_col_trajectories_colonly.rds"

# OUTPUTS
output_dir <- "/group/sms029/mnieuwenh/gbm_noise_analysis/results/03_noise_analysis"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

processed_noise_data_path <- file.path(output_dir, "noise_analysis_complete_data.csv")
hvg_plot_path <- file.path(output_dir, "hvg_identification_plot.png")
stats_outfile <- file.path(output_dir, "overall_noise_statistics_corrected.txt")


# --- 2. DATA LOADING AND PREPARATION ---

message("Loading and preparing data...")
expr_mat_log <- as.matrix(read.csv(expr_matrix_path, row.names = 1, check.names = FALSE))
expr_mat <- expm1(expr_mat_log) # Convert to linear scale
rm(expr_mat_log)

cell_metadata <- read.csv(cell_meta_path, row.names = 1, check.names = FALSE)

gene_annotations <- read.csv(gene_anno_path, stringsAsFactors = FALSE) %>%
  mutate(
    cahn_group = case_when(
      Cahn_Methylation_status == "gbM" ~ "gbM",
      Cahn_Methylation_status == "TE-like methylation" ~ "TE-like",
      Cahn_Methylation_status == "Unmethylated" ~ "Unmethylated",
      TRUE ~ NA_character_
    ),
    bewick_group = case_when(
      Bewick_Classification == "gbM" ~ "gbM",
      Bewick_Classification %in% c("mCHH", "mCHG") ~ "TE-like",
      Bewick_Classification == "Unmethylated" ~ "Unmethylated",
      TRUE ~ NA_character_
    ),
    h2az_group = case_when(
      H2AZ_Depleted == TRUE ~ "H2A.Z-Depleted",
      H2AZ_Enriched == TRUE ~ "H2A.Z-Enriched",
      TRUE ~ "Other"
    )
  ) %>%
  select(gene = Gene_ID, cahn_group, bewick_group, h2az_group) %>%
  distinct(gene, .keep_all = TRUE)

message("Data loading complete.")


# --- 3. CALCULATE RAW WITHIN-CELL-TYPE NOISE METRICS ---

message("Calculating raw noise metrics within each cell type...")
cell_types <- unique(cell_metadata$identity)
MIN_CELLS_FOR_NOISE_CALC <- 20

noise_data_list <- lapply(cell_types, function(ct) {
  if (is.na(ct) || ct == "") return(NULL)
  cells_in_ct <- rownames(cell_metadata)[cell_metadata$identity == ct]
  
  if (length(cells_in_ct) < MIN_CELLS_FOR_NOISE_CALC) return(NULL)
  
  expr_mat_ct <- expr_mat[, cells_in_ct, drop = FALSE]
  
  mean_expr_ct <- rowMeans(expr_mat_ct)
  var_expr_ct <- apply(expr_mat_ct, 1, var)
  
  data.frame(
    gene = rownames(expr_mat_ct),
    cell_type = ct,
    n_cells = length(cells_in_ct),
    mean_expr = mean_expr_ct,
    variance = var_expr_ct,
    cv2 = var_expr_ct / (mean_expr_ct^2),
    pct_expressing = rowSums(expr_mat_ct > 0) / length(cells_in_ct),
    stringsAsFactors = FALSE
  )
})

raw_noise_data <- bind_rows(noise_data_list)
message("Raw noise calculation complete.")


# --- 4. STATISTICAL HVG IDENTIFICATION & CORRECTED NOISE METRIC ---

message("Loading full Seurat object for VST modeling...")
seurat_obj <- readRDS(seurat_object_path)
DefaultAssay(seurat_obj) <- "RNA"

message("Identifying HVGs and calculating corrected noise scores using VST...")
seurat_obj <- FindVariableFeatures(seurat_obj, selection.method = "vst", nfeatures = 2000)

hvg_info <- HVFInfo(seurat_obj, assay = "RNA")
hvg_info$gene <- rownames(hvg_info)

# CORRECTED: Use 'variance.standardized' as the column name
top_10_pct_cutoff <- quantile(hvg_info$variance.standardized, 0.9, na.rm = TRUE)
hvg_info <- hvg_info %>%
  mutate(is_hvg = variance.standardized >= top_10_pct_cutoff)

# Visualize the HVG selection process
p_hvg <- VariableFeaturePlot(seurat_obj)
p_hvg_labeled <- LabelPoints(plot = p_hvg, points = head(VariableFeatures(seurat_obj), 10), repel = TRUE) +
  ggtitle("Identification of Highly Variable Genes (VST Method)")

ggsave(hvg_plot_path, plot = p_hvg_labeled, width = 10, height = 6)
message(paste("HVG plot saved to:", hvg_plot_path))
rm(seurat_obj); gc() # Free up memory


# --- 5. INTEGRATE DATA AND APPLY ROBUST FILTERING ---

message("Integrating all metrics into a final data table...")
# Merge raw noise data with the corrected VST scores and annotations
# CORRECTED: Select 'variance.standardized'
complete_noise_data <- raw_noise_data %>%
  left_join(hvg_info %>% select(gene, variance.standardized, is_hvg), by = "gene") %>%
  left_join(gene_annotations, by = "gene")

message(paste("Saving final complete (unfiltered) data to:", processed_noise_data_path))
write.csv(complete_noise_data, processed_noise_data_path, row.names = FALSE)

# Now, create the filtered version for the analysis within THIS script.
# CORRECTED: Filter using 'variance.standardized'
noise_data_filtered <- complete_noise_data %>%
  filter(
    mean_expr > 0.01,
    pct_expressing > 0.1,
    is.finite(variance.standardized)
  )


# --- 6. ANALYSIS & VISUALIZATION USING THE CORRECTED NOISE METRIC ---

message("Performing initial analysis using the corrected noise metric (VST standardized variance)...")

# --- Visualization Function ---
plot_corrected_noise <- function(data, group_col, title) {
  plot_data <- data %>% filter(!is.na(.data[[group_col]]))
  
  # === MODIFICATION: Reorder factor levels for specific plots ===
  desired_order <- c("gbM", "Unmethylated", "TE-like")
  current_levels <- unique(as.character(plot_data[[group_col]]))
  
  # Apply reordering only if the column is one of the targeted ones and contains all levels
  if (all(desired_order %in% current_levels) && (group_col == "cahn_group" || group_col == "bewick_group")) {
    plot_data[[group_col]] <- factor(plot_data[[group_col]], levels = desired_order)
    message(paste("Reordered x-axis levels for the", group_col, "plot."))
  }
  # === END MODIFICATION ===

  # CORRECTED: Use 'variance.standardized'
  y_zoom <- quantile(plot_data$variance.standardized, probs = c(0.01, 0.99), na.rm = TRUE)

  # CORRECTED: Use 'variance.standardized'
  p <- ggplot(plot_data, aes(x = .data[[group_col]], y = variance.standardized, fill = .data[[group_col]])) +
    geom_violin(trim = FALSE, scale = "width", alpha = 0.8) +
    geom_boxplot(width = 0.1, outlier.size = 0.5, fill = "white") +
    labs(
      title = title,
      subtitle = "Comparing within-cell-type noise, corrected for mean expression",
      x = "Gene Group",
      y = "Corrected Noise (Standardized Variance)"
    ) +
    coord_cartesian(ylim = y_zoom) +
    theme_bw(base_size = 14) +
    theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust = 1))
  return(p)
}

p_cahn_corrected   <- plot_corrected_noise(noise_data_filtered, "cahn_group", "Corrected Noise Comparison by Cahn Classification")
p_bewick_corrected <- plot_corrected_noise(noise_data_filtered, "bewick_group", "Corrected Noise Comparison by Bewick Classification")
p_h2az_corrected   <- plot_corrected_noise(noise_data_filtered, "h2az_group", "Corrected Noise Comparison by H2A.Z Occupancy")

# Save the primary comparison plots
ggsave(file.path(output_dir, "overall_noise_cahn_corrected.png"),   plot = p_cahn_corrected,   width = 8, height = 7)
ggsave(file.path(output_dir, "overall_noise_bewick_corrected.png"), plot = p_bewick_corrected, width = 8, height = 7)
ggsave(file.path(output_dir, "overall_noise_h2az_corrected.png"),   plot = p_h2az_corrected,   width = 8, height = 7)

# --- Statistical Testing Function ---
run_corrected_stats <- function(data, group_col) {
  stat_data <- data %>% filter(!is.na(.data[[group_col]]))
  if (n_distinct(stat_data[[group_col]]) < 2) return(NULL)
  # CORRECTED: Use 'variance.standardized'
  kruskal_test <- kruskal.test(variance.standardized ~ get(group_col), data = stat_data)
  pairwise_tests <- pairwise.wilcox.test(stat_data$variance.standardized, stat_data[[group_col]], p.adjust.method = "BH")
  return(list(kruskal = kruskal_test, pairwise = pairwise_tests))
}

stats_cahn_corrected <- run_corrected_stats(noise_data_filtered, "cahn_group")
stats_bewick_corrected <- run_corrected_stats(noise_data_filtered, "bewick_group")
stats_h2az_corrected <- run_corrected_stats(noise_data_filtered, "h2az_group")

# --- Save Statistical Results ---
sink(stats_outfile)
cat("=================================================================\n")
cat(" OVERALL NOISE ANALYSIS - STATISTICAL SUMMARY (CORRECTED METRIC)\n")
cat("=================================================================\n\n")
cat("NOTE: All tests performed on the standardized variance score.\n")
cat("Data is filtered for mean_expr > 0.01 and pct_expressing > 0.1 in each cell type.\n\n")

cat("--- Cahn Classification ---\n"); print(stats_cahn_corrected$kruskal); cat("\nPairwise (BH-adjusted):\n"); print(stats_cahn_corrected$pairwise)
cat("\n\n--- Bewick Classification ---\n"); print(stats_bewick_corrected$kruskal); cat("\nPairwise (BH-adjusted):\n"); print(stats_bewick_corrected$pairwise)
cat("\n\n--- H2A.Z Occupancy ---\n"); print(stats_h2az_corrected$kruskal); cat("\nPairwise (BH-adjusted):\n"); print(stats_h2az_corrected$pairwise)
sink()

message(paste("Statistical summary saved to:", stats_outfile))
message("Script finished successfully.")