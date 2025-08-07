# ==============================================================================
# SCRIPT: 03_calculate_expression_noise.R
#
# PURPOSE:
#   This script performs the central analysis of the study: quantifying gene
#   expression noise. To control for the known, confounding relationship
#   between mean expression and variance, we employ Seurat's Variance-
#   Stabilizing Transformation (VST). This method models the mean-variance
#   relationship and produces a standardized, corrected noise score for each
#   gene, which serves as our primary metric for all subsequent comparisons.
#
# METHODOLOGY:
#   1. Load expression data (from script 02) and epigenetic annotations
#      (from script 01).
#   2. Calculate raw noise metrics (mean, CV²) for each gene within each cell type.
#   3. Use the full Seurat object to model the mean-variance trend via VST.
#   4. Extract the corrected noise metric ('variance.standardized') for each gene.
#   5. Integrate raw metrics, corrected metrics, and epigenetic annotations.
#   6. Perform initial statistical comparisons and visualizations of the
#      CORRECTED noise metric across different epigenetic groups.
# ==============================================================================

# --- 1. SETUP: Load libraries and define paths ---

# Ensure required packages are installed and loaded for analysis and plotting
packages_to_load <- c("Seurat", "dplyr", "tidyr", "ggplot2", "patchwork")
for (pkg in packages_to_load) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, repos = "https://cloud.r-project.org/")
  }
  library(pkg, character.only = TRUE)
}

# --- Define I/O Paths ---
# INPUTS
# Data exported from script 02
expr_matrix_path <- "/group/sms029/mnieuwenh/gbm_noise_analysis/results/02_prepare_expression_data/02_expression_matrix.csv"
cell_meta_path <- "/group/sms029/mnieuwenh/gbm_noise_analysis/results/02_prepare_expression_data/02_cell_metadata.csv"
# Epigenetic annotations from script 01 (UPDATED PATH)
gene_anno_path <- "/group/sms029/mnieuwenh/gbm_noise_analysis/results/01_prepare_epigenetic_data/01_epigenetic_annotations.csv"
# The full Seurat object is required for VST modeling
seurat_object_path <- "/group/sms029/Oliva_dataset/integrated_col_trajectories_colonly.rds"

# OUTPUTS (UPDATED)
output_dir <- "/group/sms029/mnieuwenh/gbm_noise_analysis/results/03_calculate_expression_noise"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# Consistent output filenames
processed_noise_data_path <- file.path(output_dir, "03_gene_noise_metrics_complete.csv")
hvg_plot_path <- file.path(output_dir, "03_hvg_identification_plot.png")
stats_outfile <- file.path(output_dir, "03_noise_statistics_summary.txt")


# --- 2. DATA LOADING AND PREPARATION ---

message("Loading and preparing data...")
# Load log-normalized expression and convert back to linear scale for noise calculation
expr_mat_log <- as.matrix(read.csv(expr_matrix_path, row.names = 1, check.names = FALSE))
expr_mat <- expm1(expr_mat_log)
rm(expr_mat_log); gc() # Free up memory

cell_metadata <- read.csv(cell_meta_path, row.names = 1, check.names = FALSE)

# Load epigenetic annotations and create simplified grouping columns for analysis
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
  distinct(gene, .keep_all = TRUE) # Ensure one entry per gene

message("Data loading complete.")


# --- 3. CALCULATE RAW WITHIN-CELL-TYPE NOISE METRICS ---

message("Calculating raw noise metrics (mean, CV^2) within each cell type...")
cell_types <- unique(cell_metadata$identity)
MIN_CELLS_FOR_NOISE_CALC <- 20 # A sensible threshold to avoid spurious calculations

noise_data_list <- lapply(cell_types, function(ct) {
  if (is.na(ct) || ct == "") return(NULL) # Skip missing or empty cell type labels
  cells_in_ct <- rownames(cell_metadata)[cell_metadata$identity == ct]

  # Skip cell types with too few cells
  if (length(cells_in_ct) < MIN_CELLS_FOR_NOISE_CALC) return(NULL)

  expr_mat_ct <- expr_mat[, cells_in_ct, drop = FALSE]

  # Calculate metrics
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
# This function models the mean-variance relationship and calculates 'variance.standardized'
seurat_obj <- FindVariableFeatures(seurat_obj, selection.method = "vst", nfeatures = 2000)

# Extract the VST information, which includes our corrected noise metric
hvg_info <- HVFInfo(seurat_obj, assay = "RNA")
hvg_info$gene <- rownames(hvg_info)

# Define Highly Variable Genes (HVGs) as those in the top 10% of corrected noise scores
top_10_pct_cutoff <- quantile(hvg_info$variance.standardized, 0.9, na.rm = TRUE)
hvg_info <- hvg_info %>%
  mutate(is_hvg = variance.standardized >= top_10_pct_cutoff)

# Visualize the VST model and highlight the top variable genes
p_hvg <- VariableFeaturePlot(seurat_obj)
p_hvg_labeled <- LabelPoints(plot = p_hvg, points = head(VariableFeatures(seurat_obj), 10), repel = TRUE) +
  ggtitle("Identification of Highly Variable Genes (VST Method)",
          subtitle = "The y-axis represents the corrected noise score (standardized variance)")

ggsave(hvg_plot_path, plot = p_hvg_labeled, width = 10, height = 6)
message(paste("HVG plot saved to:", hvg_plot_path))
rm(seurat_obj); gc() # Free up memory


# --- 5. INTEGRATE DATA AND APPLY ROBUST FILTERING ---

message("Integrating all metrics into a final data table...")
# Merge raw noise data with the corrected VST scores and epigenetic annotations
complete_noise_data <- raw_noise_data %>%
  left_join(hvg_info %>% select(gene, variance.standardized, is_hvg), by = "gene") %>%
  left_join(gene_annotations, by = "gene")

message(paste("Saving final complete (unfiltered) data to:", processed_noise_data_path))
write.csv(complete_noise_data, processed_noise_data_path, row.names = FALSE)

# Create a filtered version for the analysis within THIS script to ensure robustness
noise_data_filtered <- complete_noise_data %>%
  filter(
    mean_expr > 0.01,       # Filter out genes with very low expression
    pct_expressing > 0.1,   # Filter out genes detected in too few cells
    is.finite(variance.standardized) # Ensure corrected noise score is valid
  )


# --- 6. ANALYSIS & VISUALIZATION USING THE CORRECTED NOISE METRIC ---

message("Performing initial analysis using the corrected noise metric (VST standardized variance)...")

# --- Reusable Visualization Function ---
plot_corrected_noise <- function(data, group_col, title) {
  plot_data <- data %>% filter(!is.na(.data[[group_col]]))

  # Reorder factor levels for methylation plots for a more logical presentation
  desired_order <- c("gbM", "Unmethylated", "TE-like")
  if (group_col %in% c("cahn_group", "bewick_group")) {
    current_levels <- unique(as.character(plot_data[[group_col]]))
    if (all(desired_order %in% current_levels)) {
      plot_data[[group_col]] <- factor(plot_data[[group_col]], levels = desired_order)
    }
  }

  # Zoom in on the core 98% of the data to avoid extreme outliers skewing the view
  y_zoom <- quantile(plot_data$variance.standardized, probs = c(0.01, 0.99), na.rm = TRUE)

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

p_cahn_corrected   <- plot_corrected_noise(noise_data_filtered, "cahn_group", "Corrected Noise by Cahn Classification")
p_bewick_corrected <- plot_corrected_noise(noise_data_filtered, "bewick_group", "Corrected Noise by Bewick Classification")
p_h2az_corrected   <- plot_corrected_noise(noise_data_filtered, "h2az_group", "Corrected Noise by H2A.Z Occupancy")

# Save the primary comparison plots with new, consistent filenames
ggsave(file.path(output_dir, "03_noise_by_cahn_group.png"),   plot = p_cahn_corrected,   width = 8, height = 7)
ggsave(file.path(output_dir, "03_noise_by_bewick_group.png"), plot = p_bewick_corrected, width = 8, height = 7)
ggsave(file.path(output_dir, "03_noise_by_h2az_group.png"),   plot = p_h2az_corrected,   width = 8, height = 7)

# --- Reusable Statistical Testing Function ---
run_corrected_stats <- function(data, group_col) {
  stat_data <- data %>% filter(!is.na(.data[[group_col]]))
  if (n_distinct(stat_data[[group_col]]) < 2) return(NULL)
  # Kruskal-Wallis for overall difference, followed by pairwise Wilcoxon for specific comparisons
  kruskal_test <- kruskal.test(variance.standardized ~ get(group_col), data = stat_data)
  pairwise_tests <- pairwise.wilcox.test(stat_data$variance.standardized, stat_data[[group_col]], p.adjust.method = "BH")
  return(list(kruskal = kruskal_test, pairwise = pairwise_tests))
}

stats_cahn_corrected <- run_corrected_stats(noise_data_filtered, "cahn_group")
stats_bewick_corrected <- run_corrected_stats(noise_data_filtered, "bewick_group")
stats_h2az_corrected <- run_corrected_stats(noise_data_filtered, "h2az_group")

# --- Save Statistical Results to a Text File ---
sink(stats_outfile)
cat("=================================================================\n")
cat(" STATISTICAL SUMMARY OF CORRECTED NOISE METRICS\n")
cat("=================================================================\n\n")
cat("NOTE: All tests performed on the VST-based 'variance.standardized' score.\n")
cat("Data is filtered for mean_expr > 0.01 and pct_expressing > 0.1 in each cell type.\n\n")

cat("--- Cahn Classification ---\n"); print(stats_cahn_corrected$kruskal); cat("\nPairwise (BH-adjusted):\n"); print(stats_cahn_corrected$pairwise)
cat("\n\n--- Bewick Classification ---\n"); print(stats_bewick_corrected$kruskal); cat("\nPairwise (BH-adjusted):\n"); print(stats_bewick_corrected$pairwise)
cat("\n\n--- H2A.Z Occupancy ---\n"); print(stats_h2az_corrected$kruskal); cat("\nPairwise (BH-adjusted):\n"); print(stats_h2az_corrected$pairwise)
sink()

message(paste("Statistical summary saved to:", stats_outfile))
message("Script finished successfully.")