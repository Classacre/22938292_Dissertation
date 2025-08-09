# ==============================================================================
# SCRIPT: 07_GO_enrichment_analysis.R
#
# PURPOSE:
#   Perform GSEA using clusterProfiler::gseGO on genes ranked by
#   VST-corrected noise, globally and per cell type.
#
# UPDATES/FIXES:
#   - Fixes quoting typo in dotplot(split = ".sign").
#   - Uses centered (z-scored) statistics so ranks contain both positive and
#     negative values; this reduces fgsea warnings and aligns with good GSEA practice.
#   - Sets eps = 0 to better estimate very small p-values.
#   - Defensive tryCatch; qualifies dplyr verbs to avoid masking issues.
#
# REFERENCES:
#   clusterProfiler/gseGO from YuLab@SMU (https://yulab-smu.top/contribution-knowledge-mining/)
# ==============================================================================

# --- 1. SETUP: Load libraries and define paths ---
pkgs <- c("clusterProfiler", "org.At.tair.db", "dplyr", "ggplot2", "patchwork", "enrichplot")
for (pkg in pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    if (pkg %in% c("clusterProfiler", "org.At.tair.db", "enrichplot")) {
      if (!requireNamespace("BiocManager", quietly = TRUE)) {
        install.packages("BiocManager", repos = "https://cloud.r-project.org/")
      }
      BiocManager::install(pkg, ask = FALSE, update = FALSE)
    } else {
      install.packages(pkg, repos = "https://cloud.r-project.org/")
    }
  }
  library(pkg, character.only = TRUE)
}

# --- Define I/O Paths ---
# INPUT
annotated_noise_data_path <- "/group/sms029/mnieuwenh/gbm_noise_analysis/results/05_analyze_responsive_genes/05_noise_data_with_responsiveness.csv"
# OUTPUT
output_dir <- "/group/sms029/mnieuwenh/gbm_noise_analysis/results/07_GO_enrichment_analysis"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

celltype_output_dir <- file.path(output_dir, "celltype_analyses")
dir.create(celltype_output_dir, showWarnings = FALSE, recursive = TRUE)

gsea_figure_path <- file.path(output_dir, "07_GSEA_summary_plots_GLOBAL.png")
gsea_noise_results_path <- file.path(output_dir, "07_GSEA_by_noise_results_GLOBAL.csv")
summary_file_path <- file.path(output_dir, "07_GSEA_summary_GLOBAL.txt")

# --- 2. LOAD AND PREPARE DATA ---
message("Loading and preparing data for GSEA...")
if (!file.exists(annotated_noise_data_path)) stop("Annotated noise data file not found!")
all_data <- read.csv(annotated_noise_data_path)

# Keep only complete rows for noise and gene ID
all_data <- all_data %>%
  dplyr::filter(is.finite(variance.standardized), !is.na(gene))

# Helper: build a centered, named ranking vector (z-scored)
make_ranked_list <- function(df) {
  # df must have columns: gene, variance.standardized
  stats <- df %>%
    dplyr::group_by(gene) %>%
    dplyr::summarise(avg_noise = mean(variance.standardized, na.rm = TRUE), .groups = "drop") %>%
    dplyr::filter(!is.na(avg_noise))
  if (nrow(stats) == 0) return(NULL)

  # Center/scale to create both positive and negative scores
  z <- as.numeric(scale(stats$avg_noise, center = TRUE, scale = TRUE))
  names(z) <- stats$gene
  # Decreasing order for fgsea/gseGO
  z[order(z, decreasing = TRUE)]
}

# --- 3. GLOBAL GSEA ---
message("Performing Global GSEA on genes ranked by average (centered) noise...")
ranked_list_noise <- make_ranked_list(all_data)
if (is.null(ranked_list_noise) || length(ranked_list_noise) < 100) {
  warning("Not enough genes for global GSEA.")
  gsea_noise <- NULL
} else {
  set.seed(123)
  gsea_noise <- tryCatch({
    clusterProfiler::gseGO(
      geneList = ranked_list_noise,
      OrgDb = org.At.tair.db,
      ont = "BP",
      keyType = "TAIR",
      minGSSize = 10,
      maxGSSize = 5000,
      pvalueCutoff = 0.05,
      pAdjustMethod = "BH",
      eps = 0,
      verbose = FALSE
    )
  }, error = function(e) {
    message("Global gseGO failed: ", conditionMessage(e))
    NULL
  })
}

if (!is.null(gsea_noise) && nrow(gsea_noise@result) > 0) {
  p_noise <- enrichplot::dotplot(gsea_noise, showCategory = 20, split = ".sign") +
    facet_grid(~.sign) +
    labs(title = "Functional Profile of Gene Expression Noise (Global)")
  ggsave(gsea_figure_path, plot = p_noise, width = 12, height = 10, dpi = 300)
  write.csv(as.data.frame(gsea_noise@result), gsea_noise_results_path, row.names = FALSE)
} else {
  message("Global GSEA returned no significant categories.")
}

# --- 4. CELL-TYPE-SPECIFIC GSEA ---
message("Performing cell-type-specific GSEA...")
unique_cell_types <- unique(na.omit(all_data$base_cell_type))

for (cell_type in unique_cell_types) {
  message(paste("... GSEA for cell type:", cell_type))
  df_ct <- all_data %>% dplyr::filter(base_cell_type == cell_type)
  ranked_list_ct <- make_ranked_list(df_ct)

  if (is.null(ranked_list_ct) || length(ranked_list_ct) < 100) {
    message("...... Skipping due to insufficient genes.")
    next
  }

  set.seed(123)
  gsea_ct <- tryCatch({
    clusterProfiler::gseGO(
      geneList = ranked_list_ct,
      OrgDb = org.At.tair.db,
      ont = "BP",
      keyType = "TAIR",
      minGSSize = 10,
      maxGSSize = 5000,
      pvalueCutoff = 0.05,
      pAdjustMethod = "BH",
      eps = 0,
      verbose = FALSE
    )
  }, error = function(e) {
    message("...... gseGO failed for ", cell_type, ": ", conditionMessage(e))
    NULL
  })

  if (!is.null(gsea_ct) && nrow(gsea_ct@result) > 0) {
    p_ct <- enrichplot::dotplot(gsea_ct, showCategory = 20, split = ".sign") +
      facet_grid(~.sign) +
      labs(title = paste("Functional Profile of Noise in:", cell_type))
    ct_plot_path <- file.path(celltype_output_dir, paste0("07_GSEA_noise_", gsub(" ", "_", cell_type), ".png"))
    ggsave(ct_plot_path, plot = p_ct, width = 12, height = 10, dpi = 300)

    ct_data_path <- file.path(celltype_output_dir, paste0("07_GSEA_noise_results_", gsub(" ", "_", cell_type), ".csv"))
    write.csv(as.data.frame(gsea_ct@result), ct_data_path, row.names = FALSE)
  } else {
    message("...... No significant enrichment found for ", cell_type, ".")
  }
}

# --- 5. SUMMARY TEXT ---
sink(summary_file_path)
cat("=================================================================\n")
cat("         GLOBAL SUMMARY: 07_GO_enrichment_analysis.R\n")
cat("=================================================================\n\n")
cat(paste("Summary generated on:", Sys.Date()), "\n\n")
cat("Method: clusterProfiler::gseGO (BP ontology), keyType = 'TAIR'.\n")
cat("Ranking: genes ordered by centered (z-scored) average VST-corrected noise,\n")
cat("which yields both positive and negative statistics as recommended for GSEA.\n")
sink()

message(paste("Script 07 finished. Cell-type-specific outputs saved in:", celltype_output_dir))