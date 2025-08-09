# ==============================================================================
# SCRIPT: 05_analyze_responsive_genes.R
#
# PURPOSE:
#   This script identifies "responsive genes" (those upregulated in ERS+ cells)
#   and tests whether this functional state is associated with higher expression noise.
#   It performs this analysis globally and also on a per-cell-type basis.
#
#   **MODIFICATIONS**:
#   1. Saves full DE results (including avg_log2FC) for more powerful modeling.
#   2. Directly tests the hypothesis that responsive genes are enriched for TATA
#      boxes, providing a mechanistic link.
#   3. Adds a loop to perform DE analysis and visualization for each cell type
#      individually, saving results to a dedicated subdirectory.
# ==============================================================================

# --- 1. SETUP: Load libraries and define paths ---
packages_to_load <- c("Seurat", "dplyr", "ggplot2", "patchwork")
for (pkg in packages_to_load) {
  if (!requireNamespace(pkg, quietly = TRUE)) install.packages(pkg, repos = "https://cloud.r-project.org/")
  library(pkg, character.only = TRUE)
}

# --- Define I/O Paths ---
# INPUTS
seurat_object_path <- "/group/sms029/Oliva_dataset/integrated_col_trajectories_colonly.rds"
noise_data_path <- "/group/sms029/mnieuwenh/gbm_noise_analysis/results/03_calculate_expression_noise/03_gene_noise_metrics_complete.csv"
# OUTPUTS
output_dir <- "/group/sms029/mnieuwenh/gbm_noise_analysis/results/05_analyze_responsive_genes"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# Create the new subdirectory for cell type specific results
celltype_output_dir <- file.path(output_dir, "celltype_analyses")
dir.create(celltype_output_dir, showWarnings = FALSE, recursive = TRUE)

output_plot_path <- file.path(output_dir, "05_responsive_genes_noise_analysis_GLOBAL.png")
output_data_path <- file.path(output_dir, "05_noise_data_with_responsiveness.csv")
summary_file_path <- file.path(output_dir, "05_responsive_genes_summary_GLOBAL.txt")


# --- 2. DATA LOADING AND PREPARATION ---
message("Loading Seurat object and noise data...")
seurat_obj <- readRDS(seurat_object_path)
noise_data <- read.csv(noise_data_path)

# --- 3. IDENTIFY RESPONSIVE GENES (GLOBAL DE ANALYSIS) ---
message("Identifying responsive genes (ERS+ vs ERS-) across all cells...")

seurat_obj@meta.data <- seurat_obj@meta.data %>%
  mutate(
    identity_char = as.character(identity),
    ers_state = case_when(
      grepl("ERS\\+$", identity_char) ~ "ERS+",
      grepl("ERS-$", identity_char) ~ "ERS-",
      TRUE ~ "NA"
    ),
    base_cell_type = trimws(gsub("ERS\\+?/?-?$", "", identity_char))
  )

Idents(seurat_obj) <- seurat_obj$ers_state
de_results <- FindMarkers(seurat_obj, ident.1 = "ERS+", ident.2 = "ERS-",
                          logfc.threshold = 0, min.pct = 0.1, test.use = "wilcox")
de_results$gene <- rownames(de_results)

LOGFC_THRESHOLD <- 0.25
PVAL_ADJ_THRESHOLD <- 0.05
responsive_genes <- de_results %>%
  filter(p_val_adj < PVAL_ADJ_THRESHOLD & avg_log2FC > LOGFC_THRESHOLD) %>%
  pull(gene)

# --- 4. ANNOTATE NOISE DATA AND PERFORM GLOBAL STATS ---
message("Annotating noise data with global responsiveness status...")
noise_data_annotated <- noise_data %>%
  left_join(de_results %>% select(gene, avg_log2FC, p_val_adj), by = "gene") %>%
  mutate(is_responsive = gene %in% responsive_genes)

filtered_annotated_data <- noise_data_annotated %>%
  filter(mean_expr > 0.01, pct_expressing > 0.1, is.finite(variance.standardized))

stats_responsive_noise <- wilcox.test(variance.standardized ~ is_responsive, data = filtered_annotated_data)
tata_responsive_summary <- filtered_annotated_data %>%
  distinct(gene, .keep_all = TRUE) %>%
  filter(!is.na(has_TATA_box))
contingency_table <- table(tata_responsive_summary$is_responsive, tata_responsive_summary$has_TATA_box)
colnames(contingency_table) <- c("TATA-less", "TATA-containing")
rownames(contingency_table) <- c("Non-responsive", "Responsive")
chi_sq_test <- chisq.test(contingency_table)

# --- 5. VISUALIZATION (GLOBAL) ---
message("Generating global visualizations...")
publication_theme <- theme_bw() + theme(plot.title = element_text(hjust = 0.5, face = "bold"))
panel_a <- ggplot(filtered_annotated_data, aes(x = is_responsive, y = variance.standardized, fill = is_responsive)) +
  geom_violin(trim = FALSE, alpha = 0.6) + geom_boxplot(width = 0.1, outlier.shape = NA, fill = "white") +
  scale_x_discrete(labels = c("FALSE" = "Non-responsive", "TRUE" = "Responsive")) +
  scale_fill_manual(values = c("FALSE" = "cornflowerblue", "TRUE" = "coral1")) +
  labs(title = "A. Responsive Genes Exhibit Higher Noise (Global)", subtitle = paste0("Wilcoxon p=", format.pval(stats_responsive_noise$p.value, digits=2)), x = "Gene State", y = "Corrected Noise (VST)") +
  publication_theme + theme(legend.position = "none")

panel_b_data <- as.data.frame(contingency_table)
colnames(panel_b_data) <- c("Responsiveness", "PromoterType", "Count")
panel_b <- ggplot(panel_b_data, aes(x = Responsiveness, y = Count, fill = PromoterType)) +
  geom_bar(stat = "identity", position = "fill") + scale_fill_manual(values = c("TATA-less" = "grey70", "TATA-containing" = "#F0E442"), name = "Promoter Type") +
  scale_y_continuous(labels = scales::percent_format()) +
  labs(title = "B. Responsive Genes Enriched for TATA Boxes (Global)", subtitle = paste0("Chi-squared p=", format.pval(chi_sq_test$p.value, digits=2)), x = "Gene State", y = "Proportion of Genes") +
  publication_theme

final_figure <- panel_a | panel_b
ggsave(output_plot_path, plot = final_figure, width = 12, height = 6)
write.csv(noise_data_annotated, output_data_path, row.names = FALSE)

sink(summary_file_path)
cat("=================================================================\n"); cat("       GLOBAL STATISTICAL SUMMARY: 05_analyze_responsive_genes.R\n"); cat("=================================================================\n\n")
cat(paste("Summary generated on:", Sys.Date()), "\n\n"); cat("--- Differential Expression Analysis (ERS+ vs ERS-) ---\n")
cat(paste("LogFC Threshold:", LOGFC_THRESHOLD), "\n"); cat(paste("Adj P-val Threshold:", PVAL_ADJ_THRESHOLD), "\n")
cat(paste("Total 'Responsive' genes:", length(responsive_genes)), "\n\n"); cat("--- Noise Analysis --- \n"); print(stats_responsive_noise)
cat("\n\n--- TATA Box Enrichment Analysis ---\n"); cat("Contingency Table:\n"); print(contingency_table); cat("\nChi-squared Test:\n"); print(chi_sq_test)
sink()

# ==============================================================================
# --- 6. CELL-TYPE-SPECIFIC ANALYSIS OF RESPONSIVE GENES ---
# ==============================================================================
message("Starting cell-type-specific analysis of responsive genes...")

# Identify base cell types that have both ERS+ and ERS- cells
cell_types_for_de <- seurat_obj@meta.data %>%
  filter(ers_state %in% c("ERS+", "ERS-")) %>%
  group_by(base_cell_type) %>%
  summarise(n_states = n_distinct(ers_state)) %>%
  filter(n_states == 2) %>%
  pull(base_cell_type)

# Loop through each valid cell type
for (cell_type in cell_types_for_de) {
  message(paste("... Processing cell type:", cell_type))
  
  # Subset Seurat object to the current cell type
  seurat_subset <- subset(seurat_obj, subset = base_cell_type == cell_type)
  
  # Perform DE analysis within this cell type
  Idents(seurat_subset) <- seurat_subset$ers_state
  de_results_ct <- FindMarkers(seurat_subset, ident.1 = "ERS+", ident.2 = "ERS-", logfc.threshold = 0, min.pct = 0.1, test.use = "wilcox")
  de_results_ct$gene <- rownames(de_results_ct)
  
  responsive_genes_ct <- de_results_ct %>%
    filter(p_val_adj < PVAL_ADJ_THRESHOLD & avg_log2FC > LOGFC_THRESHOLD) %>%
    pull(gene)
    
  # Filter the main noise data for this cell type
  noise_data_ct <- filtered_annotated_data %>%
    filter(base_cell_type == cell_type) %>%
    mutate(is_responsive_ct = gene %in% responsive_genes_ct)
    
  # Check if there's enough data to proceed
  if(nrow(noise_data_ct) < 50 || n_distinct(noise_data_ct$is_responsive_ct) < 2) {
    message(paste("...... Skipping", cell_type, "due to insufficient data for comparison."))
    next
  }

  # Perform stats for this cell type
  stats_responsive_noise_ct <- wilcox.test(variance.standardized ~ is_responsive_ct, data = noise_data_ct)
  tata_summary_ct <- noise_data_ct %>% distinct(gene, .keep_all = TRUE) %>% filter(!is.na(has_TATA_box))
  
  if(n_distinct(tata_summary_ct$is_responsive_ct) < 2) {
      chi_sq_test_ct <- list(p.value = NA) # Cannot run test
  } else {
      contingency_table_ct <- table(tata_summary_ct$is_responsive_ct, tata_summary_ct$has_TATA_box)
      chi_sq_test_ct <- chisq.test(contingency_table_ct)
  }

  # Generate plots for this cell type
  p_a_ct <- ggplot(noise_data_ct, aes(x = is_responsive_ct, y = variance.standardized, fill = is_responsive_ct)) +
    geom_violin(trim = FALSE, alpha = 0.6) + geom_boxplot(width = 0.1, outlier.shape = NA, fill = "white") +
    scale_x_discrete(labels = c("FALSE" = "Non-resp.", "TRUE" = "Responsive")) +
    scale_fill_manual(values = c("FALSE" = "cornflowerblue", "TRUE" = "coral1")) +
    labs(subtitle = paste0("Wilcoxon p=", format.pval(stats_responsive_noise_ct$p.value, digits=2)), x = "Gene State", y = "Corrected Noise (VST)") +
    publication_theme + theme(legend.position = "none")

  if(n_distinct(tata_summary_ct$is_responsive_ct) < 2) {
      p_b_ct <- ggplot() + theme_void() + labs(title = "B. TATA Enrichment", subtitle="Not enough data")
  } else {
      panel_b_data_ct <- as.data.frame(contingency_table_ct)
      colnames(panel_b_data_ct) <- c("Responsiveness", "PromoterType", "Count")
      p_b_ct <- ggplot(panel_b_data_ct, aes(x = Responsiveness, y = Count, fill = PromoterType)) +
        geom_bar(stat = "identity", position = "fill") + scale_fill_manual(values = c("TATA-less" = "grey70", "TATA-containing" = "#F0E442"), name = "Promoter") +
        scale_y_continuous(labels = scales::percent_format()) +
        scale_x_discrete(labels = c("FALSE" = "Non-resp.", "TRUE" = "Responsive")) +
        labs(subtitle = paste0("Chi-sq p=", format.pval(chi_sq_test_ct$p.value, digits=2)), x = "Gene State", y = "Proportion") +
        publication_theme
  }

  # Combine and save the cell-type-specific plot
  final_figure_ct <- (p_a_ct | p_b_ct) + plot_annotation(title = paste("Responsiveness Analysis in:", cell_type))
  ct_plot_path <- file.path(celltype_output_dir, paste0("05_responsive_analysis_", gsub(" ", "_", cell_type), ".png"))
  ggsave(ct_plot_path, plot = final_figure_ct, width = 12, height = 6)
}

message(paste("Script finished. Cell-type-specific outputs saved in:", celltype_output_dir))