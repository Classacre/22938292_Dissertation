# Compare expression noise between gbM and TE-like methylation genes (Cahn and Bewick)
library(ggplot2)
if (!requireNamespace("ggridges", quietly = TRUE)) install.packages("ggridges", repos="https://cloud.r-project.org/")
if (!requireNamespace("ComplexUpset", quietly = TRUE)) install.packages("ComplexUpset", repos="https://cloud.r-project.org/")
if (!requireNamespace("pheatmap", quietly = TRUE)) install.packages("pheatmap", repos="https://cloud.r-project.org/")
library(ggridges)
library(ComplexUpset)
library(pheatmap)
library(gridExtra) # For marrangeGrob, tableGrob
library(grid)      # For grid.newpage and grid.draw
library(dplyr)     # For data manipulation (e.g., group_by, summarise)
library(tidyr)     # For pivot_longer
library(stringr)   ## FIX: Added stringr for str_to_title function used later

# Define the output directory
output_dir <- "/group/sms029/mnieuwenh/gbm_noise_analysis/results/high_low_noise"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE) # Ensure directory exists

# Load expression matrix
dat_path <- "/group/sms029/mnieuwenh/gbm_noise_analysis/expression_matrix.csv"
expr_mat <- as.matrix(read.csv(dat_path, row.names=1, check.names=FALSE))

# Load combined methylation annotation
gene_anno <- read.csv("/group/sms029/mnieuwenh/Methylation_Data/combined_methylation_data.csv", stringsAsFactors=FALSE)
rownames(gene_anno) <- gene_anno$Gene_ID

# Annotate genes for each methylation status
genes <- rownames(expr_mat)
anno <- gene_anno[genes, ]
cahn_gbm <- !is.na(anno$Cahn_Methylation_status) & anno$Cahn_Methylation_status == "gbM"
bewick_gbm <- !is.na(anno$Bewick_Classification) & anno$Bewick_Classification == "gbM"
bewick_te_like <- !is.na(anno$Bewick_Classification) & (anno$Bewick_Classification == "mCHG" | anno$Bewick_Classification == "mCHH")

# Calculate mean, variance, CV for each gene (overall)
mean_expr <- rowMeans(expr_mat)
var_expr <- apply(expr_mat, 1, var)
cv_expr <- sqrt(var_expr) / mean_expr

# Create results data frame (overall gene properties)
results <- data.frame(
  gene = genes,
  cahn_gbm = cahn_gbm,
  bewick_gbm = bewick_gbm,
  bewick_te_like = bewick_te_like,
  mean_expr = mean_expr,
  var_expr = var_expr,
  cv_expr = cv_expr
)

# Assign methylation groups for plotting
results$methylation_group <- NA
results$methylation_group[!is.na(anno$Cahn_Methylation_status) & anno$Cahn_Methylation_status == "gbM"] <- "Cahn gbM"
results$methylation_group[!is.na(anno$Bewick_Classification) & anno$Bewick_Classification == "gbM"] <- "Bewick gbM"
results$methylation_group[!is.na(anno$Bewick_Classification) & anno$Bewick_Classification == "mCHH"] <- "CHH"
results$methylation_group[!is.na(anno$Bewick_Classification) & anno$Bewick_Classification == "mCHG"] <- "CHG"
results$methylation_group[!is.na(anno$Cahn_Methylation_status) & anno$Cahn_Methylation_status == "TE-like methylation"] <- "TE-Like Methylation"
results$methylation_group[!is.na(anno$Cahn_Methylation_status) & anno$Cahn_Methylation_status == "Unmethylated" & is.na(results$methylation_group)] <- "Unmethylated"
results$methylation_group[is.na(results$methylation_group)] <- "Non-gbM"
results$methylation_group <- factor(results$methylation_group, levels=c("Unmethylated", "Non-gbM", "Bewick gbM", "Cahn gbM", "CHH", "CHG", "TE-Like Methylation"))

# Assign Cahn and Bewick groups for plotting
results$cahn_group <- NA
results$cahn_group[!is.na(anno$Cahn_Methylation_status) & anno$Cahn_Methylation_status == "gbM"] <- "gbM"
results$cahn_group[!is.na(anno$Cahn_Methylation_status) & anno$Cahn_Methylation_status == "TE-like methylation"] <- "TE-like"
results$cahn_group[!is.na(anno$Cahn_Methylation_status) & anno$Cahn_Methylation_status == "Unmethylated"] <- "Unmethylated"
results$cahn_group <- factor(results$cahn_group, levels=c("Unmethylated", "gbM", "TE-like"))
results$bewick_group <- NA
results$bewick_group[!is.na(anno$Bewick_Classification) & anno$Bewick_Classification == "gbM"] <- "gbM"
results$bewick_group[!is.na(anno$Bewick_Classification) & anno$Bewick_Classification == "mCHH"] <- "CHH"
results$bewick_group[!is.na(anno$Bewick_Classification) & anno$Bewick_Classification == "mCHG"] <- "CHG"
results$bewick_group[!is.na(anno$Bewick_Classification) & anno$Bewick_Classification == "Unmethylated"] <- "Unmethylated"
results$bewick_group <- factor(results$bewick_group, levels=c("Unmethylated", "gbM", "CHH", "CHG"))

# Add sd_expr column before filtering for plotting
results$sd_expr <- sqrt(results$var_expr)

# Filter out non-finite CV values for plotting (overall data)
plot_results <- results[is.finite(results$cv_expr), ]
# Remove extreme outliers (e.g., above 99th percentile) for plotting
cv_cap <- quantile(plot_results$cv_expr, 0.99, na.rm=TRUE)
plot_results <- plot_results[plot_results$cv_expr <= cv_cap, ]

# --- James Lloyd Filter (Mean Expression Percentile Cutoff) ---
# Function to apply James Lloyd filter and basic cleaning
apply_jl_filter <- function(data_df, percentile_cutoff = NULL, mean_cutoff = NULL) {
  filtered_df <- data_df
  if (!is.null(percentile_cutoff)) {
    cutoff_value <- quantile(data_df$mean_expr, percentile_cutoff, na.rm=TRUE)
    filtered_df <- filtered_df[filtered_df$mean_expr >= cutoff_value, ]
  } else if (!is.null(mean_cutoff)) {
    filtered_df <- filtered_df[filtered_df$mean_expr >= mean_cutoff, ]
  }
  
  # Apply basic cleaning: remove non-finite CVs and extreme outliers
  filtered_df <- filtered_df[is.finite(filtered_df$cv_expr), ]
  if (nrow(filtered_df) > 0) { # Ensure there's data before calculating quantile
    cv_cap <- quantile(filtered_df$cv_expr, 0.99, na.rm=TRUE)
    filtered_df <- filtered_df[filtered_df$cv_expr <= cv_cap, ]
  }
  return(filtered_df)
}

# Apply filters
plot_results_jl_10_percent_filtered <- apply_jl_filter(results, percentile_cutoff = 0.10)
plot_results_jl_20_percent_filtered <- apply_jl_filter(results, percentile_cutoff = 0.20)
plot_results_mean_gt_0.5_filtered <- apply_jl_filter(results, mean_cutoff = 0.5)
plot_results_mean_gt_1_filtered <- apply_jl_filter(results, mean_cutoff = 1)

# Debugging: Check counts of Bewick groups in filtered dataframes
print("--- Debugging Bewick Group Counts ---")
print("Counts of Bewick groups in plot_results (baseline):")
print(table(plot_results$bewick_group, useNA = "ifany"))

print("Counts of Bewick groups in plot_results_jl_10_percent_filtered:")
print(table(plot_results_jl_10_percent_filtered$bewick_group, useNA = "ifany"))

print("Counts of Bewick groups in plot_results_jl_20_percent_filtered:")
print(table(plot_results_jl_20_percent_filtered$bewick_group, useNA = "ifany"))

print("Counts of Bewick groups in plot_results_mean_gt_0.5_filtered:")
print(table(plot_results_mean_gt_0.5_filtered$bewick_group, useNA = "ifany"))

print("Counts of Bewick groups in plot_results_mean_gt_1_filtered:")
print(table(plot_results_mean_gt_1_filtered$bewick_group, useNA = "ifany"))
print("-------------------------------------")


# Plots (original individual plots - not saved to PDF but kept for Rmd)
gbm_labels <- c("FALSE"="Non-gbM", "TRUE"="gbM")
te_labels <- c("FALSE"="Other", "TRUE"="TE-like")
p1 <- ggplot(results, aes(x=factor(cahn_gbm), y=cv_expr, fill=factor(cahn_gbm))) +
  geom_boxplot(outlier.size=0.5) + ## FIX: `outlier.linewidth` is deprecated, using `outlier.size`
  scale_x_discrete(labels=gbm_labels) +
  labs(title="Expression Noise (CV) by Cahn gbM Status", x="Cahn gbM", y="CV") +
  theme_bw()
p2 <- ggplot(results, aes(x=factor(bewick_gbm), y=cv_expr, fill=factor(bewick_gbm))) +
  geom_boxplot(outlier.size=0.5) + ## FIX: `outlier.linewidth` is deprecated, using `outlier.size`
  scale_x_discrete(labels=gbm_labels, drop = FALSE) + # Added drop=FALSE
  labs(title="Expression Noise (CV) by Bewick gbM Status", x="Bewick gbM", y="CV") +
  theme_bw()
p3 <- ggplot(results, aes(x=factor(bewick_te_like), y=cv_expr, fill=factor(bewick_te_like))) +
  geom_boxplot(outlier.size=0.5) + ## FIX: `outlier.linewidth` is deprecated, using `outlier.size`
  scale_x_discrete(labels=te_labels) +
  labs(title="Expression Noise (CV) by Bewick TE-like Status", x="Bewick TE-like", y="CV") +
  theme_bw()
p_combined <- ggplot(results, aes(x=methylation_group, y=cv_expr, color=methylation_group)) +
  geom_jitter(width=0.25, alpha=0.5, size=0.7) +
  labs(title="Univariant Scatter Plot: Expression Noise (CV) by Methylation Group", x="Methylation Group", y="CV") +
  theme_bw() +
  theme(axis.text.x = element_text(angle=45, hjust=1))

# Save plot objects for Rmd
save(p1, p2, p3, p_combined, plot_results, file=file.path(output_dir, "gbm_noise_plots.RData"))

# --- H2A.Z annotation ---
h2az_depleted <- !is.na(anno$H2AZ_Depleted) & anno$H2AZ_Depleted == TRUE
h2az_enriched <- !is.na(anno$H2AZ_Enriched) & anno$H2AZ_Enriched == TRUE
h2az_group <- rep(NA, length(genes))
h2az_group[h2az_depleted] <- "H2A.Z-Depleted"
h2az_group[h2az_enriched] <- "H2A.Z-Enriched"
h2az_group[is.na(h2az_group)] <- "Other/NA"
results$h2az_group <- factor(h2az_group, levels=c("H2A.Z-Depleted", "H2A.Z-Enriched", "Other/NA"))

# --- Boxplots ---
p_h2az <- ggplot(results, aes(x=h2az_group, y=cv_expr, fill=h2az_group)) +
  geom_boxplot(outlier.size=0.5) + ## FIX: `outlier.linewidth` is deprecated, using `outlier.size`
  labs(title="Expression Noise (CV) by H2A.Z Group", x="H2A.Z Group", y="CV") +
  theme_bw()
ggsave(file.path(output_dir, "gbm_noise_boxplot_h2az.png"), plot=p_h2az, width=8, height=6)

# --- Scatter plots: mean vs SD and mean vs CV by H2A.Z group ---
p_mean_sd_h2az <- ggplot(results, aes(x=mean_expr, y=sd_expr, color=h2az_group)) +
  geom_point(alpha=0.5, size=0.7) +
  labs(title="Mean Expression vs SD by H2A.Z Group", x="Mean Expression", y="Standard Deviation") +
  theme_bw()
ggsave(file.path(output_dir, "mean_vs_sd_by_h2az_group.png"), plot=p_mean_sd_h2az, width=8, height=6)

p_mean_cv_h2az <- ggplot(results, aes(x=mean_expr, y=cv_expr, color=h2az_group)) +
  geom_point(alpha=0.5, size=0.7) +
  labs(title="Mean Expression vs CV by H2A.Z Group", x="Mean Expression", y="Coefficient of Variation (CV)") +
  theme_bw()
ggsave(file.path(output_dir, "mean_vs_cv_by_h2az_group.png"), plot=p_mean_cv_h2az, width=8, height=6)

# --- Combined methylation + H2A.Z group ---
results$meth_h2az_group <- paste0(results$methylation_group, "+", results$h2az_group)
p_meth_h2az <- ggplot(results, aes(x=meth_h2az_group, y=cv_expr, fill=meth_h2az_group)) +
  geom_boxplot(outlier.size=0.5) + ## FIX: `outlier.linewidth` is deprecated, using `outlier.size`
  labs(title="Expression Noise (CV) by Methylation + H2A.Z Group", x="Methylation + H2A.Z Group", y="CV") +
  theme_bw() +
  theme(axis.text.x = element_text(angle=45, hjust=1))
ggsave(file.path(output_dir, "gbm_noise_boxplot_meth_h2az.png"), width=12, height=6)

# --- Violin plots ---
p_cahn_violin <- ggplot(plot_results[!is.na(plot_results$cahn_group),], aes(x=cahn_group, y=cv_expr, fill=cahn_group)) +
  geom_violin(trim=FALSE, scale="width") +
  geom_boxplot(width=0.1, outlier.size=0.5, outlier.shape=16, outlier.alpha=0.3, fill="white") + ## FIX: `outlier.linewidth` is deprecated, using `outlier.size`
  labs(title="Expression Noise (CV) by Cahn Methylation Group (Violin)", x="Cahn Methylation Group", y="CV") +
  theme_bw() +
  theme(axis.text.x = element_text(angle=45, hjust=1))
ggsave(file.path(output_dir, "gbm_noise_violin_cahn.png"), plot=p_cahn_violin, width=8, height=6)

p_bewick_violin <- ggplot(plot_results[!is.na(plot_results$bewick_group),], aes(x=bewick_group, y=cv_expr, fill=bewick_group)) +
  geom_violin(trim=FALSE, scale="width") +
  geom_boxplot(width=0.1, outlier.size=0.5, outlier.shape=16, outlier.alpha=0.3, fill="white") + ## FIX: `outlier.linewidth` is deprecated, using `outlier.size`
  scale_x_discrete(drop = FALSE) + # Added for Bewick plots
  labs(title="Expression Noise (CV) by Bewick Methylation Group (Violin)", x="Bewick Methylation Group", y="CV") +
  theme_bw() +
  theme(axis.text.x = element_text(angle=45, hjust=1))
ggsave(file.path(output_dir, "gbm_noise_violin_bewick.png"), plot=p_bewick_violin, width=8, height=6)

# --- Ridge plots ---
p_cahn_ridge <- ggplot(plot_results[!is.na(plot_results$cahn_group),], aes(x=cv_expr, y=cahn_group, fill=cahn_group)) +
  ggridges::geom_density_ridges(scale=2.5, alpha=0.7, rel_min_height=0.01) +
  labs(title="CV Distribution by Cahn Methylation Group (Ridge)", x="CV", y="Cahn Methylation Group") +
  theme_bw()
ggsave(file.path(output_dir, "gbm_noise_ridge_cahn.png"), plot=p_cahn_ridge, width=8, height=8)

p_bewick_ridge <- ggplot(plot_results[!is.na(plot_results$bewick_group),], aes(x=cv_expr, y=bewick_group, fill=bewick_group)) +
  ggridges::geom_density_ridges(scale=2.5, alpha=0.7, rel_min_height=0.01) +
  scale_y_discrete(drop = FALSE) + # Added for Bewick plots
  labs(title="CV Distribution by Bewick Methylation Group (Ridge)", x="CV", y="Bewick Methylation Group") +
  theme_bw()
ggsave(file.path(output_dir, "gbm_noise_ridge_bewick.png"), plot=p_bewick_ridge, width=8, height=8)

# --- UpSet plot: Overlap of high-noise, gbM, and H2A.Z-depleted genes ---
upset_data <- data.frame(
  High_Noise = results$cv_expr >= quantile(results$cv_expr, 0.9, na.rm=TRUE),
  gbM = results$cahn_gbm | results$bewick_gbm,
  H2AZ_Depleted = results$h2az_group == "H2A.Z-Depleted"
)
upset_features <- c("High_Noise","gbM","H2AZ_Depleted")
## FIX: Assign the plot object to a variable before saving with ggsave
p_upset <- ComplexUpset::upset(
  upset_data, 
  intersect = upset_features, 
  width_ratio=0.1, 
  set_sizes=ComplexUpset::upset_set_size() + theme(axis.text.x = element_text(angle=45, hjust=1, size=14))
)
ggsave(file.path(output_dir, "upset_highnoise_gbm_h2az.png"), plot = p_upset, width=14, height=7)


# --- Heatmap: Top 30 high-noise genes (CV) ---
if (exists("plot_results")) {
  topN <- 30
  top_genes <- head(plot_results[order(-plot_results$cv_expr),], topN)
  if (nrow(top_genes) > 0) {
    mat <- expr_mat[top_genes$gene, , drop=FALSE]
    pheatmap::pheatmap(mat, cluster_rows=TRUE, cluster_cols=TRUE, show_rownames=TRUE, show_colnames=FALSE, main=paste("Top", topN, "High-Noise Genes (CV)"), filename=file.path(output_dir, "heatmap_top30_highnoise.png"), width=16, height=8)
  }
}

# Make statistics available for Rmd report
stat_cahn <- wilcox.test(cv_expr ~ cahn_gbm, data=results)
stat_bewick <- wilcox.test(cv_expr ~ bewick_gbm, data=results)
stat_te_like <- wilcox.test(cv_expr ~ bewick_te_like, data=results)
stat_h2az <- wilcox.test(cv_expr ~ h2az_group, data=results, subset=h2az_group %in% c("H2A.Z-Depleted", "H2A.Z-Enriched"))

if (exists("stat_cahn") && exists("stat_bewick") && exists("stat_te_like") && exists("stat_h2az")) {
  save(stat_cahn, stat_bewick, stat_te_like, stat_h2az, file=file.path(output_dir, "gbm_noise_stats.RData"))
}

# --- Celltype and Lineage annotation ---
# Load cell metadata for celltype and lineage groupings
cell_metadata <- read.csv("/group/sms029/mnieuwenh/seurat_metadata/seurat_metadata_full.csv", row.names=1, check.names=FALSE)

# Match cells between expression matrix and metadata
common_cells <- intersect(colnames(expr_mat), rownames(cell_metadata))
expr_mat <- expr_mat[, common_cells]
cell_metadata <- cell_metadata[common_cells, ]

# Calculate mean expression per gene per celltype and lineage (robustly)
gene_means_celltype_list <- lapply(unique(cell_metadata$identity), function(ct) {
  cells <- rownames(cell_metadata)[cell_metadata$identity == ct]
  if (length(cells) < 1) return(rep(NA, nrow(expr_mat)))
  rowMeans(expr_mat[, cells, drop=FALSE])
})
gene_means_celltype <- do.call(cbind, gene_means_celltype_list)
colnames(gene_means_celltype) <- unique(cell_metadata$identity)
rownames(gene_means_celltype) <- rownames(expr_mat) # Ensure row names are genes

gene_means_lineage_list <- lapply(unique(cell_metadata$lineage), function(ln) {
  cells <- rownames(cell_metadata)[cell_metadata$lineage == ln]
  if (length(cells) < 1) return(rep(NA, nrow(expr_mat)))
  rowMeans(expr_mat[, cells, drop=FALSE])
})
gene_means_lineage <- do.call(cbind, gene_means_lineage_list)
colnames(gene_means_lineage) <- unique(cell_metadata$lineage)
rownames(gene_means_lineage) <- rownames(expr_mat) # Ensure row names are genes


# Calculate mean, variance, CV for each gene within each celltype and lineage
gene_data_by_celltype <- list()
for (ct in unique(cell_metadata$identity)) {
  if (is.na(ct) || ct == "") next
  cells_in_ct <- rownames(cell_metadata)[cell_metadata$identity == ct]
  if (length(cells_in_ct) > 1) { # Need at least 2 cells to calculate variance
    expr_mat_ct <- expr_mat[, cells_in_ct, drop=FALSE]
    mean_expr_ct <- rowMeans(expr_mat_ct)
    var_expr_ct <- apply(expr_mat_ct, 1, var)
    cv_expr_ct <- sqrt(var_expr_ct) / mean_expr_ct

    # Combine with gene annotations from 'results' (which has cahn_group, bewick_group)
    temp_df <- data.frame(
      gene = rownames(expr_mat_ct),
      mean_expr_ct = mean_expr_ct,
      var_expr_ct = var_expr_ct,
      cv_expr_ct = cv_expr_ct
    ) %>%
    left_join(results %>% dplyr::select(gene, cahn_group, bewick_group), by = "gene") %>%
    filter(is.finite(cv_expr_ct)) # Filter out non-finite CVs for this celltype

    # Remove extreme outliers for this celltype's CV distribution (99th percentile)
    cv_cap_ct <- quantile(temp_df$cv_expr_ct, 0.99, na.rm=TRUE)
    temp_df <- temp_df[temp_df$cv_expr_ct <= cv_cap_ct, ]

    gene_data_by_celltype[[ct]] <- temp_df
  }
}

gene_data_by_lineage <- list()
for (ln in unique(cell_metadata$lineage)) {
  if (is.na(ln) || ln == "") next
  cells_in_ln <- rownames(cell_metadata)[cell_metadata$lineage == ln]
  if (length(cells_in_ln) > 1) { # Need at least 2 cells to calculate variance
    expr_mat_ln <- expr_mat[, cells_in_ln, drop=FALSE]
    mean_expr_ln <- rowMeans(expr_mat_ln)
    var_expr_ln <- apply(expr_mat_ln, 1, var)
    cv_expr_ln <- sqrt(var_expr_ln) / mean_expr_ln

    # Combine with gene annotations from 'results'
    temp_df <- data.frame(
      gene = rownames(expr_mat_ln),
      mean_expr_ln = mean_expr_ln,
      var_expr_ln = var_expr_ln,
      cv_expr_ln = cv_expr_ln
    ) %>%
    left_join(results %>% dplyr::select(gene, cahn_group, bewick_group), by = "gene") %>%
    filter(is.finite(cv_expr_ln)) # Filter out non-finite CVs for this lineage

    # Remove extreme outliers for this lineage's CV distribution
    cv_cap_ln <- quantile(temp_df$cv_expr_ln, 0.99, na.rm=TRUE)
    temp_df <- temp_df[temp_df$cv_expr_ln <= cv_cap_ln, ]

    gene_data_by_lineage[[ln]] <- temp_df
  }
}

# Define celltype_levels and lineage_levels for lapply calls
celltype_levels <- unique(cell_metadata$identity)
lineage_levels <- unique(cell_metadata$lineage)

# Debugging: Check counts of Bewick groups in celltype-specific dataframes
print("--- Debugging Bewick Group Counts in Celltype Data ---")
for (ct in celltype_levels) {
  d_ct <- gene_data_by_celltype[[ct]]
  if (!is.null(d_ct) && nrow(d_ct) > 0) {
    print(paste0("Counts of Bewick groups in gene_data_by_celltype for ", ct, ":"))
    print(table(d_ct$bewick_group, useNA = "ifany"))
    
    d_filtered_jl_ct_10 <- d_ct %>% filter(gene %in% plot_results_jl_10_percent_filtered$gene)
    if (nrow(d_filtered_jl_ct_10) > 0) {
      print(paste0("Counts of Bewick groups in gene_data_by_celltype for ", ct, " (JL10 filtered):"))
      print(table(d_filtered_jl_ct_10$bewick_group, useNA = "ifany"))
    }
    
    d_filtered_jl_ct_20 <- d_ct %>% filter(gene %in% plot_results_jl_20_percent_filtered$gene)
    if (nrow(d_filtered_jl_ct_20) > 0) {
      print(paste0("Counts of Bewick groups in gene_data_by_celltype for ", ct, " (JL20 filtered):"))
      print(table(d_filtered_jl_ct_20$bewick_group, useNA = "ifany"))
    }
    
    d_filtered_mean_gt_0.5_ct <- d_ct %>% filter(gene %in% plot_results_mean_gt_0.5_filtered$gene)
    if (nrow(d_filtered_mean_gt_0.5_ct) > 0) {
      print(paste0("Counts of Bewick groups in gene_data_by_celltype for ", ct, " (Mean > 0.5 filtered):"))
      print(table(d_filtered_mean_gt_0.5_ct$bewick_group, useNA = "ifany"))
    }
    
    d_filtered_mean_gt_1_ct <- d_ct %>% filter(gene %in% plot_results_mean_gt_1_filtered$gene)
    if (nrow(d_filtered_mean_gt_1_ct) > 0) {
      print(paste0("Counts of Bewick groups in gene_data_by_celltype for ", ct, " (Mean > 1 filtered):"))
      print(table(d_filtered_mean_gt_1_ct$bewick_group, useNA = "ifany"))
    }
  }
}
print("-------------------------------------")

# Debugging: Check counts of Bewick groups in lineage-specific dataframes
print("--- Debugging Bewick Group Counts in Lineage Data ---")
for (ln in lineage_levels) {
  d_ln <- gene_data_by_lineage[[ln]]
  if (!is.null(d_ln) && nrow(d_ln) > 0) {
    print(paste0("Counts of Bewick groups in gene_data_by_lineage for ", ln, ":"))
    print(table(d_ln$bewick_group, useNA = "ifany"))

    d_filtered_jl_ln_10 <- d_ln %>% filter(gene %in% plot_results_jl_10_percent_filtered$gene)
    if (nrow(d_filtered_jl_ln_10) > 0) {
      print(paste0("Counts of Bewick groups in gene_data_by_lineage for ", ln, " (JL10 filtered):"))
      print(table(d_filtered_jl_ln_10$bewick_group, useNA = "ifany"))
    }

    d_filtered_jl_ln_20 <- d_ln %>% filter(gene %in% plot_results_jl_20_percent_filtered$gene)
    if (nrow(d_filtered_jl_ln_20) > 0) {
      print(paste0("Counts of Bewick groups in gene_data_by_lineage for ", ln, " (JL20 filtered):"))
      print(table(d_filtered_jl_ln_20$bewick_group, useNA = "ifany"))
    }

    d_filtered_mean_gt_0.5_ln <- d_ln %>% filter(gene %in% plot_results_mean_gt_0.5_filtered$gene)
    if (nrow(d_filtered_mean_gt_0.5_ln) > 0) {
      print(paste0("Counts of Bewick groups in gene_data_by_lineage for ", ln, " (Mean > 0.5 filtered):"))
      print(table(d_filtered_mean_gt_0.5_ln$bewick_group, useNA = "ifany"))
    }

    d_filtered_mean_gt_1_ln <- d_ln %>% filter(gene %in% plot_results_mean_gt_1_filtered$gene)
    if (nrow(d_filtered_mean_gt_1_ln) > 0) {
      print(paste0("Counts of Bewick groups in gene_data_by_lineage for ", ln, " (Mean > 1 filtered):"))
      print(table(d_filtered_mean_gt_1_ln$bewick_group, useNA = "ifany"))
    }
  }
}
print("-------------------------------------")


# New helper function to run specific tests and format results for a table
run_specific_group_stats <- function(data_df, value_col, group_col_name, context_name, filter_tag, specific_context) {
  results_list <- list()
  
  # Ensure the data_df is not empty and has the required columns
  if (is.null(data_df) || nrow(data_df) == 0 || !value_col %in% colnames(data_df) || !group_col_name %in% colnames(data_df)) {
    return(NULL)
  }
  
  # Filter out NA values for the group_col and value_col
  d_test <- data_df[!is.na(data_df[[group_col_name]]) & !is.na(data_df[[value_col]]),]
  
  # Define groups for comparison based on group_col_name
  base_group <- "Unmethylated"
  compare_groups <- c()
  
  if (group_col_name == "cahn_group") {
    compare_groups <- c("gbM", "TE-like")
  } else if (group_col_name == "bewick_group") {
    # For Bewick, compare against gbM, CHH, and CHG individually
    compare_groups <- c("gbM", "CHH", "CHG")
  } else {
    # If the group_col_name is not recognized for specific comparisons, skip
    return(NULL)
  }
  
  # Kruskal-Wallis test for overall significance among all relevant groups
  all_relevant_groups <- unique(d_test[[group_col_name]])[unique(d_test[[group_col_name]]) %in% c(base_group, compare_groups)]
  d_kruskal <- d_test[d_test[[group_col_name]] %in% all_relevant_groups, ]
  
  if (nrow(d_kruskal) > 0) {
    d_kruskal[[group_col_name]] <- droplevels(d_kruskal[[group_col_name]]) # Drop unused levels
  
    if (length(unique(d_kruskal[[group_col_name]])) >= 2 && all(table(d_kruskal[[group_col_name]]) >= 2)) {
      kruskal_result <- tryCatch({
        kruskal.test(d_kruskal[[value_col]] ~ d_kruskal[[group_col_name]])
      }, error = function(e) {
        warning(paste("Kruskal-Wallis error for", context_name, filter_tag, specific_context, ":", e$message))
        return(NULL)
      })
      if (!is.null(kruskal_result)) {
        results_list[[length(results_list) + 1]] <- data.frame(
          Context = context_name,
          Filter_Tag = filter_tag,
          Specific_Context = specific_context,
          Group_Type = group_col_name,
          Comparison = "Overall (Kruskal-Wallis)",
          P_value = kruskal_result$p.value,
          stringsAsFactors = FALSE
        )
      }
    }
  }


  # Pairwise Wilcoxon tests against base_group
  for (comp_group in compare_groups) {
    group1_data <- d_test[[value_col]][d_test[[group_col_name]] == base_group]
    group2_data <- d_test[[value_col]][d_test[[group_col_name]] == comp_group]
    
    if (length(group1_data) >= 2 && length(group2_data) >= 2) {
      wilcox_result <- tryCatch({
        wilcox.test(group1_data, group2_data)
      }, error = function(e) {
        warning(paste("Wilcoxon error for", base_group, "vs", comp_group, ":", e$message))
        return(NULL)
      })
      
      if (!is.null(wilcox_result)) {
        results_list[[length(results_list) + 1]] <- data.frame(
          Context = context_name,
          Filter_Tag = filter_tag,
          Specific_Context = specific_context,
          Group_Type = group_col_name,
          Comparison = paste(comp_group, "vs", base_group),
          P_value = wilcox_result$p.value,
          stringsAsFactors = FALSE
        )
      }
    }
  }
  
  if (length(results_list) > 0) {
    return(do.call(rbind, results_list))
  } else {
    return(NULL)
  }
}

# Helper function to generate and save a PDF with plots and a statistics table
save_grid_pdf_with_stats <- function(plot_list, stat_data_list, filename, title_text, nrow_grid = 2, ncol_grid = 3) {
  valid_plots <- plot_list[!sapply(plot_list, is.null)]
  valid_stats <- stat_data_list[!sapply(stat_data_list, is.null)]

  # Only open PDF if there are plots or stats to write
  if (length(valid_plots) == 0 && length(valid_stats) == 0) {
    message(paste("Skipping PDF generation for", filename, "as no plots or stats are available."))
    return(invisible(NULL))
  }

  pdf(filename, width=18, height=12)

  # Arrange plots in a grid and paginate
  if (length(valid_plots) > 0) {
    ## FIX: Use unname() to ensure marrangeGrob receives an unnamed list
    grob_pages <- marrangeGrob(unname(valid_plots), nrow=nrow_grid, ncol=ncol_grid, top=title_text)
    # Print each page
    for (page in grob_pages) {
      grid.newpage()
      grid.draw(page)
    }
  } else {
    grid.newpage()
    grid.text("No plots to display for this category.", gp = gpar(fontsize = 12))
  }

  # Generate and print statistics table
  if (length(valid_stats) > 0) {
    combined_stats_df <- bind_rows(valid_stats)

    if (nrow(combined_stats_df) > 0) {
        combined_stats_df$P_value_Formatted <- formatC(combined_stats_df$P_value, format="e", digits=4)
        table_to_display <- combined_stats_df %>% dplyr::select(-P_value)

        ## FIX: Simplified and corrected the p-value coloring logic
        # Create a color matrix based on the table that will actually be displayed
        colors <- matrix("black", nrow = nrow(table_to_display), ncol = ncol(table_to_display))
        p_col_idx <- which(colnames(table_to_display) == "P_value_Formatted")
        
        if (length(p_col_idx) > 0) {
            # Use the original numeric p-value from combined_stats_df for the logical condition
            significant_rows <- which(combined_stats_df$P_value < 0.05)
            if (length(significant_rows) > 0) {
                colors[significant_rows, p_col_idx] <- "red"
            }
        }

        my_theme <- ttheme_default(
            core = list(fg_params = list(col = colors, fontsize = 8)),
            colhead = list(fg_params = list(col = "black", fontsize = 9, fontface = "bold")),
            rowhead = list(fg_params = list(col = "black", fontsize = 8, fontface = "bold"))
        )
        
        # Paginate the table if it's too long
        rows_per_page <- 30
        num_pages <- ceiling(nrow(table_to_display) / rows_per_page)
        table_title_grob <- textGrob(paste0("Statistical Comparisons for ", title_text), gp = gpar(fontsize = 14, fontface = "bold"))

        for (p in 1:num_pages) {
            start_row <- (p - 1) * rows_per_page + 1
            end_row <- min(p * rows_per_page, nrow(table_to_display))
            
            page_data <- table_to_display[start_row:end_row, , drop = FALSE]
            page_theme <- my_theme
            page_theme$core$fg_params$col <- colors[start_row:end_row, , drop = FALSE]

            table_grob_page <- tableGrob(page_data, rows = NULL, theme = page_theme)
            
            grid.newpage()
            grid.draw(arrangeGrob(table_title_grob, table_grob_page, ncol = 1, heights = unit.c(grobHeight(table_title_grob) + unit(1, "line"), unit(1,"npc") - (grobHeight(table_title_grob) + unit(1, "line")))))
        }
    } else {
      grid.newpage()
      grid.text("No statistical results available for these plots.", gp = gpar(fontsize = 12))
    }
  } else {
    grid.newpage()
    grid.text("No statistical results available for these plots.", gp = gpar(fontsize = 12))
  }

  dev.off()
}


# --- Boxplots by celltype (Cahn group) - Mean Expression ---
plots_cahn_celltype_unfiltered <- lapply(celltype_levels, function(ct) {
  d <- plot_results # Overall unfiltered gene set
  d$cell_mean <- gene_means_celltype[d$gene, ct] # Add celltype-specific mean
  plot_data <- d[!is.na(d$cahn_group),]
  if (nrow(plot_data) == 0 || all(is.na(plot_data$cell_mean))) return(NULL) # Return NULL if no data
  ggplot(plot_data, aes(x=factor(cahn_group), y=cell_mean, fill=factor(cahn_group))) +
    geom_boxplot(outlier.size=0.5) + ## FIX: `outlier.linewidth` is deprecated, using `outlier.size`
    labs(title=paste("Cahn Group in", ct), x="Cahn Group", y="Mean Expression") +
    theme_bw() +
    theme(legend.position = "none", axis.text.x = element_text(angle=45, hjust=1))
})
save_grid_pdf_with_stats(plots_cahn_celltype_unfiltered, list(), file.path(output_dir, "grid_mean_expr_cahn_by_celltype_unfiltered.pdf"), "Mean Expression by Cahn Group in Celltypes (Unfiltered)")

plots_cahn_celltype_jl_10_percent_filtered <- lapply(celltype_levels, function(ct) {
  d <- plot_results_jl_10_percent_filtered
  d$cell_mean <- gene_means_celltype[d$gene, ct]
  plot_data <- d[!is.na(d$cahn_group),]
  if (nrow(plot_data) == 0 || all(is.na(plot_data$cell_mean))) return(NULL) # Return NULL if no data
  ggplot(plot_data, aes(x=factor(cahn_group), y=cell_mean, fill=factor(cahn_group))) +
    geom_boxplot(outlier.size=0.5) + ## FIX: `outlier.linewidth` is deprecated, using `outlier.size`
    labs(title=paste("Cahn Group in", ct), x="Cahn Group", y="Mean Expression") +
    theme_bw() +
    theme(legend.position = "none", axis.text.x = element_text(angle=45, hjust=1))
})
save_grid_pdf_with_stats(plots_cahn_celltype_jl_10_percent_filtered, list(), file.path(output_dir, "grid_mean_expr_cahn_by_celltype_jl_10_percent_filtered.pdf"), "Mean Expression by Cahn Group in Celltypes (James Lloyd 10% Filtered)")

plots_cahn_celltype_jl_20_percent_filtered <- lapply(celltype_levels, function(ct) {
  d <- plot_results_jl_20_percent_filtered
  d$cell_mean <- gene_means_celltype[d$gene, ct]
  plot_data <- d[!is.na(d$cahn_group),]
  if (nrow(plot_data) == 0 || all(is.na(plot_data$cell_mean))) return(NULL) # Return NULL if no data
  ggplot(plot_data, aes(x=factor(cahn_group), y=cell_mean, fill=factor(cahn_group))) +
    geom_boxplot(outlier.size=0.5) + ## FIX: `outlier.linewidth` is deprecated, using `outlier.size`
    labs(title=paste("Cahn Group in", ct), x="Cahn Group", y="Mean Expression") +
    theme_bw() +
    theme(legend.position = "none", axis.text.x = element_text(angle=45, hjust=1))
})
save_grid_pdf_with_stats(plots_cahn_celltype_jl_20_percent_filtered, list(), file.path(output_dir, "grid_mean_expr_cahn_by_celltype_jl_20_percent_filtered.pdf"), "Mean Expression by Cahn Group in Celltypes (James Lloyd 20% Filtered)")

plots_cahn_celltype_mean_gt_0.5_filtered <- lapply(celltype_levels, function(ct) {
  d <- plot_results_mean_gt_0.5_filtered
  d$cell_mean <- gene_means_celltype[d$gene, ct]
  plot_data <- d[!is.na(d$cahn_group),]
  if (nrow(plot_data) == 0 || all(is.na(plot_data$cell_mean))) return(NULL) # Return NULL if no data
  ggplot(plot_data, aes(x=factor(cahn_group), y=cell_mean, fill=factor(cahn_group))) +
    geom_boxplot(outlier.size=0.5) + ## FIX: `outlier.linewidth` is deprecated, using `outlier.size`
    labs(title=paste("Cahn Group in", ct), x="Cahn Group", y="Mean Expression") +
    theme_bw() +
    theme(legend.position = "none", axis.text.x = element_text(angle=45, hjust=1))
})
save_grid_pdf_with_stats(plots_cahn_celltype_mean_gt_0.5_filtered, list(), file.path(output_dir, "grid_mean_expr_cahn_by_celltype_mean_gt_0.5_filtered.pdf"), "Mean Expression by Cahn Group in Celltypes (Mean > 0.5 Filtered)")

plots_cahn_celltype_mean_gt_1_filtered <- lapply(celltype_levels, function(ct) {
  d <- plot_results_mean_gt_1_filtered
  d$cell_mean <- gene_means_celltype[d$gene, ct]
  plot_data <- d[!is.na(d$cahn_group),]
  if (nrow(plot_data) == 0 || all(is.na(plot_data$cell_mean))) return(NULL) # Return NULL if no data
  ggplot(plot_data, aes(x=factor(cahn_group), y=cell_mean, fill=factor(cahn_group))) +
    geom_boxplot(outlier.size=0.5) + ## FIX: `outlier.linewidth` is deprecated, using `outlier.size`
    labs(title=paste("Cahn Group in", ct), x="Cahn Group", y="Mean Expression") +
    theme_bw() +
    theme(legend.position = "none", axis.text.x = element_text(angle=45, hjust=1))
})
save_grid_pdf_with_stats(plots_cahn_celltype_mean_gt_1_filtered, list(), file.path(output_dir, "grid_mean_expr_cahn_by_celltype_mean_gt_1_filtered.pdf"), "Mean Expression by Cahn Group in Celltypes (Mean > 1 Filtered)")


# --- Boxplots by celltype (Bewick group) - Mean Expression ---
plots_bewick_celltype_unfiltered <- lapply(celltype_levels, function(ct) {
  d <- plot_results
  d$cell_mean <- gene_means_celltype[d$gene, ct]
  plot_data <- d[!is.na(d$bewick_group),]
  if (nrow(plot_data) == 0 || all(is.na(plot_data$cell_mean))) return(NULL) # Return NULL if no data
  ggplot(plot_data, aes(x=factor(bewick_group), y=cell_mean, fill=factor(bewick_group))) +
    geom_boxplot(outlier.size=0.5) + ## FIX: `outlier.linewidth` is deprecated, using `outlier.size`
    scale_x_discrete(drop = FALSE) + # Added for Bewick plots
    labs(title=paste("Bewick Group in", ct), x="Bewick Group", y="Mean Expression") +
    theme_bw() +
    theme(legend.position = "none", axis.text.x = element_text(angle=45, hjust=1))
})
save_grid_pdf_with_stats(plots_bewick_celltype_unfiltered, list(), file.path(output_dir, "grid_mean_expr_bewick_by_celltype_unfiltered.pdf"), "Mean Expression by Bewick Group in Celltypes (Unfiltered)")

plots_bewick_celltype_jl_10_percent_filtered <- lapply(celltype_levels, function(ct) {
  d <- plot_results_jl_10_percent_filtered
  d$cell_mean <- gene_means_celltype[d$gene, ct]
  plot_data <- d[!is.na(d$bewick_group),]
  if (nrow(plot_data) == 0 || all(is.na(plot_data$cell_mean))) return(NULL) # Return NULL if no data
  ggplot(plot_data, aes(x=factor(bewick_group), y=cell_mean, fill=factor(bewick_group))) +
    geom_boxplot(outlier.size=0.5) + ## FIX: `outlier.linewidth` is deprecated, using `outlier.size`
    scale_x_discrete(drop = FALSE) + # Added for Bewick plots
    labs(title=paste("Bewick Group in", ct), x="Bewick Group", y="Mean Expression") +
    theme_bw() +
    theme(legend.position = "none", axis.text.x = element_text(angle=45, hjust=1))
})
save_grid_pdf_with_stats(plots_bewick_celltype_jl_10_percent_filtered, list(), file.path(output_dir, "grid_mean_expr_bewick_by_celltype_jl_10_percent_filtered.pdf"), "Mean Expression by Bewick Group in Celltypes (James Lloyd 10% Filtered)")

plots_bewick_celltype_jl_20_percent_filtered <- lapply(celltype_levels, function(ct) {
  d <- plot_results_jl_20_percent_filtered
  d$cell_mean <- gene_means_celltype[d$gene, ct]
  plot_data <- d[!is.na(d$bewick_group),]
  if (nrow(plot_data) == 0 || all(is.na(plot_data$cell_mean))) return(NULL) # Return NULL if no data
  ggplot(plot_data, aes(x=factor(bewick_group), y=cell_mean, fill=factor(bewick_group))) +
    geom_boxplot(outlier.size=0.5) + ## FIX: `outlier.linewidth` is deprecated, using `outlier.size`
    scale_x_discrete(drop = FALSE) + # Added for Bewick plots
    labs(title=paste("Bewick Group in", ct), x="Bewick Group", y="Mean Expression") +
    theme_bw() +
    theme(legend.position = "none", axis.text.x = element_text(angle=45, hjust=1))
})
save_grid_pdf_with_stats(plots_bewick_celltype_jl_20_percent_filtered, list(), file.path(output_dir, "grid_mean_expr_bewick_by_celltype_jl_20_percent_filtered.pdf"), "Mean Expression by Bewick Group in Celltypes (James Lloyd 20% Filtered)")

plots_bewick_celltype_mean_gt_0.5_filtered <- lapply(celltype_levels, function(ct) {
  d <- plot_results_mean_gt_0.5_filtered
  d$cell_mean <- gene_means_celltype[d$gene, ct]
  plot_data <- d[!is.na(d$bewick_group),]
  if (nrow(plot_data) == 0 || all(is.na(plot_data$cell_mean))) return(NULL) # Return NULL if no data
  ggplot(plot_data, aes(x=factor(bewick_group), y=cell_mean, fill=factor(bewick_group))) +
    geom_boxplot(outlier.size=0.5) + ## FIX: `outlier.linewidth` is deprecated, using `outlier.size`
    scale_x_discrete(drop = FALSE) + # Added for Bewick plots
    labs(title=paste("Bewick Group in", ct), x="Bewick Group", y="Mean Expression") +
    theme_bw() +
    theme(legend.position = "none", axis.text.x = element_text(angle=45, hjust=1))
})
save_grid_pdf_with_stats(plots_bewick_celltype_mean_gt_0.5_filtered, list(), file.path(output_dir, "grid_mean_expr_bewick_by_celltype_mean_gt_0.5_filtered.pdf"), "Mean Expression by Bewick Group in Celltypes (Mean > 0.5 Filtered)")

plots_bewick_celltype_mean_gt_1_filtered <- lapply(celltype_levels, function(ct) {
  d <- plot_results_mean_gt_1_filtered
  d$cell_mean <- gene_means_celltype[d$gene, ct]
  plot_data <- d[!is.na(d$bewick_group),]
  if (nrow(plot_data) == 0 || all(is.na(plot_data$cell_mean))) return(NULL) # Return NULL if no data
  ggplot(plot_data, aes(x=factor(bewick_group), y=cell_mean, fill=factor(bewick_group))) +
    geom_boxplot(outlier.size=0.5) + ## FIX: `outlier.linewidth` is deprecated, using `outlier.size`
    scale_x_discrete(drop = FALSE) + # Added for Bewick plots
    labs(title=paste("Bewick Group in", ct), x="Bewick Group", y="Mean Expression") +
    theme_bw() +
    theme(legend.position = "none", axis.text.x = element_text(angle=45, hjust=1))
})
save_grid_pdf_with_stats(plots_bewick_celltype_mean_gt_1_filtered, list(), file.path(output_dir, "grid_mean_expr_bewick_by_celltype_mean_gt_1_filtered.pdf"), "Mean Expression by Bewick Group in Celltypes (Mean > 1 Filtered)")


# --- Boxplots by lineage (Cahn group) - Mean Expression ---
plots_cahn_lineage_unfiltered <- lapply(lineage_levels, function(ln) {
  d <- plot_results
  d$lineage_mean <- gene_means_lineage[d$gene, ln]
  plot_data <- d[!is.na(d$cahn_group),]
  if (nrow(plot_data) == 0 || all(is.na(plot_data$lineage_mean))) return(NULL)
  ggplot(plot_data, aes(x=factor(cahn_group), y=lineage_mean, fill=factor(cahn_group))) +
    geom_boxplot(outlier.size=0.5) + ## FIX: `outlier.linewidth` is deprecated, using `outlier.size`
    labs(title=paste("Cahn Group in", ln), x="Cahn Group", y="Mean Expression") +
    theme_bw() +
    theme(legend.position = "none", axis.text.x = element_text(angle=45, hjust=1))
})
save_grid_pdf_with_stats(plots_cahn_lineage_unfiltered, list(), file.path(output_dir, "grid_mean_expr_cahn_by_lineage_unfiltered.pdf"), "Mean Expression by Cahn Group in Lineages (Unfiltered)")

plots_cahn_lineage_jl_10_percent_filtered <- lapply(lineage_levels, function(ln) {
  d <- plot_results_jl_10_percent_filtered
  d$lineage_mean <- gene_means_lineage[d$gene, ln]
  plot_data <- d[!is.na(d$cahn_group),]
  if (nrow(plot_data) == 0 || all(is.na(plot_data$lineage_mean))) return(NULL)
  ggplot(plot_data, aes(x=factor(cahn_group), y=lineage_mean, fill=factor(cahn_group))) +
    geom_boxplot(outlier.size=0.5) + ## FIX: `outlier.linewidth` is deprecated, using `outlier.size`
    labs(title=paste("Cahn Group in", ln), x="Cahn Group", y="Mean Expression") +
    theme_bw() +
    theme(legend.position = "none", axis.text.x = element_text(angle=45, hjust=1))
})
save_grid_pdf_with_stats(plots_cahn_lineage_jl_10_percent_filtered, list(), file.path(output_dir, "grid_mean_expr_cahn_by_lineage_jl_10_percent_filtered.pdf"), "Mean Expression by Cahn Group in Lineages (James Lloyd 10% Filtered)")

plots_cahn_lineage_jl_20_percent_filtered <- lapply(lineage_levels, function(ln) {
  d <- plot_results_jl_20_percent_filtered
  d$lineage_mean <- gene_means_lineage[d$gene, ln]
  plot_data <- d[!is.na(d$cahn_group),]
  if (nrow(plot_data) == 0 || all(is.na(plot_data$lineage_mean))) return(NULL)
  ggplot(plot_data, aes(x=factor(cahn_group), y=lineage_mean, fill=factor(cahn_group))) +
    geom_boxplot(outlier.size=0.5) + ## FIX: `outlier.linewidth` is deprecated, using `outlier.size`
    labs(title=paste("Cahn Group in", ln), x="Cahn Group", y="Mean Expression") +
    theme_bw() +
    theme(legend.position = "none", axis.text.x = element_text(angle=45, hjust=1))
})
save_grid_pdf_with_stats(plots_cahn_lineage_jl_20_percent_filtered, list(), file.path(output_dir, "grid_mean_expr_cahn_by_lineage_jl_20_percent_filtered.pdf"), "Mean Expression by Cahn Group in Lineages (James Lloyd 20% Filtered)")

plots_cahn_lineage_mean_gt_0.5_filtered <- lapply(lineage_levels, function(ln) {
  d <- plot_results_mean_gt_0.5_filtered
  d$lineage_mean <- gene_means_lineage[d$gene, ln]
  plot_data <- d[!is.na(d$cahn_group),]
  if (nrow(plot_data) == 0 || all(is.na(plot_data$lineage_mean))) return(NULL)
  ggplot(plot_data, aes(x=factor(cahn_group), y=lineage_mean, fill=factor(cahn_group))) +
    geom_boxplot(outlier.size=0.5) + ## FIX: `outlier.linewidth` is deprecated, using `outlier.size`
    labs(title=paste("Cahn Group in", ln), x="Cahn Group", y="Mean Expression") +
    theme_bw() +
    theme(legend.position = "none", axis.text.x = element_text(angle=45, hjust=1))
})
save_grid_pdf_with_stats(plots_cahn_lineage_mean_gt_0.5_filtered, list(), file.path(output_dir, "grid_mean_expr_cahn_by_lineage_mean_gt_0.5_filtered.pdf"), "Mean Expression by Cahn Group in Lineages (Mean > 0.5 Filtered)")

plots_cahn_lineage_mean_gt_1_filtered <- lapply(lineage_levels, function(ln) {
  d <- plot_results_mean_gt_1_filtered
  d$lineage_mean <- gene_means_lineage[d$gene, ln]
  plot_data <- d[!is.na(d$cahn_group),]
  if (nrow(plot_data) == 0 || all(is.na(plot_data$lineage_mean))) return(NULL)
  ggplot(plot_data, aes(x=factor(cahn_group), y=lineage_mean, fill=factor(cahn_group))) +
    geom_boxplot(outlier.size=0.5) + ## FIX: `outlier.linewidth` is deprecated, using `outlier.size`
    labs(title=paste("Cahn Group in", ln), x="Cahn Group", y="Mean Expression") +
    theme_bw() +
    theme(legend.position = "none", axis.text.x = element_text(angle=45, hjust=1))
})
save_grid_pdf_with_stats(plots_cahn_lineage_mean_gt_1_filtered, list(), file.path(output_dir, "grid_mean_expr_cahn_by_lineage_mean_gt_1_filtered.pdf"), "Mean Expression by Cahn Group in Lineages (Mean > 1 Filtered)")


# --- Boxplots by lineage (Bewick group) - Mean Expression ---
plots_bewick_lineage_unfiltered <- lapply(lineage_levels, function(ln) {
  d <- plot_results
  d$lineage_mean <- gene_means_lineage[d$gene, ln]
  plot_data <- d[!is.na(d$bewick_group),]
  if (nrow(plot_data) == 0 || all(is.na(plot_data$lineage_mean))) return(NULL)
  ggplot(plot_data, aes(x=factor(bewick_group), y=lineage_mean, fill=factor(bewick_group))) +
    geom_boxplot(outlier.size=0.5) + ## FIX: `outlier.linewidth` is deprecated, using `outlier.size`
    scale_x_discrete(drop = FALSE) + # Added for Bewick plots
    labs(title=paste("Bewick Group in", ln), x="Bewick Group", y="Mean Expression") +
    theme_bw() +
    theme(legend.position = "none", axis.text.x = element_text(angle=45, hjust=1))
})
save_grid_pdf_with_stats(plots_bewick_lineage_unfiltered, list(), file.path(output_dir, "grid_mean_expr_bewick_by_lineage_unfiltered.pdf"), "Mean Expression by Bewick Group in Lineages (Unfiltered)")

plots_bewick_lineage_jl_10_percent_filtered <- lapply(lineage_levels, function(ln) {
  d <- plot_results_jl_10_percent_filtered
  d$lineage_mean <- gene_means_lineage[d$gene, ln]
  plot_data <- d[!is.na(d$bewick_group),]
  if (nrow(plot_data) == 0 || all(is.na(plot_data$lineage_mean))) return(NULL)
  ggplot(plot_data, aes(x=factor(bewick_group), y=lineage_mean, fill=factor(bewick_group))) +
    geom_boxplot(outlier.size=0.5) + ## FIX: `outlier.linewidth` is deprecated, using `outlier.size`
    scale_x_discrete(drop = FALSE) + # Added for Bewick plots
    labs(title=paste("Bewick Group in", ln), x="Bewick Group", y="Mean Expression") +
    theme_bw() +
    theme(legend.position = "none", axis.text.x = element_text(angle=45, hjust=1))
})
save_grid_pdf_with_stats(plots_bewick_lineage_jl_10_percent_filtered, list(), file.path(output_dir, "grid_mean_expr_bewick_by_lineage_jl_10_percent_filtered.pdf"), "Mean Expression by Bewick Group in Lineages (James Lloyd 10% Filtered)")

plots_bewick_lineage_jl_20_percent_filtered <- lapply(lineage_levels, function(ln) {
  d <- plot_results_jl_20_percent_filtered
  d$lineage_mean <- gene_means_lineage[d$gene, ln]
  plot_data <- d[!is.na(d$bewick_group),]
  if (nrow(plot_data) == 0 || all(is.na(plot_data$lineage_mean))) return(NULL)
  ggplot(plot_data, aes(x=factor(bewick_group), y=lineage_mean, fill=factor(bewick_group))) +
    geom_boxplot(outlier.size=0.5) + ## FIX: `outlier.linewidth` is deprecated, using `outlier.size`
    scale_x_discrete(drop = FALSE) + # Added for Bewick plots
    labs(title=paste("Bewick Group in", ln), x="Bewick Group", y="Mean Expression") +
    theme_bw() +
    theme(legend.position = "none", axis.text.x = element_text(angle=45, hjust=1))
})
save_grid_pdf_with_stats(plots_bewick_lineage_jl_20_percent_filtered, list(), file.path(output_dir, "grid_mean_expr_bewick_by_lineage_jl_20_percent_filtered.pdf"), "Mean Expression by Bewick Group in Lineages (James Lloyd 20% Filtered)")

plots_bewick_lineage_mean_gt_0.5_filtered <- lapply(lineage_levels, function(ln) {
  d <- plot_results_mean_gt_0.5_filtered
  d$lineage_mean <- gene_means_lineage[d$gene, ln]
  plot_data <- d[!is.na(d$bewick_group),]
  if (nrow(plot_data) == 0 || all(is.na(plot_data$lineage_mean))) return(NULL)
  ggplot(plot_data, aes(x=factor(bewick_group), y=lineage_mean, fill=factor(bewick_group))) +
    geom_boxplot(outlier.size=0.5) + ## FIX: `outlier.linewidth` is deprecated, using `outlier.size`
    scale_x_discrete(drop = FALSE) + # Added for Bewick plots
    labs(title=paste("Bewick Group in", ln), x="Bewick Group", y="Mean Expression") +
    theme_bw() +
    theme(legend.position = "none", axis.text.x = element_text(angle=45, hjust=1))
})
save_grid_pdf_with_stats(plots_bewick_lineage_mean_gt_0.5_filtered, list(), file.path(output_dir, "grid_mean_expr_bewick_by_lineage_mean_gt_0.5_filtered.pdf"), "Mean Expression by Bewick Group in Lineages (Mean > 0.5 Filtered)")

plots_bewick_lineage_mean_gt_1_filtered <- lapply(lineage_levels, function(ln) {
  d <- plot_results_mean_gt_1_filtered
  d$lineage_mean <- gene_means_lineage[d$gene, ln]
  plot_data <- d[!is.na(d$bewick_group),]
  if (nrow(plot_data) == 0 || all(is.na(plot_data$lineage_mean))) return(NULL)
  ggplot(plot_data, aes(x=factor(bewick_group), y=lineage_mean, fill=factor(bewick_group))) +
    geom_boxplot(outlier.size=0.5) + ## FIX: `outlier.linewidth` is deprecated, using `outlier.size`
    scale_x_discrete(drop = FALSE) + # Added for Bewick plots
    labs(title=paste("Bewick Group in", ln), x="Bewick Group", y="Mean Expression") +
    theme_bw() +
    theme(legend.position = "none", axis.text.x = element_text(angle=45, hjust=1))
})
save_grid_pdf_with_stats(plots_bewick_lineage_mean_gt_1_filtered, list(), file.path(output_dir, "grid_mean_expr_bewick_by_lineage_mean_gt_1_filtered.pdf"), "Mean Expression by Bewick Group in Lineages (Mean > 1 Filtered)")


# --- Boxplots by celltype (Cahn group) - Expression Noise (CV) ---
plots_cahn_celltype_cv_unfiltered <- list()
stats_cahn_celltype_cv_unfiltered <- list()
for (ct in celltype_levels) {
  d_ct <- gene_data_by_celltype[[ct]]
  if (is.null(d_ct) || nrow(d_ct) == 0) next
  
  plot_data <- d_ct[!is.na(d_ct$cahn_group),]
  if (nrow(plot_data) == 0 || all(is.na(plot_data$cv_expr_ct))) {
    plots_cahn_celltype_cv_unfiltered[[ct]] <- NULL
  } else {
    plot_title <- paste0("Cahn Group in ", ct)
    plots_cahn_celltype_cv_unfiltered[[ct]] <- ggplot(plot_data, aes(x=factor(cahn_group), y=cv_expr_ct, fill=factor(cahn_group))) +
      geom_boxplot(outlier.size=0.5) + ## FIX: `outlier.linewidth` is deprecated, using `outlier.size`
      labs(title=plot_title, x="Cahn Group", y="Expression Noise (CV)") +
      theme_bw() +
      theme(legend.position = "none", axis.text.x = element_text(angle=45, hjust=1))
  }
  stats_cahn_celltype_cv_unfiltered[[ct]] <- run_specific_group_stats(d_ct, "cv_expr_ct", "cahn_group", "CV Expr by Celltype", "Unfiltered", ct)
}
save_grid_pdf_with_stats(plots_cahn_celltype_cv_unfiltered, stats_cahn_celltype_cv_unfiltered, file.path(output_dir, "grid_cv_expr_cahn_by_celltype_unfiltered.pdf"), "Expression Noise (CV) by Cahn Group in Celltypes (Unfiltered)")


plots_cahn_celltype_cv_jl_10_percent_filtered <- list()
stats_cahn_celltype_cv_jl_10_percent_filtered <- list()
for (ct in celltype_levels) {
  d_ct <- gene_data_by_celltype[[ct]]
  if (is.null(d_ct) || nrow(d_ct) == 0) next

  d_filtered_jl_ct <- d_ct %>% filter(gene %in% plot_results_jl_10_percent_filtered$gene)
  plot_data <- d_filtered_jl_ct[!is.na(d_filtered_jl_ct$cahn_group),]
  if (nrow(plot_data) == 0 || all(is.na(plot_data$cv_expr_ct))) {
    plots_cahn_celltype_cv_jl_10_percent_filtered[[ct]] <- NULL
  } else {
    plot_title <- paste0("Cahn Group in ", ct)
    plots_cahn_celltype_cv_jl_10_percent_filtered[[ct]] <- ggplot(plot_data, aes(x=factor(cahn_group), y=cv_expr_ct, fill=factor(cahn_group))) +
      geom_boxplot(outlier.size=0.5) + ## FIX: `outlier.linewidth` is deprecated, using `outlier.size`
      labs(title=plot_title, x="Cahn Group", y="Expression Noise (CV) (James Lloyd 10% Filtered)") +
      theme_bw() +
      theme(legend.position = "none", axis.text.x = element_text(angle=45, hjust=1))
  }
  stats_cahn_celltype_cv_jl_10_percent_filtered[[ct]] <- run_specific_group_stats(d_filtered_jl_ct, "cv_expr_ct", "cahn_group", "CV Expr by Celltype", "James Lloyd 10% Filtered", ct)
}
save_grid_pdf_with_stats(plots_cahn_celltype_cv_jl_10_percent_filtered, stats_cahn_celltype_cv_jl_10_percent_filtered, file.path(output_dir, "grid_cv_expr_cahn_by_celltype_james_lloyd_10_percent_filtered.pdf"), "Expression Noise (CV) by Cahn Group in Celltypes (James Lloyd 10% Filtered)")

plots_cahn_celltype_cv_jl_20_percent_filtered <- list()
stats_cahn_celltype_cv_jl_20_percent_filtered <- list()
for (ct in celltype_levels) {
  d_ct <- gene_data_by_celltype[[ct]]
  if (is.null(d_ct) || nrow(d_ct) == 0) next
  d_filtered_jl_ct <- d_ct %>% filter(gene %in% plot_results_jl_20_percent_filtered$gene)
  plot_data <- d_filtered_jl_ct[!is.na(d_filtered_jl_ct$cahn_group),]
  if (nrow(plot_data) == 0 || all(is.na(plot_data$cv_expr_ct))) {
    plots_cahn_celltype_cv_jl_20_percent_filtered[[ct]] <- NULL
  } else {
    plot_title <- paste0("Cahn Group in ", ct)
    plots_cahn_celltype_cv_jl_20_percent_filtered[[ct]] <- ggplot(plot_data, aes(x=factor(cahn_group), y=cv_expr_ct, fill=factor(cahn_group))) +
      geom_boxplot(outlier.size=0.5) + ## FIX: `outlier.linewidth` is deprecated, using `outlier.size`
      labs(title=plot_title, x="Cahn Group", y="Expression Noise (CV) (James Lloyd 20% Filtered)") +
      theme_bw() +
      theme(legend.position = "none", axis.text.x = element_text(angle=45, hjust=1))
  }
  stats_cahn_celltype_cv_jl_20_percent_filtered[[ct]] <- run_specific_group_stats(d_filtered_jl_ct, "cv_expr_ct", "cahn_group", "CV Expr by Celltype", "James Lloyd 20% Filtered", ct)
}
save_grid_pdf_with_stats(plots_cahn_celltype_cv_jl_20_percent_filtered, stats_cahn_celltype_cv_jl_20_percent_filtered, file.path(output_dir, "grid_cv_expr_cahn_by_celltype_james_lloyd_20_percent_filtered.pdf"), "Expression Noise (CV) by Cahn Group in Celltypes (James Lloyd 20% Filtered)")

plots_cahn_celltype_cv_mean_gt_0.5_filtered <- list()
stats_cahn_celltype_cv_mean_gt_0.5_filtered <- list()
for (ct in celltype_levels) {
  d_ct <- gene_data_by_celltype[[ct]]
  if (is.null(d_ct) || nrow(d_ct) == 0) next
  d_filtered_mean_gt_0.5_ct <- d_ct %>% filter(gene %in% plot_results_mean_gt_0.5_filtered$gene)
  plot_data <- d_filtered_mean_gt_0.5_ct[!is.na(d_filtered_mean_gt_0.5_ct$cahn_group),]
  if (nrow(plot_data) == 0 || all(is.na(plot_data$cv_expr_ct))) {
    plots_cahn_celltype_cv_mean_gt_0.5_filtered[[ct]] <- NULL
  } else {
    plot_title <- paste0("Cahn Group in ", ct)
    plots_cahn_celltype_cv_mean_gt_0.5_filtered[[ct]] <- ggplot(plot_data, aes(x=factor(cahn_group), y=cv_expr_ct, fill=factor(cahn_group))) +
      geom_boxplot(outlier.size=0.5) + ## FIX: `outlier.linewidth` is deprecated, using `outlier.size`
      labs(title=plot_title, x="Cahn Group", y="Expression Noise (CV) (Mean > 0.5 Filtered)") +
      theme_bw() +
      theme(legend.position = "none", axis.text.x = element_text(angle=45, hjust=1))
  }
  stats_cahn_celltype_cv_mean_gt_0.5_filtered[[ct]] <- run_specific_group_stats(d_filtered_mean_gt_0.5_ct, "cv_expr_ct", "cahn_group", "CV Expr by Celltype", "Mean > 0.5 Filtered", ct)
}
save_grid_pdf_with_stats(plots_cahn_celltype_cv_mean_gt_0.5_filtered, stats_cahn_celltype_cv_mean_gt_0.5_filtered, file.path(output_dir, "grid_cv_expr_cahn_by_celltype_mean_gt_0.5_filtered.pdf"), "Expression Noise (CV) by Cahn Group in Celltypes (Mean > 0.5 Filtered)")

plots_cahn_celltype_cv_mean_gt_1_filtered <- list()
stats_cahn_celltype_cv_mean_gt_1_filtered <- list()
for (ct in celltype_levels) {
  d_ct <- gene_data_by_celltype[[ct]]
  if (is.null(d_ct) || nrow(d_ct) == 0) next
  d_filtered_mean_gt_1_ct <- d_ct %>% filter(gene %in% plot_results_mean_gt_1_filtered$gene)
  plot_data <- d_filtered_mean_gt_1_ct[!is.na(d_filtered_mean_gt_1_ct$cahn_group),]
  if (nrow(plot_data) == 0 || all(is.na(plot_data$cv_expr_ct))) {
    plots_cahn_celltype_cv_mean_gt_1_filtered[[ct]] <- NULL
  } else {
    plot_title <- paste0("Cahn Group in ", ct)
    plots_cahn_celltype_cv_mean_gt_1_filtered[[ct]] <- ggplot(plot_data, aes(x=factor(cahn_group), y=cv_expr_ct, fill=factor(cahn_group))) +
      geom_boxplot(outlier.size=0.5) + ## FIX: `outlier.linewidth` is deprecated, using `outlier.size`
      labs(title=plot_title, x="Cahn Group", y="Expression Noise (CV) (Mean > 1 Filtered)") +
      theme_bw() +
      theme(legend.position = "none", axis.text.x = element_text(angle=45, hjust=1))
  }
  stats_cahn_celltype_cv_mean_gt_1_filtered[[ct]] <- run_specific_group_stats(d_filtered_mean_gt_1_ct, "cv_expr_ct", "cahn_group", "CV Expr by Celltype", "Mean > 1 Filtered", ct)
}
save_grid_pdf_with_stats(plots_cahn_celltype_cv_mean_gt_1_filtered, stats_cahn_celltype_cv_mean_gt_1_filtered, file.path(output_dir, "grid_cv_expr_cahn_by_celltype_mean_gt_1_filtered.pdf"), "Expression Noise (CV) by Cahn Group in Celltypes (Mean > 1 Filtered)")


# --- Boxplots by celltype (Bewick group) - Expression Noise (CV) ---
plots_bewick_celltype_cv_unfiltered <- list()
stats_bewick_celltype_cv_unfiltered <- list()
for (ct in celltype_levels) {
  d_ct <- gene_data_by_celltype[[ct]]
  if (is.null(d_ct) || nrow(d_ct) == 0) next
  
  plot_data <- d_ct[!is.na(d_ct$bewick_group),]
  if (nrow(plot_data) == 0 || all(is.na(plot_data$cv_expr_ct))) {
    plots_bewick_celltype_cv_unfiltered[[ct]] <- NULL
  } else {
    plot_title <- paste0("Bewick Group in ", ct)
    plots_bewick_celltype_cv_unfiltered[[ct]] <- ggplot(plot_data, aes(x=factor(bewick_group), y=cv_expr_ct, fill=factor(bewick_group))) +
      geom_boxplot(outlier.size=0.5) + ## FIX: `outlier.linewidth` is deprecated, using `outlier.size`
      scale_x_discrete(drop = FALSE) + # Added for Bewick plots
      labs(title=plot_title, x="Bewick Group", y="Expression Noise (CV)") +
      theme_bw() +
      theme(legend.position = "none", axis.text.x = element_text(angle=45, hjust=1))
  }
  stats_bewick_celltype_cv_unfiltered[[ct]] <- run_specific_group_stats(d_ct, "cv_expr_ct", "bewick_group", "CV Expr by Celltype", "Unfiltered", ct)
}
save_grid_pdf_with_stats(plots_bewick_celltype_cv_unfiltered, stats_bewick_celltype_cv_unfiltered, file.path(output_dir, "grid_cv_expr_bewick_by_celltype_unfiltered.pdf"), "Expression Noise (CV) by Bewick Group in Celltypes (Unfiltered)")

plots_bewick_celltype_cv_jl_10_percent_filtered <- list()
stats_bewick_celltype_cv_jl_10_percent_filtered <- list()
for (ct in celltype_levels) {
  d_ct <- gene_data_by_celltype[[ct]]
  if (is.null(d_ct) || nrow(d_ct) == 0) next
  d_filtered_jl_ct <- d_ct %>% filter(gene %in% plot_results_jl_10_percent_filtered$gene)
  plot_data <- d_filtered_jl_ct[!is.na(d_filtered_jl_ct$bewick_group),]
  if (nrow(plot_data) == 0 || all(is.na(plot_data$cv_expr_ct))) {
    plots_bewick_celltype_cv_jl_10_percent_filtered[[ct]] <- NULL
  } else {
    plot_title <- paste0("Bewick Group in ", ct)
    plots_bewick_celltype_cv_jl_10_percent_filtered[[ct]] <- ggplot(plot_data, aes(x=factor(bewick_group), y=cv_expr_ct, fill=factor(bewick_group))) +
      geom_boxplot(outlier.size=0.5) + ## FIX: `outlier.linewidth` is deprecated, using `outlier.size`
      scale_x_discrete(drop = FALSE) + # Added for Bewick plots
      labs(title=plot_title, x="Bewick Group", y="Expression Noise (CV) (James Lloyd 10% Filtered)") +
      theme_bw() +
      theme(legend.position = "none", axis.text.x = element_text(angle=45, hjust=1))
  }
  stats_bewick_celltype_cv_jl_10_percent_filtered[[ct]] <- run_specific_group_stats(d_filtered_jl_ct, "cv_expr_ct", "bewick_group", "CV Expr by Celltype", "James Lloyd 10% Filtered", ct)
}
save_grid_pdf_with_stats(plots_bewick_celltype_cv_jl_10_percent_filtered, stats_bewick_celltype_cv_jl_10_percent_filtered, file.path(output_dir, "grid_cv_expr_bewick_by_celltype_james_lloyd_10_percent_filtered.pdf"), "Expression Noise (CV) by Bewick Group in Celltypes (James Lloyd 10% Filtered)")

plots_bewick_celltype_cv_jl_20_percent_filtered <- list()
stats_bewick_celltype_cv_jl_20_percent_filtered <- list()
for (ct in celltype_levels) {
  d_ct <- gene_data_by_celltype[[ct]]
  if (is.null(d_ct) || nrow(d_ct) == 0) next
  d_filtered_jl_ct <- d_ct %>% filter(gene %in% plot_results_jl_20_percent_filtered$gene)
  plot_data <- d_filtered_jl_ct[!is.na(d_filtered_jl_ct$bewick_group),]
  if (nrow(plot_data) == 0 || all(is.na(plot_data$cv_expr_ct))) {
    plots_bewick_celltype_cv_jl_20_percent_filtered[[ct]] <- NULL
  } else {
    plot_title <- paste0("Bewick Group in ", ct)
    plots_bewick_celltype_cv_jl_20_percent_filtered[[ct]] <- ggplot(plot_data, aes(x=factor(bewick_group), y=cv_expr_ct, fill=factor(bewick_group))) +
      geom_boxplot(outlier.size=0.5) + ## FIX: `outlier.linewidth` is deprecated, using `outlier.size`
      scale_x_discrete(drop = FALSE) + # Added for Bewick plots
      labs(title=plot_title, x="Bewick Group", y="Expression Noise (CV) (James Lloyd 20% Filtered)") +
      theme_bw() +
      theme(legend.position = "none", axis.text.x = element_text(angle=45, hjust=1))
  }
  stats_bewick_celltype_cv_jl_20_percent_filtered[[ct]] <- run_specific_group_stats(d_filtered_jl_ct, "cv_expr_ct", "bewick_group", "CV Expr by Celltype", "James Lloyd 20% Filtered", ct)
}
save_grid_pdf_with_stats(plots_bewick_celltype_cv_jl_20_percent_filtered, stats_bewick_celltype_cv_jl_20_percent_filtered, file.path(output_dir, "grid_cv_expr_bewick_by_celltype_james_lloyd_20_percent_filtered.pdf"), "Expression Noise (CV) by Bewick Group in Celltypes (James Lloyd 20% Filtered)")

plots_bewick_celltype_cv_mean_gt_0.5_filtered <- list()
stats_bewick_celltype_cv_mean_gt_0.5_filtered <- list()
for (ct in celltype_levels) {
  d_ct <- gene_data_by_celltype[[ct]]
  if (is.null(d_ct) || nrow(d_ct) == 0) next
  d_filtered_mean_gt_0.5_ct <- d_ct %>% filter(gene %in% plot_results_mean_gt_0.5_filtered$gene)
  plot_data <- d_filtered_mean_gt_0.5_ct[!is.na(d_filtered_mean_gt_0.5_ct$bewick_group),]
  if (nrow(plot_data) == 0 || all(is.na(plot_data$cv_expr_ct))) {
    plots_bewick_celltype_cv_mean_gt_0.5_filtered[[ct]] <- NULL
  } else {
    plot_title <- paste0("Bewick Group in ", ct)
    plots_bewick_celltype_cv_mean_gt_0.5_filtered[[ct]] <- ggplot(plot_data, aes(x=factor(bewick_group), y=cv_expr_ct, fill=factor(bewick_group))) +
      geom_boxplot(outlier.size=0.5) + ## FIX: `outlier.linewidth` is deprecated, using `outlier.size`
      scale_x_discrete(drop = FALSE) + # Added for Bewick plots
      labs(title=plot_title, x="Bewick Group", y="Expression Noise (CV) (Mean > 0.5 Filtered)") +
      theme_bw() +
      theme(legend.position = "none", axis.text.x = element_text(angle=45, hjust=1))
  }
  stats_bewick_celltype_cv_mean_gt_0.5_filtered[[ct]] <- run_specific_group_stats(d_filtered_mean_gt_0.5_ct, "cv_expr_ct", "bewick_group", "CV Expr by Celltype", "Mean > 0.5 Filtered", ct)
}
save_grid_pdf_with_stats(plots_bewick_celltype_cv_mean_gt_0.5_filtered, stats_bewick_celltype_cv_mean_gt_0.5_filtered, file.path(output_dir, "grid_cv_expr_bewick_by_celltype_mean_gt_0.5_filtered.pdf"), "Expression Noise (CV) by Bewick Group in Celltypes (Mean > 0.5 Filtered)")

plots_bewick_celltype_cv_mean_gt_1_filtered <- list()
stats_bewick_celltype_cv_mean_gt_1_filtered <- list()
for (ct in celltype_levels) {
  d_ct <- gene_data_by_celltype[[ct]]
  if (is.null(d_ct) || nrow(d_ct) == 0) next
  d_filtered_mean_gt_1_ct <- d_ct %>% filter(gene %in% plot_results_mean_gt_1_filtered$gene)
  plot_data <- d_filtered_mean_gt_1_ct[!is.na(d_filtered_mean_gt_1_ct$bewick_group),]
  if (nrow(plot_data) == 0 || all(is.na(plot_data$cv_expr_ct))) {
    plots_bewick_celltype_cv_mean_gt_1_filtered[[ct]] <- NULL
  } else {
    plot_title <- paste0("Bewick Group in ", ct)
    plots_bewick_celltype_cv_mean_gt_1_filtered[[ct]] <- ggplot(plot_data, aes(x=factor(bewick_group), y=cv_expr_ct, fill=factor(bewick_group))) +
      geom_boxplot(outlier.size=0.5) + ## FIX: `outlier.linewidth` is deprecated, using `outlier.size`
      scale_x_discrete(drop = FALSE) + # Added for Bewick plots
      labs(title=plot_title, x="Bewick Group", y="Expression Noise (CV) (Mean > 1 Filtered)") +
      theme_bw() +
      theme(legend.position = "none", axis.text.x = element_text(angle=45, hjust=1))
  }
  stats_bewick_celltype_cv_mean_gt_1_filtered[[ct]] <- run_specific_group_stats(d_filtered_mean_gt_1_ct, "cv_expr_ct", "bewick_group", "CV Expr by Celltype", "Mean > 1 Filtered", ct)
}
save_grid_pdf_with_stats(plots_bewick_celltype_cv_mean_gt_1_filtered, stats_bewick_celltype_cv_mean_gt_1_filtered, file.path(output_dir, "grid_cv_expr_bewick_by_celltype_mean_gt_1_filtered.pdf"), "Expression Noise (CV) by Bewick Group in Celltypes (Mean > 1 Filtered)")


# --- Boxplots by lineage (Cahn group) - Expression Noise (CV) ---
plots_cahn_lineage_cv_unfiltered <- list()
stats_cahn_lineage_cv_unfiltered <- list()
for (ln in lineage_levels) {
  d_ln <- gene_data_by_lineage[[ln]]
  if (is.null(d_ln) || nrow(d_ln) == 0) next
  
  plot_data <- d_ln[!is.na(d_ln$cahn_group),]
  if (nrow(plot_data) == 0 || all(is.na(plot_data$cv_expr_ln))) {
    plots_cahn_lineage_cv_unfiltered[[ln]] <- NULL
  } else {
    plot_title <- paste0("Cahn Group in ", ln)
    plots_cahn_lineage_cv_unfiltered[[ln]] <- ggplot(plot_data, aes(x=factor(cahn_group), y=cv_expr_ln, fill=factor(cahn_group))) +
      geom_boxplot(outlier.size=0.5) + ## FIX: `outlier.linewidth` is deprecated, using `outlier.size`
      labs(title=plot_title, x="Cahn Group", y="Expression Noise (CV)") +
      theme_bw() +
      theme(legend.position = "none", axis.text.x = element_text(angle=45, hjust=1))
  }
  stats_cahn_lineage_cv_unfiltered[[ln]] <- run_specific_group_stats(d_ln, "cv_expr_ln", "cahn_group", "CV Expr by Lineage", "Unfiltered", ln)
}
save_grid_pdf_with_stats(plots_cahn_lineage_cv_unfiltered, stats_cahn_lineage_cv_unfiltered, file.path(output_dir, "grid_cv_expr_cahn_by_lineage_unfiltered.pdf"), "Expression Noise (CV) by Cahn Group in Lineages (Unfiltered)")

plots_cahn_lineage_cv_jl_10_percent_filtered <- list()
stats_cahn_lineage_cv_jl_10_percent_filtered <- list()
for (ln in lineage_levels) {
  d_ln <- gene_data_by_lineage[[ln]]
  if (is.null(d_ln) || nrow(d_ln) == 0) next
  d_filtered_jl_ln <- d_ln %>% filter(gene %in% plot_results_jl_10_percent_filtered$gene)
  plot_data <- d_filtered_jl_ln[!is.na(d_filtered_jl_ln$cahn_group),]
  if (nrow(plot_data) == 0 || all(is.na(plot_data$cv_expr_ln))) {
    plots_cahn_lineage_cv_jl_10_percent_filtered[[ln]] <- NULL
  } else {
    plot_title <- paste0("Cahn Group in ", ln)
    plots_cahn_lineage_cv_jl_10_percent_filtered[[ln]] <- ggplot(plot_data, aes(x=factor(cahn_group), y=cv_expr_ln, fill=factor(cahn_group))) +
      geom_boxplot(outlier.size=0.5) + ## FIX: `outlier.linewidth` is deprecated, using `outlier.size`
      labs(title=plot_title, x="Cahn Group", y="Expression Noise (CV) (James Lloyd 10% Filtered)") +
      theme_bw() +
      theme(legend.position = "none", axis.text.x = element_text(angle=45, hjust=1))
  }
  stats_cahn_lineage_cv_jl_10_percent_filtered[[ln]] <- run_specific_group_stats(d_filtered_jl_ln, "cv_expr_ln", "cahn_group", "CV Expr by Lineage", "James Lloyd 10% Filtered", ln)
}
save_grid_pdf_with_stats(plots_cahn_lineage_cv_jl_10_percent_filtered, stats_cahn_lineage_cv_jl_10_percent_filtered, file.path(output_dir, "grid_cv_expr_cahn_by_lineage_james_lloyd_10_percent_filtered.pdf"), "Expression Noise (CV) by Cahn Group in Lineages (James Lloyd 10% Filtered)")

plots_cahn_lineage_cv_jl_20_percent_filtered <- list()
stats_cahn_lineage_cv_jl_20_percent_filtered <- list()
for (ln in lineage_levels) {
  d_ln <- gene_data_by_lineage[[ln]]
  if (is.null(d_ln) || nrow(d_ln) == 0) next
  d_filtered_jl_ln <- d_ln %>% filter(gene %in% plot_results_jl_20_percent_filtered$gene)
  plot_data <- d_filtered_jl_ln[!is.na(d_filtered_jl_ln$cahn_group),]
  if (nrow(plot_data) == 0 || all(is.na(plot_data$cv_expr_ln))) {
    plots_cahn_lineage_cv_jl_20_percent_filtered[[ln]] <- NULL
  } else {
    plot_title <- paste0("Cahn Group in ", ln)
    plots_cahn_lineage_cv_jl_20_percent_filtered[[ln]] <- ggplot(plot_data, aes(x=factor(cahn_group), y=cv_expr_ln, fill=factor(cahn_group))) +
      geom_boxplot(outlier.size=0.5) + ## FIX: `outlier.linewidth` is deprecated, using `outlier.size`
      labs(title=plot_title, x="Cahn Group", y="Expression Noise (CV) (James Lloyd 20% Filtered)") +
      theme_bw() +
      theme(legend.position = "none", axis.text.x = element_text(angle=45, hjust=1))
  }
  stats_cahn_lineage_cv_jl_20_percent_filtered[[ln]] <- run_specific_group_stats(d_filtered_jl_ln, "cv_expr_ln", "cahn_group", "CV Expr by Lineage", "James Lloyd 20% Filtered", ln)
}
save_grid_pdf_with_stats(plots_cahn_lineage_cv_jl_20_percent_filtered, stats_cahn_lineage_cv_jl_20_percent_filtered, file.path(output_dir, "grid_cv_expr_cahn_by_lineage_james_lloyd_20_percent_filtered.pdf"), "Expression Noise (CV) by Cahn Group in Lineages (James Lloyd 20% Filtered)")

plots_cahn_lineage_cv_mean_gt_0.5_filtered <- list()
stats_cahn_lineage_cv_mean_gt_0.5_filtered <- list()
for (ln in lineage_levels) {
  d_ln <- gene_data_by_lineage[[ln]]
  if (is.null(d_ln) || nrow(d_ln) == 0) next
  d_filtered_mean_gt_0.5_ln <- d_ln %>% filter(gene %in% plot_results_mean_gt_0.5_filtered$gene)
  plot_data <- d_filtered_mean_gt_0.5_ln[!is.na(d_filtered_mean_gt_0.5_ln$cahn_group),]
  if (nrow(plot_data) == 0 || all(is.na(plot_data$cv_expr_ln))) {
    plots_cahn_lineage_cv_mean_gt_0.5_filtered[[ln]] <- NULL
  } else {
    plot_title <- paste0("Cahn Group in ", ln)
    plots_cahn_lineage_cv_mean_gt_0.5_filtered[[ln]] <- ggplot(plot_data, aes(x=factor(cahn_group), y=cv_expr_ln, fill=factor(cahn_group))) +
      geom_boxplot(outlier.size=0.5) + ## FIX: `outlier.linewidth` is deprecated, using `outlier.size`
      labs(title=plot_title, x="Cahn Group", y="Expression Noise (CV) (Mean > 0.5 Filtered)") +
      theme_bw() +
      theme(legend.position = "none", axis.text.x = element_text(angle=45, hjust=1))
  }
  stats_cahn_lineage_cv_mean_gt_0.5_filtered[[ln]] <- run_specific_group_stats(d_filtered_mean_gt_0.5_ln, "cv_expr_ln", "cahn_group", "CV Expr by Lineage", "Mean > 0.5 Filtered", ln)
}
save_grid_pdf_with_stats(plots_cahn_lineage_cv_mean_gt_0.5_filtered, stats_cahn_lineage_cv_mean_gt_0.5_filtered, file.path(output_dir, "grid_cv_expr_cahn_by_lineage_mean_gt_0.5_filtered.pdf"), "Expression Noise (CV) by Cahn Group in Lineages (Mean > 0.5 Filtered)")

plots_cahn_lineage_cv_mean_gt_1_filtered <- list()
stats_cahn_lineage_cv_mean_gt_1_filtered <- list()
for (ln in lineage_levels) {
  d_ln <- gene_data_by_lineage[[ln]]
  if (is.null(d_ln) || nrow(d_ln) == 0) next
  d_filtered_mean_gt_1_ln <- d_ln %>% filter(gene %in% plot_results_mean_gt_1_filtered$gene)
  plot_data <- d_filtered_mean_gt_1_ln[!is.na(d_filtered_mean_gt_1_ln$cahn_group),]
  if (nrow(plot_data) == 0 || all(is.na(plot_data$cv_expr_ln))) {
    plots_cahn_lineage_cv_mean_gt_1_filtered[[ln]] <- NULL
  } else {
    plot_title <- paste0("Cahn Group in ", ln)
    plots_cahn_lineage_cv_mean_gt_1_filtered[[ln]] <- ggplot(plot_data, aes(x=factor(cahn_group), y=cv_expr_ln, fill=factor(cahn_group))) +
      geom_boxplot(outlier.size=0.5) + ## FIX: `outlier.linewidth` is deprecated, using `outlier.size`
      labs(title=plot_title, x="Cahn Group", y="Expression Noise (CV) (Mean > 1 Filtered)") +
      theme_bw() +
      theme(legend.position = "none", axis.text.x = element_text(angle=45, hjust=1))
  }
  stats_cahn_lineage_cv_mean_gt_1_filtered[[ln]] <- run_specific_group_stats(d_filtered_mean_gt_1_ln, "cv_expr_ln", "cahn_group", "CV Expr by Lineage", "Mean > 1 Filtered", ln)
}
save_grid_pdf_with_stats(plots_cahn_lineage_cv_mean_gt_1_filtered, stats_cahn_lineage_cv_mean_gt_1_filtered, file.path(output_dir, "grid_cv_expr_cahn_by_lineage_mean_gt_1_filtered.pdf"), "Expression Noise (CV) by Cahn Group in Lineages (Mean > 1 Filtered)")


# --- Boxplots by lineage (Bewick group) - Expression Noise (CV) ---
plots_bewick_lineage_cv_unfiltered <- list()
stats_bewick_lineage_cv_unfiltered <- list()
for (ln in lineage_levels) {
  d_ln <- gene_data_by_lineage[[ln]]
  if (is.null(d_ln) || nrow(d_ln) == 0) next
  
  plot_data <- d_ln[!is.na(d_ln$bewick_group),]
  if (nrow(plot_data) == 0 || all(is.na(plot_data$cv_expr_ln))) {
    plots_bewick_lineage_cv_unfiltered[[ln]] <- NULL
  } else {
    plot_title <- paste0("Bewick Group in ", ln)
    plots_bewick_lineage_cv_unfiltered[[ln]] <- ggplot(plot_data, aes(x=factor(bewick_group), y=cv_expr_ln, fill=factor(bewick_group))) +
      geom_boxplot(outlier.size=0.5) + ## FIX: `outlier.linewidth` is deprecated, using `outlier.size`
      scale_x_discrete(drop = FALSE) + # Added for Bewick plots
      labs(title=plot_title, x="Bewick Group", y="Expression Noise (CV)") +
      theme_bw() +
      theme(legend.position = "none", axis.text.x = element_text(angle=45, hjust=1))
  }
  stats_bewick_lineage_cv_unfiltered[[ln]] <- run_specific_group_stats(d_ln, "cv_expr_ln", "bewick_group", "CV Expr by Lineage", "Unfiltered", ln)
}
save_grid_pdf_with_stats(plots_bewick_lineage_cv_unfiltered, stats_bewick_lineage_cv_unfiltered, file.path(output_dir, "grid_cv_expr_bewick_by_lineage_unfiltered.pdf"), "Expression Noise (CV) by Bewick Group in Lineages (Unfiltered)")

plots_bewick_lineage_cv_jl_10_percent_filtered <- list()
stats_bewick_lineage_cv_jl_10_percent_filtered <- list()
for (ln in lineage_levels) {
  d_ln <- gene_data_by_lineage[[ln]]
  if (is.null(d_ln) || nrow(d_ln) == 0) next
  d_filtered_jl_ln <- d_ln %>% filter(gene %in% plot_results_jl_10_percent_filtered$gene)
  plot_data <- d_filtered_jl_ln[!is.na(d_filtered_jl_ln$bewick_group),]
  if (nrow(plot_data) == 0 || all(is.na(plot_data$cv_expr_ln))) {
    plots_bewick_lineage_cv_jl_10_percent_filtered[[ln]] <- NULL
  } else {
    plot_title <- paste0("Bewick Group in ", ln)
    plots_bewick_lineage_cv_jl_10_percent_filtered[[ln]] <- ggplot(plot_data, aes(x=factor(bewick_group), y=cv_expr_ln, fill=factor(bewick_group))) +
      geom_boxplot(outlier.size=0.5) + ## FIX: `outlier.linewidth` is deprecated, using `outlier.size`
      scale_x_discrete(drop = FALSE) + # Added for Bewick plots
      labs(title=plot_title, x="Bewick Group", y="Expression Noise (CV) (James Lloyd 10% Filtered)") +
      theme_bw() +
      theme(legend.position = "none", axis.text.x = element_text(angle=45, hjust=1))
  }
  stats_bewick_lineage_cv_jl_10_percent_filtered[[ln]] <- run_specific_group_stats(d_filtered_jl_ln, "cv_expr_ln", "bewick_group", "CV Expr by Lineage", "James Lloyd 10% Filtered", ln)
}
save_grid_pdf_with_stats(plots_bewick_lineage_cv_jl_10_percent_filtered, stats_bewick_lineage_cv_jl_10_percent_filtered, file.path(output_dir, "grid_cv_expr_bewick_by_lineage_james_lloyd_10_percent_filtered.pdf"), "Expression Noise (CV) by Bewick Group in Lineages (James Lloyd 10% Filtered)")

plots_bewick_lineage_cv_jl_20_percent_filtered <- list()
stats_bewick_lineage_cv_jl_20_percent_filtered <- list()
for (ln in lineage_levels) {
  d_ln <- gene_data_by_lineage[[ln]]
  if (is.null(d_ln) || nrow(d_ln) == 0) next
  d_filtered_jl_ln <- d_ln %>% filter(gene %in% plot_results_jl_20_percent_filtered$gene)
  plot_data <- d_filtered_jl_ln[!is.na(d_filtered_jl_ln$bewick_group),]
  if (nrow(plot_data) == 0 || all(is.na(plot_data$cv_expr_ln))) {
    plots_bewick_lineage_cv_jl_20_percent_filtered[[ln]] <- NULL
  } else {
    plot_title <- paste0("Bewick Group in ", ln)
    plots_bewick_lineage_cv_jl_20_percent_filtered[[ln]] <- ggplot(plot_data, aes(x=factor(bewick_group), y=cv_expr_ln, fill=factor(bewick_group))) +
      geom_boxplot(outlier.size=0.5) + ## FIX: `outlier.linewidth` is deprecated, using `outlier.size`
      scale_x_discrete(drop = FALSE) + # Added for Bewick plots
      labs(title=plot_title, x="Bewick Group", y="Expression Noise (CV) (James Lloyd 20% Filtered)") +
      theme_bw() +
      theme(legend.position = "none", axis.text.x = element_text(angle=45, hjust=1))
  }
  stats_bewick_lineage_cv_jl_20_percent_filtered[[ln]] <- run_specific_group_stats(d_filtered_jl_ln, "cv_expr_ln", "bewick_group", "CV Expr by Lineage", "James Lloyd 20% Filtered", ln)
}
save_grid_pdf_with_stats(plots_bewick_lineage_cv_jl_20_percent_filtered, stats_bewick_lineage_cv_jl_20_percent_filtered, file.path(output_dir, "grid_cv_expr_bewick_by_lineage_james_lloyd_20_percent_filtered.pdf"), "Expression Noise (CV) by Bewick Group in Lineages (James Lloyd 20% Filtered)")

plots_bewick_lineage_cv_mean_gt_0.5_filtered <- list()
stats_bewick_lineage_cv_mean_gt_0.5_filtered <- list()
for (ln in lineage_levels) {
  d_ln <- gene_data_by_lineage[[ln]]
  if (is.null(d_ln) || nrow(d_ln) == 0) next
  d_filtered_mean_gt_0.5_ln <- d_ln %>% filter(gene %in% plot_results_mean_gt_0.5_filtered$gene)
  plot_data <- d_filtered_mean_gt_0.5_ln[!is.na(d_filtered_mean_gt_0.5_ln$bewick_group),]
  if (nrow(plot_data) == 0 || all(is.na(plot_data$cv_expr_ln))) {
    plots_bewick_lineage_cv_mean_gt_0.5_filtered[[ln]] <- NULL
  } else {
    plot_title <- paste0("Bewick Group in ", ln)
    plots_bewick_lineage_cv_mean_gt_0.5_filtered[[ln]] <- ggplot(plot_data, aes(x=factor(bewick_group), y=cv_expr_ln, fill=factor(bewick_group))) +
      geom_boxplot(outlier.size=0.5) + ## FIX: `outlier.linewidth` is deprecated, using `outlier.size`
      scale_x_discrete(drop = FALSE) + # Added for Bewick plots
      labs(title=plot_title, x="Bewick Group", y="Expression Noise (CV) (Mean > 0.5 Filtered)") +
      theme_bw() +
      theme(legend.position = "none", axis.text.x = element_text(angle=45, hjust=1))
  }
  stats_bewick_lineage_cv_mean_gt_0.5_filtered[[ln]] <- run_specific_group_stats(d_filtered_mean_gt_0.5_ln, "cv_expr_ln", "bewick_group", "CV Expr by Lineage", "Mean > 0.5 Filtered", ln)
}
save_grid_pdf_with_stats(plots_bewick_lineage_cv_mean_gt_0.5_filtered, stats_bewick_lineage_cv_mean_gt_0.5_filtered, file.path(output_dir, "grid_cv_expr_bewick_by_lineage_mean_gt_0.5_filtered.pdf"), "Expression Noise (CV) by Bewick Group in Lineages (Mean > 0.5 Filtered)")

plots_bewick_lineage_cv_mean_gt_1_filtered <- list()
stats_bewick_lineage_cv_mean_gt_1_filtered <- list()
for (ln in lineage_levels) {
  d_ln <- gene_data_by_lineage[[ln]]
  if (is.null(d_ln) || nrow(d_ln) == 0) next
  d_filtered_mean_gt_1_ln <- d_ln %>% filter(gene %in% plot_results_mean_gt_1_filtered$gene)
  plot_data <- d_filtered_mean_gt_1_ln[!is.na(d_filtered_mean_gt_1_ln$bewick_group),]
  if (nrow(plot_data) == 0 || all(is.na(plot_data$cv_expr_ln))) {
    plots_bewick_lineage_cv_mean_gt_1_filtered[[ln]] <- NULL
  } else {
    plot_title <- paste0("Bewick Group in ", ln)
    plots_bewick_lineage_cv_mean_gt_1_filtered[[ln]] <- ggplot(plot_data, aes(x=factor(bewick_group), y=cv_expr_ln, fill=factor(bewick_group))) +
      geom_boxplot(outlier.size=0.5) + ## FIX: `outlier.linewidth` is deprecated, using `outlier.size`
      scale_x_discrete(drop = FALSE) + # Added for Bewick plots
      labs(title=plot_title, x="Bewick Group", y="Expression Noise (CV) (Mean > 1 Filtered)") +
      theme_bw() +
      theme(legend.position = "none", axis.text.x = element_text(angle=45, hjust=1))
  }
  stats_bewick_lineage_cv_mean_gt_1_filtered[[ln]] <- run_specific_group_stats(d_filtered_mean_gt_1_ln, "cv_expr_ln", "bewick_group", "CV Expr by Lineage", "Mean > 1 Filtered", ln)
}
save_grid_pdf_with_stats(plots_bewick_lineage_cv_mean_gt_1_filtered, stats_bewick_lineage_cv_mean_gt_1_filtered, file.path(output_dir, "grid_cv_expr_bewick_by_lineage_mean_gt_1_filtered.pdf"), "Expression Noise (CV) by Bewick Group in Lineages (Mean > 1 Filtered)")


# --- Comprehensive Filter Visualization Plot (Main PDF) ---
# Recalculate James Lloyd filter status on plot_results (the baseline data)
plot_results$jl_10_percent_valid <- plot_results$mean_expr >= quantile(results$mean_expr, 0.10, na.rm=TRUE)
plot_results$jl_20_percent_valid <- plot_results$mean_expr >= quantile(results$mean_expr, 0.20, na.rm=TRUE)
plot_results$mean_gt_0.5_valid <- plot_results$mean_expr >= 0.5
plot_results$mean_gt_1_valid <- plot_results$mean_expr >= 1


# Create a new column for categorized plotting based on filter validity
plot_results$filter_category_jl <- "Unfiltered" # Default

# Assign categories based on James Lloyd filter validity (order matters for assignment)
# Assign most stringent filters last so they override less stringent ones
plot_results$filter_category_jl[plot_results$jl_10_percent_valid] <- "JL 10% Filtered"
plot_results$filter_category_jl[plot_results$jl_20_percent_valid] <- "JL 20% Filtered"
plot_results$filter_category_jl[plot_results$mean_gt_0.5_valid] <- "Mean > 0.5 Filtered"
plot_results$filter_category_jl[plot_results$mean_gt_1_valid] <- "Mean > 1 Filtered"


# Convert to factor with desired order for legend and plotting
plot_results$filter_category_jl <- factor(plot_results$filter_category_jl,
                                       levels = c("Unfiltered",
                                                 "JL 10% Filtered",
                                                 "JL 20% Filtered",
                                                 "Mean > 0.5 Filtered",
                                                 "Mean > 1 Filtered"))

# Define colors for these categories
category_colors_jl <- c("Unfiltered" = "grey",
                     "JL 10% Filtered" = "darkgreen",
                     "JL 20% Filtered" = "purple",
                     "Mean > 0.5 Filtered" = "blue",
                     "Mean > 1 Filtered" = "red")

# Get actual cutoff values for subtitle
mean_expr_10th_percentile_val <- quantile(results$mean_expr, 0.10, na.rm=TRUE)
mean_expr_20th_percentile_val <- quantile(results$mean_expr, 0.20, na.rm=TRUE)

p_combined_filter_visualization <- ggplot(plot_results, aes(x = mean_expr, y = cv_expr)) +
  geom_point(aes(color = filter_category_jl), alpha = 0.5, size = 0.7) +
  scale_color_manual(values = category_colors_jl, name = "Gene Filter Status") +
  geom_vline(xintercept = mean_expr_10th_percentile_val, linetype = "dashed", color = "darkgreen", linewidth = 1) +
  geom_vline(xintercept = mean_expr_20th_percentile_val, linetype = "dashed", color = "purple", linewidth = 1) +
  geom_vline(xintercept = 0.5, linetype = "dashed", color = "blue", linewidth = 1) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "red", linewidth = 1) +
  scale_x_log10(labels = scales::scientific) +
  scale_y_log10(labels = scales::scientific) +
  labs(title = "Mean Expression vs CV with James Lloyd Filtering Methods",
       subtitle = paste0("10% (", formatC(mean_expr_10th_percentile_val, format="e", digits=2), "), ",
                         "20% (", formatC(mean_expr_20th_percentile_val, format="e", digits=2), "), ",
                         "Mean > 0.5, Mean > 1 (Dashed Lines)"),
       x = "Mean Expression",
       y = "Coefficient of Variation (CV)") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5),
        legend.position = "bottom",
        legend.box = "horizontal")


# --- Mean Expression vs CV plots for comparison of filtering methods ---
# Collect all plots in a list to save to a single PDF
all_plots_for_pdf <- list()

# Helper function to conditionally add plots to the list
add_plot_if_data <- function(plot_list_name, data_df, group_col, x_var, y_var, title_suffix, log_scale = TRUE) {
  # Ensure data_df has rows and relevant columns have non-NA values
  if (nrow(data_df) == 0 || all(is.na(data_df[[x_var]])) || all(is.na(data_df[[y_var]])) || all(is.na(data_df[[group_col]]))) {
    return(NULL)
  }
  
  ## FIX: `aes_string` is deprecated. Using tidy evaluation with `.data[[]]`
  p <- ggplot(data_df, aes(x = .data[[x_var]], y = .data[[y_var]], color = .data[[group_col]])) +
    geom_point(alpha=0.5, size=0.7) +
    labs(title=paste0("Mean vs CV by ", stringr::str_to_title(gsub("_group", "", group_col)), " Group (", title_suffix, ")"),
         x="Mean Expression", y="Coefficient of Variation (CV)") +
    theme_bw()
  
  if (log_scale) {
    p <- p + scale_x_log10() + scale_y_log10()
  }
  
  # Ensure all factor levels are shown on Bewick plots even if no data
  if (group_col == "bewick_group") {
      p <- p + scale_color_discrete(drop = FALSE)
  }
  return(p)
}


# Add the new combined filter visualization plot as the first page
if (nrow(plot_results) > 0) {
  all_plots_for_pdf[["combined_filter_overview"]] <- p_combined_filter_visualization
}


# Scenario 1: Baseline (only basic cleaning)
p_cahn_baseline <- add_plot_if_data("cahn_baseline", plot_results, "cahn_group", "mean_expr", "cv_expr", "Baseline - Basic Cleaning Only", log_scale = TRUE)
if (!is.null(p_cahn_baseline)) all_plots_for_pdf[["cahn_baseline"]] <- p_cahn_baseline

p_bewick_baseline <- add_plot_if_data("bewick_baseline", plot_results, "bewick_group", "mean_expr", "cv_expr", "Baseline - Basic Cleaning Only", log_scale = TRUE)
if (!is.null(p_bewick_baseline)) all_plots_for_pdf[["bewick_baseline"]] <- p_bewick_baseline


# Scenario 2: James Lloyd 10% Filter (Log Scale)
p_cahn_jl_10_percent <- add_plot_if_data("cahn_jl_10_percent", plot_results_jl_10_percent_filtered, "cahn_group", "mean_expr", "cv_expr", "James Lloyd 10% Filtered", log_scale = TRUE)
if (!is.null(p_cahn_jl_10_percent)) all_plots_for_pdf[["cahn_jl_10_percent"]] <- p_cahn_jl_10_percent

p_bewick_jl_10_percent <- add_plot_if_data("bewick_jl_10_percent", plot_results_jl_10_percent_filtered, "bewick_group", "mean_expr", "cv_expr", "James Lloyd 10% Filtered", log_scale = TRUE)
if (!is.null(p_bewick_jl_10_percent)) all_plots_for_pdf[["bewick_jl_10_percent"]] <- p_bewick_jl_10_percent


# Scenario 3: James Lloyd 20% Filter (Log Scale)
p_cahn_jl_20_percent <- add_plot_if_data("cahn_jl_20_percent", plot_results_jl_20_percent_filtered, "cahn_group", "mean_expr", "cv_expr", "James Lloyd 20% Filtered", log_scale = TRUE)
if (!is.null(p_cahn_jl_20_percent)) all_plots_for_pdf[["cahn_jl_20_percent"]] <- p_cahn_jl_20_percent

p_bewick_jl_20_percent <- add_plot_if_data("bewick_jl_20_percent", plot_results_jl_20_percent_filtered, "bewick_group", "mean_expr", "cv_expr", "James Lloyd 20% Filtered", log_scale = TRUE)
if (!is.null(p_bewick_jl_20_percent)) all_plots_for_pdf[["bewick_jl_20_percent"]] <- p_bewick_jl_20_percent

# Scenario 4: Mean > 0.5 Filter (Log Scale)
p_cahn_mean_gt_0.5 <- add_plot_if_data("cahn_mean_gt_0.5", plot_results_mean_gt_0.5_filtered, "cahn_group", "mean_expr", "cv_expr", "Mean > 0.5 Filtered", log_scale = TRUE)
if (!is.null(p_cahn_mean_gt_0.5)) all_plots_for_pdf[["cahn_mean_gt_0.5"]] <- p_cahn_mean_gt_0.5

p_bewick_mean_gt_0.5 <- add_plot_if_data("bewick_mean_gt_0.5", plot_results_mean_gt_0.5_filtered, "bewick_group", "mean_expr", "cv_expr", "Mean > 0.5 Filtered", log_scale = TRUE)
if (!is.null(p_bewick_mean_gt_0.5)) all_plots_for_pdf[["bewick_mean_gt_0.5"]] <- p_bewick_mean_gt_0.5

# Scenario 5: Mean > 1 Filter (Log Scale)
p_cahn_mean_gt_1 <- add_plot_if_data("cahn_mean_gt_1", plot_results_mean_gt_1_filtered, "cahn_group", "mean_expr", "cv_expr", "Mean > 1 Filtered", log_scale = TRUE)
if (!is.null(p_cahn_mean_gt_1)) all_plots_for_pdf[["cahn_mean_gt_1"]] <- p_cahn_mean_gt_1

p_bewick_mean_gt_1 <- add_plot_if_data("bewick_mean_gt_1", plot_results_mean_gt_1_filtered, "bewick_group", "mean_expr", "cv_expr", "Mean > 1 Filtered", log_scale = TRUE)
if (!is.null(p_bewick_mean_gt_1)) all_plots_for_pdf[["bewick_mean_gt_1"]] <- p_bewick_mean_gt_1


# --- NEW: Mean Expression vs CV plots (NO Log Scale) ---

# James Lloyd 10% Filter Only (No Log Scale)
p_cahn_jl_10_percent_nolog <- add_plot_if_data("cahn_jl_10_percent_nolog", plot_results_jl_10_percent_filtered, "cahn_group", "mean_expr", "cv_expr", "James Lloyd 10% Filtered, No Log Scale", log_scale = FALSE)
if (!is.null(p_cahn_jl_10_percent_nolog)) all_plots_for_pdf[["cahn_jl_10_percent_nolog"]] <- p_cahn_jl_10_percent_nolog

p_bewick_jl_10_percent_nolog <- add_plot_if_data("bewick_jl_10_percent_nolog", plot_results_jl_10_percent_filtered, "bewick_group", "mean_expr", "cv_expr", "James Lloyd 10% Filtered, No Log Scale", log_scale = FALSE)
if (!is.null(p_bewick_jl_10_percent_nolog)) all_plots_for_pdf[["bewick_jl_10_percent_nolog"]] <- p_bewick_jl_10_percent_nolog

# James Lloyd 20% Filter Only (No Log Scale)
p_cahn_jl_20_percent_nolog <- add_plot_if_data("cahn_jl_20_percent_nolog", plot_results_jl_20_percent_filtered, "cahn_group", "mean_expr", "cv_expr", "James Lloyd 20% Filtered, No Log Scale", log_scale = FALSE)
if (!is.null(p_cahn_jl_20_percent_nolog)) all_plots_for_pdf[["cahn_jl_20_percent_nolog"]] <- p_cahn_jl_20_percent_nolog

p_bewick_jl_20_percent_nolog <- add_plot_if_data("bewick_jl_20_percent_nolog", plot_results_jl_20_percent_filtered, "bewick_group", "mean_expr", "cv_expr", "James Lloyd 20% Filtered, No Log Scale", log_scale = FALSE)
if (!is.null(p_bewick_jl_20_percent_nolog)) all_plots_for_pdf[["bewick_jl_20_percent_nolog"]] <- p_bewick_jl_20_percent_nolog

# Mean > 0.5 Filter (No Log Scale)
p_cahn_mean_gt_0.5_nolog <- add_plot_if_data("cahn_mean_gt_0.5_nolog", plot_results_mean_gt_0.5_filtered, "cahn_group", "mean_expr", "cv_expr", "Mean > 0.5 Filtered, No Log Scale", log_scale = FALSE)
if (!is.null(p_cahn_mean_gt_0.5_nolog)) all_plots_for_pdf[["cahn_mean_gt_0.5_nolog"]] <- p_cahn_mean_gt_0.5_nolog

p_bewick_mean_gt_0.5_nolog <- add_plot_if_data("bewick_mean_gt_0.5_nolog", plot_results_mean_gt_0.5_filtered, "bewick_group", "mean_expr", "cv_expr", "Mean > 0.5 Filtered, No Log Scale", log_scale = FALSE)
if (!is.null(p_bewick_mean_gt_0.5_nolog)) all_plots_for_pdf[["bewick_mean_gt_0.5_nolog"]] <- p_bewick_mean_gt_0.5_nolog

# Mean > 1 Filter (No Log Scale)
p_cahn_mean_gt_1_nolog <- add_plot_if_data("cahn_mean_gt_1_nolog", plot_results_mean_gt_1_filtered, "cahn_group", "mean_expr", "cv_expr", "Mean > 1 Filtered, No Log Scale", log_scale = FALSE)
if (!is.null(p_cahn_mean_gt_1_nolog)) all_plots_for_pdf[["cahn_mean_gt_1_nolog"]] <- p_cahn_mean_gt_1_nolog

p_bewick_mean_gt_1_nolog <- add_plot_if_data("bewick_mean_gt_1_nolog", plot_results_mean_gt_1_filtered, "bewick_group", "mean_expr", "cv_expr", "Mean > 1 Filtered, No Log Scale", log_scale = FALSE)
if (!is.null(p_bewick_mean_gt_1_nolog)) all_plots_for_pdf[["bewick_mean_gt_1_nolog"]] <- p_bewick_mean_gt_1_nolog


# --- STATISTICAL ANALYSIS FOR FILTERING METHODS AND GRID PLOTS ---

# Initialize a list to store all statistical results for grid plots
all_grid_stats_df_rows_main_pdf <- list()
stat_counter_main_pdf <- 0

# Function to perform Kruskal-Wallis and post-hoc Wilcoxon tests for main PDF summary
perform_group_stats_main_pdf <- function(data_df, value_col, group_col_name, context_name, filter_tag, current_ct_ln) {
  # Create a copy of the dataframe to avoid modifying the original
  d_test <- data_df

  # Handle Bewick's TE-like grouping by combining CHH and CHG for statistical comparison
  if (group_col_name == "bewick_group") {
    d_test$temp_group_for_stats <- as.character(d_test[[group_col_name]])
    d_test$temp_group_for_stats[d_test$temp_group_for_stats %in% c("CHH", "CHG")] <- "TE-like"
    # Re-factor to ensure correct levels and order for comparison
    d_test$temp_group_for_stats <- factor(d_test$temp_group_for_stats, levels=c("Unmethylated", "gbM", "TE-like"))
    test_group_col <- "temp_group_for_stats"
  } else {
    test_group_col <- group_col_name
  }

  # Ensure the data_df is not empty and has the required columns
  if (is.null(d_test) || nrow(d_test) == 0 || !value_col %in% colnames(d_test) || !test_group_col %in% colnames(d_test)) {
    return(NULL)
  }

  # Filter out NA values for the group_col and value_col
  d_test <- d_test[!is.na(d_test[[test_group_col]]) & !is.na(d_test[[value_col]]),]

  # Ensure there's more than one group and at least 2 data points per group for meaningful tests
  if (length(unique(d_test[[test_group_col]])) < 2 || any(table(d_test[[test_group_col]]) < 2)) {
    return(NULL)
  }

  kruskal_result <- tryCatch({
    kruskal.test(d_test[[value_col]] ~ d_test[[test_group_col]])
  }, error = function(e) {
    warning(paste("Kruskal-Wallis error for", context_name, filter_tag, current_ct_ln, ":", e$message))
    return(NULL)
  })

  if (is.null(kruskal_result)) {
    return(NULL)
  }

  stat_counter_main_pdf <<- stat_counter_main_pdf + 1
  result_entry <- data.frame(
    ID = stat_counter_main_pdf,
    Context = context_name,
    Filter_Tag = filter_tag,
    Group_Type = group_col_name, # Keep original group name for clarity
    Value_Type = value_col,
    Specific_Context = current_ct_ln,
    Test = "Kruskal-Wallis",
    Comparison = paste(levels(as.factor(d_test[[test_group_col]])), collapse = " vs "),
    P_value = formatC(kruskal_result$p.value, format="e", digits=4), # Format for display
    Statistic = round(kruskal_result$statistic, 3),
    stringsAsFactors = FALSE
  )
  all_grid_stats_df_rows_main_pdf[[length(all_grid_stats_df_rows_main_pdf) + 1]] <<- result_entry

  if (kruskal_result$p.value < 0.05) {
    pairwise_result <- tryCatch({
      pairwise.wilcox.test(d_test[[value_col]], d_test[[test_group_col]], p.adjust.method = "BH")
    }, error = function(e) {
      warning(paste("Pairwise Wilcoxon error for", context_name, filter_tag, current_ct_ln, ":", e$message))
      return(NULL)
    })

    if (!is.null(pairwise_result)) {
      p_matrix <- pairwise_result$p.value
      # Convert p_matrix to a long format dataframe for easier processing
      p_df <- as.data.frame(p_matrix) %>%
        tibble::rownames_to_column("Group1") %>%
        pivot_longer(cols = -Group1, names_to = "Group2", values_to = "P_value") %>%
        filter(!is.na(P_value))

      for (i in 1:nrow(p_df)) {
        stat_counter_main_pdf <<- stat_counter_main_pdf + 1
        pairwise_entry <- data.frame(
          ID = stat_counter_main_pdf,
          Context = context_name,
          Filter_Tag = filter_tag,
          Group_Type = group_col_name,
          Value_Type = value_col,
          Specific_Context = current_ct_ln,
          Test = "Pairwise Wilcoxon (BH-adjusted)",
          Comparison = paste(p_df$Group2[i], "vs", p_df$Group1[i]),
          P_value = formatC(p_df$P_value[i], format="e", digits=4), # Format for display
          Statistic = NA, # Pairwise Wilcoxon doesn't directly return a single statistic
          stringsAsFactors = FALSE
        )
        all_grid_stats_df_rows_main_pdf[[length(all_grid_stats_df_rows_main_pdf) + 1]] <<- pairwise_entry
      }
    }
  }
}

# Iterate and collect stats for Mean Expression by Celltype/Lineage (for main PDF)
for (ct in celltype_levels) { 
  if (is.na(ct) || ct == "") next
  d_unfiltered <- plot_results
  d_unfiltered$cell_mean <- gene_means_celltype[d_unfiltered$gene, ct]
  perform_group_stats_main_pdf(d_unfiltered, "cell_mean", "cahn_group", "Mean Expr by Celltype", "Unfiltered", ct)
  perform_group_stats_main_pdf(d_unfiltered, "cell_mean", "bewick_group", "Mean Expr by Celltype", "Unfiltered", ct)

  d_jl_10_percent_filtered <- plot_results_jl_10_percent_filtered
  d_jl_10_percent_filtered$cell_mean <- gene_means_celltype[d_jl_10_percent_filtered$gene, ct]
  perform_group_stats_main_pdf(d_jl_10_percent_filtered, "cell_mean", "cahn_group", "Mean Expr by Celltype", "James Lloyd 10% Filtered", ct)
  perform_group_stats_main_pdf(d_jl_10_percent_filtered, "cell_mean", "bewick_group", "Mean Expr by Celltype", "James Lloyd 10% Filtered", ct)

  d_jl_20_percent_filtered <- plot_results_jl_20_percent_filtered
  d_jl_20_percent_filtered$cell_mean <- gene_means_celltype[d_jl_20_percent_filtered$gene, ct]
  perform_group_stats_main_pdf(d_jl_20_percent_filtered, "cell_mean", "cahn_group", "Mean Expr by Celltype", "James Lloyd 20% Filtered", ct)
  perform_group_stats_main_pdf(d_jl_20_percent_filtered, "cell_mean", "bewick_group", "Mean Expr by Celltype", "James Lloyd 20% Filtered", ct)

  d_mean_gt_0.5_filtered <- plot_results_mean_gt_0.5_filtered
  d_mean_gt_0.5_filtered$cell_mean <- gene_means_celltype[d_mean_gt_0.5_filtered$gene, ct]
  perform_group_stats_main_pdf(d_mean_gt_0.5_filtered, "cell_mean", "cahn_group", "Mean Expr by Celltype", "Mean > 0.5 Filtered", ct)
  perform_group_stats_main_pdf(d_mean_gt_0.5_filtered, "cell_mean", "bewick_group", "Mean Expr by Celltype", "Mean > 0.5 Filtered", ct)

  d_mean_gt_1_filtered <- plot_results_mean_gt_1_filtered
  d_mean_gt_1_filtered$cell_mean <- gene_means_celltype[d_mean_gt_1_filtered$gene, ct]
  perform_group_stats_main_pdf(d_mean_gt_1_filtered, "cell_mean", "cahn_group", "Mean Expr by Celltype", "Mean > 1 Filtered", ct)
  perform_group_stats_main_pdf(d_mean_gt_1_filtered, "cell_mean", "bewick_group", "Mean Expr by Celltype", "Mean > 1 Filtered", ct)
}

for (ln in lineage_levels) {
  if (is.na(ln) || ln == "") next
  d_unfiltered <- plot_results
  d_unfiltered$lineage_mean <- gene_means_lineage[d_unfiltered$gene, ln]
  perform_group_stats_main_pdf(d_unfiltered, "lineage_mean", "cahn_group", "Mean Expr by Lineage", "Unfiltered", ln)
  perform_group_stats_main_pdf(d_unfiltered, "lineage_mean", "bewick_group", "Mean Expr by Lineage", "Unfiltered", ln)

  d_jl_10_percent_filtered <- plot_results_jl_10_percent_filtered
  d_jl_10_percent_filtered$lineage_mean <- gene_means_lineage[d_jl_10_percent_filtered$gene, ln]
  perform_group_stats_main_pdf(d_jl_10_percent_filtered, "lineage_mean", "cahn_group", "Mean Expr by Lineage", "James Lloyd 10% Filtered", ln)
  perform_group_stats_main_pdf(d_jl_10_percent_filtered, "lineage_mean", "bewick_group", "Mean Expr by Lineage", "James Lloyd 10% Filtered", ln)

  d_jl_20_percent_filtered <- plot_results_jl_20_percent_filtered
  d_jl_20_percent_filtered$lineage_mean <- gene_means_lineage[d_jl_20_percent_filtered$gene, ln]
  perform_group_stats_main_pdf(d_jl_20_percent_filtered, "lineage_mean", "cahn_group", "Mean Expr by Lineage", "James Lloyd 20% Filtered", ln)
  perform_group_stats_main_pdf(d_jl_20_percent_filtered, "lineage_mean", "bewick_group", "Mean Expr by Lineage", "James Lloyd 20% Filtered", ln)

  d_mean_gt_0.5_filtered <- plot_results_mean_gt_0.5_filtered
  d_mean_gt_0.5_filtered$lineage_mean <- gene_means_lineage[d_mean_gt_0.5_filtered$gene, ln]
  perform_group_stats_main_pdf(d_mean_gt_0.5_filtered, "lineage_mean", "cahn_group", "Mean Expr by Lineage", "Mean > 0.5 Filtered", ln)
  perform_group_stats_main_pdf(d_mean_gt_0.5_filtered, "lineage_mean", "bewick_group", "Mean Expr by Lineage", "Mean > 0.5 Filtered", ln)

  d_mean_gt_1_filtered <- plot_results_mean_gt_1_filtered
  d_mean_gt_1_filtered$lineage_mean <- gene_means_lineage[d_mean_gt_1_filtered$gene, ln]
  perform_group_stats_main_pdf(d_mean_gt_1_filtered, "lineage_mean", "cahn_group", "Mean Expr by Lineage", "Mean > 1 Filtered", ln)
  perform_group_stats_main_pdf(d_mean_gt_1_filtered, "lineage_mean", "bewick_group", "Mean Expr by Lineage", "Mean > 1 Filtered", ln)
}

# Iterate and collect stats for CV Expression by Celltype/Lineage (for main PDF)
for (ct in celltype_levels) {
  if (is.na(ct) || ct == "") next
  d_ct_unfiltered <- gene_data_by_celltype[[ct]]
  if (!is.null(d_ct_unfiltered) && nrow(d_ct_unfiltered) > 0) {
    perform_group_stats_main_pdf(d_ct_unfiltered, "cv_expr_ct", "cahn_group", "CV Expr by Celltype", "Unfiltered", ct)
    perform_group_stats_main_pdf(d_ct_unfiltered, "cv_expr_ct", "bewick_group", "CV Expr by Celltype", "Unfiltered", ct)
  }

  d_ct_jl_10_percent_filtered <- d_ct_unfiltered %>% filter(gene %in% plot_results_jl_10_percent_filtered$gene)
  if (!is.null(d_ct_jl_10_percent_filtered) && nrow(d_ct_jl_10_percent_filtered) > 0) {
    perform_group_stats_main_pdf(d_ct_jl_10_percent_filtered, "cv_expr_ct", "cahn_group", "CV Expr by Celltype", "James Lloyd 10% Filtered", ct)
    perform_group_stats_main_pdf(d_ct_jl_10_percent_filtered, "cv_expr_ct", "bewick_group", "CV Expr by Celltype", "James Lloyd 10% Filtered", ct)
  }

  d_ct_jl_20_percent_filtered <- d_ct_unfiltered %>% filter(gene %in% plot_results_jl_20_percent_filtered$gene)
  if (!is.null(d_ct_jl_20_percent_filtered) && nrow(d_ct_jl_20_percent_filtered) > 0) {
    perform_group_stats_main_pdf(d_ct_jl_20_percent_filtered, "cv_expr_ct", "cahn_group", "CV Expr by Celltype", "James Lloyd 20% Filtered", ct)
    perform_group_stats_main_pdf(d_ct_jl_20_percent_filtered, "cv_expr_ct", "bewick_group", "CV Expr by Celltype", "James Lloyd 20% Filtered", ct)
  }

  d_ct_mean_gt_0.5_filtered <- d_ct_unfiltered %>% filter(gene %in% plot_results_mean_gt_0.5_filtered$gene)
  if (!is.null(d_ct_mean_gt_0.5_filtered) && nrow(d_ct_mean_gt_0.5_filtered) > 0) {
    perform_group_stats_main_pdf(d_ct_mean_gt_0.5_filtered, "cv_expr_ct", "cahn_group", "CV Expr by Celltype", "Mean > 0.5 Filtered", ct)
    perform_group_stats_main_pdf(d_ct_mean_gt_0.5_filtered, "cv_expr_ct", "bewick_group", "CV Expr by Celltype", "Mean > 0.5 Filtered", ct)
  }

  d_ct_mean_gt_1_filtered <- d_ct_unfiltered %>% filter(gene %in% plot_results_mean_gt_1_filtered$gene)
  if (!is.null(d_ct_mean_gt_1_filtered) && nrow(d_ct_mean_gt_1_filtered) > 0) {
    perform_group_stats_main_pdf(d_ct_mean_gt_1_filtered, "cv_expr_ct", "cahn_group", "CV Expr by Celltype", "Mean > 1 Filtered", ct)
    perform_group_stats_main_pdf(d_ct_mean_gt_1_filtered, "cv_expr_ct", "bewick_group", "CV Expr by Celltype", "Mean > 1 Filtered", ct)
  }
}

for (ln in lineage_levels) {
  if (is.na(ln) || ln == "") next
  d_ln_unfiltered <- gene_data_by_lineage[[ln]]
  if (!is.null(d_ln_unfiltered) && nrow(d_ln_unfiltered) > 0) {
    perform_group_stats_main_pdf(d_ln_unfiltered, "cv_expr_ln", "cahn_group", "CV Expr by Lineage", "Unfiltered", ln)
    perform_group_stats_main_pdf(d_ln_unfiltered, "cv_expr_ln", "bewick_group", "CV Expr by Lineage", "Unfiltered", ln)
  }

  d_ln_jl_10_percent_filtered <- d_ln_unfiltered %>% filter(gene %in% plot_results_jl_10_percent_filtered$gene)
  if (!is.null(d_ln_jl_10_percent_filtered) && nrow(d_ln_jl_10_percent_filtered) > 0) {
    perform_group_stats_main_pdf(d_ln_jl_10_percent_filtered, "cv_expr_ln", "cahn_group", "CV Expr by Lineage", "James Lloyd 10% Filtered", ln)
    perform_group_stats_main_pdf(d_ln_jl_10_percent_filtered, "cv_expr_ln", "bewick_group", "CV Expr by Lineage", "James Lloyd 10% Filtered", ln)
  }

  d_ln_jl_20_percent_filtered <- d_ln_unfiltered %>% filter(gene %in% plot_results_jl_20_percent_filtered$gene)
  if (!is.null(d_ln_jl_20_percent_filtered) && nrow(d_ln_jl_20_percent_filtered) > 0) {
    perform_group_stats_main_pdf(d_ln_jl_20_percent_filtered, "cv_expr_ln", "cahn_group", "CV Expr by Lineage", "James Lloyd 20% Filtered", ln)
    perform_group_stats_main_pdf(d_ln_jl_20_percent_filtered, "cv_expr_ln", "bewick_group", "CV Expr by Lineage", "James Lloyd 20% Filtered", ln)
  }

  d_ln_mean_gt_0.5_filtered <- d_ln_unfiltered %>% filter(gene %in% plot_results_mean_gt_0.5_filtered$gene)
  if (!is.null(d_ln_mean_gt_0.5_filtered) && nrow(d_ln_mean_gt_0.5_filtered) > 0) {
    perform_group_stats_main_pdf(d_ln_mean_gt_0.5_filtered, "cv_expr_ln", "cahn_group", "CV Expr by Lineage", "Mean > 0.5 Filtered", ln)
    perform_group_stats_main_pdf(d_ln_mean_gt_0.5_filtered, "cv_expr_ln", "bewick_group", "CV Expr by Lineage", "Mean > 0.5 Filtered", ln)
  }

  d_ln_mean_gt_1_filtered <- d_ln_unfiltered %>% filter(gene %in% plot_results_mean_gt_1_filtered$gene)
  if (!is.null(d_ln_mean_gt_1_filtered) && nrow(d_ln_mean_gt_1_filtered) > 0) {
    perform_group_stats_main_pdf(d_ln_mean_gt_1_filtered, "cv_expr_ln", "cahn_group", "CV Expr by Lineage", "Mean > 1 Filtered", ln)
    perform_group_stats_main_pdf(d_ln_mean_gt_1_filtered, "cv_expr_ln", "bewick_group", "CV Expr by Lineage", "Mean > 1 Filtered", ln)
  }
}


# Combine all collected stats into a single data frame for main PDF
final_grid_stats_df_main_pdf <- if (length(all_grid_stats_df_rows_main_pdf) > 0) {
  do.call(rbind, all_grid_stats_df_rows_main_pdf)
} else {
  data.frame(ID=numeric(0), Context=character(0), Filter_Tag=character(0), Group_Type=character(0),
             Value_Type=character(0), Specific_Context=character(0), Test=character(0),
             Comparison=character(0), P_value=character(0), Statistic=numeric(0))
}


# Summarize the counts and median values for each filter category
filter_category_summary <- plot_results %>%
  group_by(filter_category_jl) %>% # Use the new JL filter category
  summarise(
    Count = n(),
    Median_CV = median(cv_expr, na.rm = TRUE),
    Median_Mean_Expr = median(mean_expr, na.rm = TRUE)
  ) %>%
  mutate(
    Median_CV = formatC(Median_CV, format="e", digits=2),
    Median_Mean_Expr = formatC(Median_Mean_Expr, format="e", digits=2)
  )

# Original summary statistics
stat_table <- data.frame(
  Comparison = c(
    "Cahn gbM vs non-gbM",
    "Bewick gbM vs non-gbM",
    "Bewick TE-like vs other",
    "H2A.Z-Depleted vs H2A.Z-Enriched"
  ),
  P_value = c(
    formatC(stat_cahn$p.value, format="e", digits=4),
    formatC(stat_bewick$p.value, format="e", digits=4),
    formatC(stat_te_like$p.value, format="e", digits=4),
    formatC(stat_h2az$p.value, format="e", digits=4)
  ),
  Statistic = c(
    round(stat_cahn$statistic, 3),
    round(stat_bewick$statistic, 3),
    round(stat_te_like$statistic, 3),
    round(stat_h2az$statistic, 3)
  ),
  stringsAsFactors = FALSE
)

# --- Generate Tables for PDF ---
all_tables_for_pdf <- list()

# Table 1: Filter Category Summary
filter_summary_grob <- tableGrob(filter_category_summary, rows = NULL, theme = ttheme_default(base_size = 10))
filter_summary_title <- textGrob("Summary of Genes by Filter Category", gp = gpar(fontsize = 14, fontface = "bold"))
all_tables_for_pdf[[length(all_tables_for_pdf) + 1]] <- arrangeGrob(filter_summary_title, filter_summary_grob, ncol = 1, heights = c(0.1, 0.9))

# Table 2: Overall Methylation Group Statistics
overall_stats_grob <- tableGrob(stat_table, rows = NULL, theme = ttheme_default(base_size = 10))
overall_stats_title <- textGrob("Overall Methylation Group Comparison Statistics", gp = gpar(fontsize = 14, fontface = "bold"))
all_tables_for_pdf[[length(all_tables_for_pdf) + 1]] <- arrangeGrob(overall_stats_title, overall_stats_grob, ncol = 1, heights = c(0.1, 0.9))

# Table 3: Grid Plot Statistics (Paginated if necessary) - for Main PDF
if (nrow(final_grid_stats_df_main_pdf) > 0) {
  rows_per_page <- 40 # Adjust as needed for readability
  num_pages <- ceiling(nrow(final_grid_stats_df_main_pdf) / rows_per_page)

  for (p in 1:num_pages) {
    start_row <- (p - 1) * rows_per_page + 1
    end_row <- min(p * rows_per_page, nrow(final_grid_stats_df_main_pdf))
    page_data <- final_grid_stats_df_main_pdf[start_row:end_row, ]

    grid_stats_grob <- tableGrob(page_data, rows = NULL, theme = ttheme_default(base_size = 8))
    grid_stats_title <- textGrob(paste0("Grid Plot Statistics (Page ", p, " of ", num_pages, ")"), gp = gpar(fontsize = 14, fontface = "bold"))
    all_tables_for_pdf[[length(all_tables_for_pdf) + 1]] <- arrangeGrob(grid_stats_title, grid_stats_grob, ncol = 1, heights = c(0.1, 0.9))
  }
} else {
  warning("No grid plot statistics to add to main PDF.")
}


# --- Save all plots and tables to a single PDF ---
pdf_filename <- file.path(output_dir, "mean_cv_filtering_comparison.pdf")

if (length(all_plots_for_pdf) > 0 || length(all_tables_for_pdf) > 0) {
  pdf(pdf_filename, width=10, height=8)

  # Print all plots
  for (plot_obj in all_plots_for_pdf) {
    if (!is.null(plot_obj)) {
      print(plot_obj)
    }
  }

  # Print all tables
  for (table_grob_obj in all_tables_for_pdf) {
    if (!is.null(table_grob_obj)) {
      grid.newpage()
      grid.draw(table_grob_obj)
    }
  }

  dev.off()
  message(paste("All plots and statistics saved to:", pdf_filename))
} else {
  message(paste("Skipping main PDF generation for", pdf_filename, "as no plots or tables are available."))
}


# --- NEW: Save James Lloyd Filtered CV Boxplots to separate PDFs ---
# This section was moved to the end to avoid confusion with the main PDF generation.
# The plotting logic inside the loops was already corrected earlier.

# --- NEW STANDALONE PNG: Mean Expression vs CV (James Lloyd Filtered, All Cells) ---
p_jl_10_percent_filtered_all_cells_cahn <- add_plot_if_data("cahn_jl_10_percent_all_cells", plot_results_jl_10_percent_filtered, "cahn_group", "mean_expr", "cv_expr", "James Lloyd 10% Filtered, All Cells", log_scale = FALSE)
if (!is.null(p_jl_10_percent_filtered_all_cells_cahn)) {
  ggsave(file.path(output_dir, "mean_cv_cahn_james_lloyd_10_percent_filtered_all_cells.png"), plot=p_jl_10_percent_filtered_all_cells_cahn, width=10, height=8)
}

p_jl_10_percent_filtered_all_cells_bewick <- add_plot_if_data("bewick_jl_10_percent_all_cells", plot_results_jl_10_percent_filtered, "bewick_group", "mean_expr", "cv_expr", "James Lloyd 10% Filtered, All Cells", log_scale = FALSE)
if (!is.null(p_jl_10_percent_filtered_all_cells_bewick)) {
  ggsave(file.path(output_dir, "mean_cv_bewick_james_lloyd_10_percent_filtered_all_cells.png"), plot=p_jl_10_percent_filtered_all_cells_bewick, width=10, height=8)
}

p_jl_20_percent_filtered_all_cells_cahn <- add_plot_if_data("cahn_jl_20_percent_all_cells", plot_results_jl_20_percent_filtered, "cahn_group", "mean_expr", "cv_expr", "James Lloyd 20% Filtered, All Cells", log_scale = FALSE)
if (!is.null(p_jl_20_percent_filtered_all_cells_cahn)) {
  ggsave(file.path(output_dir, "mean_cv_cahn_james_lloyd_20_percent_filtered_all_cells.png"), plot=p_jl_20_percent_filtered_all_cells_cahn, width=10, height=8)
}

p_jl_20_percent_filtered_all_cells_bewick <- add_plot_if_data("bewick_jl_20_percent_all_cells", plot_results_jl_20_percent_filtered, "bewick_group", "mean_expr", "cv_expr", "James Lloyd 20% Filtered, All Cells", log_scale = FALSE)
if (!is.null(p_jl_20_percent_filtered_all_cells_bewick)) {
  ggsave(file.path(output_dir, "mean_cv_bewick_james_lloyd_20_percent_filtered_all_cells.png"), plot=p_jl_20_percent_filtered_all_cells_bewick, width=10, height=8)
}

p_mean_gt_0.5_filtered_all_cells_cahn <- add_plot_if_data("cahn_mean_gt_0.5_all_cells", plot_results_mean_gt_0.5_filtered, "cahn_group", "mean_expr", "cv_expr", "Mean > 0.5 Filtered, All Cells", log_scale = FALSE)
if (!is.null(p_mean_gt_0.5_filtered_all_cells_cahn)) {
  ggsave(file.path(output_dir, "mean_cv_cahn_mean_gt_0.5_filtered_all_cells.png"), plot=p_mean_gt_0.5_filtered_all_cells_cahn, width=10, height=8)
}

p_mean_gt_0.5_filtered_all_cells_bewick <- add_plot_if_data("bewick_mean_gt_0.5_all_cells", plot_results_mean_gt_0.5_filtered, "bewick_group", "mean_expr", "cv_expr", "Mean > 0.5 Filtered, All Cells", log_scale = FALSE)
if (!is.null(p_mean_gt_0.5_filtered_all_cells_bewick)) {
  ggsave(file.path(output_dir, "mean_cv_bewick_mean_gt_0.5_filtered_all_cells.png"), plot=p_mean_gt_0.5_filtered_all_cells_bewick, width=10, height=8)
}

p_mean_gt_1_filtered_all_cells_cahn <- add_plot_if_data("cahn_mean_gt_1_all_cells", plot_results_mean_gt_1_filtered, "cahn_group", "mean_expr", "cv_expr", "Mean > 1 Filtered, All Cells", log_scale = FALSE)
if (!is.null(p_mean_gt_1_filtered_all_cells_cahn)) {
  ggsave(file.path(output_dir, "mean_cv_cahn_mean_gt_1_filtered_all_cells.png"), plot=p_mean_gt_1_filtered_all_cells_cahn, width=10, height=8)
}

p_mean_gt_1_filtered_all_cells_bewick <- add_plot_if_data("bewick_mean_gt_1_all_cells", plot_results_mean_gt_1_filtered, "bewick_group", "mean_expr", "cv_expr", "Mean > 1 Filtered, All Cells", log_scale = FALSE)
if (!is.null(p_mean_gt_1_filtered_all_cells_bewick)) {
  ggsave(file.path(output_dir, "mean_cv_bewick_mean_gt_1_filtered_all_cells.png"), plot=p_mean_gt_1_filtered_all_cells_bewick, width=10, height=8)
}

# --- NEW STANDALONE PNG: CV Boxplots (James Lloyd Filtered, All Cells) ---
p_cv_boxplot_cahn_jl_10_percent_filtered_all_cells <- ggplot(plot_results_jl_10_percent_filtered, aes(x=cahn_group, y=cv_expr, fill=cahn_group)) +
  geom_boxplot(outlier.size=0.5) + ## FIX: `outlier.linewidth` is deprecated, using `outlier.size`
  labs(title="Expression Noise (CV) by Cahn Group (James Lloyd 10% Filtered, All Cells)",
       x="Cahn Group", y="Coefficient of Variation (CV)") +
  theme_bw()
ggsave(file.path(output_dir, "cv_boxplot_cahn_james_lloyd_10_percent_filtered_all_cells.png"), plot=p_cv_boxplot_cahn_jl_10_percent_filtered_all_cells, width=8, height=6)

p_cv_boxplot_bewick_jl_10_percent_filtered_all_cells <- ggplot(plot_results_jl_10_percent_filtered, aes(x=bewick_group, y=cv_expr, fill=bewick_group)) +
  geom_boxplot(outlier.size=0.5) + ## FIX: `outlier.linewidth` is deprecated, using `outlier.size`
  scale_x_discrete(drop = FALSE) + # Added for Bewick plots
  labs(title="Expression Noise (CV) by Bewick Group (James Lloyd 10% Filtered, All Cells)",
       x="Bewick Group", y="Coefficient of Variation (CV)") +
  theme_bw()
ggsave(file.path(output_dir, "cv_boxplot_bewick_james_lloyd_10_percent_filtered_all_cells.png"), plot=p_cv_boxplot_bewick_jl_10_percent_filtered_all_cells, width=8, height=6)

p_cv_boxplot_cahn_jl_20_percent_filtered_all_cells <- ggplot(plot_results_jl_20_percent_filtered, aes(x=cahn_group, y=cv_expr, fill=cahn_group)) +
  geom_boxplot(outlier.size=0.5) + ## FIX: `outlier.linewidth` is deprecated, using `outlier.size`
  labs(title="Expression Noise (CV) by Cahn Group (James Lloyd 20% Filtered, All Cells)",
       x="Cahn Group", y="Coefficient of Variation (CV)") +
  theme_bw()
ggsave(file.path(output_dir, "cv_boxplot_cahn_james_lloyd_20_percent_filtered_all_cells.png"), plot=p_cv_boxplot_cahn_jl_20_percent_filtered_all_cells, width=8, height=6)

p_cv_boxplot_bewick_jl_20_percent_filtered_all_cells <- ggplot(plot_results_jl_20_percent_filtered, aes(x=bewick_group, y=cv_expr, fill=bewick_group)) +
  geom_boxplot(outlier.size=0.5) + ## FIX: `outlier.linewidth` is deprecated, using `outlier.size`
  scale_x_discrete(drop = FALSE) + # Added for Bewick plots
  labs(title="Expression Noise (CV) by Bewick Group (James Lloyd 20% Filtered, All Cells)",
       x="Bewick Group", y="Coefficient of Variation (CV)") +
  theme_bw()
ggsave(file.path(output_dir, "cv_boxplot_bewick_james_lloyd_20_percent_filtered_all_cells.png"), plot=p_cv_boxplot_bewick_jl_20_percent_filtered_all_cells, width=8, height=6)

p_cv_boxplot_cahn_mean_gt_0.5_filtered_all_cells <- ggplot(plot_results_mean_gt_0.5_filtered, aes(x=cahn_group, y=cv_expr, fill=cahn_group)) +
  geom_boxplot(outlier.size=0.5) + ## FIX: `outlier.linewidth` is deprecated, using `outlier.size`
  labs(title="Expression Noise (CV) by Cahn Group (Mean > 0.5 Filtered, All Cells)",
       x="Cahn Group", y="Coefficient of Variation (CV)") +
  theme_bw()
ggsave(file.path(output_dir, "cv_boxplot_cahn_mean_gt_0.5_filtered_all_cells.png"), plot=p_cv_boxplot_cahn_mean_gt_0.5_filtered_all_cells, width=8, height=6)

p_cv_boxplot_bewick_mean_gt_0.5_filtered_all_cells <- ggplot(plot_results_mean_gt_0.5_filtered, aes(x=bewick_group, y=cv_expr, fill=bewick_group)) +
  geom_boxplot(outlier.size=0.5) + ## FIX: `outlier.linewidth` is deprecated, using `outlier.size`
  scale_x_discrete(drop = FALSE) + # Added for Bewick plots
  labs(title="Expression Noise (CV) by Bewick Group (Mean > 0.5 Filtered, All Cells)",
       x="Bewick Group", y="Coefficient of Variation (CV)") +
  theme_bw()
ggsave(file.path(output_dir, "cv_boxplot_bewick_mean_gt_0.5_filtered_all_cells.png"), plot=p_cv_boxplot_bewick_mean_gt_0.5_filtered_all_cells, width=8, height=6)

p_cv_boxplot_cahn_mean_gt_1_filtered_all_cells <- ggplot(plot_results_mean_gt_1_filtered, aes(x=cahn_group, y=cv_expr, fill=cahn_group)) +
  geom_boxplot(outlier.size=0.5) + ## FIX: `outlier.linewidth` is deprecated, using `outlier.size`
  labs(title="Expression Noise (CV) by Cahn Group (Mean > 1 Filtered, All Cells)",
       x="Cahn Group", y="Coefficient of Variation (CV)") +
  theme_bw()
ggsave(file.path(output_dir, "cv_boxplot_cahn_mean_gt_1_filtered_all_cells.png"), plot=p_cv_boxplot_cahn_mean_gt_1_filtered_all_cells, width=8, height=6)

p_cv_boxplot_bewick_mean_gt_1_filtered_all_cells <- ggplot(plot_results_mean_gt_1_filtered, aes(x=bewick_group, y=cv_expr, fill=bewick_group)) +
  geom_boxplot(outlier.size=0.5) + ## FIX: `outlier.linewidth` is deprecated, using `outlier.size`
  scale_x_discrete(drop = FALSE) + # Added for Bewick plots
  labs(title="Expression Noise (CV) by Bewick Group (Mean > 1 Filtered, All Cells)",
       x="Bewick Group", y="Coefficient of Variation (CV)") +
  theme_bw()
ggsave(file.path(output_dir, "cv_boxplot_bewick_mean_gt_1_filtered_all_cells.png"), plot=p_cv_boxplot_bewick_mean_gt_1_filtered_all_cells, width=8, height=6)