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

# Filter out non-finite CV values for plotting (overall data for Brennecke fit)
plot_results <- results[is.finite(results$cv_expr), ]
# Remove extreme outliers (e.g., above 99th percentile) for plotting
cv_cap <- quantile(plot_results$cv_expr, 0.99, na.rm=TRUE)
plot_results <- plot_results[plot_results$cv_expr <= cv_cap, ]

# Brennecke et al. (2013) noise filtering: fit mean--CV^2 relationship
plot_results$cv2 <- plot_results$cv_expr^2
fit <- nls(cv2 ~ a / mean_expr + b, data=plot_results[plot_results$mean_expr > 0,], start=list(a=1, b=0.1))
fit_vals <- predict(fit, newdata=data.frame(mean_expr=plot_results$mean_expr))
residuals <- plot_results$cv2 - fit_vals
threshold <- quantile(residuals, 0.95, na.rm=TRUE)
plot_results$biol_variable <- residuals > threshold
filtered_biol <- plot_results[plot_results$biol_variable & plot_results$mean_expr > 0, ]

# Plots (original individual plots - not saved to PDF but kept for Rmd)
gbm_labels <- c("FALSE"="Non-gbM", "TRUE"="gbM")
te_labels <- c("FALSE"="Other", "TRUE"="TE-like")
p1 <- ggplot(results, aes(x=factor(cahn_gbm), y=cv_expr, fill=factor(cahn_gbm))) +
  geom_boxplot(outlier.size=0.5) +
  scale_x_discrete(labels=gbm_labels) +
  labs(title="Expression Noise (CV) by Cahn gbM Status", x="Cahn gbM", y="CV") +
  theme_bw()
p2 <- ggplot(results, aes(x=factor(bewick_gbm), y=cv_expr, fill=factor(bewick_gbm))) +
  geom_boxplot(outlier.size=0.5) +
  scale_x_discrete(labels=gbm_labels) +
  labs(title="Expression Noise (CV) by Bewick gbM Status", x="Bewick gbM", y="CV") +
  theme_bw()
p3 <- ggplot(results, aes(x=factor(bewick_te_like), y=cv_expr, fill=factor(bewick_te_like))) +
  geom_boxplot(outlier.size=0.5) +
  scale_x_discrete(labels=te_labels) +
  labs(title="Expression Noise (CV) by Bewick TE-like Status", x="Bewick TE-like", y="CV") +
  theme_bw()
p_combined <- ggplot(results, aes(x=methylation_group, y=cv_expr, color=methylation_group)) +
  geom_jitter(width=0.25, alpha=0.5, size=0.7) +
  labs(title="Univariant Scatter Plot: Expression Noise (CV) by Methylation Group", x="Methylation Group", y="CV") +
  theme_bw() +
  theme(axis.text.x = element_text(angle=45, hjust=1))

# Save plot objects for Rmd
save(p1, p2, p3, p_combined, plot_results, filtered_biol, fit, file=file.path(output_dir, "gbm_noise_plots.RData"))

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
  geom_boxplot(outlier.size=0.5) +
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
  geom_boxplot(outlier.size=0.5) +
  labs(title="Expression Noise (CV) by Methylation + H2A.Z Group", x="Methylation + H2A.Z Group", y="CV") +
  theme_bw() +
  theme(axis.text.x = element_text(angle=45, hjust=1))
ggsave(file.path(output_dir, "gbm_noise_boxplot_meth_h2az.png"), width=12, height=6)

# --- Violin plots ---
p_cahn_violin <- ggplot(plot_results[!is.na(plot_results$cahn_group),], aes(x=cahn_group, y=cv_expr, fill=cahn_group)) +
  geom_violin(trim=FALSE, scale="width") +
  geom_boxplot(width=0.1, outlier.size=0.5, outlier.shape=16, outlier.alpha=0.3, fill="white") +
  labs(title="Expression Noise (CV) by Cahn Methylation Group (Violin)", x="Cahn Methylation Group", y="CV") +
  theme_bw() +
  theme(axis.text.x = element_text(angle=45, hjust=1))
ggsave(file.path(output_dir, "gbm_noise_violin_cahn.png"), plot=p_cahn_violin, width=8, height=6)

p_bewick_violin <- ggplot(plot_results[!is.na(plot_results$bewick_group),], aes(x=bewick_group, y=cv_expr, fill=bewick_group)) +
  geom_violin(trim=FALSE, scale="width") +
  geom_boxplot(width=0.1, outlier.size=0.5, outlier.shape=16, outlier.alpha=0.3, fill="white") +
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
  labs(title="CV Distribution by Bewick Methylation Group (Ridge)", x="CV", y="Bewick Methylation Group") +
  theme_bw()
ggsave(file.path(output_dir, "gbm_noise_ridge_bewick.png"), plot=p_bewick_ridge, width=8, height=8)

# --- UpSet plot: Overlap of high-noise, gbM, and H2A.Z-depleted genes ---
upset_data <- data.frame(
  High_Noise = results$cv_expr >= quantile(results$cv_expr, 0.9, na.rm=TRUE),
  gbM = results$cahn_gbm | results$bewick_gbm,
  H2AZ_Depleted = results$h2az_group == "H2A.Z-Depleted"
)
# Fix for the 'intersect' naming conflict
upset_features <- c("High_Noise","gbM","H2AZ_Depleted")
ComplexUpset::upset(upset_data, intersect = upset_features, width_ratio=0.1, set_sizes=ComplexUpset::upset_set_size() + theme(axis.text.x = element_text(angle=45, hjust=1, size=14)))
ggsave(file.path(output_dir, "upset_highnoise_gbm_h2az.png"), width=14, height=7)

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

# --- Boxplots by celltype (Cahn group) - Mean Expression ---
plots_cahn_celltype_unfiltered <- lapply(celltype_levels, function(ct) {
  d <- plot_results # Overall unfiltered gene set
  d$cell_mean <- gene_means_celltype[d$gene, ct] # Add celltype-specific mean
  ggplot(d[!is.na(d$cahn_group),], aes(x=factor(cahn_group), y=cell_mean, fill=factor(cahn_group))) +
    geom_boxplot(outlier.size=0.5) +
    labs(title=paste("Cahn Group in", ct), x="Cahn Group", y="Mean Expression") +
    theme_bw() +
    theme(legend.position = "none", axis.text.x = element_text(angle=45, hjust=1))
})
g_cahn_celltype_unfiltered <- marrangeGrob(plots_cahn_celltype_unfiltered, nrow=2, ncol=3, top="Mean Expression by Cahn Group in Celltypes (Unfiltered)")
pdf(file.path(output_dir, "grid_mean_expr_cahn_by_celltype_unfiltered.pdf"), width=18, height=12)
for (i in 1:length(g_cahn_celltype_unfiltered)) {
  grid.newpage()
  grid.draw(g_cahn_celltype_unfiltered[[i]])
}
dev.off()

plots_cahn_celltype_filtered <- lapply(celltype_levels, function(ct) {
  d <- filtered_biol # Overall filtered gene set
  d$cell_mean <- gene_means_celltype[d$gene, ct] # Add celltype-specific mean
  ggplot(d[!is.na(d$cahn_group),], aes(x=factor(cahn_group), y=cell_mean, fill=factor(cahn_group))) +
    geom_boxplot(outlier.size=0.5) +
    labs(title=paste("Cahn Group in", ct), x="Cahn Group", y="Mean Expression") +
    theme_bw() +
    theme(legend.position = "none", axis.text.x = element_text(angle=45, hjust=1))
})
g_cahn_celltype_filtered <- marrangeGrob(plots_cahn_celltype_filtered, nrow=2, ncol=3, top="Mean Expression by Cahn Group in Celltypes (Filtered)")
pdf(file.path(output_dir, "grid_mean_expr_cahn_by_celltype_filtered.pdf"), width=18, height=12)
for (i in 1:length(g_cahn_celltype_filtered)) {
  grid.newpage()
  grid.draw(g_cahn_celltype_filtered[[i]])
}
dev.off()

# --- Boxplots by celltype (Bewick group) - Mean Expression ---
plots_bewick_celltype_unfiltered <- lapply(celltype_levels, function(ct) {
  d <- plot_results
  d$cell_mean <- gene_means_celltype[d$gene, ct]
  ggplot(d[!is.na(d$bewick_group),], aes(x=factor(bewick_group), y=cell_mean, fill=factor(bewick_group))) +
    geom_boxplot(outlier.size=0.5) +
    labs(title=paste("Bewick Group in", ct), x="Bewick Group", y="Mean Expression") +
    theme_bw() +
    theme(legend.position = "none", axis.text.x = element_text(angle=45, hjust=1))
})
g_bewick_celltype_unfiltered <- marrangeGrob(plots_bewick_celltype_unfiltered, nrow=2, ncol=3, top="Mean Expression by Bewick Group in Celltypes (Unfiltered)")
pdf(file.path(output_dir, "grid_mean_expr_bewick_by_celltype_unfiltered.pdf"), width=18, height=12)
for (i in 1:length(g_bewick_celltype_unfiltered)) {
  grid.newpage()
  grid.draw(g_bewick_celltype_unfiltered[[i]])
}
dev.off()

plots_bewick_celltype_filtered <- lapply(celltype_levels, function(ct) {
  d <- filtered_biol
  d$cell_mean <- gene_means_celltype[d$gene, ct]
  ggplot(d[!is.na(d$bewick_group),], aes(x=factor(bewick_group), y=cell_mean, fill=factor(bewick_group))) +
    geom_boxplot(outlier.size=0.5) +
    labs(title=paste("Bewick Group in", ct), x="Bewick Group", y="Mean Expression") +
    theme_bw() +
    theme(legend.position = "none", axis.text.x = element_text(angle=45, hjust=1))
})
g_bewick_celltype_filtered <- marrangeGrob(plots_bewick_celltype_filtered, nrow=2, ncol=3, top="Mean Expression by Bewick Group in Celltypes (Filtered)")
pdf(file.path(output_dir, "grid_mean_expr_bewick_by_celltype_filtered.pdf"), width=18, height=12)
for (i in 1:length(g_bewick_celltype_filtered)) {
  grid.newpage()
  grid.draw(g_bewick_celltype_filtered[[i]])
}
dev.off()


# --- Boxplots by lineage (Cahn group) - Mean Expression ---
plots_cahn_lineage_unfiltered <- lapply(lineage_levels, function(ln) {
  d <- plot_results
  d$lineage_mean <- gene_means_lineage[d$gene, ln]
  ggplot(d[!is.na(d$cahn_group),], aes(x=factor(cahn_group), y=lineage_mean, fill=factor(cahn_group))) +
    geom_boxplot(outlier.size=0.5) +
    labs(title=paste("Cahn Group in", ln), x="Cahn Group", y="Mean Expression") +
    theme_bw() +
    theme(legend.position = "none", axis.text.x = element_text(angle=45, hjust=1))
})
g_cahn_lineage_unfiltered <- marrangeGrob(plots_cahn_lineage_unfiltered, nrow=2, ncol=3, top="Mean Expression by Cahn Group in Lineages (Unfiltered)")
pdf(file.path(output_dir, "grid_mean_expr_cahn_by_lineage_unfiltered.pdf"), width=18, height=12)
for (i in 1:length(g_cahn_lineage_unfiltered)) {
  grid.newpage()
  grid.draw(g_cahn_lineage_unfiltered[[i]])
}
dev.off()

plots_cahn_lineage_filtered <- lapply(lineage_levels, function(ln) {
  d <- filtered_biol
  d$lineage_mean <- gene_means_lineage[d$gene, ln]
  ggplot(d[!is.na(d$cahn_group),], aes(x=factor(cahn_group), y=lineage_mean, fill=factor(cahn_group))) +
    geom_boxplot(outlier.size=0.5) +
    labs(title=paste("Cahn Group in", ln), x="Cahn Group", y="Mean Expression") +
    theme_bw() +
    theme(legend.position = "none", axis.text.x = element_text(angle=45, hjust=1))
})
g_cahn_lineage_filtered <- marrangeGrob(plots_cahn_lineage_filtered, nrow=2, ncol=3, top="Mean Expression by Cahn Group in Lineages (Filtered)")
pdf(file.path(output_dir, "grid_mean_expr_cahn_by_lineage_filtered.pdf"), width=18, height=12)
for (i in 1:length(g_cahn_lineage_filtered)) {
  grid.newpage()
  grid.draw(g_cahn_lineage_filtered[[i]])
}
dev.off()

# --- Boxplots by lineage (Bewick group) - Mean Expression ---
plots_bewick_lineage_unfiltered <- lapply(lineage_levels, function(ln) {
  d <- plot_results
  d$lineage_mean <- gene_means_lineage[d$gene, ln]
  ggplot(d[!is.na(d$bewick_group),], aes(x=factor(bewick_group), y=lineage_mean, fill=factor(bewick_group))) +
    geom_boxplot(outlier.size=0.5) +
    labs(title=paste("Bewick Group in", ln), x="Bewick Group", y="Mean Expression") +
    theme_bw() +
    theme(legend.position = "none", axis.text.x = element_text(angle=45, hjust=1))
})
g_bewick_lineage_unfiltered <- marrangeGrob(plots_bewick_lineage_unfiltered, nrow=2, ncol=3, top="Mean Expression by Bewick Group in Lineages (Unfiltered)")
pdf(file.path(output_dir, "grid_mean_expr_bewick_by_lineage_unfiltered.pdf"), width=18, height=12)
for (i in 1:length(g_bewick_lineage_unfiltered)) {
  grid.newpage()
  grid.draw(g_bewick_lineage_unfiltered[[i]])
}
dev.off()

plots_bewick_lineage_filtered <- lapply(lineage_levels, function(ln) {
  d <- filtered_biol
  d$lineage_mean <- gene_means_lineage[d$gene, ln]
  ggplot(d[!is.na(d$bewick_group),], aes(x=factor(bewick_group), y=lineage_mean, fill=factor(bewick_group))) +
    geom_boxplot(outlier.size=0.5) +
    labs(title=paste("Bewick Group in", ln), x="Bewick Group", y="Mean Expression") +
    theme_bw() +
    theme(legend.position = "none", axis.text.x = element_text(angle=45, hjust=1))
})
g_bewick_lineage_filtered <- marrangeGrob(plots_bewick_lineage_filtered, nrow=2, ncol=3, top="Mean Expression by Bewick Group in Lineages (Filtered)")
pdf(file.path(output_dir, "grid_mean_expr_bewick_by_lineage_filtered.pdf"), width=18, height=12)
for (i in 1:length(g_bewick_lineage_filtered)) {
  grid.newpage()
  grid.draw(g_bewick_lineage_filtered[[i]])
}
dev.off()


# --- Boxplots by celltype (Cahn group) - Expression Noise (CV) ---
plots_cahn_celltype_cv_unfiltered <- lapply(celltype_levels, function(ct) {
  d_ct <- gene_data_by_celltype[[ct]] # Use celltype-specific data
  if (is.null(d_ct) || nrow(d_ct) == 0) return(NULL) # Handle empty data for a celltype

  ggplot(d_ct[!is.na(d_ct$cahn_group),], aes(x=factor(cahn_group), y=cv_expr_ct, fill=factor(cahn_group))) +
    geom_boxplot(outlier.size=0.5) +
    labs(title=paste("Cahn Group in", ct), x="Cahn Group", y="Expression Noise (CV)") +
    theme_bw() +
    theme(legend.position = "none", axis.text.x = element_text(angle=45, hjust=1))
})
g_cahn_celltype_cv_unfiltered <- marrangeGrob(plots_cahn_celltype_cv_unfiltered, nrow=2, ncol=3, top="Expression Noise (CV) by Cahn Group in Celltypes (Unfiltered)")
pdf(file.path(output_dir, "grid_cv_expr_cahn_by_celltype_unfiltered.pdf"), width=18, height=12)
for (i in 1:length(g_cahn_celltype_cv_unfiltered)) {
  grid.newpage()
  grid.draw(g_cahn_celltype_cv_unfiltered[[i]])
}
dev.off()

plots_cahn_celltype_cv_filtered <- lapply(celltype_levels, function(ct) {
  d_ct <- gene_data_by_celltype[[ct]] # Start with celltype-specific data
  if (is.null(d_ct) || nrow(d_ct) == 0) return(NULL)

  # Filter to include only genes that were overall Brennecke-filtered
  d_filtered_ct <- d_ct %>% filter(gene %in% filtered_biol$gene)
  if (nrow(d_filtered_ct) == 0) return(NULL)

  ggplot(d_filtered_ct[!is.na(d_filtered_ct$cahn_group),], aes(x=factor(cahn_group), y=cv_expr_ct, fill=factor(cahn_group))) +
    geom_boxplot(outlier.size=0.5) +
    labs(title=paste("Cahn Group in", ct), x="Cahn Group", y="Expression Noise (CV) (Brennecke Filtered)") +
    theme_bw() +
    theme(legend.position = "none", axis.text.x = element_text(angle=45, hjust=1))
})
g_cahn_celltype_cv_filtered <- marrangeGrob(plots_cahn_celltype_cv_filtered, nrow=2, ncol=3, top="Expression Noise (CV) by Cahn Group in Celltypes (Brennecke Filtered)")
pdf(file.path(output_dir, "grid_cv_expr_cahn_by_celltype_brennecke_filtered.pdf"), width=18, height=12)
for (i in 1:length(g_cahn_celltype_cv_filtered)) {
  grid.newpage()
  grid.draw(g_cahn_celltype_cv_filtered[[i]])
}
dev.off()

# --- Boxplots by celltype (Bewick group) - Expression Noise (CV) ---
plots_bewick_celltype_cv_unfiltered <- lapply(celltype_levels, function(ct) {
  d_ct <- gene_data_by_celltype[[ct]] # Use celltype-specific data
  if (is.null(d_ct) || nrow(d_ct) == 0) return(NULL)

  ggplot(d_ct[!is.na(d_ct$bewick_group),], aes(x=factor(bewick_group), y=cv_expr_ct, fill=factor(bewick_group))) +
    geom_boxplot(outlier.size=0.5) +
    labs(title=paste("Bewick Group in", ct), x="Bewick Group", y="Expression Noise (CV)") +
    theme_bw() +
    theme(legend.position = "none", axis.text.x = element_text(angle=45, hjust=1))
})
g_bewick_celltype_cv_unfiltered <- marrangeGrob(plots_bewick_celltype_cv_unfiltered, nrow=2, ncol=3, top="Expression Noise (CV) by Bewick Group in Celltypes (Unfiltered)")
pdf(file.path(output_dir, "grid_cv_expr_bewick_by_celltype_unfiltered.pdf"), width=18, height=12)
for (i in 1:length(g_bewick_celltype_cv_unfiltered)) {
  grid.newpage()
  grid.draw(g_bewick_celltype_cv_unfiltered[[i]])
}
dev.off()

plots_bewick_celltype_cv_filtered <- lapply(celltype_levels, function(ct) {
  d_ct <- gene_data_by_celltype[[ct]] # Start with celltype-specific data
  if (is.null(d_ct) || nrow(d_ct) == 0) return(NULL)

  # Filter to include only genes that were overall Brennecke-filtered
  d_filtered_ct <- d_ct %>% filter(gene %in% filtered_biol$gene)
  if (nrow(d_filtered_ct) == 0) return(NULL)

  ggplot(d_filtered_ct[!is.na(d_filtered_ct$bewick_group),], aes(x=factor(bewick_group), y=cv_expr_ct, fill=factor(bewick_group))) +
    geom_boxplot(outlier.size=0.5) +
    labs(title=paste("Bewick Group in", ct), x="Bewick Group", y="Expression Noise (CV) (Brennecke Filtered)") +
    theme_bw() +
    theme(legend.position = "none", axis.text.x = element_text(angle=45, hjust=1))
})
g_bewick_celltype_cv_filtered <- marrangeGrob(plots_bewick_celltype_cv_filtered, nrow=2, ncol=3, top="Expression Noise (CV) by Bewick Group in Celltypes (Brennecke Filtered)")
pdf(file.path(output_dir, "grid_cv_expr_bewick_by_celltype_brennecke_filtered.pdf"), width=18, height=12)
for (i in 1:length(g_bewick_celltype_cv_filtered)) {
  grid.newpage()
  grid.draw(g_bewick_celltype_cv_filtered[[i]])
}
dev.off()


# --- Boxplots by lineage (Cahn group) - Expression Noise (CV) ---
plots_cahn_lineage_cv_unfiltered <- lapply(lineage_levels, function(ln) {
  d_ln <- gene_data_by_lineage[[ln]] # Use lineage-specific data
  if (is.null(d_ln) || nrow(d_ln) == 0) return(NULL)

  ggplot(d_ln[!is.na(d_ln$cahn_group),], aes(x=factor(cahn_group), y=cv_expr_ln, fill=factor(cahn_group))) +
    geom_boxplot(outlier.size=0.5) +
    labs(title=paste("Cahn Group in", ln), x="Cahn Group", y="Expression Noise (CV)") +
    theme_bw() +
    theme(legend.position = "none", axis.text.x = element_text(angle=45, hjust=1))
})
g_cahn_lineage_cv_unfiltered <- marrangeGrob(plots_cahn_lineage_cv_unfiltered, nrow=2, ncol=3, top="Expression Noise (CV) by Cahn Group in Lineages (Unfiltered)")
pdf(file.path(output_dir, "grid_cv_expr_cahn_by_lineage_unfiltered.pdf"), width=18, height=12)
for (i in 1:length(g_cahn_lineage_cv_unfiltered)) {
  grid.newpage()
  grid.draw(g_cahn_lineage_cv_unfiltered[[i]])
}
dev.off()

plots_cahn_lineage_cv_filtered <- lapply(lineage_levels, function(ln) {
  d_ln <- gene_data_by_lineage[[ln]] # Start with lineage-specific data
  if (is.null(d_ln) || nrow(d_ln) == 0) return(NULL)

  # Filter to include only genes that were overall Brennecke-filtered
  d_filtered_ln <- d_ln %>% filter(gene %in% filtered_biol$gene)
  if (nrow(d_filtered_ln) == 0) return(NULL)

  ggplot(d_filtered_ln[!is.na(d_filtered_ln$cahn_group),], aes(x=factor(cahn_group), y=cv_expr_ln, fill=factor(cahn_group))) +
    geom_boxplot(outlier.size=0.5) +
    labs(title=paste("Cahn Group in", ln), x="Cahn Group", y="Expression Noise (CV)") +
    theme_bw() +
    theme(legend.position = "none", axis.text.x = element_text(angle=45, hjust=1))
})
g_cahn_lineage_cv_filtered <- marrangeGrob(plots_cahn_lineage_cv_filtered, nrow=2, ncol=3, top="Expression Noise (CV) by Cahn Group in Lineages (Filtered)")
pdf(file.path(output_dir, "grid_cv_expr_cahn_by_lineage_filtered.pdf"), width=18, height=12)
for (i in 1:length(g_cahn_lineage_cv_filtered)) {
  grid.newpage()
  grid.draw(g_cahn_lineage_cv_filtered[[i]])
}
dev.off()

# --- Boxplots by lineage (Bewick group) - Expression Noise (CV) ---
plots_bewick_lineage_cv_unfiltered <- lapply(lineage_levels, function(ln) {
  d_ln <- gene_data_by_lineage[[ln]] # Use lineage-specific data
  if (is.null(d_ln) || nrow(d_ln) == 0) return(NULL)

  ggplot(d_ln[!is.na(d_ln$bewick_group),], aes(x=factor(bewick_group), y=cv_expr_ln, fill=factor(bewick_group))) +
    geom_boxplot(outlier.size=0.5) +
    labs(title=paste("Bewick Group in", ln), x="Bewick Group", y="Expression Noise (CV)") +
    theme_bw() +
    theme(legend.position = "none", axis.text.x = element_text(angle=45, hjust=1))
})
g_bewick_lineage_cv_unfiltered <- marrangeGrob(plots_bewick_lineage_cv_unfiltered, nrow=2, ncol=3, top="Expression Noise (CV) by Bewick Group in Lineages (Unfiltered)")
pdf(file.path(output_dir, "grid_cv_expr_bewick_by_lineage_unfiltered.pdf"), width=18, height=12)
for (i in 1:length(g_bewick_lineage_cv_unfiltered)) {
  grid.newpage()
  grid.draw(g_bewick_lineage_cv_unfiltered[[i]])
}
dev.off()

plots_bewick_lineage_cv_filtered <- lapply(lineage_levels, function(ln) {
  d_ln <- gene_data_by_lineage[[ln]] # Start with lineage-specific data
  if (is.null(d_ln) || nrow(d_ln) == 0) return(NULL)

  # Filter to include only genes that were overall Brennecke-filtered
  d_filtered_ln <- d_ln %>% filter(gene %in% filtered_biol$gene)
  if (nrow(d_filtered_ln) == 0) return(NULL)

  ggplot(d_filtered_ln[!is.na(d_filtered_ln$bewick_group),], aes(x=factor(bewick_group), y=cv_expr_ln, fill=factor(bewick_group))) +
    geom_boxplot(outlier.size=0.5) +
    labs(title=paste("Bewick Group in", ln), x="Bewick Group", y="Expression Noise (CV)") +
    theme_bw() +
    theme(legend.position = "none", axis.text.x = element_text(angle=45, hjust=1))
})
g_bewick_lineage_cv_filtered <- marrangeGrob(plots_bewick_lineage_cv_filtered, nrow=2, ncol=3, top="Expression Noise (CV) by Bewick Group in Lineages (Filtered)")
pdf(file.path(output_dir, "grid_cv_expr_bewick_by_lineage_filtered.pdf"), width=18, height=12)
for (i in 1:length(g_bewick_lineage_cv_filtered)) {
  grid.newpage()
  grid.draw(g_bewick_lineage_cv_filtered[[i]])
}
dev.off()

# --- NEW SECTION: Filter by bottom 10% least expressed genes and re-plot Mean vs CV ---

# 1. Determine the 10th percentile of mean expression
mean_expr_10th_percentile <- quantile(results$mean_expr, 0.10, na.rm=TRUE)

# 2. Create new dataframes with genes below this threshold removed
#    This dataframe is used for the "James Lloyd Filter Only" scenario
results_mean_filtered_only <- results[results$mean_expr >= mean_expr_10th_percentile, ]

# 3. Apply basic cleaning to results_mean_filtered_only to get plot_results_mean_filtered
plot_results_mean_filtered <- results_mean_filtered_only[is.finite(results_mean_filtered_only$cv_expr), ]
cv_cap_mean_filtered <- quantile(plot_results_mean_filtered$cv_expr, 0.99, na.rm=TRUE)
plot_results_mean_filtered <- plot_results_mean_filtered[plot_results_mean_filtered$cv_expr <= cv_cap_mean_filtered, ]


# 4. Apply Brennecke filter to plot_results_mean_filtered to get filtered_biol_mean_filtered
#    This dataframe is used for the "James Lloyd Filter THEN Brennecke Filter" scenario
filtered_biol_mean_filtered <- data.frame() # Initialize as empty
if (nrow(plot_results_mean_filtered[plot_results_mean_filtered$mean_expr > 0,]) > 2) {
  plot_results_mean_filtered$cv2 <- plot_results_mean_filtered$cv_expr^2
  fit_mean_filtered <- nls(cv2 ~ a / mean_expr + b, data=plot_results_mean_filtered[plot_results_mean_filtered$mean_expr > 0,], start=list(a=1, b=0.1))
  fit_vals_mean_filtered <- predict(fit_mean_filtered, newdata=data.frame(mean_expr=plot_results_mean_filtered$mean_expr))
  residuals_mean_filtered <- plot_results_mean_filtered$cv2 - fit_vals_mean_filtered
  threshold_mean_filtered <- quantile(residuals_mean_filtered, 0.95, na.rm=TRUE)
  plot_results_mean_filtered$biol_variable <- residuals_mean_filtered > threshold_mean_filtered
  filtered_biol_mean_filtered <- plot_results_mean_filtered[plot_results_mean_filtered$biol_variable & plot_results_mean_filtered$mean_expr > 0, ]
} else {
  warning("Not enough data points to perform Brennecke fit after 10% mean expression filtering for combined filter.")
}


# --- Comprehensive Filter Visualization Plot (New plot for first page of PDF) ---

# Recalculate James Lloyd filter status on plot_results (the baseline data)
plot_results$james_lloyd_valid <- plot_results$mean_expr >= mean_expr_10th_percentile

# Create a new column for categorized plotting based on filter validity
plot_results$filter_category <- "Not Valid" # Default

# Apply categories based on filter validity (order matters for assignment)
# Assign most specific categories first to ensure correct overlay
plot_results$filter_category[plot_results$biol_variable & plot_results$james_lloyd_valid] <- "Valid by Both (Brennecke & JL)"
plot_results$filter_category[plot_results$biol_variable & !plot_results$james_lloyd_valid] <- "Valid by Brennecke Only"
plot_results$filter_category[!plot_results$biol_variable & plot_results$james_lloyd_valid] <- "Valid by James Lloyd Only"

# Convert to factor with desired order for legend and plotting
plot_results$filter_category <- factor(plot_results$filter_category,
                                       levels = c("Not Valid",
                                                 "Valid by James Lloyd Only",
                                                 "Valid by Brennecke Only",
                                                 "Valid by Both (Brennecke & JL)"))

# Define colors for these categories
category_colors <- c("Not Valid" = "grey",
                     "Valid by James Lloyd Only" = "darkgreen",
                     "Valid by Brennecke Only" = "blue",
                     "Valid by Both (Brennecke & JL)" = "purple")

# Generate a sequence of mean_expr values for a smooth Brennecke fit line
mean_expr_for_fit_line <- seq(min(plot_results$mean_expr[plot_results$mean_expr > 0], na.rm=TRUE),
                              max(plot_results$mean_expr, na.rm=TRUE), length.out = 100)
fit_line_data <- data.frame(mean_expr = mean_expr_for_fit_line,
                            cv2 = predict(fit, newdata = data.frame(mean_expr = mean_expr_for_fit_line)))

# Create the combined filter visualization plot
p_combined_filter_visualization <- ggplot(plot_results, aes(x = mean_expr, y = cv2)) +
  geom_point(aes(color = filter_category), alpha = 0.5, size = 0.7) +
  scale_color_manual(values = category_colors, name = "Gene Filter Status") +
  geom_line(data = fit_line_data, aes(x = mean_expr, y = cv2), color = "red", size = 1) + # Brennecke fit line
  geom_vline(xintercept = mean_expr_10th_percentile, linetype = "dashed", color = "darkorange", size = 1) + # James Lloyd line
  scale_x_log10(labels = scales::scientific) +
  scale_y_log10(labels = scales::scientific) +
  labs(title = "Mean Expression vs CV^2 with Filtering Methods",
       subtitle = paste0("Brennecke Fit (Red Line), James Lloyd Filter (Orange Dashed Line at ",
                         formatC(mean_expr_10th_percentile, format="e", digits=2), " Mean Expression)"),
       x = "Mean Expression",
       y = "CV^2") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5),
        legend.position = "bottom",
        legend.box = "horizontal") # Place legend at bottom for better readability


# --- Mean Expression vs CV plots for comparison of filtering methods ---
# Collect all plots in a list to save to a single PDF
all_plots_for_pdf <- list()

# Add the new combined filter visualization plot as the first page
all_plots_for_pdf[["combined_filter_overview"]] <- p_combined_filter_visualization


# Scenario 1: Baseline (only basic cleaning)
p_mean_cv_cahn_baseline <- ggplot(plot_results, aes(x=mean_expr, y=cv_expr, color=cahn_group)) +
  geom_point(alpha=0.5, size=0.7) +
  labs(title="Mean vs CV by Cahn Group (Baseline - Basic Cleaning Only)",
       x="Mean Expression", y="Coefficient of Variation (CV)") +
  theme_bw() +
  scale_x_log10() + scale_y_log10()
all_plots_for_pdf[["cahn_baseline"]] <- p_mean_cv_cahn_baseline

p_mean_cv_bewick_baseline <- ggplot(plot_results, aes(x=mean_expr, y=cv_expr, color=bewick_group)) +
  geom_point(alpha=0.5, size=0.7) +
  labs(title="Mean vs CV by Bewick Group (Baseline - Basic Cleaning Only)",
       x="Mean Expression", y="Coefficient of Variation (CV)") +
  theme_bw() +
  scale_x_log10() + scale_y_log10()
all_plots_for_pdf[["bewick_baseline"]] <- p_mean_cv_bewick_baseline


# Scenario 2: Brennecke Filter Only (Log Scale)
p_mean_cv_cahn_brennecke_filtered <- ggplot(filtered_biol, aes(x=mean_expr, y=cv_expr, color=cahn_group)) +
  geom_point(alpha=0.5, size=0.7) +
  labs(title="Mean vs CV by Cahn Group (Brennecke Filter Only)",
       x="Mean Expression", y="Coefficient of Variation (CV)") +
  theme_bw() +
  scale_x_log10() + scale_y_log10()
all_plots_for_pdf[["cahn_brennecke"]] <- p_mean_cv_cahn_brennecke_filtered

p_mean_cv_bewick_brennecke_filtered <- ggplot(filtered_biol, aes(x=mean_expr, y=cv_expr, color=bewick_group)) +
  geom_point(alpha=0.5, size=0.7) +
  labs(title="Mean vs CV by Bewick Group (Brennecke Filter Only)",
       x="Mean Expression", y="Coefficient of Variation (CV)") +
  theme_bw() +
  scale_x_log10() + scale_y_log10()
all_plots_for_pdf[["bewick_brennecke"]] <- p_mean_cv_bewick_brennecke_filtered


# Scenario 3: James Lloyd Filter Only (Log Scale)
p_mean_cv_cahn_james_lloyd_filtered <- ggplot(plot_results_mean_filtered, aes(x=mean_expr, y=cv_expr, color=cahn_group)) +
  geom_point(alpha=0.5, size=0.7) +
  labs(title="Mean vs CV by Cahn Group (James Lloyd Filter Only)",
       x="Mean Expression", y="Coefficient of Variation (CV)") +
  theme_bw() +
  scale_x_log10() + scale_y_log10()
all_plots_for_pdf[["cahn_james_lloyd"]] <- p_mean_cv_cahn_james_lloyd_filtered

p_mean_cv_bewick_james_lloyd_filtered <- ggplot(plot_results_mean_filtered, aes(x=mean_expr, y=cv_expr, color=bewick_group)) +
  geom_point(alpha=0.5, size=0.7) +
  labs(title="Mean vs CV by Bewick Group (James Lloyd Filter Only)",
       x="Mean Expression", y="Coefficient of Variation (CV)") +
  theme_bw() +
  scale_x_log10() + scale_y_log10()
all_plots_for_pdf[["bewick_james_lloyd"]] <- p_mean_cv_bewick_james_lloyd_filtered


# Scenario 4: James Lloyd Filter THEN Brennecke Filter (Log Scale)
if (nrow(filtered_biol_mean_filtered) > 0) {
  p_mean_cv_cahn_james_lloyd_then_brennecke_filtered <- ggplot(filtered_biol_mean_filtered, aes(x=mean_expr, y=cv_expr, color=cahn_group)) +
    geom_point(alpha=0.5, size=0.7) +
    labs(title="Mean vs CV by Cahn Group (James Lloyd Filter THEN Brennecke)",
         x="Mean Expression", y="Coefficient of Variation (CV)") +
    theme_bw() +
    scale_x_log10() + scale_y_log10()
  all_plots_for_pdf[["cahn_james_lloyd_then_brennecke"]] <- p_mean_cv_cahn_james_lloyd_then_brennecke_filtered
} else {
  warning("Skipping Cahn James Lloyd THEN Brennecke plot: filtered_biol_mean_filtered is empty.")
}

if (nrow(filtered_biol_mean_filtered) > 0) {
  p_mean_cv_bewick_james_lloyd_then_brennecke_filtered <- ggplot(filtered_biol_mean_filtered, aes(x=mean_expr, y=cv_expr, color=bewick_group)) +
    geom_point(alpha=0.5, size=0.7) +
    labs(title="Mean vs CV by Bewick Group (James Lloyd Filter THEN Brennecke)",
         x="Mean Expression", y="Coefficient of Variation (CV)") +
    theme_bw() +
    scale_x_log10() + scale_y_log10()
  all_plots_for_pdf[["bewick_james_lloyd_then_brennecke"]] <- p_mean_cv_bewick_james_lloyd_then_brennecke_filtered
} else {
  warning("Skipping Bewick James Lloyd THEN Brennecke plot: filtered_biol_mean_filtered is empty.")
}

# --- NEW: Mean Expression vs CV plots (NO Log Scale) ---

# James Lloyd Filter Only (No Log Scale)
p_mean_cv_cahn_james_lloyd_nolog <- ggplot(plot_results_mean_filtered, aes(x=mean_expr, y=cv_expr, color=cahn_group)) +
  geom_point(alpha=0.5, size=0.7) +
  labs(title="Mean vs CV by Cahn Group (James Lloyd Filter Only, No Log Scale)",
       x="Mean Expression", y="Coefficient of Variation (CV)") +
  theme_bw()
all_plots_for_pdf[["cahn_james_lloyd_nolog"]] <- p_mean_cv_cahn_james_lloyd_nolog

p_mean_cv_bewick_james_lloyd_nolog <- ggplot(plot_results_mean_filtered, aes(x=mean_expr, y=cv_expr, color=bewick_group)) +
  geom_point(alpha=0.5, size=0.7) +
  labs(title="Mean vs CV by Bewick Group (James Lloyd Filter Only, No Log Scale)",
       x="Mean Expression", y="Coefficient of Variation (CV)") +
  theme_bw()
all_plots_for_pdf[["bewick_james_lloyd_nolog"]] <- p_mean_cv_bewick_james_lloyd_nolog

# Brennecke Filter Only (No Log Scale)
p_mean_cv_cahn_brennecke_nolog <- ggplot(filtered_biol, aes(x=mean_expr, y=cv_expr, color=cahn_group)) +
  geom_point(alpha=0.5, size=0.7) +
  labs(title="Mean vs CV by Cahn Group (Brennecke Filter Only, No Log Scale)",
       x="Mean Expression", y="Coefficient of Variation (CV)") +
  theme_bw()
all_plots_for_pdf[["cahn_brennecke_nolog"]] <- p_mean_cv_cahn_brennecke_nolog

p_mean_cv_bewick_brennecke_nolog <- ggplot(filtered_biol, aes(x=mean_expr, y=cv_expr, color=bewick_group)) +
  geom_point(alpha=0.5, size=0.7) +
  labs(title="Mean vs CV by Bewick Group (Brennecke Filter Only, No Log Scale)",
       x="Mean Expression", y="Coefficient of Variation (CV)") +
  theme_bw()
all_plots_for_pdf[["bewick_brennecke_nolog"]] <- p_mean_cv_bewick_brennecke_nolog


# --- STATISTICAL ANALYSIS FOR FILTERING METHODS AND GRID PLOTS ---

# Initialize a list to store all statistical results for grid plots
all_grid_stats_df_rows <- list()
stat_counter <- 0

# Function to perform Kruskal-Wallis and post-hoc Wilcoxon tests
perform_group_stats <- function(data_df, value_col, group_col_name, context_name, filter_tag, current_ct_ln) {
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
  if (nrow(d_test) == 0 || !value_col %in% colnames(d_test) || !test_group_col %in% colnames(d_test)) {
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

  stat_counter <<- stat_counter + 1
  result_entry <- data.frame(
    ID = stat_counter,
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
  all_grid_stats_df_rows[[length(all_grid_stats_df_rows) + 1]] <<- result_entry

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
        stat_counter <<- stat_counter + 1
        pairwise_entry <- data.frame(
          ID = stat_counter,
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
        all_grid_stats_df_rows[[length(all_grid_stats_df_rows) + 1]] <<- pairwise_entry
      }
    }
  }
}

# Iterate and collect stats for Mean Expression by Celltype/Lineage
for (ct in celltype_levels) { # Use celltype_levels
  if (is.na(ct) || ct == "") next
  d_unfiltered <- plot_results
  d_unfiltered$cell_mean <- gene_means_celltype[d_unfiltered$gene, ct]
  perform_group_stats(d_unfiltered, "cell_mean", "cahn_group", "Mean Expr by Celltype", "Unfiltered", ct)
  perform_group_stats(d_unfiltered, "cell_mean", "bewick_group", "Mean Expr by Celltype", "Unfiltered", ct)

  d_filtered <- filtered_biol
  d_filtered$cell_mean <- gene_means_celltype[d_filtered$gene, ct]
  perform_group_stats(d_filtered, "cell_mean", "cahn_group", "Mean Expr by Celltype", "Brennecke Filtered", ct)
  perform_group_stats(d_filtered, "cell_mean", "bewick_group", "Mean Expr by Celltype", "Brennecke Filtered", ct)

  d_jl_filtered <- plot_results_mean_filtered # Data filtered by James Lloyd
  d_jl_filtered$cell_mean <- gene_means_celltype[d_jl_filtered$gene, ct]
  perform_group_stats(d_jl_filtered, "cell_mean", "cahn_group", "Mean Expr by Celltype", "James Lloyd Filtered", ct)
  perform_group_stats(d_jl_filtered, "cell_mean", "bewick_group", "Mean Expr by Celltype", "James Lloyd Filtered", ct)
}

for (ln in lineage_levels) { # Use lineage_levels
  if (is.na(ln) || ln == "") next
  d_unfiltered <- plot_results
  d_unfiltered$lineage_mean <- gene_means_lineage[d_unfiltered$gene, ln]
  perform_group_stats(d_unfiltered, "lineage_mean", "cahn_group", "Mean Expr by Lineage", "Unfiltered", ln)
  perform_group_stats(d_unfiltered, "lineage_mean", "bewick_group", "Mean Expr by Lineage", "Unfiltered", ln)

  d_filtered <- filtered_biol
  d_filtered$lineage_mean <- gene_means_lineage[d_filtered$gene, ln]
  perform_group_stats(d_filtered, "lineage_mean", "cahn_group", "Mean Expr by Lineage", "Brennecke Filtered", ln)
  perform_group_stats(d_filtered, "lineage_mean", "bewick_group", "Mean Expr by Lineage", "Brennecke Filtered", ln)

  d_jl_filtered <- plot_results_mean_filtered # Data filtered by James Lloyd
  d_jl_filtered$lineage_mean <- gene_means_lineage[d_jl_filtered$gene, ln]
  perform_group_stats(d_jl_filtered, "lineage_mean", "cahn_group", "Mean Expr by Lineage", "James Lloyd Filtered", ln)
  perform_group_stats(d_jl_filtered, "lineage_mean", "bewick_group", "Mean Expr by Lineage", "James Lloyd Filtered", ln)
}

# Iterate and collect stats for CV Expression by Celltype/Lineage
for (ct in celltype_levels) { # Use celltype_levels
  if (is.na(ct) || ct == "") next
  # For CV, use the celltype-specific data frame
  d_ct_unfiltered <- gene_data_by_celltype[[ct]]
  if (!is.null(d_ct_unfiltered) && nrow(d_ct_unfiltered) > 0) {
    perform_group_stats(d_ct_unfiltered, "cv_expr_ct", "cahn_group", "CV Expr by Celltype", "Unfiltered", ct)
    perform_group_stats(d_ct_unfiltered, "cv_expr_ct", "bewick_group", "CV Expr by Celltype", "Unfiltered", ct)
  }

  # Brennecke Filtered CV
  d_ct_brennecke_filtered <- d_ct_unfiltered %>% filter(gene %in% filtered_biol$gene)
  if (!is.null(d_ct_brennecke_filtered) && nrow(d_ct_brennecke_filtered) > 0) {
    perform_group_stats(d_ct_brennecke_filtered, "cv_expr_ct", "cahn_group", "CV Expr by Celltype", "Brennecke Filtered", ct)
    perform_group_stats(d_ct_brennecke_filtered, "cv_expr_ct", "bewick_group", "CV Expr by Celltype", "Brennecke Filtered", ct)
  }

  # James Lloyd Filtered CV
  d_ct_james_lloyd_filtered <- d_ct_unfiltered %>% filter(gene %in% plot_results_mean_filtered$gene)
  if (!is.null(d_ct_james_lloyd_filtered) && nrow(d_ct_james_lloyd_filtered) > 0) {
    perform_group_stats(d_ct_james_lloyd_filtered, "cv_expr_ct", "cahn_group", "CV Expr by Celltype", "James Lloyd Filtered", ct)
    perform_group_stats(d_ct_james_lloyd_filtered, "cv_expr_ct", "bewick_group", "CV Expr by Celltype", "James Lloyd Filtered", ct)
  }
}

for (ln in lineage_levels) { # Use lineage_levels
  if (is.na(ln) || ln == "") next
  # For CV, use the lineage-specific data frame
  d_ln_unfiltered <- gene_data_by_lineage[[ln]]
  if (!is.null(d_ln_unfiltered) && nrow(d_ln_unfiltered) > 0) {
    perform_group_stats(d_ln_unfiltered, "cv_expr_ln", "cahn_group", "CV Expr by Lineage", "Unfiltered", ln)
    perform_group_stats(d_ln_unfiltered, "cv_expr_ln", "bewick_group", "CV Expr by Lineage", "Unfiltered", ln)
  }

  # Brennecke Filtered CV
  d_ln_brennecke_filtered <- d_ln_unfiltered %>% filter(gene %in% filtered_biol$gene)
  if (!is.null(d_ln_brennecke_filtered) && nrow(d_ln_brennecke_filtered) > 0) {
    perform_group_stats(d_ln_brennecke_filtered, "cv_expr_ln", "cahn_group", "CV Expr by Lineage", "Brennecke Filtered", ln)
    perform_group_stats(d_ln_brennecke_filtered, "cv_expr_ln", "bewick_group", "CV Expr by Lineage", "Brennecke Filtered", ln)
  }

  # James Lloyd Filtered CV
  d_ln_james_lloyd_filtered <- d_ln_unfiltered %>% filter(gene %in% plot_results_mean_filtered$gene)
  if (!is.null(d_ln_james_lloyd_filtered) && nrow(d_ln_james_lloyd_filtered) > 0) {
    perform_group_stats(d_ln_james_lloyd_filtered, "cv_expr_ln", "cahn_group", "CV Expr by Lineage", "James Lloyd Filtered", ln)
    perform_group_stats(d_ln_james_lloyd_filtered, "cv_expr_ln", "bewick_group", "CV Expr by Lineage", "James Lloyd Filtered", ln)
  }
}


# Combine all collected stats into a single data frame
final_grid_stats_df <- if (length(all_grid_stats_df_rows) > 0) {
  do.call(rbind, all_grid_stats_df_rows)
} else {
  data.frame(ID=numeric(0), Context=character(0), Filter_Tag=character(0), Group_Type=character(0),
             Value_Type=character(0), Specific_Context=character(0), Test=character(0),
             Comparison=character(0), P_value=character(0), Statistic=numeric(0))
}


# Now for the filter comparison plot statistics.
# Summarize the counts and median values for each filter category
filter_category_summary <- plot_results %>%
  group_by(filter_category) %>%
  summarise(
    Count = n(),
    Median_CV = median(cv_expr, na.rm = TRUE),
    Median_CV2 = median(cv2, na.rm = TRUE),
    Median_Mean_Expr = median(mean_expr, na.rm = TRUE)
  ) %>%
  mutate(
    Median_CV = formatC(Median_CV, format="e", digits=2),
    Median_CV2 = formatC(Median_CV2, format="e", digits=2),
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

# Table 3: Grid Plot Statistics (Paginated if necessary)
if (nrow(final_grid_stats_df) > 0) {
  rows_per_page <- 40 # Adjust as needed for readability
  num_pages <- ceiling(nrow(final_grid_stats_df) / rows_per_page)

  for (p in 1:num_pages) {
    start_row <- (p - 1) * rows_per_page + 1
    end_row <- min(p * rows_per_page, nrow(final_grid_stats_df))
    page_data <- final_grid_stats_df[start_row:end_row, ]

    grid_stats_grob <- tableGrob(page_data, rows = NULL, theme = ttheme_default(base_size = 8))
    grid_stats_title <- textGrob(paste0("Grid Plot Statistics (Page ", p, " of ", num_pages, ")"), gp = gpar(fontsize = 14, fontface = "bold"))
    all_tables_for_pdf[[length(all_tables_for_pdf) + 1]] <- arrangeGrob(grid_stats_title, grid_stats_grob, ncol = 1, heights = c(0.1, 0.9))
  }
} else {
  warning("No grid plot statistics to add to PDF.")
}


# --- Save all plots and tables to a single PDF ---
pdf_filename <- file.path(output_dir, "mean_cv_filtering_comparison.pdf")
pdf(pdf_filename, width=10, height=8) # Set dimensions for plots and tables

# Print all plots
for (plot_obj in all_plots_for_pdf) {
  if (!is.null(plot_obj)) {
    print(plot_obj)
  }
}

# Print all tables
for (table_grob_obj in all_tables_for_pdf) {
  if (!is.null(table_grob_obj)) {
    grid.newpage() # Start a new page for each table
    grid.draw(table_grob_obj)
  }
}

dev.off()

message(paste("All plots and statistics saved to:", pdf_filename))

# --- NEW: Save James Lloyd Filtered CV Boxplots to separate PDFs ---

# Boxplots by celltype (Cahn group) - Expression Noise (CV) (James Lloyd Filtered)
plots_cahn_celltype_cv_james_lloyd_filtered <- lapply(celltype_levels, function(ct) {
  d_ct <- gene_data_by_celltype[[ct]] # Start with celltype-specific data
  if (is.null(d_ct) || nrow(d_ct) == 0) return(NULL)

  # Filter to include only genes that were overall James Lloyd-filtered
  d_filtered_jl_ct <- d_ct %>% filter(gene %in% plot_results_mean_filtered$gene)
  if (nrow(d_filtered_jl_ct) == 0) return(NULL)

  ggplot(d_filtered_jl_ct[!is.na(d_filtered_jl_ct$cahn_group),], aes(x=factor(cahn_group), y=cv_expr_ct, fill=factor(cahn_group))) +
    geom_boxplot(outlier.size=0.5) +
    labs(title=paste("Cahn Group in", ct), x="Cahn Group", y="Expression Noise (CV) (James Lloyd Filtered)") +
    theme_bw() +
    theme(legend.position = "none", axis.text.x = element_text(angle=45, hjust=1))
})
g_cahn_celltype_cv_james_lloyd_filtered <- marrangeGrob(plots_cahn_celltype_cv_james_lloyd_filtered, nrow=2, ncol=3, top="Expression Noise (CV) by Cahn Group in Celltypes (James Lloyd Filtered)")
pdf(file.path(output_dir, "grid_cv_expr_cahn_by_celltype_james_lloyd_filtered.pdf"), width=18, height=12)
for (i in 1:length(g_cahn_celltype_cv_james_lloyd_filtered)) {
  grid.newpage()
  grid.draw(g_cahn_celltype_cv_james_lloyd_filtered[[i]])
}
dev.off()

# Boxplots by celltype (Bewick group) - Expression Noise (CV) (James Lloyd Filtered)
plots_bewick_celltype_cv_james_lloyd_filtered <- lapply(celltype_levels, function(ct) {
  d_ct <- gene_data_by_celltype[[ct]] # Start with celltype-specific data
  if (is.null(d_ct) || nrow(d_ct) == 0) return(NULL)

  # Filter to include only genes that were overall James Lloyd-filtered
  d_filtered_jl_ct <- d_ct %>% filter(gene %in% plot_results_mean_filtered$gene)
  if (nrow(d_filtered_jl_ct) == 0) return(NULL)

  ggplot(d_filtered_jl_ct[!is.na(d_filtered_jl_ct$bewick_group),], aes(x=factor(bewick_group), y=cv_expr_ct, fill=factor(bewick_group))) +
    geom_boxplot(outlier.size=0.5) +
    labs(title=paste("Bewick Group in", ct), x="Bewick Group", y="Expression Noise (CV) (James Lloyd Filtered)") +
    theme_bw() +
    theme(legend.position = "none", axis.text.x = element_text(angle=45, hjust=1))
})
g_bewick_celltype_cv_james_lloyd_filtered <- marrangeGrob(plots_bewick_celltype_cv_james_lloyd_filtered, nrow=2, ncol=3, top="Expression Noise (CV) by Bewick Group in Celltypes (James Lloyd Filtered)")
pdf(file.path(output_dir, "grid_cv_expr_bewick_by_celltype_james_lloyd_filtered.pdf"), width=18, height=12)
for (i in 1:length(g_bewick_celltype_cv_james_lloyd_filtered)) {
  grid.newpage()
  grid.draw(g_bewick_celltype_cv_james_lloyd_filtered[[i]])
}
dev.off()

# Boxplots by lineage (Cahn group) - Expression Noise (CV) (James Lloyd Filtered)
plots_cahn_lineage_cv_james_lloyd_filtered <- lapply(lineage_levels, function(ln) {
  d_ln <- gene_data_by_lineage[[ln]] # Start with lineage-specific data
  if (is.null(d_ln) || nrow(d_ln) == 0) return(NULL)

  # Filter to include only genes that were overall James Lloyd-filtered
  d_filtered_jl_ln <- d_ln %>% filter(gene %in% plot_results_mean_filtered$gene)
  if (nrow(d_filtered_jl_ln) == 0) return(NULL)

  ggplot(d_filtered_jl_ln[!is.na(d_filtered_jl_ln$cahn_group),], aes(x=factor(cahn_group), y=cv_expr_ln, fill=factor(cahn_group))) +
    geom_boxplot(outlier.size=0.5) +
    labs(title=paste("Cahn Group in", ln), x="Cahn Group", y="Expression Noise (CV) (James Lloyd Filtered)") +
    theme_bw() +
    theme(legend.position = "none", axis.text.x = element_text(angle=45, hjust=1))
})
g_cahn_lineage_cv_james_lloyd_filtered <- marrangeGrob(plots_cahn_lineage_cv_james_lloyd_filtered, nrow=2, ncol=3, top="Expression Noise (CV) by Cahn Group in Lineages (James Lloyd Filtered)")
pdf(file.path(output_dir, "grid_cv_expr_cahn_by_lineage_james_lloyd_filtered.pdf"), width=18, height=12)
for (i in 1:length(g_cahn_lineage_cv_james_lloyd_filtered)) {
  grid.newpage()
  grid.draw(g_cahn_lineage_cv_james_lloyd_filtered[[i]])
}
dev.off()

# Boxplots by lineage (Bewick group) - Expression Noise (CV) (James Lloyd Filtered)
plots_bewick_lineage_cv_james_lloyd_filtered <- lapply(lineage_levels, function(ln) {
  d_ln <- gene_data_by_lineage[[ln]] # Start with lineage-specific data
  if (is.null(d_ln) || nrow(d_ln) == 0) return(NULL)

  # Filter to include only genes that were overall James Lloyd-filtered
  d_filtered_jl_ln <- d_ln %>% filter(gene %in% plot_results_mean_filtered$gene)
  if (nrow(d_filtered_jl_ln) == 0) return(NULL)

  ggplot(d_filtered_jl_ln[!is.na(d_filtered_jl_ln$bewick_group),], aes(x=factor(bewick_group), y=cv_expr_ln, fill=factor(bewick_group))) +
    geom_boxplot(outlier.size=0.5) +
    labs(title=paste("Bewick Group in", ln), x="Bewick Group", y="Expression Noise (CV) (James Lloyd Filtered)") +
    theme_bw() +
    theme(legend.position = "none", axis.text.x = element_text(angle=45, hjust=1))
})
g_bewick_lineage_cv_james_lloyd_filtered <- marrangeGrob(plots_bewick_lineage_cv_james_lloyd_filtered, nrow=2, ncol=3, top="Expression Noise (CV) by Bewick Group in Lineages (James Lloyd Filtered)")
pdf(file.path(output_dir, "grid_cv_expr_bewick_by_lineage_james_lloyd_filtered.pdf"), width=18, height=12)
for (i in 1:length(g_bewick_lineage_cv_james_lloyd_filtered)) {
  grid.newpage()
  grid.draw(g_bewick_lineage_cv_james_lloyd_filtered[[i]])
}
dev.off()