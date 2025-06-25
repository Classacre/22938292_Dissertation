# Load required libraries
library(Seurat)
library(ggplot2)

# Path to Seurat object
seurat_path <- "/group/sms029/Oliva_dataset/integrated_col_trajectories_colonly.rds"

# Load the Seurat object
seurat_obj <- readRDS(seurat_path)

# Subset to only 'col' genotype (exclude NA and others)
seurat_obj <- subset(seurat_obj, subset = genotype == "col")

# Set output directory for plots
plot_dir <- "/group/sms029/mnieuwenh/seurat_plots"
if (!dir.exists(plot_dir)) dir.create(plot_dir)

caption_text <- "Only the Col genotype was used for this analysis."

# UMAP plot colored by cluster identity
if ("umap.rpca" %in% names(seurat_obj@reductions)) {
  p1 <- DimPlot(seurat_obj, reduction = "umap.rpca", group.by = "seurat_clusters") +
    ggtitle("UMAP by Seurat Clusters") +
    labs(caption = caption_text)
  ggsave(filename = file.path(plot_dir, "umap_seurat_clusters.png"), plot = p1, width = 8, height = 6)
}

# UMAP plot colored by cell type (if 'identity' or similar exists)
if ("identity" %in% colnames(seurat_obj@meta.data)) {
  p2 <- DimPlot(seurat_obj, reduction = "umap.rpca", group.by = "identity") +
    ggtitle("UMAP by Identity") +
    labs(caption = caption_text)
  ggsave(filename = file.path(plot_dir, "umap_identity.png"), plot = p2, width = 8, height = 6)
}

# UMAP plot colored by genotype (if 'genotype' exists)
if ("genotype" %in% colnames(seurat_obj@meta.data)) {
  p3 <- DimPlot(seurat_obj, reduction = "umap.rpca", group.by = "genotype") +
    ggtitle("UMAP by Genotype") +
    labs(caption = caption_text)
  ggsave(filename = file.path(plot_dir, "umap_genotype.png"), plot = p3, width = 8, height = 6)
}

# Feature plot for a marker gene (replace 'ACT2' with a gene of interest)
marker_gene <- "ACT2"
if (marker_gene %in% rownames(seurat_obj)) {
  p4 <- FeaturePlot(seurat_obj, features = marker_gene, reduction = "umap.rpca") +
    ggtitle(paste("FeaturePlot:", marker_gene)) +
    labs(caption = caption_text)
  ggsave(filename = file.path(plot_dir, paste0("featureplot_", marker_gene, ".png")), plot = p4, width = 8, height = 6)
}

# Violin plots for QC metrics (now using cell type/identity)
qc_features <- c("nFeature_RNA", "nCount_RNA")
for (feature in qc_features) {
  if (feature %in% colnames(seurat_obj@meta.data)) {
    p <- VlnPlot(seurat_obj, features = feature, group.by = "identity") +
      ggtitle(paste("Violin plot of", feature, "by Cell Type (Identity)")) +
      labs(caption = caption_text) +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
        plot.margin = margin(10, 10, 10, 10)
      )
    ggsave(filename = file.path(plot_dir, paste0("violin_", feature, ".png")), plot = p, width = 16, height = 6)
  }
}

# Dot plot for marker genes (replace with relevant genes)
marker_genes <- c("ACT2", "SCR", "WOX5")
marker_genes_present <- marker_genes[marker_genes %in% rownames(seurat_obj)]
if (length(marker_genes_present) > 0) {
  p_dot <- DotPlot(seurat_obj, features = marker_genes_present, group.by = "identity") + RotatedAxis() +
    ggtitle("DotPlot of Marker Genes by Identity") +
    labs(caption = caption_text)
  ggsave(filename = file.path(plot_dir, "dotplot_marker_genes.png"), plot = p_dot, width = 10, height = 6)
}

# Heatmap for top 10 variable genes
if (length(seurat_obj@assays$RNA@var.features) >= 10) {
  top10 <- seurat_obj@assays$RNA@var.features[1:10]
  p_heat <- DoHeatmap(seurat_obj, features = top10, group.by = "seurat_clusters") +
    ggtitle("Heatmap of Top 10 Variable Genes") +
    labs(caption = caption_text)
  ggsave(filename = file.path(plot_dir, "heatmap_top10_variable_genes.png"), plot = p_heat, width = 10, height = 8)
}

# Ridge plot for a marker gene (replace 'ACT2' if desired)
if ("ACT2" %in% rownames(seurat_obj)) {
  p_ridge <- RidgePlot(seurat_obj, features = "ACT2", group.by = "seurat_clusters") +
    ggtitle("Ridge Plot: ACT2 by Seurat Clusters") +
    labs(caption = caption_text)
  ggsave(filename = file.path(plot_dir, "ridgeplot_ACT2.png"), plot = p_ridge, width = 10, height = 6)
}

# UMAP colored by pseudotime (if available)
if ("pseudotime" %in% colnames(seurat_obj@meta.data)) {
  seurat_obj$pseudotime_factor <- as.numeric(seurat_obj@meta.data$pseudotime)
  p_pt <- FeaturePlot(seurat_obj, features = "pseudotime_factor", reduction = "umap.rpca") +
    ggtitle("UMAP colored by Pseudotime") +
    labs(caption = caption_text)
  ggsave(filename = file.path(plot_dir, "umap_pseudotime.png"), plot = p_pt, width = 8, height = 6)
}

cat("Plots saved to:", plot_dir, "\n")