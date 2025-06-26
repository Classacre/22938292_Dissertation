# Cell count per cell type and lineage from Seurat object
library(Seurat)

seurat_path <- "/group/sms029/Oliva_dataset/integrated_col_trajectories_colonly.rds"
outfile_celltype <- "/group/sms029/mnieuwenh/seurat_utilities/cell_counts_by_celltype.csv"
outfile_lineage <- "/group/sms029/mnieuwenh/seurat_utilities/cell_counts_by_lineage.csv"

seurat_obj <- readRDS(seurat_path)

# Cell type counts
if ("identity" %in% colnames(seurat_obj@meta.data)) {
  celltype_counts <- as.data.frame(table(seurat_obj@meta.data$identity))
  colnames(celltype_counts) <- c("celltype", "cell_count")
  write.csv(celltype_counts, outfile_celltype, row.names=FALSE)
  cat("Cell type counts saved to:", outfile_celltype, "\n")
} else {
  cat("No 'identity' column found in Seurat metadata.\n")
}

# Lineage counts
if ("lineage" %in% colnames(seurat_obj@meta.data)) {
  lineage_counts <- as.data.frame(table(seurat_obj@meta.data$lineage))
  colnames(lineage_counts) <- c("lineage", "cell_count")
  write.csv(lineage_counts, outfile_lineage, row.names=FALSE)
  cat("Lineage counts saved to:", outfile_lineage, "\n")
} else {
  cat("No 'lineage' column found in Seurat metadata.\n")
}
