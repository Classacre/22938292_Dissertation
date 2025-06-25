# Export full Seurat metadata to CSV
library(Seurat)

seurat_path <- "/group/sms029/Oliva_dataset/integrated_col_trajectories_colonly.rds"
outfile <- "/group/sms029/mnieuwenh/seurat_metadata/seurat_metadata_full.csv"

seurat_obj <- readRDS(seurat_path)
write.csv(seurat_obj@meta.data, file=outfile, row.names=TRUE)
cat("Full metadata saved to:", outfile, "\n")
