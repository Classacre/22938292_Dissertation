# Check Seurat object RNA assay structure and available slots
library(Seurat)

seurat_path <- "/group/sms029/Oliva_dataset/integrated_col_trajectories.rds"

seurat_obj <- readRDS(seurat_path)

cat("Assay class:", class(seurat_obj@assays$RNA), "\n\n")
cat("Assay slots:\n")
print(slotNames(seurat_obj@assays$RNA))
cat("\nAssay names:\n")
print(names(seurat_obj@assays$RNA))
cat("\nAssay structure:\n")
str(seurat_obj@assays$RNA)
cat("\nAvailable layers (if any):\n")
if ("layers" %in% slotNames(seurat_obj@assays$RNA)) {
  print(names(seurat_obj@assays$RNA@layers))
}
cat("\nSeurat object version:", seurat_obj@version, "\n")
cat("Matrix check complete.\n")
