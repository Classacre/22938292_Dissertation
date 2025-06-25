# Script to remove NA genotypes from a Seurat object and save the filtered object
library(Seurat)

input_path <- "/group/sms029/Oliva_dataset/integrated_col_trajectories_colonly.rds"
output_path <- "/group/sms029/Oliva_dataset/integrated_col_trajectories_colonly.rds"

seurat_obj <- readRDS(input_path)

# Subset to only cells with genotype == 'col' (removes NA and all others)
seurat_obj_col <- subset(seurat_obj, subset = genotype == "col")

saveRDS(seurat_obj_col, output_path)
cat("Filtered Seurat object saved to:", output_path, "\n")
