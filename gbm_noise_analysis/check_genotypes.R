# Check genotype diversity in Seurat metadata
seurat_obj <- readRDS("/group/sms029/Oliva_dataset/integrated_col_trajectories_colonly.rds")
genotype_table <- table(seurat_obj@meta.data$genotype, useNA = "ifany")
write.csv(as.data.frame(genotype_table), "/group/sms029/mnieuwenh/gbm_noise_analysis/genotype_counts.csv", row.names=FALSE)
cat("Genotype counts saved to gbm_noise_analysis/genotype_counts.csv\n")
