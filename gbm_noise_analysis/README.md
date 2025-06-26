# Pipeline for Comparing Gene Expression Noise and Methylation in Seurat Data

This workflow analyzes gene expression noise in Arabidopsis root single-cell RNA-seq data, focusing on the effects of gene body methylation (gbM) and TE-like methylation, using a large Seurat object on an HPC cluster with SLURM.

## Files and Inputs
- **Seurat object:** `/group/sms029/Oliva_dataset/integrated_col_trajectories_colonly.rds`
- **Combined methylation annotation:** `/group/sms029/mnieuwenh/Methylation_Data/combined_methylation_data.csv` (from Cahn et al. and Bewick et al.)
- **Output directory:** `/group/sms029/mnieuwenh/gbm_noise_analysis/`

## Data Availability
- The methylation annotation includes public data from Cahn et al. ([PubMed](https://pubmed.ncbi.nlm.nih.gov/39632087/)) and Bewick et al. ([PNAS](https://www.pnas.org/doi/10.1073/pnas.1604666113)).

## Main Steps
1. **Export gene expression matrix from Seurat**
   - `export_expression_matrix.R`/`.sh`: Exports log-normalized gene expression matrix for all cells.
2. **Calculate expression noise and compare methylation groups**
   - `gbm_noise_analysis.R`/`.sh`: Calculates mean, variance, standard deviation, and coefficient of variation (CV) for each gene. Annotates each gene with methylation status from both Cahn and Bewick studies.
   - Compares noise (CV) between:
     - Cahn: gbM, TE-like methylation, Unmethylated
     - Bewick: gbM, mCHG, mCHH, Unmethylated
   - Generates boxplots and scatter plots for each group.
3. **Responsive gene analysis**
   - `responsive_gene_analysis.R`/`.sh`: Identifies genes highly expressed in one cell type but low in others (responsive genes).
4. **Genomic feature annotation**
   - `genomic_feature_analysis.R`/`.sh`: Annotates genes with noise, responsiveness, and methylation features for downstream analysis.

## How to Run
1. Export the expression matrix:
   ```sh
   sbatch export_expression_matrix.sh
   ```
2. Run the noise and methylation analysis:
   ```sh
   sbatch gbm_noise_analysis.sh
   ```
3. (Optional) Run responsive gene and feature annotation analyses:
   ```sh
   sbatch responsive_gene_analysis.sh
   sbatch genomic_feature_analysis.sh
   ```

## Output Files
- **Expression matrix:** `expression_matrix.csv`
- **Noise analysis results:** `results/high_low_noise/gbm_noise_comparison.csv`
- **Boxplots:**
  - `gbm_noise_boxplot_cahn.png` (Cahn: gbM, TE-like, Unmethylated)
  - `gbm_noise_boxplot_bewick.png` (Bewick: gbM, mCHG, mCHH, Unmethylated)
  - `gbm_noise_boxplot_combined.png` (all groups)
  - `gbm_noise_boxplot_te_like.png`, `gbm_noise_boxplot_by_celltype.png`, `gbm_noise_boxplot_by_lineage.png`, etc.
  - `gbm_noise_boxplot_h2az.png` (CV by H2A.Z group)
  - `gbm_noise_boxplot_meth_h2az.png` (CV by methylation + H2A.Z group)
- **Scatter plots:**
  - `gbm_noise_scatter_cahn.png`, `gbm_noise_scatter_bewick.png`, `gbm_noise_scatter_combined.png`
  - `mean_vs_sd_by_cahn_group.png`, `mean_vs_cv_by_cahn_group.png`
  - `mean_vs_sd_by_bewick_group.png`, `mean_vs_cv_by_bewick_group.png`
  - `mean_vs_sd_by_methylation_group.png`, `mean_vs_cv_by_methylation_group.png`
  - `mean_vs_sd_by_h2az_group.png`, `mean_vs_cv_by_h2az_group.png`
  - `mean_vs_sd_by_meth_h2az_group.png`, `mean_vs_cv_by_meth_h2az_group.png`
- **Statistical results:**
  - `gbm_noise_stats.txt` (Wilcoxon test p-values for group comparisons)
  - `gbm_noise_report.txt` (summary report)
  - `fisher_h2az.txt` (Fisher's test for H2A.Z enrichment)
- **Responsive gene results:**
  - `results/responsive_genes/responsive_genes.csv`, `responsive_noise_stats.txt`, `responsive_vs_nonresponsive_noise.png`
  - `prop_responsive_by_h2az.csv` (proportion of responsive genes by H2A.Z group)
  - `cv_by_h2az_and_responsive.png` (CV by H2A.Z group and responsiveness)
  - `cv_by_meth_h2az_and_responsive.png` (CV by methylation + H2A.Z group and responsiveness)
- **Logs:**
  - All SLURM output and error logs are in the `logs/` subfolder.

## Batch Scripts
- All batch scripts (`*.sh`) are provided for SLURM job submission and reproducibility.
- Adjust memory and CPU in the sbatch scripts if needed (default: 120G RAM, 8 CPUs, 24h for large jobs).
- All scripts assume the conda environment at `/group/sms029/conda_environment/R`.

---

## Notes
- The pipeline now supports flexible comparison of gene expression noise and responsiveness across both methylation and H2A.Z categories, including their combinations.
- All code and workflow steps are provided for transparency and future reproducibility.
- You may need to edit gene names if your Seurat object uses a different format than the methylation annotation.
