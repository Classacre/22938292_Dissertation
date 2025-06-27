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
- **Violin plots:**
  - `gbm_noise_violin_cahn.png`, `gbm_noise_violin_bewick.png` (CV by methylation group, violin)
  - `results/responsive_genes/responsive_vs_nonresponsive_violin.png` (CV for responsive vs non-responsive genes, violin)
- **Ridge plots:**
  - `gbm_noise_ridge_cahn.png`, `gbm_noise_ridge_bewick.png` (CV by methylation group, ridge)
  - `results/responsive_genes/responsive_vs_nonresponsive_ridge.png` (CV for responsive vs non-responsive genes, ridge)
- **UpSet plots:**
  - `upset_highnoise_gbm_h2az.png` (Overlap of high-noise, gbM, and H2A.Z-depleted genes)
  - `results/responsive_genes/upset_responsive_gbm_h2az.png` (Overlap of responsive, gbM, and H2A.Z-depleted genes)
- **Heatmaps:**
  - `heatmap_top100_highnoise.png` (Top 100 high-noise genes, wide canvas)
  - `results/responsive_genes/heatmap_top100_responsive.png` (Top 100 responsive genes, wide canvas)
- **Statistical results:**
  - `gbm_noise_stats_summary.csv` (table of p-values for all group comparisons)
  - `gbm_noise_stats_summary.txt` (plain-language summary of findings)
- **Responsive gene results:**
  - `results/responsive_genes/responsive_stats_summary.csv` (table of p-values for all group comparisons)
  - `results/responsive_genes/responsive_stats_summary.txt` (plain-language summary of findings)
- **Logs:**
  - All SLURM output and error logs are in the `logs/` subfolder.

## Batch Scripts
- All batch scripts (`*.sh`) are provided for SLURM job submission and reproducibility.
- Adjust memory and CPU in the sbatch scripts if needed (default: 120G RAM, 8 CPUs, 24h for large jobs).
- All scripts assume the conda environment at `/group/sms029/conda_environment/R`.

---

## Notes
- All statistical results are now summarized in a single table and a single summary text file per analysis folder for clarity. Old .txt files with individual p-values have been removed for a cleaner output directory.
- The pipeline now supports flexible comparison of gene expression noise and responsiveness across both methylation and H2A.Z categories, including their combinations.
- All code and workflow steps are provided for transparency and future reproducibility.
- You may need to edit gene names if your Seurat object uses a different format than the methylation annotation.

## Noise Filtering Approach (Brennecke et al., 2013)
To distinguish technical noise from biological variability, we use the method of Brennecke et al. (2013), which models the expected technical noise as a function of mean expression:

\[
\text{CV}^2 = a_0 + \frac{a_1}{\mu}
\]

where:
- \(\text{CV}^2\) is the squared coefficient of variation for each gene,
- \(\mu\) is the mean expression of the gene,
- \(a_0\) and \(a_1\) are parameters estimated from the data by fitting the curve to the majority of genes (assumed to be dominated by technical noise).

Genes with observed \(\text{CV}^2\) significantly above this fitted curve are considered to have biological variability (not just technical noise). Typically, genes above the 95% confidence interval of the fit are retained for downstream noise analysis.

This approach ensures that only genes with evidence of true biological variability are included in noise analyses, reducing the impact of technical artifacts from lowly expressed genes.

- See `mean_vs_cv2_brennecke.png` for a visualization of the fitted curve and selected genes.

# In-depth Explanation: Noise Filtering and Its Importance

## Why Filter for Technical Noise?
Single-cell RNA-seq data is inherently noisy, especially for genes with low expression. This noise can arise from technical limitations (such as inefficient RNA capture or amplification) and does not reflect true biological differences. If not accounted for, technical noise can obscure real biological patterns and lead to misleading conclusions.

## Step-by-Step: How We Distinguish Technical Noise from Biological Variability

1. **Calculate Mean and Variance for Each Gene**
   - For every gene, we compute the mean expression (average across all cells) and the variance (how much the expression varies between cells).
   - These are basic statistical measures: mean tells us the typical expression, variance tells us how much it fluctuates.

2. **Compute the Coefficient of Variation (CV)**
   - The CV is the standard deviation divided by the mean. It tells us how variable a gene is, relative to its average expression.
   - We use the squared CV (CV²) for mathematical convenience and to match the Brennecke method.

3. **Plot Mean Expression vs. CV² for All Genes**
   - Most genes, especially those with low mean expression, show high CV² due to technical noise.
   - By plotting all genes, we can see the general trend: as mean expression increases, technical noise (CV²) decreases.

4. **Fit the Brennecke et al. (2013) Model**
   - The Brennecke model describes the expected technical noise as:
     
     \[
     \text{CV}^2 = a_0 + \frac{a_1}{\mu}
     \]
     - Here, \(\mu\) is the mean expression, and \(a_0\), \(a_1\) are fitted parameters.
     - This formula captures the fact that technical noise is higher for lowly expressed genes.

5. **Identify Genes with True Biological Variability**
   - Genes whose observed CV² is significantly above the fitted curve are likely to be truly variable due to biology, not just technical error.
   - We typically use a statistical threshold (e.g., 95% confidence interval) to select these genes.

6. **Filter Out Genes Dominated by Technical Noise**
   - Only genes above the threshold are kept for downstream analysis. This ensures our results reflect real biological differences, not artifacts.

## Why This Matters
- **Accurate Results:** By removing genes dominated by technical noise, we avoid false positives and focus on genes that are genuinely variable between cells.
- **Reproducibility:** This method is widely used and accepted in the single-cell field, making our results robust and comparable to other studies.
- **Transparency:** The approach is visualized in `mean_vs_cv2_brennecke.png`, so anyone can see how the threshold was set and which genes were selected.

## Mathematical Significance (in Simple Terms)
- The formula \(\text{CV}^2 = a_0 + \frac{a_1}{\mu}\) means that for genes with low average expression (small \(\mu\)), the expected noise is high. As the average expression increases, the expected noise drops.
- By fitting this curve to the data, we can tell which genes are more variable than expected by chance (technical noise), and which are likely to be biologically interesting.

## For Non-Experts
- Think of technical noise like static on a radio: it’s louder when the signal is weak (low expression), but fades as the signal gets stronger (high expression).
- Our method finds the point where the “static” is no longer the main thing we hear, so we can focus on the real “music” (biological signal).

---
