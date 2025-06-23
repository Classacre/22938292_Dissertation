# Pipeline for Comparing GBM Gene Expression Noise in Seurat Data

This workflow compares expression noise between gene-body methylated (GBM) and non-GBM genes using a large Seurat object on an HPC cluster with SLURM.

## Files and Inputs
- **Seurat object:** `/group/sms029/Oliva_dataset/integrated_col_trajectories.rds`
- **GBM gene list:** `/group/sms029/mnieuwenh/gbM_data/unique_gbm_genes.csv` (column: `Gene`)
- **Output directory:** `/group/sms029/mnieuwenh/gbm_noise_analysis/`

## Data Availability
- The GBM gene list is not included in this repository as it is the property of Cahn et al ([PubMed](https://pubmed.ncbi.nlm.nih.gov/39632087/)) and Bewick et al ([PNAS](https://www.pnas.org/doi/10.1073/pnas.1604666113)).

## Steps
1. **Export gene expression matrix from Seurat (R script)**
2. **Calculate expression noise and compare GBM vs. non-GBM (R script)**
3. **Run each step with SLURM using provided bash scripts**

## How to Run
1. Submit the export job:
   ```sh
   sbatch export_expression_matrix.sh
   ```
2. When that finishes, submit the analysis job:
   ```sh
   sbatch gbm_noise_analysis.sh
   ```

## Additional Analysis Steps

### 4. Check Genotype Diversity
To check which genotypes are present in your Seurat object, run:
```sh
sbatch batch_scripts/check_genotypes.sh
```
This will output a CSV file `gbm_noise_analysis/genotype_counts.csv` listing all unique genotype values and their counts.

### 5. Identify High- and Low-Noise Genes
- Calculates coefficient of variation (CV) for each gene across all cells.
- Classifies top 10% as "high-noise" and bottom 10% as "low-noise" genes.

### 6. Identify Responsive Genes
- Finds genes highly expressed in one cell type (identity) but low in others.
- Outputs a list of responsive genes.

### 7. Annotate Genes with gbM and Compare Features
- Annotates genes as gbM/non-gbM, high/low noise, and responsive/non-responsive.
- Compares gbM frequency among these groups.

### 8. Statistical Analysis
- Performs Fisher's exact test to compare gbM frequency in high/low noise and responsive/non-responsive genes.

## How to Run All Steps
1. Export expression matrix:
   ```sh
   sbatch batch_scripts/export_expression_matrix.sh
   ```
2. Check genotype diversity:
   ```sh
   sbatch batch_scripts/check_genotypes.sh
   ```
3. Run noise and responsiveness analysis:
   ```sh
   sbatch batch_scripts/gbm_noise_analysis.sh
   ```

## Output Files
- Expression matrix: `gbm_noise_analysis/expression_matrix.csv`
- Genotype counts: `gbm_noise_analysis/genotype_counts.csv`
- Noise analysis results: `gbm_noise_analysis/gbm_noise_comparison.csv`, `gbm_noise_analysis/gbm_noise_boxplot.png`, `gbm_noise_analysis/gbm_noise_stats.txt`
- Additional gene annotation and responsiveness results will be added as scripts are extended.

## Batch Scripts Location
All batch scripts and their output/error logs are in `mnieuwenh/batch_scripts/`.

---

## Notes
- Adjust memory and CPU in the sbatch scripts if needed (default: 120G RAM, 8 CPUs, 24h for large jobs).
- All scripts assume the conda environment at `/group/sms029/conda_environment/R`.
- You may need to edit gene names if your Seurat object uses a different format than the GBM list.
