#!/bin/bash
#SBATCH --job-name=gbm_noise_pipeline
#SBATCH --output=/group/sms029/mnieuwenh/gbm_noise_analysis/logs/gbm_noise_pipeline.out
#SBATCH --error=/group/sms029/mnieuwenh/gbm_noise_analysis/logs/gbm_noise_pipeline.err
#SBATCH --partition=work
#SBATCH --mem=128G
#SBATCH --cpus-per-task=8
#SBATCH --time=48:00:00
#SBATCH --mail-user=22938292@student.uwa.edu.au
#SBATCH --mail-type=END,FAIL

# --- ENVIRONMENT SETUP ---
# Exit immediately if a command exits with a non-zero status.
# This is crucial for preventing the pipeline from continuing if a step fails.
set -e

# Load necessary modules and activate the Conda environment
module load Anaconda3/2024.06
source ~/.bashrc
conda activate /group/sms029/conda_environment/R

echo "Environment activated successfully. Starting pipeline execution..."
echo "----------------------------------------------------------------"

# --- PIPELINE EXECUTION ---
# # The 'echo' statements will print progress to your .out file.
# # The 'set -e' command ensures that the job will stop immediately if any script fails.

# echo "--- Starting Step 0: Data Overview and Setup ---"
# Rscript /group/sms029/mnieuwenh/gbm_noise_analysis/00_data_overview_and_setup.R
# echo "--- Finished Step 0 ---"

# echo "--- Starting Step 1: Prepare Epigenetic Data ---"
# Rscript /group/sms029/mnieuwenh/gbm_noise_analysis/01_prepare_epigenetic_data.R
# echo "--- Finished Step 1 ---"

# echo "--- Starting Step 2: Prepare Expression Data ---"
# Rscript /group/sms029/mnieuwenh/gbm_noise_analysis/02_prepare_expression_data.R
# echo "--- Finished Step 2 ---"

# echo "--- Starting Step 3: Calculate Expression Noise ---"
# Rscript /group/sms029/mnieuwenh/gbm_noise_analysis/03_calculate_expression_noise.R
# echo "--- Finished Step 3 ---"

# echo "--- Starting Step 4: Validate Gene Filtering ---"
# Rscript /group/sms029/mnieuwenh/gbm_noise_analysis/04_validate_gene_filtering.R
# echo "--- Finished Step 4 ---"

# echo "--- Starting Step 5: Analyze Responsive Genes ---"
# Rscript /group/sms029/mnieuwenh/gbm_noise_analysis/05_analyze_responsive_genes.R
# echo "--- Finished Step 5 ---"

# echo "--- Starting Step 6: Generate Summary Figures ---"
# Rscript /group/sms029/mnieuwenh/gbm_noise_analysis/06_generate_summary_figures.R
# echo "--- Finished Step 6 ---"

echo "--- Starting Step 7: Perform GO Enrichment Analysis ---"
Rscript /group/sms029/mnieuwenh/gbm_noise_analysis/07_perform_go_enrichment.R
echo "--- Finished Step 7 ---"

echo "----------------------------------------------------------------"
echo "Pipeline completed successfully."