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

# --- SCRIPT SETUP ---
set -e

echo "================================================================"
echo "Job starting at: $(date)"
echo "================================================================"

# --- ENVIRONMENT SETUP ---
# This is now extremely fast as the environment is 100% pre-built.
echo "--- Loading modules and activating final Conda environment ---"
module load Anaconda3/2024.06
source ~/.bashrc
# Activate the final, correct environment
conda activate R_final_env
echo "Conda environment activated: $(which R)"
echo "----------------------------------------------------------------"

# --- NO INSTALLATION STEPS ARE NEEDED ---

# # --- PIPELINE EXECUTION ---
echo "--- Starting main analysis pipeline ---"
Rscript /group/sms029/mnieuwenh/gbm_noise_analysis/00_data_overview_and_setup.R
echo "--- Finished Step 0 ---"

Rscript /group/sms029/mnieuwenh/gbm_noise_analysis/01_prepare_epigenetic_data.R
echo "--- Finished Step 1 ---"

Rscript /group/sms029/mnieuwenh/gbm_noise_analysis/02_prepare_expression_data.R
echo "--- Finished Step 2 ---"

Rscript /group/sms029/mnieuwenh/gbm_noise_analysis/03_calculate_expression_noise.R
echo "--- Finished Step 3 ---"

Rscript /group/sms029/mnieuwenh/gbm_noise_analysis/04_validate_gene_filtering.R
echo "--- Finished Step 4 ---"

Rscript /group/sms029/mnieuwenh/gbm_noise_analysis/05_analyze_responsive_genes.R
echo "--- Finished Step 5 ---"

Rscript /group/sms029/mnieuwenh/gbm_noise_analysis/06_generate_summary_figures.R
echo "--- Finished Step 6 ---"

Rscript /group/sms029/mnieuwenh/gbm_noise_analysis/07_GO_enrichment_analysis.R
echo "--- Finished Step 7 ---"

Rscript /group/sms029/mnieuwenh/gbm_noise_analysis/08_bootstrap_analysis.R
echo "--- Finished Step 8 ---"

echo "================================================================"
echo "Pipeline completed successfully at: $(date)"
echo "================================================================"