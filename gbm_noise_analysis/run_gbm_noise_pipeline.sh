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

module load Anaconda3/2024.06
source ~/.bashrc
conda activate /group/sms029/conda_environment/R

echo "Environment activated successfully."

# --- PIPELINE EXECUTION ---
# The 'echo' statements will print progress to your .out file.
# If a script fails, the job will now stop immediately.

echo "--- Starting Step 1: Calculate Noise ---"
Rscript /group/sms029/mnieuwenh/gbm_noise_analysis/01_calculate_noise.R
echo "--- Finished Step 1 ---"

#echo "--- Starting Step 2: Analyze Responsive Genes ---"
#Rscript /group/sms029/mnieuwenh/gbm_noise_analysis/02_analyze_responsive_genes.R
#echo "--- Finished Step 2 ---"

#echo "--- Starting Step 3: Create Master Summary ---"
#Rscript /group/sms029/mnieuwenh/gbm_noise_analysis/03_create_master_gene_summary.R
#echo "--- Finished Step 3 ---"

#echo "--- Starting Step 4: Subsampling Analysis ---"
#Rscript /group/sms029/mnieuwenh/gbm_noise_analysis/04_subsampling_validation.R
#echo "--- Finished Step 4 ---"

#echo "--- Starting Step 5: Visualize Filtering ---"
#Rscript /group/sms029/mnieuwenh/gbm_noise_analysis/05_visualize_filtering.R
#echo "--- Finished Step 5 ---"

echo "--- Starting Step 6: Visualize Filtering ---"
Rscript /group/sms029/mnieuwenh/gbm_noise_analysis/06_supervisor_requests.R
echo "--- Finished Step 6 ---"

echo "Pipeline completed successfully."