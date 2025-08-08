#!/bin/bash
#SBATCH --job-name=seurat_diagnostic
#SBATCH --output=/group/sms029/mnieuwenh/gbm_noise_analysis/logs/seurat_diagnostic.out
#SBATCH --error=/group/sms029/mnieuwenh/gbm_noise_analysis/logs/seurat_diagnostic.err
#SBATCH --partition=work
#SBATCH --mem=64G
#SBATCH --cpus-per-task=4
#SBATCH --time=01:00:00
#SBATCH --mail-user=22938292@student.uwa.edu.au
#SBATCH --mail-type=END,FAIL

# --- ENVIRONMENT SETUP ---
set -e
module load Anaconda3/2024.06
source ~/.bashrc
conda activate /group/sms029/conda_environment/R

echo "Environment activated. Running Seurat diagnostic script..."
echo "---------------------------------------------------------"

# --- EXECUTE THE DIAGNOSTIC SCRIPT ---
Rscript /group/sms029/mnieuwenh/test/diagnose_seurat_object.R

echo "---------------------------------------------------------"
echo "Diagnostic script finished."