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

# Run scripts in order
Rscript /group/sms029/mnieuwenh/gbm_noise_analysis/01_calculate_noise.R
Rscript /group/sms029/mnieuwenh/gbm_noise_analysis/02_analyze_responsive_genes.R
Rscript /group/sms029/mnieuwenh/gbm_noise_analysis/03_create_master_gene_summary.R
Rscript /group/sms029/mnieuwenh/gbm_noise_analysis/04_subsampling_analysis.R
Rscript /group/sms029/mnieuwenh/gbm_noise_analysis/05_visualize_filtering.R