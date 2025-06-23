#!/bin/bash
#SBATCH --job-name=export_expr
#SBATCH --output=/group/sms029/mnieuwenh/gbm_noise_analysis/logs/export_expr.out
#SBATCH --error=/group/sms029/mnieuwenh/gbm_noise_analysis/logs/export_expr.err
#SBATCH --partition=work
#SBATCH --mem=120G
#SBATCH --cpus-per-task=8
#SBATCH --time=24:00:00
#SBATCH --mail-user=22938292@uwa.edu.au
#SBATCH --mail-type=END,FAIL

module load Anaconda3/2024.06
source ~/.bashrc
conda activate /group/sms029/conda_environment/R

Rscript /group/sms029/mnieuwenh/gbm_noise_analysis/export_expression_matrix.R
