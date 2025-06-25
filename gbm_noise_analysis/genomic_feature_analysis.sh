#!/bin/bash
#SBATCH --job-name=genomic_feature_analysis
#SBATCH --output=/group/sms029/mnieuwenh/gbm_noise_analysis/results/genomic_features/genomic_feature_analysis.out
#SBATCH --error=/group/sms029/mnieuwenh/gbm_noise_analysis/results/genomic_features/genomic_feature_analysis.err
#SBATCH --partition=work
#SBATCH --mem=16G
#SBATCH --cpus-per-task=1
#SBATCH --time=02:00:00
#SBATCH --mail-user=22938292@student.uwa.edu.au
#SBATCH --mail-type=END,FAIL

module load Anaconda3/2024.06
source ~/.bashrc
conda activate /group/sms029/conda_environment/R

Rscript /group/sms029/mnieuwenh/gbm_noise_analysis/genomic_feature_analysis.R
