#!/bin/bash
#SBATCH --job-name=check_genotypes
#SBATCH --output=/group/sms029/mnieuwenh/gbm_noise_analysis/logs/check_genotypes.out
#SBATCH --error=/group/sms029/mnieuwenh/gbm_noise_analysis/logs/check_genotypes.err
#SBATCH --partition=work
#SBATCH --mem=16G
#SBATCH --cpus-per-task=2
#SBATCH --time=01:00:00
#SBATCH --mail-user=22938292@student.uwa.edu.au
#SBATCH --mail-type=END,FAIL

module load Anaconda3/2024.06
source ~/.bashrc
conda activate /group/sms029/conda_environment/R

Rscript /group/sms029/mnieuwenh/gbm_noise_analysis/check_genotypes.R
