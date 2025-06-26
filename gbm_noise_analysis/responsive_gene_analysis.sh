#!/bin/bash
#SBATCH --job-name=responsive_gene_analysis
#SBATCH --output=/group/sms029/mnieuwenh/gbm_noise_analysis/logs/responsive_gene_analysis.out
#SBATCH --error=/group/sms029/mnieuwenh/gbm_noise_analysis/logs/responsive_gene_analysis.err
#SBATCH --partition=work
#SBATCH --mem=32G
#SBATCH --cpus-per-task=2
#SBATCH --time=04:00:00
#SBATCH --mail-user=22938292@student.uwa.edu.au
#SBATCH --mail-type=END,FAIL

module load Anaconda3/2024.06
source ~/.bashrc
conda activate /group/sms029/conda_environment/R

# Install required R packages via conda if not present
conda install -y -c conda-forge r-ggridges r-complexupset r-pheatmap

Rscript /group/sms029/mnieuwenh/gbm_noise_analysis/responsive_gene_analysis.R
