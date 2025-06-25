#!/bin/bash
#SBATCH --job-name=remove_na_genotype
#SBATCH --output=/group/sms029/mnieuwenh/seurat_utilities/remove_na_genotype.out
#SBATCH --error=/group/sms029/mnieuwenh/seurat_utilities/remove_na_genotype.err
#SBATCH --partition=work
#SBATCH --mem=120G
#SBATCH --cpus-per-task=8
#SBATCH --time=08:00:00
#SBATCH --mail-user=22938292@student.uwa.edu.au
#SBATCH --mail-type=END,FAIL

module load Anaconda3/2024.06
source ~/.bashrc
conda activate /group/sms029/conda_environment/R

Rscript /group/sms029/mnieuwenh/seurat_utilities/remove_na_genotype.R
