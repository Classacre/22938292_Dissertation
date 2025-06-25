#!/bin/bash
#SBATCH --job-name=export_full_metadata
#SBATCH --output=/group/sms029/mnieuwenh/seurat_metadata/logs/export_full_metadata.out
#SBATCH --error=/group/sms029/mnieuwenh/seurat_metadata/logs/export_full_metadata.err
#SBATCH --partition=work
#SBATCH --mem=16G
#SBATCH --cpus-per-task=1
#SBATCH --time=01:00:00
#SBATCH --mail-user=22938292@student.uwa.edu.au
#SBATCH --mail-type=END,FAIL

module load Anaconda3/2024.06
source ~/.bashrc
conda activate /group/sms029/conda_environment/R

Rscript /group/sms029/mnieuwenh/seurat_metadata/export_full_metadata.R
