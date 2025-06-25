#!/bin/bash
#SBATCH --job-name=check_seurat_structure
#SBATCH --output=/group/sms029/mnieuwenh/test/logs/seurat_summary.out
#SBATCH --error=/group/sms029/mnieuwenh/test/logs/seurat_summary.err
#SBATCH --partition=work
#SBATCH --mem=120G
#SBATCH --cpus-per-task=8
#SBATCH --time=24:00:00
#SBATCH --mail-user=22938292@student.uwa.edu.au
#SBATCH --mail-type=END,FAIL

module load Anaconda3/2024.06
source ~/.bashrc
conda activate /group/sms029/conda_environment/R

Rscript /group/sms029/mnieuwenh/test/check_seurat_structure.R
