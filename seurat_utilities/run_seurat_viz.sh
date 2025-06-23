#!/bin/bash
#SBATCH --job-name=seurat_viz
#SBATCH --output=/group/sms029/mnieuwenh/batch_scripts/logs/seurat_viz.out
#SBATCH --error=/group/sms029/mnieuwenh/batch_scripts/logs/seurat_viz.err
#SBATCH --partition=work
#SBATCH --mem=120G
#SBATCH --cpus-per-task=8
#SBATCH --time=24:00:00

module load Anaconda3/2024.06
source ~/.bashrc
conda activate /group/sms029/conda_environment/R

Rscript /group/sms029/mnieuwenh/seurat_utilities/Seurat_Visualization.R