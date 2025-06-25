#!/bin/bash
#SBATCH --job-name=check_cellname_overlap
#SBATCH --output=/group/sms029/mnieuwenh/test/logs/check_cellname_overlap.out
#SBATCH --error=/group/sms029/mnieuwenh/test/logs/check_cellname_overlap.err
#SBATCH --partition=work
#SBATCH --mem=16G
#SBATCH --cpus-per-task=1
#SBATCH --time=01:00:00
#SBATCH --mail-user=22938292@student.uwa.edu.au
#SBATCH --mail-type=END,FAIL

module load Anaconda3/2024.06
source ~/.bashrc
conda activate /group/sms029/conda_environment/R

Rscript /group/sms029/mnieuwenh/test/check_cellname_overlap.R
