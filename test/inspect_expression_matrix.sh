#!/bin/bash
#SBATCH --job-name=inspect_expr_matrix
#SBATCH --output=/group/sms029/mnieuwenh/test/logs/inspect_expression_matrix.out
#SBATCH --error=/group/sms029/mnieuwenh/test/logs/inspect_expression_matrix.err
#SBATCH --partition=work
#SBATCH --mem=120G
#SBATCH --cpus-per-task=2
#SBATCH --time=00:30:00
#SBATCH --mail-user=22938292@student.uwa.edu.au
#SBATCH --mail-type=END,FAIL

module load Anaconda3/2024.06
source ~/.bashrc
conda activate /group/sms029/conda_environment/R

Rscript /group/sms029/mnieuwenh/test/inspect_expression_matrix.R
