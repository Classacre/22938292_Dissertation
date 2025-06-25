#!/bin/bash
#SBATCH --job-name=replace_seurat_path
#SBATCH --output=/group/sms029/mnieuwenh/test/replace_seurat_path.out
#SBATCH --error=/group/sms029/mnieuwenh/test/replace_seurat_path.err
#SBATCH --partition=work
#SBATCH --mem=4G
#SBATCH --cpus-per-task=1
#SBATCH --time=00:10:00
#SBATCH --mail-user=22938292@student.uwa.edu.au
#SBATCH --mail-type=END,FAIL

bash /group/sms029/mnieuwenh/test/replace_seurat_path.sh
