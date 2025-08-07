#!/bin/bash

#==============================================================================
# SLURM SCRIPT FOR SUPERVISOR'S REQUESTED PLOTS (SCRIPT 06)
#==============================================================================

# --- SLURM JOB CONFIGURATION ---
# Job name for easy identification in the queue
#SBATCH --job-name=supervisor_plots

# Paths for the output and error logs
#SBATCH --output=/group/sms029/mnieuwenh/gbm_noise_analysis/logs/supervisor_plots.out
#SBATCH --error=/group/sms029/mnieuwenh/gbm_noise_analysis/logs/supervisor_plots.err

# Resource allocation
#SBATCH --partition=work         # The partition (queue) to run on
#SBATCH --mem=128G               # Memory request (generous to handle large data files)
#SBATCH --cpus-per-task=8        # Number of CPU cores
#SBATCH --time=02:00:00          # Time limit (2 hours should be plenty)

# Email notifications for job status
#SBATCH --mail-user=22938292@student.uwa.edu.au
#SBATCH --mail-type=END,FAIL     # Notify on job end or failure


# --- ENVIRONMENT SETUP ---
echo "========================================================"
echo "Job started on $(date)"
echo "Running on node: $(hostname)"
echo "========================================================"

echo "Loading Anaconda module..."
module load Anaconda3/2024.06

echo "Sourcing .bashrc to initialize Conda..."
source ~/.bashrc

echo "Activating Conda environment: /group/sms029/conda_environment/R"
conda activate /group/sms029/conda_environment/R

# Check if activation was successful
if [ $? -eq 0 ]; then
    echo "Environment activated successfully."
    echo "Using R from: $(which R)"
else
    echo "ERROR: Failed to activate Conda environment. Halting job."
    exit 1
fi

# --- IMPORTANT PRE-REQUISITE ---
# This script assumes that the required R packages have already been installed
# via Conda. If not, cancel this job and run this command in your terminal:
#
conda install -c conda-forge r-dplyr r-tidyr r-ggplot2 r-patchwork r-ggpubr r-readr
#

# --- SCRIPT EXECUTION ---
echo ""
echo "--- Starting Step 6: Generate Supervisor's Requested Plots ---"
Rscript /group/sms029/mnieuwenh/gbm_noise_analysis/06_supervisor_requests.R

# Check the exit code of the R script
if [ $? -eq 0 ]; then
    echo "--- Finished Step 6 successfully. ---"
else
    echo "--- ERROR: Step 6 failed. Check the supervisor_plots.err log file. ---"
    exit 1 # Exit with an error code
fi

echo ""
echo "========================================================"
echo "Job completed successfully on $(date)."
echo "========================================================"