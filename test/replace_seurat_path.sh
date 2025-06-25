#!/bin/bash
# Script to replace all references to the old Seurat object with the new col-only version
# Searches all files in /group/sms029/mnieuwenh recursively

find /group/sms029/mnieuwenh -type f -exec sed -i 's|/group/sms029/Oliva_dataset/integrated_col_trajectories_colonly.rds|/group/sms029/Oliva_dataset/integrated_col_trajectories_colonly.rds|g' {} +

echo "All references updated to use integrated_col_trajectories_colonly.rds."
