#!/bin/bash
# Script to move and rename analysis scripts for project consistency
# Also updates references in batch scripts and README files

set -e

# Move and rename R script.R to seurat_metadata_summary.R
if [ -f "misc_files/R script.R" ]; then
  mv "misc_files/R script.R" "seurat_metadata/seurat_metadata_summary.R"
  echo "Moved: misc_files/R script.R -> seurat_metadata/seurat_metadata_summary.R"
fi

# Move Seurat_Visualization.R to seurat_utilities
if [ -f "misc_files/Seurat_Visualization.R" ]; then
  mv "misc_files/Seurat_Visualization.R" "seurat_utilities/Seurat_Visualization.R"
  echo "Moved: misc_files/Seurat_Visualization.R -> seurat_utilities/Seurat_Visualization.R"
fi

# Update references in run_seurat_summary.sh
sed -i 's|R script.R|seurat_metadata/seurat_metadata_summary.R|g' seurat_utilities/run_seurat_summary.sh

# Update references in README.md and README_CLEANUP.txt
sed -i 's|R script.R|seurat_metadata/seurat_metadata_summary.R|g' README.md
if [ -f "misc_files/README_CLEANUP.txt" ]; then
  sed -i 's|R script.R|seurat_metadata/seurat_metadata_summary.R|g' misc_files/README_CLEANUP.txt
fi

# Update references in seurat_metadata_sample.csv output path in seurat_metadata_summary.R
sed -i 's|/group/sms029/mnieuwenh/seurat_metadata_sample.csv|/group/sms029/mnieuwenh/seurat_metadata/seurat_metadata_sample.csv|g' seurat_metadata/seurat_metadata_summary.R

# Update references in seurat_summary.txt output path in seurat_metadata_summary.R
sed -i 's|/group/sms029/mnieuwenh/seurat_summary.txt|/group/sms029/mnieuwenh/seurat_metadata/seurat_summary.txt|g' seurat_metadata/seurat_metadata_summary.R

# Update references in run_seurat_viz.sh
sed -i 's|Seurat_Visualization.R|seurat_utilities/Seurat_Visualization.R|g' seurat_utilities/run_seurat_viz.sh

# Done
echo "All moves and reference updates complete."
