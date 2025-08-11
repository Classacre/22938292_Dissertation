#!/usr/bin/env bash
#SBATCH --job-name=smk_driver
#SBATCH --partition=work
#SBATCH --cpus-per-task=2
#SBATCH --mem=32G
#SBATCH --time=72:00:00
#SBATCH --mail-user=22938292@student.uwa.edu.au
#SBATCH --mail-type=END,FAIL
#SBATCH --output=/group/sms029/mnieuwenh/workflow2/logs/run_snakemake.out
#SBATCH --error=/group/sms029/mnieuwenh/workflow2/logs/run_snakemake.err

set -euo pipefail

cd /group/sms029/mnieuwenh/workflow2
mkdir -p logs .conda_pkgs

export CONDARC="/group/sms029/mnieuwenh/workflow2/.condarc"
export CONDA_CHANNEL_PRIORITY=strict
export CONDA_PKGS_DIRS="/group/sms029/mnieuwenh/workflow2/.conda_pkgs"

# Load base conda
source /etc/profile.d/modules.sh
module load Anaconda3/2024.06 || true
source "$(conda info --base)/etc/profile.d/conda.sh"

# Ensure a dedicated snakemake env exists
if ! conda env list | awk '{print $1}' | grep -qx "snakemake"; then
  echo "Creating main 'snakemake' environment..."
  conda create -y -n snakemake -c conda-forge -c bioconda python=3.11 "snakemake-minimal>=8" pip
fi
conda activate snakemake

# Ensure conda >= 24.7.1 and libmamba solver inside this env (does not touch base)
echo "Checking conda version inside 'snakemake' environment..."
if ! python - <<'PY'
import sys
try:
    import conda
    from packaging.version import Version
    sys.exit(0 if Version(conda.__version__) >= Version("24.7.1") else 1)
except Exception:
    sys.exit(1)
PY
then
  echo "Upgrading conda and installing libmamba solver..."
  conda install -y -c conda-forge "conda>=24.7.1" conda-libmamba-solver
  conda config --env --set solver libmamba
fi

# Install the Snakemake cluster executor plugin (cluster-generic) if missing
if ! python - <<'PY'
import importlib.util, sys
sys.exit(0 if importlib.util.find_spec("snakemake_executor_plugin_cluster_generic") else 1)
PY
then
  echo "Installing snakemake-executor-plugin-cluster-generic..."
  conda install -y -c conda-forge snakemake-executor-plugin-cluster-generic
fi

# Where Snakemake stores rule envs
export SNAKEMAKE_CONDA_PREFIX="/group/sms029/mnieuwenh/workflow2/.snakemake/conda"

# Clean up any leftover, incomplete env directories
if [ -d "$SNAKEMAKE_CONDA_PREFIX" ]; then
  echo "Cleaning up potentially broken conda environments..."
  while IFS= read -r -d '' d; do
    [ -d "$d/conda-meta" ] || rm -rf "$d"
  done < <(find "$SNAKEMAKE_CONDA_PREFIX" -mindepth 1 -maxdepth 1 -type d -print0)
fi

# Step 1: Pre-create all rule environments locally (no SLURM submission)
echo "Pre-creating rule conda environments locally..."
snakemake \
  --executor local \
  --cores 1 \
  --use-conda \
  --conda-prefix "$SNAKEMAKE_CONDA_PREFIX" \
  --conda-create-envs-only

# Step 2: Run the workflow on SLURM using the profile
SMK_J=100
echo "Starting Snakemake on SLURM..."
snakemake \
  --profile config/slurm_profile \
  --rerun-incomplete \
  -j "$SMK_J" \
  "$@"

echo "Snakemake execution finished."