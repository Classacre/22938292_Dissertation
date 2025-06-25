#!/bin/bash
#SBATCH --job-name=bfg_clean
#SBATCH --output=/group/sms029/mnieuwenh/test/bfg_clean.out
#SBATCH --error=/group/sms029/mnieuwenh/test/bfg_clean.err
#SBATCH --partition=work
#SBATCH --mem=8G
#SBATCH --cpus-per-task=1
#SBATCH --time=01:00:00
#SBATCH --mail-user=22938292@student.uwa.edu.au
#SBATCH --mail-type=END,FAIL

# Set variables
REPO_DIR="/group/sms029/mnieuwenh"
BFG_JAR="bfg-1.14.0.jar"
BFG_URL="https://repo1.maven.org/maven2/com/madgag/bfg/1.14.0/bfg-1.14.0.jar"

cd "$REPO_DIR"

# Download BFG Repo-Cleaner if not present
if [ ! -f "$BFG_JAR" ]; then
  wget "$BFG_URL" -O "$BFG_JAR"
fi

# Run BFG to remove all files over 2GB from git history
java -jar "$BFG_JAR" --strip-blobs-bigger-than 2G

# Clean and update the repo
cd "$REPO_DIR"
git reflog expire --expire=now --all
git gc --prune=now --aggressive

git push --force

# Remove BFG jar and this script after completion
rm -f "$BFG_JAR"
rm -f "$0"
