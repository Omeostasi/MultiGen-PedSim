#!/bin/bash
#SBATCH --job-name=founder_benchmark
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --mem=100G
#SBATCH --time=12:00:00
#SBATCH --output=logs/founder_benchmark_%j.log
#SBATCH --error=logs/founder_benchmark_%j.err

# ── Setup ─────────────────────────────────────────────────────────────────────
set -euo pipefail

mkdir -p logs results

echo "=========================================="
echo "Job ID      : ${SLURM_JOB_ID}"
echo "Node        : ${SLURMD_NODENAME}"
echo "CPUs        : ${SLURM_CPUS_PER_TASK}"
echo "Start time  : $(date '+%Y-%m-%dT%H:%M:%S')"
echo "=========================================="

# ── Conda environment ─────────────────────────────────────────────────────────
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate alphasim

# ── Run ───────────────────────────────────────────────────────────────────────
Rscript time_to_generate_founders.R

EXIT_CODE=$?

echo "=========================================="
echo "End time    : $(date '+%Y-%m-%dT%H:%M:%S')"
echo "Exit code   : ${EXIT_CODE}"
echo "=========================================="

exit ${EXIT_CODE}
