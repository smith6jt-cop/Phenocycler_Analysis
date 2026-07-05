#!/bin/bash
# Step 3 — RESTORE normalization on the REDSEA-corrected cohort.
# The RESTORE fit pools all images (one job); the per-cell gate/normalize APPLY is
# per-donor independent and runs across --jobs processes inside this job.
#SBATCH --job-name=pheno_restore
#SBATCH --output=logs/restore_%j.out
#SBATCH --error=logs/restore_%j.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=128gb
#SBATCH --time=24:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=your_email@ufl.edu
#SBATCH --qos=your_qos
#SBATCH --account=your_account

set -euo pipefail
module load conda 2>/dev/null || true
conda activate phenocycler_analysis
mkdir -p logs

# Reads data/cells_redsea by default; writes data/restore_gated_redsea + thresholds.
python -m phenocycler.restore --jobs "${SLURM_CPUS_PER_TASK:-16}"
