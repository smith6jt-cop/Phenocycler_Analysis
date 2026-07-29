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
# `conda activate` needs conda's shell hook, which a non-interactive SLURM shell does
# not source. Without this the job dies immediately under `set -e` with "shell has not
# been properly configured" — before any Python runs, so the error is opaque.
module load conda 2>/dev/null || true
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate phenocycler_analysis
mkdir -p logs

# Reads data/cells_redsea by default; writes data/restore_gated_redsea + thresholds.
python -m phenocycler.restore --jobs "${SLURM_CPUS_PER_TASK:-16}"

# The extra pass (CD99 / B3TUBB / MPO -> restore_gated_redsea_extra/). Separate on purpose so
# the validated 10-marker gates stay byte-identical, but NOT optional: lineage.py raises
# FileNotFoundError("extra-gated markers missing for donor ...") without it, so omitting this
# meant step 4 died for every donor after step 3 had already burned its 24 h.
python -m phenocycler.restore --extra --jobs "${SLURM_CPUS_PER_TASK:-16}"
