#!/bin/bash
# Step 1 — build the donor-partitioned cells parquet from the QuPath CSV (DuckDB).
# One job; DuckDB streams the ~138 GB CSV out-of-core with N threads.
#SBATCH --job-name=pheno_cells
#SBATCH --output=logs/cells_%j.out
#SBATCH --error=logs/cells_%j.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64gb
#SBATCH --time=12:00:00
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

# --threads matches --cpus-per-task; --memory-limit keeps DuckDB inside the alloc.
python -m phenocycler.cells_parquet --threads "${SLURM_CPUS_PER_TASK:-8}" --memory-limit 56GB
