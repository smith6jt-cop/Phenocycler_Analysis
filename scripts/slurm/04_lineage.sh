#!/bin/bash
# Step 4 (+5) — broad lineage assignment (parallel over donors) + QuPath export.
#SBATCH --job-name=pheno_lineage
#SBATCH --output=logs/lineage_%j.out
#SBATCH --error=logs/lineage_%j.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64gb
#SBATCH --time=04:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=your_email@ufl.edu
#SBATCH --qos=your_qos
#SBATCH --account=your_account

set -euo pipefail
module load conda 2>/dev/null || true
conda activate phenocycler_analysis
mkdir -p logs

python -m phenocycler.lineage --jobs "${SLURM_CPUS_PER_TASK:-8}"
python -m phenocycler.qupath_export
