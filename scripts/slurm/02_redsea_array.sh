#!/bin/bash
# Step 2 — REDSEA spillover correction, ONE DONOR PER ARRAY TASK.
# Every donor is fully independent (own GeoJSON, qptiff, output parquet), so a
# SLURM array is the cleanest cohort-scale parallelism. Set --array=0-(N-1) where
# N = number of donors under data/cells/donor_id=* (e.g. --array=0-14 for 15).
#
#   sbatch --array=0-14 scripts/slurm/02_redsea_array.sh
#
#SBATCH --job-name=pheno_redsea
#SBATCH --output=logs/redsea_%A_%a.out
#SBATCH --error=logs/redsea_%A_%a.err
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
module load conda 2>/dev/null || true
conda activate phenocycler_analysis
mkdir -p logs

# Map this array index to a donor id (sorted, matches cfg.discover_donors()).
DONOR=$(python - <<'PY'
import os
from phenocycler import load_config
donors = load_config().discover_donors()
i = int(os.environ["SLURM_ARRAY_TASK_ID"])
if i >= len(donors):
    raise SystemExit(f"array index {i} >= {len(donors)} donors; shrink --array range")
print(donors[i])
PY
)
echo "[array ${SLURM_ARRAY_TASK_ID}] donor=${DONOR}"

# One donor per task -> serial within the task (the array IS the parallelism).
python -m phenocycler.redsea --donor "${DONOR}" --keep-mask --save-intermediates
