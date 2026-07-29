#!/bin/bash
# Step 2 — REDSEA spillover correction, ONE SECTION PER ARRAY TASK.
#
# A SECTION, not a donor. The work unit is one image: one qptiff, one GeoJSON, one
# coordinate frame. A donor with a pancreas and a lymph node contributes TWO — `6539`
# and `6539pln` — so indexing donors would leave every pLN section uncorrected while
# the job still exited 0. Size the array off `discover_sections()`:
#
#   N=$(python -c "from phenocycler import load_config; print(len(load_config().discover_sections()))")
#   sbatch --array=0-$((N-1)) scripts/slurm/02_redsea_array.sh
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
# `conda activate` needs conda's shell hook, which a non-interactive SLURM shell does
# not source. Without this the job dies immediately under `set -e` with "shell has not
# been properly configured" — before any Python runs, so the error is opaque.
module load conda 2>/dev/null || true
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate phenocycler_analysis
mkdir -p logs

# Map this array index to a SECTION key (sorted, matches cfg.discover_sections(), which is
# also what `redsea --donor` resolves against).
SECTION=$(python - <<'PY'
import os
from phenocycler import load_config
sections = load_config().discover_sections()
i = int(os.environ["SLURM_ARRAY_TASK_ID"])
if i >= len(sections):
    raise SystemExit(f"array index {i} >= {len(sections)} sections; shrink --array range")
print(sections[i])
PY
)
echo "[array ${SLURM_ARRAY_TASK_ID}] section=${SECTION}"

# One section per task -> serial within the task (the array IS the parallelism).
# --keep-mask/--save-intermediates cost disk but make `--from-intermediates` possible, which
# turns an alpha / comp-mode sweep from a full qptiff re-read into seconds. Drop them if you
# are not sweeping REDSEA parameters.
python -m phenocycler.redsea --donor "${SECTION}" --keep-mask --save-intermediates
