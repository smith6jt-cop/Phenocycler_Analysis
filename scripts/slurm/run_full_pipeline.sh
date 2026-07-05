#!/bin/bash
# Full raw-data -> broad-lineage pipeline as a dependency chain on HiPerGator.
# Submits: cells -> redsea (array) -> restore -> lineage, each waiting on the last.
# Edit --array range to match your donor count, and set the SBATCH creds in each
# child script (or export SBATCH_ACCOUNT / SBATCH_QOS before running this).
#
#   bash scripts/slurm/run_full_pipeline.sh 15      # 15 donors -> array 0-14
#
set -euo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")/../.."   # repo root
mkdir -p logs

NDONORS="${1:-15}"
LAST=$((NDONORS - 1))

cells=$(sbatch --parsable scripts/slurm/01_cells_parquet.sh)
echo "cells   : $cells"

redsea=$(sbatch --parsable --dependency=afterok:${cells} \
                --array=0-${LAST} scripts/slurm/02_redsea_array.sh)
echo "redsea  : $redsea (array 0-${LAST})"

restore=$(sbatch --parsable --dependency=afterok:${redsea} scripts/slurm/03_restore.sh)
echo "restore : $restore"

lineage=$(sbatch --parsable --dependency=afterok:${restore} scripts/slurm/04_lineage.sh)
echo "lineage : $lineage"

echo "Submitted chain. Monitor with:  squeue -u \$USER"
echo "Status after completion:  python -m phenocycler.pipeline --status"
