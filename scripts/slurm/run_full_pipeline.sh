#!/bin/bash
# Full raw-data -> broad-lineage pipeline as a dependency chain on HiPerGator.
# Submits: cells -> redsea (array) -> restore(+extra) -> floor/lineage/qupath/figures,
# each waiting on the last.
#
# The array is sized in SECTIONS, not donors: the work unit is one image, so a donor with a
# pancreas and a lymph node contributes two (`6539` and `6539pln`). Pass the count, or omit it
# and the count is read from data/cells/donor_id=* once step 1 has run.
#
#   bash scripts/slurm/run_full_pipeline.sh 40      # 40 sections -> array 0-39
#
# SBATCH credentials: either edit the three placeholders in each child script, or export
# SBATCH_ACCOUNT / SBATCH_QOS before running this (sbatch input environment variables take
# precedence over in-script #SBATCH directives). There is no environment equivalent for
# --mail-user, so that line still needs editing or deleting.
set -euo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")/../.."   # repo root
mkdir -p logs

# Sections, not donors. If not given, ask the config — but only step 1 can answer, so an
# explicit count is required on a cold start.
NSECTIONS="${1:-}"
if [[ -z "$NSECTIONS" ]]; then
  NSECTIONS=$(python -c "from phenocycler import load_config; print(len(load_config().discover_sections()))" 2>/dev/null || echo 0)
  if [[ "$NSECTIONS" -eq 0 ]]; then
    echo "No sections found under data/cells/donor_id=* — pass the section count explicitly:" >&2
    echo "  bash scripts/slurm/run_full_pipeline.sh 40" >&2
    exit 1
  fi
  echo "Discovered ${NSECTIONS} sections under data/cells/"
fi
LAST=$((NSECTIONS - 1))

cells=$(sbatch --parsable scripts/slurm/01_cells_parquet.sh)
echo "cells   : $cells"

redsea=$(sbatch --parsable --dependency=afterok:${cells} \
                --array=0-${LAST} scripts/slurm/02_redsea_array.sh)
echo "redsea  : $redsea (array 0-${LAST})"

restore=$(sbatch --parsable --dependency=afterok:${redsea} scripts/slurm/03_restore.sh)
echo "restore : $restore"

lineage=$(sbatch --parsable --dependency=afterok:${restore} scripts/slurm/04_lineage.sh)
echo "lineage : $lineage  (floor -> lineage -> qupath -> figures -> reassess_diag)"

echo "Submitted chain. Monitor with:  squeue -u \$USER"
echo "Status after completion:  python -m phenocycler.pipeline --status"
echo ""
echo "NOTE: afterok only checks the job's exit code. Per-section failures inside REDSEA and"
echo "RESTORE are LOGGED, not raised (parallel.map_donors defaults to on_error=log), and"
echo "--status only checks that a stage wrote at least one partition. Count partitions"
echo "between stages before trusting the chain:"
echo "  for d in cells cells_redsea restore_gated_redsea restore_gated_redsea_extra; do"
echo "    echo \"\$d \$(ls -d data/\$d/donor_id=* 2>/dev/null | wc -l)\"; done"
