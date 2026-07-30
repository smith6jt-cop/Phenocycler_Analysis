#!/bin/bash
# Steps 5-8 — the norm floor, broad lineage, QuPath export and identity figures.
#
# Driven through the orchestrator rather than the module CLIs so the hormone_floor stage gets
# its pre-floor backup (restore_gated_redsea.pre_hormonefloor). That backup is what makes a
# later threshold change reproducible from clean RESTORE gates instead of re-flooring an
# already-floored table.
#
# The floor is NOT optional: lineage.py's endocrine and immune _pos calls assume it has run,
# and an earlier version of this chain went straight from restore to lineage, which silently
# reinstated the false-endocrine and false-immune over-calling the floor exists to remove.
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
# `conda activate` needs conda's shell hook, which a non-interactive SLURM shell does
# not source. Without this the job dies immediately under `set -e` with "shell has not
# been properly configured" — before any Python runs, so the error is opaque.
module load conda 2>/dev/null || true
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate phenocycler_analysis
mkdir -p logs

# Re-running this job is safe: hormone_floor is idempotent (it recomputes _pos from _norm,
# which it never touches), and the other three skip when their outputs exist. After CHANGING
# a threshold, add --force — `lineage` is judged current by column shape, not by the K it was
# computed at, so it would otherwise be skipped and the change would appear to do nothing.
python -m phenocycler.pipeline --only hormone_floor lineage qupath figures \
                               --jobs "${SLURM_CPUS_PER_TASK:-8}"

# The acceptance yardstick for any positivity change. Read-only, no stage depends on it, and
# it is the only thing in the repo that says whether a threshold change helped or just moved
# cells around.
python -m phenocycler.reassess_diag || echo "[warn] reassess_diag failed (non-fatal)"
