#!/bin/bash
# One manifest-driven, content-addressed workflow job.
#
# Resource/account policy belongs on the submission command, for example:
#   sbatch --account=<account> --qos=<qos> --cpus-per-task=8 --mem=120G \
#     --time=2-00:00:00 scripts/slurm/run_full_pipeline.sh config.ini
#
# The cohort manifest referenced by config.ini supplies the exact donor set and
# immutable QuPath inputs. Stage manifests provide restart/current-state logic;
# no hand-maintained donor count, array range, or dependent legacy scripts are
# involved.
set -euo pipefail

REPOSITORY_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
cd "${REPOSITORY_DIR}"

CONFIG_PATH="${1:-config.ini}"
JOBS="${PHENOCYCLER_JOBS:-${SLURM_CPUS_PER_TASK:-1}}"
if [[ ! -f "${CONFIG_PATH}" ]]; then
  echo "ERROR: config file not found: ${CONFIG_PATH}" >&2
  exit 2
fi

exec python -m phenocycler.pipeline run \
  --config "${CONFIG_PATH}" \
  --jobs "${JOBS}"
