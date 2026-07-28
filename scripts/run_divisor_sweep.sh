#!/usr/bin/env bash
# Divisor-policy sweep driver — one `restore_validation evaluate-locked` run per
# (control_definition, divisor_statistic) arm, into an isolated output directory.
#
# WHY THIS EXISTS: the sweep arms under data/restore_pair_validation/sweep_qc/ and sweep_v10/ were
# produced by ad-hoc CLI invocations that were never recorded anywhere — `grep -rIn sweep_qc` over the
# repo returns nothing, and ~/.bash_history predates the runs. The arms are the evidence a Gate-4
# divisor decision rests on, so they have to be reproducible.
#
# The screen writes NO normalized values, NO per-cell calls and NO phenotype output. It only emits the
# per-(donor, pair) evaluation tables. Canonical data/ is never written.
#
#   ./scripts/run_divisor_sweep.sh <output-root> [--pair-set SET] [--donor ID]... [--arm ARM]...
#
# Arms are "<control_definition>/<divisor_statistic>", e.g. reference_only/p99. With no --arm the
# default matrix is the four that define the policy question plus the two robust quantiles.
#
# Example — reproduce data/restore_pair_validation/sweep_qc:
#   ./scripts/run_divisor_sweep.sh data/restore_pair_validation/sweep_qc_repro \
#       --donor 115 --donor 6380 --donor 6436 --donor 6450 --donor 6476 --donor 6539 --donor 6591
set -euo pipefail

# print the header comment block (line 2 up to the first non-comment line)
usage() { sed -n '2,${/^[^#]/q;s/^# \{0,1\}//;p}' "$0"; exit "${1:-0}"; }
[[ $# -ge 1 ]] || usage 1
case "$1" in -h|--help) usage 0 ;; esac

OUT_ROOT="$1"; shift

# Resolve the repo root from this script's location, so the relative data paths below work from
# anywhere. scripts/ lives inside the Phenocycler_Analysis submodule; data/ lives in its parent.
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/../.." && pwd)"

PAIR_SET="accepted"
DONORS=()
ARMS=()
# Default is the locked (Membrane, Cell). Markers added after the compartment sample was extracted have
# only whole-cell values and need "Cell" alone -- see restore_validation.load_validation_sample.
FEATURE_COMPARTMENTS=(Membrane Cell)
while [[ $# -gt 0 ]]; do
  case "$1" in
    --pair-set) PAIR_SET="$2"; shift 2 ;;
    --donor)    DONORS+=("$2"); shift 2 ;;
    --arm)      ARMS+=("$2"); shift 2 ;;
    --feature-compartments) shift; FEATURE_COMPARTMENTS=()
        while [[ $# -gt 0 && "$1" != --* ]]; do FEATURE_COMPARTMENTS+=("$1"); shift; done ;;
    *) echo "unknown argument: $1" >&2; usage 1 ;;
  esac
done

# Default arm matrix. reference_and_target/max reproduces v9 (the control was additionally capped at
# the target's own Otsu floor, which made the divisor BE that floor); the reference_only rows walk the
# divisor statistic from the manuscript maximum down through the robust quantiles.
if [[ ${#ARMS[@]} -eq 0 ]]; then
  ARMS=(reference_and_target/max reference_only/max reference_only/p999 reference_only/p99 reference_only/p95)
fi
# The 7 Gate-3 pilot donors (6539 and 6476 by standing decision, plus the weakest case per pair).
if [[ ${#DONORS[@]} -eq 0 ]]; then
  DONORS=(115 6380 6436 6450 6476 6539 6591)
fi

SAMPLE="${REPO_ROOT}/data/restore_pair_validation/qupath_compartment_full_20donor"
RAW_CELLS="${REPO_ROOT}/data/cells"
REDSEA_CELLS="${REPO_ROOT}/data/cells_redsea"
PAIR_REVIEWS="${REPO_ROOT}/data/restore_pair_validation/pair_reviews_accepted.csv"
# The SAME expression reviews the production apply reads (config.restore_expression_reviews_csv).
# Omitting this does not merely lose provenance -- it changes states: without it donor 6566 INS returns
# MODEL_UNSTABLE instead of the NO_TARGET_EXPRESSION_CONFIRMED that production records, so a sweep arm
# would silently disagree with the run it is meant to explain.
EXPRESSION_REVIEWS="${REPO_ROOT}/data/restore_pair_validation/expression_reviews.csv"

for p in "$SAMPLE" "$RAW_CELLS" "$REDSEA_CELLS" "$PAIR_REVIEWS" "$EXPRESSION_REVIEWS"; do
  [[ -e "$p" ]] || { echo "[err] missing required input: $p" >&2; exit 1; }
done

DONOR_ARGS=(); for d in "${DONORS[@]}"; do DONOR_ARGS+=(--donor "$d"); done
mkdir -p "$OUT_ROOT"

echo "[sweep] pair-set=${PAIR_SET} donors=${DONORS[*]}"
echo "[sweep] arms=${ARMS[*]}"
echo "[sweep] output root=${OUT_ROOT}"

for arm in "${ARMS[@]}"; do
  ctrl="${arm%%/*}"; stat="${arm##*/}"
  out="${OUT_ROOT}/ctrl-${ctrl}_div-${stat}"
  if [[ -d "$out" ]]; then
    echo "[skip] ${arm} — ${out} already exists"
    continue
  fi
  # An aborted run leaves the module's atomic-write staging dir behind; it is not a result.
  rm -rf "${out}.tmp"
  echo "[run ] ${arm} -> ${out}"
  python -m phenocycler.restore_validation evaluate-locked \
    --sample "$SAMPLE" \
    --raw-cells "$RAW_CELLS" \
    --redsea-cells "$REDSEA_CELLS" \
    --pair-set "$PAIR_SET" \
    --pair-reviews "$PAIR_REVIEWS" \
    --expression-reviews "$EXPRESSION_REVIEWS" \
    --control-definition "$ctrl" \
    --divisor-statistic "$stat" \
    --feature-compartments "${FEATURE_COMPARTMENTS[@]}" \
    --output "$out" \
    "${DONOR_ARGS[@]}"
done

echo "[sweep] done. Arms under ${OUT_ROOT}:"
ls -1 "${OUT_ROOT}"
