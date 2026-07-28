#!/bin/bash
# Step 6 — PhenoCycler <-> Xenium integration.
#
# Runs in the phenocycler_integration environment (environment-integration.yml), which is the
# union of the core stack and the Xenium/spatial one — export_pheno reads hive parquet while
# import_xenium reads h5ad/zarr, so both have to be importable in one process.
#
# The whole thing is idempotent, so a re-submit after a partial run picks up where it stopped.
# Per-donor work inside register/match is not embarrassingly parallel the way REDSEA is (the
# expensive part is a handful of image registrations), so this is a single task rather than an
# array; use --donor to split it across jobs if the cohort grows.
#SBATCH --job-name=pheno_xen_integration
#SBATCH --output=logs/integration_%j.out
#SBATCH --error=logs/integration_%j.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=128gb
#SBATCH --time=12:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=your_email@ufl.edu
#SBATCH --qos=your_qos
#SBATCH --account=your_account

set -euo pipefail
module load conda 2>/dev/null || true
conda activate phenocycler_integration
mkdir -p logs

MODE="${MODE:-sequential}"
ROI="${ROI:-panc}"

# Always print the pairing table first — it is where donor/roi problems surface, and it costs
# nothing. Donors reported as `donor_unknown` need a row in data/integration/donor_overrides.csv.
python -m phenocycler.integration.manifest --report

python -m phenocycler.integration.pipeline --mode "$MODE" --roi "$ROI" "$@"

# The report is the thing to read afterwards: per-donor QC verdicts plus the cross-modal
# concordance summary.
cat data/integration/qc/integration_report.txt 2>/dev/null || true
