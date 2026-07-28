#!/bin/bash
# Quick-start setup for the Phenocycler_Analysis pipeline
# (raw data -> broad lineage via REDSEA + RESTORE).
set -e

echo "=========================================="
echo "Phenocycler_Analysis — setup"
echo "=========================================="

# 1) Vendored submodules: RESTORE (normalization) + XeniumPanelExplorer (the curated
#    identity_core / panel_roles CSVs the integration crosswalk is derived from).
echo "[1/4] Fetching vendored submodules (RESTORE, XeniumPanelExplorer) ..."
git submodule update --init --recursive

# 2) Conda environment
if ! command -v conda &> /dev/null; then
  echo "ERROR: conda not found. Install Miniconda/Anaconda first:"
  echo "  https://docs.conda.io/en/latest/miniconda.html"
  exit 1
fi
echo "[2/4] Creating conda env 'phenocycler_analysis' ..."
if conda env list | grep -q "^phenocycler_analysis "; then
  echo "  env already exists (skipping; 'conda env update -f environment.yml' to refresh)"
else
  conda env create -f environment.yml
fi

# 3) Verify the pipeline + vendored RESTORE import
echo "[3/4] Verifying installation ..."
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate phenocycler_analysis
python - <<'PY'
import phenocycler
from phenocycler import load_config
cfg = load_config()
print("  phenocycler", phenocycler.__version__)
for m in ("cells_parquet", "redsea", "restore", "lineage", "qupath_export", "pipeline"):
    __import__(f"phenocycler.{m}")
print("  all pipeline modules import OK")
import sys; sys.path.insert(0, str(cfg.restore_vendor))
try:
    from RESTORE import Normalization  # noqa
    print("  vendored RESTORE import OK")
except Exception as e:
    print(f"  NOTE: RESTORE import needs 'spams' (pip install spams-bin): {e}")
PY

# 4) Create the data directory skeleton
echo "[4/4] Creating data directories ..."
mkdir -p data/raw data/redsea_scratch/geojson logs

echo ""
echo "=========================================="
echo "Setup complete."
echo "  1. Edit config.ini [paths] to point at your qptiff images, Cellmeasurements.csv,"
echo "     and donor metadata (see DATA_README.md)."
echo "  2. conda activate phenocycler_analysis"
echo "  3. Run the pipeline:  python -m phenocycler.pipeline --status"
echo "     or open notebooks/00_run_full_pipeline.ipynb"
echo "  4. Run the tests:     pytest tests/"
echo ""
echo "  For PhenoCycler <-> Xenium integration (a separate, heavier environment):"
echo "     conda env create -f environment-integration.yml"
echo "     conda activate phenocycler_integration && pip install -e ."
echo "     python -m phenocycler.integration.manifest --report"
echo "=========================================="
