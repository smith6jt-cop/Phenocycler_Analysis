#!/bin/bash
# Create/update the focused environment and install the current package.
set -euo pipefail

REPOSITORY_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "${REPOSITORY_DIR}"
ENVIRONMENT_NAME="phenocycler-analysis"

if ! command -v conda >/dev/null 2>&1; then
  echo "ERROR: conda is required. Install Miniconda or Miniforge first." >&2
  exit 1
fi

source "$(conda info --base)/etc/profile.d/conda.sh"
if conda env list | awk '{print $1}' | grep -Fxq "${ENVIRONMENT_NAME}"; then
  echo "[setup] updating ${ENVIRONMENT_NAME}"
  conda env update --name "${ENVIRONMENT_NAME}" --file environment.yml
else
  echo "[setup] creating ${ENVIRONMENT_NAME}"
  conda env create --file environment.yml
fi

conda activate "${ENVIRONMENT_NAME}"
python -m pip install --editable ".[test]"

echo "[setup] verifying package data and operational imports"
python - <<'PY'
from importlib.resources import files

import phenocycler
from phenocycler.hierarchical_typing import TypingRegistry
from phenocycler.marker_registry import load_registry
from phenocycler.pipeline import STAGE_ORDER

package = files("phenocycler")
for name in ("marker_registry.json", "typing_rules.json"):
    resource = package.joinpath(name)
    if not resource.is_file():
        raise SystemExit(f"missing installed package data: {name}")
registry = load_registry()
TypingRegistry.from_marker_registry(registry)
print(f"  phenocycler {phenocycler.__version__}")
print(f"  registry {registry.registry_version} ({len(registry.markers)} markers)")
print(f"  stages: {' -> '.join(STAGE_ORDER)}")
PY

mkdir -p data logs

echo
echo "Setup complete."
echo "  conda activate ${ENVIRONMENT_NAME}"
echo "  phenocycler manifest template --out data/qupath_manifest_spec.json"
echo "  # Edit the specification, then:"
echo "  phenocycler manifest create --spec data/qupath_manifest_spec.json --out data/qupath_manifest.json"
echo "  phenocycler run --config config.ini"
echo "  phenocycler status --config config.ini"
echo "  python -m pytest tests"
