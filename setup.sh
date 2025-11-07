#!/bin/bash
# Quick Start Script for Phenocycler Analysis Pipeline
# This script helps you get started with the Phenocycler analysis pipeline

set -e  # Exit on error

echo "=========================================="
echo "Phenocycler Analysis Pipeline - Quick Start"
echo "=========================================="
echo ""

# Check if conda is available
if ! command -v conda &> /dev/null; then
    echo "ERROR: conda is not installed or not in PATH"
    echo "Please install Miniconda or Anaconda first"
    echo "Visit: https://docs.conda.io/en/latest/miniconda.html"
    exit 1
fi

echo "Step 1: Creating conda environment..."
echo "This may take several minutes..."
if conda env list | grep -q "^phenocycler_analysis "; then
    echo "Environment 'phenocycler_analysis' already exists."
    read -p "Do you want to recreate it? (y/n) " -n 1 -r
    echo
    if [[ $REPLY =~ ^[Yy]$ ]]; then
        conda env remove -n phenocycler_analysis -y
        conda env create -f environment.yml
    fi
else
    conda env create -f environment.yml
fi

echo ""
echo "Step 2: Activating environment..."
source $(conda info --base)/etc/profile.d/conda.sh
conda activate phenocycler_analysis

echo ""
echo "Step 3: Verifying installation..."
python << 'PYEOF'
import sys
print(f"Python version: {sys.version}")

try:
    import scanpy as sc
    print(f"✓ scanpy: {sc.__version__}")
except ImportError as e:
    print(f"✗ scanpy: {e}")

try:
    import squidpy as sq
    print(f"✓ squidpy: {sq.__version__}")
except ImportError as e:
    print(f"✗ squidpy: {e}")

try:
    import scvi
    print(f"✓ scvi-tools: {scvi.__version__}")
except ImportError as e:
    print(f"✗ scvi-tools: {e}")

try:
    import numpy as np
    print(f"✓ numpy: {np.__version__}")
except ImportError as e:
    print(f"✗ numpy: {e}")

try:
    import pandas as pd
    print(f"✓ pandas: {pd.__version__}")
except ImportError as e:
    print(f"✗ pandas: {e}")

try:
    import tifffile
    print(f"✓ tifffile: {tifffile.__version__}")
except ImportError as e:
    print(f"✗ tifffile: {e}")

print("\nAll core packages loaded successfully!")
PYEOF

echo ""
echo "Step 4: Creating directory structure..."
mkdir -p data/raw
mkdir -p data/processed
mkdir -p data/exports
mkdir -p figures
mkdir -p logs
mkdir -p notebooks
mkdir -p scripts/R
mkdir -p scripts/slurm
mkdir -p utils

echo ""
echo "=========================================="
echo "Setup Complete!"
echo "=========================================="
echo ""
echo "Next steps:"
echo "1. Place your Phenocycler qptiff or CSV files in: data/raw/"
echo "2. Activate the environment:"
echo "   conda activate phenocycler_analysis"
echo "3. Launch Jupyter Lab:"
echo "   jupyter lab"
echo "4. Open and run: notebooks/01_preprocessing.ipynb"
echo ""
echo "For HiPerGator users:"
echo "1. Edit SLURM scripts in scripts/slurm/"
echo "2. Update email, QOS, and account settings"
echo "3. Submit jobs:"
echo "   sbatch scripts/slurm/01_run_preprocessing.sh"
echo ""
echo "For help, see README.md and DATA_README.md"
echo "=========================================="
