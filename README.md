# Phenocycler Spatial Proteomics Analysis Pipeline

A comprehensive, jupyter-notebook-based analysis pipeline for Phenocycler spatial proteomics data. Optimized for HiPerGator HPC cluster and includes scvi-tools, scanpy, squidpy, and R integration.

## Table of Contents
- [Features](#features)
- [Installation](#installation)
- [Quick Start](#quick-start)
- [Workflow](#workflow)
- [HiPerGator Usage](#hipergator-usage)
- [Directory Structure](#directory-structure)
- [Notebooks](#notebooks)
- [Scripts](#scripts)
- [Citation](#citation)

## Features

### Analysis Capabilities
- ✅ **Preprocessing**: QC, filtering, normalization from qptiff or CSV files (QuPath export)
- ✅ **Phenotyping**: Cell type annotation with marker proteins and scVI
- ✅ **Spatial Analysis**: Neighborhood enrichment, cell-cell interactions, spatial statistics
- ✅ **Group Comparisons**: Differential expression and composition analysis
- ✅ **Tissue Comparisons**: Multi-sample integration and comparison

### Technologies
- **Python**: scanpy, squidpy, scvi-tools, tifffile
- **R**: Seurat, ggplot2, rpy2 integration
- **HPC**: SLURM scripts for HiPerGator
- **Visualization**: matplotlib, seaborn, plotly, napari

## Installation

### 1. Clone Repository

```bash
git clone https://github.com/smith6jt-cop/Phenocycler_Analysis.git
cd Phenocycler_Analysis
```

### 2. Create Conda Environment

```bash
conda env create -f environment.yml
conda activate phenocycler_analysis
```

This creates an environment with all required packages including:
- Python 3.10
- scanpy, squidpy, scvi-tools
- R and Seurat
- Jupyter Lab
- All dependencies

### 3. Verify Installation

```bash
python -c "import scanpy as sc; import squidpy as sq; import scvi; print('All packages loaded successfully!')"
```

## Quick Start

### 1. Prepare Your Data

Place your Phenocycler data in the `data/raw/` directory:

**Option 1: qptiff files**
```bash
data/raw/
├── phenocycler_sample_01.qptiff
├── phenocycler_sample_02.qptiff
└── ...
```

**Option 2: CSV files from QuPath**
```bash
data/raw/
├── sample_01_expression.csv
├── sample_01_metadata.csv
├── sample_01_coordinates.csv
└── ...
```

### 2. Run Analysis Locally

```bash
# Activate environment
conda activate phenocycler_analysis

# Launch Jupyter Lab
jupyter lab

# Open and run notebooks in order:
# 01_preprocessing.ipynb
# 02_phenotyping.ipynb
# 03_spatial_analysis.ipynb
# ...
```

### 3. Run on HiPerGator

```bash
# Submit individual jobs
sbatch scripts/slurm/01_run_preprocessing.sh

# Or run complete pipeline
sbatch scripts/slurm/run_full_pipeline.sh
```

## Workflow

The analysis pipeline consists of 5 sequential notebooks:

```
Raw qptiff/CSV files
     ↓
[01_preprocessing.ipynb]
     ↓
Filtered & normalized data
     ↓
[02_phenotyping.ipynb]
     ↓
Annotated cell types
     ↓
[03_spatial_analysis.ipynb]
     ↓
Spatial metrics & interactions
     ↓
[04_group_comparisons.ipynb]
     ↓
Differential expression
     ↓
[05_tissue_comparisons.ipynb]
     ↓
Final integrated dataset
```

## HiPerGator Usage

### Setup on HiPerGator

1. **Load conda module**:
```bash
module load conda
```

2. **Create environment**:
```bash
conda env create -f environment.yml
```

3. **Update SLURM scripts**:

Edit the following in each SLURM script:
- `--mail-user`: Your UFL email
- `--qos`: Your QOS allocation
- `--account`: Your account/group

### Submit Jobs

**Individual notebooks**:
```bash
# Create logs directory
mkdir -p logs

# Submit jobs
sbatch scripts/slurm/01_run_preprocessing.sh
sbatch scripts/slurm/02_run_phenotyping.sh
sbatch scripts/slurm/03_run_spatial_analysis.sh
```

**Full pipeline**:
```bash
sbatch scripts/slurm/run_full_pipeline.sh
```

### Monitor Jobs

```bash
# Check job status
squeue -u $USER

# Check job output
tail -f logs/preprocess_JOBID.out

# Cancel job
scancel JOBID
```

### Resource Requirements

| Step | CPUs | Memory | Time | Notes |
|------|------|--------|------|-------|
| Preprocessing | 8 | 64GB | 24h | Data loading and QC |
| Phenotyping | 8 | 64GB | 24h | Includes scVI training |
| Spatial Analysis | 16 | 128GB | 48h | Computationally intensive |
| Full Pipeline | 16 | 128GB | 96h | All steps sequentially |

## Directory Structure

```
Phenocycler_Analysis/
├── README.md                          # This file
├── environment.yml                    # Conda environment
├── notebooks/                         # Analysis notebooks
│   ├── 01_preprocessing.ipynb
│   ├── 02_phenotyping.ipynb
│   ├── 03_spatial_analysis.ipynb
│   ├── 04_group_comparisons.ipynb
│   └── 05_tissue_comparisons.ipynb
├── scripts/
│   ├── R/                            # R analysis scripts
│   │   └── phenocycler_analysis.R
│   └── slurm/                        # HiPerGator SLURM scripts
│       ├── 01_run_preprocessing.sh
│       ├── 02_run_phenotyping.sh
│       ├── 03_run_spatial_analysis.sh
│       └── run_full_pipeline.sh
├── utils/                            # Utility functions
│   └── analysis_utils.py
├── data/
│   ├── raw/                          # Raw qptiff/CSV files (not tracked)
│   ├── processed/                    # Processed data (not tracked)
│   └── exports/                      # Exported results (not tracked)
└── figures/                          # Generated figures (not tracked)
    ├── 01_preprocessing/
    ├── 02_phenotyping/
    ├── 03_spatial_analysis/
    ├── 04_group_comparisons/
    └── 05_tissue_comparisons/
```

## Notebooks

### 01_preprocessing.ipynb
- Load Phenocycler data from qptiff or CSV files (QuPath export)
- Quality control metrics
- Cell and marker filtering
- Normalization and transformation
- Dimensionality reduction (PCA, UMAP)
- Initial clustering

**Input**: Raw qptiff files or CSV files from QuPath  
**Output**: `*_preprocessed.h5ad`

### 02_phenotyping.ipynb
- Marker protein analysis
- Cell type scoring
- Automated annotation
- scVI-based clustering
- Manual curation support
- Final cell type assignment

**Input**: `*_preprocessed.h5ad`  
**Output**: `*_annotated.h5ad`

### 03_spatial_analysis.ipynb
- Spatial neighborhood graph
- Neighborhood enrichment
- Co-occurrence analysis
- Spatial autocorrelation (Moran's I)
- Spatial domains
- Ligand-receptor interactions
- Ripley's statistics

**Input**: `*_annotated.h5ad`  
**Output**: `*_spatial_analysis.h5ad`

### 04_group_comparisons.ipynb
- Cell type composition analysis
- Differential expression between groups
- Cell type-specific DE
- Volcano plots
- Statistical testing

**Input**: `*_spatial_analysis.h5ad`  
**Output**: DE results, composition tables

### 05_tissue_comparisons.ipynb
- Multi-tissue integration
- Batch correction with scVI
- Cross-tissue differential expression
- Tissue-specific spatial patterns
- Comparative analysis

**Input**: Multiple `*_annotated.h5ad` files  
**Output**: `integrated_tissues.h5ad`

## Scripts

### R Integration (`scripts/R/phenocycler_analysis.R`)

Advanced statistical analysis and visualization:

```r
source("scripts/R/phenocycler_analysis.R")

# Load data
seurat_obj <- load_phenocycler_h5ad("data/processed/sample_annotated.h5ad")

# Differential expression
de_results <- run_de_analysis(seurat_obj, "condition", "Treatment", "Control")

# Visualizations
plot_volcano(de_results)
plot_spatial_features(seurat_obj, c("CD3", "CD8"))
plot_composition(seurat_obj, "celltype", "condition")
```

### SLURM Scripts

All SLURM scripts support:
- Email notifications
- Resource allocation
- Error logging
- Automatic output organization

Modify these variables in each script:
```bash
#SBATCH --mail-user=your_email@ufl.edu
#SBATCH --qos=your_qos
#SBATCH --account=your_account
```

## Utility Functions

The `utils/analysis_utils.py` module provides reusable functions:

```python
from utils.analysis_utils import *

# Load data from qptiff
adata = load_qptiff_data(qptiff_path, "sample_01")

# Load data from QuPath CSV
adata = load_qupath_csv(csv_dir, "sample_01")

# Calculate QC metrics
adata = calculate_qc_metrics(adata)

# Filter
adata = filter_cells_and_markers(adata)

# Normalize and identify highly variable proteins
adata = normalize_and_hvg(adata)

# Run standard workflow
adata = run_standard_workflow(adata)

# Export results
export_to_csv(adata, output_dir, "sample_01")
```

## Customization

### Cell Type Markers

Edit marker proteins in `02_phenotyping.ipynb`:

```python
marker_proteins = {
    'T cells': ['CD3', 'CD8', 'CD4'],
    'B cells': ['CD19', 'CD20'],
    'Macrophages': ['CD68', 'CD163'],
    'Epithelial': ['PanCK', 'EpCAM'],
    # Add your tissue-specific markers
}
```

### QC Thresholds

Adjust in `01_preprocessing.ipynb`:

```python
MIN_MARKERS = 5
MIN_CELLS = 3
MAX_BACKGROUND = 0.2
MIN_INTENSITY = 10
```

### Spatial Parameters

Modify in `03_spatial_analysis.ipynb`:

```python
# Spatial neighbors
sq.gr.spatial_neighbors(adata, n_neighs=10)

# Spatial domains
n_domains = 5
```

## Troubleshooting

### Common Issues

**1. Memory errors on HiPerGator**
```bash
# Increase memory allocation
#SBATCH --mem=256gb
```

**2. SLURM job fails immediately**
```bash
# Check account and QOS
sacctmgr show assoc user=$USER
```

**3. Conda environment issues**
```bash
# Remove and recreate
conda env remove -n phenocycler_analysis
conda env create -f environment.yml
```

**4. Module not found errors**
```bash
# Verify environment is activated
conda activate phenocycler_analysis
which python
```

## Best Practices

1. **Start with small test dataset** to verify pipeline
2. **Save intermediate results** after each major step
3. **Use version control** for custom modifications
4. **Document parameters** used for reproducibility
5. **Keep raw data separate** from processed outputs
6. **Regular backups** of analysis results

## Data Organization

Recommended file naming:
```
<project>_<tissue>_<condition>_<replicate>_<step>.h5ad

Examples:
study1_liver_control_rep1_preprocessed.h5ad
study1_liver_treatment_rep1_annotated.h5ad
```

## Citation

If you use this pipeline, please cite:

- **Scanpy**: Wolf et al., Genome Biology (2018)
- **Squidpy**: Palla et al., Nature Methods (2022)
- **scvi-tools**: Lopez et al., Nature Methods (2018)
- **Seurat**: Hao et al., Cell (2021)

## Support

For questions and issues:
- Open an issue on GitHub
- Check HiPerGator documentation: https://help.rc.ufl.edu/
- Scanpy tutorials: https://scanpy-tutorials.readthedocs.io/

## License

This pipeline is provided as-is for research purposes.

## Acknowledgments

- University of Florida Research Computing
- HiPerGator HPC cluster
- Scanpy, Squidpy, and scvi-tools development teams
