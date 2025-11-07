# Data Directory Structure

This document describes the expected data structure for the Phenocycler analysis pipeline.

## Directory Overview

```
data/
├── raw/              # Raw Phenocycler files (qptiff or CSV)
├── processed/        # Intermediate and final processed files
└── exports/          # Exported results and tables
```

## Raw Data (`data/raw/`)

### Phenocycler Data Formats

The pipeline supports two input formats:

**Option 1: qptiff files**
```
data/raw/
└── phenocycler_sample_01.qptiff
```

Note: qptiff loading requires custom implementation based on your specific format. We recommend using QuPath to export to CSV format.

**Option 2: CSV files from QuPath (Recommended)**
```
data/raw/
├── sample_01_expression.csv      # Cell x protein MFI matrix
├── sample_01_coordinates.csv     # Spatial coordinates (x, y)
└── sample_01_metadata.csv        # Optional cell metadata
```

### CSV Format Specifications

If providing CSV files from QuPath, use the following structure:

**expression.csv**
```csv
cell_id,CD3,CD4,CD8,CD19,CD68,PanCK,...
cell_0001,125.3,82.1,15.2,5.5,12.1,450.2,...
cell_0002,12.1,8.5,180.7,3.3,11.8,8.5,...
...
```

**coordinates.csv**
```csv
cell_id,x,y
cell_0001,1523.4,2847.6
cell_0002,1528.9,2851.2
...
```

**metadata.csv** (optional)
```csv
cell_id,cluster,region,sample_id
cell_0001,cluster_1,tumor,sample_01
cell_0002,cluster_2,stroma,sample_01
...
```

### Exporting from QuPath

To export data from QuPath:

1. Open your qptiff file in QuPath
2. Run cell detection/segmentation
3. Measure cell intensities for all channels
4. Export measurements as CSV:
   - File → Export → Export as CSV
5. Organize exported files according to the structure above

### File Naming Convention

Use descriptive, consistent naming:
```
<project>_<tissue>_<condition>_<replicate>_<type>.csv

Examples:
study1_liver_control_rep1_expression.csv
study1_liver_control_rep1_coordinates.csv
study1_liver_treatment_rep1_expression.csv
```

## Processed Data (`data/processed/`)

The pipeline automatically generates processed files at each step:

```
data/processed/
├── sample_01_preprocessed.h5ad          # After 01_preprocessing.ipynb
├── sample_01_annotated.h5ad             # After 02_phenotyping.ipynb
├── sample_01_spatial_analysis.h5ad      # After 03_spatial_analysis.ipynb
├── sample_01_celltypes.csv              # Cell type annotations
├── sample_01_cluster_markers.csv        # Marker proteins per cluster
├── sample_01_de_results.csv             # Differential expression results
└── integrated_tissues.h5ad              # After multi-tissue integration
```

## Exports (`data/exports/`)

Export directory for results:

```
data/exports/
├── sample_01_expression.csv       # Exported expression matrix
├── sample_01_metadata.csv         # Exported cell metadata
├── sample_01_coordinates.csv      # Exported spatial coordinates
└── sample_01_marker_info.csv      # Marker protein information
```

## Data Requirements

### Minimum Requirements

For Phenocycler data:
- Protein expression matrix (cells × proteins)
- Spatial coordinates (x, y) for each cell
- Cell IDs

### Recommended Additional Data

- Cell segmentation boundaries
- Channel/marker names
- Sample metadata (condition, replicate, etc.)
- QC metrics from imaging platform

## Converting qptiff to CSV

If you have qptiff files but need CSV format:

### Using QuPath

1. Install QuPath: https://qupath.github.io/
2. Open your qptiff file
3. Run cell detection: Analyze → Cell detection → Cell detection
4. Measure intensities: Measure → Show detection measurements
5. Export: File → Export → Export as CSV

### Using Python (Custom)

```python
import tifffile
import pandas as pd
import numpy as np

# Load qptiff
with tifffile.TiffFile('sample.qptiff') as tif:
    images = tif.asarray()
    
# Perform cell segmentation (using your preferred method)
# Extract cell-level features
# Create CSV files

# Example structure
expression_df.to_csv('sample_expression.csv')
coordinates_df.to_csv('sample_coordinates.csv')
```

## Storage Considerations

### Disk Space Requirements

Typical file sizes:
- Raw qptiff: 500MB - 5GB per sample
- CSV exports: 50MB - 500MB per sample
- Preprocessed h5ad: 100MB - 1GB per sample
- Figures: 10-50 MB per sample
- Total per sample: ~1-7 GB

### Backup Strategy

Important data to backup:
1. ✅ Raw data files (always keep originals)
2. ✅ Final processed files (annotated, integrated)
3. ✅ Analysis reports and figures
4. ⚠️ Intermediate files (can be regenerated)

## Data Privacy and Compliance

### For Human Data

- Ensure proper IRB approval
- De-identify patient data
- Follow HIPAA/GDPR guidelines
- Secure data storage
- Document data sharing agreements

### For Collaborative Projects

- Use consistent naming conventions
- Document metadata thoroughly
- Version control analysis parameters
- Share preprocessing scripts
- Track data provenance

## Troubleshooting

### Common Issues

**1. File not found errors**
```bash
# Check file exists
ls -lh data/raw/

# Check file permissions
chmod 644 data/raw/your_file.csv
```

**2. CSV format errors**
```python
# Check CSV structure
import pandas as pd
df = pd.read_csv('data/raw/expression.csv', index_col=0)
print(df.shape)
print(df.head())
```

**3. Memory errors with large files**
```python
# Use chunked reading
for chunk in pd.read_csv('large_file.csv', chunksize=10000):
    process(chunk)
```

**4. Coordinate mismatch**
```python
# Verify cell IDs match
expr = pd.read_csv('expression.csv', index_col=0)
coords = pd.read_csv('coordinates.csv', index_col=0)
assert expr.index.equals(coords.index), "Cell IDs don't match!"
```

## Example Datasets

For testing the pipeline, you can use:

1. **CODEX Mouse Spleen Data** (public)
   - https://data.mendeley.com/datasets/zjnpwh8m5b/1
   - Similar spatial proteomics data

2. **CyCIF Tissue Atlas** (public)
   - https://www.tissue-atlas.org/
   - High-plex imaging data

3. **Generate Synthetic Data** (for testing)

```python
import numpy as np
import pandas as pd

# Create synthetic data for testing
n_cells = 10000
n_markers = 30

# Random expression
expr = np.random.lognormal(3, 1, (n_cells, n_markers))
expr_df = pd.DataFrame(
    expr,
    index=[f'cell_{i:05d}' for i in range(n_cells)],
    columns=[f'Marker_{i}' for i in range(n_markers)]
)

# Random coordinates
coords_df = pd.DataFrame({
    'x': np.random.uniform(0, 1000, n_cells),
    'y': np.random.uniform(0, 1000, n_cells)
}, index=expr_df.index)

# Save
expr_df.to_csv('data/raw/test_expression.csv')
coords_df.to_csv('data/raw/test_coordinates.csv')
```

## Contact

For data format questions:
- Check README.md for general documentation
- See example notebooks for data loading patterns
- Open an issue on GitHub
