"""
Phenocycler Analysis Utility Functions

This module provides reusable functions for Phenocycler spatial proteomics data analysis.
Includes data loading, QC, preprocessing, and export utilities.
"""

import numpy as np
import pandas as pd
import scanpy as sc
import squidpy as sq
import anndata as ad
from pathlib import Path
from typing import Optional, Union, List, Dict
import warnings
import tifffile
import h5py
from scipy import sparse

# Suppress warnings
warnings.filterwarnings('ignore')


def load_qptiff_data(
    qptiff_path: Union[str, Path],
    sample_name: str,
    channel_names: Optional[List[str]] = None
) -> ad.AnnData:
    """
    Load Phenocycler data from qptiff file.
    
    NOTE: This function requires custom implementation for your specific qptiff format.
    We recommend using QuPath to export cell measurements as CSV instead.
    
    To implement qptiff loading:
    1. Load multi-channel TIFF with tifffile
    2. Perform cell segmentation (e.g., with cellpose, stardist, or QuPath)
    3. Extract mean fluorescence intensity (MFI) per channel per cell
    4. Extract cell centroid coordinates
    5. Create AnnData object with expression matrix and spatial coords
    
    Example implementation:
        from cellpose import models
        # Segment cells
        model = models.Cellpose(model_type='cyto')
        masks, flows, styles = model.eval(images[0])  # Use DAPI/nuclear channel
        
        # Extract intensities per cell per channel
        from skimage.measure import regionprops_table
        for channel_idx, channel_name in enumerate(channel_names):
            props = regionprops_table(masks, images[channel_idx], 
                                     properties=['mean_intensity', 'centroid'])
    
    Parameters
    ----------
    qptiff_path : str or Path
        Path to qptiff file
    sample_name : str
        Name of the sample
    channel_names : list of str, optional
        Names of protein channels
    
    Returns
    -------
    adata : AnnData
        AnnData object with expression data and spatial coordinates
    """
    print(f"Loading qptiff file: {qptiff_path}")
    
    # Load TIFF file
    with tifffile.TiffFile(qptiff_path) as tif:
        # Read image data
        images = tif.asarray()
        
        # Get metadata if available
        if tif.pages[0].description:
            print(f"Image metadata: {tif.pages[0].description[:100]}...")
    
    raise NotImplementedError(
        "qptiff loading requires custom cell segmentation implementation.\n\n"
        "RECOMMENDED APPROACH: Use QuPath to export cell measurements as CSV:\n"
        "  1. Open qptiff in QuPath (https://qupath.github.io/)\n"
        "  2. Run: Analyze → Cell detection → Cell detection\n"
        "  3. Measure: Measure → Show detection measurements\n"
        "  4. Export: File → Export → Export as CSV\n"
        "  5. Use load_qupath_csv() to load the exported data\n\n"
        "ALTERNATIVE: Implement custom segmentation (see docstring for example)\n"
        "For help, see DATA_README.md section 'Converting qptiff to CSV'"
    )


def load_qupath_csv(
    csv_dir: Union[str, Path],
    sample_name: str,
    expression_file: str = "expression.csv",
    metadata_file: Optional[str] = "metadata.csv",
    coordinates_file: Optional[str] = "coordinates.csv"
) -> ad.AnnData:
    """
    Load Phenocycler data from QuPath CSV export.
    
    QuPath export should contain:
    - expression.csv: Cell x protein expression matrix
    - coordinates.csv: Cell spatial coordinates (x, y)
    - metadata.csv (optional): Additional cell metadata
    
    Parameters
    ----------
    csv_dir : str or Path
        Directory containing CSV files
    sample_name : str
        Name of the sample
    expression_file : str
        Name of expression CSV file
    metadata_file : str, optional
        Name of metadata CSV file
    coordinates_file : str, optional
        Name of coordinates CSV file
    
    Returns
    -------
    adata : AnnData
        AnnData object with expression data and spatial coordinates
    """
    csv_dir = Path(csv_dir)
    print(f"Loading QuPath CSV data from: {csv_dir}")
    
    # Load expression data
    expr_path = csv_dir / expression_file
    if not expr_path.exists():
        raise FileNotFoundError(f"Expression file not found: {expr_path}")
    
    expr_df = pd.read_csv(expr_path, index_col=0)
    print(f"Loaded expression matrix: {expr_df.shape}")
    
    # Create AnnData object
    adata = ad.AnnData(X=expr_df.values)
    adata.obs_names = expr_df.index.astype(str)
    adata.var_names = expr_df.columns
    adata.obs['sample'] = sample_name
    
    # Load spatial coordinates
    if coordinates_file:
        coords_path = csv_dir / coordinates_file
        if coords_path.exists():
            coords_df = pd.read_csv(coords_path, index_col=0)
            # Ensure same order as expression data
            coords_df = coords_df.loc[adata.obs_names]
            adata.obsm['spatial'] = coords_df[['x', 'y']].values
            print(f"Added spatial coordinates: {adata.obsm['spatial'].shape}")
        else:
            print(f"Warning: Coordinates file not found: {coords_path}")
    
    # Load metadata
    if metadata_file:
        meta_path = csv_dir / metadata_file
        if meta_path.exists():
            meta_df = pd.read_csv(meta_path, index_col=0)
            # Ensure same order
            meta_df = meta_df.loc[adata.obs_names]
            # Add metadata columns to obs
            for col in meta_df.columns:
                adata.obs[col] = meta_df[col].values
            print(f"Added metadata: {list(meta_df.columns)}")
        else:
            print(f"Warning: Metadata file not found: {meta_path}")
    
    return adata


def calculate_qc_metrics(
    adata: ad.AnnData,
    background_threshold: float = 0.1
) -> ad.AnnData:
    """
    Calculate quality control metrics for Phenocycler data.
    
    Parameters
    ----------
    adata : AnnData
        AnnData object
    background_threshold : float
        Threshold for identifying background markers
    
    Returns
    -------
    adata : AnnData
        Updated AnnData with QC metrics
    """
    print("Calculating QC metrics...")
    
    # Calculate basic metrics
    adata.obs['n_markers'] = (adata.X > 0).sum(axis=1)
    adata.obs['total_intensity'] = adata.X.sum(axis=1)
    adata.obs['mean_intensity'] = adata.X.mean(axis=1)
    adata.obs['max_intensity'] = adata.X.max(axis=1)
    
    # Marker-level metrics
    adata.var['n_cells'] = (adata.X > 0).sum(axis=0).A1 if sparse.issparse(adata.X) else (adata.X > 0).sum(axis=0)
    adata.var['mean_expression'] = adata.X.mean(axis=0).A1 if sparse.issparse(adata.X) else adata.X.mean(axis=0)
    adata.var['max_expression'] = adata.X.max(axis=0).A1 if sparse.issparse(adata.X) else adata.X.max(axis=0)
    
    # Identify potential background markers
    adata.var['is_background'] = adata.var['mean_expression'] < background_threshold
    
    print(f"QC metrics calculated:")
    print(f"  - Cells: {adata.n_obs}")
    print(f"  - Markers: {adata.n_vars}")
    print(f"  - Mean markers/cell: {adata.obs['n_markers'].mean():.1f}")
    print(f"  - Background markers: {adata.var['is_background'].sum()}")
    
    return adata


def filter_cells_and_markers(
    adata: ad.AnnData,
    min_markers: int = 5,
    min_cells: int = 3,
    min_intensity: float = 10,
    max_background_ratio: float = 0.5
) -> ad.AnnData:
    """
    Filter low-quality cells and markers.
    
    Parameters
    ----------
    adata : AnnData
        AnnData object
    min_markers : int
        Minimum number of markers per cell
    min_cells : int
        Minimum number of cells per marker
    min_intensity : float
        Minimum total intensity per cell
    max_background_ratio : float
        Maximum fraction of background markers
    
    Returns
    -------
    adata : AnnData
        Filtered AnnData object
    """
    print(f"Filtering cells and markers...")
    print(f"Starting with {adata.n_obs} cells and {adata.n_vars} markers")
    
    # Filter cells
    cell_filter = (
        (adata.obs['n_markers'] >= min_markers) &
        (adata.obs['total_intensity'] >= min_intensity)
    )
    print(f"Removing {(~cell_filter).sum()} low-quality cells")
    adata = adata[cell_filter, :].copy()
    
    # Filter markers
    marker_filter = adata.var['n_cells'] >= min_cells
    print(f"Removing {(~marker_filter).sum()} low-count markers")
    adata = adata[:, marker_filter].copy()
    
    # Optionally remove background markers
    if 'is_background' in adata.var.columns:
        background_markers = adata.var['is_background'].sum()
        if background_markers / adata.n_vars > max_background_ratio:
            print(f"Warning: {background_markers} background markers detected")
    
    print(f"After filtering: {adata.n_obs} cells and {adata.n_vars} markers")
    
    return adata


def normalize_and_hvg(
    adata: ad.AnnData,
    method: str = 'zscore',
    n_top_markers: int = 30,
    target_sum: float = 10000
) -> ad.AnnData:
    """
    Normalize data and identify highly variable markers.
    
    Parameters
    ----------
    adata : AnnData
        AnnData object
    method : str
        Normalization method: 'zscore', 'log', or 'arcsinh'
    n_top_markers : int
        Number of highly variable markers to identify
    target_sum : float
        Target sum for normalization (if using log method)
    
    Returns
    -------
    adata : AnnData
        Normalized AnnData object
    """
    print(f"Normalizing data using {method} method...")
    
    # Store raw data
    adata.raw = adata.copy()
    
    if method == 'zscore':
        # Z-score normalization (standard for protein data)
        from sklearn.preprocessing import StandardScaler
        scaler = StandardScaler()
        adata.X = scaler.fit_transform(adata.X)
    
    elif method == 'log':
        # Log normalization (similar to RNA-seq)
        sc.pp.normalize_total(adata, target_sum=target_sum)
        sc.pp.log1p(adata)
    
    elif method == 'arcsinh':
        # Arcsinh transformation (common for CyTOF/CODEX)
        cofactor = 5
        adata.X = np.arcsinh(adata.X / cofactor)
    
    else:
        raise ValueError(f"Unknown normalization method: {method}")
    
    # Identify highly variable markers
    if n_top_markers and n_top_markers < adata.n_vars:
        sc.pp.highly_variable_genes(
            adata,
            n_top_genes=n_top_markers,
            flavor='seurat_v3'
        )
        print(f"Identified {adata.var['highly_variable'].sum()} highly variable markers")
    
    return adata


def run_standard_workflow(
    adata: ad.AnnData,
    n_pcs: int = 30,
    n_neighbors: int = 15,
    resolution: float = 1.0
) -> ad.AnnData:
    """
    Run standard dimensionality reduction and clustering workflow.
    
    Parameters
    ----------
    adata : AnnData
        Normalized AnnData object
    n_pcs : int
        Number of principal components
    n_neighbors : int
        Number of neighbors for UMAP
    resolution : float
        Clustering resolution
    
    Returns
    -------
    adata : AnnData
        AnnData with PCA, UMAP, and clusters
    """
    print("Running standard workflow...")
    
    # PCA
    print(f"Computing PCA with {n_pcs} components...")
    sc.tl.pca(adata, n_comps=n_pcs, svd_solver='arpack')
    
    # Neighborhood graph
    print(f"Computing neighborhood graph with {n_neighbors} neighbors...")
    sc.pp.neighbors(adata, n_neighbors=n_neighbors, n_pcs=n_pcs)
    
    # UMAP
    print("Computing UMAP...")
    sc.tl.umap(adata)
    
    # Leiden clustering
    print(f"Running Leiden clustering with resolution {resolution}...")
    sc.tl.leiden(adata, resolution=resolution)
    
    print("Standard workflow complete!")
    
    return adata


def export_to_csv(
    adata: ad.AnnData,
    output_dir: Union[str, Path],
    sample_name: str
) -> None:
    """
    Export AnnData to CSV files.
    
    Parameters
    ----------
    adata : AnnData
        AnnData object to export
    output_dir : str or Path
        Output directory
    sample_name : str
        Sample name for file naming
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    print(f"Exporting data to {output_dir}...")
    
    # Export expression matrix
    expr_df = pd.DataFrame(
        adata.X if not sparse.issparse(adata.X) else adata.X.toarray(),
        index=adata.obs_names,
        columns=adata.var_names
    )
    expr_path = output_dir / f"{sample_name}_expression.csv"
    expr_df.to_csv(expr_path)
    print(f"Exported expression: {expr_path}")
    
    # Export metadata
    meta_path = output_dir / f"{sample_name}_metadata.csv"
    adata.obs.to_csv(meta_path)
    print(f"Exported metadata: {meta_path}")
    
    # Export spatial coordinates if available
    if 'spatial' in adata.obsm:
        coords_df = pd.DataFrame(
            adata.obsm['spatial'],
            index=adata.obs_names,
            columns=['x', 'y']
        )
        coords_path = output_dir / f"{sample_name}_coordinates.csv"
        coords_df.to_csv(coords_path)
        print(f"Exported coordinates: {coords_path}")
    
    # Export marker info
    marker_path = output_dir / f"{sample_name}_marker_info.csv"
    adata.var.to_csv(marker_path)
    print(f"Exported marker info: {marker_path}")


def create_summary_report(
    adata: ad.AnnData,
    sample_name: str,
    output_path: Union[str, Path]
) -> pd.DataFrame:
    """
    Create a summary report of the dataset.
    
    Parameters
    ----------
    adata : AnnData
        AnnData object
    sample_name : str
        Sample name
    output_path : str or Path
        Path to save summary CSV
    
    Returns
    -------
    summary : DataFrame
        Summary statistics
    """
    summary_data = {
        'Sample': [sample_name],
        'n_cells': [adata.n_obs],
        'n_markers': [adata.n_vars],
        'mean_markers_per_cell': [adata.obs['n_markers'].mean()],
        'median_markers_per_cell': [adata.obs['n_markers'].median()],
        'mean_total_intensity': [adata.obs['total_intensity'].mean()],
        'has_spatial_coords': ['spatial' in adata.obsm],
        'has_clusters': ['leiden' in adata.obs.columns or 'louvain' in adata.obs.columns],
    }
    
    # Add cell type info if available
    if 'celltype' in adata.obs.columns:
        summary_data['n_celltypes'] = [adata.obs['celltype'].nunique()]
    
    summary_df = pd.DataFrame(summary_data)
    
    # Save to file
    summary_df.to_csv(output_path, index=False)
    print(f"Summary report saved to: {output_path}")
    
    return summary_df
