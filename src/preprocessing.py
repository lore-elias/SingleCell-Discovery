"""
Preprocessing module for single-cell RNA-seq data.

Functions for data loading, quality control, normalization, and gene selection.
"""

from typing import Optional
import scanpy as sc
import numpy as np
import pandas as pd
from anndata import AnnData

def load_h5ad(file_path: str) -> AnnData:
    """Load an H5AD file into an AnnData object.
    
    Parameters
    ----------
    file_path : str
        Path to the H5AD file.
        
    Returns
    -------
    AnnData
        Annotated data matrix.
        
    Raises
    ------
    FileNotFoundError
        If the file does not exist.
    """
    try:
        adata = sc.read_h5ad(file_path)
        print(f"Loaded data: {adata.shape[0]} cells, {adata.shape[1]} genes")
        return adata
    except FileNotFoundError:
        raise FileNotFoundError(f"H5AD file not found: {file_path}")
    except Exception as e:
        raise RuntimeError(f"Error loading H5AD file: {e}")

def filter_cells(
    adata: AnnData,
    min_genes: int = 200,
    max_genes: int = 2500,
    max_mito: float = 0.05
) -> AnnData:
    """Filter cells based on quality control metrics.
    
    Parameters
    ----------
    adata : AnnData
        Annotated data matrix.
    min_genes : int, optional
        Minimum number of genes detected per cell (default: 200).
    max_genes : int, optional
        Maximum number of genes detected per cell (default: 2500).
    max_mito : float, optional
        Maximum mitochondrial percentage (default: 0.05).
        
    Returns
    -------
    AnnData
        Filtered annotated data matrix.
    """
    # Calculate QC metrics
    adata.obs['n_genes'] = (adata.X > 0).sum(axis=1).A1 if hasattr(adata.X, 'A1') else (adata.X > 0).sum(axis=1)
    adata.obs['n_counts'] = np.asarray(adata.X.sum(axis=1)).flatten()
    
    # Calculate mitochondrial percentage
    mito_genes = adata.var_names.str.startswith('MT-')
    if mito_genes.sum() > 0:
        mito_counts = np.asarray(adata[:, mito_genes].X.sum(axis=1)).flatten()
        adata.obs['mito_percent'] = mito_counts / adata.obs['n_counts']
    else:
        adata.obs['mito_percent'] = 0
    
    # Filter cells
    initial_cells = adata.shape[0]
    sc.pp.filter_cells(adata, min_genes=min_genes)
    sc.pp.filter_cells(adata, max_genes=max_genes)
    adata = adata[adata.obs['mito_percent'] < max_mito]
    
    print(f"Filtered: {initial_cells} → {adata.shape[0]} cells")
    return adata

def normalize_and_log_transform(adata: AnnData, target_sum: float = 1e4) -> AnnData:
    """Normalize counts and log-transform.
    
    Parameters
    ----------
    adata : AnnData
        Annotated data matrix.
    target_sum : float, optional
        Target sum for normalization (default: 1e4).
        
    Returns
    -------
    AnnData
        Normalized and log-transformed annotated data matrix.
    """
    sc.pp.normalize_total(adata, target_sum=target_sum)
    sc.pp.log1p(adata)
    return adata

def select_variable_genes(adata: AnnData, n_top_genes: int = 2000) -> AnnData:
    """Select highly variable genes.
    
    Parameters
    ----------
    adata : AnnData
        Annotated data matrix.
    n_top_genes : int, optional
        Number of highly variable genes to select (default: 2000).
        
    Returns
    -------
    AnnData
        Annotated data matrix with only highly variable genes.
    """
    initial_genes = adata.shape[1]
    sc.pp.highly_variable_genes(adata, n_top_genes=n_top_genes)
    adata = adata[:, adata.var['highly_variable']]
    print(f"Gene selection: {initial_genes} → {adata.shape[1]} genes")
    return adata


def preprocess_pipeline(
    adata: AnnData,
    min_genes: int = 200,
    max_genes: int = 2500,
    max_mito: float = 0.05,
    n_top_genes: int = 2000,
    target_sum: float = 1e4
) -> AnnData:
    """Complete preprocessing pipeline.
    
    Parameters
    ----------
    adata : AnnData
        Annotated data matrix.
    min_genes : int, optional
        Minimum genes per cell (default: 200).
    max_genes : int, optional
        Maximum genes per cell (default: 2500).
    max_mito : float, optional
        Maximum mitochondrial percentage (default: 0.05).
    n_top_genes : int, optional
        Number of variable genes (default: 2000).
    target_sum : float, optional
        Normalization target sum (default: 1e4).
        
    Returns
    -------
    AnnData
        Preprocessed annotated data matrix.
    """
    print("Starting preprocessing pipeline...")
    adata = filter_cells(adata, min_genes, max_genes, max_mito)
    adata = normalize_and_log_transform(adata, target_sum)
    adata = select_variable_genes(adata, n_top_genes)
    print("Preprocessing complete!")
    return adata


