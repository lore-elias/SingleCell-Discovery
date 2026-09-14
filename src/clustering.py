"""
Clustering module for single-cell RNA-seq data.

Functions for PCA reduction, UMAP embedding, and Leiden clustering.
"""

from typing import Optional, Tuple
import scanpy as sc
from anndata import AnnData

def compute_pca(adata: AnnData, n_pcs: int = 50) -> AnnData:
    """Compute PCA on the data.
    
    Parameters
    ----------
    adata : AnnData
        Annotated data matrix.
    n_pcs : int, optional
        Number of principal components (default: 50).
        
    Returns
    -------
    AnnData
        Annotated data matrix with PCA coordinates in obsm['X_pca'].
    """
    sc.tl.pca(adata, n_comps=n_pcs)
    print(f"PCA computed: {n_pcs} components")
    return adata


def compute_neighbors(adata: AnnData, n_neighbors: int = 15, n_pcs: int = 50) -> AnnData:
    """Compute neighborhood graph.
    
    Parameters
    ----------
    adata : AnnData
        Annotated data matrix with PCA coordinates.
    n_neighbors : int, optional
        Number of neighbors (default: 15).
    n_pcs : int, optional
        Number of PCA dimensions to use (default: 50).
        
    Returns
    -------
    AnnData
        Annotated data matrix with neighbor graph.
    """
    sc.pp.neighbors(adata, n_neighbors=n_neighbors, n_pcs=n_pcs)
    print(f"Neighbors computed: k={n_neighbors}")
    return adata


def compute_umap(adata: AnnData, min_dist: float = 0.1) -> AnnData:
    """Compute UMAP embedding.
    
    Parameters
    ----------
    adata : AnnData
        Annotated data matrix with neighbor graph.
    min_dist : float, optional
        Minimum distance for UMAP (default: 0.1).
        
    Returns
    -------
    AnnData
        Annotated data matrix with UMAP coordinates in obsm['X_umap'].
    """
    sc.tl.umap(adata, min_dist=min_dist, random_state=42)
    print("UMAP computed")
    return adata


def leiden_clustering(adata: AnnData, resolution: float = 1.0) -> AnnData:
    """Leiden clustering.
    
    Parameters
    ----------
    adata : AnnData
        Annotated data matrix with neighbor graph.
    resolution : float, optional
        Resolution parameter (default: 1.0). Higher values lead to more clusters.
        
    Returns
    -------
    AnnData
        Annotated data matrix with cluster assignments in obs['leiden'].
    """
    sc.tl.leiden(adata, resolution=resolution, random_state=42)
    n_clusters = len(adata.obs['leiden'].unique())
    print(f"Leiden clustering: {n_clusters} clusters found")
    return adata


def clustering_pipeline(
    adata: AnnData,
    n_pcs: int = 50,
    n_neighbors: int = 15,
    min_dist: float = 0.1,
    resolution: float = 1.0
) -> AnnData:
    """Complete clustering pipeline.
    
    Parameters
    ----------
    adata : AnnData
        Annotated data matrix (preprocessed).
    n_pcs : int, optional
        Number of PCA components (default: 50).
    n_neighbors : int, optional
        Number of neighbors (default: 15).
    min_dist : float, optional
        UMAP minimum distance (default: 0.1).
    resolution : float, optional
        Leiden resolution (default: 1.0).
        
    Returns
    -------
    AnnData
        Annotated data matrix with PCA, UMAP, and cluster assignments.
    """
    print("Starting clustering pipeline...")
    adata = compute_pca(adata, n_pcs)
    adata = compute_neighbors(adata, n_neighbors, n_pcs)
    adata = compute_umap(adata, min_dist)
    adata = leiden_clustering(adata, resolution)
    print("Clustering complete!")
    return adata

