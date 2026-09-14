"""
Unit tests for clustering module.
"""

import pytest
import numpy as np
from anndata import AnnData
import scanpy as sc
from src.clustering import (
    compute_pca,
    compute_neighbors,
    compute_umap,
    leiden_clustering,
    clustering_pipeline
)


@pytest.fixture
def preprocessed_adata():
    """Create preprocessed AnnData for clustering tests."""
    np.random.seed(42)
    n_cells, n_genes = 100, 200
    
    X = np.random.negative_binomial(5, 0.3, size=(n_cells, n_genes)).astype(np.float32)
    
    adata = AnnData(X=X)
    adata.var_names = [f"GENE_{i}" for i in range(n_genes)]
    adata.obs['cell_id'] = [f"cell_{i}" for i in range(n_cells)]
    
    # Preprocess
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    
    return adata


def test_compute_pca(preprocessed_adata):
    """Test PCA computation."""
    adata = preprocessed_adata.copy()
    
    adata_pca = compute_pca(adata, n_pcs=30)
    
    assert 'X_pca' in adata_pca.obsm
    assert adata_pca.obsm['X_pca'].shape == (adata.shape[0], 30)


def test_compute_neighbors(preprocessed_adata):
    """Test neighbor graph computation."""
    adata = preprocessed_adata.copy()
    adata = compute_pca(adata, n_pcs=20)
    
    adata_neighbors = compute_neighbors(adata, n_neighbors=15, n_pcs=20)
    
    assert 'distances' in adata_neighbors.obsp
    assert 'connectivities' in adata_neighbors.obsp
    assert 'distances' in adata_neighbors.obsp 
    assert adata_neighbors.obsp['distances'].shape[0] == adata.shape[0]


def test_compute_umap(preprocessed_adata):
    """Test UMAP computation."""
    adata = preprocessed_adata.copy()
    adata = compute_pca(adata, n_pcs=20)
    adata = compute_neighbors(adata, n_neighbors=15, n_pcs=20)
    
    adata_umap = compute_umap(adata, min_dist=0.1)
    
    assert 'X_umap' in adata_umap.obsm
    assert adata_umap.obsm['X_umap'].shape == (adata.shape[0], 2)


def test_leiden_clustering(preprocessed_adata):
    """Test Leiden clustering."""
    adata = preprocessed_adata.copy()
    adata = compute_pca(adata, n_pcs=20)
    adata = compute_neighbors(adata, n_neighbors=15, n_pcs=20)
    
    adata_clustered = leiden_clustering(adata, resolution=0.5)
    
    assert 'leiden' in adata_clustered.obs
    assert adata_clustered.obs['leiden'].notna().all()


def test_clustering_pipeline(preprocessed_adata):
    """Test complete clustering pipeline."""
    adata = preprocessed_adata.copy()
    
    adata_clustered = clustering_pipeline(
        adata,
        n_pcs=20,
        n_neighbors=10,
        min_dist=0.1,
        resolution=0.5
    )
    
    # Check all results are present
    assert 'X_pca' in adata_clustered.obsm
    assert 'X_umap' in adata_clustered.obsm
    assert 'leiden' in adata_clustered.obs.columns
    assert 'distances' in adata_clustered.obsp
    assert 'connectivities' in adata_clustered.obsp
    assert 'neighbors' in adata_clustered.uns


def test_clustering_consistency(preprocessed_adata):
    """Test that clustering is reproducible with same seed."""
    adata1 = preprocessed_adata.copy()
    adata2 = preprocessed_adata.copy()
    
    result1 = clustering_pipeline(adata1, resolution=0.5)
    result2 = clustering_pipeline(adata2, resolution=0.5)
    
    # Results should be identical (UMAP may vary slightly due to randomness)
    assert result1.obs['leiden'].equals(result2.obs['leiden'])


def test_leiden_resolution_effect(preprocessed_adata):
    """Test that higher resolution produces more clusters."""
    adata_low = preprocessed_adata.copy()
    adata_high = preprocessed_adata.copy()
    
    adata_low = clustering_pipeline(adata_low, resolution=0.1)
    adata_high = clustering_pipeline(adata_high, resolution=2.0)
    
    n_clusters_low = len(adata_low.obs['leiden'].unique())
    n_clusters_high = len(adata_high.obs['leiden'].unique())
    
    # Higher resolution should give more clusters
    assert n_clusters_high >= n_clusters_low


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
