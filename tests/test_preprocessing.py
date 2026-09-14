"""
Unit tests for preprocessing module.
"""

import pytest
import numpy as np
import pandas as pd
from anndata import AnnData
from src.preprocessing import (
    load_h5ad,
    filter_cells,
    normalize_and_log_transform,
    select_variable_genes,
    preprocess_pipeline
)


@pytest.fixture
def sample_adata():
    """Create a sample AnnData object for testing."""
    np.random.seed(42)
    n_cells, n_genes = 100, 200
    
    X = np.random.negative_binomial(5, 0.3, size=(n_cells, n_genes)).astype(np.float32)
    
    adata = AnnData(X=X)
    adata.obs['cell_id'] = [f"cell_{i}" for i in range(n_cells)]
    
    var_names = [f"GENE_{i}" for i in range(n_genes)]

    for i in range(10):
        var_names[i] = f"MT-{i}"

    adata.var_names = var_names
    
    return adata


def test_filter_cells(sample_adata):
    """Test cell filtering."""
    adata = sample_adata.copy()
    initial_cells = adata.shape[0]
    
    adata_filtered = filter_cells(adata, min_genes=50, max_genes=1000, max_mito=0.5)
    
    assert adata_filtered.shape[0] <= initial_cells
    assert 'n_genes' in adata_filtered.obs.columns
    assert 'n_counts' in adata_filtered.obs.columns
    assert 'mito_percent' in adata_filtered.obs.columns


def test_normalize_and_log_transform(sample_adata):
    """Test normalization and log transformation."""
    adata = sample_adata.copy()
    
    adata_norm = normalize_and_log_transform(adata, target_sum=1e4)
    
    assert adata_norm.shape == sample_adata.shape
    # Check that data has been transformed (values should be different)
    assert not np.allclose(adata_norm.X, sample_adata.X)


def test_select_variable_genes(sample_adata):
    """Test highly variable gene selection."""
    adata = sample_adata.copy()
    initial_genes = adata.shape[1]
    
    # First normalize
    adata = normalize_and_log_transform(adata)
    
    # Select variable genes
    adata_hvg = select_variable_genes(adata, n_top_genes=50)
    
    assert adata_hvg.shape[1] == 50
    assert adata_hvg.shape[1] < initial_genes
    assert 'highly_variable' in adata_hvg.var.columns


def test_preprocess_pipeline(sample_adata):
    """Test complete preprocessing pipeline."""
    adata = sample_adata.copy()
    
    adata_processed = preprocess_pipeline(
        adata,
        min_genes=30,
        max_genes=500,
        max_mito=0.5,
        n_top_genes=50
    )
    
    # Check pipeline effects
    assert adata_processed.shape[1] == 50  # Should have 50 variable genes
    assert adata_processed.shape[0] <= adata.shape[0]  # Cells filtered
    assert 'highly_variable' in adata_processed.var.columns


def test_filter_cells_edge_cases(sample_adata):
    """Test edge cases in cell filtering."""
    adata = sample_adata.copy()
    
    # Very strict filtering
    adata_filtered = filter_cells(adata, min_genes=10000, max_genes=20000)
    assert adata_filtered.shape[0] == 0  # Should remove all cells
    
    # Very lenient filtering
    adata_filtered = filter_cells(adata, min_genes=1, max_genes=100000, max_mito=1.0)
    assert adata_filtered.shape[0] > 0  # Should keep all cells


def test_normalize_preserves_shape(sample_adata):
    """Test that normalization preserves data shape."""
    adata = sample_adata.copy()
    original_shape = adata.shape
    
    adata_norm = normalize_and_log_transform(adata)
    
    assert adata_norm.shape == original_shape


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
