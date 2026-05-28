"""
Unit tests for marker gene detection module.
"""

import pytest
import numpy as np
import pandas as pd
from anndata import AnnData
import scanpy as sc
import tempfile
import os
from src.markers import (
    rank_marker_genes,
    get_marker_genes,
    export_marker_genes,
    marker_gene_pipeline
)


@pytest.fixture
def clustered_adata():
    """Create clustered AnnData for marker tests."""
    np.random.seed(42)
    n_cells, n_genes = 100, 200
    
    X = np.random.negative_binomial(5, 0.3, size=(n_cells, n_genes)).astype(np.float32)
    
    adata = AnnData(X=X)
    adata.var_names = [f"GENE_{i}" for i in range(n_genes)]
    adata.obs['cell_id'] = [f"cell_{i}" for i in range(n_cells)]
    
    # Preprocess
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    
    # Cluster
    sc.tl.pca(adata, n_comps=20)
    sc.pp.neighbors(adata, n_neighbors=15, n_pcs=20)
    sc.tl.leiden(adata, resolution=0.5)
    
    return adata


def test_rank_marker_genes(clustered_adata):
    """Test marker gene ranking."""
    adata = clustered_adata.copy()
    
    adata_ranked = rank_marker_genes(adata, groupby='leiden', method='wilcoxon')
    
    assert 'rank_genes_groups' in adata_ranked.uns
    ranking = adata_ranked.uns['rank_genes_groups']
    assert 'names' in ranking
    assert 'scores' in ranking
    assert 'logfoldchanges' in ranking
    assert 'pvals_adj' in ranking


def test_get_marker_genes(clustered_adata):
    """Test marker gene extraction."""
    adata = clustered_adata.copy()
    adata = rank_marker_genes(adata, groupby='leiden')
    
    markers = get_marker_genes(adata, n_genes=10)
    
    assert isinstance(markers, dict)
    assert len(markers) > 0
    
    # Check each cluster has a dataframe
    for cluster, df in markers.items():
        assert isinstance(df, pd.DataFrame)
        assert len(df) == 10
        assert 'gene' in df.columns
        assert 'score' in df.columns


def test_export_marker_genes(clustered_adata):
    """Test marker gene export."""
    adata = clustered_adata.copy()
    adata = rank_marker_genes(adata)
    markers = get_marker_genes(adata, n_genes=5)
    
    with tempfile.TemporaryDirectory() as tmpdir:
        export_marker_genes(markers, tmpdir)
        
        # Check files exist
        files = os.listdir(tmpdir)
        assert len(files) > 0
        assert any(f.startswith('markers_cluster_') for f in files)


def test_marker_gene_pipeline(clustered_adata):
    """Test complete marker gene pipeline."""
    adata = clustered_adata.copy()
    
    adata_result, markers = marker_gene_pipeline(
        adata,
        groupby='leiden',
        method='wilcoxon',
        n_genes=10
    )
    
    assert 'rank_genes_groups' in adata_result.uns
    assert isinstance(markers, dict)
    assert len(markers) > 0


def test_marker_pipeline_with_export(clustered_adata):
    """Test marker pipeline with export."""
    adata = clustered_adata.copy()
    
    with tempfile.TemporaryDirectory() as tmpdir:
        adata_result, markers = marker_gene_pipeline(
            adata,
            n_genes=10,
            output_dir=tmpdir
        )
        
        # Check files
        files = os.listdir(tmpdir)
        assert len(files) > 0


def test_marker_genes_content(clustered_adata):
    """Test that marker genes have expected columns."""
    adata = clustered_adata.copy()
    adata = rank_marker_genes(adata)
    markers = get_marker_genes(adata, n_genes=5)
    
    for cluster, df in markers.items():
        assert 'gene' in df.columns
        assert 'score' in df.columns
        assert 'logfoldchange' in df.columns
        assert 'pval_adj' in df.columns
        
        # Check values are numeric
        assert df['score'].dtype in [np.float32, np.float64]
        assert df['logfoldchange'].dtype in [np.float32, np.float64]


def test_invalid_method_raises_error(clustered_adata):
    """Test that invalid method raises error."""
    adata = clustered_adata.copy()
    
    with pytest.raises(ValueError):
        rank_marker_genes(adata, method='invalid_method')


def test_missing_key_error(clustered_adata):
    """Test error when ranking key not found."""
    adata = clustered_adata.copy()
    
    with pytest.raises(ValueError):
        get_marker_genes(adata, key='nonexistent_key')


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
