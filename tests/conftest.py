"""
Pytest configuration and shared fixtures.
"""

import pytest
import numpy as np
from anndata import AnnData


@pytest.fixture(scope="session")
def random_seed():
    """Set random seed for reproducibility."""
    np.random.seed(42)
    return 42


@pytest.fixture
def basic_adata():
    """Create a basic AnnData object for testing."""
    np.random.seed(42)
    n_cells, n_genes = 50, 100
    
    X = np.random.poisson(5, size=(n_cells, n_genes)).astype(np.float32)
    
    adata = AnnData(X=X)
    adata.var_names = [f"GENE_{i:03d}" for i in range(n_genes)]
    adata.var.index.name = "gene_names"
    adata.obs['cell_id'] = [f"cell_{i:04d}" for i in range(n_cells)]
    adata.obs.index.name = "cell_names"
    
    return adata


@pytest.fixture
def sample_adata():
    """Create a sample AnnData object for annotation and cluster tests."""
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
