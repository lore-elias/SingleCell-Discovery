"""
Download datasets from CellxGene (scanpy/CELLxGENE API)

Load H5AD or count matrices

Filter low-quality cells (QC metrics: n_genes, n_counts, mitochondrial %) Normalize & log-transform counts

Variable gene selection
"""

import scanpy as sc
import numpy as np

def load_h5ad(file_path):
    """Load an H5AD file into an AnnData object"""
    adata = sc.read_h5ad(file_path)
    return adata

def filter_cells(adata, min_genes = 200, max_genes = 2500, max_mito = 0.05):
    """Filter cells based on QC metrics"""
    adata.obs['n_genes'] = (adata.X > 0).sum(axis=1)
    adata.obs['n_counts'] = adata.X.sum(axis=1)
    adata.obs['mito_percent'] = np.sum(adata[:, adata.var_names.str.startswith('MT-')].X, axis=1) / adata.obs['n_counts']
    
    sc.pp.filter_cells(adata, min_genes=min_genes)
    sc.pp.filter_cells(adata, max_genes=max_genes)
    adata = adata[adata.obs['mito_percent'] < max_mito]
    
    return adata

def normalize_and_log_transform(adata):
    """Normalize counts and log-transform"""
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    return adata

def select_variable_genes(adata, n_top_genes=2000):
    """Select highly variable genes"""
    sc.pp.highly_variable_genes(adata, n_top_genes=n_top_genes)
    adata= adata[:, adata.var['highly_variable']]
    return adata


