"""
Marker gene detection module for single-cell RNA-seq data.

Functions for differential expression analysis and marker gene identification.
"""

from typing import Dict, Optional, Tuple
import numpy as np
import pandas as pd
import scanpy as sc
from anndata import AnnData

def rank_marker_genes(
    adata: AnnData,
    groupby: str = 'leiden',
    method: str = 'wilcoxon',
    key_added: str = 'rank_genes_groups'
) -> AnnData:
    """Identify marker genes per cluster using ranking.
    
    Parameters
    ----------
    adata : AnnData
        Annotated data matrix with cluster assignments.
    groupby : str, optional
        Column in obs containing cluster labels (default: 'leiden').
    method : str, optional
        Statistical test ('wilcoxon' or 't-test', default: 'wilcoxon').
    key_added : str, optional
        Key to store results (default: 'rank_genes_groups').
        
    Returns
    -------
    AnnData
        Annotated data matrix with marker genes in .uns[key_added].
    """
    valid_methods = ['wilcoxon', 't-test']
    if method not in valid_methods:
        raise ValueError(f"Method must be one of {valid_methods}, got {method}")
    
    sc.tl.rank_genes_groups(adata, groupby=groupby, method=method, key_added=key_added)
    print(f"Marker genes ranked ({method})")
    return adata

def get_marker_genes(
    adata: AnnData,
    n_genes: int = 10,
    key: str = 'rank_genes_groups'
) -> Dict[str, pd.DataFrame]:
    """Extract top marker genes per cluster.
    
    Parameters
    ----------
    adata : AnnData
        Annotated data matrix with ranked marker genes.
    n_genes : int, optional
        Number of top genes per cluster (default: 10).
    key : str, optional
        Key containing ranking results (default: 'rank_genes_groups').
        
    Returns
    -------
    dict
        Dictionary mapping cluster names to DataFrames of top marker genes.
    """
    if key not in adata.uns:
        raise ValueError(f"Key '{key}' not found in adata.uns. Run rank_marker_genes first.")
    
    marker_dict = {}
    ranking = adata.uns[key]
    clusters = ranking['names'].dtype.names
    
    for cluster in clusters:
        genes = ranking['names'][cluster][:n_genes]
        scores = ranking['scores'][cluster][:n_genes]
        logfolds = np.full(len(genes), np.nan, dtype=float)
        if 'logfoldchanges' in ranking:
            logfolds = ranking['logfoldchanges'][cluster][:n_genes]
        if 'pvals_adj' in ranking:
            pvals = ranking['pvals_adj'][cluster][:n_genes]
        else:
            pvals = np.full(len(genes), np.nan, dtype=float)
        
        marker_dict[cluster] = pd.DataFrame({
            'gene': genes,
            'score': scores,
            'logfoldchange': logfolds,
            'pval_adj': pvals
        })
    
    return marker_dict

def export_marker_genes(
    marker_dict: Dict[str, pd.DataFrame],
    output_dir: str = 'results/'
) -> None:
    """Export marker genes to CSV files.
    
    Parameters
    ----------
    marker_dict : dict
        Dictionary of marker gene DataFrames per cluster.
    output_dir : str, optional
        Output directory for CSV files (default: 'results/').
    """
    import os
    os.makedirs(output_dir, exist_ok=True)
    
    for cluster, df in marker_dict.items():
        filepath = os.path.join(output_dir, f'markers_cluster_{cluster}.csv')
        df.to_csv(filepath, index=False)
        print(f"Exported: {filepath}")


def marker_gene_pipeline(
    adata: AnnData,
    groupby: str = 'leiden',
    method: str = 'wilcoxon',
    n_genes: int = 10,
    output_dir: Optional[str] = None
) -> Tuple[AnnData, Dict[str, pd.DataFrame]]:
    """Complete marker gene detection pipeline.
    
    Parameters
    ----------
    adata : AnnData
        Annotated data matrix with cluster assignments.
    groupby : str, optional
        Column containing cluster labels (default: 'leiden').
    method : str, optional
        Statistical test method (default: 'wilcoxon').
    n_genes : int, optional
        Number of top genes per cluster (default: 10).
    output_dir : str, optional
        Directory to export marker genes. If None, no export.
        
    Returns
    -------
    tuple
        Updated AnnData object and dictionary of marker gene DataFrames.
    """
    print("Starting marker gene detection pipeline...")
    adata = rank_marker_genes(adata, groupby=groupby, method=method)
    marker_dict = get_marker_genes(adata, n_genes=n_genes)
    
    if output_dir:
        export_marker_genes(marker_dict, output_dir)
    
    print(f"Marker genes detected for {len(marker_dict)} clusters")
    return adata, marker_dict

