"""Heatmap visualization for marker gene expression."""

from typing import List, Optional
import matplotlib.pyplot as plt
import plotly.graph_objects as go
import plotly.express as px
import numpy as np
import pandas as pd
import seaborn as sns
import scanpy as sc
from anndata import AnnData

def plot_heatmap(
    adata: AnnData,
    marker_genes: List[str],
    groupby: str,
    figsize: tuple = (12, 8),
    cmap: str = 'viridis',
    save_path: Optional[str] = None
) -> None:
    """Plot a heatmap of marker gene expression using matplotlib.
    
    Parameters
    ----------
    adata : AnnData
        The annotated data matrix.
    marker_genes : list of str
        List of marker genes to plot.
    groupby : str
        The key in adata.obs to group by (e.g., 'leiden', 'cell_type').
    figsize : tuple, optional
        Size of the figure (default: (12, 8)).
    cmap : str, optional
        Colormap to use (default: 'viridis').
    save_path : str, optional
        Path to save the plot. If None, plot is displayed.
    """
    try:
        # Use scanpy's plotting function
        sc.pl.heatmap(
            adata,
            marker_genes,
            groupby=groupby,
            cmap=cmap,
            figsize=figsize,
            show=False,
            save=save_path
        )
        if save_path:
            print(f"Heatmap saved: {save_path}")
        else:
            plt.show()
    except Exception as e:
        print(f"Error creating heatmap: {e}")


def plot_heatmap_interactive(
    adata: AnnData,
    marker_genes: List[str],
    groupby: str,
    save_path: Optional[str] = None
) -> Optional[go.Figure]:
    """Plot interactive heatmap using Plotly.
    
    Parameters
    ----------
    adata : AnnData
        The annotated data matrix.
    marker_genes : list of str
        List of marker genes to plot.
    groupby : str
        The key in adata.obs to group by.
    save_path : str, optional
        Path to save as HTML. If None, figure is returned.
        
    Returns
    -------
    go.Figure or None
        Plotly figure object.
    """
    try:
        # Extract and prepare data
        subset_adata = adata[:, marker_genes].copy()
        groups = adata.obs[groupby]
        
        # Create a DataFrame
        data = subset_adata.X.toarray() if hasattr(subset_adata.X, 'toarray') else subset_adata.X
        df = pd.DataFrame(data, columns=marker_genes)
        df[groupby] = groups.values
        
        # Compute mean expression per group
        mean_expr = df.groupby(groupby).mean()
        
        # Create heatmap
        fig = go.Figure(data=go.Heatmap(
            z=mean_expr.values,
            x=mean_expr.columns,
            y=mean_expr.index,
            colorscale='Viridis'
        ))
        
        fig.update_layout(
            title=f"Marker Gene Expression by {groupby}",
            xaxis_title="Genes",
            yaxis_title=groupby,
            height=600,
            width=1200
        )
        
        if save_path:
            fig.write_html(save_path)
            print(f"Interactive heatmap saved: {save_path}")
        
        return fig
    except Exception as e:
        print(f"Error creating interactive heatmap: {e}")
        return None


