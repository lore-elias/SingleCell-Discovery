"""UMAP visualization plots for single-cell RNA-seq data."""

from typing import Optional
import matplotlib.pyplot as plt
import plotly.graph_objects as go
import plotly.express as px
import scanpy as sc
from anndata import AnnData

def plot_umap(
    adata: AnnData,
    color_by: str,
    title: Optional[str] = None,
    save_path: Optional[str] = None,
    figsize: tuple = (10, 8)
) -> None:
    """Plot UMAP colored by cluster/metadata using matplotlib.
    
    Parameters
    ----------
    adata : AnnData
        Annotated data matrix with UMAP coordinates.
    color_by : str
        Column name in adata.obs to color by.
    title : str, optional
        Title of the plot. If None, no title is shown.
    save_path : str, optional
        Path to save the plot. If None, plot is displayed.
    figsize : tuple, optional
        Figure size (default: (10, 8)).
    """
    try:
        sc.pl.umap(adata, color=color_by, title=title, size=50, show=False)
        plt.figure(figsize=figsize)
        if title:
            plt.title(title)
        if save_path:
            plt.savefig(save_path, dpi=300, bbox_inches='tight')
            print(f"Plot saved: {save_path}")
            plt.close()
        else:
            plt.show()
    except Exception as e:
        print(f"Error creating UMAP plot: {e}")


def plot_umap_interactive(
    adata: AnnData,
    color_by: str,
    title: Optional[str] = None,
    save_path: Optional[str] = None
) -> go.Figure:
    """Plot interactive UMAP using Plotly.
    
    Parameters
    ----------
    adata : AnnData
        Annotated data matrix with UMAP coordinates.
    color_by : str
        Column name in adata.obs to color by.
    title : str, optional
        Title of the plot. If None, a default title is used.
    save_path : str, optional
        Path to save the plot as HTML. If None, plot is returned.
        
    Returns
    -------
    go.Figure
        Plotly figure object.
    """
    try:
        umap_coords = adata.obsm['X_umap']
        colors = adata.obs[color_by].astype(str)
        
        fig = px.scatter(
            x=umap_coords[:, 0],
            y=umap_coords[:, 1],
            color=colors,
            title=title or f"UMAP colored by {color_by}",
            labels={'x': 'UMAP 1', 'y': 'UMAP 2'},
            hover_name=adata.obs.index
        )
        
        fig.update_traces(marker=dict(size=8, opacity=0.7))
        
        if save_path:
            fig.write_html(save_path)
            print(f"Interactive plot saved: {save_path}")
        
        return fig
    except Exception as e:
        print(f"Error creating interactive UMAP: {e}")
        return None
