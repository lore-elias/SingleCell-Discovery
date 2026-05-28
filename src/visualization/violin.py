"""Violin and dot plots for gene expression visualization."""

from typing import Optional, List
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
import scanpy as sc
from anndata import AnnData

def violin_plot(
    adata: AnnData,
    genes: List[str],
    groupby: Optional[str] = None,
    figsize: tuple = (12, 6),
    save_path: Optional[str] = None
) -> None:
    """Violin plot for gene expression.
    
    Parameters
    ----------
    adata : AnnData
        Annotated data matrix.
    genes : list of str
        Gene names to plot.
    groupby : str, optional
        Column in adata.obs to group by. If None, all cells are plotted.
    figsize : tuple, optional
        Figure size (default: (12, 6)).
    save_path : str, optional
        Path to save the plot. If None, plot is displayed.
    """
    try:
        sc.pl.violin(
            adata,
            genes,
            groupby=groupby,
            show=False,
            figsize=figsize
        )
        if save_path:
            plt.savefig(save_path, dpi=300, bbox_inches='tight')
            print(f"Violin plot saved: {save_path}")
            plt.close()
        else:
            plt.show()
    except Exception as e:
        print(f"Error creating violin plot: {e}")


def dot_plot(
    adata: AnnData,
    genes: List[str],
    groupby: Optional[str] = None,
    figsize: tuple = (12, 6),
    save_path: Optional[str] = None
) -> None:
    """Dot plot for gene expression.
    
    Parameters
    ----------
    adata : AnnData
        Annotated data matrix.
    genes : list of str
        Gene names to plot.
    groupby : str, optional
        Column in adata.obs to group by. If None, all cells plotted together.
    figsize : tuple, optional
        Figure size (default: (12, 6)).
    save_path : str, optional
        Path to save the plot. If None, plot is displayed.
    """
    try:
        sc.pl.dotplot(
            adata,
            genes,
            groupby=groupby,
            show=False,
            figsize=figsize
        )
        if save_path:
            plt.savefig(save_path, dpi=300, bbox_inches='tight')
            print(f"Dot plot saved: {save_path}")
            plt.close()
        else:
            plt.show()
    except Exception as e:
        print(f"Error creating dot plot: {e}")

