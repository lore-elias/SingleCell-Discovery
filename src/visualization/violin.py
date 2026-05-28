"""Violin/dot plots for gene expression"""

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

def violin_plot(adata, gene, groupby=None, log=False, figsize=(8, 6), **kwargs):
    """Violin plot for gene expression

    Parameters
    ----------
    adata : AnnData
        Annotated data matrix.
    gene : str
        Gene name to plot.
    groupby : str, optional
        Column name in adata.obs to group by. If None, all cells are plotted together.
    log : bool, optional
        Whether to log-transform the expression values. Default is False.
    figsize : tuple, optional
        Figure size. Default is (8, 6).
    **kwargs
        Additional keyword arguments passed to sns.violinplot.

    Returns
    -------
    matplotlib.axes.Axes
        The axes object with the violin plot.
    """
    if groupby is not None:
        data = adata.obs[[groupby]].copy()
        data[gene] = adata[:, gene].X.toarray().flatten() if hasattr(adata[:, gene].X, 'toarray') else adata[:, gene].X.flatten()
    else:
        data = pd.DataFrame({gene: adata[:, gene].X.toarray().flatten() if hasattr(adata[:, gene].X, 'toarray') else adata[:, gene].X.flatten()})
    
    if log:
        data[gene] = np.log1p(data[gene])
    
    plt.figure(figsize=figsize)
    ax = sns.violinplot(x=groupby, y=gene, data=data, **kwargs)
    ax.set_title(f'Violin plot of {gene} expression')
    
    return ax

def dot_plot(adata, gene, groupby=None, log=False, figsize=(8, 6), **kwargs):
    """Dot plot for gene expression

    Parameters
    ----------
    adata : AnnData
        Annotated data matrix.
    gene : str
        Gene name to plot.
    groupby : str, optional
        Column name in adata.obs to group by. If None, all cells are plotted together.
    log : bool, optional
        Whether to log-transform the expression values. Default is False.
    figsize : tuple, optional
        Figure size. Default is (8, 6).
    **kwargs
        Additional keyword arguments passed to sns.stripplot.

    Returns
    -------
    matplotlib.axes.Axes
        The axes object with the dot plot.
    """
    if groupby is not None:
        data = adata.obs[[groupby]].copy()
        data[gene] = adata[:, gene].X.toarray().flatten() if hasattr(adata[:, gene].X, 'toarray') else adata[:, gene].X.flatten()
    else:
        data = pd.DataFrame({gene: adata[:, gene].X.toarray().flatten() if hasattr(adata[:, gene].X, 'toarray') else adata[:, gene].X.flatten()})
    
    if log:
        data[gene] = np.log1p(data[gene])
    
    plt.figure(figsize=figsize)
    ax = sns.stripplot(x=groupby, y=gene, data=data, jitter=True, **kwargs)
    ax.set_title(f'Dot plot of {gene} expression')
    
    return ax

