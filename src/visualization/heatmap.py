"""Marker gene heatmaps"""

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

def plot_heatmap(adata, marker_genes, groupby, figsize=(10, 10), cmap='viridis'):
    """Plot a heatmap of marker gene expression.

    Parameters
    ----------
    adata : AnnData
        The annotated data matrix.
    marker_genes : list of str
        List of marker genes to plot.
    groupby : str
        The key in adata.obs to group by (e.g., cell type).
    figsize : tuple, optional
        Size of the figure (default: (10, 10)).
    cmap : str, optional
        Colormap to use for the heatmap (default: 'viridis').

    Returns
    -------
    None
        Displays the heatmap.
    """
    # Extract the relevant data
    data = adata[:, marker_genes].X.toarray() if hasattr(adata[:, marker_genes].X, 'toarray') else adata[:, marker_genes].X
    groups = adata.obs[groupby]

    # Create a DataFrame for plotting
    df = pd.DataFrame(data, columns=marker_genes)
    df[groupby] = groups.values

    # Compute mean expression for each group
    mean_expression = df.groupby(groupby).mean()

    # Plot the heatmap
    plt.figure(figsize=figsize)
    sns.heatmap(mean_expression.T, cmap=cmap)
    plt.title('Marker Gene Expression Heatmap')
    plt.xlabel(groupby)
    plt.ylabel('Marker Genes')
    plt.show()


