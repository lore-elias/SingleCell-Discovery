"""UMAP plots colored by cluster/metadata"""

import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns

def plot_umap(adata, color_by, title = None, save_path = None):
    """Plot UMAP colored by cluster/metadata

    Parameters
    ----------
    adata : AnnData
        Annotated data matrix.
    color_by : str
        Column name in adata.obs to color by.
    title : str, optional
        Title of the plot, by default None.
    save_path : str, optional
        Path to save the plot, by default None.

    Returns
    -------
    None
    """
    # Get UMAP coordinates and metadata
    umap_coords = adata.obsm['X_umap']
    metadata = adata.obs[color_by]

    # Create a DataFrame for plotting
    plot_data = np.hstack((umap_coords, metadata.values[:, None]))
    
    # Create a scatter plot
    plt.figure(figsize=(10, 8))
    sns.scatterplot(x=plot_data[:, 0], y=plot_data[:, 1], hue=plot_data[:, 2], palette='tab10', s=50)
    
    # Set title and labels
    if title:
        plt.title(title)
    plt.xlabel('UMAP 1')
    plt.ylabel('UMAP 2')
    
    # Save or show the plot
    if save_path:
        plt.savefig(save_path)
        plt.close()
    else:
        plt.show(
        )
