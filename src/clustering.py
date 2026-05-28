"""
PCA for dimensionality reduction

UMAP embedding

Leiden/Louvain clustering with multiple resolution options

Return cluster assignments and embeddings
"""

import numpy as np
import umap
import scanpy as sc
import leidenalg
import igraph as ig
from sklearn.decomposition import PCA

def perform_clustering(data, n_pcs=50, n_neighbors=15, min_dist=0.1, resolution=1.0):
    # Step 1: PCA for dimensionality reduction
    pca = PCA(n_components=n_pcs)
    pca_result = pca.fit_transform(data)

    # Step 2: UMAP embedding
    umap_embedder = umap.UMAP(n_neighbors=n_neighbors, min_dist=min_dist)
    umap_embedding = umap_embedder.fit_transform(pca_result)

    # Step 3: Leiden clustering
    # Create a k-nearest neighbor graph
    knn_graph = sc.pp.neighbors(pca_result, n_neighbors=n_neighbors, return_graph=True)
    
    # Convert to igraph format
    g = ig.Graph.Adjacency((knn_graph > 0).tolist())
    
    # Perform Leiden clustering
    partition = leidenalg.find_partition(g, leidenalg.RBConfigurationVertexPartition, resolution_parameter=resolution)
    
    # Get cluster assignments
    cluster_assignments = np.array(partition.membership)

    return cluster_assignments, umap_embedding

# Example Usage

# data = pd.read_csv('your_data.csv')  # Load your data here
# cluster_assignments, umap_embedding = perform_clustering(data.values)

