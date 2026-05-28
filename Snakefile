# Snakemake workflow for single-cell RNA-seq discovery pipeline
# Usage:
#   snakemake -n                    # Dry run
#   snakemake --cores 4             # Run with 4 cores
#   snakemake --profile slurm       # Run on SLURM cluster

import os

configfile: "config.yaml"

# Configuration
DATA_DIR = config.get("data_dir", "data/")
OUTPUT_DIR = config.get("output_dir", "results/")
INPUT_FILE = config.get("input_file", "sample.h5ad")

# Preprocessing parameters
MIN_GENES = config.get("min_genes", 200)
MAX_GENES = config.get("max_genes", 2500)
MAX_MITO = config.get("max_mito", 0.05)
N_TOP_GENES = config.get("n_top_genes", 2000)

# Clustering parameters
N_PCS = config.get("n_pcs", 50)
N_NEIGHBORS = config.get("n_neighbors", 15)
RESOLUTION = config.get("resolution", 1.0)

# Marker genes
N_MARKERS = config.get("n_markers", 10)

# Create output directories
os.makedirs(OUTPUT_DIR, exist_ok=True)
os.makedirs(os.path.join(OUTPUT_DIR, "plots"), exist_ok=True)
os.makedirs(os.path.join(OUTPUT_DIR, "markers"), exist_ok=True)

rule all:
    input:
        h5ad = f"{OUTPUT_DIR}/processed_data.h5ad",
        markers = f"{OUTPUT_DIR}/markers/markers_complete.txt",
        plots = f"{OUTPUT_DIR}/plots/umap_clusters.html"

rule preprocess:
    input:
        f"{DATA_DIR}/{INPUT_FILE}"
    output:
        f"{OUTPUT_DIR}/preprocessed.h5ad"
    params:
        min_genes = MIN_GENES,
        max_genes = MAX_GENES,
        max_mito = MAX_MITO,
        n_top_genes = N_TOP_GENES
    shell:
        """python -c "
from src.preprocessing import preprocess_pipeline
from src import load_h5ad

adata = load_h5ad('{input}')
adata = preprocess_pipeline(adata, min_genes={params.min_genes}, max_genes={params.max_genes}, max_mito={params.max_mito}, n_top_genes={params.n_top_genes})
adata.write('{output}')
"
        """

rule clustering:
    input:
        f"{OUTPUT_DIR}/preprocessed.h5ad"
    output:
        f"{OUTPUT_DIR}/clustered.h5ad"
    params:
        n_pcs = N_PCS,
        n_neighbors = N_NEIGHBORS,
        resolution = RESOLUTION
    shell:
        """python -c "
from src import clustering_pipeline
import scanpy as sc

adata = sc.read_h5ad('{input}')
adata = clustering_pipeline(adata, n_pcs={params.n_pcs}, n_neighbors={params.n_neighbors}, resolution={params.resolution})
adata.write('{output}')
"
        """

rule markers:
    input:
        f"{OUTPUT_DIR}/clustered.h5ad"
    output:
        h5ad = f"{OUTPUT_DIR}/processed_data.h5ad",
        summary = f"{OUTPUT_DIR}/markers/markers_complete.txt"
    params:
        n_genes = N_MARKERS,
        output_dir = f"{OUTPUT_DIR}/markers"
    shell:
        """python -c "
from src import marker_gene_pipeline
import scanpy as sc

adata = sc.read_h5ad('{input}')
adata, markers = marker_gene_pipeline(adata, n_genes={params.n_genes}, output_dir='{params.output_dir}')
adata.write('{output.h5ad}')

with open('{output.summary}', 'w') as f:
    f.write(f'Marker genes detected for {{len(markers)}} clusters\n')
    for cluster, df in markers.items():
        f.write(f'Cluster {{cluster}}: {{len(df)}} genes\n')
"
        """

rule visualize:
    input:
        f"{OUTPUT_DIR}/processed_data.h5ad"
    output:
        f"{OUTPUT_DIR}/plots/umap_clusters.html"
    shell:
        """python -c "
from src.visualization import plot_umap_interactive
import scanpy as sc

adata = sc.read_h5ad('{input}')
fig = plot_umap_interactive(adata, color_by='leiden', title='UMAP: Leiden Clusters')
if fig:
    fig.write_html('{output}')
"
        """
