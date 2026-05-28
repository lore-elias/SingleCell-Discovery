"""
Snakemake workflow for single-cell RNA-seq discovery pipeline.

Usage:
    snakemake -n                          # Dry run
    snakemake --cores 4                   # Run with 4 cores
    snakemake --profile slurm             # Run on SLURM cluster
"""

configfile: "config.yaml"

import os
from pathlib import Path

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

# Create output directory
os.makedirs(OUTPUT_DIR, exist_ok=True)


rule all:
    """Final target."""
    input:
        h5ad=f"{OUTPUT_DIR}/processed_data.h5ad",
        markers=f"{OUTPUT_DIR}/markers/markers_complete.txt",
        plots=f"{OUTPUT_DIR}/plots/umap_clusters.html"


rule preprocess:
    """Preprocessing step: QC filtering, normalization, gene selection."""
    input:
        f"{DATA_DIR}/{INPUT_FILE}"
    output:
        f"{OUTPUT_DIR}/preprocessed.h5ad"
    params:
        min_genes=MIN_GENES,
        max_genes=MAX_GENES,
        max_mito=MAX_MITO,
        n_top_genes=N_TOP_GENES
    script:
        "scripts/preprocess.py"


rule clustering:
    """Clustering step: PCA, UMAP, Leiden clustering."""
    input:
        f"{OUTPUT_DIR}/preprocessed.h5ad"
    output:
        f"{OUTPUT_DIR}/clustered.h5ad"
    params:
        n_pcs=N_PCS,
        n_neighbors=N_NEIGHBORS,
        resolution=RESOLUTION
    script:
        "scripts/clustering.py"


rule markers:
    """Marker gene detection."""
    input:
        f"{OUTPUT_DIR}/clustered.h5ad"
    output:
        h5ad=f"{OUTPUT_DIR}/processed_data.h5ad",
        summary=f"{OUTPUT_DIR}/markers/markers_complete.txt"
    params:
        n_genes=N_MARKERS,
        output_dir=f"{OUTPUT_DIR}/markers"
    script:
        "scripts/markers.py"


rule visualize:
    """Generate visualizations."""
    input:
        f"{OUTPUT_DIR}/processed_data.h5ad"
    output:
        f"{OUTPUT_DIR}/plots/umap_clusters.html"
    params:
        output_dir=f"{OUTPUT_DIR}/plots"
    script:
        "scripts/visualize.py"
