# Quick Start Guide

## Installation

```bash
# Clone or navigate to the project
cd singlecell-discovery

# Create a virtual environment (recommended)
conda create -n scrna python=3.10
conda activate scrna

# Install dependencies
pip install -r requirements.txt
```

## Basic Usage

### Option 1: One-Command Pipeline (Easiest)

```bash
python run_pipeline.py \
    --input data/sample.h5ad \
    --output results/ \
    --resolution 1.0
```

All parameters with defaults:
```
--input              Input H5AD file (required)
--output             Output directory (default: results/)
--min-genes          Min genes per cell (default: 200)
--max-genes          Max genes per cell (default: 2500)
--max-mito           Max mitochondrial % (default: 0.05)
--n-genes            Highly variable genes (default: 2000)
--n-pcs              PCA components (default: 50)
--resolution         Leiden resolution (default: 1.0)
--n-markers          Top markers per cluster (default: 10)
```

### Option 2: Python API (Most Flexible)

```python
from src import (
    load_h5ad,
    preprocess_pipeline,
    clustering_pipeline,
    marker_gene_pipeline
)
from src.visualization import plot_umap_interactive

# Load data
adata = load_h5ad("data/sample.h5ad")

# Preprocessing
adata = preprocess_pipeline(adata)

# Clustering
adata = clustering_pipeline(adata)

# Find markers
adata, markers = marker_gene_pipeline(adata, n_genes=10)

# Visualize
fig = plot_umap_interactive(adata, color_by='leiden')
fig.write_html("results/umap.html")

# Save results
adata.write("results/processed.h5ad")
```

### Option 3: Snakemake Workflow (Most Professional)

```bash
# Check configuration
cat config.yaml

# Dry run to see what will execute
snakemake -n

# Run workflow with 4 cores
snakemake --cores 4

# Run on SLURM cluster
snakemake --profile slurm --cores 40
```

---

## Output Structure

```
results/
├── processed_data.h5ad      # Final processed data
├── preprocessed.h5ad        # After filtering/normalization
├── clustered.h5ad           # With clustering results
├── plots/
│   ├── umap_clusters.html   # Interactive UMAP plot
│   └── ...
└── markers/
    ├── markers_cluster_0.csv
    ├── markers_cluster_1.csv
    └── ...
```

---

## Common Workflows

### Explore Clusters
```python
from src.visualization import plot_umap_interactive

# Interactive UMAP
fig = plot_umap_interactive(adata, color_by='leiden')
fig.show()
```

### View Marker Genes
```python
markers_df = markers['0']  # For cluster 0
print(markers_df.head(10))
```

### Plot Gene Expression
```python
from src.visualization import violin_plot

# Violin plot
violin_plot(adata, genes=['CD4', 'CD8'], groupby='leiden')
```

### Heatmap of Markers
```python
from src.visualization import plot_heatmap

genes = markers['0']['gene'].head(5).tolist()
plot_heatmap(adata, genes, groupby='leiden')
```

---

## Troubleshooting

### Memory Issues
- Reduce `--n-genes` or `--n-pcs`
- Filter data first with stricter `--min-genes`

### Slow Computation
- Use fewer PCs: `--n-pcs 30`
- Use Snakemake parallel: `snakemake --cores 8`

### Missing Gene Annotations
- Ensure gene names are in `adata.var_names`
- Check for case sensitivity (CD4 vs cd4)

### Plot Not Showing
- Use Jupyter notebook instead of terminal
- Or save to HTML: `fig.write_html("plot.html")`

---

## Data Format

Input must be HDF5-AnnData format (.h5ad):
```python
# If you have a CSV file
import scanpy as sc
import pandas as pd

df = pd.read_csv("my_data.csv", index_col=0)
adata = sc.AnnData(df.T)  # Transpose: genes x cells
adata.write("data.h5ad")
```

---

## Example Analysis

From raw data to results in one script:

```python
import scanpy as sc
from src import *
from src.visualization import *

# Load and process
adata = load_h5ad("sample.h5ad")
adata = preprocess_pipeline(adata, n_top_genes=3000)
adata = clustering_pipeline(adata, resolution=0.8)
adata, markers = marker_gene_pipeline(adata, n_genes=20)

# Visualize
plot_umap_interactive(adata, 'leiden').write_html("umap.html")
plot_heatmap_interactive(
    adata,
    markers['0']['gene'].head(10).tolist(),
    groupby='leiden'
).write_html("heatmap.html")

# Save
adata.write("results.h5ad")
```

---

## Next Steps

1. **Prepare your data** in H5AD format
2. **Run the pipeline** with appropriate parameters
3. **Explore results** in `results/` directory
4. **Fine-tune parameters** if needed
5. **Generate publication figures** using saved plots

---

## Support

- See `IMPROVEMENTS.md` for detailed changes
- Check docstrings in source code: `help(preprocessing.preprocess_pipeline)`
- Scanpy documentation: https://scanpy.readthedocs.io/
- AnnData documentation: https://anndata.readthedocs.io/
