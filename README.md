# Single-cell RNA-seq Discovery Platform

A portfolio project that demonstrates a complete scRNA-seq analysis workflow: quality control, normalization, dimensionality reduction, clustering, marker discovery and interactive visualization.

This repository is designed as a clean, reproducible prototype for single-cell discovery workflows using AnnData and Scanpy.

## :gem: What the project does

- Loads `.h5ad` data
- Filters low-quality cells
- Normalizes and log-transforms counts
- Selects highly variable genes
- Computes PCA, neighbors, and UMAP
- Runs Leiden clustering
- Identifies cluster marker genes
- Exports results and interactive plots
- Provides a Streamlit dashboard for inspection

## :house: Project structure

```text
singlecell-discovery/
├── app/
│   └── streamlit_app.py
├── data/
│   ├── generate_sample_data.py
│   └── sample.h5ad
├── src/
│   ├── annotation.py
│   ├── clustering.py
│   ├── markers.py
│   ├── preprocessing.py
│   └── visualization/
├── tests/
├── config.yaml
├── QUICKSTART.md
├── Snakefile
├── run_pipeline.py
├── requirements.txt
└── README.md
```

## :gear: Tech stack

- Python
- Scanpy
- AnnData
- UMAP
- Plotly
- Streamlit
- Snakemake

## :facepunch: Quick start

### 1) Create the environment

```bash
cd singlecell-discovery
conda create -n scrna python=3.9
conda activate scrna
pip install -r requirements.txt
```

### 2) Generate sample data

```bash
python data/generate_sample_data.py
```

This creates a demo `.h5ad` file in `data/sample.h5ad`.

### 3) Run the pipeline

```bash
python run_pipeline.py \
  --input data/sample.h5ad \
  --output results/demo \
  --resolution 1.0
```

### 4) Explore outputs

The run creates an output folder such as:

```text
results/demo/
├── processed_data.h5ad
├── plots/
│   └── umap_clusters.html
├── markers/
│   ├── markers_cluster_0.csv
│   └── ...
└── ...
```

## :globe_with_meridians: Dashboard

You can launch the interactive dashboard with:

```bash
streamlit run app/streamlit_app.py
```

The app allows uploading an H5AD file and exploring clustering and marker-gene results interactively.

## :computer: Example workflow

1. Load scRNA-seq data
2. Filter poor-quality cells
3. Normalize counts
4. Select variable genes
5. Cluster the cells
6. Rank marker genes per cluster
7. Visualize with UMAP and expression plots

## :orange_book: Notes

This repository is a prototype and portfolio project for demonstrating scRNA-seq analysis workflows. It is best suited for clean demo datasets and educational exploration rather than fully production-scale biomedical pipelines.

## :rocket: Future improvements

- stronger validation for different input schema variants
- more explicit support for metadata columns
- better cell-type annotation workflows
- publication-ready plotting and report generation
- more robust CLI handling for edge cases

## :balance_scale:  License

MIT