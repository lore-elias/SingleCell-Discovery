# Quick start guide

This guide is meant to be the shortest route from a fresh clone to a working demo run.

## 1) Install dependencies

```bash
cd singlecell-discovery
conda create -n scrna python=3.9
conda activate scrna
pip install -r requirements.txt
```

## 2) Generate demo data

```bash
python data/generate_sample_data.py
```

This creates a synthetic AnnData dataset in `data/sample.h5ad`.

## 3) Run the pipeline

```bash
python run_pipeline.py \
  --input data/sample.h5ad \
  --output results/demo \
  --resolution 1.0
```

## 4) Open the results

The generated results include:

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

## 5) Launch the dashboard

```bash
streamlit run app/streamlit_app.py
```

## Common parameter examples

```bash
python run_pipeline.py --input data/sample.h5ad --output results/demo --n-pcs 30 --resolution 0.8
python run_pipeline.py --input data/sample.h5ad --output results/demo --min-genes 100 --max-genes 4000
```

## Python API example

```python
from src import load_h5ad, preprocess_pipeline, clustering_pipeline, marker_gene_pipeline

adata = load_h5ad("data/sample.h5ad")
adata = preprocess_pipeline(adata)
adata = clustering_pipeline(adata, resolution=0.8)
adata, markers = marker_gene_pipeline(adata, n_genes=10)

print(list(markers.keys())[:5])
adata.write("results/demo_processed.h5ad")
```

## Data format expected

The input should be a standard single-cell AnnData file in `.h5ad` format.

## Troubleshooting

- If the pipeline fails, make sure the environment is activated and dependencies are installed.
- If `sample.h5ad` is missing, run `python data/generate_sample_data.py` first.
- If the results seem too sparse, lower the filtering thresholds or use a smaller PCA setting.

## Notes

This project is intended as a clean demonstration workflow for a portfolio, not as a final production-grade biomedicine pipeline.
- Scanpy documentation: https://scanpy.readthedocs.io/
- AnnData documentation: https://anndata.readthedocs.io/
