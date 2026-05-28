# Single-cell RNA-seq Discovery Platform

## Workflow 

### Input 

Public single-cell RNA-seq dataset (from CellxGene)

Adapt later to user's own data

### Pipeline

1. Download Dataset
2. Clean/process cells
3. Normalize counts
4. Cluster cells
5. Visualize clusters
6. Detect marker genes
7. Predict cell types
8. Compare healthy vs diseased tissue
9. Generate biological insights automatically

### Output

- UMAP plots
- discovered cell populations
- differenntial expression tables
- automated report
- interactive dashboard

Saved in `output/` directory

## How to use

**1) One command:**

```bash
python run_pipeline.py
``` 

and everything runs

**2) Dashboard:**

User uploads data, sees:

- clusters
- UMAP 
- Marker genes


**Includes: Docker**

**Discovery platform:**
- "Novel immune states in breast cancer?"
- "Cell-type changes in Alzheimer's disease?"
- "Inflammatory signatures in autoimmune disorders?"