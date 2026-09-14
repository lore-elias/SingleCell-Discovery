"""
Which genes differ between conditions? Healthy and Disease
"""

import scanpy as sc
import pandas as pd

def differential_expression(adata, groupby = "condition", reference = "Healthy", method='wilcoxon'):
    sc.tl.rank_genes_groups(adata, groupby=groupby, reference=reference, method=method, key_added = "diffrential_expression")
    return adata

def get_results(adata, group = "Disease", key = "diffrential_expression"):
    return sc.get.rank_genes_groups_df(adata, group=group, key=key)

def differential_expression_pipeline(
    adata,
    output_file=None
):

    adata = differential_expression(adata)

    df = get_results(adata)

    if output_file:
        df.to_csv(output_file, index=False)

    return adata, df