"""
Differential expression analysis (Wilcoxon, t-test)

Identify top marker genes per cluster

Statistical significance filtering

Generate marker gene tables
"""

import numpy as np
import pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests

def differential_expression(data, labels, method='wilcoxon'):
    """
    Perform differential expression analysis between clusters.

    Parameters:
    data (pd.DataFrame): Gene expression data (genes x samples).
    labels (pd.Series): Cluster labels for each sample.
    method (str): Statistical test to use ('wilcoxon' or 't-test').

    Returns:
    pd.DataFrame: DataFrame containing marker genes and statistics.
    """
    results = []
    unique_clusters = labels.unique()
    
    for cluster in unique_clusters:
        cluster_samples = data.columns[labels == cluster]
        other_samples = data.columns[labels != cluster]
        
        for gene in data.index:
            if method == 'wilcoxon':
                stat, p_value = stats.ranksums(data.loc[gene, cluster_samples], data.loc[gene, other_samples])
            elif method == 't-test':
                stat, p_value = stats.ttest_ind(data.loc[gene, cluster_samples], data.loc[gene, other_samples])
            else:
                raise ValueError("Method must be 'wilcoxon' or 't-test'")
            
            results.append({'gene': gene, 'cluster': cluster, 'statistic': stat, 'p_value': p_value})
    
    results_df = pd.DataFrame(results)
    
    # Adjust p-values for multiple testing
    results_df['adjusted_p_value'] = multipletests(results_df['p_value'], method='fdr_bh')[1]
    
    return results_df

def filter_markers(results_df, p_value, threshold=0.05):
    """
    Filter marker genes based on adjusted p-value.

    Parameters:
    results_df (pd.DataFrame): DataFrame containing marker genes and statistics.
    p_value_threshold (float): Threshold for adjusted p-value.

    Returns:
    pd.DataFrame: Filtered DataFrame containing significant marker genes.
    """
    filtered_df = results_df[results_df['adjusted_p_value'] < threshold]
    return filtered_df.sort_values(by='adjusted_p_value')

def generate_marker_tables(filtered_df):
    """
    Generate marker gene tables for each cluster.

    Parameters:
    filtered_df (pd.DataFrame): DataFrame containing significant marker genes.

    Returns:
    dict: Dictionary of DataFrames, one for each cluster.
    """
    marker_tables = {}
    for cluster in filtered_df['cluster'].unique():
        marker_tables[cluster] = filtered_df[filtered_df['cluster'] == cluster].sort_values(by='adjusted_p_value')
    
    return marker_tables

