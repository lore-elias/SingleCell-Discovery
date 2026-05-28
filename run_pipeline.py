#!/usr/bin/env python3
"""
Main pipeline runner for single-cell RNA-seq discovery platform.

Usage:
    python run_pipeline.py --input data.h5ad --output results/
"""

import argparse
import sys
from pathlib import Path

from src import (
    load_h5ad,
    preprocess_pipeline,
    clustering_pipeline,
    marker_gene_pipeline,
    visualization
)


def main():
    """Run the complete pipeline."""
    parser = argparse.ArgumentParser(
        description="Single-cell RNA-seq discovery platform"
    )
    parser.add_argument(
        '--input',
        type=str,
        required=True,
        help='Path to H5AD input file'
    )
    parser.add_argument(
        '--output',
        type=str,
        default='results/',
        help='Output directory (default: results/)'
    )
    parser.add_argument(
        '--min-genes',
        type=int,
        default=200,
        help='Minimum genes per cell (default: 200)'
    )
    parser.add_argument(
        '--max-genes',
        type=int,
        default=2500,
        help='Maximum genes per cell (default: 2500)'
    )
    parser.add_argument(
        '--max-mito',
        type=float,
        default=0.05,
        help='Maximum mitochondrial percentage (default: 0.05)'
    )
    parser.add_argument(
        '--n-genes',
        type=int,
        default=2000,
        help='Number of variable genes (default: 2000)'
    )
    parser.add_argument(
        '--n-pcs',
        type=int,
        default=50,
        help='Number of PCA components (default: 50)'
    )
    parser.add_argument(
        '--resolution',
        type=float,
        default=1.0,
        help='Leiden resolution (default: 1.0)'
    )
    parser.add_argument(
        '--n-markers',
        type=int,
        default=10,
        help='Top marker genes per cluster (default: 10)'
    )
    
    args = parser.parse_args()
    
    # Validate inputs
    input_path = Path(args.input)
    if not input_path.exists():
        print(f"Error: Input file not found: {args.input}")
        sys.exit(1)
    
    output_dir = Path(args.output)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    print("=" * 60)
    print("Single-cell RNA-seq Discovery Platform")
    print("=" * 60)
    
    # Load data
    print("\n[1/4] Loading data...")
    adata = load_h5ad(str(input_path))
    
    # Preprocessing
    print("\n[2/4] Preprocessing...")
    adata = preprocess_pipeline(
        adata,
        min_genes=args.min_genes,
        max_genes=args.max_genes,
        max_mito=args.max_mito,
        n_top_genes=args.n_genes
    )
    
    # Clustering
    print("\n[3/4] Clustering...")
    adata = clustering_pipeline(
        adata,
        n_pcs=args.n_pcs,
        resolution=args.resolution
    )
    
    # Marker gene detection
    print("\n[4/4] Detecting marker genes...")
    adata, markers = marker_gene_pipeline(
        adata,
        n_genes=args.n_markers,
        output_dir=str(output_dir / 'markers')
    )
    
    # Save results
    print("\nSaving results...")
    adata.write(output_dir / 'processed_data.h5ad')
    print(f"Saved: {output_dir}/processed_data.h5ad")
    
    # Generate plots
    print("\nGenerating visualizations...")
    viz_dir = output_dir / 'plots'
    viz_dir.mkdir(exist_ok=True)
    
    # UMAP by clusters
    try:
        umap_fig = visualization.plot_umap_interactive(
            adata,
            color_by='leiden',
            title='UMAP: Leiden Clusters'
        )
        if umap_fig:
            umap_fig.write_html(viz_dir / 'umap_clusters.html')
            print(f"Saved: {viz_dir}/umap_clusters.html")
    except Exception as e:
        print(f"Warning: Could not create UMAP plot: {e}")
    
    print("\n" + "=" * 60)
    print("Pipeline complete!")
    print(f"Results saved to: {output_dir}")
    print("=" * 60)


if __name__ == '__main__':
    main()
