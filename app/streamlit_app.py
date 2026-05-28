#!/usr/bin/env python3
"""
Interactive Streamlit dashboard for single-cell RNA-seq discovery platform.

Run with: streamlit run app/streamlit_app.py
"""

import streamlit as st
import pandas as pd
import numpy as np
import plotly.graph_objects as go
from pathlib import Path
import tempfile
import os

# Set page config
st.set_page_config(
    page_title="scRNA-seq Discovery Platform",
    page_icon="🧬",
    layout="wide",
    initial_sidebar_state="expanded"
)

# Import pipeline functions
try:
    from src import (
        load_h5ad,
        preprocess_pipeline,
        clustering_pipeline,
        marker_gene_pipeline
    )
    from src.visualization import plot_umap_interactive, plot_heatmap_interactive
    import scanpy as sc
    HAS_IMPORTS = True
except ImportError as e:
    HAS_IMPORTS = False
    import_error = str(e)


def main():
    """Main Streamlit app."""
    
    st.title("🧬 Single-cell RNA-seq Discovery Platform")
    st.markdown("Interactive analysis and visualization of scRNA-seq data")
    
    # Check imports
    if not HAS_IMPORTS:
        st.error(f"Missing imports: {import_error}")
        st.info("Install required packages: `pip install -r requirements.txt`")
        return
    
    # Sidebar navigation
    st.sidebar.header("Navigation")
    page = st.sidebar.radio(
        "Select page:",
        ["Home", "Upload & Process", "Analysis", "Visualization", "Results"]
    )
    
    if page == "Home":
        show_home()
    elif page == "Upload & Process":
        show_upload()
    elif page == "Analysis":
        show_analysis()
    elif page == "Visualization":
        show_visualization()
    elif page == "Results":
        show_results()


def show_home():
    """Home page."""
    st.header("Welcome to scRNA-seq Discovery")
    
    col1, col2 = st.columns(2)
    
    with col1:
        st.subheader("🚀 Quick Start")
        st.markdown("""
        1. **Upload Data** - Load your H5AD file
        2. **Preprocess** - QC filtering and normalization
        3. **Cluster** - PCA, UMAP, and Leiden clustering
        4. **Discover** - Identify marker genes
        5. **Visualize** - Explore results interactively
        """)
    
    with col2:
        st.subheader("📊 Features")
        st.markdown("""
        - Quality control filtering
        - Normalization & gene selection
        - PCA dimensionality reduction
        - UMAP visualization
        - Leiden clustering
        - Marker gene detection
        - Interactive plots
        - Data export
        """)
    
    st.divider()
    st.subheader("Example Usage")
    st.code("""
# Load and process data
adata = load_h5ad("data.h5ad")
adata = preprocess_pipeline(adata)
adata = clustering_pipeline(adata)
adata, markers = marker_gene_pipeline(adata)
    """, language="python")


def show_upload():
    """Upload and preprocessing page."""
    st.header("📤 Upload & Preprocessing")
    
    col1, col2 = st.columns(2)
    
    with col1:
        st.subheader("Upload Data")
        uploaded_file = st.file_uploader(
            "Choose H5AD file",
            type=["h5ad"],
            help="Single-cell RNA-seq data in H5AD format"
        )
        
        if uploaded_file:
            st.success(f"✅ File uploaded: {uploaded_file.name}")
            
            # Save to temp file
            with tempfile.NamedTemporaryFile(suffix=".h5ad", delete=False) as tmp:
                tmp.write(uploaded_file.getbuffer())
                tmp_path = tmp.name
            
            # Load data
            try:
                adata = load_h5ad(tmp_path)
                st.session_state.adata = adata
                
                st.info(f"""
                **Data Loaded:**
                - Cells: {adata.shape[0]:,}
                - Genes: {adata.shape[1]:,}
                """)
            except Exception as e:
                st.error(f"Error loading file: {e}")
            finally:
                os.unlink(tmp_path)
    
    with col2:
        st.subheader("Preprocessing Parameters")
        
        min_genes = st.slider(
            "Min genes per cell",
            min_value=50,
            max_value=500,
            value=200,
            help="Cells with fewer genes are filtered out"
        )
        
        max_genes = st.slider(
            "Max genes per cell",
            min_value=1000,
            max_value=10000,
            value=2500,
            help="Cells with more genes are filtered out (removes doublets)"
        )
        
        max_mito = st.slider(
            "Max mitochondrial %",
            min_value=0.01,
            max_value=0.2,
            value=0.05,
            step=0.01,
            help="Cells with higher mitochondrial content are filtered out"
        )
        
        n_top_genes = st.slider(
            "Highly variable genes",
            min_value=500,
            max_value=5000,
            value=2000,
            step=100,
            help="Number of genes to keep after HVG selection"
        )
    
    st.divider()
    
    # Run preprocessing
    if "adata" in st.session_state and st.button("▶️ Run Preprocessing", type="primary"):
        with st.spinner("Preprocessing data..."):
            try:
                adata = st.session_state.adata.copy()
                adata = preprocess_pipeline(
                    adata,
                    min_genes=min_genes,
                    max_genes=max_genes,
                    max_mito=max_mito,
                    n_top_genes=n_top_genes
                )
                st.session_state.adata_preprocessed = adata
                
                st.success("✅ Preprocessing complete!")
                st.info(f"""
                **Preprocessed Data:**
                - Cells: {adata.shape[0]:,}
                - Genes: {adata.shape[1]:,}
                """)
            except Exception as e:
                st.error(f"Error in preprocessing: {e}")


def show_analysis():
    """Analysis page (clustering)."""
    st.header("🔬 Clustering Analysis")
    
    if "adata_preprocessed" not in st.session_state:
        st.warning("⚠️ Please preprocess data first on the 'Upload & Process' page")
        return
    
    col1, col2 = st.columns(2)
    
    with col1:
        n_pcs = st.slider(
            "PCA components",
            min_value=10,
            max_value=100,
            value=50,
            step=5
        )
        
        n_neighbors = st.slider(
            "Number of neighbors",
            min_value=5,
            max_value=50,
            value=15,
            step=1
        )
    
    with col2:
        min_dist = st.slider(
            "UMAP min distance",
            min_value=0.01,
            max_value=0.5,
            value=0.1,
            step=0.01
        )
        
        resolution = st.slider(
            "Leiden resolution",
            min_value=0.1,
            max_value=3.0,
            value=1.0,
            step=0.1,
            help="Higher values produce more clusters"
        )
    
    st.divider()
    
    if st.button("▶️ Run Clustering", type="primary"):
        with st.spinner("Performing clustering..."):
            try:
                adata = st.session_state.adata_preprocessed.copy()
                adata = clustering_pipeline(
                    adata,
                    n_pcs=n_pcs,
                    n_neighbors=n_neighbors,
                    min_dist=min_dist,
                    resolution=resolution
                )
                st.session_state.adata_clustered = adata
                
                n_clusters = len(adata.obs['leiden'].unique())
                st.success(f"✅ Clustering complete! Found {n_clusters} clusters")
                
            except Exception as e:
                st.error(f"Error in clustering: {e}")


def show_visualization():
    """Visualization page."""
    st.header("📊 Visualizations")
    
    if "adata_clustered" not in st.session_state:
        st.warning("⚠️ Please complete clustering first")
        return
    
    adata = st.session_state.adata_clustered
    
    tab1, tab2, tab3 = st.tabs(["UMAP", "Heatmap", "Gene Expression"])
    
    with tab1:
        st.subheader("UMAP Embedding")
        
        color_by = st.selectbox(
            "Color by:",
            ["leiden"] + [col for col in adata.obs.columns if col != "leiden"]
        )
        
        try:
            fig = plot_umap_interactive(adata, color_by=color_by)
            if fig:
                st.plotly_chart(fig, use_container_width=True)
        except Exception as e:
            st.error(f"Error creating UMAP: {e}")
    
    with tab2:
        st.subheader("Marker Gene Heatmap")
        
        # Detect markers if not done
        if "markers" not in st.session_state:
            if st.button("Detect markers"):
                with st.spinner("Detecting marker genes..."):
                    try:
                        adata_temp, markers = marker_gene_pipeline(adata, n_genes=10)
                        st.session_state.adata_clustered = adata_temp
                        st.session_state.markers = markers
                        st.rerun()
                    except Exception as e:
                        st.error(f"Error detecting markers: {e}")
        
        if "markers" in st.session_state:
            markers = st.session_state.markers
            cluster = st.selectbox("Select cluster:", list(markers.keys()))
            
            if cluster:
                marker_genes = markers[cluster]['gene'].head(15).tolist()
                
                try:
                    fig = plot_heatmap_interactive(adata, marker_genes, groupby='leiden')
                    if fig:
                        st.plotly_chart(fig, use_container_width=True)
                except Exception as e:
                    st.error(f"Error creating heatmap: {e}")
    
    with tab3:
        st.subheader("Gene Expression")
        
        gene = st.selectbox("Select gene:", adata.var_names[:100].tolist())
        groupby = st.checkbox("Group by leiden", value=True)
        
        if gene:
            try:
                expr = adata[:, gene].X
                if hasattr(expr, 'toarray'):
                    expr = expr.toarray().flatten()
                else:
                    expr = np.asarray(expr).flatten()
                
                df = pd.DataFrame({
                    'Expression': expr,
                    'Leiden': adata.obs['leiden'].values if groupby else 'All'
                })
                
                fig = go.Figure()
                if groupby:
                    for cluster in sorted(adata.obs['leiden'].unique()):
                        cluster_expr = df[df['Leiden'] == cluster]['Expression']
                        fig.add_trace(go.Box(y=cluster_expr, name=str(cluster)))
                else:
                    fig.add_trace(go.Histogram(x=expr, name=gene))
                
                fig.update_layout(title=f"{gene} Expression", height=500)
                st.plotly_chart(fig, use_container_width=True)
                
            except Exception as e:
                st.error(f"Error plotting gene: {e}")


def show_results():
    """Results and export page."""
    st.header("💾 Results & Export")
    
    if "adata_clustered" not in st.session_state:
        st.warning("⚠️ No analysis results yet")
        return
    
    adata = st.session_state.adata_clustered
    
    col1, col2, col3 = st.columns(3)
    
    with col1:
        st.metric("Cells", f"{adata.shape[0]:,}")
    
    with col2:
        st.metric("Genes", f"{adata.shape[1]:,}")
    
    with col3:
        n_clusters = len(adata.obs['leiden'].unique())
        st.metric("Clusters", n_clusters)
    
    st.divider()
    
    st.subheader("Export Options")
    
    col1, col2, col3 = st.columns(3)
    
    with col1:
        if st.button("📥 Download H5AD", use_container_width=True):
            # Save to temporary file
            with tempfile.NamedTemporaryFile(suffix=".h5ad", delete=False) as tmp:
                adata.write(tmp.name)
                with open(tmp.name, "rb") as f:
                    st.download_button(
                        label="Download",
                        data=f.read(),
                        file_name="results.h5ad",
                        mime="application/octet-stream"
                    )
                os.unlink(tmp.name)
    
    with col2:
        if st.button("📊 Download Markers (CSV)", use_container_width=True):
            if "markers" in st.session_state:
                markers = st.session_state.markers
                all_markers = pd.concat(
                    [df.assign(cluster=c) for c, df in markers.items()]
                )
                csv = all_markers.to_csv(index=False)
                st.download_button(
                    label="Download",
                    data=csv,
                    file_name="markers.csv",
                    mime="text/csv"
                )
    
    with col3:
        if st.button("📈 Download Summary", use_container_width=True):
            summary = f"""
Single-cell RNA-seq Analysis Summary
====================================

Data:
- Cells: {adata.shape[0]:,}
- Genes: {adata.shape[1]:,}
- Clusters: {len(adata.obs['leiden'].unique())}

Parameters:
- PCA components: {adata.n_obs if hasattr(adata, 'n_obs') else 'N/A'}
- Resolution: {adata.obs['leiden'].max() if 'leiden' in adata.obs.columns else 'N/A'}
            """
            st.download_button(
                label="Download",
                data=summary,
                file_name="summary.txt",
                mime="text/plain"
            )


if __name__ == "__main__":
    main()
