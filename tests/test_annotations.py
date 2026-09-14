import pytest
from src.annotation import (
    annotate_cell_types,
    annotate_disease_association,
    annotation_pipeline
)

def test_annotate_cell_types(sample_adata):
    """Test cell type annotation."""
    adata = sample_adata.copy()
    
    adata_annotated = annotate_cell_types(adata)
    
    assert 'cell_type' in adata_annotated.obs
    assert adata_annotated.obs['cell_type'].notna().all()

def test_annotate_disease_association(sample_adata):
    """Test disease association annotation."""
    adata = sample_adata.copy()
    
    adata_disease = annotate_disease_association(adata)
    
    assert 'disease_association' in adata_disease.obs
    assert adata_disease.obs['disease_association'].notna().all()

def test_annotation_pipeline(sample_adata):
    """Test complete annotation pipeline."""
    adata = sample_adata.copy()
    
    adata_annotated = annotation_pipeline(adata)
    
    assert 'cell_type' in adata_annotated.obs
    assert 'disease_association' in adata_annotated.obs
    assert adata_annotated.obs['cell_type'].notna().all()
    assert adata_annotated.obs['disease_association'].notna().all()

def test_annotation_consistency(sample_adata):
    """Test consistency of annotations."""
    adata = sample_adata.copy()
    
    adata_annotated = annotation_pipeline(adata)
    
    # Check that cell types are consistent with known markers
    cell_type_markers = {
        'T_cell': ['CD3D', 'CD3E', 'IL7R'],
        'B_cell': ['MS4A1', 'CD79A'],
        'Monocyte': ['LYZ', 'FCN1', 'S100A8'],
        'NK_cell': ['NKG7', 'GNLY']
    }
    
    for cell_type, markers in cell_type_markers.items():
        cells_of_type = adata_annotated.obs['cell_type'] == cell_type
        for marker in markers:
            assert adata_annotated[cells_of_type, marker].X.mean() > adata_annotated[~cells_of_type, marker].X.mean()

