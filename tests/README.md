# Testing

This directory contains unit tests for the single-cell RNA-seq discovery platform.

## Test Structure

- `test_preprocessing.py` - Tests for data preprocessing (QC, normalization, gene selection)
- `test_clustering.py` - Tests for clustering analysis (PCA, UMAP, Leiden)
- `test_markers.py` - Tests for marker gene detection
- `conftest.py` - Pytest configuration and shared fixtures

## Running Tests

### Run all tests
```bash
pytest tests/ -v
```

### Run specific test file
```bash
pytest tests/test_preprocessing.py -v
```

### Run specific test function
```bash
pytest tests/test_preprocessing.py::test_filter_cells -v
```

### Run with coverage report
```bash
pytest tests/ --cov=src --cov-report=html
```

### Run with verbose output
```bash
pytest tests/ -vv
```

## Test Coverage

Current test coverage:
- **preprocessing.py**: 5 tests covering filtering, normalization, gene selection
- **clustering.py**: 5 tests covering PCA, neighbors, UMAP, Leiden
- **markers.py**: 7 tests covering ranking, extraction, export

## Fixtures

### `sample_adata`
Basic AnnData object with random count data (100 cells × 200 genes)

### `preprocessed_adata`
AnnData object that has been normalized and log-transformed

### `clustered_adata`
AnnData object with clustering results (leiden clusters)

### `basic_adata` (session scope)
Minimal AnnData for quick tests (50 cells × 100 genes)

## Sample Data

Generate sample data for manual testing:
```bash
python data/generate_sample_data.py
```

This creates `data/sample.h5ad` with synthetic scRNA-seq data suitable for testing the full pipeline.

## CI/CD Integration

Tests can be integrated with CI/CD systems:

```yaml
# Example GitHub Actions
- name: Run tests
  run: pytest tests/ --cov=src
```

## Contributing

When adding new functionality:
1. Write tests first (TDD approach)
2. Ensure all tests pass: `pytest tests/ -v`
3. Check coverage: `pytest tests/ --cov=src`
4. Maintain >80% code coverage
