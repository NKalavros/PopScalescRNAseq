#!/usr/bin/env python
"""
Create test dataset for foundation model testing
"""

import scanpy as sc
import numpy as np

def create_test_data():
    """Create a test .h5ad file using scanpy PBMC data"""
    print("Creating test dataset...")
    
    # Load PBMC data from scanpy
    adata = sc.datasets.pbmc68k_reduced()
    
    # Basic preprocessing
    sc.pp.filter_cells(adata, min_genes=200)
    sc.pp.filter_genes(adata, min_cells=3)
    
    # Add gene names if missing
    if "gene_name" not in adata.var.columns:
        adata.var["gene_name"] = adata.var.index.tolist()
    
    print(f"Test data shape: {adata.shape}")
    print(f"Gene names: {adata.var_names[:5].tolist()}")
    
    # Save test data
    output_path = "test_data.h5ad"
    adata.write_h5ad(output_path)
    print(f"✅ Test data saved to {output_path}")
    
    return output_path

if __name__ == "__main__":
    create_test_data()