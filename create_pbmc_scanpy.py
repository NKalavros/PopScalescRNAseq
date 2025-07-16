#!/usr/bin/env python
"""
Create PBMC dataset using scanpy for foundation model testing
"""

import scanpy as sc
import pandas as pd

def create_pbmc_data():
    """Create PBMC dataset using scanpy"""
    print("Loading PBMC dataset from scanpy...")
    
    try:
        # Try to load PBMC 3k dataset (smaller, more manageable)
        adata = sc.datasets.pbmc3k_processed()
        print(f"Loaded PBMC3k dataset: {adata.shape}")
        
    except Exception as e:
        print(f"Error loading PBMC3k: {e}")
        try:
            # Fallback to PBMC68k reduced
            adata = sc.datasets.pbmc68k_reduced()  
            print(f"Loaded PBMC68k reduced dataset: {adata.shape}")
        except Exception as e2:
            print(f"Error loading PBMC68k: {e2}")
            raise e2
    
    # Make sure gene names are in var.index
    if 'gene_symbols' in adata.var.columns:
        adata.var_names = adata.var['gene_symbols']
        adata.var_names_make_unique()
    
    print(f"PBMC data shape: {adata.shape}")
    print(f"Sample gene names: {adata.var_names[:10].tolist()}")
    print(f"Expression range: {adata.X.min():.2f} - {adata.X.max():.2f}")
    
    # Save as h5ad
    output_path = "pbmc_scanpy_data.h5ad"
    adata.write_h5ad(output_path)
    print(f"✅ PBMC data saved to {output_path}")
    
    return output_path

if __name__ == "__main__":
    create_pbmc_data()