#!/usr/bin/env python
"""
Download PBMC10k dataset for foundation model testing
Uses 10x Genomics PBMC10k dataset with real gene symbols
"""

import numpy as np
import pandas as pd
import anndata as ad
import requests
import gzip
import os
from scipy.sparse import csr_matrix
from pathlib import Path

def download_file(url, filename):
    """Download a file from URL"""
    print(f"Downloading {filename}...")
    response = requests.get(url, stream=True)
    response.raise_for_status()
    
    with open(filename, 'wb') as f:
        for chunk in response.iter_content(chunk_size=8192):
            f.write(chunk)
    print(f"✅ Downloaded {filename}")

def download_pbmc10k():
    """Download and prepare PBMC10k dataset"""
    print("Downloading PBMC10k dataset from 10x Genomics...")
    
    # Create data directory
    data_dir = Path("data/pbmc10k")
    data_dir.mkdir(parents=True, exist_ok=True)
    
    # URLs for PBMC10k filtered data
    base_url = "https://cf.10xgenomics.com/samples/cell-exp/3.0.0/pbmc_10k_v3/pbmc_10k_v3_filtered_feature_bc_matrix"
    
    files = {
        "barcodes.tsv.gz": f"{base_url}/barcodes.tsv.gz",
        "features.tsv.gz": f"{base_url}/features.tsv.gz", 
        "matrix.mtx.gz": f"{base_url}/matrix.mtx.gz"
    }
    
    # Download files
    for filename, url in files.items():
        filepath = data_dir / filename
        if not filepath.exists():
            download_file(url, filepath)
        else:
            print(f"✅ {filename} already exists")
    
    # Read the data
    print("Reading matrix data...")
    
    # Read matrix
    from scipy.io import mmread
    with gzip.open(data_dir / "matrix.mtx.gz", 'rt') as f:
        X = mmread(f).T.tocsr()  # Transpose to get cells x genes
    
    # Read barcodes (cell names)
    with gzip.open(data_dir / "barcodes.tsv.gz", 'rt') as f:
        barcodes = [line.strip() for line in f]
    
    # Read features (gene info)
    with gzip.open(data_dir / "features.tsv.gz", 'rt') as f:
        features = []
        for line in f:
            parts = line.strip().split('\t')
            features.append({
                'gene_id': parts[0],
                'gene_symbol': parts[1],
                'feature_type': parts[2]
            })
    
    features_df = pd.DataFrame(features)
    
    # Create AnnData object
    print("Creating AnnData object...")
    adata = ad.AnnData(
        X=X,
        obs=pd.DataFrame(index=barcodes),
        var=pd.DataFrame(features_df, index=features_df['gene_symbol'])
    )
    
    # Filter to only include gene expression data (not antibody features)
    gene_mask = adata.var['feature_type'] == 'Gene Expression'
    adata = adata[:, gene_mask].copy()
    
    # Make gene symbols unique (some may be duplicated)
    adata.var_names_make_unique()
    
    print(f"PBMC10k data shape: {adata.shape}")
    print(f"Sample gene names: {adata.var_names[:10].tolist()}")
    print(f"Expression range: {adata.X.min():.2f} - {adata.X.max():.2f}")
    
    # Save as h5ad
    output_path = "pbmc10k_data.h5ad"
    adata.write_h5ad(output_path)
    print(f"✅ PBMC10k data saved to {output_path}")
    
    return output_path

if __name__ == "__main__":
    download_pbmc10k()