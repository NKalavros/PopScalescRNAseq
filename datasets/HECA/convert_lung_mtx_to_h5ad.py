#!/usr/bin/env python
"""
Convert lung study MTX files to AnnData h5ad format
Processes data from MTX count matrices + CSV metadata
"""

import os
import sys
import pandas as pd
import numpy as np
import anndata as ad
import scipy.io as io
from pathlib import Path

def load_mtx_data(mtx_path, genes_path, cells_path):
    """
    Load MTX data with genes and cells metadata

    Parameters:
    -----------
    mtx_path : str
        Path to .mtx file
    genes_path : str
        Path to genes.txt file
    cells_path : str
        Path to cells.csv file

    Returns:
    --------
    adata : AnnData object
    """
    print(f"Loading MTX from {mtx_path}")

    # Load count matrix
    X = io.mmread(mtx_path).T.tocsr()  # Transpose to (cells, genes)

    # Load gene names
    print(f"Loading genes from {genes_path}")
    genes = pd.read_csv(genes_path, header=None, sep='\t', names=['gene_name'])
    gene_names = genes['gene_name'].values

    # Load cell metadata
    print(f"Loading cells from {cells_path}")
    cells_df = pd.read_csv(cells_path)
    cell_names = cells_df['cell_name'].values if 'cell_name' in cells_df.columns else cells_df.iloc[:, 0].values

    print(f"Data shape: {X.shape}")
    print(f"Genes: {len(gene_names)}, Cells: {len(cell_names)}")

    # Create AnnData object
    adata = ad.AnnData(
        X=X,
        obs=cells_df.set_index(cells_df.columns[0]) if len(cells_df.columns) > 0 else pd.DataFrame(index=cell_names),
        var=pd.DataFrame(index=gene_names)
    )

    adata.obs_names.name = None
    adata.var_names.name = None

    return adata

def process_lung_study(study_name, base_path):
    """
    Process a single lung study from MTX to h5ad

    Parameters:
    -----------
    study_name : str
        Name of the study (e.g., 'Laughney2020')
    base_path : str
        Base path to BigPurple lung studies
    """
    # Define paths based on study structure
    study_path = os.path.join(base_path, study_name)
    data_dir = os.path.join(study_path, 'data', f'Data_{study_name}_Lung')
    output_path = os.path.join(study_path, 'data.h5ad')

    # Check if data directory exists
    if not os.path.exists(data_dir):
        print(f"WARNING: Data directory not found for {study_name}: {data_dir}")
        return False

    # Check if output already exists
    if os.path.exists(output_path):
        print(f"Output file already exists for {study_name}: {output_path}")
        return True

    try:
        # Find MTX file (could be different names)
        mtx_files = [f for f in os.listdir(data_dir) if f.endswith('.mtx')]
        if not mtx_files:
            print(f"ERROR: No MTX file found in {data_dir}")
            return False
        mtx_file = mtx_files[0]
        mtx_path = os.path.join(data_dir, mtx_file)

        # Find genes file
        genes_file = 'Genes.txt'
        genes_path = os.path.join(data_dir, genes_file)
        if not os.path.exists(genes_path):
            print(f"ERROR: Genes file not found at {genes_path}")
            return False

        # Find cells file
        cells_file = 'Cells.csv'
        cells_path = os.path.join(data_dir, cells_file)
        if not os.path.exists(cells_path):
            print(f"ERROR: Cells file not found at {cells_path}")
            return False

        print(f"\n{'='*60}")
        print(f"Processing {study_name}")
        print(f"{'='*60}")

        # Load data
        adata = load_mtx_data(mtx_path, genes_path, cells_path)

        # Create output directory if needed
        os.makedirs(os.path.dirname(output_path), exist_ok=True)

        # Save to h5ad
        print(f"Saving to {output_path}")
        adata.write_h5ad(output_path)

        print(f"Successfully created {output_path}")
        return True

    except Exception as e:
        print(f"ERROR processing {study_name}: {str(e)}")
        import traceback
        traceback.print_exc()
        return False

def main():
    """Process all lung studies"""

    # BigPurple path
    base_path = '/gpfs/scratch/nk4167/LungAtlas'

    # Studies to process (excluding Qian2020 which has metadata only)
    studies = [
        'Laughney2020',
        'Maynard2020',
        'Song2019',
        'Xing2021',
        'Zilionis2019'
    ]

    print(f"Processing lung studies from {base_path}\n")

    results = {}
    for study in studies:
        success = process_lung_study(study, base_path)
        results[study] = 'SUCCESS' if success else 'FAILED'

    print(f"\n{'='*60}")
    print("Summary:")
    print(f"{'='*60}")
    for study, status in results.items():
        print(f"{study}: {status}")

if __name__ == '__main__':
    main()
