#!/usr/bin/env python
"""
Convert lung study MTX files to AnnData h5ad format
Processes data from MTX count matrices + CSV metadata

Handles:
- Checking for existing h5ad files (skips if present)
- Nested folder structures: {Study}/data/Data_{Study}_Lung/
- Case-insensitive gene file names (genes.txt vs Genes.txt)
- Various MTX file types (UMIcounts, TPM, counts)
"""

import os
import sys
import pandas as pd
import numpy as np
import anndata as ad
import scipy.io as io
from pathlib import Path
import glob

def find_file_case_insensitive(directory, filename):
    """
    Find a file in a directory, ignoring case sensitivity.

    Parameters:
    -----------
    directory : str
        Directory to search
    filename : str
        Filename to find (case-insensitive)

    Returns:
    --------
    str or None : Full path to file if found, None otherwise
    """
    try:
        files = os.listdir(directory)
        for f in files:
            if f.lower() == filename.lower():
                return os.path.join(directory, f)
    except Exception as e:
        print(f"Error searching directory {directory}: {e}")
    return None

def load_mtx_data(mtx_path, genes_path, cells_path):
    """
    Load MTX data with genes and cells metadata

    Parameters:
    -----------
    mtx_path : str
        Path to .mtx file
    genes_path : str
        Path to genes.txt file (case-insensitive)
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

    # Check dimensions match
    if X.shape[1] != len(gene_names):
        print(f"WARNING: Gene count mismatch! Matrix has {X.shape[1]} genes but {len(gene_names)} gene names loaded")
    if X.shape[0] != len(cell_names):
        print(f"WARNING: Cell count mismatch! Matrix has {X.shape[0]} cells but {len(cell_names)} cell names loaded")

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
        Base path to BigPurple lung studies (/gpfs/scratch/nk4167/LungAtlas)

    Returns:
    --------
    bool : True if successful or already exists, False if failed
    """
    # Define paths
    study_path = os.path.join(base_path, study_name)
    output_path = os.path.join(study_path, 'data.h5ad')

    # Check if output already exists
    if os.path.exists(output_path):
        print(f"✓ {study_name}: data.h5ad already exists, skipping")
        return True

    # Find the actual data directory (handles nested structure)
    # Possible paths: study_path/data/Data_{study_name}_Lung/ or similar variations
    data_search_paths = [
        os.path.join(study_path, 'data', f'Data_{study_name}_Lung'),
        os.path.join(study_path, 'data'),
    ]

    data_dir = None
    for candidate_path in data_search_paths:
        if os.path.isdir(candidate_path):
            # Check if this directory has the actual data files
            files = os.listdir(candidate_path)
            if any(f.endswith('.mtx') for f in files):
                data_dir = candidate_path
                break

    if data_dir is None:
        print(f"✗ {study_name}: Could not find data directory with MTX files")
        print(f"  Checked paths: {data_search_paths}")
        return False

    print(f"\n{'='*60}")
    print(f"Processing {study_name}")
    print(f"Data directory: {data_dir}")
    print(f"{'='*60}")

    try:
        # Find MTX file (could be Exp_data_UMIcounts.mtx, Exp_data_TPM.mtx, Exp_data_counts.mtx, etc.)
        mtx_files = [f for f in os.listdir(data_dir) if f.endswith('.mtx')]
        if not mtx_files:
            print(f"✗ {study_name}: No MTX file found")
            return False

        mtx_file = mtx_files[0]
        if len(mtx_files) > 1:
            print(f"  WARNING: Found multiple MTX files, using {mtx_file}")
        mtx_path = os.path.join(data_dir, mtx_file)
        print(f"  Using MTX file: {mtx_file}")

        # Find genes file (case-insensitive: genes.txt or Genes.txt)
        genes_path = find_file_case_insensitive(data_dir, 'genes.txt')
        if genes_path is None:
            print(f"✗ {study_name}: Genes file not found (genes.txt or Genes.txt)")
            return False
        print(f"  Using genes file: {os.path.basename(genes_path)}")

        # Find cells file (case-insensitive)
        cells_path = find_file_case_insensitive(data_dir, 'cells.csv')
        if cells_path is None:
            print(f"✗ {study_name}: Cells.csv file not found")
            return False
        print(f"  Using cells file: {os.path.basename(cells_path)}")

        # Load data
        adata = load_mtx_data(mtx_path, genes_path, cells_path)

        # Create output directory if needed
        os.makedirs(os.path.dirname(output_path), exist_ok=True)

        # Save to h5ad
        print(f"Saving to {output_path}")
        adata.write_h5ad(output_path)

        print(f"✓ Successfully created {output_path}")
        return True

    except Exception as e:
        print(f"✗ {study_name}: Error during processing: {str(e)}")
        import traceback
        traceback.print_exc()
        return False

def main():
    """Process all lung studies"""

    # BigPurple path
    base_path = '/gpfs/scratch/nk4167/LungAtlas'

    # All studies with data files
    studies = [
        'Bischoff2021',
        'Chan2021',
        'Guo2018',
        'Ireland2020',
        'Kim2020',
        'Laughney2020',
        'Maynard2020',
        'Song2019',
        'Xing2021',
        'Zilionis2019',
        'Qian2020',  # Has data now
    ]

    print(f"Processing lung studies from {base_path}\n")

    results = {}
    skipped = 0
    processed = 0
    failed = 0

    for study in studies:
        success = process_lung_study(study, base_path)
        if success:
            # Check if we actually created it or skipped it
            study_path = os.path.join(base_path, study)
            if os.path.exists(os.path.join(study_path, 'data.h5ad')):
                results[study] = 'SUCCESS'
                processed += 1
            else:
                results[study] = 'SKIPPED (already exists)'
                skipped += 1
        else:
            results[study] = 'FAILED'
            failed += 1

    print(f"\n{'='*60}")
    print("Summary:")
    print(f"{'='*60}")
    for study, status in sorted(results.items()):
        print(f"{study:20s}: {status}")

    print(f"\nTotal: {len(studies)} studies")
    print(f"  Created new: {processed}")
    print(f"  Skipped (already exist): {skipped}")
    print(f"  Failed: {failed}")
    print(f"{'='*60}")

if __name__ == '__main__':
    main()
