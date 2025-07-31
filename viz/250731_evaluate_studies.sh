
import os
import scanpy as sc
import pandas as pd
import numpy as np
import anndata as ad
import scib
import sklearn
import time
from scipy import sparse
from scib_metrics.benchmark import Benchmarker, BioConservation, BatchCorrection
import scvi
# Import subprocess
import subprocess
# Set base directory (relative to script location)
BASE_DIR = '/gpfs/scratch/nk4167/KidneyAtlas'  # Change this to your base directory


N_OBS = 50000  # Number of observations to subsample for speed
EMBEDDING_METHODS = ['scgpt', 'scimilarity', 'geneformer', 'scfoundation', 'uce']  # Embedding methods to evaluate
DIRECTORIES_TO_EVALUATE = ['lake_scrna','lake_snrna', 'Abedini', 'SCP1288', 'Krishna', 'Braun']  # Directories to evaluate
RUN_SCIB = True  # Whether to run scib evaluation
RUN_SUBSET = True  # Whether to subsample the data
RUN_UMAP = True  # Whether to run UMAP
all_files = {}
for dir in DIRECTORIES_TO_EVALUATE:
    print(f"Processing directory: {dir}")
    figures_dir = os.path.join(BASE_DIR, dir, 'figures')
    umaps_dir = os.path.join(BASE_DIR, dir, 'umaps')
    if not os.path.exists(figures_dir):
        os.makedirs(figures_dir)
    if not os.path.exists(umaps_dir):
        os.makedirs(umaps_dir)
    # Get all subdirectories that end with 'embeddings' (not just containing 'embeddings')
    dir_path = os.path.join(BASE_DIR, dir)
    subdirs = [d for d in os.listdir(dir_path)
              if os.path.isdir(os.path.join(dir_path, d)) and d.endswith('embeddings') and not d.startswith('.')]
    # Gather all the embedding files for each subdirectory
    for subdir in subdirs:
        print(f"Processing subdirectory: {subdir}")
        embedding_files = {}
        embedding_files[subdir] = []
        for method in EMBEDDING_METHODS:
            # Geneformer generates geneformer.csv, so we need to check for that
            if method.lower() == 'geneformer':
                file_path = os.path.join(dir_path, subdir, 'geneformer.csv')
            else:
                file_path = os.path.join(dir_path, subdir, f'{method}.h5ad')
            if os.path.exists(file_path):
                embedding_files[subdir].append(file_path)
        all_files[dir + "_" + subdir] = embedding_files
# Let's go through each files now and evaluate them
for method in EMBEDDING_METHODS:
    for study in all_files:
        print(f"Evaluating {method} for {study}")
        embedding_files = all_files[study]['embeddings']
        # Check for each method in the keys for the embedding files with an in (because it's not the actual name of the method but method.h5ad or method.csv)
        # This is a bit tricky because the keys are not just the method names but also the subdirectory names
        # So we need to check if the method is in the keys of the embedding_files
        # and then process the files accordingly
        # If the method is in the embedding_files, we process those files
        method_file = f"{method}.h5ad" if method != 'geneformer' else 'geneformer.csv'
        embedding_file = [f for f in embedding_files if method_file in f]
        file_path = embedding_file[0] if embedding_file else None
        if file_path:
            print(f"Processing file: {file_path}")
            # Load the data
            if file_path.endswith('.h5ad'):
                adata = sc.read_h5ad(file_path)
            elif file_path.endswith('.csv'):
                adata = sc.read_csv(file_path)

                adata.obsm['X_geneformer'] = adata.X.copy()
            else:
                print(f"Unsupported file format for {file_path}")
                continue
            print(f"Loaded {file_path} with shape {adata.shape}")
            # Print existing obsms
            print(f"Available obsm keys: {adata.obsm.keys()}")
            # Rename all representations to lowercase
            adata.obsm = {k.lower(): v for k, v in adata.obsm.items()}
            # Print all obs columns we have
            # Subsample if needed
            if RUN_SUBSET and adata.n_obs > N_OBS:
                adata = adata[np.random.choice(adata.n_obs, N_OBS, replace=False), :]
            # Get the name of the desired representation
            representation_key = f"x_{method.lower()}" if method != 'geneformer' else 'X_geneformer'
            if representation_key not in adata.obsm:
                print(f"Representation {representation_key} not found in obsm. Skipping {method} for {study}.")
                continue
            # Run UMAP if needed
            if RUN_UMAP:
                sc.pp.neighbors(adata, n_neighbors=15, use_rep=representation_key)
                sc.tl.umap(adata)
            # Plot UMAP and save to file in the correct figures directory
            # Determine the correct figures_dir for this study/subdir
            study_dir = study.replace('_embeddings', '')  # Remove '_embeddings' from the study name

            figures_dir = os.path.join(BASE_DIR, study_dir, 'figures')
            if not os.path.exists(figures_dir):
                os.makedirs(figures_dir)
            save_path = os.path.join(figures_dir, f'{study_dir}_{method}.png')
            sc.pl.embedding(adata, basis='umap', color=['cell_type'], save=study_dir + '_' +  method, show=False)
            # Upload with rclone to GDrive
            umaps_gdrive = f'GDrive:KidneyAtlas/{study_dir}/umaps/'
            subprocess.run(['rclone', 'copy', 'figures/umap' + study_dir + '_' + method + '.pdf', umaps_gdrive])
            # Run preprocessing and PCA (unintegrated)
            sc.pp.normalize_total(adata, target_sum=1e4)
            sc.pp.log1p(adata)
            # Highly variable genes
            sc.pp.highly_variable_genes(adata, n_top_genes=2000, subset=True)
            sc.pp.scale(adata, max_value=10)
            sc.pp.pca(adata, n_comps=50, use_highly_variable=True)
            # Run quick harmony for comparison 
            sc.external.pp.harmony_integrate(adata, key='library_id', basis='X_pca', max_iter_harmony=20, verbose=False)
            # Run quick scVI for comparison
            if 'X_scvi' not in adata.obsm:
            # Run scib evaluation
            if RUN_SCIB:
                biocons = BioConservation(isolated_labels=False)
                bm = Benchmarker(
                    adata,
                    batch_key="library_id",
                    label_key="cell_type",
                    embedding_obsm_keys=["X_umap"],
                    pre_integrated_embedding_obsm_key="X_pca",
                    bio_conservation_metrics=biocons,
                    batch_correction_metrics=BatchCorrection(),
                    n_jobs=-1,
                )
                bm.prepare()
                bm.benchmark()
                bm.plot_results_table(show=False, min_max_scale=False)