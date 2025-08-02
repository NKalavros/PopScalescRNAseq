
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

# Refactored: Merge all model embeddings into a single AnnData object per study

for study in all_files:
    print(f"Merging embeddings for {study}")
    embedding_files = all_files[study]['embeddings']
    # Always load scGPT first
    scgpt_file = [f for f in embedding_files if 'scgpt.h5ad' in f]
    if not scgpt_file:
        print(f"scGPT embedding not found for {study}, skipping.")
        continue
    base_adata = sc.read_h5ad(scgpt_file[0])
    # Subsample only once, on scGPT
    if RUN_SUBSET and base_adata.n_obs > N_OBS:
        idx = np.random.choice(base_adata.n_obs, N_OBS, replace=False)
        base_adata = base_adata[idx, :]
    print(f"Set base AnnData for {study} with shape {base_adata.shape}")
    # Lowercase all obsm keys
    base_adata.obsm = {k.lower(): v for k, v in base_adata.obsm.items()}
    # Now add other embeddings, aligning by index
    for method in EMBEDDING_METHODS:
        if method == 'scgpt':
            continue  # already loaded
        method_file = f"{method}.h5ad" if method != 'geneformer' else 'geneformer.csv'
        embedding_file = [f for f in embedding_files if method_file in f]
        file_path = embedding_file[0] if embedding_file else None
        if file_path:
            print(f"Loading embedding: {file_path}")
            if file_path.endswith('.h5ad'):
                adata_tmp = sc.read_h5ad(file_path)
            elif file_path.endswith('.csv'):
                adata_tmp = sc.read_csv(file_path)
                adata_tmp.obsm['X_geneformer'] = adata_tmp.X.copy()
            else:
                print(f"Unsupported file format for {file_path}")
                continue
            adata_tmp.obsm = {k.lower(): v for k, v in adata_tmp.obsm.items()}
            # Align cells by index (assume same order, or use intersection if needed)
            # Here, we use the intersection of cell names (index)
            common_idx = base_adata.obs_names.intersection(adata_tmp.obs_names)
            if len(common_idx) == 0:
                print(f"No overlapping cells between base AnnData and {file_path}, skipping.")
                del adata_tmp
                continue
            rep_key = f"x_{method.lower()}" if method != 'geneformer' else 'X_geneformer'
            if rep_key in adata_tmp.obsm:
                arr = adata_tmp[common_idx, :].obsm[rep_key]
                # Only add if the intersection matches the base AnnData's obs_names exactly
                if list(common_idx) == list(base_adata.obs_names):
                    base_adata.obsm[rep_key] = arr
                    print(f"Added {rep_key} to base AnnData.obsm for {len(common_idx)} cells.")
                else:
                    print(f"Cell indices do not match exactly for {rep_key} in {file_path}. Skipping this embedding.")
            else:
                print(f"Representation {rep_key} not found in obsm for {file_path}")
            del adata_tmp
    print(f"Merged AnnData for {study} with representations: {list(base_adata.obsm.keys())}")

    # --- Batch Correction: PCA (unintegrated) ---
    print(f"Running PCA (unintegrated) for {study}")

    sc.pp.normalize_total(base_adata, target_sum=1e4)
    sc.pp.log1p(base_adata)
    # Save log-normalized counts in a layer only once
    base_adata.layers['lognorm'] = base_adata.X.copy()
    sc.pp.highly_variable_genes(base_adata, n_top_genes=2000, subset=True)
    sc.pp.scale(base_adata, max_value=10)
    sc.pp.pca(base_adata, n_comps=50, use_highly_variable=True)
    print(f"PCA complete. Representation 'X_pca' added to obsm.")


    # --- Define all integration keys and readable names ---
    # Add foundation model embeddings if present
    foundation_methods = [
        ("x_scgpt", "scgpt"),
        ("x_scimilarity", "scimilarity"),
        ("x_geneformer", "geneformer"),
        ("x_scfoundation", "scfoundation"),
        ("x_uce", "uce"),
    ]
    integration_methods = [
        ("X_pca", "pca"),
        ("X_pca_harmony", "harmony"),
        ("X_scvi", "scvi"),
        ("X_pca_fastmnn", "fastmnn"),
        ("X_scanorama", "scanorama"),
    ]
    # Add foundation model embeddings if present
    for obsm_key, method_name in foundation_methods:
        if obsm_key in base_adata.obsm:
            integration_methods.append((obsm_key, method_name))
    # bbknn does not create a new embedding, but updates neighbors graph
    # We'll treat it separately
    # --- Batch Correction: Combat (PCA on ComBat-corrected data) ---
    print(f"Running Combat integration for {study}")
    try:
        # Set X to lognorm for Combat
        base_adata.X = base_adata.layers['lognorm'].copy()
        sc.pp.combat(base_adata, key='library_id')
        # After Combat, run PCA on corrected X
        sc.pp.scale(base_adata, max_value=10)
        sc.pp.pca(base_adata, n_comps=50)
        # Save PCA as X_pca_combat (explicitly as batch correction embedding)
        base_adata.obsm['X_pca_combat'] = base_adata.obsm['X_pca'].copy()
        print(f"Combat batch correction complete. Representation 'X_pca_combat' added to obsm.")
        # Restore X to lognorm for downstream steps
        base_adata.X = base_adata.layers['lognorm'].copy()
    except Exception as e:
        print(f"Combat integration failed for {study}: {e}")

    # Always include ComBat PCA as a batch correction embedding
    if ("X_pca_combat", "combat") not in integration_methods:
        integration_methods.append(("X_pca_combat", "combat"))


    # --- For each integration, compute UMAP, plot, and rclone upload ---
    available_obsms = []
    for obsm_key, method_name in integration_methods:
        if obsm_key in base_adata.obsm:
            available_obsms.append(obsm_key)
            print(f"Running UMAP for {method_name} integration for {study}")
            try:
                sc.pp.neighbors(base_adata, use_rep=obsm_key)
                sc.tl.umap(base_adata)
                # Plot UMAP
                study_dir = study.replace('_embeddings', '')
                figures_dir = os.path.join(BASE_DIR, study_dir, 'figures')
                if not os.path.exists(figures_dir):
                    os.makedirs(figures_dir)
                save_path = os.path.join(figures_dir, f'{study_dir}_{method_name}.png')
                sc.pl.embedding(base_adata, basis='umap', color=['cell_type'], save=f"_{study_dir}_{method_name}.png", show=False)
                # rclone upload
                umaps_gdrive = f'GDrive:KidneyAtlas/{study_dir}/umaps/'
                local_umap_path = os.path.join('figures', f'umap_{study_dir}_{method_name}.png')
                subprocess.run(['rclone', 'copy', local_umap_path, umaps_gdrive])
                print(f"UMAP for {method_name} saved and uploaded.")
            except Exception as e:
                print(f"UMAP/plot/rclone failed for {method_name} in {study}: {e}")

    # --- bbKNN: neighbors/UMAP/plot/rclone/scIB ---

    print(f"Running UMAP for bbknn integration for {study}")
    try:
        sc.tl.umap(base_adata)
        study_dir = study.replace('_embeddings', '')
        figures_dir = os.path.join(BASE_DIR, study_dir, 'figures')
        if not os.path.exists(figures_dir):
            os.makedirs(figures_dir)
        save_path = os.path.join(figures_dir, f'{study_dir}_bbknn.png')
        sc.pl.embedding(base_adata, basis='umap', color=['cell_type'], save=f"_{study_dir}_bbknn.png", show=False)
        umaps_gdrive = f'GDrive:KidneyAtlas/{study_dir}/umaps/'
        local_umap_path = os.path.join('figures', f'umap_{study_dir}_bbknn.png')
        subprocess.run(['rclone', 'copy', local_umap_path, umaps_gdrive])
        print(f"UMAP for bbknn saved and uploaded.")
        available_obsms.append("X_umap")
    except Exception as e:
        print(f"UMAP/plot/rclone failed for bbknn in {study}: {e}")

    # --- Single scIB evaluation for all embeddings ---
    if RUN_SCIB and available_obsms:
        try:
            biocons = BioConservation(isolated_labels=False)
            bm = Benchmarker(
                base_adata,
                batch_key="library_id",
                label_key="cell_type",
                embedding_obsm_keys=available_obsms,
                pre_integrated_embedding_obsm_key="X_pca",
                bio_conservation_metrics=biocons,
                batch_correction_metrics=BatchCorrection(),
                n_jobs=-1,
            )
            bm.prepare()
            bm.benchmark()
            bm.plot_results_table(show=False, min_max_scale=False)
            print(f"scIB evaluation complete for all integrations in {study}.")
        except Exception as e:
            print(f"scIB evaluation failed for all integrations in {study}: {e}")


    # --- Batch Correction: Harmony ---
    print(f"Running Harmony integration for {study}")
    try:
        sc.external.pp.harmony_integrate(base_adata, key='library_id', basis='X_pca', max_iter_harmony=20, verbose=False)
        print(f"Harmony integration complete. Representation 'X_pca_harmony' added to obsm.")
    except Exception as e:
        print(f"Harmony integration failed for {study}: {e}")

    # --- Batch Correction: SCVI ---
    print(f"Running SCVI integration for {study}")
    try:
        import scvi
        scvi.model.SCVI.setup_anndata(base_adata, batch_key="library_id")
        scvi_model = scvi.model.SCVI(base_adata)
        scvi_model.train()
        base_adata.obsm["X_scvi"] = scvi_model.get_latents()
        print(f"SCVI integration complete. Representation 'X_scvi' added to obsm.")
    except Exception as e:
        print(f"SCVI integration failed for {study}: {e}")

    # --- Batch Correction: FastMNN ---
    print(f"Running FastMNN integration for {study}")
    try:
        import scanpy.external as sce
        # Split AnnData by batch
        batches = [base_adata[base_adata.obs['library_id'] == b].copy() for b in base_adata.obs['library_id'].unique()]
        mnn_corrected, *_ = sce.pp.mnn_correct(*batches, batch_key="library_id")
        # mnn_corrected is a concatenated AnnData object; extract the corrected PCA
        if hasattr(mnn_corrected, 'obsm') and 'X_pca' in mnn_corrected.obsm:
            base_adata.obsm['X_pca_fastmnn'] = mnn_corrected.obsm['X_pca']
            print(f"FastMNN integration complete. Representation 'X_pca_fastmnn' added to obsm.")
        else:
            print(f"FastMNN did not return expected 'X_pca' in obsm.")
    except Exception as e:
        print(f"FastMNN integration failed for {study}: {e}")

    # --- Batch Correction: bbKNN ---
    print(f"Running bbKNN integration for {study}")
    try:
        import bbknn
        bbknn.bbknn(base_adata, batch_key='library_id')
        print(f"bbKNN integration complete. Neighbors graph updated for batch correction.")
    except Exception as e:
        print(f"bbKNN integration failed for {study}: {e}")

    # --- Batch Correction: Scanorama ---
    print(f"Running Scanorama integration for {study}")
    try:
        import scanpy.external as sce
        sce.pp.scanorama_integrate(base_adata, 'library_id', basis='X_pca', adjusted_basis='X_scanorama')
        print(f"Scanorama integration complete. Representation 'X_scanorama' added to obsm.")
    except Exception as e:
        print(f"Scanorama integration failed for {study}: {e}")