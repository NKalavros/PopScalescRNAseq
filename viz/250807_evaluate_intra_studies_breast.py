
import os
import scanpy as sc # type: ignore[import]
import scib_metrics # type: ignore[import]
import pandas as pd # type: ignore[import]
import numpy as np
import anndata as ad # type: ignore[import]
import scib # type: ignore[import]
import sklearn # type: ignore[import]
import time 
from scipy import sparse # type: ignore[import]
from scib_metrics.benchmark import Benchmarker, BioConservation, BatchCorrection # type: ignore[import]
import subprocess
import gc
# Set base directory (relative to script location)
BASE_DIR = '/gpfs/scratch/nk4167/BreastAtlas'  # Change this to your base directory
N_OBS = 10000  # Number of observations to subsample for speed
EMBEDDING_METHODS = ['scgpt', 'scimilarity', 'geneformer', 'scfoundation', 'uce']  # Embedding methods to evaluate
DIRECTORIES_TO_EVALUATE = [ 'Bassez', 'Wu', 'Gray', 'Klughammer_snRNA','Klughammer_scRNA','Tietscher']  # Directories to evaluate
RUN_SCIB = True  # Whether to run scib evaluation
RUN_SUBSET = True  # Whether to subsample the data
RUN_UMAP = True  # Whether to run UMAP
LIBRARY_KEY = 'library_id'  # Key for batch/library information in obs (biosample_id,orig.ident, Sample,SampleID, library_id)
CELL_TYPE_KEY = 'cell_type'  # Key for cell type information in obs (FinalCellType, subclass.l1, subclass.l2 Major.subtype, cell_type, celltype_major,ClusterName_AllCells,cluster,Cluster_Idents)
all_files = {}
for dir in DIRECTORIES_TO_EVALUATE:
    print(f"Processing directory: {dir}")
    figures_dir = os.path.join(BASE_DIR, dir, 'figures')
    umaps_dir = os.path.join(BASE_DIR, dir, 'umaps')
    if not os.path.exists(figures_dir):
        os.makedirs(figures_dir)
    if not os.path.exists(umaps_dir):
        os.makedirs(umaps_dir)
    # Get all subdirectoriesrun_ that end with 'embeddings' (not just containing 'embeddings')
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
                file_path = os.path.join(dir_path, subdir, 'geneformer.h5ad')
            else:
                file_path = os.path.join(dir_path, subdir, f'{method}.h5ad')
            if os.path.exists(file_path):
                embedding_files[subdir].append(file_path)
        all_files[dir] = embedding_files


# --- Utility Functions ---
def run_combat(adata, lognorm_layer, batch_key='library_id', n_comps=50):
    try:
        adata.X = adata.layers[lognorm_layer].copy()
        sc.pp.combat(adata, key=batch_key)
        sc.pp.scale(adata, max_value=10)
        sc.pp.pca(adata, n_comps=n_comps)
        adata.obsm['X_pca_combat'] = adata.obsm['X_pca'].copy()
        adata.X = adata.layers[lognorm_layer].copy()
        print(f"Combat batch correction complete. Representation 'X_pca_combat' added to obsm.")
        return True
    except Exception as e:
        print(f"Combat integration failed: {e}")
        return False

def run_harmony(adata, batch_key='library_id', basis='X_pca'):
    try:
        sc.external.pp.harmony_integrate(adata, key=batch_key, basis=basis, max_iter_harmony=20, verbose=False)
        print(f"Harmony integration complete. Representation 'X_pca_harmony' added to obsm.")
        return True
    except Exception as e:
        print(f"Harmony integration failed: {e}")
        return False

def run_scvi(adata, batch_key='library_id', counts_layer='counts'):
    try:
        import scvi # type: ignore[import]
        # Restore raw counts for scVI
        if counts_layer in adata.layers:
            adata.X = adata.layers[counts_layer].copy()
        else:
            print(f"Warning: counts layer '{counts_layer}' not found. Using current adata.X for scVI.")
        scvi.model.SCVI.setup_anndata(adata, batch_key=batch_key)
        scvi_model = scvi.model.SCVI(adata)
        scvi_model.train()
        adata.obsm["X_scvi"] = scvi_model.get_latent_representation()
        print(f"SCVI integration complete. Representation 'X_scvi' added to obsm.")
        return True
    except Exception as e:
        print(f"SCVI integration failed: {e}")
        return False

def run_fastmnn(adata, batch_key='library_id'):
    try:
        import scanpy.external as sce # type: ignore[import]
        batches = [adata[adata.obs[batch_key] == b].copy() for b in adata.obs[batch_key].unique()]
        mnn_corrected, *_ = sce.pp.mnn_correct(*batches, batch_key=batch_key)
        if hasattr(mnn_corrected, 'obsm') and 'X_pca' in mnn_corrected.obsm:
            adata.obsm['X_pca_fastmnn'] = mnn_corrected.obsm['X_pca']
            print(f"FastMNN integration complete. Representation 'X_pca_fastmnn' added to obsm.")
            return True
        else:
            print(f"FastMNN did not return expected 'X_pca' in obsm.")
            return False
    except Exception as e:
        print(f"FastMNN integration failed: {e}")
        return False


def run_scanorama(adata, batch_key='library_id', basis='X_pca', adjusted_basis='X_scanorama'):
    try:
        import scanpy.external as sce # type: ignore[import]
        # We must store cells by batch order for scanorama
        if batch_key not in adata.obs:
            raise ValueError(f"Batch key '{batch_key}' not found in adata.obs.")
        batch_idx_order = adata.obs[batch_key].astype('category').cat.codes
        adata.obs['batch_idx_order'] = batch_idx_order
        # Sort by batch_idx_order to ensure scanorama processes cells in the correct order
        adata = adata[adata.obs['batch_idx_order'].argsort()]
        sce.pp.scanorama_integrate(adata, batch_key, basis=basis, adjusted_basis=adjusted_basis)
        print(adata)
        # Weird patchh to allow scanorama to run within the function with a global object (it does not do in place mods)
        adata.obsm[adjusted_basis] = adata.obsm[adjusted_basis].copy()
        print(f"Scanorama integration complete. Representation '{adjusted_basis}' added to obsm.")
        return(adata)
    except Exception as e:
        print(f"Scanorama integration failed: {e}")
        return adata

def run_umap_and_plot(adata, obsm_key, method_name, study, base_dir, color=['cell_type']):
    try:
        sc.pp.neighbors(adata, use_rep=obsm_key)
        print('Neighbors computed for UMAP.')
        sc.tl.umap(adata)
        print(f"UMAP computed for {method_name} in {study}.")
        study_dir = study.replace('_embeddings', '')
        figures_dir = os.path.join(base_dir, study_dir, 'figures')
        if not os.path.exists(figures_dir):
            os.makedirs(figures_dir)
        save_path = os.path.join(figures_dir, f'{study_dir}_{method_name}.png')
        sc.pl.embedding(adata, basis='umap', color=color, save=f"_{study_dir}_{method_name}.png", show=False)
        umaps_gdrive = f'GDrive:BreastAtlas/{study_dir}/umaps/'
        local_umap_path = os.path.join('figures', f'umap_{study_dir}_{method_name}.png')
        subprocess.run(['rclone', 'copy', local_umap_path, umaps_gdrive])
        print(f"UMAP for {method_name} saved and uploaded.")
        return True
    except Exception as e:
        print(f"UMAP/plot/rclone failed for {method_name} in {study}: {e}")
        return False
    

# Load the whole anndata object for each study and merge embeddings
print("Loading and merging embeddings for all studies...")
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
        stored_obs_names = base_adata.obs_names.copy() #Store to assign to geneformer later
        base_adata = base_adata[idx, :]
    print(f"Set base AnnData for {study} with shape {base_adata.shape}")
    # Lowercase all obsm keys
    base_adata.obsm = {k.lower(): v for k, v in base_adata.obsm.items()}
    # Now add other embeddings, aligning by index
    for method in EMBEDDING_METHODS:
        if method == 'scgpt':
            continue  # already loaded
        method_file = f"{method}.h5ad"
        embedding_file = [f for f in embedding_files if method_file in f]
        file_path = embedding_file[0] if embedding_file else None
        if file_path:
            print(f"Loading embedding: {file_path}")
            if file_path.endswith('.h5ad'):
                adata_tmp = sc.read_h5ad(file_path)
            # This branch will never actually happen now
            elif file_path.endswith('.csv'):
                adata_tmp = sc.read_csv(file_path)
                adata_tmp = adata_tmp[1:, :]  # Assuming the first row is not a cell but the index of the CSV
                adata_tmp.obsm['X_geneformer'] = adata_tmp.X.copy()
                # Ensure that the number of observations between base_adata and adata_tmp match
                print(f"adata_tmp shape: {adata_tmp.shape}, base_adata shape: {base_adata.shape}")
                if adata_tmp.shape[0] ==  len(stored_obs_names):
                    adata_tmp.obs_names = stored_obs_names
                else:
                    print(f"adata_tmp shape {adata_tmp.shape} does not match base_adata shape {base_adata.shape}, skipping.")
                    del adata_tmp
                    continue

            else:
                print(f"Unsupported file format for {file_path}")
                continue
            adata_tmp.obsm = {k.lower(): v for k, v in adata_tmp.obsm.items()}
            # Align cells by index (assume same order, or use intersection if needed)
            # Here, we use the intersection of cell names (index)
            common_idx = base_adata.obs_names.intersection(adata_tmp.obs_names)
            # Reorder base adata to make sure we'll add the correct embeddings to the correct places
            base_adata = base_adata[common_idx, :]
            if len(common_idx) == 0:
                print(f"No overlapping cells between base AnnData and {file_path}, skipping.")
                del adata_tmp
                continue
            # scfoundation and geneformer have specific representation keys (scfoundation_embeddings and X_geneformer)
            if method == 'geneformer':
                rep_key = 'x_geneformer'
            elif method == 'scfoundation':
                rep_key = 'scfoundation_embeddings'
            else:
                rep_key = f"x_{method.lower()}"
            if rep_key in adata_tmp.obsm:
                arr = adata_tmp[common_idx, :].obsm[rep_key]
                # Only add if the intersection matches the base AnnData's obs_names exactly
                if list(common_idx) == list(base_adata.obs_names):
                    # Introducing a minor name change here for extra consistency.
                    if rep_key != 'scfoundation_embeddings':
                        rep_key_new = rep_key
                    else:
                        rep_key_new = 'x_scfoundation'
                    base_adata.obsm[rep_key_new] = arr
                    print(f"Added {rep_key} to base AnnData.obsm for {len(common_idx)} cells.")
                else:
                    print(f"Cell indices do not match exactly for {rep_key} in {file_path}. Skipping this embedding.")
            else:
                print(f"Representation {rep_key} not found in obsm for {file_path}")
            del adata_tmp
            gc.collect()
    print(f"Merged AnnData for {study} with representations: {list(base_adata.obsm.keys())}")
    print(f'Metadata columns in adata.obs: {base_adata.obs.columns.tolist()}')
    if LIBRARY_KEY not in base_adata.obs:
        print(f"Warning: {LIBRARY_KEY} not found in adata.obs, using potential searches.")
        potential_keys = ['biosample_id', 'orig.ident', 'Sample', 'SampleID', 'library_id','orig_ident','sample','patient_id','sampleid']
        for key in potential_keys:
            if key in base_adata.obs:
                # Check if there are more than one unique values
                if base_adata.obs[key].nunique() > 1:
                    print(f"Found {key} in adata.obs, using it as LIBRARY_KEY.")
                    base_adata.obs[LIBRARY_KEY] = base_adata.obs[key]
                    break
                else:
                    print(f"Found {key} in adata.obs, but it has only one unique value, skipping.")
                    continue
        else:
            print(f"Error: No valid LIBRARY_KEY found in adata.obs. Exiting.")
            continue
    if CELL_TYPE_KEY not in base_adata.obs:
        print(f"Warning: {CELL_TYPE_KEY} not found in adata.obs, using potential searches.")
        potential_keys = ['FinalCellType', 'subclass.l1', 'subclass.l2', 'Major.subtype', 'cell_type', 'celltype_major', 'ClusterName_AllCells', 'cluster', 'Cluster_Idents','cellType']
        for key in potential_keys:
            if key in base_adata.obs:
                # Check if there are more than one unique values
                if base_adata.obs[key].nunique() > 1:
                    print(f"Found {key} in adata.obs, using it as CELL_TYPE_KEY.")
                    base_adata.obs[CELL_TYPE_KEY] = base_adata.obs[key]
                    break
                else:
                    print(f"Found {key} in adata.obs, but it has only one unique value, skipping.")
                    continue
        else:
            print(f"Error: No valid CELL_TYPE_KEY found in adata.obs. Exiting.")
            continue
    # Ensure the library_id and cell_type columns are categorical
    base_adata.obs[LIBRARY_KEY] = base_adata.obs[LIBRARY_KEY].astype('str').astype('category')
    base_adata.obs[CELL_TYPE_KEY] = base_adata.obs[CELL_TYPE_KEY].astype('str').astype('category')
    print(f"Final AnnData for {study} has {base_adata.n_obs} cells and {base_adata.n_vars} genes.")
    print(f"Library key '{LIBRARY_KEY}' has {base_adata.obs[LIBRARY_KEY].nunique()} unique values.")
    print(f"Cell type key '{CELL_TYPE_KEY}' has {base_adata.obs[CELL_TYPE_KEY].nunique()} unique values.")

    # --- Batch Correction and Integration Methods ---
    print(f"Running PCA (unintegrated) for {study}")
    # Save raw counts before normalization/log1p for scVI
    sc.pp.filter_genes(base_adata, min_cells=10)
    sc.pp.filter_cells(base_adata, min_genes=300)
    sc.pp.filter_cells(base_adata, min_counts=600)
    base_adata.layers['counts'] = base_adata.X.copy()
    sc.pp.normalize_total(base_adata, target_sum=1e4)
    sc.pp.log1p(base_adata)
    base_adata.layers['lognorm'] = base_adata.X.copy()
    sc.pp.highly_variable_genes(base_adata, n_top_genes=2000, subset=True)
    sc.pp.scale(base_adata, max_value=10)
    sc.pp.pca(base_adata, n_comps=50, use_highly_variable=True)
    print(f"PCA complete. Representation 'X_pca' added to obsm.")

    # Define all integration keys and readable names
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
        ("X_scanorama", "scanorama"),
        ('X_pca_combat', 'combat'),
    ]
    # Remaking library key to remove non-existent categories
    base_adata.obs[LIBRARY_KEY] = base_adata.obs[LIBRARY_KEY].astype('str').astype('category')
    # Run batch correction methods using utility functions
    print(f"Running Combat integration for {study}")
    run_combat(base_adata, lognorm_layer='lognorm', batch_key=LIBRARY_KEY, n_comps=50)
    if ("X_pca_combat", "combat") not in integration_methods:
        integration_methods.append(("X_pca_combat", "combat"))
    print(f"Running Harmony integration for {study}")
    run_harmony(base_adata, batch_key=LIBRARY_KEY, basis='X_pca')

    print(f"Running SCVI integration for {study}")
    run_scvi(base_adata, batch_key=LIBRARY_KEY, counts_layer='counts')

    print(f"Not Running FastMNN integration for {study} as it does not exist for Python 3.12")
    #run_fastmnn(base_adata, batch_key=LIBRARY_KEY)

    print(f"Running Scanorama integration for {study}")
    base_adata = run_scanorama(base_adata, batch_key=LIBRARY_KEY, basis='X_pca', adjusted_basis='X_scanorama')

    for obsm_key, method_name in foundation_methods:
        if obsm_key in base_adata.obsm:
            integration_methods.append((obsm_key, method_name))
    # --- For each integration, compute UMAP, plot, and rclone upload ---
    available_obsms = []
    for obsm_key, method_name in integration_methods:
        if obsm_key in base_adata.obsm:
            available_obsms.append(obsm_key)
            print(f"Running UMAP for {method_name} integration for {study}")
            run_umap_and_plot(base_adata, obsm_key, method_name, study, BASE_DIR, color=['cell_type'])


    # --- Single scIB evaluation for all embeddings ---
    if RUN_SCIB and available_obsms:
        try:
            biocons = BioConservation(isolated_labels=False)
            bm = Benchmarker(
                base_adata,
                batch_key=LIBRARY_KEY,
                label_key=CELL_TYPE_KEY,
                embedding_obsm_keys=available_obsms,
                pre_integrated_embedding_obsm_key="X_pca",
                bio_conservation_metrics=biocons,
                batch_correction_metrics=BatchCorrection(),
                n_jobs=-1,
            )
            bm.prepare()
            bm.benchmark()
            bm.plot_results_table(show=False, min_max_scale=False,save_dir=os.path.join(BASE_DIR, study, 'figures'))
            os.system(f"rclone copy --progress {os.path.join(BASE_DIR, study, 'figures')} GDrive:BreastAtlas/{study}/figures/")
            print(f"scIB evaluation complete for all integrations in {study}.")
            del bm
            del biocons
            gc.collect()
        except Exception as e:
            print(f"scIB evaluation failed for all integrations in {study}: {e}")
        del available_obsms
        gc.collect()

    # --- Save the final AnnData object with all embeddings ---`
    final_save_path = os.path.join(BASE_DIR, study.replace('_embeddings',''), f"{study.replace('_embeddings','')}_final_embeddings.h5ad")
    base_adata.write(final_save_path)
    print(f"Final AnnData object saved for {study.replace('_embeddings','')} at {final_save_path}")
    del base_adata
    gc.collect()