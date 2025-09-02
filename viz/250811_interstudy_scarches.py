"""Inter-study integration & label transfer using scArches (SCANVI surgery).

This module trains a SCANVI (scArches) reference model on the largest study
then maps (surgery) each remaining study as a query. It extracts a unified
latent embedding and performs KNN-based label transfer from the reference to
each query. Results are aggregated into metrics and the joint embedding is
stored back into the combined AnnData.

Key additions (in-place on provided combined AnnData):
  obsm['X_scarches']  -> unified latent embedding (all cells)
  obs['cell_type_scarches_pred']  -> predicted / transferred labels for queries
  obs['cell_type_scarches_uncert'] -> uncertainty scores (if available)

Returned metrics DataFrame columns:
  ['query_study','n_cells','accuracy','macro_f1','weighted_f1'] plus overall row.

References:
  HLCA mapping tutorial (scArches surgery):
    https://docs.scarches.org/en/latest/hlca_map_classify.html

Simplifications:
  - Uses largest study as reference (can be parameterized later).
  - Single label level (cell_type_key) instead of multi-level HLCA annotations.
  - Saves model to a temporary directory within BASE_DIR for query preparation.
"""
from __future__ import annotations
import os
import shutil
import tempfile
from typing import Dict, List, Tuple, Optional

import numpy as np
import pandas as pd
import anndata as ad  # type: ignore

DEFAULT_STUDY_KEY = 'study'
DEFAULT_BATCH_KEY = 'library_id'
DEFAULT_CELL_TYPE_KEY = 'cell_type'

try:  # scArches / scvi-tools
    import scarches as sca  # type: ignore
except Exception as e:  # pragma: no cover
    sca = None  # type: ignore
    _SCARCHES_IMPORT_ERROR = e
else:
    _SCARCHES_IMPORT_ERROR = None

try:
    from sklearn.metrics import f1_score, accuracy_score  # type: ignore
except Exception as e:  # pragma: no cover
    f1_score = None  # type: ignore
    accuracy_score = None  # type: ignore


def _check_requirements():
    if sca is None:
        raise ImportError(f"scArches (scarches) not available: {_SCARCHES_IMPORT_ERROR}")
    if f1_score is None or accuracy_score is None:
        raise ImportError("scikit-learn metrics not available for accuracy/F1 computation.")


def _align_genes_between_studies(adata: ad.AnnData, study_key: str) -> ad.AnnData:
    """Align genes across studies to common intersection for scArches compatibility.
    
    scArches requires all studies (reference and queries) to have identical gene sets.
    This function finds the intersection of genes across all studies and subsets 
    each study to this common gene set.
    
    Returns a copy of adata with aligned genes.
    """
    print("[DEBUG][scArches] Starting gene alignment across studies...")
    
    # Get all unique studies
    studies = adata.obs[study_key].unique()
    print(f"[DEBUG][scArches] Studies to align: {list(studies)}")
    
    # Find intersection of genes across all studies
    common_genes = None
    study_gene_counts = {}
    
    for study in studies:
        study_mask = adata.obs[study_key] == study
        study_adata = adata[study_mask]
        study_genes = set(study_adata.var_names)
        study_gene_counts[study] = len(study_genes)
        
        if common_genes is None:
            common_genes = study_genes
        else:
            common_genes = common_genes.intersection(study_genes)
    
    print(f"[DEBUG][scArches] Gene counts per study: {study_gene_counts}")
    
    if common_genes is None or len(common_genes) == 0:
        raise ValueError("No common genes found across studies")
    
    print(f"[DEBUG][scArches] Common genes across all studies: {len(common_genes)}")
    
    # Convert to sorted list for consistent ordering
    common_genes_list = sorted(list(common_genes))
    print(f"[DEBUG][scArches] Using {len(common_genes_list)} common genes for scArches")
    
    # Subset to common genes
    aligned_adata = adata[:, common_genes_list].copy()
    
    # Verify gene alignment
    for study in studies:
        study_mask = aligned_adata.obs[study_key] == study
        study_subset = aligned_adata[study_mask]
        print(f"[DEBUG][scArches] Study '{study}' after alignment: {study_subset.shape}")
    
    return aligned_adata


def _subset_by_study(adata: ad.AnnData, study_key: str, study_name: str) -> ad.AnnData:
    return adata[adata.obs[study_key] == study_name].copy()


def _evaluate(true_labels: List[str], pred_labels: List[str]) -> Dict[str, float]:
    return {
        'accuracy': float(accuracy_score(true_labels, pred_labels)),  # type: ignore[arg-type]
        'macro_f1': float(f1_score(true_labels, pred_labels, average='macro', zero_division=0)),  # type: ignore[arg-type]
        'weighted_f1': float(f1_score(true_labels, pred_labels, average='weighted', zero_division=0)),  # type: ignore[arg-type]
    }


def run_scarches_interstudy(
    adata: ad.AnnData,
    base_dir: Optional[str] = None,
    study_key: str = DEFAULT_STUDY_KEY,
    batch_key: Optional[str] = DEFAULT_BATCH_KEY,
    cell_type_key: str = DEFAULT_CELL_TYPE_KEY,
    embedding_key: str = 'X_scarches',
    pred_suffix: str = '_scarches_pred',
    uncert_suffix: str = '_scarches_uncert',
    latent_dim: int = 30,
    reference_max_epochs: int = 80,
    query_max_epochs: int = 30,
    unlabeled_category: str = 'unlabeled',
    early_stopping_patience: int = 15,
    n_neighbors_transfer: int = 50,
    seed: int = 0,
) -> Tuple[pd.DataFrame, str]:
    """Run scArches (SCANVI surgery) inter-study mapping.

    Parameters
    ----------
    adata : AnnData
        Combined multi-study AnnData (modified in-place).
    base_dir : Optional[str]
        Directory to persist temporary reference model; if None uses temp dir.
    latent_dim : int
        Latent dimensionality for SCANVI.
    reference_max_epochs / query_max_epochs : int
        Epoch counts for reference training & query surgery.
    unlabeled_category : str
        Category name used for unlabeled cells in SCANVI.
    n_neighbors_transfer : int
        K for weighted KNN label transfer from reference to query.

    Returns
    -------
    (metrics_df, reference_study_name)
    """
    _check_requirements()
    if study_key not in adata.obs:
        raise KeyError(f"Study key '{study_key}' not found in adata.obs")
    if cell_type_key not in adata.obs:
        raise KeyError(f"Cell type key '{cell_type_key}' not found in adata.obs")

    # Align genes across studies to prevent dimension mismatch errors
    print("[INFO][scArches] Aligning genes across studies...")
    adata_aligned = _align_genes_between_studies(adata, study_key)
    
    # CRITICAL: Convert sparse matrices to dense to prevent index corruption issues
    # scArches/scvi-tools can suffer from the same sparse matrix corruption as scPoli
    if hasattr(adata_aligned.X, 'toarray'):
        print(f"[DEBUG][scArches] Converting sparse matrix to dense to avoid indexing errors")
        adata_aligned.X = adata_aligned.X.toarray().astype(np.float32)
    else:
        print(f"[DEBUG][scArches] Working with dense matrix - ensuring contiguous memory layout")
        adata_aligned.X = np.ascontiguousarray(adata_aligned.X)
    
    # Choose reference (largest study)
    counts = adata_aligned.obs[study_key].value_counts()
    reference_study = counts.idxmax()
    print(f"[INFO][scArches] Using '{reference_study}' as SCANVI reference (n={counts.max()} cells)")
    ref_adata = _subset_by_study(adata_aligned, study_key, reference_study)
    
    # CRITICAL: Ensure reference data is dense and contiguous
    if hasattr(ref_adata.X, 'toarray'):
        print(f"[DEBUG][scArches] Converting reference sparse matrix to dense to avoid indexing errors")
        ref_adata.X = ref_adata.X.toarray().astype(np.float32)
    else:
        print(f"[DEBUG][scArches] Reference working with dense matrix - ensuring contiguous memory layout")
        ref_adata.X = np.ascontiguousarray(ref_adata.X)

    # Ensure unlabeled category exists
    if ref_adata.obs[cell_type_key].dtype.name == 'category':
        if unlabeled_category not in ref_adata.obs[cell_type_key].cat.categories:
            ref_adata.obs[cell_type_key] = ref_adata.obs[cell_type_key].cat.add_categories([unlabeled_category])
    else:
        ref_adata.obs[cell_type_key] = ref_adata.obs[cell_type_key].astype('category')

    print(f"[DEBUG][scArches] Reference study '{reference_study}' aligned shape: {ref_adata.shape}")
    print(f"[DEBUG][scArches] Unique cell types in reference: {len(ref_adata.obs[cell_type_key].unique())}")

    setup_kwargs = dict(labels_key=cell_type_key, unlabeled_category=unlabeled_category)
    if batch_key and batch_key in ref_adata.obs:
        setup_kwargs['batch_key'] = batch_key
    sca.models.SCANVI.setup_anndata(ref_adata, **setup_kwargs)  # type: ignore[attr-defined]

    print('[INFO][scArches] Initializing SCANVI reference model...')
    ref_model = sca.models.SCANVI(ref_adata, n_latent=latent_dim)  # type: ignore[attr-defined]

    early_stop_ref = dict(
        early_stopping_monitor='elbo_train',
        early_stopping_patience=early_stopping_patience,
        early_stopping_min_delta=0.001,
        plan_kwargs={'weight_decay': 0.0},
    )
    ref_model.train(max_epochs=reference_max_epochs, **early_stop_ref)

    # Save reference to disk
    if base_dir is None:
        tmp_dir = tempfile.mkdtemp(prefix='scarches_ref_')
        cleanup_dir = True
    else:
        tmp_dir = os.path.join(base_dir, 'scarches_reference_model')
        os.makedirs(tmp_dir, exist_ok=True)
        cleanup_dir = False
    ref_model.save(tmp_dir, overwrite=True)

    ref_latent = ref_model.get_latent_representation(ref_adata)
    ref_latent_map: Dict[str, np.ndarray] = {cid: row.astype(np.float32, copy=False) for cid, row in zip(ref_adata.obs_names, ref_latent)}

    metrics_rows: List[Dict[str, object]] = []
    query_latents: Dict[str, np.ndarray] = {}
    query_preds: Dict[str, str] = {}
    query_uncerts: Dict[str, float] = {}

    # KNN trainer on reference latent
    try:
        ref_emb_adata = ad.AnnData(ref_latent)
        ref_emb_adata.obs = ref_adata.obs[[cell_type_key]].copy()
        ref_emb_adata.X = ref_latent
        
        # CRITICAL: Ensure reference latent embedding is dense and contiguous for KNN operations
        if hasattr(ref_emb_adata.X, 'toarray'):
            print(f"[DEBUG][scArches] Converting reference latent sparse matrix to dense for KNN")
            ref_emb_adata.X = ref_emb_adata.X.toarray().astype(np.float32)
        else:
            print(f"[DEBUG][scArches] Reference latent working with dense matrix - ensuring contiguous memory layout")
            ref_emb_adata.X = np.ascontiguousarray(ref_emb_adata.X.astype(np.float32))
        
        # ULTRA-CRITICAL: Fix numerical stability issues that cause KNN index corruption
        print(f"[DEBUG][scArches] Fixing numerical stability issues in reference latent embedding")
        
        # 1. Check for and fix any NaN/inf values
        if np.any(~np.isfinite(ref_emb_adata.X)):
            print(f"[WARN][scArches] Found NaN/inf values in reference latent, replacing with zeros")
            ref_emb_adata.X = np.nan_to_num(ref_emb_adata.X, nan=0.0, posinf=0.0, neginf=0.0)
        
        # 2. First normalize to reasonable range to prevent numerical overflow
        print(f"[DEBUG][scArches] Normalizing reference embedding to prevent numerical overflow")
        max_abs_val = np.max(np.abs(ref_emb_adata.X))
        if max_abs_val > 1e-8:  # Only normalize if not already tiny
            ref_emb_adata.X = ref_emb_adata.X / max_abs_val
        
        # 3. Add noise to prevent zero standard deviations AFTER normalization
        print(f"[DEBUG][scArches] Adding noise to prevent zero standard deviations")
        noise_scale = 1e-6  # Larger noise scale to be more aggressive
        ref_emb_adata.X += np.random.RandomState(seed).normal(0, noise_scale, ref_emb_adata.X.shape).astype(np.float32)
        
        # 4. Check standard deviations and aggressively fix any that are still too small
        stds = np.std(ref_emb_adata.X, axis=0)
        problematic_cols = stds < 1e-6  # More aggressive threshold
        if np.any(problematic_cols):
            print(f"[WARN][scArches] Found {np.sum(problematic_cols)} columns with problematic std, fixing aggressively")
            # Replace problematic columns with random normal values
            ref_emb_adata.X[:, problematic_cols] = np.random.RandomState(seed+2).normal(0, 1e-4, (ref_emb_adata.X.shape[0], np.sum(problematic_cols))).astype(np.float32)
        
        # 5. Final check: ensure ALL standard deviations are reasonable
        final_stds = np.std(ref_emb_adata.X, axis=0)
        min_std = np.min(final_stds)
        print(f"[DEBUG][scArches] Final reference std check: min={min_std:.8f}, max={np.max(final_stds):.8f}")
        if min_std < 1e-8:
            print(f"[ERROR][scArches] Still have near-zero std after fixing! Adding global noise...")
            ref_emb_adata.X += np.random.RandomState(seed+3).normal(0, 1e-5, ref_emb_adata.X.shape).astype(np.float32)
        
        print(f"[DEBUG][scArches] Reference embedding shape: {ref_emb_adata.X.shape}, dtype: {ref_emb_adata.X.dtype}")
        print(f"[DEBUG][scArches] Reference embedding stats: min={np.min(ref_emb_adata.X):.6f}, max={np.max(ref_emb_adata.X):.6f}, std={np.std(ref_emb_adata.X):.6f}")
        
        knn_trainer = sca.utils.knn.weighted_knn_trainer(  # type: ignore[attr-defined]
            train_adata=ref_emb_adata,
            train_adata_emb='X',
            n_neighbors=n_neighbors_transfer,
        )
    except Exception as e:  # pragma: no cover
        print(f"[WARN][scArches] Failed to initialize weighted KNN trainer: {e}")
        knn_trainer = None

    for study_name in counts.index:
        if study_name == reference_study:
            continue
        query_adata = _subset_by_study(adata_aligned, study_key, study_name)
        print(f"[DEBUG][scArches] Processing query study '{study_name}' with shape {query_adata.shape}")
        
        # CRITICAL: Ensure query data is dense and contiguous to prevent index corruption
        if hasattr(query_adata.X, 'toarray'):
            print(f"[DEBUG][scArches] Converting query sparse matrix to dense to avoid indexing errors")
            query_adata.X = query_adata.X.toarray().astype(np.float32)
        else:
            print(f"[DEBUG][scArches] Query working with dense matrix - ensuring contiguous memory layout")
            query_adata.X = np.ascontiguousarray(query_adata.X)
        
        try:
            q_prep = sca.models.SCANVI.prepare_query_anndata(adata=query_adata, reference_model=tmp_dir, inplace=False)  # type: ignore[attr-defined]
        except Exception as e:
            print(f"[ERROR][scArches] prepare_query_anndata failed for study {study_name}: {e}")
            print(f"[ERROR][scArches] This indicates gene alignment issues between reference and query")
            raise RuntimeError(f"scArches prepare_query_anndata failed for study {study_name}: {e}")
        
        # CRITICAL: Ensure prepared query data is also dense and contiguous
        if hasattr(q_prep.X, 'toarray'):
            print(f"[DEBUG][scArches] Converting prepared query sparse matrix to dense to avoid indexing errors")
            q_prep.X = q_prep.X.toarray().astype(np.float32)
        else:
            print(f"[DEBUG][scArches] Prepared query working with dense matrix - ensuring contiguous memory layout")
            q_prep.X = np.ascontiguousarray(q_prep.X)
        if cell_type_key in q_prep.obs:
            q_prep.obs[cell_type_key] = unlabeled_category
        try:
            q_model = sca.models.SCANVI.load_query_data(q_prep, tmp_dir, freeze_dropout=True)  # type: ignore[attr-defined]
            q_model.train(max_epochs=query_max_epochs, **early_stop_ref)
            q_latent = q_model.get_latent_representation(q_prep)
        except Exception as e:
            print(f"[ERROR][scArches] Surgery failed for study {study_name}: {e}")
            print(f"[ERROR][scArches] This indicates gene alignment issues in scArches surgery")
            raise RuntimeError(f"scArches surgery failed for study {study_name}: {e}")
        try:
            if knn_trainer is not None:
                q_emb_adata = ad.AnnData(q_latent)
                q_emb_adata.X = q_latent
                
                # CRITICAL: Ensure latent embeddings are dense and contiguous for KNN operations
                if hasattr(q_emb_adata.X, 'toarray'):
                    print(f"[DEBUG][scArches] Converting query latent sparse matrix to dense for KNN")
                    q_emb_adata.X = q_emb_adata.X.toarray().astype(np.float32)
                else:
                    print(f"[DEBUG][scArches] Query latent working with dense matrix - ensuring contiguous memory layout")
                    q_emb_adata.X = np.ascontiguousarray(q_emb_adata.X.astype(np.float32))
                
                # ULTRA-CRITICAL: Fix numerical stability issues that cause KNN index corruption
                print(f"[DEBUG][scArches] Fixing numerical stability issues in query latent embedding")
                
                # 1. Check for and fix any NaN/inf values
                if np.any(~np.isfinite(q_emb_adata.X)):
                    print(f"[WARN][scArches] Found NaN/inf values in query latent, replacing with zeros")
                    q_emb_adata.X = np.nan_to_num(q_emb_adata.X, nan=0.0, posinf=0.0, neginf=0.0)
                
                # 2. First normalize to reasonable range to prevent numerical overflow
                print(f"[DEBUG][scArches] Normalizing query embedding to prevent numerical overflow")
                max_abs_val = np.max(np.abs(q_emb_adata.X))
                if max_abs_val > 1e-8:  # Only normalize if not already tiny
                    q_emb_adata.X = q_emb_adata.X / max_abs_val
                
                # 3. Add noise to prevent zero standard deviations AFTER normalization
                print(f"[DEBUG][scArches] Adding noise to prevent zero standard deviations")
                noise_scale = 1e-6  # Larger noise scale to be more aggressive
                # Ensure seed is within valid range [0, 2**32-1]
                query_seed = abs(hash(study_name)) % (2**32 - 1)
                q_emb_adata.X += np.random.RandomState(query_seed).normal(0, noise_scale, q_emb_adata.X.shape).astype(np.float32)
                
                # 4. Check standard deviations and aggressively fix any that are still too small
                stds = np.std(q_emb_adata.X, axis=0)
                problematic_cols = stds < 1e-6  # More aggressive threshold
                if np.any(problematic_cols):
                    print(f"[WARN][scArches] Found {np.sum(problematic_cols)} columns with problematic std, fixing aggressively")
                    # Replace problematic columns with random normal values
                    q_emb_adata.X[:, problematic_cols] = np.random.RandomState(query_seed+2).normal(0, 1e-4, (q_emb_adata.X.shape[0], np.sum(problematic_cols))).astype(np.float32)
                
                # 5. Final check: ensure ALL standard deviations are reasonable
                final_stds = np.std(q_emb_adata.X, axis=0)
                min_std = np.min(final_stds)
                print(f"[DEBUG][scArches] Final query std check: min={min_std:.8f}, max={np.max(final_stds):.8f}")
                if min_std < 1e-8:
                    print(f"[ERROR][scArches] Still have near-zero std after fixing! Adding global noise...")
                    q_emb_adata.X += np.random.RandomState(query_seed+3).normal(0, 1e-5, q_emb_adata.X.shape).astype(np.float32)
                
                print(f"[DEBUG][scArches] Query embedding shape: {q_emb_adata.X.shape}, dtype: {q_emb_adata.X.dtype}")
                print(f"[DEBUG][scArches] Query embedding stats: min={np.min(q_emb_adata.X):.6f}, max={np.max(q_emb_adata.X):.6f}, std={np.std(q_emb_adata.X):.6f}")
                
                # Debug information about shapes before KNN transfer
                print(f"[DEBUG][scArches] Reference embedding shape: {ref_emb_adata.shape}")
                print(f"[DEBUG][scArches] Query embedding shape: {q_emb_adata.shape}")
                print(f"[DEBUG][scArches] Query latent shape: {q_latent.shape}")
                print(f"[DEBUG][scArches] Reference has {len(ref_emb_adata.obs[cell_type_key].unique())} cell types")
                
                labels_df, uncert_df = sca.utils.knn.weighted_knn_transfer(  # type: ignore[attr-defined]
                    query_adata=q_emb_adata,
                    query_adata_emb='X',
                    label_keys=cell_type_key,
                    knn_model=knn_trainer,
                    ref_adata_obs=ref_emb_adata.obs,
                )
                preds_list = labels_df[cell_type_key].astype(str).tolist()
                uncert_list = uncert_df[cell_type_key].astype(float).tolist()
            else:
                raise RuntimeError('KNN trainer unavailable - gene alignment may have failed')
        except Exception as e:
            print(f"[ERROR][scArches] Label transfer failed for {study_name}: {e}")
            print(f"[ERROR][scArches] This indicates gene alignment issues between reference and query")
            raise RuntimeError(f"scArches failed for study {study_name} due to gene alignment: {e}")
        
        true_labels = query_adata.obs[cell_type_key].astype(str).tolist()
        metrics = _evaluate(true_labels, preds_list)
        metrics_rows.append({
            'query_study': study_name,
            'n_cells': q_prep.n_obs,
            **metrics,
        })
        for cid, row, pl, ul in zip(q_prep.obs_names, q_latent, preds_list, uncert_list):
            query_latents[cid] = row.astype(np.float32, copy=False)
            query_preds[cid] = pl
            query_uncerts[cid] = float(ul)

    metrics_df = pd.DataFrame(metrics_rows).sort_values('query_study')
    if not metrics_df.empty:
        overall = {
            'query_study': '__OVERALL_MACRO__',
            'n_cells': metrics_df['n_cells'].sum(),
            'accuracy': (metrics_df['accuracy'] * metrics_df['n_cells']).sum() / metrics_df['n_cells'].sum(),
            'macro_f1': metrics_df['macro_f1'].mean(),
            'weighted_f1': (metrics_df['weighted_f1'] * metrics_df['n_cells']).sum() / metrics_df['n_cells'].sum(),
        }
        metrics_df = pd.concat([metrics_df, pd.DataFrame([overall])], ignore_index=True)

    try:
        dim = next(iter(ref_latent_map.values())).shape[0]
        mat = np.zeros((adata.n_obs, dim), dtype=np.float32)
        preds_series: List[Optional[str]] = []
        uncert_series: List[Optional[float]] = []
        missing = 0
        for idx, cid in enumerate(adata.obs_names):
            if cid in ref_latent_map:
                mat[idx] = ref_latent_map[cid]
                preds_series.append(None)
                uncert_series.append(None)
            elif cid in query_latents:
                mat[idx] = query_latents[cid]
                preds_series.append(query_preds.get(cid))
                uncert_series.append(query_uncerts.get(cid))
            else:
                missing += 1
                preds_series.append(None)
                uncert_series.append(None)
        if missing:
            print(f"[WARN][scArches] Missing scArches latent for {missing} cells.")
        adata.obsm[embedding_key] = mat
        adata.obs[f"{cell_type_key}{pred_suffix}"] = preds_series
        adata.obs[f"{cell_type_key}{uncert_suffix}"] = uncert_series
        print(f"[INFO][scArches] Stored scArches latent embedding in obsm['{embedding_key}'] shape {mat.shape}")
    except Exception as e:  # pragma: no cover
        print(f"[WARN][scArches] Failed to attach scArches embedding: {e}")

    if 'cleanup_dir' in locals() and cleanup_dir:
        try:
            shutil.rmtree(tmp_dir, ignore_errors=True)
        except Exception:
            pass

    return metrics_df, reference_study

__all__ = [
    'run_scarches_interstudy',
]
