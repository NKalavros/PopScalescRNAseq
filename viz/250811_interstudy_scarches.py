"""Inter-study integration & label transfer using scArches (SCANVI surgery).

This module follows the official scArches SCANVI surgery pipeline to train a 
reference model on the largest study, then perform surgery on each remaining 
study as queries. It uses the built-in SCANVI classifier for predictions
instead of custom KNN transfer to avoid numerical stability issues.

Key additions (in-place on provided combined AnnData):
  obsm['X_scarches']  -> unified latent embedding (all cells)
  obs['cell_type_scarches_pred']  -> predicted labels from SCANVI classifier

Returned metrics DataFrame columns:
  ['query_study','n_cells','accuracy','macro_f1','weighted_f1'] plus overall row.

References:
  Official scArches SCANVI surgery pipeline:
    https://docs.scarches.org/en/latest/scanvi_surgery_pipeline.html

Implementation follows the official tutorial exactly to avoid numerical issues.
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
    from scarches.dataset.trvae.data_handling import remove_sparsity  # type: ignore
except Exception as e:  # pragma: no cover
    sca = None  # type: ignore
    remove_sparsity = None  # type: ignore
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
    if remove_sparsity is None:
        raise ImportError(f"remove_sparsity function not available: {_SCARCHES_IMPORT_ERROR}")
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
    latent_dim: int = 30,
    reference_max_epochs: int = 80,
    query_max_epochs: int = 30,
    unlabeled_category: str = 'Unknown',
    early_stopping_patience: int = 15,
    seed: int = 0,
) -> Tuple[pd.DataFrame, str]:
    """Run scArches (SCANVI surgery) inter-study mapping following official pipeline.

    This implementation follows the official scArches SCANVI surgery pipeline exactly:
    https://docs.scarches.org/en/latest/scanvi_surgery_pipeline.html

    Parameters
    ----------
    adata : AnnData
        Combined multi-study AnnData (modified in-place).
    base_dir : Optional[str]
        Directory to persist reference model; if None uses temp dir.
    latent_dim : int
        Latent dimensionality for SCANVI.
    reference_max_epochs / query_max_epochs : int
        Epoch counts for reference training & query surgery.
    unlabeled_category : str
        Category name used for unlabeled cells in SCANVI.

    Returns
    -------
    (metrics_df, reference_study_name)
    """
    _check_requirements()
    if study_key not in adata.obs:
        raise KeyError(f"Study key '{study_key}' not found in adata.obs")
    if cell_type_key not in adata.obs:
        raise KeyError(f"Cell type key '{cell_type_key}' not found in adata.obs")

    # Align genes across studies and ensure proper data format
    print("[INFO][scArches] Aligning genes across studies...")
    adata_aligned = _align_genes_between_studies(adata, study_key)
    
    # CRITICAL: Following official pipeline - remove sparsity and ensure count data in X
    print("[INFO][scArches] Removing sparsity and ensuring count data format...")
    assert remove_sparsity is not None, "remove_sparsity function not available"
    adata_aligned = remove_sparsity(adata_aligned)
    
    # Choose reference (largest study) 
    counts = adata_aligned.obs[study_key].value_counts()
    reference_study = counts.idxmax()
    print(f"[INFO][scArches] Using '{reference_study}' as SCANVI reference (n={counts.max()} cells)")
    
    # Split into reference and target datasets following official pipeline
    source_adata = adata_aligned[adata_aligned.obs[study_key] == reference_study].copy()
    print(f"[DEBUG][scArches] Reference dataset shape: {source_adata.shape}")
    print(f"[DEBUG][scArches] Unique cell types in reference: {len(source_adata.obs[cell_type_key].unique())}")

    # Setup reference model following official pipeline
    setup_kwargs = dict(labels_key=cell_type_key, unlabeled_category=unlabeled_category)
    if batch_key and batch_key in source_adata.obs:
        setup_kwargs['batch_key'] = batch_key
    
    print("[INFO][scArches] Setting up SCANVI reference model...")
    assert sca is not None, "scArches not available"
    sca.models.SCANVI.setup_anndata(source_adata, **setup_kwargs)

    # Create and train SCANVI model following official pipeline
    print("[INFO][scArches] Creating SCANVI reference model...")
    assert sca is not None, "scArches not available"
    scanvae = sca.models.SCANVI(
        source_adata,
        n_latent=latent_dim,
        unlabeled_category=unlabeled_category
    )
    
    print("[INFO][scArches] Training SCANVI reference model...")
    scanvae.train(
        max_epochs=reference_max_epochs,
        early_stopping_patience=early_stopping_patience,
        early_stopping_monitor='elbo_train',
        early_stopping_min_delta=0.001,
        plan_kwargs={'weight_decay': 0.0},
    )

    # Save reference model
    if base_dir is None:
        tmp_dir = tempfile.mkdtemp(prefix='scarches_ref_')
        cleanup_dir = True
    else:
        tmp_dir = os.path.join(base_dir, 'scarches_reference_model')
        os.makedirs(tmp_dir, exist_ok=True)
        cleanup_dir = False
    
    print(f"[INFO][scArches] Saving reference model to {tmp_dir}")
    scanvae.save(tmp_dir, overwrite=True)

    # Get reference latent representation
    ref_latent = scanvae.get_latent_representation()
    ref_latent_map: Dict[str, np.ndarray] = {
        cid: row.astype(np.float32, copy=False) 
        for cid, row in zip(source_adata.obs_names, ref_latent)
    }

    metrics_rows: List[Dict[str, object]] = []
    query_latents: Dict[str, np.ndarray] = {}
    query_preds: Dict[str, str] = {}

    # Process each query study following official surgery pipeline
    for study_name in counts.index:
        if study_name == reference_study:
            continue
            
        print(f"[INFO][scArches] Processing query study '{study_name}'...")
        target_adata = adata_aligned[adata_aligned.obs[study_key] == study_name].copy()
        print(f"[DEBUG][scArches] Query dataset shape: {target_adata.shape}")
        
        # Following official pipeline: set all query labels to unlabeled category
        # This is critical for proper SCANVI surgery behavior
        target_adata.obs['orig_cell_types'] = target_adata.obs[cell_type_key].copy()
        target_adata.obs[cell_type_key] = unlabeled_category
        
        try:
            # Load query data using official surgery pipeline
            print(f"[DEBUG][scArches] Loading query data for surgery...")
            assert sca is not None, "scArches not available"
            model = sca.models.SCANVI.load_query_data(
                target_adata,
                tmp_dir,
                freeze_dropout=True,
            )
            
            # Set indices for unsupervised surgery following official pipeline
            model._unlabeled_indices = np.arange(target_adata.n_obs)
            model._labeled_indices = []
            print(f"[DEBUG][scArches] Unlabeled indices: {len(model._unlabeled_indices)}")
            
            # Train query model
            print(f"[DEBUG][scArches] Training query model for surgery...")
            model.train(
                max_epochs=query_max_epochs,
                plan_kwargs={'weight_decay': 0.0},
                early_stopping_patience=early_stopping_patience,
                early_stopping_monitor='elbo_train',
                early_stopping_min_delta=0.001,
            )
            
            # Get latent representation and predictions using official method
            q_latent = model.get_latent_representation()
            predictions = model.predict()
            
        except Exception as e:
            print(f"[ERROR][scArches] Surgery failed for study {study_name}: {e}")
            raise RuntimeError(f"scArches surgery failed for study {study_name}: {e}")
        
        # Evaluate predictions against original true labels
        true_labels = target_adata.obs['orig_cell_types'].astype(str).tolist()
        pred_labels = predictions.astype(str).tolist()
        
        metrics = _evaluate(true_labels, pred_labels)
        metrics_rows.append({
            'query_study': study_name,
            'n_cells': target_adata.n_obs,
            **metrics,
        })
        
        # Store results
        for cid, row, pred in zip(target_adata.obs_names, q_latent, predictions):
            query_latents[cid] = row.astype(np.float32, copy=False)
            query_preds[cid] = str(pred)

    # Aggregate metrics
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

    # Store results in original adata
    try:
        dim = next(iter(ref_latent_map.values())).shape[0]
        mat = np.zeros((adata.n_obs, dim), dtype=np.float32)
        preds_series: List[Optional[str]] = []
        missing = 0
        
        for idx, cid in enumerate(adata.obs_names):
            if cid in ref_latent_map:
                mat[idx] = ref_latent_map[cid]
                preds_series.append(None)  # Reference has no predictions
            elif cid in query_latents:
                mat[idx] = query_latents[cid]
                preds_series.append(query_preds.get(cid))
            else:
                missing += 1
                preds_series.append(None)
                
        if missing:
            print(f"[WARN][scArches] Missing scArches latent for {missing} cells.")
            
        adata.obsm[embedding_key] = mat
        adata.obs[f"{cell_type_key}{pred_suffix}"] = preds_series
        print(f"[INFO][scArches] Stored scArches latent embedding in obsm['{embedding_key}'] shape {mat.shape}")
        
    except Exception as e:  # pragma: no cover
        print(f"[WARN][scArches] Failed to attach scArches embedding: {e}")

    # Cleanup
    if 'cleanup_dir' in locals() and cleanup_dir:
        try:
            shutil.rmtree(tmp_dir, ignore_errors=True)
        except Exception:
            pass

    return metrics_df, reference_study

__all__ = [
    'run_scarches_interstudy',
]
