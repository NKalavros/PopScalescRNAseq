"""Inter-study integration and label transfer using scPoli (scArches).

This module mirrors the Symphony workflow but uses scPoli to build a reference
on the largest study (by cell count) and incrementally map the remaining
studies as queries. For each query we fine-tune via surgery and perform label
transfer. A unified latent embedding is stored in the original combined AnnData.

Key outputs stored in-place (where possible):
- obsm['X_scpoli'] : consolidated latent embedding (cells x latent_dim)
- obs['cell_type_scpoli_pred'] : transferred / predicted cell types (per cell)
- obs['cell_type_scpoli_uncert'] : uncertainty scores (if available)

Returns a metrics DataFrame similar to Symphony: per-query accuracy, macro F1,
weighted F1 plus an overall row.

NOTE: This is a simplified pipeline intended for benchmarking. Hyperparameters
are minimal and may need tuning for production use. Training epochs are kept
modest for speed.
"""
from __future__ import annotations
import os
from typing import Dict, List, Tuple, Optional

import numpy as np
import pandas as pd
import anndata as ad  # type: ignore
import scanpy as sc  # type: ignore

DEFAULT_STUDY_KEY = 'study'
DEFAULT_BATCH_KEY = 'library_id'  # (optional, scPoli condition)
DEFAULT_CELL_TYPE_KEY = 'cell_type'

try:
    from scarches.models.scpoli import scPoli  # type: ignore
except Exception as e:  # pragma: no cover
    scPoli = None  # type: ignore
    _SCPOLI_IMPORT_ERROR = e
else:
    _SCPOLI_IMPORT_ERROR = None

try:
    from sklearn.metrics import f1_score, accuracy_score  # type: ignore
except Exception as e:  # pragma: no cover
    f1_score = None  # type: ignore
    accuracy_score = None  # type: ignore


def _check_requirements():
    if scPoli is None:
        raise ImportError(f"scPoli not available (scarches.models.scpoli import failed): {_SCPOLI_IMPORT_ERROR}")
    if f1_score is None or accuracy_score is None:
        raise ImportError("scikit-learn metrics not available for accuracy/F1 computation.")


def _subset_by_study(adata: ad.AnnData, study_key: str, study_name: str) -> ad.AnnData:
    return adata[adata.obs[study_key] == study_name].copy()


def _evaluate(true_labels: List[str], pred_labels: List[str]) -> Dict[str, float]:
    return {
        'accuracy': float(accuracy_score(true_labels, pred_labels)),  # type: ignore[call-arg]
        'macro_f1': float(f1_score(true_labels, pred_labels, average='macro', zero_division=0)),  # type: ignore[call-arg]
        'weighted_f1': float(f1_score(true_labels, pred_labels, average='weighted', zero_division=0)),  # type: ignore[call-arg]
    }


def run_scpoli_interstudy(
    adata: ad.AnnData,
    study_key: str = DEFAULT_STUDY_KEY,
    batch_key: Optional[str] = DEFAULT_BATCH_KEY,
    cell_type_key: str = DEFAULT_CELL_TYPE_KEY,
    embedding_key: str = 'X_scpoli',
    pred_suffix: str = '_scpoli_pred',
    uncert_suffix: str = '_scpoli_uncert',
    reference_fraction: float = 1.0,
    latent_dim: int = 30,
    pretraining_epochs: int = 40,
    total_epochs: int = 50,
    eta: int = 5,
    max_queries: Optional[int] = None,
    seed: int = 0,
) -> Tuple[pd.DataFrame, str]:
    if scPoli is None:
        raise ImportError(f"scPoli not available: {_SCPOLI_IMPORT_ERROR}")
    _check_requirements()
    if study_key not in adata.obs:
        raise KeyError(f"Study key '{study_key}' not found in adata.obs")
    if cell_type_key not in adata.obs:
        raise KeyError(f"Cell type key '{cell_type_key}' not found in adata.obs")

    counts = adata.obs[study_key].value_counts()
    reference_study = counts.idxmax()
    print(f"[INFO][scPoli] Using '{reference_study}' as reference (n={counts.max()} cells)")

    ref_adata = _subset_by_study(adata, study_key, reference_study)
    if 0 < reference_fraction < 1.0:
        ref_indices = np.random.default_rng(seed).choice(ref_adata.obs_names, size=int(len(ref_adata)*reference_fraction), replace=False)
        ref_adata = ref_adata[ref_indices].copy()
        print(f"[INFO][scPoli] Subsampled reference to {ref_adata.n_obs} cells (fraction={reference_fraction})")

    # Basic filtering to remove genes with all zeros
    sc.pp.filter_genes(ref_adata, min_counts=1)

    condition_key = study_key  # treat study as condition for integration
    condition_keys = condition_key

    print("[INFO][scPoli] Initializing reference scPoli model...")
    model = scPoli(
        adata=ref_adata,
        condition_keys=condition_keys,
        cell_type_keys=cell_type_key,
        embedding_dims=latent_dim,
        recon_loss='nb',
    )

    early_stopping_kwargs = {
        'early_stopping_metric': 'val_prototype_loss',
        'mode': 'min',
        'threshold': 0,
        'patience': 10,
        'reduce_lr': True,
        'lr_patience': 7,
        'lr_factor': 0.2,
    }

    model.train(
        n_epochs=total_epochs,
        pretraining_epochs=pretraining_epochs,
        early_stopping_kwargs=early_stopping_kwargs,
        eta=eta,
    )

    # Collect embeddings + predictions for reference cells
    model.model.eval()
    ref_latent = model.get_latent(ref_adata, mean=True)
    ref_pred = [None] * ref_adata.n_obs  # predictions only for queries; keep placeholder
    ref_uncert = [None] * ref_adata.n_obs

    metrics_rows: List[Dict[str, object]] = []

    # Iterate query studies
    processed_queries = 0
    query_latents: Dict[str, np.ndarray] = {}
    query_preds: Dict[str, str] = {}
    query_uncerts: Dict[str, float] = {}

    for study_name in counts.index:
        if study_name == reference_study:
            continue
        if max_queries is not None and processed_queries >= max_queries:
            break
        query_adata = _subset_by_study(adata, study_key, study_name)
        # Align genes intersection
        common = ref_adata.var_names.intersection(query_adata.var_names)
        if len(common) == 0:
            print(f"[WARN][scPoli] No overlapping genes with query {study_name}; skipping.")
            continue
        query_adata = query_adata[:, common].copy()
        # Some preprocessing to ensure gene order matches reference (subset)
        ref_sub = ref_adata[:, common].copy()

        print(f"[INFO][scPoli] Loading query data for study '{study_name}' (n={query_adata.n_obs})")
        query_model = scPoli.load_query_data(
            adata=query_adata,
            reference_model=model,
            labeled_indices=[],  # treat as unlabeled for adaptation
        )
        query_model.train(
            n_epochs=total_epochs,
            pretraining_epochs=pretraining_epochs,
            eta=eta,
        )
        results = query_model.classify(query_adata, scale_uncertainties=True)
        preds = results[cell_type_key]['preds']
        uncert = results[cell_type_key]['uncert']
        true_labels = query_adata.obs[cell_type_key].astype(str).tolist()
        pred_labels = preds.tolist()
        metrics = _evaluate(true_labels, pred_labels)
        metrics_rows.append({
            'query_study': study_name,
            'n_cells': len(query_adata),
            **metrics,
        })

        # Latent
        query_latent = query_model.get_latent(query_adata, mean=True)
        # Store per-cell
        for cid, row, pl, ul in zip(query_adata.obs_names, query_latent, pred_labels, uncert.tolist()):
            query_latents[cid] = row.astype(np.float32, copy=False)
            query_preds[cid] = pl
            query_uncerts[cid] = float(ul)

        processed_queries += 1

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

    # Assemble full embedding matrix in combined order
    try:
        dim = ref_latent.shape[1]
        full_mat = np.zeros((adata.n_obs, dim), dtype=np.float32)
        # Reference cells
        ref_map = {cid: row.astype(np.float32, copy=False) for cid, row in zip(ref_adata.obs_names, ref_latent)}
        missing = 0
        preds_series = []
        uncert_series = []
        for cid in adata.obs_names:
            if cid in ref_map:
                full_mat[adata.obs_names.get_loc(cid)] = ref_map[cid]
                preds_series.append(None)
                uncert_series.append(None)
            elif cid in query_latents:
                full_mat[adata.obs_names.get_loc(cid)] = query_latents[cid]
                preds_series.append(query_preds.get(cid))
                uncert_series.append(query_uncerts.get(cid))
            else:
                missing += 1
                preds_series.append(None)
                uncert_series.append(None)
        if missing:
            print(f"[WARN][scPoli] Missing latent for {missing} cells (likely not processed)")
        adata.obsm[embedding_key] = full_mat
        adata.obs[f"{cell_type_key}{pred_suffix}"] = preds_series
        adata.obs[f"{cell_type_key}{uncert_suffix}"] = uncert_series
        print(f"[INFO][scPoli] Stored scPoli latent embedding in obsm['{embedding_key}'] shape {full_mat.shape}")
    except Exception as e:  # pragma: no cover
        print(f"[WARN][scPoli] Failed to attach scPoli embedding: {e}")

    return metrics_df, reference_study

__all__ = [
    'run_scpoli_interstudy',
]
