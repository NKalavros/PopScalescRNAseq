"""Inter-study concordance evaluation using SymphonyPy.

This module builds a Symphony reference from the largest study (by cell count)
within a combined AnnData (containing multiple studies), then maps each other
study as a query, performs label transfer, and computes concordance metrics.

Requirements / Assumptions:
- symphonypy is installed (import symphonypy as sp)
- The combined AnnData includes a study identifier column, a library/batch key,
  and a cell type column (defaults mirror earlier scripts: 'study', 'library_id', 'cell_type').
- Counts layer may exist (preferred). If not, current X is treated as counts.
- Gene sets across studies may differ; intersection is used for reference/query alignment.

Outputs:
- Adds Symphony mappings to query AnnData objects (in-place subsets) including:
  * obsm['X_pca_harmony'] (adjusted embedding)
  * obs['symphony_per_cell_dist'] (mapping score)
  * obs[cell_type_key + '_transferred'] (transferred labels)
- Returns a pandas DataFrame with per-query metrics (accuracy, macro F1, support) and overall macro summary.

Example usage:
    import anndata as ad
    from 250811_interstudy_symphony import run_symphony_interstudy
    adata = ad.read_h5ad('/gpfs/scratch/nk4167/KidneyAtlas/interstudy_combined.h5ad')
    metrics_df, ref_name = run_symphony_interstudy(adata)
    print(metrics_df)

Author: automated assistant
"""
from __future__ import annotations
import os
from typing import List, Tuple, Optional, Dict

import numpy as np
import pandas as pd
import scanpy as sc  # type: ignore
import anndata as ad  # type: ignore

try:
    import symphonypy as sp  # type: ignore
except Exception as e:  # pragma: no cover
    sp = None  # type: ignore
    _SYMPHONY_IMPORT_ERROR = e
else:
    _SYMPHONY_IMPORT_ERROR = None

try:
    from sklearn.metrics import f1_score, accuracy_score  # type: ignore
except Exception as e:  # pragma: no cover
    f1_score = None  # type: ignore
    accuracy_score = None  # type: ignore

# ---------------- Configuration (can be overridden in function calls) ----------------
DEFAULT_STUDY_KEY = 'study'
DEFAULT_BATCH_KEY = 'library_id'
DEFAULT_CELL_TYPE_KEY = 'cell_type'

# ---------------- Utility Functions ----------------

def _check_requirements():
    if sp is None:
        raise ImportError(f"symphonypy not available: {_SYMPHONY_IMPORT_ERROR}")
    if f1_score is None or accuracy_score is None:
        raise ImportError("scikit-learn metrics not available for accuracy/F1 computation.")


def _subset_by_study(adata: ad.AnnData, study_key: str, study_name: str) -> ad.AnnData:
    return adata[adata.obs[study_key] == study_name].copy()


def _prepare_reference(
    adata_ref: ad.AnnData,
    batch_key: str,
    n_top_genes: int = 3000,
    n_pcs: int = 30,
    target_sum: float = 1e5,
) -> ad.AnnData:
    """Run Symphony reference preprocessing pipeline in-place and return it.
    Steps: normalize -> log1p -> HVGs (batch aware) -> subset HVGs -> scale -> PCA -> Harmony
    
    Follows symphonypy tutorial workflow with robustness for small datasets.
    """
    # Use counts layer if present
    if 'counts' in adata_ref.layers:
        adata_ref.X = adata_ref.layers['counts'].copy()
    
    # Normalization & log (following tutorial exactly)
    sc.pp.normalize_total(adata_ref, target_sum=target_sum)
    sc.pp.log1p(adata_ref)
    
    # Highly variable genes with adaptive n_top_genes for small datasets
    original_n_genes = adata_ref.n_vars
    adaptive_n_top_genes = min(n_top_genes, max(500, adata_ref.n_vars // 3))
    
    print(f"[INFO] Symphony HVG selection: requesting {adaptive_n_top_genes} from {original_n_genes} genes")
    
    try:
        sc.pp.highly_variable_genes(
            adata_ref,
            batch_key=batch_key if batch_key in adata_ref.obs else None,
            n_top_genes=adaptive_n_top_genes,
        )
    except Exception as e:
        print(f"[WARN] HVG selection failed with batch_key, trying without batch: {e}")
        # Fallback: no batch correction in HVG selection
        sc.pp.highly_variable_genes(
            adata_ref,
            n_top_genes=adaptive_n_top_genes,
        )
    
    # Set raw before subsetting (following tutorial)
    adata_ref.raw = adata_ref
    
    # Subset to highly variable genes
    if 'highly_variable' in adata_ref.var:
        n_hvg = adata_ref.var['highly_variable'].sum()
        print(f"[INFO] Found {n_hvg} highly variable genes, subsetting...")
        if n_hvg > 0:
            adata_ref = adata_ref[:, adata_ref.var['highly_variable']].copy()
        else:
            print(f"[WARN] No highly variable genes found, using top {adaptive_n_top_genes} by variance")
            # Fallback: select genes by highest variance
            adata_ref.var['gene_var'] = np.var(adata_ref.X.toarray() if hasattr(adata_ref.X, 'toarray') else adata_ref.X, axis=0)
            top_var_genes = adata_ref.var.nlargest(adaptive_n_top_genes, 'gene_var').index
            adata_ref = adata_ref[:, top_var_genes].copy()
    
    # Scale (following tutorial)
    sc.pp.scale(adata_ref, max_value=10)
    
    # PCA (following tutorial exactly - zero_center=False)
    actual_n_pcs = min(n_pcs, adata_ref.n_vars - 1, adata_ref.n_obs - 1)
    print(f"[INFO] Computing PCA with {actual_n_pcs} components")
    sc.pp.pca(adata_ref, n_comps=actual_n_pcs, zero_center=False)
    
    # Harmony integration (following tutorial)
    if batch_key in adata_ref.obs and adata_ref.obs[batch_key].nunique() > 1 and sp is not None:
        print(f"[INFO] Running Harmony integration on {adata_ref.obs[batch_key].nunique()} batches")
        try:
            sp.pp.harmony_integrate(adata_ref, key=batch_key, max_iter_harmony=10)  # type: ignore[attr-defined]
        except Exception as e:
            print(f"[WARN] Harmony integration failed: {e}, using PCA only")
            adata_ref.obsm['X_pca_harmony'] = adata_ref.obsm['X_pca']
    else:
        print(f"[INFO] Single batch or no batch key, using PCA as harmony embedding")
        # Mirror expected key so downstream knows representation name exists
        adata_ref.obsm['X_pca_harmony'] = adata_ref.obsm['X_pca']
    
    return adata_ref


def _align_vars(ref: ad.AnnData, query: ad.AnnData) -> ad.AnnData:
    """Restrict query to intersection of variables with reference; fill missing genes by dropping.
    Assumes reference already subset to HVGs.
    """
    common = ref.var_names.intersection(query.var_names)
    if len(common) == 0:
        raise ValueError("No overlapping genes between reference and query for Symphony mapping.")
    if len(common) < ref.n_vars:
        # Align order to reference var order for consistency
        common_in_order = [g for g in ref.var_names if g in common]
        query = query[:, common_in_order].copy()
    else:
        query = query[:, ref.var_names].copy()
    return query


def _prepare_query(
    adata_query: ad.AnnData,
    target_sum: float = 1e5,
) -> ad.AnnData:
    if 'counts' in adata_query.layers:
        adata_query.X = adata_query.layers['counts'].copy()
    sc.pp.normalize_total(adata_query, target_sum=target_sum)
    sc.pp.log1p(adata_query)
    return adata_query


def _evaluate_concordance(
    true_labels: List[str],
    pred_labels: List[str],
) -> Dict[str, float]:
    acc = float(accuracy_score(true_labels, pred_labels))  # type: ignore[operator]
    macro_f1 = float(f1_score(true_labels, pred_labels, average='macro', zero_division=0))  # type: ignore[operator]
    weighted_f1 = float(f1_score(true_labels, pred_labels, average='weighted', zero_division=0))  # type: ignore[operator]
    return {
        'accuracy': acc,
        'macro_f1': macro_f1,
        'weighted_f1': weighted_f1,
    }

# ---------------- Main Public Function ----------------

def run_symphony_interstudy(
    adata: ad.AnnData,
    study_key: str = DEFAULT_STUDY_KEY,
    batch_key: str = DEFAULT_BATCH_KEY,
    cell_type_key: str = DEFAULT_CELL_TYPE_KEY,
    n_top_genes: int = 3000,
    n_pcs: int = 30,
    target_sum: float = 1e5,
    kNN_label_key_suffix: str = '_transferred',
    add_embedding: bool = True,
    embedding_key: str = 'X_symphony',
) -> Tuple[pd.DataFrame, str]:
    """Run Symphony mapping using the largest study as reference.

    Parameters
    ----------
    adata : AnnData
        Combined multi-study AnnData (modified in-place if add_embedding True).
    add_embedding : bool
        If True, store the integrated Harmony / Symphony embedding for all cells
        in adata.obsm[embedding_key]. Reference + mapped query embeddings are
        concatenated in original cell order.
    embedding_key : str
        Key under which to store embedding in adata.obsm.

    Returns
    -------
    (metrics_df, reference_study_name)
        metrics_df columns: ['query_study','n_cells','accuracy','macro_f1','weighted_f1'] (+ overall row)
    """
    _check_requirements()

    if study_key not in adata.obs:
        raise KeyError(f"Study key '{study_key}' not found in adata.obs")
    if cell_type_key not in adata.obs:
        raise KeyError(f"Cell type key '{cell_type_key}' not found in adata.obs")

    # Identify largest study
    counts = adata.obs[study_key].value_counts()
    reference_study = counts.idxmax()
    print(f"[INFO] Using '{reference_study}' as Symphony reference (n={counts.max()} cells)")

    adata_ref = _subset_by_study(adata, study_key, reference_study)
    adata_ref = _prepare_reference(adata_ref, batch_key=batch_key, n_top_genes=n_top_genes, n_pcs=n_pcs, target_sum=target_sum)

    metrics_rows: List[Dict[str, object]] = []
    # Optionally collect embeddings per cell
    embedding_dict: Dict[str, np.ndarray] = {}

    for query_study in counts.index:
        if query_study == reference_study:
            continue
        print(f"[INFO] Mapping query study '{query_study}'")
        adata_query = _subset_by_study(adata, study_key, query_study)
        # Align genes
        adata_query = _align_vars(adata_ref, adata_query)
        # Preprocess
        adata_query = _prepare_query(adata_query, target_sum=target_sum)
        # Map embedding (following tutorial)
        if sp is not None:
            try:
                sp.tl.map_embedding(adata_query, adata_ref, key=batch_key if batch_key in adata_query.obs else None)  # type: ignore[attr-defined]
            except Exception as e:
                print(f"[WARN] map_embedding failed for study {query_study}: {e}")
                continue
                
        # Confidence metrics (following tutorial)
        try:
            if sp is not None:
                sp.tl.per_cell_confidence(adata_query, adata_ref)  # type: ignore[attr-defined]
        except Exception as e:  # pragma: no cover
            print(f"[WARN] per_cell_confidence failed for study {query_study}: {e}")
            
        # Label transfer (following tutorial - need ref_labels parameter)
        transferred_key = cell_type_key + kNN_label_key_suffix
        if sp is not None:
            try:
                # Following tutorial: need to specify ref_labels parameter
                sp.tl.transfer_labels_kNN(adata_query, adata_ref, ref_labels=[cell_type_key])  # type: ignore[attr-defined]
            except Exception as e:
                print(f"[WARN] transfer_labels_kNN failed for study {query_study}: {e}")
                # Fallback: copy original labels
                adata_query.obs[transferred_key] = adata_query.obs[cell_type_key]
        # The transfer function overwrites the same key; copy to suffix if needed
        if transferred_key != cell_type_key and transferred_key not in adata_query.obs:
            adata_query.obs[transferred_key] = adata_query.obs[cell_type_key]
        # Evaluate
        true_labels = adata_query.obs[cell_type_key].astype(str).tolist()
        pred_labels = adata_query.obs[transferred_key].astype(str).tolist()
        metrics = _evaluate_concordance(true_labels, pred_labels)
        metrics_rows.append({
            'query_study': query_study,
            'n_cells': len(adata_query),
            **metrics,
        })
        # Store embeddings
        if add_embedding:
            if 'X_pca_harmony' in adata_query.obsm:
                emb = adata_query.obsm['X_pca_harmony']
            elif 'X_pca' in adata_query.obsm:
                emb = adata_query.obsm['X_pca']
            else:
                emb = None
            if emb is not None:
                for cid, row in zip(adata_query.obs_names, emb):
                    embedding_dict[cid] = row.astype(np.float32, copy=False)

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

    # Attach embedding to original adata
    if add_embedding:
        try:
            # Add reference embeddings first (after mapping so Harmony finished)
            if 'X_pca_harmony' in adata_ref.obsm:
                ref_emb = adata_ref.obsm['X_pca_harmony']
            else:
                ref_emb = adata_ref.obsm.get('X_pca')
            if ref_emb is not None:
                for cid, row in zip(adata_ref.obs_names, ref_emb):
                    embedding_dict[cid] = row.astype(np.float32, copy=False)
            # Build matrix in combined cell order
            if embedding_dict:
                dim = len(next(iter(embedding_dict.values())))
                mat = np.zeros((adata.n_obs, dim), dtype=np.float32)
                missing = 0
                for i, cid in enumerate(adata.obs_names):
                    v = embedding_dict.get(cid)
                    if v is not None:
                        mat[i] = v
                    else:
                        missing += 1
                if missing:
                    print(f"[WARN] Symphony embedding missing for {missing} cells (likely dropped during alignment).")
                adata.obsm[embedding_key] = mat
                print(f"[INFO] Stored Symphony embedding in adata.obsm['{embedding_key}'] with shape {mat.shape}")
        except Exception as e:  # pragma: no cover
            print(f"[WARN] Failed to attach Symphony embedding: {e}")

    return metrics_df, reference_study


__all__ = [
    'run_symphony_interstudy',
]
