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

def _reset_categorical_variables(adata: ad.AnnData, keys: List[str]) -> None:
    """Reset categorical variables to prevent index corruption issues."""
    for key in keys:
        if key in adata.obs:
            # Convert to string, remove any NaN, then back to clean categorical
            series = adata.obs[key].astype(str)
            series = series.replace(['nan', 'None', 'null'], 'Unknown')
            adata.obs[key] = pd.Categorical(series)
            print(f"[DEBUG][scPoli] Reset categorical '{key}': {len(adata.obs[key].cat.categories)} categories")


def _align_genes_between_studies(adata: ad.AnnData, study_key: str) -> ad.AnnData:
    """Align genes across studies to common intersection for scPoli compatibility."""
    print("[DEBUG][scPoli] Starting gene alignment across studies...")
    studies = adata.obs[study_key].unique()
    common_genes = None
    for st in studies:
        sub = adata[adata.obs[study_key] == st]
        genes = set(sub.var_names)
        common_genes = genes if common_genes is None else common_genes.intersection(genes)
    if not common_genes:
        raise ValueError("No common genes found across studies for scPoli")
    common_list = sorted(common_genes)
    print(f"[DEBUG][scPoli] Using {len(common_list)} common genes for alignment")
    
    # CRITICAL: Create copy and preserve categorical variables during gene subsetting
    aligned_adata = adata[:, common_list].copy()
    
    # Ensure all categorical variables are preserved correctly
    for col in adata.obs.columns:
        if isinstance(adata.obs[col].dtype, pd.CategoricalDtype):
            print(f"[DEBUG][scPoli] Preserving categorical variable '{col}' during gene alignment")
            aligned_adata.obs[col] = adata.obs[col].copy()
    
    return aligned_adata


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
    # Align genes across all studies to avoid dimension mismatches
    print("[INFO][scPoli] Aligning genes across studies before reference building...")
    adata = _align_genes_between_studies(adata, study_key)
    if study_key not in adata.obs:
        raise KeyError(f"Study key '{study_key}' not found in adata.obs")
    if cell_type_key not in adata.obs:
        raise KeyError(f"Cell type key '{cell_type_key}' not found in adata.obs")

    counts = adata.obs[study_key].value_counts()
    
    # Check cell type coverage across studies
    print("[INFO][scPoli] Analyzing cell type coverage across studies...")
    all_cell_types = set(adata.obs[cell_type_key].unique())
    study_cell_types = {}
    
    for study in counts.index:
        study_adata = _subset_by_study(adata, study_key, study)
        study_cell_types[study] = set(study_adata.obs[cell_type_key].unique())
    
    # Find the study with the most comprehensive cell type coverage
    coverage_scores = {}
    for study in counts.index:
        coverage = len(study_cell_types[study])
        cell_count = counts[study]
        # Score combines cell type diversity and total cell count
        coverage_scores[study] = (coverage, cell_count)
    
    # Sort by cell type coverage first, then by cell count
    best_study = max(coverage_scores.keys(), key=lambda x: coverage_scores[x])
    reference_study = best_study
    
    ref_cell_types = study_cell_types[reference_study]
    missing_types = all_cell_types - ref_cell_types
    
    if missing_types:
        print(f"[INFO][scPoli] Reference study '{reference_study}' missing {len(missing_types)} cell types: {missing_types}")
        print("[WARN][scPoli] Skipping reference augmentation to avoid categorical corruption issues")
        print("[WARN][scPoli] This may reduce label transfer accuracy for missing cell types")
        
        # Skip augmentation to avoid categorical corruption that causes the indexing error
        ref_adata = _subset_by_study(adata, study_key, reference_study)
        
        # # Original augmentation code (disabled for debugging):
        # print("[INFO][scPoli] Will augment reference with samples from other studies to ensure full coverage")
        # 
        # # Create augmented reference by adding cells from other studies for missing types
        # ref_adata = _subset_by_study(adata, study_key, reference_study)
        # 
        # # Add representative cells for missing cell types from other studies
        # augment_cells = []
        # for missing_type in missing_types:
        #     # Find studies that have this cell type
        #     donor_studies = [s for s in counts.index if missing_type in study_cell_types[s] and s != reference_study]
        #     if donor_studies:
        #         # Take from the study with most cells of this type
        #         best_donor = max(donor_studies, key=lambda x: counts[x])
        #         donor_adata = _subset_by_study(adata, study_key, best_donor)
        #         type_cells = donor_adata[donor_adata.obs[cell_type_key] == missing_type]
        #         
        #         # Take a reasonable sample (up to 50 cells per missing type)
        #         n_sample = min(50, len(type_cells))
        #         if n_sample > 0:
        #             sampled = type_cells[:n_sample].copy()
        #             augment_cells.append(sampled)
        #             print(f"[INFO][scPoli] Adding {n_sample} '{missing_type}' cells from study '{best_donor}'")
        # 
        # # Combine reference with augmented cells
        # if augment_cells:
        #     # Debug: Check categorical variables before concatenation
        #     print(f"[DEBUG][scPoli] Before concat - ref study categories: {ref_adata.obs[study_key].unique()}")
        #     print(f"[DEBUG][scPoli] Before concat - ref library categories: {ref_adata.obs.get(batch_key, pd.Series()).unique()}")
        #     
        #     for i, aug_cell in enumerate(augment_cells):
        #         print(f"[DEBUG][scPoli] Augment {i} - study cats: {aug_cell.obs[study_key].unique()}")
        #         print(f"[DEBUG][scPoli] Augment {i} - library cats: {aug_cell.obs.get(batch_key, pd.Series()).unique()}")
        #     
        #     # Use concatenate with batch_key to avoid categorical corruption
        #     augmented_ref = ad.concat([ref_adata] + augment_cells, join='outer', merge='unique')
        #     ref_adata = augmented_ref
        #     
        #     # Debug: Check categorical variables after concatenation
        #     print(f"[DEBUG][scPoli] After concat - study categories: {ref_adata.obs[study_key].unique()}")
        #     print(f"[DEBUG][scPoli] After concat - library categories: {ref_adata.obs.get(batch_key, pd.Series()).unique()}")
        #     
        #     # Ensure categorical variables are properly encoded as strings then back to categories
        #     ref_adata.obs[study_key] = ref_adata.obs[study_key].astype(str).astype('category')
        #     if batch_key and batch_key in ref_adata.obs:
        #         ref_adata.obs[batch_key] = ref_adata.obs[batch_key].astype(str).astype('category')
        #     ref_adata.obs[cell_type_key] = ref_adata.obs[cell_type_key].astype(str).astype('category')
        #     
        #     print(f"[INFO][scPoli] Augmented reference from {counts[reference_study]} to {ref_adata.n_obs} cells with full cell type coverage")
    else:
        ref_adata = _subset_by_study(adata, study_key, reference_study)
        print(f"[INFO][scPoli] Reference study '{reference_study}' has complete cell type coverage")
    
    # Debug: show final reference cell types
    ref_cell_types_final = set(ref_adata.obs[cell_type_key].unique())
    print(f"[DEBUG][scPoli] Reference cell types after augmentation: {len(ref_cell_types_final)}")
    print(f"[DEBUG][scPoli] Reference types: {sorted(ref_cell_types_final)}")
    
    print(f"[INFO][scPoli] Using '{reference_study}' as reference base (n={counts[reference_study]} cells)")
    
    # Debug: show cell types being used by scPoli
    all_cell_types_in_data = set(adata.obs[cell_type_key].unique())
    print(f"[DEBUG][scPoli] Total unique cell types in dataset: {len(all_cell_types_in_data)}")
    print(f"[DEBUG][scPoli] Cell types: {sorted(all_cell_types_in_data)}")
    
    # Find common genes across ALL datasets BEFORE any filtering to ensure consistency
    print("[INFO][scPoli] Finding common genes across all studies before any filtering...")
    all_studies = list(counts.index)
    all_gene_sets = []
    
    # Collect genes from each study without any filtering first
    for study in all_studies:
        study_adata = _subset_by_study(adata, study_key, study)
        all_gene_sets.append(set(study_adata.var_names))
    
    # Find intersection of all gene sets - genes present in ALL studies
    common_genes = set.intersection(*all_gene_sets)
    common_genes = sorted(list(common_genes))  # Sort for consistent ordering
    
    if len(common_genes) == 0:
        raise ValueError("[ERROR][scPoli] No genes common to all studies")
    
    print(f"[INFO][scPoli] Found {len(common_genes)} genes common to all studies")
    
    # Now subset reference to common genes BEFORE any filtering
    ref_adata = ref_adata[:, common_genes].copy()
    
    # Apply filtering only AFTER gene alignment to maintain consistency
    print("[INFO][scPoli] Applying gene filtering to reference after alignment...")
    sc.pp.filter_genes(ref_adata, min_cells=3)
    
    # Get the final filtered gene set from reference
    final_genes = sorted(list(ref_adata.var_names))
    print(f"[INFO][scPoli] After filtering, reference has {len(final_genes)} genes")
    
    if 0 < reference_fraction < 1.0:
        ref_indices = np.random.default_rng(seed).choice(ref_adata.obs_names, size=int(len(ref_adata)*reference_fraction), replace=False)
        ref_adata = ref_adata[ref_indices].copy()
        print(f"[INFO][scPoli] Subsampled reference to {ref_adata.n_obs} cells (fraction={reference_fraction})")
    
    # Ensure reference data is float32 to avoid dtype mismatches
    if hasattr(ref_adata.X, 'dtype') and ref_adata.X.dtype != np.float32:
        ref_adata.X = ref_adata.X.astype(np.float32)
    
    # CRITICAL: Convert sparse matrices to dense if needed to avoid indexing issues
    # The axis 1 index error often occurs with sparse matrix corruption
    if hasattr(ref_adata.X, 'toarray'):
        print(f"[DEBUG][scPoli] Converting sparse matrix to dense to avoid indexing errors")
        ref_adata.X = ref_adata.X.toarray().astype(np.float32)
    
    # Clean cell type labels to remove any problematic characters/formats
    ref_adata.obs[cell_type_key] = ref_adata.obs[cell_type_key].astype(str)
    
    # Setup condition keys - use study as the main condition
    condition_keys = [study_key]
    if batch_key and batch_key in ref_adata.obs:
        condition_keys.append(batch_key)

    # CRITICAL: Reset all categorical variables to prevent corruption
    print("[DEBUG][scPoli] Resetting categorical variables to prevent index corruption...")
    _reset_categorical_variables(ref_adata, condition_keys + [cell_type_key])

    # CRITICAL: Debug categorical variables before scPoli initialization
    print(f"[DEBUG][scPoli] FINAL CHECK before model init:")
    print(f"[DEBUG][scPoli] Reference dimensions: {ref_adata.shape}")
    print(f"[DEBUG][scPoli] Gene set length: {len(final_genes)}")
    print(f"[DEBUG][scPoli] Condition keys: {condition_keys}")
    
    for key in condition_keys + [cell_type_key]:
        if key in ref_adata.obs:
            cats = ref_adata.obs[key].unique()
            print(f"[DEBUG][scPoli] Key '{key}': {len(cats)} categories = {sorted(cats)}")
            # Check for any NaN or problematic values
            nan_count = ref_adata.obs[key].isnull().sum()
            if nan_count > 0:
                print(f"[ERROR][scPoli] Key '{key}' has {nan_count} NaN values!")
                # Fill NaN values with string representation
                ref_adata.obs[key] = ref_adata.obs[key].fillna('Unknown').astype(str).astype('category')
        else:
            print(f"[ERROR][scPoli] Key '{key}' missing from obs!")
    
    # Additional check for X data integrity
    print(f"[DEBUG][scPoli] X matrix info: shape={ref_adata.X.shape}, dtype={ref_adata.X.dtype}")
    if hasattr(ref_adata.X, 'nnz'):
        print(f"[DEBUG][scPoli] X is sparse, nnz={ref_adata.X.nnz}")
    
    print(f"[DEBUG][scPoli] Reference dimensions before model init: {ref_adata.shape}")
    print(f"[DEBUG][scPoli] Reference gene set length: {len(final_genes)}")
    print(f"[DEBUG][scPoli] Reference condition keys: {condition_keys}")

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
        
        # Use the exact same final gene set as the reference (after filtering)
        # This ensures perfect alignment between reference and query
        missing_genes = set(final_genes) - set(query_adata.var_names)
        if missing_genes:
            raise ValueError(f"[ERROR][scPoli] Query study '{study_name}' missing genes after filtering: {missing_genes}")
        
        query_adata = query_adata[:, final_genes].copy()
        
        print(f"[INFO][scPoli] Query '{study_name}' aligned to {len(final_genes)} genes (same as reference)")
        print(f"[DEBUG][scPoli] Query dimensions: {query_adata.shape}")
        print(f"[DEBUG][scPoli] Query gene names match reference: {list(query_adata.var_names) == list(ref_adata.var_names)}")
        
        # Ensure data is float32 to avoid dtype mismatches
        if hasattr(query_adata.X, 'dtype') and query_adata.X.dtype != np.float32:
            query_adata.X = query_adata.X.astype(np.float32)
        
        # CRITICAL: Convert sparse matrices to dense if needed to avoid indexing issues
        if hasattr(query_adata.X, 'toarray'):
            print(f"[DEBUG][scPoli] Converting query sparse matrix to dense to avoid indexing errors")
            query_adata.X = query_adata.X.toarray().astype(np.float32)
        
        # Clean query cell type labels
        query_adata.obs[cell_type_key] = query_adata.obs[cell_type_key].astype(str)
        
        # CRITICAL: Ensure query categorical variables match reference encoding
        print(f"[DEBUG][scPoli] Query '{study_name}' categorical debugging:")
        _reset_categorical_variables(query_adata, condition_keys + [cell_type_key])
        
        for key in condition_keys + [cell_type_key]:
            if key in query_adata.obs:
                cats = query_adata.obs[key].unique()
                print(f"[DEBUG][scPoli] Query key '{key}': {len(cats)} categories = {sorted(cats)}")
                # Check for NaN values
                nan_count = query_adata.obs[key].isnull().sum()
                if nan_count > 0:
                    print(f"[ERROR][scPoli] Query key '{key}' has {nan_count} NaN values!")
                    query_adata.obs[key] = query_adata.obs[key].fillna('Unknown').astype(str).astype('category')
                # Ensure consistency with reference
                query_adata.obs[key] = query_adata.obs[key].astype(str).astype('category')
            else:
                print(f"[ERROR][scPoli] Query key '{key}' missing from obs!")
        
        # Verify gene order consistency
        if not np.array_equal(query_adata.var_names, ref_adata.var_names):
            print(f"[ERROR][scPoli] Gene order mismatch between reference and query '{study_name}'")
            print(f"Reference genes[:5]: {list(ref_adata.var_names[:5])}")
            print(f"Query genes[:5]: {list(query_adata.var_names[:5])}")
            raise ValueError(f"Gene order mismatch for query '{study_name}'")
        
        print(f"[INFO][scPoli] Loading query data for study '{study_name}' (n={query_adata.n_obs}, genes={query_adata.n_vars})")
        
        try:
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
            
        except Exception as e:
            print(f"[WARN][scPoli] Failed to process query {study_name}: {e}")
            continue

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
