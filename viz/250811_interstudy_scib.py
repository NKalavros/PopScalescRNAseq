"""scIB benchmarking for multiple embeddings (foundation + integration methods).

This module evaluates a set of embedding representations stored in
`adata.obsm` using scIB metrics for one or more batch keys. It compares each
embedding (treated as the integrated data) against the original unintegrated
counts (the provided `adata`).

Approach:
  For each embedding key K in embedding_keys present in adata.obsm:
    1. Construct a lightweight AnnData `adata_int` whose X matrix is the
       embedding (dense) and obs copied from original.
    2. Run scib.metrics.metrics with specified batch_key and label_key.
    3. Collect metric outputs into a long-form DataFrame:
         columns: embedding_key, batch_key, metric, value.

Notes:
  - Some scIB metrics need gene-level info or HVG overlap; since embeddings do
    not retain gene expression we pass only the embedding as X; metrics that
    rely on expression-based features may be NA or fail and are skipped.
  - For efficiency we disable very slow metrics by default (kBET, LISI) unless
    `full=True` is passed.

Configuration parameters allow selection of fast / slim / all metrics.
"""
from __future__ import annotations
from typing import List, Dict, Optional
import numpy as np
import pandas as pd
import anndata as ad  # type: ignore

try:
    import scib  # type: ignore
    from scib.metrics import metrics as scib_metrics  # type: ignore
    import scanpy as sc  # type: ignore
except Exception as e:  # pragma: no cover
    scib = None  # type: ignore
    scib_metrics = None  # type: ignore
    sc = None  # type: ignore
    _SCIB_IMPORT_ERROR = e
else:
    _SCIB_IMPORT_ERROR = None


def _check_requirements():
    if scib is None or scib_metrics is None or sc is None:
        raise ImportError(f"scIB or scanpy not available: {_SCIB_IMPORT_ERROR}")


def run_scib_benchmark(
    adata: ad.AnnData,
    embedding_keys: List[str],
    batch_keys: List[str],
    label_key: str,
    fast: bool = True,
    slim: bool = False,
    full: bool = False,
    min_cells_per_label: int = 5,
) -> pd.DataFrame:
    """Run scIB metrics for each embedding and batch key.

    Parameters
    ----------
    adata : AnnData
        Original (unintegrated) AnnData with expression matrix.
    embedding_keys : list[str]
        obsm keys containing embedding matrices to benchmark.
    batch_keys : list[str]
        Batch columns to evaluate (each separately).
    label_key : str
        Cell type / biological label column in obs.
    fast / slim / full : bool
        Select scope of scIB metrics. Only one should typically be True; priority full>slim>fast.

    Returns
    -------
    DataFrame long-form: embedding_key, batch_key, metric, value
    """
    _check_requirements()
    results_rows: List[Dict[str, object]] = []

    if full:
        preset = 'all'
    elif slim:
        preset = 'slim'
    else:
        preset = 'fast'

    METRIC_FLAGS_BASE = dict(
        ari_=True,
        nmi_=True,
        silhouette_=True,  # silhouette label and batch ASW
        graph_conn_=False,
        pcr_=True,
        ilisi_=full,  # only in full
        clisi_=full,  # only in full
        kBET_=False,  # Always disable kBET (very slow)
        isolated_labels_f1_=True,  # Disable isolated labels (very slow clustering optimization)
        isolated_labels_asw_=False,  # Disable isolated labels (very slow clustering optimization)
        hvg_score_=False,
        cell_cycle_=False,
        trajectory_=False,
    )

    if label_key in adata.obs:
        vc = adata.obs[label_key].value_counts()
        rare = set(vc[vc < min_cells_per_label].index)
        if rare:
            tmp_label_key = f"{label_key}__scibtmp"
            adata.obs[tmp_label_key] = adata.obs[label_key].astype(str).where(~adata.obs[label_key].isin(rare), 'Rare')
            label_to_use = tmp_label_key
        else:
            label_to_use = label_key
    else:
        raise KeyError(f"Label key '{label_key}' not in adata.obs")

    for emb_key in embedding_keys:
        if emb_key not in adata.obsm:
            continue
        emb = adata.obsm[emb_key]
        if not isinstance(emb, np.ndarray):
            try:
                emb = emb.A
            except Exception:
                emb = np.asarray(emb)
        adata_int = ad.AnnData(X=emb.copy(), obs=adata.obs.copy())
        
        # Compute PCA and kNN graph on the embedding for scIB metrics that require them
        # Some metrics require X_pca to be available
        # For embeddings, we compute PCA on the embedding itself
        try:
            sc.tl.pca(adata_int, n_comps=min(50, emb.shape[1] - 1, emb.shape[0] - 1))
            # Also compute neighbors as some metrics may need it
            sc.pp.neighbors(adata_int, use_rep='X_pca')
        except Exception:
            # If preprocessing fails, try without PCA (will limit available metrics)
            try:
                sc.pp.neighbors(adata_int, use_rep='X')
            except Exception:
                pass  # Some metrics may still work without neighbors
            
        for bkey in batch_keys:
            if bkey not in adata.obs:
                continue
            try:
                scores = scib_metrics(
                    adata,
                    adata_int,
                    batch_key=bkey,
                    label_key=label_to_use,
                    **METRIC_FLAGS_BASE,
                )
                for col, val in scores.iloc[0].items():
                    results_rows.append({
                        'embedding_key': emb_key,
                        'batch_key': bkey,
                        'metric': col,
                        'value': float(val) if isinstance(val, (int, float, np.floating)) else val,
                        'preset': preset,
                    })
            except Exception as e:  # pragma: no cover
                results_rows.append({
                    'embedding_key': emb_key,
                    'batch_key': bkey,
                    'metric': 'ERROR',
                    'value': str(e),
                    'preset': preset,
                })

    return pd.DataFrame(results_rows)

__all__ = [
    'run_scib_benchmark',
]
