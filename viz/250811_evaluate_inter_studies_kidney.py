"""
Inter-study evaluation script.

Goals (initial version):
1. Load final embedding AnnData files from multiple study directories.
2. Concatenate them into a single AnnData object for inter-study benchmarking.
3. Add two batch-related metadata columns:
   - Existing library-level batch key (LIBRARY_KEY) already present in per-study objects.
   - New Study key derived from directory / study name.
4. Provide an optional LLM-powered helper to standardize (harmonize) cell type labels across studies.
   (Not executed by default; requires OPENAI_API_KEY and network access.)
5. Save the combined AnnData plus (optionally) a JSON mapping of harmonized labels.

DEBUG MODE:
- Set DEBUG_MODE=True to subsample each study to DEBUG_N_CELLS cells for faster testing
- Debug mode automatically adds '_debug' suffix to all output files
- Uses random seed 42 for reproducible subsampling

Future (not yet implemented):
- Run cross-study batch correction metrics (scIB) using both batch keys (study vs library) in separate benchmark runs.
- Allow selecting subset of embeddings / dimensionality reductions for evaluation.

This script is intentionally modular and side-effect light so it can be imported.
"""
from __future__ import annotations
import os
import json
import time
from typing import Dict, List, Optional, Tuple

import scanpy as sc  # type: ignore
import anndata as ad  # type: ignore
import pandas as pd  # type: ignore
import numpy as np  # type: ignore

# ------------------------------- Configuration ---------------------------------
BASE_DIR = '/gpfs/scratch/nk4167/KidneyAtlas'  # Adjust as needed
STUDIES: List[str] = [
    'lake_scrna', 'lake_snrna', 'Abedini', 'SCP1288', 'Krishna', 'Braun'
]
FINAL_SUFFIX = '_final_embeddings.h5ad'
LIBRARY_KEY = 'library_id'
CELL_TYPE_KEY = 'cell_type'
STUDY_KEY = 'study'  # New key we add
DEBUG_MODE = True  # Set True to use only subset of cells for faster testing
DEBUG_N_CELLS = 1000  # Number of cells to subsample per study in debug mode
RUN_LLM_STANDARDIZATION = True  # Set True to attempt LLM-based label harmonization
LLM_MODEL = 'gpt-4o'  # Change if desired
LLM_TEMPERATURE = 0.0
LLM_MAX_RETRIES = 3
LLM_SLEEP_BETWEEN_RETRIES = 5
LLM_USE_STRUCTURED = True  # Attempt JSON schema structured output (newer SDKs)
OUTPUT_COMBINED_FILENAME = 'interstudy_combined_debug.h5ad' if DEBUG_MODE else 'interstudy_combined.h5ad'
LABEL_MAPPING_JSON = 'interstudy_celltype_mapping_debug.json' if DEBUG_MODE else 'interstudy_celltype_mapping.json'
MAX_STANDARD_LABELS = 30  # Hard cap for standardized labels
RUN_SYMPHONY = True  # Run Symphony inter-study mapping & label transfer evaluation
SYMPHONY_MODULE_FILENAME = '250811_interstudy_symphony.py'  # File containing run_symphony_interstudy
SYMPHONY_METRICS_CSV = 'symphony_interstudy_metrics_debug.csv' if DEBUG_MODE else 'symphony_interstudy_metrics.csv'
RUN_SCPOLI = True  # Run scPoli inter-study integration & label transfer
SCPOLI_MODULE_FILENAME = '250811_interstudy_scpoli.py'
SCPOLI_METRICS_CSV = 'scpoli_interstudy_metrics_debug.csv' if DEBUG_MODE else 'scpoli_interstudy_metrics.csv'
RUN_SCARCHES =True  # Run scArches (SCANVI surgery) inter-study integration
SCARCHES_MODULE_FILENAME = '250811_interstudy_scarches.py'
SCARCHES_METRICS_CSV = 'scarches_interstudy_metrics_debug.csv' if DEBUG_MODE else 'scarches_interstudy_metrics.csv'
RUN_SCIB_BENCHMARKER = True  # Run scIB benchmarker evaluation like intra-studies script
RUN_SCIB = True  # Run scIB benchmarking across embeddings
SCIB_MODULE_FILENAME = '250811_interstudy_scib.py'
SCIB_METRICS_CSV = 'scib_metrics_long_debug.csv' if DEBUG_MODE else 'scib_metrics_long.csv'
SCIB_EMBEDDINGS: List[str] = [
    'X_symphony', 'X_scpoli', 'X_scarches',  # integration methods
    # Potential foundation model embeddings (add if present):
    'x_geneformer_helical', 'x_scgpt_helical', 'x_uce', 'x_scfoundation', 'x_transcriptformer','x_scimilarity'
]
SCIB_BATCH_KEYS: List[str] = [STUDY_KEY, LIBRARY_KEY]

# ------------------------------- Utilities -------------------------------------

def clean_cell_type_labels(adata: ad.AnnData, cell_type_key: str = CELL_TYPE_KEY) -> None:
    """Clean cell type labels by removing problematic suffixes like .1, .2, etc.
    
    This helps with label harmonization across studies that may have added unique
    suffixes to prevent duplicate labels during processing.
    """
    if cell_type_key not in adata.obs:
        print(f"[WARN] '{cell_type_key}' not found in adata.obs; skipping label cleaning.")
        return
    
    import re
    original_labels = adata.obs[cell_type_key].astype(str)
    
    # Remove suffixes like .1, .2, .10, etc. but preserve meaningful ones like CD4+, CD8+, etc.
    cleaned_labels = original_labels.str.replace(r'\.(\d+)$', '', regex=True)
    
    # Count how many labels were changed
    n_changed = (original_labels != cleaned_labels).sum()
    if n_changed > 0:
        print(f"[INFO] Cleaned {n_changed} cell type labels (removed numeric suffixes)")
        adata.obs[cell_type_key] = cleaned_labels.astype('category')
    else:
        print(f"[INFO] No cell type labels needed cleaning")


def find_final_file(study: str, base_dir: str = BASE_DIR, suffix: str = FINAL_SUFFIX) -> Optional[str]:
    """Return path to the expected final embeddings file if it exists."""
    candidate = os.path.join(base_dir, study, f"{study}{suffix}")
    return candidate if os.path.isfile(candidate) else None


def load_studies(studies: List[str]) -> List[ad.AnnData]:
    """Load available studies' final embeddings AnnData objects, adding study key.

    Skips studies whose final file is missing. Warns if required metadata keys are missing.
    If DEBUG_MODE is True, subsamples each study to DEBUG_N_CELLS cells.
    """
    adatas: List[ad.AnnData] = []
    for s in studies:
        path = find_final_file(s)
        if not path:
            print(f"[WARN] Final embeddings file not found for study '{s}', skipping.")
            continue
        print(f"[INFO] Loading {path}")
        a = sc.read_h5ad(path)
        
        # Debug mode: subsample cells
        if DEBUG_MODE and a.n_obs > DEBUG_N_CELLS:
            print(f"[DEBUG] Subsampling study '{s}' from {a.n_obs} to {DEBUG_N_CELLS} cells")
            # Use scanpy's random subsampling to maintain diversity
            sc.pp.subsample(a, n_obs=DEBUG_N_CELLS, random_state=42)
        elif DEBUG_MODE:
            print(f"[DEBUG] Study '{s}' has {a.n_obs} cells (≤ {DEBUG_N_CELLS}), keeping all")
        
        # add study key
        a.obs[STUDY_KEY] = s
        
        # Clean cell type labels (remove .1, .2, etc. suffixes)
        clean_cell_type_labels(a, CELL_TYPE_KEY)
        
        # Ensure categorical types
        for key in [LIBRARY_KEY, CELL_TYPE_KEY, STUDY_KEY]:
            if key in a.obs:
                a.obs[key] = a.obs[key].astype(str).astype('category')
            else:
                print(f"[WARN] Key '{key}' missing in study '{s}'.")
        adatas.append(a)
    return adatas


def concatenate(adatas: List[ad.AnnData]) -> ad.AnnData:
    """Concatenate multiple AnnData objects aligning variables (outer join)."""
    if not adatas:
        raise ValueError("No AnnData objects provided for concatenation.")
    print(f"[INFO] Concatenating {len(adatas)} studies...")
    combined = ad.concat(adatas, join='outer', label=STUDY_KEY, keys=[a.obs[STUDY_KEY].unique()[0] for a in adatas])
    # Reinforce categorical types
    for key in [LIBRARY_KEY, CELL_TYPE_KEY, STUDY_KEY]:
        if key in combined.obs:
            combined.obs[key] = combined.obs[key].astype(str).astype('category')
    print(f"[INFO] Combined shape: {combined.shape}")
    return combined

# -------------------- LLM-based label harmonization helper ---------------------
LLM_SYSTEM_PROMPT = (
    "You are a domain expert harmonizing single-cell RNA-seq cell type labels across datasets. "
    "Given a list of raw labels, produce a concise JSON object mapping each ORIGINAL label to a STANDARD label. "
    "Rules: (1) Preserve biologically meaningful granularity (e.g., 'CD4+ T cell' not just 'T cell' unless overly specific). "
    "(2) Use consistent capitalization (Title Case except common markers like CD4+). "
    "(3) If two labels are clear synonyms, map both to a single standard form. "
    "(4) Return ONLY valid JSON (no markdown) with keys as original labels and values as standardized labels. "
    f"(5) Limit the TOTAL number of UNIQUE standardized labels to at most {MAX_STANDARD_LABELS} by merging related or rare originals into broader canonical categories (e.g., merge granular T cell activation states into 'CD4+ T Cell' or 'CD8+ T Cell'). "
    f"(6) Prefer biologically meaningful parent types over 'Other'; only use 'Other' if absolutely necessary to stay within {MAX_STANDARD_LABELS}."
)


def _call_openai(
    messages: List[Dict[str, str]],
    model: str = LLM_MODEL,
    temperature: float = LLM_TEMPERATURE,
    response_format: Optional[dict] = None,
) -> str:
    """Internal: call OpenAI Chat Completions API (current stable) with fallbacks.

    We first attempt the modern client.chat.completions.create pathway. If unavailable, try legacy
    openai.ChatCompletion.create. Expects 'messages' list of dicts.
    """
    api_key = os.getenv('OPENAI_API_KEY')
    if not api_key:
        raise RuntimeError("OPENAI_API_KEY not set in environment; cannot call LLM.")
    try:
        # New style
        from openai import OpenAI  # type: ignore
        client = OpenAI(api_key=api_key)
        try:
            kwargs = dict(model=model, messages=messages, temperature=temperature)
            if response_format is not None:
                # Newer SDKs: response_format expects a dict schema; add defensively.
                # type: ignore[arg-type]
                kwargs["response_format"] = response_format  # type: ignore
            resp = client.chat.completions.create(**kwargs)
            # Standard shape: resp.choices[0].message.content
            if hasattr(resp, 'choices') and resp.choices:
                choice0 = resp.choices[0]
                # Different SDK versions may store content under message or directly
                content = None
                if hasattr(choice0, 'message') and getattr(choice0.message, 'content', None):
                    content = choice0.message.content
                elif hasattr(choice0, 'text'):
                    content = choice0.text
                if not content:
                    content = str(resp)
                return content
        except AttributeError:
            # Fallback to legacy import path
            pass
    except Exception:
        # Legacy global namespace fallback
        pass
    # Legacy fallback: openai.ChatCompletion
    try:  # pragma: no cover - network/SDK dependent
        import openai  # type: ignore
        openai.api_key = api_key  # type: ignore[attr-defined]
        # Legacy interface ignores response_format
        resp = openai.ChatCompletion.create(  # type: ignore[attr-defined]
            model=model,
            messages=messages,
            temperature=temperature,
        )
        if 'choices' in resp and resp['choices']:
            return resp['choices'][0]['message']['content']
        return str(resp)
    except Exception as e:  # pragma: no cover
        raise RuntimeError(f"All OpenAI chat completion attempts failed: {e}")


def _build_mapping_schema() -> dict:
    """Return JSON schema dict for structured cell type mapping output.

    The model should return an object with one property 'mappings': array of {original, standard} objects.
    """
    return {
        "type": "json_schema",
        "json_schema": {
            "name": "cell_type_mapping",
            "schema": {
                "type": "object",
                "properties": {
                    "mappings": {
                        "type": "array",
                        "items": {
                            "type": "object",
                            "properties": {
                                "original": {"type": "string"},
                                "standard": {"type": "string"},
                            },
                            "required": ["original", "standard"],
                            "additionalProperties": False,
                        },
                        "minItems": 1,
                    }
                },
                "required": ["mappings"],
                "additionalProperties": False,
            },
        },
    }


def generate_cell_type_mapping(raw_labels: List[str], model: str = LLM_MODEL) -> Dict[str, str]:
    """Call LLM to produce mapping original_label -> standardized_label.

    Returns mapping. Uses best available mapping even if it exceeds MAX_STANDARD_LABELS.
    """
    unique_labels = sorted(set(l for l in raw_labels if l and l.lower() != 'nan'))
    user_prompt = (
        "Here are the cell type labels to harmonize (one per line):\n" + "\n".join(unique_labels)
    )
    messages = [
        {"role": "system", "content": LLM_SYSTEM_PROMPT + (" Return structured JSON per schema if provided." if LLM_USE_STRUCTURED else "")},
        {"role": "user", "content": user_prompt},
    ]

    last_err: Optional[Exception] = None
    adaptive_messages = list(messages)
    best_mapping: Optional[Dict[str, str]] = None
    best_unique_count = float('inf')
    
    for attempt in range(1, LLM_MAX_RETRIES + 1):
        try:
            response_format = _build_mapping_schema() if LLM_USE_STRUCTURED else None
            raw = _call_openai(adaptive_messages, model=model, response_format=response_format)
            # First try structured schema pathway (object with mappings array)
            parsed = json.loads(raw)
            mapping: Dict[str, str] = {}
            if isinstance(parsed, dict) and 'mappings' in parsed and isinstance(parsed['mappings'], list):
                for item in parsed['mappings']:
                    if isinstance(item, dict) and 'original' in item and 'standard' in item:
                        o = str(item['original']).strip()
                        s = str(item['standard']).strip()
                        if o:
                            mapping[o] = s or o
            elif isinstance(parsed, dict):
                # Fallback: assume direct mapping original->standard
                mapping = {str(k): str(v) for k, v in parsed.items()}
            else:
                raise ValueError("Unexpected JSON shape for mapping.")
            # Identity fill-ins
            for lbl in unique_labels:
                mapping.setdefault(lbl, lbl)
            unique_standard = list(dict.fromkeys(mapping.values()))
            
            # Track the best mapping regardless of whether it meets threshold
            if len(unique_standard) < best_unique_count:
                best_mapping = mapping
                best_unique_count = len(unique_standard)
            
            # If we meet the threshold, return immediately
            if len(unique_standard) <= MAX_STANDARD_LABELS:
                print(f"[INFO] LLM produced {len(unique_standard)} standardized labels (≤ {MAX_STANDARD_LABELS})")
                return mapping
            
            print(f"[WARN] Attempt {attempt}: LLM produced {len(unique_standard)} standardized labels (> {MAX_STANDARD_LABELS}). Retrying with stricter instructions.")
            if attempt < LLM_MAX_RETRIES:  # Only add reduction hint if we have more attempts
                overflow_sample = unique_standard[MAX_STANDARD_LABELS:MAX_STANDARD_LABELS+10]
                reduction_hint = (
                    f"You produced {len(unique_standard)} unique standardized labels which exceeds the limit of {MAX_STANDARD_LABELS}. "
                    f"Please merge semantically similar or rare subtypes so the final total is <= {MAX_STANDARD_LABELS}. "
                    f"Overflow examples: {overflow_sample}. Return ONLY the final mapping JSON."
                )
                adaptive_messages.append({"role": "user", "content": reduction_hint})
            continue
        except Exception as e:  # pragma: no cover - network dependent
            last_err = e
            print(f"[WARN] LLM attempt {attempt} failed: {e}")
            if attempt < LLM_MAX_RETRIES:
                print(f"[INFO] Sleeping {LLM_SLEEP_BETWEEN_RETRIES}s before retry...")
                time.sleep(LLM_SLEEP_BETWEEN_RETRIES)
    
    # If we couldn't meet the threshold, use the best mapping we got
    if best_mapping is not None:
        print(f"[WARN] Using best available mapping with {best_unique_count} standardized labels (> {MAX_STANDARD_LABELS} threshold)")
        return best_mapping
    
    # Final fallback: if all attempts failed, create identity mapping
    print(f"[ERROR] All LLM attempts failed, using identity mapping. Last error: {last_err}")
    return {lbl: lbl for lbl in unique_labels}


def apply_cell_type_mapping(adata: ad.AnnData, mapping: Dict[str, str], cell_type_key: str = CELL_TYPE_KEY) -> None:
    """In-place replacement of cell type labels using mapping."""
    if cell_type_key not in adata.obs:
        raise KeyError(f"'{cell_type_key}' not in adata.obs")
    adata.obs[cell_type_key] = adata.obs[cell_type_key].astype(str).map(mapping).astype('category')

# -------------------- DataFrame sanitation (obs) before writing ---------------
def sanitize_obs(adata: ad.AnnData) -> None:
    """Sanitize adata.obs columns to avoid h5py write errors.

    Converts object dtyped columns to either numeric (if all scalar numeric) or string.
    Ensures no python list/dict objects remain. Updates in place.
    """
    to_convert: List[str] = []
    for col in adata.obs.columns:
        ser = adata.obs[col]
        if ser.dtype == 'O':
            # Check if all entries are scalar numerics / NaNs
            def _is_scalar_numeric(x):
                return (isinstance(x, (int, float)) and not isinstance(x, bool)) or pd.isna(x)
            if ser.map(_is_scalar_numeric).all():
                try:
                    adata.obs[col] = pd.to_numeric(ser, errors='coerce')
                    print(f"[SANITIZE] Converted object numeric column '{col}' -> float")
                    continue
                except Exception:
                    pass
            # Fallback: stringify complex objects (lists, dicts, sets, tuples, mixed types)
            adata.obs[col] = ser.map(lambda v: 'NA' if pd.isna(v) else str(v))
            print(f"[SANITIZE] Stringified object column '{col}'")
        # Re-cast categories to avoid stale categories referencing objects
        if isinstance(adata.obs[col].dtype, pd.CategoricalDtype):
            adata.obs[col] = adata.obs[col].astype(str).astype('category')
    # Specific safeguard for commonly problematic count columns
    for special in ['nCount_RNA', 'nFeature_RNA']:
        if special in adata.obs:
            try:
                adata.obs[special] = pd.to_numeric(adata.obs[special], errors='coerce')
            except Exception:
                pass

# ------------------------------- Main Flow -------------------------------------

def main():
    # Set random seed for reproducibility
    np.random.seed(42)
    
    if DEBUG_MODE:
        print(f"[DEBUG] Running in DEBUG mode - using max {DEBUG_N_CELLS} cells per study")
        print(f"[DEBUG] Output files will have debug suffix")
        print(f"[DEBUG] Using smaller parameters for faster processing")
    else:
        print("[INFO] Running in FULL mode - using all cells")
    
    adatas = load_studies(STUDIES)
    if not adatas:
        print("[ERROR] No studies loaded; exiting.")
        return
    combined = concatenate(adatas)
    print(combined.obsm.keys())
    # Clean cell type labels after concatenation to ensure consistency
    print("[INFO] Final cell type label cleaning after concatenation...")
    clean_cell_type_labels(combined, CELL_TYPE_KEY)
    
    # Debug: show unique cell types after cleaning
    if DEBUG_MODE:
        unique_types = combined.obs[CELL_TYPE_KEY].unique()
        print(f"[DEBUG] Unique cell types after cleaning ({len(unique_types)}): {sorted(unique_types)}")

    mapping: Optional[Dict[str, str]] = None
    if RUN_LLM_STANDARDIZATION and CELL_TYPE_KEY in combined.obs:
        print("[INFO] Generating LLM-based cell type mapping...")
        try:
            original_labels = list(combined.obs[CELL_TYPE_KEY].astype(str).unique())
            print(f"[DEBUG] Original cell types ({len(original_labels)}): {sorted(original_labels)}")
            
            mapping = generate_cell_type_mapping(original_labels)
            print(f"[INFO] Mapping size: {len(mapping)}")
            
            # Apply the mapping
            original_obs = combined.obs[CELL_TYPE_KEY].copy()
            apply_cell_type_mapping(combined, mapping, CELL_TYPE_KEY)
            
            # Show detailed debugging information
            harmonized_labels = list(combined.obs[CELL_TYPE_KEY].astype(str).unique())
            print(f"[DEBUG] Harmonized cell types ({len(harmonized_labels)}): {sorted(harmonized_labels)}")
            
            # Show the actual mappings that were applied
            applied_mappings = {}
            for orig in original_labels:
                if orig in mapping:
                    applied_mappings[orig] = mapping[orig]
            
            print(f"[DEBUG] LLM mapping applied:")
            for orig, std in sorted(applied_mappings.items()):
                if orig != std:  # Only show mappings that actually changed
                    print(f"  '{orig}' → '{std}'")
            
            unchanged_count = sum(1 for orig, std in applied_mappings.items() if orig == std)
            changed_count = len(applied_mappings) - unchanged_count
            print(f"[DEBUG] {changed_count} labels harmonized, {unchanged_count} unchanged")
            
        except Exception as e:
            print(f"[WARN] LLM standardization skipped due to error: {e}")
            mapping = None
    else:
        if not RUN_LLM_STANDARDIZATION:
            print("[INFO] RUN_LLM_STANDARDIZATION is False; skipping label harmonization.")

    # Save artifacts
    out_path = os.path.join(BASE_DIR, OUTPUT_COMBINED_FILENAME)
    # Sanitize obs to prevent h5py TypeError (e.g., mixed object dtypes like 'nCount_RNA')
    sanitize_obs(combined)

    # Final debug information before running integration methods
    if DEBUG_MODE:
        final_cell_types = sorted(combined.obs[CELL_TYPE_KEY].astype(str).unique())
        print(f"[DEBUG] Final cell types ready for integration ({len(final_cell_types)}): {final_cell_types}")
        if mapping:
            print(f"[DEBUG] Using LLM-harmonized labels for all downstream integration methods")
        else:
            print(f"[DEBUG] Using original labels for all downstream integration methods")

    # Optional: Symphony inter-study evaluation (largest study as reference)
    if RUN_SYMPHONY:
        try:
            import importlib.util, importlib.machinery
            symphony_path = os.path.join(os.path.dirname(__file__), SYMPHONY_MODULE_FILENAME)
            if os.path.isfile(symphony_path):
                spec = importlib.util.spec_from_file_location('interstudy_symphony_mod', symphony_path)
                if spec and spec.loader:
                    symphony_mod = importlib.util.module_from_spec(spec)
                    spec.loader.exec_module(symphony_mod)  # type: ignore[attr-defined]
                    if hasattr(symphony_mod, 'run_symphony_interstudy'):
                        print('[INFO] Running Symphony inter-study mapping...')
                        # Adapt parameters for debug mode (smaller datasets need fewer HVGs)
                        n_top_genes = 1000 if DEBUG_MODE else 3000
                        n_pcs = 20 if DEBUG_MODE else 30
                        print(f"[INFO] Symphony params: n_top_genes={n_top_genes}, n_pcs={n_pcs}")
                        metrics_df, ref_study = symphony_mod.run_symphony_interstudy(
                            combined,
                            study_key=STUDY_KEY,
                            batch_key=LIBRARY_KEY,
                            cell_type_key=CELL_TYPE_KEY,
                            add_embedding=True,
                            embedding_key='X_symphony',
                            n_top_genes=n_top_genes,
                            n_pcs=n_pcs,
                        )
                        metrics_csv_path = os.path.join(BASE_DIR, SYMPHONY_METRICS_CSV)
                        metrics_df.to_csv(metrics_csv_path, index=False)
                        print(f"[INFO] Symphony metrics saved to {metrics_csv_path} (reference={ref_study})")
                    else:
                        print('[WARN] run_symphony_interstudy not found in Symphony module; skipping.')
                else:
                    print('[WARN] Could not load Symphony module spec; skipping.')
            else:
                print('[WARN] Symphony module file not found; skipping inter-study mapping.')
        except Exception as e:  # pragma: no cover
            print(f"[WARN] Symphony inter-study evaluation failed: {e}")

    # Optional: scPoli inter-study evaluation (errors will propagate)
    if RUN_SCPOLI:
        try:
            import importlib.util
            scpoli_path = os.path.join(os.path.dirname(__file__), SCPOLI_MODULE_FILENAME)
            if os.path.isfile(scpoli_path):
                spec = importlib.util.spec_from_file_location('interstudy_scpoli_mod', scpoli_path)
                if spec and spec.loader:
                    scpoli_mod = importlib.util.module_from_spec(spec)
                    spec.loader.exec_module(scpoli_mod)  # type: ignore[attr-defined]
                    if hasattr(scpoli_mod, 'run_scpoli_interstudy'):
                        print('[INFO] Running scPoli inter-study integration...')
                        latent_dim = 20 if DEBUG_MODE else 30
                        pretraining_epochs = 10 if DEBUG_MODE else 40
                        total_epochs = 20 if DEBUG_MODE else 50
                        print(f"[INFO] scPoli params: latent_dim={latent_dim}, epochs={total_epochs}")
                        
                        # Check obsm keys before scPoli
                        obsm_before = set(combined.obsm.keys())
                        print(f"[DEBUG] obsm keys before scPoli: {sorted(obsm_before)}")
                        
                        scpoli_metrics_df, scpoli_ref = scpoli_mod.run_scpoli_interstudy(
                            combined,
                            study_key=STUDY_KEY,
                            batch_key=LIBRARY_KEY,
                            cell_type_key=CELL_TYPE_KEY,
                            embedding_key='X_scpoli',
                            latent_dim=latent_dim,
                            pretraining_epochs=pretraining_epochs,
                            total_epochs=total_epochs,
                        )
                        
                        # Check obsm keys after scPoli
                        obsm_after = set(combined.obsm.keys())
                        new_keys = obsm_after - obsm_before
                        print(f"[DEBUG] obsm keys after scPoli: {sorted(obsm_after)}")
                        if new_keys:
                            print(f"[DEBUG] New obsm keys added by scPoli: {sorted(new_keys)}")
                        else:
                            print("[WARN] No new obsm keys added by scPoli!")
                        
                        scpoli_csv = os.path.join(BASE_DIR, SCPOLI_METRICS_CSV)
                        scpoli_metrics_df.to_csv(scpoli_csv, index=False)
                        print(f"[INFO] scPoli metrics saved to {scpoli_csv} (reference={scpoli_ref})")
                        
                        # Debug: Check if scPoli embedding was actually added
                        if 'X_scpoli' in combined.obsm:
                            print(f"[DEBUG] scPoli embedding successfully added to obsm with shape {combined.obsm['X_scpoli'].shape}")
                        else:
                            print("[ERROR] scPoli embedding 'X_scpoli' not found in obsm after execution!")
                            print(f"[DEBUG] Available obsm keys after scPoli: {list(combined.obsm.keys())}")
                            # Check for any scpoli-related keys
                            scpoli_keys = [k for k in combined.obsm.keys() if 'scpoli' in k.lower()]
                            if scpoli_keys:
                                print(f"[DEBUG] Found alternative scPoli keys: {scpoli_keys}")
                            else:
                                print("[DEBUG] No scPoli-related keys found")
                    else:
                        raise AttributeError('run_scpoli_interstudy not found in scPoli module')
                else:
                    raise ImportError('Could not load scPoli module spec')
            else:
                raise FileNotFoundError(f'scPoli module file not found: {scpoli_path}')
        except Exception as e:
            print(f"[ERROR] scPoli inter-study evaluation failed: {e}")
            import traceback
            traceback.print_exc()
            print("[WARN] Continuing without scPoli embedding")

    # Optional: scArches (SCANVI surgery) inter-study evaluation (errors will propagate)
    if RUN_SCARCHES:
        import importlib.util
        scarches_path = os.path.join(os.path.dirname(__file__), SCARCHES_MODULE_FILENAME)
        if os.path.isfile(scarches_path):
            spec = importlib.util.spec_from_file_location('interstudy_scarches_mod', scarches_path)
            if spec and spec.loader:
                scarches_mod = importlib.util.module_from_spec(spec)
                spec.loader.exec_module(scarches_mod)  # type: ignore[attr-defined]
                if hasattr(scarches_mod, 'run_scarches_interstudy'):
                    print('[INFO] Running scArches (SCANVI surgery) inter-study integration...')
                    n_epochs = 20 if DEBUG_MODE else 50
                    n_latent = 20 if DEBUG_MODE else 30
                    print(f"[INFO] scArches params: n_epochs={n_epochs}, n_latent={n_latent}")
                    scarches_metrics_df, scarches_ref = scarches_mod.run_scarches_interstudy(
                        combined,
                        base_dir=BASE_DIR,
                        study_key=STUDY_KEY,
                        batch_key=LIBRARY_KEY,
                        cell_type_key=CELL_TYPE_KEY,
                        embedding_key='X_scarches',
                        reference_max_epochs=n_epochs,
                        latent_dim=n_latent,
                    )
                    scarches_csv = os.path.join(BASE_DIR, SCARCHES_METRICS_CSV)
                    scarches_metrics_df.to_csv(scarches_csv, index=False)
                    print(f"[INFO] scArches metrics saved to {scarches_csv} (reference={scarches_ref})")
                    
                    # Debug: Check if scArches embedding was actually added
                    if 'X_scarches' in combined.obsm:
                        print(f"[DEBUG] scArches embedding successfully added to obsm with shape {combined.obsm['X_scarches'].shape}")
                    else:
                        print("[ERROR] scArches embedding 'X_scarches' not found in obsm after execution!")
                        print(f"[DEBUG] Available obsm keys after scArches: {list(combined.obsm.keys())}")
                else:
                    raise AttributeError('run_scarches_interstudy not found in scArches module')
            else:
                raise ImportError('Could not load scArches module spec')
        else:
            raise FileNotFoundError(f'scArches module file not found: {scarches_path}')

    if mapping:
        mapping_path = os.path.join(BASE_DIR, LABEL_MAPPING_JSON)
        print(f"[INFO] Writing label mapping JSON to {mapping_path}")
        with open(mapping_path, 'w') as f:
            json.dump(mapping, f, indent=2)

    # scIB Benchmarker evaluation (like intra-studies script)
    if RUN_SCIB_BENCHMARKER:
        try:
            import gc
            from scib_metrics.benchmark import Benchmarker, BioConservation, BatchCorrection  # type: ignore[import]
            
            print('[INFO] Running scIB Benchmarker evaluation...')
            
            # Debug: Comprehensive check of integration method results
            print("[DEBUG] Checking integration method results before scIB evaluation:")
            expected_embeddings = {
                'X_symphony': RUN_SYMPHONY,
                'X_scpoli': RUN_SCPOLI, 
                'X_scarches': RUN_SCARCHES
            }
            
            for embedding_key, was_enabled in expected_embeddings.items():
                if was_enabled:
                    if embedding_key in combined.obsm:
                        print(f"[DEBUG] ✓ {embedding_key}: Present (shape: {combined.obsm[embedding_key].shape})")
                    else:
                        print(f"[DEBUG] ✗ {embedding_key}: MISSING despite being enabled!")
                else:
                    print(f"[DEBUG] - {embedding_key}: Disabled (RUN flag was False)")
            
            print(f"[DEBUG] All obsm keys: {sorted(combined.obsm.keys())}")
            
            # Collect available embeddings from integration methods
            available_embeddings = []
            integration_methods = [
                ('X_symphony', 'symphony'),
                ('X_scpoli', 'scpoli'), 
                ('X_scarches', 'scarches'),
                # Foundation model embeddings if present
                ('x_geneformer_helical', 'geneformer_helical'),
                ('x_scgpt_helical', 'scgpt_helical'),
                ('x_uce', 'uce'),
                ('x_scfoundation', 'scfoundation'),
                ('x_transcriptformer', 'transcriptformer'),
                ('x_scimilarity', 'scimilarity'),
                # Standard integration methods
                ('X_pca', 'pca_unintegrated'),
            ]
            
            available_obsms = []
            for obsm_key, method_name in integration_methods:
                if obsm_key in combined.obsm:
                    available_obsms.append(obsm_key)
                    print(f"[INFO] Found embedding: {obsm_key} ({method_name})")
                else:
                    print(f"[DEBUG] Missing embedding: {obsm_key} ({method_name})")
            
            print('Here are all available obsm keys:', list(combined.obsm.keys()))
            
            # Special debugging for scPoli
            if 'X_scpoli' not in available_obsms:
                print("[WARN] X_scpoli embedding not found!")
                print("[DEBUG] Checking if scPoli ran successfully...")
                if RUN_SCPOLI:
                    print("[DEBUG] RUN_SCPOLI is True - scPoli should have run")
                    # Check if there are any scpoli-related keys
                    scpoli_keys = [k for k in combined.obsm.keys() if 'scpoli' in k.lower()]
                    if scpoli_keys:
                        print(f"[DEBUG] Found scPoli-related keys: {scpoli_keys}")
                        # Add the first scpoli key we find and update integration_methods
                        for scpoli_key in scpoli_keys:
                            if scpoli_key not in available_obsms:
                                available_obsms.append(scpoli_key)
                                print(f"[INFO] Added {scpoli_key} as scPoli embedding")
                                # Also add to integration_methods for proper tracking
                                integration_methods.append((scpoli_key, f'scpoli_{scpoli_key}'))
                    else:
                        print("[DEBUG] No scPoli-related keys found in obsm")
                        print("[DEBUG] This suggests scPoli integration failed or didn't store results")
                else:
                    print("[DEBUG] RUN_SCPOLI is False - scPoli was skipped")
            
            # Also check for any other integration method embeddings that might be missing
            for embed_key, method_name in [('X_symphony', 'symphony'), ('X_scarches', 'scarches')]:
                if embed_key not in available_obsms and embed_key in combined.obsm:
                    print(f"[DEBUG] Found {embed_key} in obsm but missed in initial detection")
                    available_obsms.append(embed_key)
            if not available_obsms:
                print('[WARN] No embeddings found for scIB evaluation; skipping.')
            else:
                print(f'[INFO] Evaluating {len(available_obsms)} embeddings with scIB...')
                
                # Check and clean embeddings for NaN values before scIB evaluation
                print("[INFO] Checking embeddings for NaN values...")
                clean_obsms = []
                for obsm_key in available_obsms:
                    embedding = combined.obsm[obsm_key]
                    n_nan = np.isnan(embedding).sum()
                    n_inf = np.isinf(embedding).sum()
                    
                    if n_nan > 0 or n_inf > 0:
                        print(f"[WARN] {obsm_key}: Found {n_nan} NaN and {n_inf} infinite values, cleaning...")
                        # Replace NaN and infinite values with zeros
                        embedding_clean = np.nan_to_num(embedding, nan=0.0, posinf=0.0, neginf=0.0)
                        
                        # Ensure float64 dtype for stability
                        embedding_clean = embedding_clean.astype(np.float64)
                        
                        # Final verification
                        if np.isnan(embedding_clean).sum() == 0 and np.isinf(embedding_clean).sum() == 0:
                            combined.obsm[obsm_key] = embedding_clean
                            clean_obsms.append(obsm_key)
                            print(f"[INFO] {obsm_key}: Successfully cleaned problematic values")
                        else:
                            print(f"[ERROR] {obsm_key}: Failed to clean problematic values, excluding from evaluation")
                    else:
                        print(f"[INFO] {obsm_key}: No problematic values detected")
                        # Still ensure consistent dtype
                        combined.obsm[obsm_key] = embedding.astype(np.float64)
                        clean_obsms.append(obsm_key)
                
                if not clean_obsms:
                    print('[ERROR] No clean embeddings available for scIB evaluation; skipping.')
                else:
                    print(f'[INFO] Using {len(clean_obsms)} clean embeddings for scIB evaluation: {clean_obsms}')
                    
                    # Additional data validation
                    print("[INFO] Validating AnnData object for scIB compatibility...")
                    
                    # Ensure we have PCA as baseline if it's missing
                    if 'X_pca' not in combined.obsm and 'X_pca' not in clean_obsms:
                        print("[INFO] X_pca missing, computing PCA for baseline comparison...")
                        try:
                            # Ensure we have the raw or log-normalized data for PCA
                            if combined.X is not None:
                                sc.tl.pca(combined, n_comps=50)
                                if 'X_pca' in combined.obsm:
                                    clean_obsms.append('X_pca')
                                    print("[INFO] Successfully computed X_pca baseline")
                                else:
                                    print("[WARN] PCA computation did not produce X_pca")
                            else:
                                print("[WARN] No .X data available for PCA computation")
                        except Exception as e:
                            print(f"[WARN] PCA computation failed: {e}")
                    
                    # Check for missing observations in key columns
                    for key in [STUDY_KEY, CELL_TYPE_KEY]:
                        if key in combined.obs:
                            n_na = combined.obs[key].isna().sum()
                            if n_na > 0:
                                print(f"[WARN] Found {n_na} missing values in {key}, filling with 'Unknown'")
                                combined.obs[key] = combined.obs[key].fillna('Unknown').astype('category')
                            else:
                                print(f"[INFO] {key}: No missing values")
                                combined.obs[key] = combined.obs[key].astype('category')
                    
                    print("[INFO] Setting up scIB Benchmarker...")
                    
                    # Force garbage collection before intensive operations
                    gc.collect()
                    
                    # Setup scIB metrics
                    biocons = BioConservation(isolated_labels=False)
                    
                    print(f"[INFO] Final embeddings for scIB evaluation: {clean_obsms}")
                    print(f"[INFO] Batch key: {STUDY_KEY}, Label key: {CELL_TYPE_KEY}")
                    print(f"[INFO] Pre-integrated baseline: {'X_pca' if 'X_pca' in clean_obsms else clean_obsms[0]}")
                    
                    bm = Benchmarker(
                        combined,
                        batch_key=STUDY_KEY,  # Use study as batch key for inter-study evaluation
                        label_key=CELL_TYPE_KEY,
                        embedding_obsm_keys=clean_obsms,  # Use cleaned embeddings
                        pre_integrated_embedding_obsm_key="X_pca" if "X_pca" in clean_obsms else clean_obsms[0],
                        bio_conservation_metrics=biocons,
                        batch_correction_metrics=BatchCorrection(),
                        n_jobs=-1,
                    )
                    
                    print("[INFO] Preparing scIB benchmarker...")
                    try:
                        bm.prepare()
                        print("[INFO] scIB benchmarker preparation successful")
                    except Exception as e:
                        print(f"[ERROR] scIB benchmarker preparation failed: {e}")
                        print("[INFO] Attempting to continue with benchmark...")
                    
                    print("[INFO] Running scIB benchmark...")
                    try:
                        bm.benchmark()
                        print("[INFO] scIB benchmark execution successful")
                    except Exception as e:
                        print(f"[ERROR] scIB benchmark execution failed: {e}")
                        raise
                    
                    print("[INFO] Generating scIB results...")
                    gc.collect()
                    
                    # Create figures directory
                    figures_dir = os.path.join(BASE_DIR, 'interstudy_figures')
                    os.makedirs(figures_dir, exist_ok=True)
                    
                    # Plot results table
                    bm.plot_results_table(
                        show=False, 
                        min_max_scale=False,
                        save_dir=figures_dir
                    )
                    
                    # Save results to CSV
                    results_df = bm.get_results()
                    results_csv_path = os.path.join(figures_dir, f'scib_results_interstudy_{"debug" if DEBUG_MODE else "full"}.csv')
                    results_df.to_csv(results_csv_path)
                    
                    print(f"[INFO] scIB evaluation complete")
                    print(f"[INFO] Results saved to: {results_csv_path}")
                    print(f"[INFO] Figures saved to: {figures_dir}")
                    
                    # Print summary of results
                    print(f"[INFO] scIB Results Summary:")
                    print(results_df.round(3))
                    
                    # Upload results if configured
                    try:
                        import subprocess
                        if os.path.exists(figures_dir):
                            gdrive_path = 'GDrive:KidneyAtlas/interstudy_figures/'
                            subprocess.run(['rclone', 'copy', '--progress', figures_dir, gdrive_path])
                            print(f"[INFO] Results uploaded to {gdrive_path}")
                    except Exception as e:
                        print(f"[WARN] Upload failed: {e}")
                        
                    gc.collect()
                
        except Exception as e:
            print(f"[ERROR] scIB Benchmarker evaluation failed: {e}")
            import traceback
            traceback.print_exc()

    # Optional: scIB benchmarking
    if RUN_SCIB:
        try:
            import importlib.util
            scib_path = os.path.join(os.path.dirname(__file__), SCIB_MODULE_FILENAME)
            if os.path.isfile(scib_path):
                spec = importlib.util.spec_from_file_location('interstudy_scib_mod', scib_path)
                if spec and spec.loader:
                    scib_mod = importlib.util.module_from_spec(spec)
                    spec.loader.exec_module(scib_mod)  # type: ignore[attr-defined]
                    if hasattr(scib_mod, 'run_scib_benchmark'):
                        print('[INFO] Running scIB benchmarking over embeddings...')
                        present_embeddings = [k for k in SCIB_EMBEDDINGS if k in combined.obsm]
                        if not present_embeddings:
                            print('[WARN] No specified scIB embeddings present in adata.obsm; skipping scIB.')
                        else:
                            scib_df = scib_mod.run_scib_benchmark(
                                combined,
                                embedding_keys=present_embeddings,
                                batch_keys=SCIB_BATCH_KEYS,
                                label_key=CELL_TYPE_KEY,
                                fast=True,
                                slim=False,
                                full=False,
                            )
                            scib_csv = os.path.join(BASE_DIR, SCIB_METRICS_CSV)
                            scib_df.to_csv(scib_csv, index=False)
                            print(f"[INFO] scIB metrics saved to {scib_csv}")
                    else:
                        print('[WARN] run_scib_benchmark not found in scIB module; skipping.')
                else:
                    print('[WARN] Could not load scIB module spec; skipping scIB benchmarking.')
            else:
                print('[WARN] scIB module file not found; skipping scIB benchmarking.')
        except Exception as e:  # pragma: no cover
            print(f"[WARN] scIB benchmarking failed: {e}")
    print(f"[INFO] Writing combined AnnData to {out_path}")
    combined.write(out_path)
    print("[DONE] Inter-study combination complete.")

if __name__ == '__main__':  # pragma: no cover
    main()
    adata = ad.read_h5ad(os.path.join(BASE_DIR, OUTPUT_COMBINED_FILENAME))
    print(f"[FINAL] Combined AnnData loaded with shape {adata.shape} and obs keys: {', '.join(adata.obs.keys())}")
    
    # Always show final cell types regardless of whether mapping was used
    adata.obs[CELL_TYPE_KEY] = adata.obs[CELL_TYPE_KEY].astype(str)
    unique_cell_types = sorted(adata.obs[CELL_TYPE_KEY].unique())
    print(f"[FINAL] Unique cell types in combined data: {len(unique_cell_types)}")
    print(f"[FINAL] Cell types: {', '.join(unique_cell_types)}")
    
    if LABEL_MAPPING_JSON in os.listdir(BASE_DIR):
        with open(os.path.join(BASE_DIR, LABEL_MAPPING_JSON), 'r') as f:
            mapping = json.load(f)
        print(f"[FINAL] Loaded cell type mapping with {len(mapping)} entries from JSON")
        
        # Show how many original labels were harmonized
        changed_mappings = {k: v for k, v in mapping.items() if k != v}
        if changed_mappings:
            print(f"[FINAL] {len(changed_mappings)} labels were harmonized by LLM:")
            for orig, harmonized in sorted(changed_mappings.items()):
                print(f"  '{orig}' → '{harmonized}'")
        else:
            print(f"[FINAL] All labels remained unchanged (identity mapping)")
    else:
        print(f"[FINAL] No label mapping JSON found - original labels were used")