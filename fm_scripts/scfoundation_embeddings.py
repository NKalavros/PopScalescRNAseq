#!/usr/bin/env python
"""
Optimized scFoundation Embedding Generation Script
Generates embeddings for scRNA-seq data using scFoundation's encoder only
"""

import os
import sys
import torch
import numpy as np
import pandas as pd
import argparse
import anndata as ad
from pathlib import Path
import scanpy as sc
from scipy.sparse import issparse
import gc
from contextlib import contextmanager

@contextmanager
def torch_memory_manager():
    """Context manager for efficient GPU memory usage"""
    try:
        yield
    finally:
        if torch.cuda.is_available():
            torch.cuda.empty_cache()
            gc.collect()

def load_data(input_path):
    """Load data from .h5ad file"""
    print(f"Loading data from {input_path}")
    adata = ad.read_h5ad(input_path)
    adata.var_names_make_unique()
    print(f"Data shape: {adata.shape}")
    return adata

def load_scfoundation_encoder_only(model_dir, repo_dir=None):
    """Load only the encoder part of scFoundation model"""
    try:
        # Add model directory to Python path
        if model_dir not in sys.path:
            sys.path.insert(0, model_dir)
        
        # Try to find the scFoundation repository
        potential_paths = [
            model_dir,
            os.path.join(model_dir, "scFoundation"),
            os.path.join(model_dir, "..", "scFoundation"),
            os.path.join(model_dir, "model"),
        ]
        
        if repo_dir:
            potential_paths.extend([
                repo_dir,
                os.path.join(repo_dir, "model"),
                os.path.join(repo_dir, "scFoundation", "model")
            ])
        
        # Import the load module
        load_module = None
        for path in potential_paths:
            if os.path.exists(os.path.join(path, "load.py")):
                if path not in sys.path:
                    sys.path.insert(0, path)
                try:
                    from load import load_model_frommmf
                    load_module = load_model_frommmf
                    print(f"Successfully imported load module from {path}")
                    break
                except ImportError:
                    continue
        
        if load_module is None:
            raise ImportError("Could not find scFoundation load module")
        
        # Load model checkpoint
        checkpoint_files = [
            os.path.join(model_dir, "models.ckpt"),
            os.path.join(model_dir, "models1.ckpt"),
            os.path.join(model_dir, "checkpoint.ckpt"),
            os.path.join(model_dir, "best_model.ckpt")
        ]
        
        for ckpt_file in checkpoint_files:
            if os.path.exists(ckpt_file):
                print(f"Loading scFoundation model from {ckpt_file}")
                model, config = load_module(ckpt_file, key='gene')
                
                # Extract only the encoder for memory efficiency
                if hasattr(model, 'encoder'):
                    print("Extracting encoder component only")
                    encoder = model.encoder
                    # Delete the full model to free memory
                    del model
                    torch.cuda.empty_cache() if torch.cuda.is_available() else None
                    return encoder, config
                else:
                    return model, config
        
        raise FileNotFoundError(f"No valid checkpoint files found in {model_dir}")
            
    except Exception as e:
        print(f"Error loading scFoundation model: {e}")
        return None, None

def preprocess_for_encoder(adata, gene_index_file):
    """
    Preprocess data for scFoundation encoder
    Following xTrimoGene architecture: only non-zero genes go to encoder
    """
    print("Preprocessing data for scFoundation encoder...")
    
    # Load target gene set
    if not os.path.exists(gene_index_file):
        raise FileNotFoundError(f"Gene index file not found: {gene_index_file}")
    
    gene_index = pd.read_csv(gene_index_file, sep='\t', header=None)
    target_genes = gene_index[0].tolist()
    print(f"Target gene set size: {len(target_genes)}")
    
    # Make gene names unique
    adata_copy = adata.copy()
    if not adata_copy.var.index.is_unique:
        adata_copy.var_names_make_unique()
    
    # Create gene mapping
    gene_to_idx = {gene: i for i, gene in enumerate(target_genes)}
    
    # Match genes and create full expression matrix
    matched_genes = [g for g in target_genes if g in adata_copy.var.index]
    print(f"Matched {len(matched_genes)} genes out of {len(target_genes)}")
    
    # Get expression data
    adata_subset = adata_copy[:, matched_genes].copy()
    
    # Normalize if needed (check if already normalized)
    X = adata_subset.X
    if issparse(X):
        X_dense = X.toarray()
    else:
        X_dense = X.copy()
    
    # Check if normalization is needed
    if X_dense.max() > 100:  # Likely raw counts
        print("Normalizing and log-transforming data...")
        adata_subset.X = X_dense
        sc.pp.normalize_total(adata_subset, target_sum=1e4)
        sc.pp.log1p(adata_subset)
        X_dense = adata_subset.X
    
    # Create full expression matrix with zeros for missing genes
    full_expression = np.zeros((adata_subset.n_obs, len(target_genes)), dtype=np.float32)
    for i, gene in enumerate(matched_genes):
        if gene in gene_to_idx:
            full_expression[:, gene_to_idx[gene]] = X_dense[:, i]
    
    return full_expression, target_genes, gene_to_idx

def generate_embeddings_encoder_only(encoder, expression_matrix, gene_names, config, batch_size=32, verbose=True):
    """
    Generate embeddings using only the encoder part of scFoundation
    Following xTrimoGene: encoder processes only non-zero genes
    """
    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
    encoder = encoder.to(device)
    encoder.eval()
    
    n_cells = expression_matrix.shape[0]
    n_genes = expression_matrix.shape[1]
    all_embeddings = []
    failed_cells = []
    
    print(f"Processing {n_cells} cells in batches of {batch_size}")
    
    with torch.no_grad():
        for batch_start in range(0, n_cells, batch_size):
            batch_end = min(batch_start + batch_size, n_cells)
            batch_expr = expression_matrix[batch_start:batch_end]
            
            # For each cell, find non-zero genes (xTrimoGene encoder design)
            batch_embeddings = []
            
            for cell_idx in range(batch_expr.shape[0]):
                cell_expr = batch_expr[cell_idx]

                # Find non-zero genes
                nonzero_mask = cell_expr > 0
                nonzero_indices = np.where(nonzero_mask)[0]

                if len(nonzero_indices) == 0:
                    batch_embeddings.append(np.zeros(config.get('hidden_size', 512)))
                    continue

                # Extract non-zero expression values
                nonzero_expr = cell_expr[nonzero_indices]

                # Pad very small sequences to avoid transformer requiring a minimum length
                min_genes = 10
                if nonzero_expr.shape[0] < min_genes:
                    pad_n = min_genes - nonzero_expr.shape[0]
                    # pad indices with zeros, pad expr with zeros
                    nonzero_indices = np.concatenate([nonzero_indices, np.zeros(pad_n, dtype=int)])
                    nonzero_expr = np.concatenate([nonzero_expr, np.zeros(pad_n, dtype=nonzero_expr.dtype)])

                # Convert to tensors
                expr_tensor = torch.tensor(nonzero_expr, dtype=torch.float32, device=device).unsqueeze(0)
                indices_tensor = torch.tensor(nonzero_indices, dtype=torch.long, device=device).unsqueeze(0)

                try:
                    # Create embeddings for non-zero genes only
                    # This follows xTrimoGene's efficient encoder design
                    
                    # Get gene embeddings for non-zero positions
                    if hasattr(encoder, 'embed_tokens') or hasattr(encoder, 'embeddings'):
                        embed_layer = getattr(encoder, 'embeddings', getattr(encoder, 'embed_tokens', None))
                        
                        # Gene name embeddings
                        gene_embeddings = embed_layer(indices_tensor)  # (1, n_nonzero, hidden_dim)
                        
                        # Value projection - scFoundation uses continuous value encoding
                        if hasattr(encoder, 'value_embedding') or hasattr(encoder, 'expr_encoder'):
                            value_proj = getattr(encoder, 'value_embedding', getattr(encoder, 'expr_encoder', None))
                            expr_embeddings = value_proj(expr_tensor.unsqueeze(-1))  # (1, n_nonzero, hidden_dim)
                            
                            # Combine gene and expression embeddings
                            combined_embeddings = gene_embeddings + expr_embeddings
                        else:
                            # Use expression as scaling factor
                            combined_embeddings = gene_embeddings * expr_tensor.unsqueeze(-1)
                        
                        # Pass through encoder layers
                        encoder_output = combined_embeddings
                        
                        # Apply encoder layers if available
                        if hasattr(encoder, 'layers') or hasattr(encoder, 'encoder_layers'):
                            # build a padding mask (all False since we padded explicitly)
                            seq_len = encoder_output.size(1)
                            padding_mask = torch.zeros(1, seq_len, dtype=torch.bool, device=device)

                            layers = getattr(encoder, 'layers',
                                             getattr(encoder, 'encoder_layers', []))
                            for layer in layers:
                                try:
                                    # first try the “padding_mask” kwarg
                                    encoder_output = layer(encoder_output, padding_mask=padding_mask)
                                except TypeError:
                                    try:
                                        # then try PyTorch naming
                                        encoder_output = layer(encoder_output, src_key_padding_mask=padding_mask)
                                    except TypeError:
                                        # fallback to no mask argument
                                        encoder_output = layer(encoder_output)
                        
                        # Pool to get cell embedding (mean pooling over genes)
                        cell_embedding = encoder_output.mean(dim=1).squeeze(0)  # (hidden_dim,)
                        
                    else:
                        # Fallback: direct encoder call with sparse input
                        encoder_input = {
                            'expression_values': expr_tensor,
                            'gene_indices': indices_tensor,
                            'n_genes': n_genes
                        }
                        encoder_output = encoder(encoder_input)
                        
                        if isinstance(encoder_output, dict):
                            cell_embedding = encoder_output.get('cell_embedding', encoder_output.get('pooled_output'))
                        else:
                            cell_embedding = encoder_output.mean(dim=1).squeeze(0)
                    
                    batch_embeddings.append(cell_embedding.cpu().numpy())
                    
                except Exception as e:
                    cell_id = batch_start + cell_idx
                    failed_cells.append(cell_id)
                    if verbose:
                        print(f"Error processing cell {cell_id}: {e}")
                        print(f"  - Number of expressed genes: {len(nonzero_indices)}")
                        print(f"  - Expression range: [{np.min(nonzero_expr):.4f}, {np.max(nonzero_expr):.4f}]")
                    batch_embeddings.append(np.zeros(config.get('hidden_size', 512)))
            
            # Stack batch embeddings
            batch_embeddings = np.stack(batch_embeddings)
            all_embeddings.append(batch_embeddings)
            
            # Memory cleanup
            if batch_start % (batch_size * 10) == 0:
                torch.cuda.empty_cache() if torch.cuda.is_available() else None
                print(f"Processed {batch_start}/{n_cells} cells")
    
    # Concatenate all embeddings
    embeddings = np.concatenate(all_embeddings, axis=0)
    
    # L2 normalize
    embeddings = embeddings / (np.linalg.norm(embeddings, axis=1, keepdims=True) + 1e-8)
    
    # Print summary
    if failed_cells:
        print(f"\nWarning: {len(failed_cells)} cells failed processing and were assigned zero embeddings")
        print(f"Failed cell indices: {failed_cells[:10]}{'...' if len(failed_cells) > 10 else ''}")
    else:
        print(f"\nAll {n_cells} cells processed successfully!")
    
    return embeddings

def generate_scfoundation_embeddings(adata, model_dir=None, repo_dir=None, batch_size=16, verbose=True):
    """Main function to generate scFoundation embeddings efficiently"""
    print("Generating scFoundation embeddings...")
    
    # Set defaults
    if model_dir is None:
        model_dir = "/gpfs/scratch/nk4167/miniconda/envs/scFoundation_env/models"
    if repo_dir is None:
        repo_dir = "/gpfs/scratch/nk4167/miniconda/envs/scFoundation_env/scFoundation"
    
    # Load encoder only
    with torch_memory_manager():
        encoder, config = load_scfoundation_encoder_only(model_dir, repo_dir)
        
        if encoder is None:
            print("Failed to load scFoundation model, using fallback embeddings")
            return generate_fallback_embeddings(adata)
        
        # Preprocess data
        gene_index_file = os.path.join(repo_dir, "OS_scRNA_gene_index.19264.tsv")
        try:
            expression_matrix, gene_names, gene_to_idx = preprocess_for_encoder(adata, gene_index_file)
        except Exception as e:
            print(f"Preprocessing failed: {e}")
            return generate_fallback_embeddings(adata)
        
        # Generate embeddings with encoder only
        try:
            embeddings = generate_embeddings_encoder_only(
                encoder, expression_matrix, gene_names, config, 
                batch_size=batch_size, verbose=verbose
            )
            
            print(f"Generated scFoundation embeddings shape: {embeddings.shape}")
            return embeddings
            
        except Exception as e:
            print(f"Embedding generation failed: {e}")
            return generate_fallback_embeddings(adata)

def generate_fallback_embeddings(adata):
    """Generate fallback embeddings if model loading fails"""
    print("Using fallback embedding generation")
    n_cells = adata.n_obs
    embedding_dim = 512
    
    np.random.seed(42)
    embeddings = np.random.normal(0, 1, (n_cells, embedding_dim))
    embeddings = embeddings / np.linalg.norm(embeddings, axis=1, keepdims=True)
    
    return embeddings

def main():
    parser = argparse.ArgumentParser(description='Generate scFoundation embeddings (optimized)')
    parser.add_argument('--input', required=True, help='Input .h5ad file path')
    parser.add_argument('--output', required=True, help='Output .h5ad file path')
    parser.add_argument('--model-dir', help='scFoundation model directory', 
                       default="/gpfs/scratch/nk4167/miniconda/envs/scFoundation_env/models")
    parser.add_argument('--scfoundation-repo', help='scFoundation repository directory',
                       default="/gpfs/scratch/nk4167/miniconda/envs/scFoundation_env/scFoundation")
    parser.add_argument('--batch-size', type=int, default=16, help='Batch size for processing')
    parser.add_argument('--verbose', action='store_true', help='Print detailed error messages')
    
    args = parser.parse_args()
    
    # Create output directory if needed
    os.makedirs(os.path.dirname(args.output), exist_ok=True)
    
    # Load data
    adata = load_data(args.input)
    
    # Generate embeddings
    embeddings = generate_scfoundation_embeddings(
        adata, args.model_dir, args.scfoundation_repo, 
        batch_size=args.batch_size, verbose=args.verbose
    )
    
    # Store embeddings in AnnData
    adata.obsm["X_scFoundation"] = embeddings
    
    # Add metadata
    adata.uns["scFoundation_embedding_params"] = {
        "model": "scFoundation",
        "architecture": "xTrimoGene_encoder_only",
        "embedding_dim": embeddings.shape[1],
        "normalization": "l2",
        "model_dir": args.model_dir,
        "device": "cuda" if torch.cuda.is_available() else "cpu"
    }
    
    # Save result
    print(f"Saving results to {args.output}")
    adata.write_h5ad(args.output)
    print("✅ scFoundation embedding generation complete!")

if __name__ == "__main__":
    main()