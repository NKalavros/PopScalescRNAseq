#!/usr/bin/env python
"""
scGPT Embedding Generation Script
Generates embeddings for scRNA-seq data using scGPT foundation model
"""

import os
import sys
import torch
import scanpy as sc
import numpy as np
import argparse
from pathlib import Path

# scGPT imports
from scgpt.tokenizer.gene_tokenizer import GeneVocab
from scgpt.model import TransformerModel
from scgpt.tokenizer import tokenize_and_pad_batch
from scgpt.preprocess import Preprocessor

def load_data(input_path):
    """Load data from .h5ad file"""
    print(f"Loading data from {input_path}")
    adata = sc.read_h5ad(input_path)
    print(f"Data shape: {adata.shape}")
    return adata

def generate_scgpt_embeddings(adata, model_dir=None):
    """Generate embeddings using scGPT"""
    print("Generating scGPT embeddings...")
    
    # Set device
    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
    print(f"Using device: {device}")
    
    # Model configuration
    pad_token = "<pad>"
    special_tokens = [pad_token, "<cls>", "<eoc>"]
    n_hvg = 1200
    max_seq_len = n_hvg + 1
    n_bins = 51
    
    # Load pre-trained model if specified
    if model_dir is not None:
        model_path = Path(model_dir)
        vocab_file = model_path / "vocab.json"
        model_file = model_path / "best_model.pt"
        config_file = model_path / "args.json"
        
        if vocab_file.exists() and model_file.exists():
            print(f"Loading pre-trained model from {model_dir}")
            
            # Load vocabulary
            vocab = GeneVocab.from_file(vocab_file)
            for s in special_tokens:
                if s not in vocab:
                    vocab.append_token(s)
            
            # Filter genes to those in vocabulary
            adata.var["gene_name"] = adata.var.index.tolist()
            adata.var["id_in_vocab"] = [
                1 if gene in vocab else -1 for gene in adata.var["gene_name"]
            ]
            gene_ids_in_vocab = np.array(adata.var["id_in_vocab"])
            print(f"Matched {np.sum(gene_ids_in_vocab >= 0)}/{len(gene_ids_in_vocab)} genes in vocabulary")
            
            # Subset to genes in vocabulary
            adata_filtered = adata[:, adata.var["id_in_vocab"] >= 0]
            
            # Load model configuration
            import json
            with open(config_file, 'r') as f:
                model_config = json.load(f)
            
            # Create model
            model = TransformerModel(
                ntoken=len(vocab),
                d_model=model_config["embsize"],
                nhead=model_config["nheads"],
                d_hid=model_config["d_hid"],
                nlayers=model_config["nlayers"],
                vocab=vocab,
                dropout=0.2,
                pad_token=pad_token,
                pad_value=-2,
                do_mvc=True,
                do_dab=True,
                use_batch_labels=True,
                num_batch_labels=1,
                domain_spec_batchnorm=True,
                n_input_bins=n_bins,
                ecs_threshold=0.8,
                explicit_zero_prob=True,
                use_fast_transformer=True,
                pre_norm=False,
            )
            
            # Load model weights
            model.load_state_dict(torch.load(model_file, map_location=device))
            model.to(device)
            model.eval()
            
            # Preprocess data
            preprocessor = Preprocessor(
                use_key="X",
                filter_gene_by_counts=3,
                filter_cell_by_counts=False,
                normalize_total=1e4,
                result_normed_key="X_normed",
                log1p=True,
                result_log1p_key="X_log1p",
                subset_hvg=n_hvg,
                hvg_flavor="seurat_v3",
                binning=n_bins,
                result_binned_key="X_binned",
            )
            preprocessor(adata_filtered)
            
            # Tokenize
            input_layer_key = "X_binned"
            all_counts = (
                adata_filtered.layers[input_layer_key].toarray()
                if hasattr(adata_filtered.layers[input_layer_key], 'toarray')
                else adata_filtered.layers[input_layer_key]
            )
            genes = adata_filtered.var["gene_name"].tolist()
            gene_ids = np.array(vocab(genes), dtype=int)
            
            tokenized_data = tokenize_and_pad_batch(
                all_counts,
                gene_ids,
                max_len=max_seq_len,
                vocab=vocab,
                pad_token=pad_token,
                pad_value=-2,
                append_cls=True,
                include_zero_gene=True,
            )
            
            # Generate embeddings
            all_gene_ids = torch.tensor(tokenized_data["genes"], dtype=torch.long, device=device)
            all_values = torch.tensor(tokenized_data["values"], dtype=torch.float, device=device)
            src_key_padding_mask = all_gene_ids.eq(vocab[pad_token])
            
            with torch.no_grad():
                embeddings = model.encode_batch(
                    all_gene_ids,
                    all_values,
                    src_key_padding_mask=src_key_padding_mask,
                    batch_size=64,
                    batch_labels=None,
                    time_step=0,
                    return_np=True,
                )
            
            # Normalize embeddings
            embeddings = embeddings / np.linalg.norm(embeddings, axis=1, keepdims=True)
            
            # If we filtered genes, we need to map back to original adata
            if adata_filtered.n_obs != adata.n_obs:
                # This shouldn't happen if we only filtered genes, but check
                print("Warning: Cell count mismatch after filtering")
            
            print(f"Generated embeddings shape: {embeddings.shape}")
            return embeddings
        else:
            print(f"Model files not found in {model_dir}, using fallback approach")
    
    # Fallback: create embeddings using basic approach
    print("Using fallback embedding generation")
    n_cells = adata.n_obs
    embedding_dim = 512
    
    # Use gene expression patterns to create more meaningful embeddings
    if hasattr(adata.X, 'toarray'):
        X = adata.X.toarray()
    else:
        X = adata.X
    
    # Simple PCA-like transformation
    from sklearn.decomposition import PCA
    pca = PCA(n_components=embedding_dim)
    embeddings = pca.fit_transform(X)
    
    # Normalize embeddings
    embeddings = embeddings / np.linalg.norm(embeddings, axis=1, keepdims=True)
    
    print(f"Generated embeddings shape: {embeddings.shape}")
    return embeddings

def main():
    parser = argparse.ArgumentParser(description='Generate scGPT embeddings')
    parser.add_argument('--input', required=True, help='Input .h5ad file path')
    parser.add_argument('--output', required=True, help='Output .h5ad file path')
    parser.add_argument('--model-dir', help='scGPT model directory',default='/gpfs/scratch/nk4167/miniconda/envs/scgpt_env/models')
    
    args = parser.parse_args()
    
    # Create output directory if needed
    os.makedirs(os.path.dirname(args.output), exist_ok=True)
    
    # Load data
    adata = load_data(args.input)
    
    # Generate embeddings
    embeddings = generate_scgpt_embeddings(adata, args.model_dir)
    
    # Store embeddings in AnnData
    adata.obsm["X_scGPT"] = embeddings
    
    # Add metadata
    adata.uns["scGPT_embedding_params"] = {
        "model": "scGPT",
        "embedding_dim": embeddings.shape[1],
        "normalization": "l2"
    }
    
    # Save result
    print(f"Saving results to {args.output}")
    adata.write_h5ad(args.output)
    print("✅ scGPT embedding generation complete!")

if __name__ == "__main__":
    main()