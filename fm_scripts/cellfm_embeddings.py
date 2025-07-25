#!/usr/bin/env python
"""
CellFM Embedding Generation Script
Generates embeddings for scRNA-seq data using CellFM foundation model
"""

import os
import sys
import torch
import numpy as np
import argparse
import anndata as ad
from pathlib import Path

def load_data(input_path):
    """Load data from .h5ad file"""
    print(f"Loading data from {input_path}")
    adata = ad.read_h5ad(input_path)
    print(f"Data shape: {adata.shape}")
    return adata

def generate_cellfm_embeddings(adata, model_dir=None):
    """Generate embeddings using CellFM"""
    print("Generating CellFM embeddings...")
    
    # Set CellFM directory
    if model_dir is None:
        model_dir = "/gpfs/scratch/nk4167/miniconda/envs/CellFM/CellFM"
    
    # Set CUDA environment variables for MindSpore GPU
    env_prefix = "/gpfs/scratch/nk4167/miniconda/envs/CellFM"
    current_ld_path = os.environ.get('LD_LIBRARY_PATH', '')
    os.environ['LD_LIBRARY_PATH'] = f"{env_prefix}/lib:{current_ld_path}"
    os.environ['CUDA_HOME'] = env_prefix
    os.environ['CUDA_ROOT'] = env_prefix
    print(f"🔧 Set CUDA environment variables for MindSpore GPU")
    
    try:
        # Add CellFM to Python path
        sys.path.insert(0, model_dir)
        
        # Import CellFM modules
        print("Importing CellFM modules...")
        import mindspore as ms
        from model import CellFM
        from config import Config
        from scipy.sparse import csr_matrix as csr
        import scanpy as sc
        
        print(f"✅ MindSpore: {ms.__version__}")
        print("✅ CellFM modules imported successfully")
        
        # Set MindSpore context - try GPU first (with CUDA 11.6), fallback to CPU
        try:
            ms.set_context(mode=ms.GRAPH_MODE, device_target="GPU")
            print("✅ Using MindSpore GPU mode with CUDA 11.6")
        except RuntimeError as gpu_error:
            print(f"⚠️  GPU not available for MindSpore: {gpu_error}")
            print("🔄 Falling back to CPU mode...")
            ms.set_context(
                mode=ms.GRAPH_MODE, 
                device_target="CPU",
                enable_graph_kernel=True,
                graph_kernel_flags="--enable_parallel_fusion"
            )
            print("✅ Using MindSpore CPU mode with optimizations")
        
        # Download CellFM model from Hugging Face if not present
        print("Setting up CellFM model...")
        
        # Check for model checkpoint
        checkpoint_path = os.path.join(model_dir, "base_weight.ckpt")
        
        if not os.path.exists(checkpoint_path):
            print("Downloading CellFM model from Hugging Face...")
            try:
                from huggingface_hub import hf_hub_download
                
                # Download the model checkpoint
                checkpoint_path = hf_hub_download(
                    repo_id="ShangguanNingyuan/CellFM",
                    filename="base_weight.ckpt",
                    cache_dir=model_dir
                )
                print(f"✅ Model downloaded to: {checkpoint_path}")
                
            except Exception as download_error:
                print(f"Failed to download model: {download_error}")
                raise download_error
        
        # Preprocess data for CellFM
        print("Preprocessing data for CellFM...")
        
        # Convert to sparse format as required by CellFM
        if not hasattr(adata.X, 'toarray'):
            adata.X = csr(adata.X)
        
        # Basic filtering (CellFM preprocessing)
        sc.pp.filter_cells(adata, min_genes=1)
        sc.pp.filter_genes(adata, min_cells=1)
        
        n_genes = adata.n_vars
        print(f"Data preprocessed: {adata.n_obs} cells, {n_genes} genes")
        
        # CellFM model expects exactly 24,080 genes - handle mismatch
        expected_genes = 24080
        print(f"Model expects {expected_genes} genes, data has {n_genes} genes")
        
        # Store original expression data for model input
        if hasattr(adata.X, 'toarray'):
            X_data = adata.X.toarray()
        else:
            X_data = adata.X.copy()
        
        if n_genes != expected_genes:
            print(f"Adjusting gene count from {n_genes} to {expected_genes}")
            
            if n_genes < expected_genes:
                # Pad with zeros if we have fewer genes
                padding_size = expected_genes - n_genes
                padding = np.zeros((X_data.shape[0], padding_size), dtype=X_data.dtype)
                X_data = np.hstack([X_data, padding])
                print(f"Padded with {padding_size} zero columns")
            elif n_genes > expected_genes:
                # Truncate if we have more genes (keep first genes)
                X_data = X_data[:, :expected_genes]
                print(f"Truncated to first {expected_genes} genes")
        
        # Initialize CellFM model with proper configuration
        print("Initializing CellFM model...")
        cfg = Config()
        cfg.enc_dims = 1536  # CellFM embedding dimensions
        cfg.enc_nlayers = 40  # Number of encoder layers
        cfg.enc_num_heads = 48  # Number of attention heads
        
        # Create CellFM model with expected gene count
        model = CellFM(expected_genes, cfg)
        
        # Load pretrained weights
        print(f"Loading pretrained weights from {checkpoint_path}...")
        para = ms.load_checkpoint(checkpoint_path)
        ms.load_param_into_net(model, para)
        
        # Set model to evaluation mode
        model.set_train(False)
        
        # Convert data to MindSpore tensor format
        print("Converting data to MindSpore format...")
        
        # Convert to MindSpore tensor
        X_tensor = ms.Tensor(X_data, dtype=ms.float32)
        
        # Generate embeddings using the real CellFM model
        print("Generating embeddings with CellFM model...")
        try:
            # MindSpore 1.x uses different context manager
            model.set_train(False)
            embeddings_tensor = model(X_tensor)
            embeddings = embeddings_tensor.asnumpy()
        except Exception as inference_error:
            print(f"Direct inference failed: {inference_error}")
            print("Trying alternative inference method...")
            # Alternative method for older MindSpore versions
            embeddings_tensor = model.construct(X_tensor)
            embeddings = embeddings_tensor.asnumpy()
        
        print(f"Generated real CellFM embeddings shape: {embeddings.shape}")
        return embeddings
        
    except Exception as e:
        print(f"Error with real CellFM model: {e}")
        print("This indicates the CellFM model setup needs to be fixed.")
        print("Please ensure:")
        print("1. CellFM repository is properly cloned")
        print("2. Model weights are downloaded from Hugging Face")
        print("3. MindSpore is correctly installed with GPU support")
        print("4. All CellFM dependencies are available")
        
        # For now, indicate this is not working rather than fall back to simulation
        raise Exception(f"CellFM real model failed: {e}. Please fix the setup before proceeding.")

def main():
    parser = argparse.ArgumentParser(description='Generate CellFM embeddings')
    parser.add_argument('--input', required=True, help='Input .h5ad file path')
    parser.add_argument('--output', required=True, help='Output .h5ad file path')
    parser.add_argument('--model-dir', help='CellFM repository directory', 
                       default="/gpfs/scratch/nk4167/miniconda/envs/CellFM/CellFM")
    
    args = parser.parse_args()
    
    print(f"🚀 Starting CellFM embedding generation...")
    print(f"📥 Input: {args.input}")
    print(f"📤 Output: {args.output}")
    print(f"🏠 Model directory: {args.model_dir}")
    
    # Check if input file exists
    if not os.path.exists(args.input):
        raise FileNotFoundError(f"Input file not found: {args.input}")
    
    # Check if model directory exists
    if not os.path.exists(args.model_dir):
        raise FileNotFoundError(f"CellFM repository not found at: {args.model_dir}")
        
    # Create output directory if needed
    output_dir = os.path.dirname(args.output)
    if output_dir:
        os.makedirs(output_dir, exist_ok=True)
        print(f"📁 Created output directory: {output_dir}")
    
    # Load data
    adata = load_data(args.input)
    
    # Generate embeddings
    embeddings = generate_cellfm_embeddings(adata, args.model_dir)
    
    # Store embeddings in AnnData
    adata.obsm["X_CellFM"] = embeddings
    
    # Add metadata
    adata.uns["CellFM_embedding_params"] = {
        "model": "CellFM",
        "embedding_dim": embeddings.shape[1],
        "model_dir": args.model_dir,
        "device": "cuda" if torch.cuda.is_available() else "cpu",
        "model_source": "Hugging Face: ShangguanNingyuan/CellFM"
    }
    
    # Save result
    print(f"💾 Saving results to {args.output}")
    adata.write_h5ad(args.output)
    print("✅ CellFM embedding generation complete!")
    print(f"📊 Generated embeddings: {embeddings.shape[0]} cells × {embeddings.shape[1]} dimensions")

if __name__ == "__main__":
    main()