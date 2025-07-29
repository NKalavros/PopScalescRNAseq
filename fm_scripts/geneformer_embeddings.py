#!/usr/bin/env python
"""
Geneformer Embedding Generation Script
Generates embeddings for scRNA-seq data using Geneformer foundation model
"""

import os
import sys
import torch
import numpy as np
import argparse
import anndata as ad
from pathlib import Path
import scanpy as sc
from geneformer import TranscriptomeTokenizer
import mygene
from geneformer import EmbExtractor
import huggingface_hub

# Add flash attention as default for torch


def parse_args():
    parser = argparse.ArgumentParser(description="Generate Geneformer embeddings for scRNA-seq data.")
    parser.add_argument('--input', type=str, required=True, help='Path to input h5ad file with raw data.',default='pbmc_scanpy_data.h5ad')
    parser.add_argument('--output', type=str, required=True, help='Path to output h5ad file with embeddings.',default='pbmc3k_embeddings.h5ad')
    parser.add_argument('--model_version', type=str, default="V2", choices=["V1", "V2"], help='Version of Geneformer model to use.')
    return parser.parse_args()
# Check if in ipython
if 'ipython' in sys.argv[0] or 'jupyter' in sys.argv[0]:
    print("Running in Jupyter Notebook or IPython environment.")
    # Get default args (parse_args fails in ipython)
    args = {
        "input": "pbmc_scanpy_data.h5ad",
        "output": "pbmc3k_embeddings.h5ad",
        "model_version": "V2"
    }
    input = args['input']
    output = args['output']
    model_version = args['model_version']
    print(f"Using input: {input}, output: {output}, model_version: {model_version}")
else:
    print("Running in standard Python environment.")
    args = parse_args()
    input = args.input
    output = args.output
    model_version = args.model_version
    print(f"Using input: {input}, output: {output}, model_version: {model_version}")
# read in anndata file
if 'pbmc' in input:
    print("Reading PBMC data...")
    adata = sc.datasets.pbmc3k()
else:
    print(f"Reading data from {input}...")
    adata = sc.read_h5ad(input)
# Convert gene names to ensembl
adata.var_names_make_unique()
sc.pp.filter_genes(adata, min_cells=10)  # Filter genes expressed in at least 3 cells
# Convert gene names to Ensembl IDs
mg = mygene.MyGeneInfo()

# Get gene symbols from the current var_names
gene_symbols = adata.var_names.tolist()

# Query mygene to get Ensembl IDs
result = mg.querymany(gene_symbols, scopes='symbol', fields='ensembl.gene', species='human')

# Create a mapping dictionary
ensembl_mapping = {}
for item in result:
    if 'ensembl' in item and 'gene' in item['ensembl']:
        ensembl_id = item['ensembl']['gene']
        if isinstance(ensembl_id, list):
            ensembl_id = ensembl_id[0]  # Take first if multiple
        ensembl_mapping[item['query']] = ensembl_id
    else:
        ensembl_mapping[item['query']] = None

# Add ensembl_id column to adata.var
adata.var['ensembl_id'] = [ensembl_mapping.get(gene, None) for gene in adata.var_names]

# Filter out genes without valid Ensembl IDs
valid_mask = adata.var['ensembl_id'].notna() & (adata.var['ensembl_id'] != "")
adata = adata[:, valid_mask].copy()
adata.var = adata.var.loc[valid_mask].copy()
adata.var_names = adata.var['ensembl_id']
adata.var_names_make_unique()  # Ensure unique names after conversion

print(f"Genes with Ensembl IDs: {adata.n_vars}/{len(gene_symbols)}")
# Add n_counts to adata.obs
adata.obs['n_counts'] = adata.X.sum(axis=1)
adata.var_names_make_unique()  # Ensure unique names after conversion

# Ensure all columns in adata.var and adata.obs are string type (per column)
for col in adata.var.columns:
    adata.var[col] = adata.var[col].astype(str)
for col in adata.obs.columns:
    adata.obs[col] = adata.obs[col].astype(str)

# Make a directory for the embeddings if it doesn't exist
os.makedirs("geneformer_tokenizer", exist_ok=True)
print('Tokenizing data...')
print('Data shape:', adata.shape)
# Write file there
adata.write_h5ad("geneformer_tokenizer/geneformer_tokenized.h5ad")
# Save the new anndata there
tk = TranscriptomeTokenizer(nproc=int(os.environ.get('SLURM_CPUS_PER_TASK', '1')), model_version=model_version)
tk.tokenize_data("geneformer_tokenizer/", 
                 "geneformer_tokenizer", 
                 "geneformer_tokenized", 
                 file_format="h5ad")
# Initialize an emb_extractor


torch.cuda.empty_cache()

n_cells = adata.n_obs
# 0 for last layer, -1 for second to last
layer = -1
model_version='V2' # Hardcoding this
# initiate EmbExtractor
embex = EmbExtractor(model_type="Pretrained",
                     num_classes=0,
                     max_ncells=None,
                     emb_mode='cell',
                     emb_layer=layer,
                     forward_batch_size=20,
                     nproc=int(os.environ.get('SLURM_CPUS_PER_TASK', '4')),
                      )
huggingface_hub.login(token=os.environ.get('HUGGINGFACE_TOKEN', ''))

# Set model path based on version
if model_version == "V2":
    model_path = "ctheodoris/Geneformer"
elif model_version == "V2_small":
    model_path_huggingface="https://huggingface.co/ctheodoris/Geneformer/tree/main/Geneformer-V2-104M"
    # Download the model from Hugging Face
    model_path = huggingface_hub.snapshot_download(repo_id=model_path_huggingface,
                                                  local_dir="/gpfs/scratch/nk4167/miniconda/envs/geneformer_env/",
                                                  local_dir_use_symlinks=False)
else:
    model_path = "ctheodoris/Geneformer-V1"

# extracts embedding from input data
# input data is tokenized rank value encodings generated by Geneformer tokenizer (see tokenizing_scRNAseq_data.ipynb)
embs = embex.extract_embs(model_path,
                          "geneformer_tokenizer/geneformer_tokenized.dataset",
                          ".",
                          "geneformer_embedded")