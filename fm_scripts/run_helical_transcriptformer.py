

from sklearn.neighbors import KNeighborsClassifier
from sklearn.metrics import classification_report,confusion_matrix,ConfusionMatrixDisplay
import warnings
import seaborn as sns
import matplotlib.pyplot as plt
import umap
import pandas as pd
import numpy as np
import logging
import torch
from helical.utils import get_anndata_from_hf_dataset
from datasets import load_dataset


import argparse
from helical.models.transcriptformer import TranscriptFormer, TranscriptFormerConfig
import anndata as ad #type: ignore
import os
import torch #type: ignore

def main():
    parser = argparse.ArgumentParser(description="Run Geneformer embedding generation.")
    parser.add_argument('--workdir', type=str, default='.', help='Working directory to change to')
    parser.add_argument('--batch_size', type=int, default=20, help='Batch size for Geneformer')
    parser.add_argument('--input_file', type=str, default='data.h5ad', required=False, help='Input h5ad file')
    parser.add_argument('--output_file', type=str, default='embeddings/geneformer_helical.h5ad', help='Output h5ad file')
    parser.add_argument('--model_name', type=str, default='gf-18L-316M-i4096',help='Name of the model to use')
    args = parser.parse_args()

    print(f"Arguments: {args}")
    os.chdir(args.workdir)
    model_config = TranscriptFormerConfig(model_name="tf_sapiens", batch_size=args.batch_size)
    transcriptformer = TranscriptFormer(configurer=model_config)

    adata = ad.read_h5ad(args.input_file)
    dataset = transcriptformer.process_data(adata)
    embeddings = transcriptformer.get_embeddings(dataset)
    adata.obsm['X_transcriptformer_helical'] = embeddings

    # Sanitize obs and var columns
    adata.obs.columns = adata.obs.columns.astype(str)
    adata.var.columns = adata.var.columns.astype(str)
    # Convert object dtype columns to string
    for col in adata.obs.columns:
        if adata.obs[col].dtype == 'object':
            adata.obs[col] = adata.obs[col].astype(str)
    for col in adata.var.columns:
        if adata.var[col].dtype == 'object':
            adata.var[col] = adata.var[col].astype(str)
    # Remove problematic _index column if present
    if '_index' in adata.obs.columns:
        adata.obs.drop('_index', axis=1, inplace=True)
    if '_index' in adata.var.columns:
        adata.var.drop('_index', axis=1, inplace=True)
    # Drop adata.var index
    adata.var.reset_index(drop=True, inplace=True)

    adata.write(args.output_file)

    print("Base model embeddings shape:", embeddings.shape)

if __name__ == "__main__":
    main()
