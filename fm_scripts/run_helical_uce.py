
import argparse
from helical.models.uce import UCE, UCEConfig # type: ignore
import anndata as ad #type: ignore
import os
import numpy as np # type: ignore
import torch #type: ignore

def main():
    parser = argparse.ArgumentParser(description="Run Geneformer embedding generation.")
    parser.add_argument('--workdir', type=str, default='.', help='Working directory to change to')
    parser.add_argument('--batch_size', type=int, default=20, help='Batch size for Geneformer')
    parser.add_argument('--input_file', type=str, default='data.h5ad', required=False, help='Input h5ad file')
    parser.add_argument('--output_file', type=str, default='embeddings/uce_helical.h5ad', help='Output h5ad file')
    parser.add_argument('--model_name', type=str, default='33l_8ep_1024t_1280', help='Model name')
    args = parser.parse_args()

    print(f"Arguments: {args}")
    os.chdir(args.workdir)
    model_config = UCEConfig(model_name= args.model_name,batch_size=args.batch_size, device= "cuda" if torch.cuda.is_available() else "cpu")
    uce = UCE(configurer = model_config)

    adata = ad.read_h5ad(args.input_file)
    adata.var_names_make_unique()
    batch_size = args.batch_size
    # Batched embedding generation for lower RAM usage
    # Initialize a list to store embeddings from each batch
    all_embeddings = []
    # If batch size is 0, just do all the dataset
    if batch_size <= 0 or batch_size > adata.shape[0]:
        batch_size = adata.shape[0]
    print(f"Using batch size: {batch_size}")
    # Iterate over the data in batches
    for start in range(0, adata.shape[0], batch_size):
        end = min(start + batch_size, adata.shape[0])
        ann_data_batch = adata[start:end].to_memory()

        dataset_batch = uce.process_data(ann_data_batch)
        embeddings_batch = uce.get_embeddings(dataset_batch)

        all_embeddings.append(embeddings_batch)

    # Concatenate the embeddings from each batch
    all_embeddings = np.concatenate(all_embeddings, axis=0)
    #dataset = uce.process_data(adata)
    #embeddings = uce.get_embeddings(dataset)
    adata.obsm['X_geneformer_uce'] = all_embeddings

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

    print("Base model embeddings shape:", all_embeddings.shape)

if __name__ == "__main__":
    main()