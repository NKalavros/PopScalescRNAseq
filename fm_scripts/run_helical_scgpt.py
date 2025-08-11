
import argparse
from helical.models.scgpt import scGPT, scGPTConfig # type: ignore
import anndata as ad #type: ignore
import os
import torch #type: ignore

def main():
    parser = argparse.ArgumentParser(description="Run Geneformer embedding generation.")
    parser.add_argument('--workdir', type=str, default='.', help='Working directory to change to')
    parser.add_argument('--batch_size', type=int, default=20, help='Batch size for Geneformer')
    parser.add_argument('--input_file', type=str, default='data.h5ad', required=False, help='Input h5ad file')
    parser.add_argument('--output_file', type=str, default='embeddings/geneformer_scgpt.h5ad', help='Output h5ad file')
    args = parser.parse_args()

    print(f"Arguments: {args}")
    os.chdir(args.workdir)
    model_config = scGPTConfig(batch_size=args.batch_size, device= "cuda" if torch.cuda.is_available() else "cpu")
    scgpt = scGPT(configurer = model_config)

    adata = ad.read_h5ad(args.input_file)
    dataset = scgpt.process_data(adata)
    embeddings = scgpt.get_embeddings(dataset)
    adata.obsm['X_geneformer_scgpt'] = embeddings

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