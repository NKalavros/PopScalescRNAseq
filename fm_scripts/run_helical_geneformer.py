
import argparse
from helical.models.geneformer import Geneformer, GeneformerConfig #type: ignore   
import anndata as ad #type: ignore
import os
import torch #type: ignore

def main():
    parser = argparse.ArgumentParser(description="Run Geneformer embedding generation.")
    parser.add_argument('--workdir', type=str, default='.', help='Working directory to change to')
    parser.add_argument('--batch_size', type=int, default=20, help='Batch size for Geneformer')
    parser.add_argument('--input_file', type=str, required=True, help='Input h5ad file')
    args = parser.parse_args()

    print(f"Arguments: {args}")
    os.chdir(args.workdir)
    model_config = GeneformerConfig(batch_size=args.batch_size, device= "cuda" if torch.cuda.is_available() else "cpu")
    geneformer_v2 = Geneformer(model_config)

    ann_data = ad.read_h5ad(args.input_file)
    dataset = geneformer_v2.process_data(ann_data)
    embeddings = geneformer_v2.get_embeddings(dataset)
    print("Base model embeddings shape:", embeddings.shape)

if __name__ == "__main__":
    main()