import os
import argparse

def parse_args():
    parser = argparse.ArgumentParser(description='Generate CellFM embeddings for scRNA-seq data.')
    parser.add_argument('--input', type=str, required=True, help='Path to input .h5ad file with raw data.')
    parser.add_argument('--output', type=str, required=True, help='Path to output .h5ad file with embeddings.')
    parser.add_argument('--model-dir', type=str, default="/gpfs/scratch/nk4167/miniconda/envs/uce_env/UCE/", 
                        help='Directory containing the CellFM model repository.')
    parser.add_argument('--full_env_path', type=str, default="/gpfs/scratch/nk4167/miniconda/envs/uce_env/UCE/",
                        help='Full path to the environment where the model is set up.')
    parser.add_argument('--species', type=str, default='human', help='Species for the model (e.g., human, mouse).')
    parser.add_argument('--batch_size', type=int, default=64, help='Batch size for processing.')
    return parser.parse_args()

if __name__ == "__main__":
    args = parse_args()
    FULL_ENV_PATH=args.full_env_path
    INPUT=args.input
    OUTPUT=args.output
    MODEL_DIR=args.model_dir
    SPECIES = args.species
    BATCH_SIZE = args.batch_size
    # Since the output is usually a file, just fix it as that files directory (e.g. scrna_embeddings/uce.h5ad) should become scrna_embeddings
    OUTPUT_DIR = os.path.dirname(OUTPUT)
    #python eval_single_anndata.py --adata_path {path_to_anndata} --dir {output_dir} --species {species} --model_loc {model_loc} --batch_size {batch_size}
    command = f"python {FULL_ENV_PATH}/eval_single_anndata.py --adata_path {INPUT} --dir {OUTPUT_DIR} --species {SPECIES} --model_loc {MODEL_DIR} --batch_size {BATCH_SIZE}"
    print(f"Running command: {command}")
    os.system(command)
    print(f"Embeddings saved to {OUTPUT}")
    # Note: The actual embedding extraction logic is handled in eval_single_anndata.py