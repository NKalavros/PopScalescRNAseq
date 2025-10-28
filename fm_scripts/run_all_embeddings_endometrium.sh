#!/bin/bash
# Embedding generation for HECA Endometrium datasets (cells and nuclei)
# Uses standard (non-Helical) embedding generation scripts
# Requires A100 GPUs for optimal performance

# scGPT embeddings
conda activate /gpfs/scratch/nk4167/miniconda/envs/scgpt_env
cd /gpfs/scratch/nk4167/EndometriumAtlas/HECA_cells
mkdir -p embeddings
if [ ! -f "embeddings/scgpt.h5ad" ]; then
    python /gpfs/scratch/nk4167/PopScalescRNAseq/fm_scripts/scgpt_embeddings.py  --input data_cells.h5ad --output embeddings/scgpt.h5ad
fi
cd /gpfs/scratch/nk4167/EndometriumAtlas/HECA_nuclei
mkdir -p embeddings
if [ ! -f "embeddings/scgpt.h5ad" ]; then
    python /gpfs/scratch/nk4167/PopScalescRNAseq/fm_scripts/scgpt_embeddings.py  --input data_nuclei.h5ad --output embeddings/scgpt.h5ad
fi


# Geneformer embeddings
conda activate /gpfs/scratch/nk4167/miniconda/envs/geneformer_env
cd /gpfs/scratch/nk4167/EndometriumAtlas/HECA_cells
mkdir -p embeddings
if [ ! -f "embeddings/geneformer.h5ad" ]; then
    python /gpfs/scratch/nk4167/PopScalescRNAseq/fm_scripts/geneformer_embeddings.py  --input data_cells.h5ad --output embeddings/geneformer.h5ad
fi
cd /gpfs/scratch/nk4167/EndometriumAtlas/HECA_nuclei
mkdir -p embeddings
if [ ! -f "embeddings/geneformer.h5ad" ]; then
    python /gpfs/scratch/nk4167/PopScalescRNAseq/fm_scripts/geneformer_embeddings.py  --input data_nuclei.h5ad --output embeddings/geneformer.h5ad
fi


# scFoundation embeddings
conda activate /gpfs/scratch/nk4167/miniconda/envs/scFoundation_env
cd /gpfs/scratch/nk4167/EndometriumAtlas/HECA_cells
mkdir -p embeddings
if [ ! -f "embeddings/scfoundation.h5ad" ]; then
    python /gpfs/scratch/nk4167/PopScalescRNAseq/fm_scripts/scfoundation_embeddings.py  --input data_cells.h5ad --output embeddings/scfoundation.h5ad
fi
cd /gpfs/scratch/nk4167/EndometriumAtlas/HECA_nuclei
mkdir -p embeddings
if [ ! -f "embeddings/scfoundation.h5ad" ]; then
    python /gpfs/scratch/nk4167/PopScalescRNAseq/fm_scripts/scfoundation_embeddings.py  --input data_nuclei.h5ad --output embeddings/scfoundation.h5ad
fi


# scimilarity embeddings
conda activate /gpfs/scratch/nk4167/miniconda/envs/scimilarity_env
cd /gpfs/scratch/nk4167/EndometriumAtlas/HECA_cells
mkdir -p embeddings
if [ ! -f "embeddings/scimilarity.h5ad" ]; then
    python /gpfs/scratch/nk4167/PopScalescRNAseq/fm_scripts/scimilarity_embeddings.py  --input data_cells.h5ad --output embeddings/scimilarity.h5ad
fi
cd /gpfs/scratch/nk4167/EndometriumAtlas/HECA_nuclei
mkdir -p embeddings
if [ ! -f "embeddings/scimilarity.h5ad" ]; then
    python /gpfs/scratch/nk4167/PopScalescRNAseq/fm_scripts/scimilarity_embeddings.py  --input data_nuclei.h5ad --output embeddings/scimilarity.h5ad
fi


# UCE embeddings
conda activate /gpfs/scratch/nk4167/miniconda/envs/uce_env
FULL_ENV_PATH='/gpfs/scratch/nk4167/miniconda/envs/uce_env'
cd /gpfs/scratch/nk4167/EndometriumAtlas/HECA_cells
mkdir -p embeddings
if [ ! -f "embeddings/uce.h5ad" ]; then
    rm $FULL_ENV_PATH/UCE/*embeddings*
    python /gpfs/scratch/nk4167/PopScalescRNAseq/fm_scripts/uce_embeddings.py  --input $(pwd)/data_cells.h5ad --output embeddings/uce.h5ad --batch_size 256
    mv $FULL_ENV_PATH/UCE/embeddingsdata_cells_uce_adata.h5ad embeddings/uce.h5ad
fi
cd /gpfs/scratch/nk4167/EndometriumAtlas/HECA_nuclei
mkdir -p embeddings
if [ ! -f "embeddings/uce.h5ad" ]; then
    rm $FULL_ENV_PATH/UCE/*embeddings*
    python /gpfs/scratch/nk4167/PopScalescRNAseq/fm_scripts/uce_embeddings.py  --input $(pwd)/data_nuclei.h5ad --output embeddings/uce.h5ad --batch_size 256
    mv $FULL_ENV_PATH/UCE/embeddingsdata_nuclei_uce_adata.h5ad embeddings/uce.h5ad
fi


# TranscriptFormer embeddings
cd /gpfs/scratch/nk4167/EndometriumAtlas/HECA_cells
mkdir -p embeddings
if [ ! -f "embeddings/transcriptformer.h5ad" ]; then
    bash /gpfs/scratch/nk4167/PopScalescRNAseq/fm_scripts/run_transcriptformer.sh
fi
cd /gpfs/scratch/nk4167/EndometriumAtlas/HECA_nuclei
mkdir -p embeddings
if [ ! -f "embeddings/transcriptformer.h5ad" ]; then
    bash /gpfs/scratch/nk4167/PopScalescRNAseq/fm_scripts/run_transcriptformer.sh
fi
