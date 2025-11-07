#!/bin/bash
# Embedding generation for Lung Atlas datasets
# Uses standard (non-Helical) embedding generation scripts
# Requires A100 GPUs for optimal performance

# Define studies to process (only those with successfully converted h5ad files)
# Failed conversions: Ireland2020, Qian2020
studies=("Bischoff2021" "Chan2021" "Guo2018" "Kim2020" "Laughney2020" "Maynard2020" "Song2019" "Xing2021" "Zilionis2019")

# scGPT embeddings
conda activate /gpfs/scratch/nk4167/miniconda/envs/scgpt_env
for study in "${studies[@]}"; do
    cd /gpfs/scratch/nk4167/LungAtlas/$study
    mkdir -p embeddings
    if [ ! -f "embeddings/scgpt.h5ad" ]; then
        python /gpfs/scratch/nk4167/PopScalescRNAseq/fm_scripts/scgpt_embeddings.py  --input data.h5ad --output embeddings/scgpt.h5ad
    fi
done


# Geneformer embeddings
conda activate /gpfs/scratch/nk4167/miniconda/envs/geneformer_env
for study in "${studies[@]}"; do
    cd /gpfs/scratch/nk4167/LungAtlas/$study
    mkdir -p embeddings
    if [ ! -f "embeddings/geneformer.h5ad" ]; then
        python /gpfs/scratch/nk4167/PopScalescRNAseq/fm_scripts/geneformer_embeddings.py  --input data.h5ad --output embeddings/geneformer.h5ad
    fi
done


# scFoundation embeddings
conda activate /gpfs/scratch/nk4167/miniconda/envs/scFoundation_env
for study in "${studies[@]}"; do
    cd /gpfs/scratch/nk4167/LungAtlas/$study
    mkdir -p embeddings
    if [ ! -f "embeddings/scfoundation.h5ad" ]; then
        python /gpfs/scratch/nk4167/PopScalescRNAseq/fm_scripts/scfoundation_embeddings.py  --input data.h5ad --output embeddings/scfoundation.h5ad
    fi
done


# scimilarity embeddings
conda activate /gpfs/scratch/nk4167/miniconda/envs/scimilarity_env
for study in "${studies[@]}"; do
    cd /gpfs/scratch/nk4167/LungAtlas/$study
    mkdir -p embeddings
    if [ ! -f "embeddings/scimilarity.h5ad" ]; then
        python /gpfs/scratch/nk4167/PopScalescRNAseq/fm_scripts/scimilarity_embeddings.py  --input data.h5ad --output embeddings/scimilarity.h5ad
    fi
done


# UCE embeddings
conda activate /gpfs/scratch/nk4167/miniconda/envs/uce_env
FULL_ENV_PATH='/gpfs/scratch/nk4167/miniconda/envs/uce_env'
for study in "${studies[@]}"; do
    cd /gpfs/scratch/nk4167/LungAtlas/$study
    mkdir -p embeddings
    if [ ! -f "embeddings/uce.h5ad" ]; then
        rm $FULL_ENV_PATH/UCE/*embeddings*
        python /gpfs/scratch/nk4167/PopScalescRNAseq/fm_scripts/uce_embeddings.py  --input $(pwd)/data.h5ad --output embeddings/uce.h5ad --batch_size 256
        mv $FULL_ENV_PATH/UCE/embeddingsdata_uce_adata.h5ad embeddings/uce.h5ad
    fi
done


# TranscriptFormer embeddings
source ~/.bashrc
for study in "${studies[@]}"; do
    cd /gpfs/scratch/nk4167/LungAtlas/$study
    mkdir -p embeddings
    if [ ! -f "embeddings/transcriptformer.h5ad" ]; then
        bash /gpfs/scratch/nk4167/PopScalescRNAseq/fm_scripts/run_transcriptformer.sh data.h5ad
    fi
done
