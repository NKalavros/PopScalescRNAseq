#!/bin/bash
# Embedding generation for Lung Atlas datasets
# Uses Helical-wrapped models for optimized A100 GPU performance

# Define studies to process (only those with successfully converted h5ad files)
# Failed conversions: Ireland2020, Qian2020
studies=("Bischoff2021" "Chan2021" "Guo2018" "Kim2020" "Laughney2020" "Maynard2020" "Song2019" "Xing2021" "Zilionis2019")

# Geneformer through Helical (Lung Atlas)
conda activate /gpfs/scratch/nk4167/miniconda/envs/helical-package
for study in "${studies[@]}"; do
    cd /gpfs/scratch/nk4167/LungAtlas/$study
    mkdir -p embeddings
    if [ ! -f "embeddings/geneformer_helical.h5ad" ]; then
        python /gpfs/scratch/nk4167/PopScalescRNAseq/fm_scripts/run_helical_geneformer.py  --input data.h5ad --output embeddings/geneformer_helical.h5ad
    fi
done


# scGPT through Helical (Lung Atlas)
conda activate /gpfs/scratch/nk4167/miniconda/envs/helical-package
for study in "${studies[@]}"; do
    cd /gpfs/scratch/nk4167/LungAtlas/$study
    mkdir -p embeddings
    if [ ! -f "embeddings/scgpt_helical.h5ad" ]; then
        python /gpfs/scratch/nk4167/PopScalescRNAseq/fm_scripts/run_helical_scgpt.py  --input data.h5ad --output embeddings/scgpt_helical.h5ad
    fi
done


# UCE through Helical (Lung Atlas)
conda activate /gpfs/scratch/nk4167/miniconda/envs/helical-package
for study in "${studies[@]}"; do
    cd /gpfs/scratch/nk4167/LungAtlas/$study
    mkdir -p embeddings
    if [ ! -f "embeddings/uce_helical.h5ad" ]; then
        python /gpfs/scratch/nk4167/PopScalescRNAseq/fm_scripts/run_helical_uce.py  --input data.h5ad --output embeddings/uce_helical.h5ad
    fi
done
