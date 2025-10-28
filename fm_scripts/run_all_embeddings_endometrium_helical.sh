#!/bin/bash
# Embedding generation for HECA Endometrium datasets (cells and nuclei)
# Uses Helical-wrapped models for optimized A100 GPU performance

# Geneformer through Helical (HECA Endometrium)
conda activate /gpfs/scratch/nk4167/miniconda/envs/helical-package
cd /gpfs/scratch/nk4167/EndometriumAtlas/HECA_cells
mkdir -p embeddings
if [ ! -f "embeddings/geneformer_helical.h5ad" ]; then
    python /gpfs/scratch/nk4167/PopScalescRNAseq/fm_scripts/run_helical_geneformer.py  --input data_cells.h5ad --output embeddings/geneformer_helical.h5ad
fi
cd /gpfs/scratch/nk4167/EndometriumAtlas/HECA_nuclei
mkdir -p embeddings
if [ ! -f "embeddings/geneformer_helical.h5ad" ]; then
    python /gpfs/scratch/nk4167/PopScalescRNAseq/fm_scripts/run_helical_geneformer.py  --input data_nuclei.h5ad --output embeddings/geneformer_helical.h5ad
fi


# scGPT through Helical (HECA Endometrium)
conda activate /gpfs/scratch/nk4167/miniconda/envs/helical-package
cd /gpfs/scratch/nk4167/EndometriumAtlas/HECA_cells
mkdir -p embeddings
if [ ! -f "embeddings/scgpt_helical.h5ad" ]; then
    python /gpfs/scratch/nk4167/PopScalescRNAseq/fm_scripts/run_helical_scgpt.py  --input data_cells.h5ad --output embeddings/scgpt_helical.h5ad
fi
cd /gpfs/scratch/nk4167/EndometriumAtlas/HECA_nuclei
mkdir -p embeddings
if [ ! -f "embeddings/scgpt_helical.h5ad" ]; then
    python /gpfs/scratch/nk4167/PopScalescRNAseq/fm_scripts/run_helical_scgpt.py  --input data_nuclei.h5ad --output embeddings/scgpt_helical.h5ad
fi


# UCE through Helical (HECA Endometrium)
conda activate /gpfs/scratch/nk4167/miniconda/envs/helical-package
cd /gpfs/scratch/nk4167/EndometriumAtlas/HECA_cells
mkdir -p embeddings
if [ ! -f "embeddings/uce_helical.h5ad" ]; then
    python /gpfs/scratch/nk4167/PopScalescRNAseq/fm_scripts/run_helical_uce.py  --input data_cells.h5ad --output embeddings/uce_helical.h5ad
fi
cd /gpfs/scratch/nk4167/EndometriumAtlas/HECA_nuclei
mkdir -p embeddings
if [ ! -f "embeddings/uce_helical.h5ad" ]; then
    python /gpfs/scratch/nk4167/PopScalescRNAseq/fm_scripts/run_helical_uce.py  --input data_nuclei.h5ad --output embeddings/uce_helical.h5ad
fi
