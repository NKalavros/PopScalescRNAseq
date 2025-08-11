# NOTE: This code is slightly outdates because I moved the lake files to two actual directories (lake_scrna and lake_snrna)
conda activate /gpfs/scratch/nk4167/miniconda/envs/scgpt_env
cd /gpfs/scratch/nk4167/BreastAtlas/Bassez
if [ ! -f "embeddings/scgpt.h5ad" ]; then
    python /gpfs/scratch/nk4167/PopScalescRNAseq/fm_scripts/scgpt_embeddings.py  --input data.h5ad --output embeddings/scgpt.h5ad
fi
cd /gpfs/scratch/nk4167/BreastAtlas/Wu
if [ ! -f "embeddings/scgpt.h5ad" ]; then
    python /gpfs/scratch/nk4167/PopScalescRNAseq/fm_scripts/scgpt_embeddings.py  --input data.h5ad --output embeddings/scgpt.h5ad
fi
cd /gpfs/scratch/nk4167/BreastAtlas/Pal
if [ ! -f "embeddings/scgpt.h5ad" ]; then
    python /gpfs/scratch/nk4167/PopScalescRNAseq/fm_scripts/scgpt_embeddings.py  --input data.h5ad --output embeddings/scgpt.h5ad
fi
cd /gpfs/scratch/nk4167/BreastAtlas/Gray

if [ ! -f "embeddings/scgpt.h5ad" ]; then
    python /gpfs/scratch/nk4167/PopScalescRNAseq/fm_scripts/scgpt_embeddings.py  --input data.h5ad --output embeddings/scgpt.h5ad
fi
cd /gpfs/scratch/nk4167/BreastAtlas/Klughammer_snRNA

if [ ! -f "embeddings/scgpt.h5ad" ]; then
    python /gpfs/scratch/nk4167/PopScalescRNAseq/fm_scripts/scgpt_embeddings.py  --input data.h5ad --output embeddings/scgpt.h5ad
fi
cd /gpfs/scratch/nk4167/BreastAtlas/Klughammer_scRNA
if [ ! -f "embeddings/scgpt.h5ad" ]; then
    python /gpfs/scratch/nk4167/PopScalescRNAseq/fm_scripts/scgpt_embeddings.py  --input data.h5ad --output embeddings/scgpt.h5ad
fi


# Now let's do some more models - Geneformer

# Now let's do some more models - Geneformer (BreastAtlas)
conda activate /gpfs/scratch/nk4167/miniconda/envs/geneformer_env
cd /gpfs/scratch/nk4167/BreastAtlas/Bassez
mkdir embeddings
if [ ! -f "embeddings/geneformer.h5ad" ]; then
    python /gpfs/scratch/nk4167/PopScalescRNAseq/fm_scripts/geneformer_embeddings.py  --input data.h5ad --output embeddings/geneformer.h5ad
fi
cd /gpfs/scratch/nk4167/BreastAtlas/Wu
mkdir embeddings
if [ ! -f "embeddings/geneformer.h5ad" ]; then
    python /gpfs/scratch/nk4167/PopScalescRNAseq/fm_scripts/geneformer_embeddings.py  --input data.h5ad --output embeddings/geneformer.h5ad
fi
cd /gpfs/scratch/nk4167/BreastAtlas/Pal
mkdir embeddings
if [ ! -f "embeddings/geneformer.h5ad" ]; then
    python /gpfs/scratch/nk4167/PopScalescRNAseq/fm_scripts/geneformer_embeddings.py  --input data.h5ad --output embeddings/geneformer.h5ad
fi
cd /gpfs/scratch/nk4167/BreastAtlas/Gray
mkdir embeddings
if [ ! -f "embeddings/geneformer.h5ad" ]; then
    python /gpfs/scratch/nk4167/PopScalescRNAseq/fm_scripts/geneformer_embeddings.py  --input data.h5ad --output embeddings/geneformer.h5ad
fi
cd /gpfs/scratch/nk4167/BreastAtlas/Klughammer_snRNA
mkdir embeddings
if [ ! -f "embeddings/geneformer.h5ad" ]; then
    python /gpfs/scratch/nk4167/PopScalescRNAseq/fm_scripts/geneformer_embeddings.py  --input data.h5ad --output embeddings/geneformer.h5ad
fi
cd /gpfs/scratch/nk4167/BreastAtlas/Klughammer_scRNA
mkdir embeddings
if [ ! -f "embeddings/geneformer.h5ad" ]; then
    python /gpfs/scratch/nk4167/PopScalescRNAseq/fm_scripts/geneformer_embeddings.py  --input data.h5ad --output embeddings/geneformer.h5ad
fi


# Now let's do some more models - scFoundation

# Now let's do some more models - scFoundation (BreastAtlas)
conda activate /gpfs/scratch/nk4167/miniconda/envs/scFoundation_env
cd /gpfs/scratch/nk4167/BreastAtlas/Bassez
mkdir embeddings
if [ ! -f "embeddings/scfoundation.h5ad" ]; then
    python /gpfs/scratch/nk4167/PopScalescRNAseq/fm_scripts/scfoundation_embeddings.py  --input data.h5ad --output embeddings/scfoundation.h5ad
fi
cd /gpfs/scratch/nk4167/BreastAtlas/Wu
mkdir embeddings
if [ ! -f "embeddings/scfoundation.h5ad" ]; then
    python /gpfs/scratch/nk4167/PopScalescRNAseq/fm_scripts/scfoundation_embeddings.py  --input data.h5ad --output embeddings/scfoundation.h5ad
fi
cd /gpfs/scratch/nk4167/BreastAtlas/Pal
mkdir embeddings
if [ ! -f "embeddings/scfoundation.h5ad" ]; then
    python /gpfs/scratch/nk4167/PopScalescRNAseq/fm_scripts/scfoundation_embeddings.py  --input data.h5ad --output embeddings/scfoundation.h5ad
fi
cd /gpfs/scratch/nk4167/BreastAtlas/Gray
mkdir embeddings
if [ ! -f "embeddings/scfoundation.h5ad" ]; then
    python /gpfs/scratch/nk4167/PopScalescRNAseq/fm_scripts/scfoundation_embeddings.py  --input data.h5ad --output embeddings/scfoundation.h5ad
fi
cd /gpfs/scratch/nk4167/BreastAtlas/Klughammer_snRNA
mkdir embeddings
if [ ! -f "embeddings/scfoundation.h5ad" ]; then
    python /gpfs/scratch/nk4167/PopScalescRNAseq/fm_scripts/scfoundation_embeddings.py  --input data.h5ad --output embeddings/scfoundation.h5ad
fi
cd /gpfs/scratch/nk4167/BreastAtlas/Klughammer_scRNA
mkdir embeddings
if [ ! -f "embeddings/scfoundation.h5ad" ]; then
    python /gpfs/scratch/nk4167/PopScalescRNAseq/fm_scripts/scfoundation_embeddings.py  --input data.h5ad --output embeddings/scfoundation.h5ad
fi

# Now let's do some more models - SCIMILARITY

# Now let's do some more models - SCIMILARITY (BreastAtlas)
conda activate /gpfs/scratch/nk4167/miniconda/envs/scimilarity_env
cd /gpfs/scratch/nk4167/BreastAtlas/Bassez
mkdir embeddings
if [ ! -f "embeddings/scimilarity.h5ad" ]; then
    python /gpfs/scratch/nk4167/PopScalescRNAseq/fm_scripts/scimilarity_embeddings.py  --input data.h5ad --output embeddings/scimilarity.h5ad
fi
cd /gpfs/scratch/nk4167/BreastAtlas/Wu
mkdir embeddings
if [ ! -f "embeddings/scimilarity.h5ad" ]; then
    python /gpfs/scratch/nk4167/PopScalescRNAseq/fm_scripts/scimilarity_embeddings.py  --input data.h5ad --output embeddings/scimilarity.h5ad
fi
cd /gpfs/scratch/nk4167/BreastAtlas/Pal
mkdir embeddings
if [ ! -f "embeddings/scimilarity.h5ad" ]; then
    python /gpfs/scratch/nk4167/PopScalescRNAseq/fm_scripts/scimilarity_embeddings.py  --input data.h5ad --output embeddings/scimilarity.h5ad
fi
cd /gpfs/scratch/nk4167/BreastAtlas/Gray
mkdir embeddings
if [ ! -f "embeddings/scimilarity.h5ad" ]; then
    python /gpfs/scratch/nk4167/PopScalescRNAseq/fm_scripts/scimilarity_embeddings.py  --input data.h5ad --output embeddings/scimilarity.h5ad
fi
cd /gpfs/scratch/nk4167/BreastAtlas/Klughammer_snRNA
mkdir embeddings
if [ ! -f "embeddings/scimilarity.h5ad" ]; then
    python /gpfs/scratch/nk4167/PopScalescRNAseq/fm_scripts/scimilarity_embeddings.py  --input data.h5ad --output embeddings/scimilarity.h5ad
fi
cd /gpfs/scratch/nk4167/BreastAtlas/Klughammer_scRNA
mkdir embeddings
if [ ! -f "embeddings/scimilarity.h5ad" ]; then
    python /gpfs/scratch/nk4167/PopScalescRNAseq/fm_scripts/scimilarity_embeddings.py  --input data.h5ad --output embeddings/scimilarity.h5ad
fi


# Now let's do some more models - UCE

# Now let's do some more models - UCE (BreastAtlas)
conda activate /gpfs/scratch/nk4167/miniconda/envs/uce_env
FULL_ENV_PATH='/gpfs/scratch/nk4167/miniconda/envs/uce_env'
cd /gpfs/scratch/nk4167/BreastAtlas/Bassez
mkdir embeddings
if [ ! -f "embeddings/uce.h5ad" ]; then
    rm $FULL_ENV_PATH/UCE/*embeddings*
    python /gpfs/scratch/nk4167/PopScalescRNAseq/fm_scripts/uce_embeddings.py  --input $(pwd)/data.h5ad --output embeddings/uce.h5ad --batch_size 256
    mv $FULL_ENV_PATH/UCE/embeddingsdata_uce_adata.h5ad embeddings/uce.h5ad
fi
cd /gpfs/scratch/nk4167/BreastAtlas/Wu
mkdir embeddings
if [ ! -f "embeddings/uce.h5ad" ]; then
    rm $FULL_ENV_PATH/UCE/*embeddings*
    python /gpfs/scratch/nk4167/PopScalescRNAseq/fm_scripts/uce_embeddings.py  --input $(pwd)/data.h5ad --output embeddings/uce.h5ad --batch_size 256
    mv $FULL_ENV_PATH/UCE/embeddingsdata_uce_adata.h5ad embeddings/uce.h5ad
fi
cd /gpfs/scratch/nk4167/BreastAtlas/Pal
mkdir embeddings
if [ ! -f "embeddings/uce.h5ad" ]; then
    rm $FULL_ENV_PATH/UCE/*embeddings*
    python /gpfs/scratch/nk4167/PopScalescRNAseq/fm_scripts/uce_embeddings.py  --input $(pwd)/data.h5ad --output embeddings/uce.h5ad --batch_size 256
    mv $FULL_ENV_PATH/UCE/embeddingsdata_uce_adata.h5ad embeddings/uce.h5ad
fi
cd /gpfs/scratch/nk4167/PopScalescRNAseq/fm_scripts/Gray
mkdir embeddings
if [ ! -f "embeddings/uce.h5ad" ]; then
    rm $FULL_ENV_PATH/UCE/*embeddings*
    python /gpfs/scratch/nk4167/PopScalescRNAseq/fm_scripts/uce_embeddings.py  --input $(pwd)/data.h5ad --output embeddings/uce.h5ad --batch_size 256
    mv $FULL_ENV_PATH/UCE/embeddingsdata_uce_adata.h5ad embeddings/uce.h5ad
fi
cd /gpfs/scratch/nk4167/BreastAtlas/Klughammer_snRNA
mkdir embeddings
if [ ! -f "embeddings/uce.h5ad" ]; then
    rm $FULL_ENV_PATH/UCE/*embeddings*
    python /gpfs/scratch/nk4167/PopScalescRNAseq/fm_scripts/uce_embeddings.py  --input $(pwd)/data.h5ad --output embeddings/uce.h5ad --batch_size 256
    mv $FULL_ENV_PATH/UCE/embeddingsdata_uce_adata.h5ad embeddings/uce.h5ad
fi
cd /gpfs/scratch/nk4167/BreastAtlas/Klughammer_scRNA
mkdir embeddings
if [ ! -f "embeddings/uce.h5ad" ]; then
    rm $FULL_ENV_PATH/UCE/*embeddings*
    python /gpfs/scratch/nk4167/PopScalescRNAseq/fm_scripts/uce_embeddings.py  --input $(pwd)/data.h5ad --output embeddings/uce.h5ad --batch_size 256
    mv $FULL_ENV_PATH/UCE/embeddingsdata_uce_adata.h5ad embeddings/uce.h5ad
fi