# Now let's do some more models - UCE
conda activate /gpfs/scratch/nk4167/miniconda/envs/uce_env
FULL_ENV_PATH='/gpfs/scratch/nk4167/miniconda/envs/uce_env'
cd /gpfs/scratch/nk4167/KidneyAtlas/lake_scrna
mkdir -p embeddings
mkdir -p scrna_embeddings
if [ ! -f "embeddings/uce.h5ad" ]; then
    python /gpfs/scratch/nk4167/PopScalescRNAseq/fm_scripts/uce_embeddings.py  --input $(pwd)/lake_scrna_15_genenames.h5ad --output embeddings/uce.h5ad --batch_size 256
fi
mv $FULL_ENV_PATH/UCE/scrna_embeddingslake_scrna_15_genenames_uce_adata.h5ad scrna_embeddings/uce.h5ad
cd /gpfs/scratch/nk4167/KidneyAtlas/lake_snrna
mkdir -p embeddings
mkdir -p snrna_embeddings
if [ ! -f "embeddings/uce.h5ad" ]; then
    python /gpfs/scratch/nk4167/PopScalescRNAseq/fm_scripts/uce_embeddings.py  --input $(pwd)/lake_snrna_16_genenames.h5ad --output embeddings/uce.h5ad --batch_size 128
    mv $FULL_ENV_PATH/UCE/snrna_embeddingslake_snrna_16_genenames_uce_adata.h5ad snrna_embeddings/uce.h5ad
fi
cd /gpfs/scratch/nk4167/KidneyAtlas/Abedini
mkdir -p embeddings
if [ ! -f "embeddings/uce.h5ad" ]; then
    python /gpfs/scratch/nk4167/PopScalescRNAseq/fm_scripts/uce_embeddings.py  --input $(pwd)/GSE211785_Susztak_SC_SN_ATAC_merged_PostSCVI_final_counts.h5ad --output embeddings/uce.h5ad
    mv $FULL_ENV_PATH/UCE/embeddingsGSE211785_Susztak_SC_SN_ATAC_merged_PostSCVI_final_counts_uce_adata.h5ad embeddings/uce.h5ad
fi
cd /gpfs/scratch/nk4167/KidneyAtlas/SCP1288
mkdir -p embeddings
if [ ! -f "embeddings/uce.h5ad" ]; then
    python /gpfs/scratch/nk4167/PopScalescRNAseq/fm_scripts/uce_embeddings.py  --input $(pwd)/data.h5ad --output embeddings/uce.h5ad
    mv $FULL_ENV_PATH/UCE/embeddingsdata_uce_adata.h5ad embeddings/uce.h5ad
fi
cd /gpfs/scratch/nk4167/KidneyAtlas/Krishna
mkdir -p embeddings
if [ ! -f "embeddings/uce.h5ad" ]; then
    python /gpfs/scratch/nk4167/PopScalescRNAseq/fm_scripts/uce_embeddings.py  --input $(pwd)/data.h5ad --output embeddings/uce.h5ad
    mv $FULL_ENV_PATH/UCE/embeddingsdata_uce_adata.h5ad embeddings/uce.h5ad
fi
cd /gpfs/scratch/nk4167/KidneyAtlas/Braun
mkdir -p embeddings
if [ ! -f "embeddings/uce.h5ad" ]; then
    python /gpfs/scratch/nk4167/PopScalescRNAseq/fm_scripts/uce_embeddings.py  --input $(pwd)/data.h5ad --output embeddings/uce.h5ad
    mv $FULL_ENV_PATH/UCE/embeddingsdata_uce_adata.h5ad embeddings/uce.h5ad
fi


# Now let's do some more models - TranscriptFormer
cd /gpfs/scratch/nk4167/KidneyAtlas/lake_scrna
mkdir -p embeddings
if [ ! -f "embeddings/transcriptformer.h5ad" ]; then
    bash /gpfs/scratch/nk4167/PopScalescRNAseq/fm_scripts/run_transcriptformer.sh
fi

cd /gpfs/scratch/nk4167/KidneyAtlas/lake_snrna
mkdir -p embeddings
if [ ! -f "embeddings/transcriptformer.h5ad" ]; then
    bash /gpfs/scratch/nk4167/PopScalescRNAseq/fm_scripts/run_transcriptformer.sh
fi

cd /gpfs/scratch/nk4167/KidneyAtlas/Abedini
mkdir -p embeddings
if [ ! -f "embeddings/transcriptformer.h5ad" ]; then
    bash /gpfs/scratch/nk4167/PopScalescRNAseq/fm_scripts/run_transcriptformer.sh
fi

cd /gpfs/scratch/nk4167/KidneyAtlas/SCP1288
mkdir -p embeddings
if [ ! -f "embeddings/transcriptformer.h5ad" ]; then
    bash /gpfs/scratch/nk4167/PopScalescRNAseq/fm_scripts/run_transcriptformer.sh
fi

cd /gpfs/scratch/nk4167/KidneyAtlas/Krishna
mkdir -p embeddings
if [ ! -f "embeddings/transcriptformer.h5ad" ]; then
    bash /gpfs/scratch/nk4167/PopScalescRNAseq/fm_scripts/run_transcriptformer.sh
fi
cd /gpfs/scratch/nk4167/KidneyAtlas/Braun
mkdir -p embeddings
if [ ! -f "embeddings/transcriptformer.h5ad" ]; then
    bash /gpfs/scratch/nk4167/PopScalescRNAseq/fm_scripts/run_transcriptformer.sh
fi


# Now let's do some more models - scGPT
cd /gpfs/scratch/nk4167/KidneyAtlas/lake_scrna
mkdir -p embeddings
if [ ! -f "embeddings/scgpt_helical.h5ad" ]; then
    python /gpfs/scratch/nk4167/PopScalescRNAseq/fm_scripts/run_helical_scgpt.py --input_file data_with_ensembl.h5ad --output_file embeddings/scgpt_helical.h5ad
fi

cd /gpfs/scratch/nk4167/KidneyAtlas/lake_snrna
mkdir -p embeddings
if [ ! -f "embeddings/scgpt_helical.h5ad" ]; then
    python /gpfs/scratch/nk4167/PopScalescRNAseq/fm_scripts/run_helical_scgpt.py --input_file data_with_ensembl.h5ad --output_file embeddings/scgpt_helical.h5ad
fi

cd /gpfs/scratch/nk4167/KidneyAtlas/Abedini
mkdir -p embeddings
if [ ! -f "embeddings/scgpt_helical.h5ad" ]; then
    python /gpfs/scratch/nk4167/PopScalescRNAseq/fm_scripts/run_helical_scgpt.py
fi

cd /gpfs/scratch/nk4167/KidneyAtlas/SCP1288
mkdir -p embeddings
if [ ! -f "embeddings/scgpt_helical.h5ad" ]; then
    python /gpfs/scratch/nk4167/PopScalescRNAseq/fm_scripts/run_helical_scgpt.py
fi

cd /gpfs/scratch/nk4167/KidneyAtlas/Krishna
mkdir -p embeddings
if [ ! -f "embeddings/scgpt_helical.h5ad" ]; then
    python /gpfs/scratch/nk4167/PopScalescRNAseq/fm_scripts/run_helical_scgpt.py
fi
cd /gpfs/scratch/nk4167/KidneyAtlas/Braun
mkdir -p embeddings
if [ ! -f "embeddings/scgpt_helical.h5ad" ]; then
    python /gpfs/scratch/nk4167/PopScalescRNAseq/fm_scripts/run_helical_scgpt.py
fi

