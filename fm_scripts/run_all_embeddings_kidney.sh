# NOTE: This code is slightly outdates because I moved the lake files to two actual directories (lake_scrna and lake_snrna)
conda activate /gpfs/scratch/nk4167/miniconda/envs/scgpt_env
cd /gpfs/scratch/nk4167/KidneyAtlas/lake_scrna
mkdir -p embeddings
if [ ! -f "embeddings/scgpt.h5ad" ]; then
    python /gpfs/scratch/nk4167/PopScalescRNAseq/fm_scripts/scgpt_embeddings.py  --input lake_scrna_15_genenames.h5ad --output embeddings/scgpt.h5ad
fi
cd /gpfs/scratch/nk4167/KidneyAtlas/lake_snrna
mkdir -p embeddings
if [ ! -f "embeddings/scgpt.h5ad" ]; then
    python /gpfs/scratch/nk4167/PopScalescRNAseq/fm_scripts/scgpt_embeddings.py  --input lake_snrna_16_genenames.h5ad --output embeddings/scgpt.h5ad
fi
cd /gpfs/scratch/nk4167/KidneyAtlas/Abedini
mkdir -p embeddings
if [ ! -f "embeddings/scgpt.h5ad" ]; then
    python /gpfs/scratch/nk4167/PopScalescRNAseq/fm_scripts/scgpt_embeddings.py  --input GSE211785_Susztak_SC_SN_ATAC_merged_PostSCVI_final_counts.h5ad --output embeddings/scgpt.h5ad
fi
cd /gpfs/scratch/nk4167/KidneyAtlas/SCP1288
mkdir -p embeddings
if [ ! -f "embeddings/scgpt.h5ad" ]; then
    python /gpfs/scratch/nk4167/PopScalescRNAseq/fm_scripts/scgpt_embeddings.py  --input data.h5ad --output embeddings/scgpt.h5ad
fi
cd /gpfs/scratch/nk4167/KidneyAtlas/Krishna
mkdir -p embeddings
if [ ! -f "embeddings/scgpt.h5ad" ]; then
    python /gpfs/scratch/nk4167/PopScalescRNAseq/fm_scripts/scgpt_embeddings.py  --input data.h5ad --output embeddings/scgpt.h5ad
fi
cd /gpfs/scratch/nk4167/KidneyAtlas/Braun
mkdir -p embeddings
if [ ! -f "embeddings/scgpt.h5ad" ]; then
    python /gpfs/scratch/nk4167/PopScalescRNAseq/fm_scripts/scgpt_embeddings.py  --input data.h5ad --output embeddings/scgpt.h5ad
fi

# Now let's do some more models - Geneformer
conda activate /gpfs/scratch/nk4167/miniconda/envs/geneformer_env
cd /gpfs/scratch/nk4167/KidneyAtlas/lake_scrna
mkdir -p embeddings
if [ ! -f "embeddings/geneformer.h5ad" ]; then
    python /gpfs/scratch/nk4167/PopScalescRNAseq/fm_scripts/geneformer_embeddings.py  --input lake_scrna_15_genenames.h5ad --output embeddings/geneformer.h5ad
fi
cd /gpfs/scratch/nk4167/KidneyAtlas/lake_snrna
mkdir -p embeddings
if [ ! -f "embeddings/geneformer.h5ad" ]; then
    python /gpfs/scratch/nk4167/PopScalescRNAseq/fm_scripts/geneformer_embeddings.py  --input lake_snrna_16_genenames.h5ad --output embeddings/geneformer.h5ad
fi
cd /gpfs/scratch/nk4167/KidneyAtlas/Abedini
mkdir -p embeddings
if [ ! -f "embeddings/geneformer.h5ad" ]; then
    python /gpfs/scratch/nk4167/PopScalescRNAseq/fm_scripts/geneformer_embeddings.py  --input GSE211785_Susztak_SC_SN_ATAC_merged_PostSCVI_final_counts.h5ad --output embeddings/geneformer.h5ad
fi
cd /gpfs/scratch/nk4167/KidneyAtlas/SCP1288
mkdir -p embeddings
if [ ! -f "embeddings/geneformer.h5ad" ]; then
    python /gpfs/scratch/nk4167/PopScalescRNAseq/fm_scripts/geneformer_embeddings.py  --input data.h5ad --output embeddings/geneformer.h5ad
fi
cd /gpfs/scratch/nk4167/KidneyAtlas/Krishna
mkdir -p embeddings
if [ ! -f "embeddings/geneformer.h5ad" ]; then
    python /gpfs/scratch/nk4167/PopScalescRNAseq/fm_scripts/geneformer_embeddings.py  --input data.h5ad --output embeddings/geneformer.h5ad
fi
cd /gpfs/scratch/nk4167/KidneyAtlas/Braun
mkdir -p embeddings
if [ ! -f "embeddings/geneformer.h5ad" ]; then
    python /gpfs/scratch/nk4167/PopScalescRNAseq/fm_scripts/geneformer_embeddings.py  --input data.h5ad --output embeddings/geneformer.h5ad
fi


# Now let's do some more models - scFoundation
conda activate /gpfs/scratch/nk4167/miniconda/envs/scFoundation_env
cd /gpfs/scratch/nk4167/KidneyAtlas/lake_scrna
mkdir -p embeddings
# check if the file output file exists
if [ ! -f "embeddings/scfoundation.h5ad" ]; then
    python /gpfs/scratch/nk4167/PopScalescRNAseq/fm_scripts/scfoundation_embeddings.py  --input lake_scrna_15_genenames.h5ad --output embeddings/scfoundation.h5ad
fi
cd /gpfs/scratch/nk4167/KidneyAtlas/lake_snrna
mkdir -p embeddings
if [ ! -f "embeddings/scfoundation.h5ad" ]; then
    python /gpfs/scratch/nk4167/PopScalescRNAseq/fm_scripts/scfoundation_embeddings.py  --input lake_snrna_16_genenames.h5ad --output embeddings/scfoundation.h5ad
fi
cd /gpfs/scratch/nk4167/KidneyAtlas/Abedini
mkdir -p embeddings
if [ ! -f "embeddings/scfoundation.h5ad" ]; then
    python /gpfs/scratch/nk4167/PopScalescRNAseq/fm_scripts/scfoundation_embeddings.py  --input GSE211785_Susztak_SC_SN_ATAC_merged_PostSCVI_final_counts.h5ad --output embeddings/scfoundation.h5ad
fi
cd /gpfs/scratch/nk4167/KidneyAtlas/SCP1288
mkdir -p embeddings
if [ ! -f "embeddings/scfoundation.h5ad" ]; then
    python /gpfs/scratch/nk4167/PopScalescRNAseq/fm_scripts/scfoundation_embeddings.py  --input data.h5ad --output embeddings/scfoundation.h5ad
fi
cd /gpfs/scratch/nk4167/KidneyAtlas/Krishna
mkdir -p embeddings
if [ ! -f "embeddings/scfoundation.h5ad" ]; then
    python /gpfs/scratch/nk4167/PopScalescRNAseq/fm_scripts/scfoundation_embeddings.py  --input data.h5ad --output embeddings/scfoundation.h5ad
fi
cd /gpfs/scratch/nk4167/KidneyAtlas/Braun
mkdir -p embeddings
if [ ! -f "embeddings/scfoundation.h5ad" ]; then
    python /gpfs/scratch/nk4167/PopScalescRNAseq/fm_scripts/scfoundation_embeddings.py  --input data.h5ad --output embeddings/scfoundation.h5ad
fi

# Now let's do some more models - SCIMILARITY
conda activate /gpfs/scratch/nk4167/miniconda/envs/scimilarity_env
cd /gpfs/scratch/nk4167/KidneyAtlas/lake_scrna
mkdir -p embeddings
if [ ! -f "embeddings/scimilarity.h5ad" ]; then
    python /gpfs/scratch/nk4167/PopScalescRNAseq/fm_scripts/scimilarity_embeddings.py  --input lake_scrna_15_genenames.h5ad --output embeddings/scimilarity.h5ad
fi
cd /gpfs/scratch/nk4167/KidneyAtlas/lake_snrna
mkdir -p embeddings
if [ ! -f "embeddings/scimilarity.h5ad" ]; then
    python /gpfs/scratch/nk4167/PopScalescRNAseq/fm_scripts/scimilarity_embeddings.py  --input lake_snrna_16_genenames.h5ad --output embeddings/scimilarity.h5ad
fi
cd /gpfs/scratch/nk4167/KidneyAtlas/Abedini
mkdir -p embeddings
if [ ! -f "embeddings/scimilarity.h5ad" ]; then
    python /gpfs/scratch/nk4167/PopScalescRNAseq/fm_scripts/scimilarity_embeddings.py  --input GSE211785_Susztak_SC_SN_ATAC_merged_PostSCVI_final_counts.h5ad --output embeddings/scimilarity.h5ad
fi
cd /gpfs/scratch/nk4167/KidneyAtlas/SCP1288
mkdir -p embeddings
if [ ! -f "embeddings/scimilarity.h5ad" ]; then
    python /gpfs/scratch/nk4167/PopScalescRNAseq/fm_scripts/scimilarity_embeddings.py  --input data.h5ad --output embeddings/scimilarity.h5ad
fi
cd /gpfs/scratch/nk4167/KidneyAtlas/Krishna
mkdir -p embeddings
if [ ! -f "embeddings/scimilarity.h5ad" ]; then
    python /gpfs/scratch/nk4167/PopScalescRNAseq/fm_scripts/scimilarity_embeddings.py  --input data.h5ad --output embeddings/scimilarity.h5ad
fi
cd /gpfs/scratch/nk4167/KidneyAtlas/Braun
mkdir -p embeddings
if [ ! -f "embeddings/scimilarity.h5ad" ]; then
    python /gpfs/scratch/nk4167/PopScalescRNAseq/fm_scripts/scimilarity_embeddings.py  --input data.h5ad --output embeddings/scimilarity.h5ad
fi


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
    bash /gpfs/scratch/nk4167/PopScalescRNAseq/fm_scripts/run_transcriptformer.sh data.h5ad
fi

cd /gpfs/scratch/nk4167/KidneyAtlas/lake_snrna
mkdir -p embeddings
if [ ! -f "embeddings/transcriptformer.h5ad" ]; then
    bash /gpfs/scratch/nk4167/PopScalescRNAseq/fm_scripts/run_transcriptformer.sh data.h5ad
fi

cd /gpfs/scratch/nk4167/KidneyAtlas/Abedini
mkdir -p embeddings
if [ ! -f "embeddings/transcriptformer.h5ad" ]; then
    bash /gpfs/scratch/nk4167/PopScalescRNAseq/fm_scripts/run_transcriptformer.sh data.h5ad
fi

cd /gpfs/scratch/nk4167/KidneyAtlas/SCP1288
mkdir -p embeddings
if [ ! -f "embeddings/transcriptformer.h5ad" ]; then
    bash /gpfs/scratch/nk4167/PopScalescRNAseq/fm_scripts/run_transcriptformer.sh data.h5ad
fi

cd /gpfs/scratch/nk4167/KidneyAtlas/Krishna
mkdir -p embeddings
if [ ! -f "embeddings/transcriptformer.h5ad" ]; then
    bash /gpfs/scratch/nk4167/PopScalescRNAseq/fm_scripts/run_transcriptformer.sh data.h5ad
fi

cd /gpfs/scratch/nk4167/KidneyAtlas/Braun
mkdir -p embeddings
if [ ! -f "embeddings/transcriptformer.h5ad" ]; then
    bash /gpfs/scratch/nk4167/PopScalescRNAseq/fm_scripts/run_transcriptformer.sh data.h5ad
fi
