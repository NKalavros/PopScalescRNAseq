conda activate /gpfs/scratch/nk4167/miniconda/envs/scgpt_env

if [ ! -f "scrna_embeddings/scgpt.h5ad" ]; then
    python /gpfs/scratch/nk4167/PopScalescRNAseq/fm_scripts/scgpt_embeddings.py  --input lake_scrna_15_genenames.h5ad --output scrna_embeddings/scgpt.h5ad
fi
cd /gpfs/scratch/nk4167/KidneyAtlas/lake

if [ ! -f "snrna_embeddings/scgpt.h5ad" ]; then
    python /gpfs/scratch/nk4167/PopScalescRNAseq/fm_scripts/scgpt_embeddings.py  --input lake_snrna_15_genenames.h5ad --output snrna_embeddings/scgpt.h5ad
fi
cd /gpfs/scratch/nk4167/KidneyAtlas/Abedini

if [ ! -f "embeddings/scgpt.h5ad" ]; then
    python /gpfs/scratch/nk4167/PopScalescRNAseq/fm_scripts/scgpt_embeddings.py  --input GSE211785_Susztak_SC_SN_ATAC_merged_PostSCVI_final_counts.h5ad --output embeddings/scgpt.h5ad
fi
cd /gpfs/scratch/nk4167/KidneyAtlas/SCP1288

if [ ! -f "embeddings/scgpt.h5ad" ]; then
    python /gpfs/scratch/nk4167/PopScalescRNAseq/fm_scripts/scgpt_embeddings.py  --input data.h5ad --output embeddings/scgpt.h5ad
fi
cd /gpfs/scratch/nk4167/KidneyAtlas/Krishna

if [ ! -f "embeddings/scgpt.h5ad" ]; then
    python /gpfs/scratch/nk4167/PopScalescRNAseq/fm_scripts/scgpt_embeddings.py  --input data.h5ad --output embeddings/scgpt.h5ad
fi
cd /gpfs/scratch/nk4167/KidneyAtlas/Braun
rm embeddings/scgpt.h5ad
if [ ! -f "embeddings/scgpt.h5ad" ]; then
    python /gpfs/scratch/nk4167/PopScalescRNAseq/fm_scripts/scgpt_embeddings.py  --input data.h5ad --output embeddings/scgpt.h5ad
fi

# Now let's do some more models - Geneformer
conda activate /gpfs/scratch/nk4167/miniconda/envs/geneformer_env
cd /gpfs/scratch/nk4167/KidneyAtlas/lake
if [ ! -f "scrna_embeddings/geneformer.h5ad" ]; then
    python /gpfs/scratch/nk4167/PopScalescRNAseq/fm_scripts/geneformer_embeddings.py  --input lake_scrna_15_genenames.h5ad --output scrna_embeddings/geneformer.h5ad
fi
cd /gpfs/scratch/nk4167/KidneyAtlas/lake
if [ ! -f "snrna_embeddings/geneformer.h5ad" ]; then
    python /gpfs/scratch/nk4167/PopScalescRNAseq/fm_scripts/geneformer_embeddings.py  --input lake_snrna_16_genenames.h5ad --output snrna_embeddings/geneformer.h5ad
fi
cd /gpfs/scratch/nk4167/KidneyAtlas/Abedini
if [ ! -f "embeddings/geneformer.h5ad" ]; then
    python /gpfs/scratch/nk4167/PopScalescRNAseq/fm_scripts/geneformer_embeddings.py  --input GSE211785_Susztak_SC_SN_ATAC_merged_PostSCVI_final_counts.h5ad --output embeddings/geneformer.h5ad
fi
cd /gpfs/scratch/nk4167/KidneyAtlas/SCP1288
if [ ! -f "embeddings/geneformer.h5ad" ]; then
    python /gpfs/scratch/nk4167/PopScalescRNAseq/fm_scripts/geneformer_embeddings.py  --input data.h5ad --output embeddings/geneformer.h5ad
fi
cd /gpfs/scratch/nk4167/KidneyAtlas/Krishna
if [ ! -f "embeddings/geneformer.h5ad" ]; then
    python /gpfs/scratch/nk4167/PopScalescRNAseq/fm_scripts/geneformer_embeddings.py  --input data.h5ad --output embeddings/geneformer.h5ad
fi
cd /gpfs/scratch/nk4167/KidneyAtlas/Braun
if [ ! -f "embeddings/geneformer.h5ad" ]; then
    python /gpfs/scratch/nk4167/PopScalescRNAseq/fm_scripts/geneformer_embeddings.py  --input data.h5ad --output embeddings/geneformer.h5ad
fi


# Now let's do some more models - scFoundation
conda activate /gpfs/scratch/nk4167/miniconda/envs/scFoundation_env
cd /gpfs/scratch/nk4167/KidneyAtlas/lake
# check if the file output file exists
if [ ! -f "scrna_embeddings/scfoundation.h5ad" ]; then
    python /gpfs/scratch/nk4167/PopScalescRNAseq/fm_scripts/scfoundation_embeddings.py  --input lake_scrna_15_genenames.h5ad --output scrna_embeddings/scfoundation.h5ad
fi
cd /gpfs/scratch/nk4167/KidneyAtlas/lake
if [ ! -f "snrna_embeddings/scfoundation.h5ad" ]; then
    python /gpfs/scratch/nk4167/PopScalescRNAseq/fm_scripts/scfoundation_embeddings.py  --input lake_snrna_16_genenames.h5ad --output snrna_embeddings/scfoundation.h5ad
fi
cd /gpfs/scratch/nk4167/KidneyAtlas/Abedini
if [ ! -f "embeddings/scfoundation.h5ad" ]; then
    python /gpfs/scratch/nk4167/PopScalescRNAseq/fm_scripts/scfoundation_embeddings.py  --input GSE211785_Susztak_SC_SN_ATAC_merged_PostSCVI_final_counts.h5ad --output embeddings/scfoundation.h5ad
fi
cd /gpfs/scratch/nk4167/KidneyAtlas/SCP1288
if [ ! -f "embeddings/scfoundation.h5ad" ]; then
    python /gpfs/scratch/nk4167/PopScalescRNAseq/fm_scripts/scfoundation_embeddings.py  --input data.h5ad --output embeddings/scfoundation.h5ad
fi
cd /gpfs/scratch/nk4167/KidneyAtlas/Krishna
if [ ! -f "embeddings/scfoundation.h5ad" ]; then
    python /gpfs/scratch/nk4167/PopScalescRNAseq/fm_scripts/scfoundation_embeddings.py  --input data.h5ad --output embeddings/scfoundation.h5ad
fi
cd /gpfs/scratch/nk4167/KidneyAtlas/Braun
if [ ! -f "embeddings/scfoundation.h5ad" ]; then
    python /gpfs/scratch/nk4167/PopScalescRNAseq/fm_scripts/scfoundation_embeddings.py  --input data.h5ad --output embeddings/scfoundation.h5ad
fi

# Now let's do some more models - SCIMILARITY
conda activate /gpfs/scratch/nk4167/miniconda/envs/scimilarity_env
cd /gpfs/scratch/nk4167/KidneyAtlas/lake
if [ ! -f "scrna_embeddings/scimilarity.h5ad" ]; then
    python /gpfs/scratch/nk4167/PopScalescRNAseq/fm_scripts/scimilarity_embeddings.py  --input lake_scrna_15_genenames.h5ad --output scrna_embeddings/scimilarity.h5ad
fi
cd /gpfs/scratch/nk4167/KidneyAtlas/lake
if [ ! -f "snrna_embeddings/scimilarity.h5ad" ]; then
    python /gpfs/scratch/nk4167/PopScalescRNAseq/fm_scripts/scimilarity_embeddings.py  --input lake_snrna_16_genenames.h5ad --output snrna_embeddings/scimilarity.h5ad
fi
cd /gpfs/scratch/nk4167/KidneyAtlas/Abedini
if [ ! -f "embeddings/scimilarity.h5ad" ]; then
    python /gpfs/scratch/nk4167/PopScalescRNAseq/fm_scripts/scimilarity_embeddings.py  --input GSE211785_Susztak_SC_SN_ATAC_merged_PostSCVI_final_counts.h5ad --output embeddings/scimilarity.h5ad
fi
cd /gpfs/scratch/nk4167/KidneyAtlas/SCP1288
if [ ! -f "embeddings/scimilarity.h5ad" ]; then
    python /gpfs/scratch/nk4167/PopScalescRNAseq/fm_scripts/scimilarity_embeddings.py  --input data.h5ad --output embeddings/scimilarity.h5ad
fi
cd /gpfs/scratch/nk4167/KidneyAtlas/Krishna
if [ ! -f "embeddings/scimilarity.h5ad" ]; then
    python /gpfs/scratch/nk4167/PopScalescRNAseq/fm_scripts/scimilarity_embeddings.py  --input data.h5ad --output embeddings/scimilarity.h5ad
fi
cd /gpfs/scratch/nk4167/KidneyAtlas/Braun
if [ ! -f "embeddings/scimilarity.h5ad" ]; then
    python /gpfs/scratch/nk4167/PopScalescRNAseq/fm_scripts/scimilarity_embeddings.py  --input data.h5ad --output embeddings/scimilarity.h5ad
fi


# Now let's do some more models - UCE
conda activate /gpfs/scratch/nk4167/miniconda/envs/uce_env
cd /gpfs/scratch/nk4167/KidneyAtlas/lake
if [ ! -f "scrna_embeddings/uce.h5ad" ]; then
    python /gpfs/scratch/nk4167/PopScalescRNAseq/fm_scripts/uce_embeddings.py  --input $(pwd)/lake_scrna_15_genenames.h5ad --output scrna_embeddings/uce.h5ad --batch_size 256
fi
cd /gpfs/scratch/nk4167/KidneyAtlas/lake
if [ ! -f "snrna_embeddings/uce.h5ad" ]; then
    python /gpfs/scratch/nk4167/PopScalescRNAseq/fm_scripts/uce_embeddings.py  --input $(pwd)/lake_snrna_16_genenames.h5ad --output snrna_embeddings/uce.h5ad --batch-size 128
fi
cd /gpfs/scratch/nk4167/KidneyAtlas/Abedini
if [ ! -f "embeddings/uce.h5ad" ]; then
    python /gpfs/scratch/nk4167/PopScalescRNAseq/fm_scripts/uce_embeddings.py  --input $(pwd)/GSE211785_Susztak_SC_SN_ATAC_merged_PostSCVI_final_counts.h5ad --output embeddings/uce.h5ad
fi  
cd /gpfs/scratch/nk4167/KidneyAtlas/SCP1288
if [ ! -f "embeddings/uce.h5ad" ]; then
    python /gpfs/scratch/nk4167/PopScalescRNAseq/fm_scripts/uce_embeddings.py  --input $(pwd)/data.h5ad --output embeddings/uce.h5ad
fi
cd /gpfs/scratch/nk4167/KidneyAtlas/Krishna
if [ ! -f "embeddings/uce.h5ad" ]; then
    python /gpfs/scratch/nk4167/PopScalescRNAseq/fm_scripts/uce_embeddings.py  --input $(pwd)/data.h5ad --output embeddings/uce.h5ad
fi
cd /gpfs/scratch/nk4167/KidneyAtlas/Braun
if [ ! -f "embeddings/uce.h5ad" ]; then
    python /gpfs/scratch/nk4167/PopScalescRNAseq/fm_scripts/uce_embeddings.py  --input $(pwd)/data.h5ad --output embeddings/uce.h5ad
fi  