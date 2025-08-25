#
conda activate /gpfs/scratch/nk4167/miniconda/envs/helical-package
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
    python /gpfs/scratch/nk4167/PopScalescRNAseq/fm_scripts/run_helical_scgpt.py --output_file embeddings/scgpt_helical.h5ad
fi

cd /gpfs/scratch/nk4167/KidneyAtlas/SCP1288
mkdir -p embeddings
if [ ! -f "embeddings/scgpt_helical.h5ad" ]; then
    python /gpfs/scratch/nk4167/PopScalescRNAseq/fm_scripts/run_helical_scgpt.py --output_file embeddings/scgpt_helical.h5ad
fi

cd /gpfs/scratch/nk4167/KidneyAtlas/Krishna
mkdir -p embeddings
if [ ! -f "embeddings/scgpt_helical.h5ad" ]; then
    python /gpfs/scratch/nk4167/PopScalescRNAseq/fm_scripts/run_helical_scgpt.py  --input data.h5ad --output_file embeddings/scgpt_helical.h5ad
fi
cd /gpfs/scratch/nk4167/KidneyAtlas/Braun
mkdir -p embeddings
if [ ! -f "embeddings/scgpt_helical.h5ad" ]; then
    python /gpfs/scratch/nk4167/PopScalescRNAseq/fm_scripts/run_helical_scgpt.py  --output_file embeddings/scgpt_helical.h5ad
fi

# Add Helical UCE
conda activate /gpfs/scratch/nk4167/miniconda/envs/helical-package
cd /gpfs/scratch/nk4167/KidneyAtlas/lake_scrna
mkdir -p embeddings
if [ ! -f "embeddings/helical_uce.h5ad" ]; then
    python /gpfs/scratch/nk4167/PopScalescRNAseq/fm_scripts/run_helical_uce.py --input_file data_with_ensembl.h5ad --output_file embeddings/helical_uce.h5ad --batch_size 256
fi
cd /gpfs/scratch/nk4167/KidneyAtlas/lake_snrna
mkdir -p embeddings
if [ ! -f "embeddings/helical_uce.h5ad" ]; then
    python /gpfs/scratch/nk4167/PopScalescRNAseq/fm_scripts/run_helical_uce.py --input_file data_with_ensembl.h5ad --output_file embeddings/helical_uce.h5ad --batch_size 256
fi
cd /gpfs/scratch/nk4167/KidneyAtlas/Abedini
mkdir -p embeddings
if [ ! -f "embeddings/helical_uce.h5ad" ]; then
    python /gpfs/scratch/nk4167/PopScalescRNAseq/fm_scripts/run_helical_uce.py --input_file GSE211785_Susztak_SC_SN_ATAC_merged_PostSCVI_final_counts.h5ad --output_file embeddings/helical_uce.h5ad  --batch_size 256
fi
cd /gpfs/scratch/nk4167/KidneyAtlas/SCP1288
mkdir -p embeddings
if [ ! -f "embeddings/helical_uce.h5ad" ]; then
    python /gpfs/scratch/nk4167/PopScalescRNAseq/fm_scripts/run_helical_uce.py --input_file data.h5ad --output_file embeddings/helical_uce.h5ad --batch_size 256
fi
cd /gpfs/scratch/nk4167/KidneyAtlas/Krishna
mkdir -p embeddings
if [ ! -f "embeddings/helical_uce.h5ad" ]; then
    python /gpfs/scratch/nk4167/PopScalescRNAseq/fm_scripts/run_helical_uce.py --input_file data.h5ad --output_file embeddings/helical_uce.h5ad --batch_size 256
fi
cd /gpfs/scratch/nk4167/KidneyAtlas/Braun
mkdir -p embeddings
if [ ! -f "embeddings/helical_uce.h5ad" ]; then
    python /gpfs/scratch/nk4167/PopScalescRNAseq/fm_scripts/run_helical_uce.py --input_file data.h5ad --output_file embeddings/helical_uce.h5ad --batch_size 256
fi


# Now let's do some more models - Geneformer (KidneyAtlas paths)
conda activate /gpfs/scratch/nk4167/miniconda/envs/helical-package

cd /gpfs/scratch/nk4167/KidneyAtlas/lake_scrna
mkdir -p embeddings
if [ ! -f "embeddings/geneformer_helical.h5ad" ]; then
    python /gpfs/scratch/nk4167/PopScalescRNAseq/fm_scripts/run_helical_geneformer.py --input data_with_ensembl.h5ad --output embeddings/geneformer_helical.h5ad
fi

cd /gpfs/scratch/nk4167/KidneyAtlas/lake_snrna
mkdir -p embeddings
if [ ! -f "embeddings/geneformer_helical.h5ad" ]; then
    python /gpfs/scratch/nk4167/PopScalescRNAseq/fm_scripts/run_helical_geneformer.py --input data_with_ensembl.h5ad --output embeddings/geneformer_helical.h5ad
fi

cd /gpfs/scratch/nk4167/KidneyAtlas/Abedini
mkdir -p embeddings
if [ ! -f "embeddings/geneformer_helical.h5ad" ]; then
    python /gpfs/scratch/nk4167/PopScalescRNAseq/fm_scripts/run_helical_geneformer.py --input GSE211785_Susztak_SC_SN_ATAC_merged_PostSCVI_final_counts.h5ad --output embeddings/geneformer_helical.h5ad
fi

cd /gpfs/scratch/nk4167/KidneyAtlas/SCP1288
mkdir -p embeddings
if [ ! -f "embeddings/geneformer_helical.h5ad" ]; then
    python /gpfs/scratch/nk4167/PopScalescRNAseq/fm_scripts/run_helical_geneformer.py --input data.h5ad --output embeddings/geneformer_helical.h5ad
fi

cd /gpfs/scratch/nk4167/KidneyAtlas/Krishna
mkdir -p embeddings
if [ ! -f "embeddings/geneformer_helical.h5ad" ]; then
    python /gpfs/scratch/nk4167/PopScalescRNAseq/fm_scripts/run_helical_geneformer.py --input data.h5ad --output embeddings/geneformer_helical.h5ad
fi

cd /gpfs/scratch/nk4167/KidneyAtlas/Braun
mkdir -p embeddings
if [ ! -f "embeddings/geneformer_helical.h5ad" ]; then
    python /gpfs/scratch/nk4167/PopScalescRNAseq/fm_scripts/run_helical_geneformer.py --input data.h5ad --output embeddings/geneformer_helical.h5ad
fi
