# Now let's do some more models - Geneformer (BreastAtlas)
conda activate /gpfs/scratch/nk4167/miniconda/envs/helical-package
cd /gpfs/scratch/nk4167/BreastAtlas/Bassez
mkdir embeddings
if [ ! -f "embeddings/geneformer_helical.h5ad" ]; then
    python /gpfs/scratch/nk4167/PopScalescRNAseq/fm_scripts/run_helical_geneformer.py  --input data.h5ad --output embeddings/geneformer_helical.h5ad
fi
cd /gpfs/scratch/nk4167/BreastAtlas/Wu
mkdir embeddings
if [ ! -f "embeddings/geneformer_helical.h5ad" ]; then
    python /gpfs/scratch/nk4167/PopScalescRNAseq/fm_scripts/run_helical_geneformer.py  --input data.h5ad --output embeddings/geneformer_helical.h5ad
fi
cd /gpfs/scratch/nk4167/BreastAtlas/Pal
mkdir embeddings
if [ ! -f "embeddings/geneformer_helical.h5ad" ]; then
    python /gpfs/scratch/nk4167/PopScalescRNAseq/fm_scripts/run_helical_geneformer.py  --input data.h5ad --output embeddings/geneformer_helical.h5ad
fi
cd /gpfs/scratch/nk4167/BreastAtlas/Gray
mkdir embeddings
if [ ! -f "embeddings/geneformer_helical.h5ad" ]; then
    python /gpfs/scratch/nk4167/PopScalescRNAseq/fm_scripts/run_helical_geneformer.py  --input data.h5ad --output embeddings/geneformer_helical.h5ad
fi
cd /gpfs/scratch/nk4167/BreastAtlas/Klughammer_snRNA
mkdir embeddings
if [ ! -f "embeddings/geneformer_helical.h5ad" ]; then
    python /gpfs/scratch/nk4167/PopScalescRNAseq/fm_scripts/run_helical_geneformer.py  --input data.h5ad --output embeddings/geneformer_helical.h5ad
fi
cd /gpfs/scratch/nk4167/BreastAtlas/Klughammer_scRNA
mkdir embeddings
if [ ! -f "embeddings/geneformer_helical.h5ad" ]; then
    python /gpfs/scratch/nk4167/PopScalescRNAseq/fm_scripts/run_helical_geneformer.py  --input data.h5ad --output embeddings/geneformer_helical.h5ad
fi
