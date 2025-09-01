# NOTE: This code is slightly outdates because I moved the lake files to two actual directories (lake_scrna and lake_snrna)
conda activate /gpfs/scratch/nk4167/miniconda/envs/scgpt_env
cd /gpfs/scratch/nk4167/PancreasAtlas/
# Loop over each directory in PancreasAtlas
for dir in */ ; do
    echo "Processing directory: $dir"
    # Enter the directory
    cd "$dir"
    # Create embeddings directory if it doesn't exist
    mkdir -p embeddings
    # Check if the embeddings file already exists
    if [ ! -f "embeddings/scgpt.h5ad" ]; then
        mkdir embeddings
        python /gpfs/scratch/nk4167/PopScalescRNAseq/fm_scripts/scgpt_embeddings.py  --input data.h5ad --output embeddings/scgpt.h5ad
    fi
    cd ..
done
# Now let's do some more models - scFoundation
conda activate /gpfs/scratch/nk4167/miniconda/envs/scFoundation_env
cd /gpfs/scratch/nk4167/PancreasAtlas/
for dir in */ ; do
    echo "Processing directory: $dir"
    # Enter the directory
    cd "$dir"
    # Create embeddings directory if it doesn't exist
    mkdir -p embeddings
    # Check if the embeddings file already exists
    if [ ! -f "embeddings/scfoundation.h5ad" ]; then
        mkdir embeddings
        python /gpfs/scratch/nk4167/PopScalescRNAseq/fm_scripts/scfoundation_embeddings.py  --input data.h5ad --output embeddings/scfoundation.h5ad
    fi
    cd ..
done


# Now let's do some more models - SCIMILARITY (BreastAtlas)
conda activate /gpfs/scratch/nk4167/miniconda/envs/scimilarity_env
cd /gpfs/scratch/nk4167/PancreasAtlas/
for dir in */ ; do
    echo "Processing directory: $dir"
    # Enter the directory
    cd "$dir"
    # Create embeddings directory if it doesn't exist
    mkdir -p embeddings
    # Check if the embeddings file already exists
    if [ ! -f "embeddings/scimilarity.h5ad" ]; then
        mkdir embeddings
        python /gpfs/scratch/nk4167/PopScalescRNAseq/fm_scripts/scimilarity_embeddings.py  --input data.h5ad --output embeddings/scimilarity.h5ad
    fi
    cd ..
done

# Let's move on to transcriptformer
source ~/.bashrc #Base environment
cd /gpfs/scratch/nk4167/PancreasAtlas/
for dir in */ ; do
    echo "Processing directory: $dir"
    # Enter the directory
    cd "$dir"
    # Create embeddings directory if it doesn't exist
    mkdir -p embeddings
    # Check if the embeddings file already exists
    if [ ! -f "embeddings/transcriptformer.h5ad" ]; then
        bash /gpfs/scratch/nk4167/PopScalescRNAseq/fm_scripts/run_transcriptformer.sh
    fi
    cd ..
done
# Now let's do some more models - UCE (BreastAtlas)
conda activate /gpfs/scratch/nk4167/miniconda/envs/uce_env
FULL_ENV_PATH='/gpfs/scratch/nk4167/miniconda/envs/uce_env'
cd /gpfs/scratch/nk4167/PancreasAtlas/
for dir in */ ; do
    echo "Processing directory: $dir"
    # Enter the directory
    cd "$dir"
    # Create embeddings directory if it doesn't exist
    mkdir -p embeddings
    # Check if the embeddings file already exists
    if [ ! -f "embeddings/uce.h5ad" ]; then
        rm $FULL_ENV_PATH/UCE/*embeddings*
        python /gpfs/scratch/nk4167/PopScalescRNAseq/fm_scripts/uce_embeddings.py  --input $(pwd)/data.h5ad --output embeddings/uce.h5ad --batch_size 256
        mv $(pwd)/embeddingsdata_uce_adata.h5ad embeddings/uce.h5ad
    fi
    cd ..
done

# scGPT through helical
conda activate /gpfs/scratch/nk4167/miniconda/envs/helical-package
cd /gpfs/scratch/nk4167/PancreasAtlas/
for dir in */ ; do
    echo "Processing directory: $dir"
    # Enter the directory
    cd "$dir"
    # Create embeddings directory if it doesn't exist
    mkdir -p embeddings
    # Check if the embeddings file already exists
    if [ ! -f "embeddings/scgpt_helical.h5ad" ]; then
        python /gpfs/scratch/nk4167/PopScalescRNAseq/fm_scripts/run_helical_scgpt.py --output_file embeddings/scgpt_helical.h5ad
    fi
    cd ..
done

# Geneformer through helical (only one that works)
conda activate /gpfs/scratch/nk4167/miniconda/envs/helical-package
cd /gpfs/scratch/nk4167/PancreasAtlas/
for dir in */ ; do
    echo "Processing directory: $dir"
    # Enter the directory
    cd "$dir"
    # Create embeddings directory if it doesn't exist
    mkdir -p embeddings
    # Check if the embeddings file already exists
    if [ ! -f "embeddings/geneformer_helical.h5ad" ]; then
        python /gpfs/scratch/nk4167/PopScalescRNAseq/fm_scripts/run_helical_geneformer.py --output_file embeddings/geneformer_helical.h5ad
    fi
    cd ..
done


# Run them in reverse too

# Geneformer through helical (only one that works)
conda activate /gpfs/scratch/nk4167/miniconda/envs/helical-package
cd /gpfs/scratch/nk4167/PancreasAtlas/
for dir in $(ls -d */ | sort -r) ; do
    echo "Processing directory: $dir"
    # Enter the directory
    cd "$dir"
    # Create embeddings directory if it doesn't exist
    mkdir -p embeddings
    # Check if the embeddings file already exists
    if [ ! -f "embeddings/geneformer_helical.h5ad" ]; then
        python /gpfs/scratch/nk4167/PopScalescRNAseq/fm_scripts/run_helical_geneformer.py --output_file embeddings/geneformer_helical.h5ad
    fi
    cd ..
done
