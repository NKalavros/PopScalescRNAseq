BASE_DIR=/gpfs/scratch/nk4167/KidneyAtlas
cd "$BASE_DIR"
OUTPUT_DIR="./lake"
mkdir -p "$OUTPUT_DIR"
cd "$OUTPUT_DIR"
curl -L https://datasets.cellxgene.cziscience.com/1568e555-7c47-4b32-9f29-cda5717e9186.h5ad -o lake_scrna_15.h5ad
curl -L https://datasets.cellxgene.cziscience.com/f5efcb4c-99ca-4afb-9f22-af975525e42f.h5ad -o lake_snrna_16.h5ad

conda activate /gpfs/scratch/nk4167/miniconda/envs/scGPT