BASE_DIR=/gpfs/scratch/nk4167/KidneyAtlas
cd "$BASE_DIR"
OUTPUT_DIR="./Stewart"
mkdir -p "$OUTPUT_DIR"
cd "$OUTPUT_DIR"
curl -L -o Stewart.h5ad \
https://cellgeni.cog.sanger.ac.uk/kidneycellatlas/Mature_Full_v3.h5ad