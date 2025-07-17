BASE_DIR=/gpfs/scratch/nk4167/KidneyAtlas
cd "$BASE_DIR"
OUTPUT_DIR="./Krishna"
mkdir -p "$OUTPUT_DIR"
cd "$OUTPUT_DIR"
link="https://sra-download.be-md.ncbi.nlm.nih.gov/vast/sra01/SRZ/000190/SRZ190804/ccRCC_6pat_Seurat"
curl -L "$link" -o ccRCC_6pat_Seurat.rds