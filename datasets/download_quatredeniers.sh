# Downloading from figshare
BASE_DIR=/gpfs/scratch/nk4167/KidneyAtlas
cd "$BASE_DIR"
OUTPUT_DIR="./Quatredeniers"
mkdir -p "$OUTPUT_DIR"
cd "$OUTPUT_DIR"
# Download the dataset from figshare
url=https://figshare.com/ndownloader/files/38263512
curl -L -H "User-Agent: Mozilla/5.0 (compatible; curl)" "$url" -o sc_integrated_Annot.rds