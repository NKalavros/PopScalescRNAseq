#NOTE THIS ONE I ACTUALLY DID MYSELF
BASE_DIR=/gpfs/scratch/nk4167/KidneyAtlas
cd "$BASE_DIR"
OUTPUT_DIR="./Braun"
mkdir -p "$OUTPUT_DIR"
cd "$OUTPUT_DIR"
curl -O https://pmc.ncbi.nlm.nih.gov/articles/instance/8138872/bin/NIHMS1692222-supplement-supplementary_Data_S1.zip
curl -O https://pmc.ncbi.nlm.nih.gov/articles/instance/8138872/bin/NIHMS1692222-supplement-supplementary_Data_S2.csv