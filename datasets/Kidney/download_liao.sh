BASE_DIR=/gpfs/scratch/nk4167/KidneyAtlas
cd "$BASE_DIR"
OUTPUT_DIR="./Liao"
mkdir -p "$OUTPUT_DIR"
cd "$OUTPUT_DIR"
GEO_ACCESSION=GSE131685
Rscript -e "
library(GEOquery)
options(timeout = 9999999)
options(download.file.method.GEOquery = 'curl')
getGEOSuppFiles('$GEO_ACCESSION', makeDirectory = FALSE)
untar('GSE131685_RAW.tar')"