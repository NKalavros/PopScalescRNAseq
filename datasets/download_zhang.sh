#RScript with GEOQuery to download GSE159115 supplementary files
BASE_DIR=/gpfs/scratch/nk4167/KidneyAtlas
cd "$BASE_DIR"
GSE_ID="GSE159115"
OUTPUT_DIR="./GSE159115_data"
mkdir -p "$OUTPUT_DIR"
cd "$OUTPUT_DIR"
Rscript -e "
library(GEOquery)
options(timeout = 9999999)
options(download.file.method.GEOquery = 'curl')
getGEOSuppFiles('$GSE_ID', makeDirectory = FALSE)
untar('GSE159115_RAW.tar')"