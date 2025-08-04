#https://figshare.com/ndownloader/files/53286068
# Download the Breast dataset from Chen et al.
BASE_DIR=/gpfs/scratch/nk4167/BreastAtlas
OUTPUT_DIR="./Chen"

echo "Changing to base directory: $BASE_DIR"
cd "$BASE_DIR" || { echo "Error: Could not change to $BASE_DIR"; exit 1; }

echo "Creating output directory: $OUTPUT_DIR"
mkdir -p "$OUTPUT_DIR"
cd "$OUTPUT_DIR" || { echo "Error: Could not change to $OUTPUT_DIR"; exit 1; }

echo "Downloading Chen et al. Breast dataset from Figshare..."
curl -L --retry 5 --retry-delay 5 --continue-at - \
  -A "Mozilla/5.0 (Macintosh; Intel Mac OS X 10_15_7) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/115.0.0.0 Safari/537.36" \
  -o breast.zip \
  https://figshare.com/ndownloader/files/53286068
if [ $? -ne 0 ]; then
    echo "Error: Download failed."
    exit 1
fi

echo "Download complete: breast.zip"
unzip breast.zip