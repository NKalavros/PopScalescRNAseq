# Download the Breast dataset from Bassez et al.
# Note that you need to login for that.
BASE_DIR=/gpfs/scratch/nk4167/BreastAtlas
OUTPUT_DIR="./Klughammer_scRNA"

echo "Changing to base directory: $BASE_DIR"
cd "$BASE_DIR" || { echo "Error: Could not change to $BASE_DIR"; exit 1; }

echo "Creating output directory: $OUTPUT_DIR"
mkdir -p "$OUTPUT_DIR"
cd "$OUTPUT_DIR" || { echo "Error: Could not change to $OUTPUT_DIR"; exit 1; }


DL_LINK="https://datasets.cellxgene.cziscience.com/487ed451-5945-4cc8-a696-2afeaf1caaf3.h5ad"
curl $DL_LINK --output klughammer_scrna.h5ad

python -c "
import scanpy as sc
import anndata as ad
import os
from pathlib import Path
adata = sc.read_h5ad('klughammer_scrna.h5ad')
adata.var_names = adata.var['feature_name'].tolist()
# Remake the adata object with the raw.X as X
adata = adata.raw.to_adata()
adata.var_names = adata.var['feature_name'].tolist()
adata.write_h5ad('data.h5ad')
"