BASE_DIR=/gpfs/scratch/nk4167/KidneyAtlas
cd "$BASE_DIR"
OUTPUT_DIR="./lake"
mkdir -p "$OUTPUT_DIR"
cd "$OUTPUT_DIR"
curl -L https://datasets.cellxgene.cziscience.com/1568e555-7c47-4b32-9f29-cda5717e9186.h5ad -o lake_scrna_15.h5ad
curl -L https://datasets.cellxgene.cziscience.com/f5efcb4c-99ca-4afb-9f22-af975525e42f.h5ad -o lake_snrna_16.h5ad
# Load the objects and change the names to genenames
python -c "
import scanpy as sc
import anndata as ad
import os
from pathlib import Path
adata = sc.read_h5ad('lake_scrna_15.h5ad')
adata.var_names = adata.var['feature_name'].tolist()
# Remake the adata object with the raw.X as X
adata = adata.raw.to_adata()
adata.var_names = adata.var['feature_name'].tolist()
adata.write_h5ad('lake_scrna_15_genenames.h5ad')
adata = sc.read_h5ad('lake_snrna_16.h5ad')
adata = adata.raw.to_adata()
adata.var_names = adata.var['feature_name'].tolist()
adata.write_h5ad('lake_snrna_16_genenames.h5ad')"

