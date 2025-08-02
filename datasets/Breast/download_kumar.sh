cd /gpfs/scratch/nk4167/BreastAtlas
OUTPUT_DIR="./Kumar"
mkdir -p "$OUTPUT_DIR"
cd "$OUTPUT_DIR"
if [ ! -f "kumar.h5ad" ]; then
    curl -L https://datasets.cellxgene.cziscience.com/07d2b187-26d9-4265-bb16-9f04e8aca56a.h5ad -o "kumar.h5ad"
fi

python -c "
import scanpy as sc
import anndata as ad
import os
from pathlib import Path
adata = sc.read_h5ad('kumar.h5ad')
adata.var_names = adata.var['feature_name'].tolist()
# Remake the adata object with the raw.X as X
adata = adata.raw.to_adata()
print(adata.X.max(), adata.X.min())
adata.var_names = adata.var['feature_name'].tolist()
# Print the first genenames
print(adata.var_names[:10])

adata.write_h5ad('kumar_geneneames.h5ad')"