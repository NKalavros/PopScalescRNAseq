# Downloading from figshare
BASE_DIR=/gpfs/scratch/nk4167/KidneyAtlas
cd "$BASE_DIR"
OUTPUT_DIR="./Quatredeniers"
mkdir -p "$OUTPUT_DIR"
cd "$OUTPUT_DIR"
# Download the dataset from figshare
url=https://figshare.com/ndownloader/files/38263512
curl -L -H "User-Agent: Mozilla/5.0 (compatible; curl)" "$url" -o sc_integrated_Annot.rds
# Convert the RDS file to h5ad format
Rscript -e "
library(Seurat)
library(fastMatMR)
library(Matrix)
# Load the RDS file
seurat_obj <- readRDS('sc_integrated_Annot.rds')
seurat_obj <- UpdateSeuratObject(seurat_obj)
seur_logccounts <- GetAssayData(seurat_obj, slot = 'data')
# Save with fastMatMR
write_fmm(seur_logccounts, 'logcounts.mtx')
write.csv(seurat_obj@meta.data, 'metadata.csv')
"
# Convert to h5ad format
echo "Converting RDS to h5ad format..."
python -c "
import scanpy as sc
import anndata as ad
import pandas as pd
import scipy.io
import fast_matrix_market as fmm
# Load the matrix
X = fmm.mmread('logcounts.mtx')
# Convert to CSR format for h5ad compatibility
X = X.tocsr()
# Load the metadata
metadata = pd.read_csv('metadata.csv', index_col=0)
# Create an AnnData object
adata = ad.AnnData(X=X.T, obs=metadata)
# Save as h5ad
adata.write('data.h5ad')