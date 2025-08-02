BASE_DIR=/gpfs/scratch/nk4167/KidneyAtlas
cd "$BASE_DIR"
OUTPUT_DIR="./McEvoy"
mkdir -p "$OUTPUT_DIR"
cd "$OUTPUT_DIR"
expr_matrix_link=https://cells.ucsc.edu/living-donor-kidney/exprMatrix.tsv.gz
curl -L "$expr_matrix_link" -o exprMatrix.tsv.gz
metadata_link=https://cells.ucsc.edu/living-donor-kidney/meta.tsv
curl -L "$metadata_link" -o meta.tsv
umap_coords_link=https://cells.ucsc.edu/living-donor-kidney/UMAP.coords.tsv.gz
curl -L "$umap_coords_link" -o UMAP.coords.tsv.gz
# make it into anndata
python -c "
import pandas as pd
import anndata as ad
import scipy.io
import datatable as dt
expr_matrix = dt.fread('exprMatrix.tsv.gz', fill=True).to_pandas()
# First column is gene names, set it as index and drop it
expr_matrix.set_index(expr_matrix.columns[0], inplace=True)
# Load the metadata
metadata = pd.read_csv('meta.tsv', index_col=0, sep='\\t')
# Create an AnnData object
adata = ad.AnnData(X=expr_matrix.T, obs=metadata) #this is already log counts
# Save as h5ad
adata.write('data.h5ad')