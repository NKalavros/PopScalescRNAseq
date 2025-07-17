BASE_DIR=/gpfs/scratch/nk4167/KidneyAtlas
cd "$BASE_DIR"
PROJECT_ID="SCP1288"
mkdir -p "$PROJECT_ID"
cd $PROJECT_ID
curl -k "https://singlecell.broadinstitute.org/single_cell/api/v1/bulk_download/generate_curl_config?accessions=SCP1288&auth_code=XLUWmoLs&directory=all&context=study"  -o cfg.txt; curl -K cfg.txt && rm cfg.txt

# Let's make an h5ad file from these
# ├── cluster
# │   └── Final_SCP_ClusterFile.txt
# ├── expression
# │   ├── 60c76a18771a5b0ba10ea91b
# │   │   ├── barcodes.tsv
# │   │   ├── genes.tsv
# │   │   └── matrix.mtx
# │   └── ccRCC_scRNASeq_NormalizedCounts.txt.gz
# └── metadata
#     └── Final_SCP_Metadata.txt

python -c "
import scanpy as sc
import pandas as pd
import os

# Read the matrix market format data
adata = sc.read_10x_mtx(
    'expression/60c76a18771a5b0ba10ea91b',
    var_names='gene_symbols',
    cache=True
)

# Make variable names unique
adata.var_names_make_unique()

# Read metadata
metadata = pd.read_csv('metadata/Final_SCP_Metadata.txt', sep='\t', index_col=0)
adata.obs = adata.obs.join(metadata, how='left')


# Convert all non-string columns to strings to avoid h5py errors
for col in adata.obs.columns:
    if adata.obs[col].dtype == 'object':
        adata.obs[col] = adata.obs[col].astype(str)
    elif adata.obs[col].dtype in ['int64', 'float64']:
        adata.obs[col] = adata.obs[col].astype(str)
sc.pp.filter_cells(adata, min_genes=200)  # Filter cells with at least 200 genes
sc.pp.filter_cells(adata, min_counts=400)  # Filter cells with at least 200 genes
sc.pp.filter_genes(adata, min_cells=3)  # Filter genes expressed in
sc.pp.normalize_total(adata, target_sum=1e4)  # Normalize counts
sc.pp.log1p(adata)  # Log-transform the data
# Save as h5ad
adata.write('data.h5ad')
print(f'Saved h5ad file: data.h5ad')
print(f'Shape: {adata.shape}')
"