BASE_DIR=/gpfs/scratch/nk4167/KidneyAtlas
cd "$BASE_DIR"
OUTPUT_DIR="./Krishna"
mkdir -p "$OUTPUT_DIR"
cd "$OUTPUT_DIR"
link="https://sra-download.be-md.ncbi.nlm.nih.gov/vast/sra01/SRZ/000190/SRZ190804/ccRCC_6pat_Seurat"
curl -L "$link" -o ccRCC_6pat_Seurat.rds
Rscript -e "
library(Seurat)
library(fastMatMR)
library(Matrix)
# Load the RDS file
seurat_obj <- readRDS('ccRCC_6pat_Seurat.rds')
seurat_obj <- UpdateSeuratObject(seurat_obj)
seur_counts <- GetAssayData(seurat_obj, slot = 'counts')
# Save with fastMatMR
write_fmm(seur_counts, 'counts.mtx')
# Write rownames and colnames
write.table(rownames(seur_counts), 'genes.txt', row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(colnames(seur_counts), 'barcodes.txt', row.names = FALSE, col.names = FALSE, quote = FALSE)
# Save metadata
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
X = fmm.mmread('counts.mtx')
# Convert to CSR format for h5ad compatibility
X = X.tocsr()
# Load the metadata
metadata = pd.read_csv('metadata.csv', index_col=0)
# Create an AnnData object
adata = ad.AnnData(X=X.T, obs=metadata)
# Set var_names from genes.txt
with open('genes.txt') as f:
    adata.var_names = [line.strip() for line in f]
# Set obs_names from barcodes.txt
with open('barcodes.txt') as f:
    adata.obs_names = [line.strip() for line in f]

sc.pp.filter_cells(adata, min_genes=200)  # Filter cells with at least 200 genes
sc.pp.filter_cells(adata, min_counts=400)  # Filter cells with at least 400 counts
sc.pp.filter_genes(adata, min_cells=3)  # Filter genes expressed in at least 3 cells
# Fix the names with mygene from ensembl IDs to symbols
adata.var_names_make_unique()
import mygene as mg
mg_client = mg.MyGeneInfo()
gene_info = mg_client.querymany(adata.var_names, scopes='ensembl.gene', fields='symbol', species='human')
# Create a mapping from Ensembl IDs to gene symbols
ensembl_to_symbol = {item['query']: item.get('symbol', None) for item in gene_info if 'query' in item and 'symbol' in item}
# Map the gene names
adata.var_names = [ensembl_to_symbol.get(gene, gene) for gene in adata.var_names]
# Remove genes without a symbol
adata = adata[:, adata.var_names.notna()].copy()
# Save the AnnData object
adata.write_h5ad('data.h5ad')
"