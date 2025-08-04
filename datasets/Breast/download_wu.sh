BASE_DIR=/gpfs/scratch/nk4167/BreastAtlas
cd "$BASE_DIR"
OUTPUT_DIR="./Wu"
mkdir -p "$OUTPUT_DIR"
cd "$OUTPUT_DIR"
GEO_ACCESSION=GSE176078
Rscript -e "
library(GEOquery)
options(timeout = 9999999)
options(download.file.method.GEOquery = 'curl')
getGEOSuppFiles('$GEO_ACCESSION', makeDirectory = FALSE)
untar('GSE176078_RAW.tar') # This part fails but we don't care"
tar -zxvf GSE176078_Wu_etal_2021_BRCA_scRNASeq.tar.gz

# Process these into anndata
if [ ! -f "wu.h5ad" ]; then
    python -c "
import scanpy as sc
import pandas as pd
import os
import scipy.io
from pathlib import Path
import numpy as np
# Load the gene expression data
#In [2]: os.listdir('Wu_etal_2021_BRCA_scRNASeq')
#Out[2]: 
#['count_matrix_genes.tsv',
# 'metadata.csv',
# 'count_matrix_sparse.mtx',
# 'count_matrix_barcodes.tsv']
rna_counts = scipy.io.mmread('Wu_etal_2021_BRCA_scRNASeq/count_matrix_sparse.mtx').T
genenames = pd.read_csv('Wu_etal_2021_BRCA_scRNASeq/count_matrix_genes.tsv', header=None, sep='\t')
cellnames = pd.read_csv('Wu_etal_2021_BRCA_scRNASeq/count_matrix_barcodes.tsv', header=None, sep='\t')
rna_counts = pd.DataFrame(rna_counts.todense(), index=cellnames[0], columns=genenames[0])
metadata = pd.read_csv('Wu_etal_2021_BRCA_scRNASeq/metadata.csv', index_col=0)
adata = sc.AnnData(X=rna_counts.values,
                   obs = metadata,
                    var=pd.DataFrame(index=rna_counts.columns))

print( adata.X.max())
# Set var_names to gene names
adata.var_names = adata.var.index.tolist()
# Set obs_names to cell IDs
adata.obs_names = adata.obs.index.tolist()
# Ensure var_names are unique
adata.var_names_make_unique()
# Save the AnnData object
adata.write_h5ad('data.h5ad')
"