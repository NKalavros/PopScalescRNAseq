BASE_DIR=/gpfs/scratch/nk4167/BreastAtlas
cd "$BASE_DIR"
OUTPUT_DIR="./Gray"
mkdir -p "$OUTPUT_DIR"
cd "$OUTPUT_DIR"
GEO_ACCESSION=GSE180878
Rscript -e "
library(GEOquery)
options(timeout = 9999999)
options(download.file.method.GEOquery = 'curl')
getGEOSuppFiles('$GEO_ACCESSION', makeDirectory = FALSE)
untar('GSE180878_RAW.tar') # This part fails but we don't care"

#drwxrwxr-x 5 nk4167 nk4167 4.0K Aug  2 14:58 ..
#-rw-rw-r-- 1 nk4167 nk4167 103M Aug  2 14:59 GSE180878_Li_Brugge_10XscRNAseq_GeneCellMatrix_RNAcounts_human.csv.gz
#drwxrwxr-x 2 nk4167 nk4167 4.0K Aug  2 14:59 .
#-rw-rw-r-- 1 nk4167 nk4167 313K Aug  2 14:59 GSE180878_Li_Brugge_10XscRNAseq_Metadata_human.csv.gz

# Process these into anndata
if [ ! -f "gray.h5ad" ]; then
    python -c "
import scanpy as sc
import pandas as pd
import os
from pathlib import Path
import numpy as np
# Load the gene expression data
rna_counts = pd.read_csv('GSE180878_Li_Brugge_10XscRNAseq_GeneCellMatrix_RNAcounts_human.csv.gz', 
                          index_col=0, 
                          compression='gzip')
rna_counts = rna_counts.T  # Transpose to have genes as rows and cells as columns
# Load the metadata
metadata = pd.read_csv('GSE180878_Li_Brugge_10XscRNAseq_Metadata_human.csv.gz', 
                       index_col=0, 
                       compression='gzip')
# Create an AnnData object
adata = sc.AnnData(X=rna_counts.values,
                   obs=metadata,
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
