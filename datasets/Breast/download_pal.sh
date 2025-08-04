BASE_DIR=/gpfs/scratch/nk4167/BreastAtlas
cd "$BASE_DIR"
OUTPUT_DIR="./Pal"
mkdir -p "$OUTPUT_DIR"
cd "$OUTPUT_DIR"
GEO_ACCESSION=GSE161529
Rscript -e "
library(GEOquery)
options(timeout = 9999999)
options(download.file.method.GEOquery = 'curl')
getGEOSuppFiles('$GEO_ACCESSION', makeDirectory = FALSE)
untar('GSE161529_RAW.tar') # This part fails but we don't care"

# Download from figshare https://figshare.com/ndownloader/articles/17058077/versions/1
curl 'https://figshare.com/ndownloader/articles/17058077/versions/1' \
  -H 'User-Agent: Mozilla/5.0 (Macintosh; Intel Mac OS X 10.15; rv:140.0) Gecko/20100101 Firefox/140.0' \
  -H 'Accept: text/html,application/xhtml+xml,application/xml;q=0.9,*/*;q=0.8' \
  -H 'Accept-Language: en-US,en;q=0.5' \
  -H 'Accept-Encoding: gzip, deflate, br, zstd' \
  -H 'Referer: https://figshare.com/' \
  -H 'Connection: keep-alive' \
  -H 'Cookie: fig_tracker_client=df4b334a-f924-4a91-b40e-71670531b555; GLOBAL_FIGSHARE_SESSION_KEY=82bfad5303938edd1934b2bd1ed744cd661ce626515b29fd47544b08b83e41047d10956d; FIGINSTWEBIDCD=82bfad5303938edd1934b2bd1ed744cd661ce626515b29fd47544b08b83e41047d10956d; figshare-cookies-essential=true; figshare-cookies-performance=true' \
  -H 'Upgrade-Insecure-Requests: 1' \
  -H 'Sec-Fetch-Dest: document' \
  -H 'Sec-Fetch-Mode: navigate' \
  -H 'Sec-Fetch-Site: same-origin' \
  -H 'Sec-Fetch-User: ?1' \
  -H 'Priority: u=0, i' \
  -H 'TE: trailers' \
  --output pal.zip

export  UNZIP_DISABLE_ZIPBOMB_DETECTION=TRUE
unzip -o pal.zip
unzip HumanBreast10X.zip
# Read all objects in the directory (RDS)
Rscript -e "
library(Seurat)
library(dplyr)
library(Matrix)
library(fastMatMR)
# Iterate over RDS files in the directory
rds_files <- list.files(pattern = '\\.rds$', full.names = TRUE)
objects = lapply(rds_files, readRDS)
lapply(objects,class)
for(i in seq_along(objects)) {
    try({
        objects[[i]] <- UpdateSeuratObject(objects[[i]])
    })
}
# Keep only the ones that class are  'Seurat'
objects <- objects[sapply(objects, function(x) class(x) == 'Seurat')]
# Calculate total number of cells
total_cells <- sum(sapply(objects, function(x) ncol(x)))
# Calculate total number of genes
total_genes <- sum(sapply(objects, function(x) nrow(x)))
# Print the total number of cells and genes
cat('Total number of cells:', total_cells, '\\n')
cat('Total number of genes:', total_genes, '\\n')
# Rename cells using the filenames
for (i in seq_along(objects)) {
    objects[[i]] <- RenameCells(objects[[i]], new.names = paste(gsub('./','',rds_files[i],fixed=T), colnames(objects[[i]]),sep='_'))
}
# Merge all Seurat objects into 
seur_obj <- Reduce(merge, objects)
write_fmm(LayerData(seur_obj,layer='counts',assay='RNA'), 'counts.mtx')
# Write rownames and colnames
DefaultAssay(seur_obj) <- 'RNA'
write.table(rownames(seur_obj), 'genes.txt', row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(colnames(seur_obj), 'barcodes.txt', row.names = FALSE, col.names = FALSE, quote = FALSE)
# Save metadata
write.csv(seur_obj@meta.data, 'metadata.csv')
saveRDS(seur_obj, 'seurat_obj.rds')
"

# Now make it into a Python h5ad file
python -c "
import scanpy as sc
import h5py
import numpy as np
import pandas as pd
import os
import sys
import matplotlib.pyplot as plt
import seaborn as sns
import scipy.sparse as sp
import scipy.io
import scipy
import fast_matrix_market as fmm

mtx = scipy.io.mmread('counts.mtx')
adata = sc.AnnData(X=mtx.T)
adata.obs_names = pd.read_csv('barcodes.txt', header=None).squeeze().tolist()
adata.var_names = pd.read_csv('genes.txt', header=None).squeeze().tolist()
adata.obs = pd.read_csv('metadata.csv', index_col=0)
adata.var_names_make_unique()
print(adata.X.max())
# Coo to CSR matrix
adata.X = sp.csr_matrix(adata.X)
adata.write_h5ad('data.h5ad')
"