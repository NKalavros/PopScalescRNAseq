# Download the Breast dataset from Bassez et al.
# Note that you need to login for that.
BASE_DIR=/gpfs/scratch/nk4167/BreastAtlas
OUTPUT_DIR="./Bassez"

echo "Changing to base directory: $BASE_DIR"
cd "$BASE_DIR" || { echo "Error: Could not change to $BASE_DIR"; exit 1; }

echo "Creating output directory: $OUTPUT_DIR"
mkdir -p "$OUTPUT_DIR"
cd "$OUTPUT_DIR" || { echo "Error: Could not change to $OUTPUT_DIR"; exit 1; }


echo "Downloading counts_cells_cohort1.rds..."
curl -L --retry 5 --retry-delay 10 --continue-at - \
  -A "Mozilla/5.0 (Macintosh; Intel Mac OS X 10.15; rv:140.0) Gecko/20100101 Firefox/140.0" \
  -H 'Accept: application/json, text/plain, */*' \
  -H 'Accept-Language: en-US,en;q=0.5' \
  -H 'Accept-Encoding: gzip, deflate, br, zstd' \
  -H 'ContentType: application/json' \
  -H 'Origin: https://lambrechtslab.sites.vib.be' \
  -H 'Connection: keep-alive' \
  -H 'Referer: https://lambrechtslab.sites.vib.be/' \
  -H 'Sec-Fetch-Dest: empty' \
  -H 'Sec-Fetch-Mode: cors' \
  -H 'Sec-Fetch-Site: same-site' \
  -H 'Priority: u=0' \
  -H 'TE: trailers' \
  -o counts_cells_cohort1.rds \
  https://forms-api.vib.be/api/v1/files/1863

echo "Downloading metadata1.csv..."
curl -L --retry 5 --retry-delay 10 --continue-at - \
  -A "Mozilla/5.0 (Macintosh; Intel Mac OS X 10.15; rv:140.0) Gecko/20100101 Firefox/140.0" \
  -H 'Accept: application/json, text/plain, */*' \
  -H 'Accept-Language: en-US,en;q=0.5' \
  -H 'Accept-Encoding: gzip, deflate, br, zstd' \
  -H 'ContentType: application/json' \
  -H 'Origin: https://lambrechtslab.sites.vib.be' \
  -H 'Connection: keep-alive' \
  -H 'Referer: https://lambrechtslab.sites.vib.be/' \
  -H 'Sec-Fetch-Dest: empty' \
  -H 'Sec-Fetch-Mode: cors' \
  -H 'Sec-Fetch-Site: same-site' \
  -H 'Priority: u=0' \
  -H 'TE: trailers' \
  -o metadata1.csv \
  https://forms-api.vib.be/api/v1/files/1872

echo "Downloading counts_cells_cohort2.rds..."
curl -L --retry 5 --retry-delay 10 --continue-at - \
  -A "Mozilla/5.0 (Macintosh; Intel Mac OS X 10.15; rv:140.0) Gecko/20100101 Firefox/140.0" \
  -H 'Accept: application/json, text/plain, */*' \
  -H 'Accept-Language: en-US,en;q=0.5' \
  -H 'Accept-Encoding: gzip, deflate, br, zstd' \
  -H 'ContentType: application/json' \
  -H 'Origin: https://lambrechtslab.sites.vib.be' \
  -H 'Connection: keep-alive' \
  -H 'Referer: https://lambrechtslab.sites.vib.be/' \
  -H 'Sec-Fetch-Dest: empty' \
  -H 'Sec-Fetch-Mode: cors' \
  -H 'Sec-Fetch-Site: same-site' \
  -H 'Priority: u=0' \
  -H 'TE: trailers' \
  -o counts_cells_cohort2.rds \
  https://forms-api.vib.be/api/v1/files/1867

echo "Downloading metadata2.csv..."
curl -L --retry 5 --retry-delay 10 --continue-at - \
  -A "Mozilla/5.0 (Macintosh; Intel Mac OS X 10.15; rv:140.0) Gecko/20100101 Firefox/140.0" \
  -H 'Accept: application/json, text/plain, */*' \
  -H 'Accept-Language: en-US,en;q=0.5' \
  -H 'Accept-Encoding: gzip, deflate, br, zstd' \
  -H 'ContentType: application/json' \
  -H 'Origin: https://lambrechtslab.sites.vib.be' \
  -H 'Connection: keep-alive' \
  -H 'Referer: https://lambrechtslab.sites.vib.be/' \
  -H 'Sec-Fetch-Dest: empty' \
  -H 'Sec-Fetch-Mode: cors' \
  -H 'Sec-Fetch-Site: same-site' \
  -H 'Priority: u=0' \
  -H 'TE: trailers' \
  -o metadata2.csv \
  https://forms-api.vib.be/api/v1/files/1871


echo "\n--- DEBUGGING CURL REQUEST ---"
curl -v -w "\nHTTP status: %{http_code}\n" 'https://forms-api.vib.be/api/v1/files/1871' \
  -H 'User-Agent: Mozilla/5.0 (Macintosh; Intel Mac OS X 10.15; rv:140.0) Gecko/20100101 Firefox/140.0' \
  -H 'Accept: application/json, text/plain, */*' \
  -H 'Accept-Language: en-US,en;q=0.5' \
  -H 'Accept-Encoding: gzip, deflate, br, zstd' \
  -H 'ContentType: application/json' \
  -H 'Authorization: Bearer eyJhbGciOiJSUzI1NiIsImtpZCI6IjcxQzM4RjNDRkU3MUU1REUzMDZDMDg5NjMzNkQwRjRDMkU5NjNFQzdSUzI1NiIsIng1dCI6ImNjT1BQUDV4NWQ0d2JBaVdNMjBQVEM2V1BzYyIsInR5cCI6ImF0K2p3dCJ9.eyJpc3MiOiJodHRwczovL3NlcnZpY2VzLnZpYi5iZSIsIm5iZiI6MTc1NDE4NTYyNiwiaWF0IjoxNzU0MTg1NjI2LCJleHAiOjE3NTQxODkyMjYsImF1ZCI6ImFwaS5mb3JtcyIsInNjb3BlIjpbIm9wZW5pZCIsInByb2ZpbGUiLCJ1c2Vycm9sZXMiLCJhcGkuZm9ybXMiXSwiYW1yIjpbImV4dGVybmFsIl0sImNsaWVudF9pZCI6InZpYi1mb3Jtcy1qcyIsInN1YiI6IjMwNTUzIiwiYXV0aF90aW1lIjoxNzU0MTYwNDQ4LCJpZHAiOiJHb29nbGUiLCJlbWFpbCI6ImthbGF2cm9zbmlrb2xhb3NAZ21haWwuY29tIiwibmFtZSI6Ik5pa29sYXMiLCJmYW1pbHlfbmFtZSI6IkthbGF2cm9zIiwicGkiOiJmYWxzZSIsIm5pY2tOYW1lIjoiIiwidXBuIjoia2FsYXZyb3NuaWtvbGFvc0BnbWFpbC5jb20iLCJ2aWJfY2VudGVycyI6IiIsInZpYl9sYWJzIjoiIiwicm9sZSI6IkluZHVzdHJ5Iiwic2lkIjoiOUM1NTM1MDIxMjY2MkRDM0RBNTEwOEVBRDI0Q0MyQUMifQ.OJ_oGnZ5J1orimfnLHmZVxu74-uiHeWOiWfgcijoqwjyHNZrciaT3Z-c6_hFAbqWqlMZwsLl-ZZZm4F5EECEEiprGRYgV1w6vdy2fJeVrtfRmYMRBgrR_9hTty9Og_g0v-nb9HdsHsYrMKtrfLpJJ0GuTcd9vENplnG1z9cF4dz0lqo5zZeIVtqauSRVaqJgpMwhePYveNNXgsfVbsNyOjnGUmGVJ1HK1vECPb_6g5uNMskezSAhYcVAhqS5Q0RAK49OeSqaBo0rADBeQgcIGuDeoDlxrnf18dQXo9-XwexPwkyNUdS8ArULPBZIxW9gPaabRKwIrBgoFFB92AJFMIuTUum6WurIo5vQ1RoOllY_uNqamZaceRKbNEnLfLlD6WiiBR8YShIMM_SUgBSKf59lLVKoRXuxinF6dbZwTcV6W78Hyh4oT-q9-X0bywS-UKXMerpD8pasflLbOGACiK2SK9aiQnDlqhuzlTix57Ky7Q_W6fj31EPcI5egP8tghFCZMCLmavVbe1ZE6fslN_PTneNCH7nRdklpaJ2RAgqYoNV12gB5_ZiRweWXXjUjThq87a3Ae5a2ZCVaVxGXcCa4P-5y5456ZMKsblyXpAqrwkREAQf-OS-KwZ9hPTZ70UKnOOYgJiA-uhoMVaDwzRGx4Lj60IXdZqKC7HmaCzU' \
  -H 'Origin: https://lambrechtslab.sites.vib.be' \
  -H 'Connection: keep-alive' \
  -H 'Referer: https://lambrechtslab.sites.vib.be/' \
  -H 'Sec-Fetch-Dest: empty' \
  -H 'Sec-Fetch-Mode: cors' \
  -H 'Sec-Fetch-Site: same-site' \
  -H 'Priority: u=0' \
  -H 'TE: trailers'

#1863-counts_cells_cohort1.rds  1871-BIOKEY_metaData_cohort2_web.xls
#1867-counts_cells_cohort2.rds  1872-BIOKEY_metaData_cohort1_web.xls

Rscript -e "
library(Seurat)
library(fastMatMR)
library(openxlsx)
library(readxl)
# Load the data
counts_cells_cohort1 <- readRDS('1863-counts_cells_cohort1.rds')
counts_cells_cohort2 <- readRDS('1867-counts_cells_cohort2.rds')
# Load the metadata
metadata2 <- read.csv('1871-BIOKEY_metaData_cohort2_web.xls')
metadata1 <- read.csv('1872-BIOKEY_metaData_cohort1_web.xls')
# Create Seurat objects
seurat1 <- CreateSeuratObject(counts = counts_cells_cohort1, meta.data = metadata1)
seurat2 <- CreateSeuratObject(counts = counts_cells_cohort2, meta.data = metadata2)
seurat1$Cohort = 'cohort1'
seurat2$Cohort = 'cohort2'
# Combine the Seurat objects
seur_obj <- Reduce(merge, list(seurat1, seurat2))
seur_obj <- JoinLayers(seur_obj)
write_fmm(LayerData(seur_obj,layer='counts',assay='RNA'), 'counts.mtx')
# Write rownames and colnames
DefaultAssay(seur_obj) <- 'RNA'
write.table(rownames(seur_obj), 'genes.txt', row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(colnames(seur_obj), 'barcodes.txt', row.names = FALSE, col.names = FALSE, quote = FALSE)
# Save metadata
write.csv(seur_obj@meta.data, 'metadata.csv')
saveRDS(seur_obj, 'seurat_obj.rds')"


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