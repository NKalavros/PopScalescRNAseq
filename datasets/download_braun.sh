#NOTE THIS ONE I ACTUALLY DID MYSELF
BASE_DIR=/gpfs/scratch/nk4167/KidneyAtlas
cd "$BASE_DIR"
OUTPUT_DIR="./Braun"
mkdir -p "$OUTPUT_DIR"
cd "$OUTPUT_DIR"
curl -O https://pmc.ncbi.nlm.nih.gov/articles/instance/8138872/bin/NIHMS1692222-supplement-supplementary_Data_S1.zip
unzip NIHMS1692222-supplement-supplementary_Data_S1.zip
curl -O https://pmc.ncbi.nlm.nih.gov/articles/instance/8138872/bin/NIHMS1692222-supplement-supplementary_Data_S2.csv
curl -O NIHMS1692222-supplement-Table_S3.xlsx # annotation xlsx file
# Let's make an h5ad file from these
python -c "
import datatable as dt
import pandas as pd
import anndata as ad
import scanpy as sc

# Load the data
df1 = dt.fread('Data S1.csv').to_pandas()
df2 = dt.fread('NIHMS1692222-supplement-supplementary_Data_S2.csv',fill=True).to_pandas()
# Set first column as index and drop it
df1.set_index(df1.columns[0], inplace=True)
df2.set_index(df2.columns[0], inplace=True)
# Load the annotation file. There is a header row, so we need to make it colnames
annotations = pd.read_excel('NIHMS1692222-supplement-Table_S2.xlsx', header=1)
# Set the annotation index to the first column, no need to drop it
annotations.set_index(annotations.columns[0], inplace=True)
# Concatenate df1 and df2 (counts)
df1 = pd.concat([df1, df2], axis=1)
# Transpose from R to Python style
df1 = df1.T
# Turn nans to 0
df1.fillna(0, inplace=True)
# Ensure they are the same
if not df1.index.equals(annotations.index):
    raise ValueError('Index of counts data does not match index of annotations')
# Convert object columns to strings to avoid h5ad writing issues
for col in annotations.columns:
    if annotations[col].dtype == 'object':
        annotations[col] = annotations[col].astype(str)
# Create an AnnData object
adata = ad.AnnData(X=df1, obs=annotations)
sc.pp.filter_cells(adata, min_genes=200)  # Filter cells with at least 200 genes
sc.pp.filter_cells(adata, min_counts=400)  # Filter cells with at least 200 genes
sc.pp.filter_genes(adata, min_cells=3)  # Filter genes expressed in
# Save the AnnData object
adata.write_h5ad('data.h5ad')
"