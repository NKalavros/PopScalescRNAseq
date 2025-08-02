BASE_DIR=/gpfs/scratch/nk4167/KidneyAtlas
cd "$BASE_DIR"
OUTPUT_DIR="./Abedini"
mkdir -p "$OUTPUT_DIR"
cd "$OUTPUT_DIR"
GEO_ACCESSION=GSE211785
Rscript -e "
library(GEOquery)
options(timeout = 9999999)
options(download.file.method.GEOquery = 'curl')
getGEOSuppFiles('$GEO_ACCESSION', makeDirectory = FALSE)
untar('GSE211785_RAW.tar')"
# Set .raw.X as the default matrix
python -c "
import anndata
adata = anndata.read_h5ad('GSE211785_Susztak_SC_SN_ATAC_merged_PostSCVI_final.h5ad')
# Keep the obs though!
adata = anndata.AnnData(X=adata.raw.X, obs=adata.obs, var=adata.raw.var)
sc.pp.filter_genes(adata, min_cells=3)
sc.pp.filter_cells(adata, min_genes=200)
sc.pp.filter_cells(adata, min_counts=400)
adata.write('GSE211785_Susztak_SC_SN_ATAC_merged_PostSCVI_final_counts.h5ad')
"