import scanpy as sc
import os
panpdac_obj_path = '/gpfs/data/tsirigoslab/home/camerd03/panPDAC/analysis-min1200umis-min500genes-max15mito/out-combined/seurat_obj_counts.h5ad'
outdir = '/gpfs/scratch/nk4167/PancreasAtlas/'
# Create output directory if it doesn't exist
os.makedirs(outdir, exist_ok=True)
# Load the AnnData object
adata = sc.read_h5ad(panpdac_obj_path)
print(adata.obs['Study'].value_counts())
# Split the AnnData object by 'Study'
study_groups = adata.obs['Study'].unique()
for study in study_groups:
    print(f'Processing study: {study}')
    adata_subset = adata[adata.obs['Study'] == study].copy()
    # Make a subdirectory by study name
    study_dir = os.path.join(outdir, study)
    os.makedirs(study_dir, exist_ok=True)
    # Save each subset to a separate .h5ad file named data.h5ad
    subset_path = os.path.join(study_dir, 'data.h5ad')
    adata_subset.write_h5ad(subset_path)
    print(f'Saved {study} subset with {adata_subset.n_obs} cells to {subset_path}')
    # PRint the sameple counts from orig.ident
    print(adata_subset.obs['orig.ident'].value_counts())


adata = sc.read_h5ad('/gpfs/data/tsirigoslab/home/camerd03/panPDAC/analysis-min1200umis-min500genes-max15mito/out-combined/seurat_obj_scANVI_5k.h5ad')
# Save the metadata
adata.obs.to_csv(os.path.join('/gpfs/scratch/nk4167/PancreasAtlas/metadata.csv'))