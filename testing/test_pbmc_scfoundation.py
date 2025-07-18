import scanpy as sc
import subprocess
adata = sc.read_h5ad('scrna_embeddings/scfoundation.h5ad')
# Run PCA on the embeddings
print("Running PCA on scfoundation embeddings...")
from sklearn.decomposition import PCA
pca = PCA(n_components=50)
adata.obsm['X_scfoundation_PCA'] = pca.fit_transform(adata.obsm['scFoundation_embeddings'])
# Do neighbors on PCA
print("Computing neighbors on PCA...")
sc.pp.neighbors(adata, n_neighbors=15, use_rep='X_scfoundation_PCA')
# Save new UMAP in other slot
sc.tl.umap(adata)
adata.obsm['X_umap_scfoundation'] = adata.obsm['X_umap'].copy()
# Plot and upload UMAP
sc.pl.umap(adata, color=[ 'louvain'], save='_scfoundation_umap.png')
#RClone to GDrive
subprocess.run(['rclone', 'copy', 'figures/umap_scfoundation_umap.png', 'GDrive:KidneyAtlas/test_pbmc/'])