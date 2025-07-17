import os
import scanpy as sc
import pandas as pd
import numpy as np
import scib
import sklearn
import time
from scib_metrics.benchmark import Benchmarker, BioConservation, BatchCorrection
os.chdir('/gpfs/scratch/nk4167/KidneyAtlas/lake')
scrna_embeddings = '/gpfs/scratch/nk4167/KidneyAtlas/lake/scrna_embeddings/scgpt.h5ad'
adata = sc.read_h5ad(scrna_embeddings)
# Subsample 20k for speed
adata = sc.pp.subsample(adata, n_obs=20000, copy=True)
# For those specific cells, load the scimilarity embeddings
scimilarity_embeddings = '/gpfs/scratch/nk4167/KidneyAtlas/lake/scrna_embeddings/scimilarity.h5ad'
adata.obsm["X_scimilarity"] = sc.read_h5ad(scimilarity_embeddings)[adata.obs_names,:].obsm["X_scimilarity"]
# Same for geneformer
geneformer_embeddings = '/gpfs/scratch/nk4167/KidneyAtlas/lake/scrna_embeddings/geneformer.h5ad'
adata.obsm["X_Geneformer"] = sc.read_h5ad(geneformer_embeddings)[adata.obs_names,:].obsm["X_Geneformer"]
# Do PCA on X_scimilarity embeddings with sklearn
print("Running PCA on scimilarity embeddings...")
from sklearn.decomposition import PCA
pca = PCA(n_components=50)
adata.obsm['X_scimilarity_PCA'] = pca.fit_transform(adata.obsm['X_scimilarity'])
# Do neighbors on PCA
print("Computing neighbors on PCA...")
sc.pp.neighbors(adata, n_neighbors=15, use_rep='X_scimilarity_PCA')
# Save new UMAP in other slot
sc.tl.umap(adata, neighbors_key='neighbors', key_added='X_umap_scimilarity')
# Do PCA on X_Geneformer embeddings with sklearn
print("Running PCA on Geneformer embeddings...")
pca = PCA(n_components=50)
adata.obsm['X_Geneformer_PCA'] = pca.fit_transform(adata.obsm['X_Geneformer'])
# Do neighbors on PCA
print("Computing neighbors on PCA...")
sc.pp.neighbors(adata, n_neighbors=15, use_rep='X_Geneformer_PCA')
# Save new UMAP in other slot
sc.tl.umap(adata, neighbors_key='neighbors', key_added='X_umap_Geneformer')  
# Do PCA on X_scGPT embeddings with sklearn
print("Running PCA on scGPT embeddings...")
from sklearn.decomposition import PCA
pca = PCA(n_components=50)
adata.obsm['X_scGPT_PCA'] = pca.fit_transform(adata.obsm['X_scGPT'])
# Do neighbors on PCA
print("Computing neighbors on PCA...")
sc.pp.neighbors(adata, n_neighbors=15, use_rep='X_scGPT_PCA')
# Save new UMAP in other slot
sc.tl.umap(adata, neighbors_key='neighbors', key_added='X_umap_scGPT')

sc.pp.highly_variable_genes(adata, n_top_genes=2000, subset=True)
sc.pp.scale(adata, max_value=10)
sc.pp.pca(adata, n_comps=50, use_highly_variable=True)

adata.obsm["Unintegrated"] = adata.obsm["X_pca"]  # For scib compatibility
# Run harmony for quick comparison
print("Running Harmony...")
sc.external.pp.harmony_integrate(adata, key='library_id', basis='X_pca', max_iter_harmony=20, verbose=False)
# Quick scib evaluation
print("Running scib evaluation...")
import time



biocons = BioConservation(isolated_labels=False)

start = time.time()
bm = Benchmarker(
    adata,
    batch_key="library_id",
    label_key="cell_type",
    embedding_obsm_keys=["Unintegrated",'X_pca_harmony', "X_scGPT_PCA",'X_scimilarity_PCA', "X_Geneformer_PCA"],
    pre_integrated_embedding_obsm_key="X_pca",
    bio_conservation_metrics=biocons,
    batch_correction_metrics=BatchCorrection(),
    n_jobs=-1,
)
bm.prepare()
bm.benchmark()
end = time.time()
print(f"Time: {int((end - start) / 60)} min {int((end - start) % 60)} sec")
# Save this plot
bm.plot_results_table(show=False, min_max_scale=False,save_dir='figures')

# Plot UMAPs side by side and save
import matplotlib.pyplot as plt
sc.pl.embedding(adata, basis='X_umap_scGPT', color=['cell_type'], save='_scGPT.png', show=False)
plt.close()
# Same for regular UMAP (already exists in object)
sc.pl.embedding(adata, basis='X_umap', color=['cell_type'], save='_regular.png', show=False)
plt.close()
# Same for scimilarity UMAP
sc.pl.embedding(adata, basis='X_umap_scimilarity', color=['cell_type'], save='_scimilarity.png', show=False)
plt.close()
# Same for Geneformer UMAP
sc.pl.embedding(adata, basis='X_umap_Geneformer', color=['cell_type'], save='_Geneformer.png', show=False)
plt.close()
# Upload both pngs with rclone to GDrive (case sensitive)
import subprocess
subprocess.run(['rclone', 'copy', 'figures/X_umap_scGPT_scGPT.png', 'GDrive:KidneyAtlas/lake/umaps/'])
subprocess.run(['rclone', 'copy', 'figures/X_umap_regular.png', 'GDrive:KidneyAtlas/lake/umaps/'])
subprocess.run(['rclone', 'copy', 'figures/X_umap_scimilarity_scimilarity.png', 'GDrive:KidneyAtlas/lake/umaps/'])
subprocess.run(['rclone', 'copy', 'figures/X_umap_Geneformer_Geneformer.png', 'GDrive:KidneyAtlas/lake/umaps/'])
#