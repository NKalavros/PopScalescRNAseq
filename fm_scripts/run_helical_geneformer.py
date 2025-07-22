from helical.models.geneformer import Geneformer, GeneformerConfig
import anndata as ad
import os

os.chdir('/gpfs/scratch/nk4167/KidneyAtlas/lake')
# Example configuration
model_config = GeneformerConfig(model_name="gf-12L-95M-i4096", batch_size=10)
geneformer_v2 = Geneformer(model_config)

# Example usage for base pretrained model
ann_data = ad.read_h5ad("lake_scrna_15_genenames.h5ad")
dataset = geneformer_v2.process_data(ann_data)
embeddings = geneformer_v2.get_embeddings(dataset)
print("Base model embeddings shape:", embeddings.shape)
x