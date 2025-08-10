# You can do this in base, it needs very few packages
pip install hydra-core
pip install transcriptformer --no-deps
BASE_DIR='/gpfs/scratch/nk4167'
cd $BASE_DIR
transcriptformer download tf-sapiens

cd /gpfs/scratch/nk4167/KidneyAtlas/SCP1288

python -c "
import sys
import mygene
import anndata as ad
adata_path = sys.argv[1] if len(sys.argv) > 1 else 'data.h5ad'
adata = ad.read_h5ad(adata_path)
mg = mygene.MyGeneInfo()

# Get gene symbols from the current var_names
gene_symbols = adata.var_names.tolist()

# Query mygene to get Ensembl IDs
result = mg.querymany(gene_symbols, scopes='symbol', fields='ensembl.gene', species='human')

# Create a mapping dictionary
ensembl_mapping = {}
for item in result:
    if 'ensembl' in item and 'gene' in item['ensembl']:
        ensembl_id = item['ensembl']['gene']
        if isinstance(ensembl_id, list):
            ensembl_id = ensembl_id[0]  # Take first if multiple
        ensembl_mapping[item['query']] = ensembl_id
    else:
        ensembl_mapping[item['query']] = None

# Add ensembl_id column to adata.var
adata.var['ensembl_id'] = [ensembl_mapping.get(gene, None) for gene in adata.var_names]
adata.var_names_make_unique()  # Ensure unique var_names
adata.obs['assay'] = 'SCP1288'  # Add assay column
adata.write_h5ad('data_with_ensembl.h5ad')"


transcriptformer inference \
  --checkpoint-path /gpfs/scratch/nk4167/checkpoints/tf_sapiens \
  --data-file  data_with_ensembl.h5ad \
  --output-path ./embeddings \
  --batch-size 16 \
  --remove-duplicate-genes

mv embeddings/embeddings.h5ad embeddings/transcriptformer.h5ad