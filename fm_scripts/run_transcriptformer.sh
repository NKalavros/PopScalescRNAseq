#!/bin/bash
set -euo pipefail

DATA_FILE="${1:-data.h5ad}"
PROCESSED_FILE="data_with_ensembl.h5ad"

if [ ! -f "$DATA_FILE" ]; then
    echo "Input AnnData file '$DATA_FILE' not found in $(pwd)" >&2
    exit 1
fi

python - "$DATA_FILE" "$PROCESSED_FILE" <<'PYCODE'
import sys
import mygene
import anndata as ad

adata_path, output_path = sys.argv[1:3]
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
adata.obs['index'] = adata.obs_names.astype(str)  # Ensure obs_names are strings
adata.write_h5ad(output_path)
PYCODE

transcriptformer inference \
  --checkpoint-path /gpfs/scratch/nk4167/checkpoints/tf_sapiens \
  --data-file "$PROCESSED_FILE" \
  --output-path ./embeddings \
  --batch-size 16 \
  --remove-duplicate-genes

mv embeddings/embeddings.h5ad embeddings/transcriptformer.h5ad
