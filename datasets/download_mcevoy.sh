BASE_DIR=/gpfs/scratch/nk4167/KidneyAtlas
cd "$BASE_DIR"
OUTPUT_DIR="./McEvoy"
mkdir -p "$OUTPUT_DIR"
cd "$OUTPUT_DIR"
expr_matrix_link=https://cells.ucsc.edu/living-donor-kidney/exprMatrix.tsv.gz
curl -L "$expr_matrix_link" -o exprMatrix.tsv.gz
metadata_link=https://cells.ucsc.edu/living-donor-kidney/meta.tsv
curl -L "$metadata_link" -o meta.tsv
umap_coords_link=https://cells.ucsc.edu/living-donor-kidney/UMAP.coords.tsv.gz
curl -L "$umap_coords_link" -o UMAP.coords.tsv.gz