BASE_DIR=/gpfs/scratch/nk4167/KidneyAtlas
cd "$BASE_DIR"
PROJECT_ID="SCP1288"
mkdir -p "$PROJECT_ID"
cd $PROJECT_ID
curl -k "https://singlecell.broadinstitute.org/single_cell/api/v1/bulk_download/generate_curl_config?accessions=SCP1288&auth_code=XLUWmoLs&directory=all&context=study"  -o cfg.txt; curl -K cfg.txt && rm cfg.txt