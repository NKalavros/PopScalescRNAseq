# Context Findings

## CUDA 12 Compatibility Analysis
- **scGPT**: Currently uses CUDA 11.7, needs update to CUDA 12.1
- **scFoundation, Geneformer, CellFM**: Already use CUDA 12.1 
- **Standard pattern**: `pip3 install torch torchvision torchaudio --index-url https://download.pytorch.org/whl/cu121`

## H5AD Integration Patterns
- **scGPT demo**: Uses scvi.data.pbmc_dataset() as reference
- **Preprocessing**: Standardized through Preprocessor class with binning
- **Output format**: `adata.obsm["X_scGPT"]` with L2 normalization

## Environment Structure
- **Prefix**: `/gpfs/scratch/nk4167/miniconda/envs/`
- **Pattern**: `{model_name}_env` naming convention
- **Activation**: Full path activation required

## Model Download Verification Needed
- **scFoundation**: SharePoint downloads (models.ckpt, models1.ckpt)
- **scimilarity**: Zenodo tar.gz extraction to model_v1.1/
- **Geneformer**: HuggingFace repository clone
- **CellFM**: No explicit model downloads in script

## Implementation Gaps Identified
- Missing fm_scripts/ directory for embedding generation
- Need individual demo scripts per model
- Require automated testing framework
- Need standardized embedding output format

## User Preference Update
- Test environments one at a time with manual output review
- User will run commands and provide output for validation
- Iterative approach preferred over fully automated