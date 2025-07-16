# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

This repository evaluates popular methods for atlas integration of single-cell RNA-seq data, comparing them against zero-shot scRNA-seq foundation models. The project is structured around three main components:

1. **Datasets**: Collections from multiple disease conditions + healthy samples across 3 tissues (Breast, Kidney, Endometrium)
2. **Conventional Methods**: Seurat Sketch, Symphony, scVI, scPoli, CellHint
3. **Foundation Model Methods**: scGPT, scFoundation, CellFM, Geneformer

## Directory Structure

```
environments/          # Conda environment setup scripts for BigPurple HPC
demo_scripts/          # Example scripts demonstrating foundation model usage
datasets/              # Data download and processing scripts (minimal deps: cURL, scanpy)
fm_scripts/            # Foundation model embedding generation scripts (zero-shot)
method_scripts/        # Conventional method embedding generation scripts
```

## Environment Setup Commands

The project uses isolated conda environments for each foundation model. All environments target BigPurple HPC with path prefix `/gpfs/scratch/nk4167/miniconda/envs/`:

### scGPT Environment (Most Comprehensive)
```bash
# Primary method - uses comprehensive YAML configuration
conda env create -f environments/scgpt_env.yml
conda activate /gpfs/scratch/nk4167/miniconda/envs/scgpt_env

# Alternative - automated installer with GPU detection
bash environments/install_sgpt.sh

# macOS variant (CPU-only)
bash environments/scgpt_env_macos.sh
```

### Foundation Model Environments
```bash
# scFoundation (includes model downloads from SharePoint)
bash environments/scFoundation_env.sh

# scimilarity (downloads from Zenodo)
bash environments/scimilarity_env.sh

# Geneformer (clones from Hugging Face)
bash environments/geneformer_env.sh

# CellFM (uses MindSpore framework)
bash environments/cellFM_env.sh
```

## Key Dependencies by Environment

### scGPT (Python 3.9)
- **PyTorch 2.x** with CUDA 11.7 support
- **R integration** (R-base 4.2.2, rpy2) for some dependencies
- **wandb** for experiment tracking and hyperparameter logging
- **flash-attn<1.0.5** for attention optimization
- **scgpt** package and full scientific stack (scanpy, anndata, scvi)

### scFoundation (Python 3.10)  
- **PyTorch 2.5.0** with CUDA 12.1
- **einops, local_attention** for transformer operations
- Model weights downloaded from SharePoint (requires authentication)

### scimilarity (Python 3.10)
- **scimilarity** package (lightweight installation)
- Model downloaded from Zenodo (model_v1.1.tar.gz)

### Geneformer (Python 3.10)
- **Git LFS** for large model file handling
- Full repository clone from Hugging Face (ctheodoris/Geneformer)

### CellFM (Python 3.9)
- **MindSpore 2.2.10** framework (instead of PyTorch)
- **scanpy 1.10, scib 1.1.5** for single-cell analysis

## Common Development Commands

### Environment Management
```bash
# List available environments
ls environments/

# Activate specific environment
conda activate /gpfs/scratch/nk4167/miniconda/envs/[env_name]

# Check environment setup
conda list
```

### Running Foundation Models
```bash
# scGPT demo with hyperparameter configuration and wandb logging
python demo_scripts/scgpt_demo.py

# Check model file locations after environment setup
ls /gpfs/scratch/nk4167/miniconda/envs/*/models/
ls /gpfs/scratch/nk4167/miniconda/envs/*/model_v1.1/
```

### Development Workflow
```bash
# Basic dependency check (works in all environments)
python -c "import scanpy as sc; import anndata; print('Core deps OK')"

# Check CUDA availability for torch-based models
python -c "import torch; print(f'CUDA available: {torch.cuda.is_available()}')"
```

## Architecture and Model Integration

### Core Design Pattern
Each foundation model uses isolated environments with specific dependency chains:
- **Self-contained setup scripts** handle model downloads and dependency installation
- **Zero-shot embedding generation** workflow for comparing models
- **AnnData objects** as standard data format across all models

### Key Implementation Files
- **demo_scripts/scgpt_demo.py:712**: Complete training pipeline with hyperparameter defaults, data preprocessing, model training, and evaluation metrics using wandb
- **environments/install_sgpt.sh**: 360-line automated installer with GPU detection and error handling
- Environment scripts contain model-specific download procedures (SharePoint, Zenodo, Hugging Face, Google Drive)

### Model Download Strategies
- **scFoundation**: Authenticated SharePoint downloads (models.ckpt, models1.ckpt)
- **scimilarity**: Direct Zenodo download with tar extraction
- **Geneformer**: Git LFS clone from Hugging Face repository
- **scGPT**: Google Drive downloads for macOS, package installation for Linux

### Data Processing Architecture  
- Uses standard **scanpy/anndata** workflows for single-cell RNA-seq analysis
- Rapids GPU acceleration supported for faster processing when available
- No standard testing framework - research-focused development approach