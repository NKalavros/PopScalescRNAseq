#!/bin/bash

# scGPT Installation Script with Isolated Conda Environment
# This script creates a complete conda environment for scGPT with all dependencies

set -e  # Exit on error

# Color codes for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

# Configuration variables
ENV_NAME="scgpt_env"
PYTHON_VERSION="3.9"
R_VERSION="4.2.2"
CUDA_VERSION="11.7"

# Function to print colored output
print_status() {
    echo -e "${GREEN}[INFO]${NC} $1"
}

print_error() {
    echo -e "${RED}[ERROR]${NC} $1"
}

print_warning() {
    echo -e "${YELLOW}[WARNING]${NC} $1"
}

# Function to check if conda is installed
check_conda() {
    if ! command -v conda &> /dev/null; then
        print_error "Conda is not installed. Please install Miniconda or Anaconda first."
        echo "Visit: https://docs.conda.io/en/latest/miniconda.html"
        exit 1
    fi
    print_status "Conda is installed: $(conda --version)"
}

# Function to check CUDA availability
check_cuda() {
    if command -v nvidia-smi &> /dev/null; then
        print_status "NVIDIA GPU detected"
        nvidia-smi --query-gpu=name --format=csv,noheader
        return 0
    else
        print_warning "No NVIDIA GPU detected. Will install CPU-only version."
        return 1
    fi
}

# Parse command line arguments
INSTALL_MODE="full"  # default to full installation
FORCE_REINSTALL=false

while [[ $# -gt 0 ]]; do
    case $1 in
        --basic)
            INSTALL_MODE="basic"
            shift
            ;;
        --env-name)
            ENV_NAME="$2"
            shift 2
            ;;
        --force)
            FORCE_REINSTALL=true
            shift
            ;;
        --help)
            echo "Usage: $0 [OPTIONS]"
            echo "Options:"
            echo "  --basic          Install basic version without flash-attn"
            echo "  --env-name NAME  Specify conda environment name (default: scgpt_env)"
            echo "  --force          Force reinstall if environment exists"
            echo "  --help           Show this help message"
            exit 0
            ;;
        *)
            print_error "Unknown option: $1"
            exit 1
            ;;
    esac
done

# Main installation process
print_status "Starting scGPT installation..."
print_status "Installation mode: $INSTALL_MODE"

# Check conda installation
check_conda

# Check if environment already exists
if conda env list | grep -q "^${ENV_NAME} "; then
    if [ "$FORCE_REINSTALL" = true ]; then
        print_warning "Environment '$ENV_NAME' exists. Removing and recreating..."
        conda deactivate 2>/dev/null || true
        conda env remove -n "$ENV_NAME" -y
    else
        print_error "Environment '$ENV_NAME' already exists. Use --force to reinstall."
        exit 1
    fi
fi

# Create conda environment configuration file
print_status "Creating environment configuration file..."

cat > scgpt_environment.yml << EOF
name: ${ENV_NAME}
channels:
  - pytorch
  - nvidia
  - conda-forge
  - bioconda
  - r
  - defaults
dependencies:
  # Python and core dependencies
  - python=${PYTHON_VERSION}
  - pip
  - numpy
  - scipy
  - pandas
  - matplotlib
  - seaborn
  - jupyter
  - ipykernel
  - faiss-cpu
  - louvain
  - wandb
  
  # R and R packages (required for some dependencies)
  - r-base=${R_VERSION}
  - r-essentials
  - r-devtools
  - rpy2
  
  # PyTorch and CUDA (if GPU available)
EOF

# Add PyTorch dependencies based on CUDA availability
if check_cuda; then
    cat >> scgpt_environment.yml << EOF
  - pytorch=*=*cuda${CUDA_VERSION}*
  - pytorch-cuda=${CUDA_VERSION}
  - cudatoolkit=${CUDA_VERSION}
EOF
else
    cat >> scgpt_environment.yml << EOF
  - pytorch
  - cpuonly
EOF
fi

# Add additional scientific computing packages
cat >> scgpt_environment.yml << EOF
  
  # Scientific computing
  - scikit-learn
  - scanpy
  - anndata
  - umap-learn
  - leidenalg
  - louvain
  
  # Development tools
  - git
  - gcc_linux-64  # Required for compiling some packages
  - gxx_linux-64
EOF

# Create the conda environment
print_status "Creating conda environment '$ENV_NAME'..."
conda env create -f scgpt_environment.yml

# Activate the environment
print_status "Activating environment..."
eval "$(conda shell.bash hook)"
conda activate "$ENV_NAME"

# Install additional pip packages
print_status "Installing scGPT and additional dependencies via pip..."

if [ "$INSTALL_MODE" = "full" ]; then
    # Full installation with flash-attn
    print_status "Installing full version with flash-attn..."
    
    # Install flash-attn first (version < 1.0.5 as recommended)
    pip install "flash-attn" --no-build-isolation
    
    # Install scGPT
    pip install scgpt
else
    # Basic installation without flash-attn
    print_status "Installing basic version without flash-attn..."
    pip install scgpt
fi

# Install additional recommended packages
pip install wandb  # For logging and visualization
pip install orbax==0.1.7  # Specific version to avoid compatibility issues

# Create a test script
print_status "Creating test script..."
cat > test_scgpt.py << 'EOF'
#!/usr/bin/env python
"""Test script to verify scGPT installation"""

import sys
print("Python version:", sys.version)

try:
    import torch
    print(f"PyTorch version: {torch.__version__}")
    print(f"CUDA available: {torch.cuda.is_available()}")
    if torch.cuda.is_available():
        print(f"CUDA version: {torch.version.cuda}")
        print(f"GPU: {torch.cuda.get_device_name(0)}")
except ImportError as e:
    print(f"Error importing PyTorch: {e}")

try:
    import scgpt
    print("scGPT imported successfully!")
    
    # Try to import key modules
    from scgpt.model import TransformerModel
    from scgpt.tokenizer import GeneVocab
    print("Key scGPT modules imported successfully!")
    
except ImportError as e:
    print(f"Error importing scGPT: {e}")

try:
    import scanpy as sc
    import anndata
    print("Scanpy and AnnData imported successfully!")
except ImportError as e:
    print(f"Error importing scanpy/anndata: {e}")

print("\nInstallation test completed!")
EOF

# Run the test script
print_status "Running installation test..."
python test_scgpt.py

# Create activation script
print_status "Creating activation script..."
cat > activate_scgpt.sh << EOF
#!/bin/bash
# Activation script for scGPT environment

conda activate ${ENV_NAME}
echo "scGPT environment activated!"
echo "Python: \$(python --version)"
echo "PyTorch: \$(python -c 'import torch; print(torch.__version__)')"
EOF
chmod +x activate_scgpt.sh

# Create example notebook
print_status "Creating example notebook..."
cat > scgpt_example.ipynb << 'EOF'
{
 "cells": [
  {
   "cell_type": "markdown",
   "metadata": {},
   "source": [
    "# scGPT Example Notebook\n",
    "This notebook demonstrates basic usage of scGPT"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": null,
   "metadata": {},
   "outputs": [],
   "source": [
    "# Import required libraries\n",
    "import scgpt\n",
    "import scanpy as sc\n",
    "import numpy as np\n",
    "import pandas as pd\n",
    "import torch\n",
    "\n",
    "print(f\"scGPT version: {scgpt.__version__ if hasattr(scgpt, '__version__') else 'Unknown'}\")\n",
    "print(f\"PyTorch version: {torch.__version__}\")\n",
    "print(f\"CUDA available: {torch.cuda.is_available()}\")"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": null,
   "metadata": {},
   "outputs": [],
   "source": [
    "# Example: Load and preprocess data\n",
    "# adata = sc.read_h5ad('your_data.h5ad')\n",
    "# sc.pp.normalize_total(adata)\n",
    "# sc.pp.log1p(adata)"
   ]
  },
  {
   "cell_type": "markdown",
   "metadata": {},
   "source": [
    "## Next Steps\n",
    "1. Download pretrained models from: https://github.com/bowang-lab/scGPT#pretrained-scGPT-checkpoints\n",
    "2. Follow tutorials at: https://github.com/bowang-lab/scGPT/tree/main/tutorials"
   ]
  }
 ],
 "metadata": {
  "kernelspec": {
   "display_name": "Python 3",
   "language": "python",
   "name": "python3"
  }
 },
 "nbformat": 4,
 "nbformat_minor": 4
}
EOF

# Print final instructions
print_status "Installation completed successfully!"
echo ""
echo "=========================================="
echo "scGPT has been installed in conda environment: ${ENV_NAME}"
echo ""
echo "To activate the environment, run:"
echo "  conda activate ${ENV_NAME}"
echo "  # or"
echo "  source ./activate_scgpt.sh"
echo ""
echo "To test the installation:"
echo "  python test_scgpt.py"
echo ""
echo "To start working with scGPT:"
echo "  jupyter notebook scgpt_example.ipynb"
echo ""
echo "Download pretrained models from:"
echo "  https://github.com/bowang-lab/scGPT#pretrained-scGPT-checkpoints"
echo ""
echo "For tutorials and documentation:"
echo "  https://github.com/bowang-lab/scGPT/tree/main/tutorials"
echo "  https://scgpt.readthedocs.io/"
echo "=========================================="

# Clean up
rm -f scgpt_environment.yml

# Deactivate environment
conda deactivate

print_status "Setup complete!"