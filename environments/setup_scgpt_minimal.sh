#!/bin/bash
# Minimal scGPT environment setup script for BigPurple HPC
# Creates base environment then installs packages in correct order

set -e  # Exit on any error

ENV_PREFIX="/gpfs/scratch/nk4167/miniconda/envs"
ENV_NAME="scgpt_env"
FULL_ENV_PATH="$ENV_PREFIX/$ENV_NAME"

echo "🚀 Setting up scGPT environment (minimal approach)..."

# Remove existing environment if it exists
if [ -d "$FULL_ENV_PATH" ]; then
    echo "Removing existing environment..."
    mamba env remove --prefix "$FULL_ENV_PATH" -y || true
fi

# Create minimal base environment
echo "Creating minimal base environment..."
mamba create --prefix "$FULL_ENV_PATH" python=3.9 pip git -y

# Activate environment
echo "Activating environment..."
source activate "$FULL_ENV_PATH"

# Install PyTorch first (most important)
echo "Installing PyTorch with CUDA 11.7..."
pip install torch==2.0.1+cu117 torchvision torchaudio --index-url https://download.pytorch.org/whl/cu117

# Install core scientific packages
echo "Installing core scientific packages..."
pip install numpy==1.24.3 pandas==1.5.3 scipy==1.9.3

# Install compatibility packages
echo "Installing compatibility packages..."
pip install numba==0.56.4 llvmlite==0.39.1 matplotlib==3.5.3 scikit-learn==1.1.3

# Install single-cell packages
echo "Installing single-cell packages..."
pip install anndata==0.8.0 scanpy==1.8.2

# Install other required packages
echo "Installing other packages..."
pip install torchtext==0.15.2 mudata==0.2.3 scvi-tools==0.20.3 wandb orbax==0.1.7

# Install scGPT last
echo "Installing scGPT..."
pip install scgpt

# Test the installation
echo "Testing installation..."
python -c "
import torch
print(f'✅ PyTorch: {torch.__version__}, CUDA: {torch.cuda.is_available()}')

import scanpy as sc
print('✅ Scanpy imported successfully')

import scgpt
print('✅ scGPT imported successfully')

print('🎉 scGPT environment setup COMPLETE!')
"

echo "✅ scGPT environment ready at: $FULL_ENV_PATH"
echo "To activate: mamba activate $FULL_ENV_PATH"