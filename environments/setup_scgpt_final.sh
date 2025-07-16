#!/bin/bash
# Complete scGPT environment setup script for BigPurple HPC
# This script creates a working scGPT environment with all compatibility fixes

set -e  # Exit on any error

ENV_PREFIX="/gpfs/scratch/nk4167/miniconda/envs"
ENV_NAME="scgpt_env"
FULL_ENV_PATH="$ENV_PREFIX/$ENV_NAME"

echo "🚀 Setting up scGPT environment..."

# Remove existing environment if it exists
if [ -d "$FULL_ENV_PATH" ]; then
    echo "Removing existing environment..."
    mamba env remove --prefix "$FULL_ENV_PATH" -y || true
fi

# Create environment from YAML with correct prefix
echo "Creating environment from YAML..."
mamba env create -f environments/scgpt_env.yml --prefix "$FULL_ENV_PATH"

# Activate environment
echo "Activating environment..."
source activate "$FULL_ENV_PATH"

# Install PyTorch with specific CUDA version (must be done after environment creation)
echo "Installing PyTorch with CUDA 11.7..."
pip install torch==2.0.1+cu117 torchvision torchaudio --index-url https://download.pytorch.org/whl/cu117

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