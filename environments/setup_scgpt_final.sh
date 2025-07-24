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
pip install torchtext==0.15.2 scgpt scvi-tools==0.20.3 orbax==0.1.7 torch==2.1.1+cu117 torchvision==0.16.1+cu117 torchaudio==2.1.1+cu117 --index-url https://download.pytorch.org/whl/cu117

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

cd /gpfs/scratch/nk4167/miniconda/envs/scgpt_env
mkdir -p models
# Gdrive link: https://drive.google.com/drive/folders/1oWh_-ZRdhtoGQ2Fw24HP41FgLoomVo-y?usp=sharing
echo "Downloading scGPT model..."
gdown --id 1oWh_-ZRdhtoGQ2Fw24HP41FgLoomVo-y -O models/scgpt_model.pth
# Fix up libstdc++ issue
echo "Fixing libstdc++ issue..."
mamba install -c conda-forge libstdcxx-ng -y

# Add the library to the LD_LIBRARY_PATH when running this env
echo "Adding libstdc++ to LD_LIBRARY_PATH..."
mkdir -p "$FULL_ENV_PATH/etc/conda/activate.d"
echo "export LD_LIBRARY_PATH=\"$FULL_ENV_PATH/lib:\$LD_LIBRARY_PATH\"" >> "$FULL_ENV_PATH/etc/conda/activate.d/env_vars.sh"
# Install nvcc
echo "Installing nvcc..."
mamba install -c conda-forge cudatoolkit-dev -y
# Installing gcc
echo "Installing gcc..."
 mamba install -c conda-forge gcc_linux-64==11* gxx_linux-64==11* -y# Installing flash-attn
pip install "flash-attn<1.0.5"
echo "✅ scGPT environment ready at: $FULL_ENV_PATH"
echo "To activate: mamba activate $FULL_ENV_PATH"