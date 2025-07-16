#!/bin/bash
# Fixed scFoundation environment setup script for BigPurple HPC
# Compatible with CUDA 12 and HPC environment

set -e  # Exit on any error

ENV_PREFIX="/gpfs/scratch/nk4167/miniconda/envs"
ENV_NAME="scFoundation_env"
FULL_ENV_PATH="$ENV_PREFIX/$ENV_NAME"

echo "🚀 Setting up scFoundation environment (fixed)..."

# Remove existing environment if it exists
if [ -d "$FULL_ENV_PATH" ]; then
    echo "Removing existing environment..."
    mamba env remove --prefix "$FULL_ENV_PATH" -y || true
fi

# Create minimal base environment
echo "Creating minimal base environment..."
mamba create --prefix "$FULL_ENV_PATH" python=3.10 pip git -y

# Activate environment
echo "Activating environment..."
source activate "$FULL_ENV_PATH"

# Install PyTorch with CUDA 12.1 (compatible with V100/A100)
echo "Installing PyTorch with CUDA 12.1..."
pip install torch==2.1.0+cu121 torchvision torchaudio --index-url https://download.pytorch.org/whl/cu121

# Install core scientific packages (compatible versions)
echo "Installing core scientific packages..."
pip install numpy==1.24.3 pandas==1.5.3 scipy==1.9.3

# Install scFoundation-specific packages
echo "Installing scFoundation dependencies..."
pip install einops local_attention

# Install single-cell packages (compatible versions)
echo "Installing single-cell packages..."
pip install anndata==0.8.0 scanpy==1.8.2

# Create models directory
echo "Creating models directory..."
mkdir -p "$FULL_ENV_PATH/models"

echo "⚠️  Model download requires manual authentication."
echo "Please download the following files manually and place them in:"
echo "  $FULL_ENV_PATH/models/"
echo ""
echo "Model files needed:"
echo "  - models.ckpt"
echo "  - models1.ckpt"
echo ""
echo "Download from: https://hopebio2020.sharepoint.com/sites/PublicSharedfiles/Shared%20Documents/Forms/AllItems.aspx"
echo ""

# Test the installation
echo "Testing installation..."
python -c "
import torch
print(f'✅ PyTorch: {torch.__version__}, CUDA: {torch.cuda.is_available()}')

import einops
print('✅ einops imported successfully')

import local_attention
print('✅ local_attention imported successfully')

import scanpy as sc
print('✅ scanpy imported successfully')

print('🎉 scFoundation environment setup COMPLETE!')
print('⚠️  Remember to download model files manually!')
"

echo "✅ scFoundation environment ready at: $FULL_ENV_PATH"
echo "To activate: mamba activate $FULL_ENV_PATH"
echo ""
echo "🔥 IMPORTANT: You need to manually download the model files!"
echo "   Place models.ckpt and models1.ckpt in: $FULL_ENV_PATH/models/"