#!/bin/bash
# Complete CellFM environment setup script for BigPurple HPC
# Based on official CellFM repository requirements

set -e  # Exit on any error

ENV_PREFIX="/gpfs/scratch/nk4167/miniconda/envs"
ENV_NAME="CellFM"
FULL_ENV_PATH="$ENV_PREFIX/$ENV_NAME"

echo "🚀 Setting up CellFM environment (complete)..."

# Remove existing environment if it exists
if [ -d "$FULL_ENV_PATH" ]; then
    echo "Removing existing environment..."
    mamba env remove --prefix "$FULL_ENV_PATH" -y || true
fi

# Create base environment (Python 3.9 as required by CellFM)
echo "Creating CellFM environment with Python 3.9..."
mamba create --prefix "$FULL_ENV_PATH" python=3.9 pip git -y

# Use conda run to execute commands in the environment (this works without shell activation)
echo "Installing packages in CellFM environment..."

# Install PyTorch with CUDA 12.1 (optional dependency but recommended)
echo "Installing PyTorch with CUDA 12.1..."
mamba run --prefix "$FULL_ENV_PATH" pip install torch==2.1.0+cu121 torchvision torchaudio --index-url https://download.pytorch.org/whl/cu121

# Install core scientific packages (compatible versions for Python 3.9)
echo "Installing core scientific packages..."
mamba run --prefix "$FULL_ENV_PATH" pip install numpy==1.24.3 pandas==1.5.3 scipy==1.9.3 matplotlib==3.5.3

# Install MindSpore (CellFM's main requirement) - try GPU version first
echo "Installing MindSpore..."
echo "Trying GPU version first..."
mamba run --prefix "$FULL_ENV_PATH" pip install mindspore-gpu==2.2.14 || {
    echo "GPU version failed, installing CPU version..."
    mamba run --prefix "$FULL_ENV_PATH" pip install mindspore==2.2.14
}

# Install single-cell analysis packages (exact versions from CellFM requirements)
echo "Installing single-cell packages..."
mamba run --prefix "$FULL_ENV_PATH" pip install scanpy==1.10 scib==1.1.5 anndata==0.8.0

# Install additional dependencies
echo "Installing additional dependencies..."
mamba run --prefix "$FULL_ENV_PATH" pip install transformers datasets accelerate huggingface_hub

# Clone CellFM repository
echo "Cloning CellFM repository..."
cd "$FULL_ENV_PATH"
git clone https://github.com/biomed-AI/CellFM.git
cd CellFM

# Install CellFM package in development mode
echo "Installing CellFM package in development mode..."
mamba run --prefix "$FULL_ENV_PATH" pip install -e .

# Test the installation
echo "Testing installation..."
mamba run --prefix "$FULL_ENV_PATH" python -c "
import torch
print(f'✅ PyTorch: {torch.__version__}, CUDA: {torch.cuda.is_available()}')

import mindspore as ms
print(f'✅ MindSpore: {ms.__version__}')

import scanpy as sc
print('✅ scanpy imported successfully')

import transformers
print(f'✅ transformers: {transformers.__version__}')

import huggingface_hub
print(f'✅ huggingface_hub: {huggingface_hub.__version__}')

try:
    # Test CellFM module imports
    import sys
    sys.path.append('.')
    from model import CellFM
    from config import Config
    print('✅ CellFM modules imported successfully')
except Exception as e:
    print(f'⚠️  CellFM modules: {e}')

print('🎉 CellFM environment setup COMPLETE!')
"

# Create results directory
echo "Creating results directory..."
mkdir -p "$FULL_ENV_PATH/../../results"

# Test embedding script
echo "Testing CellFM embedding script..."
echo "To test embeddings, run:"
echo "mamba activate $FULL_ENV_PATH"
echo "python fm_scripts/cellfm_embeddings.py --input pbmc3k_raw.h5ad --output results/cellfm_embeddings.h5ad"

echo "✅ CellFM environment ready at: $FULL_ENV_PATH"
echo "To activate: mamba activate $FULL_ENV_PATH"
echo "CellFM repository cloned to: $FULL_ENV_PATH/CellFM"
echo ""
echo "📝 Note: CellFM models available at:"
echo "   https://huggingface.co/ShangguanNingyuan/CellFM"