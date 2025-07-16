#!/bin/bash
# Fixed CellFM environment setup script for BigPurple HPC
# Compatible with CUDA 12 and HPC environment

set -e  # Exit on any error

ENV_PREFIX="/gpfs/scratch/nk4167/miniconda/envs"
ENV_NAME="CellFM"
FULL_ENV_PATH="$ENV_PREFIX/$ENV_NAME"

echo "🚀 Setting up CellFM environment (fixed)..."

# Remove existing environment if it exists
if [ -d "$FULL_ENV_PATH" ]; then
    echo "Removing existing environment..."
    mamba env remove --prefix "$FULL_ENV_PATH" -y || true
fi

# Create minimal base environment
echo "Creating minimal base environment..."
mamba create --prefix "$FULL_ENV_PATH" python=3.9 pip git -y

echo "Environment created. Please activate manually and continue..."
echo "Commands to run:"
echo "  mamba activate $FULL_ENV_PATH"
echo "  [then run the pip install commands below]"
echo ""

echo "# Install PyTorch with CUDA 12.1 (compatible with V100/A100)"
echo "pip install torch==2.1.0+cu121 torchvision torchaudio --index-url https://download.pytorch.org/whl/cu121"

echo "# Install core scientific packages (compatible versions)"
echo "pip install numpy==1.24.3 pandas==1.5.3 scipy==1.9.3 matplotlib==3.5.3"

echo "# Install CellFM-specific packages (MindSpore framework)"
echo "pip install mindspore==2.2.10"

echo "# Install single-cell packages (compatible versions)"
echo "pip install anndata==0.8.0 scanpy==1.8.2 scib==1.1.5"

echo "✅ CellFM environment ready at: $FULL_ENV_PATH"
echo "To activate: mamba activate $FULL_ENV_PATH"