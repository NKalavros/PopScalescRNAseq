#!/bin/bash
# CellFM Environment Setup with Real Model Support

ENV_PREFIX="/gpfs/scratch/nk4167/miniconda/envs"
ENV_NAME="CellFM"
FULL_ENV_PATH="$ENV_PREFIX/$ENV_NAME"

echo "🚀 Setting up CellFM environment with real model support..."

# Create environment
mamba create --prefix "$FULL_ENV_PATH" python=3.9 pip git -y

# Install core packages
mamba run --prefix "$FULL_ENV_PATH" pip install mindspore==2.2.10 scanpy==1.10 scib==1.1.5 anndata==0.8.0

# Install model download dependencies  
mamba run --prefix "$FULL_ENV_PATH" pip install huggingface_hub transformers

# Install PyTorch (optional, for compatibility)
mamba run --prefix "$FULL_ENV_PATH" pip install torch torchvision torchaudio --index-url https://download.pytorch.org/whl/cu121

# Clone CellFM repository
cd "$FULL_ENV_PATH"
git clone https://github.com/biomed-AI/CellFM.git

echo "✅ CellFM environment ready with real model support!"
echo "To activate: mamba activate $FULL_ENV_PATH"
