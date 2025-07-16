#!/bin/bash
# CellFM Environment Setup with Real Model Support

ENV_PREFIX="/gpfs/scratch/nk4167/miniconda/envs"
ENV_NAME="CellFM"
FULL_ENV_PATH="$ENV_PREFIX/$ENV_NAME"

echo "🚀 Setting up CellFM environment with real model support..."

# Create environment with CUDA 11.6 for MindSpore GPU
mamba create --prefix "$FULL_ENV_PATH" python=3.9 pip git -y

# Install CUDA 11.6 toolkit for MindSpore compatibility
echo "Installing CUDA 11.6 for MindSpore GPU compatibility..."
mamba run --prefix "$FULL_ENV_PATH" conda install cudatoolkit=11.6 cudnn=8.4.1 -c conda-forge -y

# Install MindSpore GPU for CUDA 11.6
echo "Installing MindSpore GPU for CUDA 11.6..."
mamba run --prefix "$FULL_ENV_PATH" pip install mindspore-gpu==1.10.0 -f https://ms-release.obs.cn-north-4.myhuaweicloud.com/1.10.0/MindSpore/gpu/x86_64/cuda-11.6/

# Install core packages
mamba run --prefix "$FULL_ENV_PATH" pip install scanpy==1.10 scib==1.1.5 anndata==0.8.0

# Install model download dependencies  
mamba run --prefix "$FULL_ENV_PATH" pip install huggingface_hub transformers

# Install PyTorch with CUDA 11.6 (latest available version)
mamba run --prefix "$FULL_ENV_PATH" pip install torch==1.13.1+cu116 torchvision==0.14.1+cu116 torchaudio==0.13.1+cu116 --index-url https://download.pytorch.org/whl/cu116

# Clone CellFM repository
cd "$FULL_ENV_PATH"
git clone https://github.com/biomed-AI/CellFM.git

echo "✅ CellFM environment ready with real model support!"
echo "To activate: mamba activate $FULL_ENV_PATH"
