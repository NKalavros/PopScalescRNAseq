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

# Install CUDA 11.6 toolkit for MindSpore GPU compatibility  
echo "Installing CUDA 11.6 toolkit for MindSpore GPU compatibility..."
mamba run --prefix "$FULL_ENV_PATH" conda install cudatoolkit=11.6 cudnn=8.4.1 cuda-nvcc==11.6.* -c nvidia -c conda-forge -y

# Install PyTorch with CUDA 11.6 (latest available version)
echo "Installing PyTorch with CUDA 11.6..."
mamba run --prefix "$FULL_ENV_PATH" pip install torch==1.13.1+cu116 torchvision==0.14.1+cu116 torchaudio==0.13.1+cu116 --index-url https://download.pytorch.org/whl/cu116

# Install core scientific packages (compatible versions for Python 3.9)
echo "Installing core scientific packages..."
mamba run --prefix "$FULL_ENV_PATH" pip install numpy==1.24.3 pandas==1.5.3 scipy==1.9.3 matplotlib==3.5.3

# Install compatible Pillow version for MindSpore
echo "Installing compatible Pillow version..."
mamba run --prefix "$FULL_ENV_PATH" pip install Pillow==9.5.0

# Install MindSpore GPU for CUDA 11.6
echo "Installing MindSpore GPU for CUDA 11.6..."
wget https://ms-release.obs.cn-north-4.myhuaweicloud.com/2.2.14/MindSpore/unified/x86_64/mindspore-2.2.14-cp39-cp39-linux_x86_64.whl
# Install single-cell analysis packages (exact versions from CellFM requirements)
echo "Installing single-cell packages..."
mamba run --prefix "$FULL_ENV_PATH" pip install scanpy==1.10 scib==1.1.5 anndata==0.8.0

# Install additional dependencies
echo "Installing additional dependencies..."
mamba run --prefix "$FULL_ENV_PATH" pip install transformers datasets accelerate huggingface_hub
mamba run --prefix "$FULL_ENV_PATH" pip install mindspore-2.2.14-cp39-cp39-linux_x86_64.whl

# Clone CellFM repository
echo "Cloning CellFM repository..."
cd "$FULL_ENV_PATH"
git clone https://github.com/biomed-AI/CellFM.git
cd CellFM

# Install CellFM package in development mode
echo "Installing CellFM package in development mode..."
#mamba run --prefix "$FULL_ENV_PATH" pip install -e .

# Create environment activation script with CUDA library paths
echo "Creating environment activation script..."
mkdir -p "$FULL_ENV_PATH/etc/conda/activate.d"
cat > "$FULL_ENV_PATH/etc/conda/activate.d/cuda_env.sh" << 'EOF'
#!/bin/bash
# Set CUDA library paths for MindSpore GPU
export LD_LIBRARY_PATH="$CONDA_PREFIX/lib:$LD_LIBRARY_PATH"
export CUDA_HOME="$CONDA_PREFIX"
export CUDA_ROOT="$CONDA_PREFIX"
export PATH="$CONDA_PREFIX/bin:$PATH"
echo "🔧 CUDA environment variables set for MindSpore GPU"
EOF

# Make activation script executable
chmod +x "$FULL_ENV_PATH/etc/conda/activate.d/cuda_env.sh"


mamba run --prefix "$FULL_ENV_PATH" pip install torch==1.13.1+cu116 torchvision==0.14.1+cu116 torchaudio==0.13.1+cu116 --index-url https://download.pytorch.org/whl/cu116

# Test the installation
echo "Testing installation..."
mamba run --prefix "$FULL_ENV_PATH" bash -c "
export LD_LIBRARY_PATH='$FULL_ENV_PATH/lib:\$LD_LIBRARY_PATH'
export CUDA_HOME='$FULL_ENV_PATH'
export CUDA_ROOT='$FULL_ENV_PATH'

python -c \"
import torch
print(f'✅ PyTorch: {torch.__version__}, CUDA available: {torch.cuda.is_available()}')
if torch.cuda.is_available():
    print(f'   CUDA version: {torch.version.cuda}')
    print(f'   GPU device: {torch.cuda.get_device_name(0)}')

import mindspore as ms
print(f'✅ MindSpore: {ms.__version__}')

# Test MindSpore GPU context
try:
    ms.set_context(mode=ms.GRAPH_MODE, device_target='GPU')
    print('✅ MindSpore GPU context set successfully')
    ms.set_context(device_target='CPU')  # Reset to CPU for safety
except Exception as e:
    print(f'⚠️  MindSpore GPU context failed: {e}')

import scanpy as sc
print('✅ scanpy imported successfully')

import transformers
print(f'✅ transformers: {transformers.__version__}')

import huggingface_hub
print(f'✅ huggingface_hub: {huggingface_hub.__version__}')

print('🎉 CellFM environment setup COMPLETE!')
\"
"

# Create results directory
echo "Creating results directory..."
mkdir -p "$FULL_ENV_PATH/../../results"

# Test embedding script
echo "Testing CellFM embedding script..."
echo "To test embeddings, run:"
echo "mamba activate $FULL_ENV_PATH"
echo "# CUDA environment variables are automatically set on activation"
echo "python fm_scripts/cellfm_embeddings.py --input pbmc3k_raw.h5ad --output results/cellfm_embeddings.h5ad"

echo "✅ CellFM environment ready at: $FULL_ENV_PATH"
echo "To activate: mamba activate $FULL_ENV_PATH"
echo "CellFM repository cloned to: $FULL_ENV_PATH/CellFM"
echo ""
echo "📝 Note: CellFM models available at:"
echo "   https://huggingface.co/ShangguanNingyuan/CellFM"