#!/bin/bash
# Fixed Geneformer environment setup script for BigPurple HPC
# Compatible with CUDA 12 and HPC environment

set -e  # Exit on any error

ENV_PREFIX="/gpfs/scratch/nk4167/miniconda/envs"
ENV_NAME="geneformer_env"
FULL_ENV_PATH="$ENV_PREFIX/$ENV_NAME"

echo "🚀 Setting up Geneformer environment (fixed)..."

# Remove existing environment if it exists
if [ -d "$FULL_ENV_PATH" ]; then
    echo "Removing existing environment..."
    mamba env remove --prefix "$FULL_ENV_PATH" -y || true
fi

# Create minimal base environment
echo "Creating minimal base environment..."
mamba create --prefix "$FULL_ENV_PATH" python=3.10 pip git git-lfs -y

# Note: Activate environment manually with: mamba activate $FULL_ENV_PATH
echo "Environment created. Please activate manually and continue..."
echo "Commands to run:"
echo "  mamba activate $FULL_ENV_PATH"
echo "  [then run the pip install commands below]"
echo ""


# Add the libraries to the conda env setup
echo "Adding libstdc++ to LD_LIBRARY_PATH..."
mkdir -p "$FULL_ENV_PATH/etc/conda/activate.d"
echo "export LD_LIBRARY_PATH=\"$FULL_ENV_PATH/lib:\$LD_LIBRARY_PATH\"" >> "$FULL_ENV_PATH/etc/conda/activate.d/env_vars.sh"



# Install PyTorch with CUDA 12.1 (compatible with V100/A100)
echo "Installing PyTorch with CUDA 12.1..."
pip install torch==2.1.0+cu121 torchvision torchaudio --index-url https://download.pytorch.org/whl/cu121

# Install core scientific packages (compatible versions)
echo "Installing core scientific packages..."
pip install numpy==1.24.3 pandas==1.5.3 scipy==1.9.3 matplotlib==3.5.3

# Install machine learning packages
echo "Installing ML packages..."
pip install transformers datasets accelerate hf_transfer

# Install single-cell packages (compatible versions)
echo "Installing single-cell packages..."
pip install anndata==0.8.0 scanpy==1.8.2

pip install ipython mygene

# Clone and install Geneformer
echo "Cloning Geneformer repository..."
cd "$FULL_ENV_PATH"
git clone https://huggingface.co/ctheodoris/Geneformer
cd Geneformer
# Install gcc g++ and libstdc++-devel for compilation
echo "Installing system dependencies..."
mamba install -c conda-forge gcc_linux-64 gxx_linux-64 libstdcxx-ng -y
echo "Installing Geneformer package..."
pip install .
# Initialize git-lfs
echo "Initializing git-lfs..."
git lfs install
# Test the installation
echo "Testing installation..."
# Weird pandas mismatch
pip install pandas==2.2.2
pip install --upgrade scanpy
python -c "
import torch
print(f'✅ PyTorch: {torch.__version__}, CUDA: {torch.cuda.is_available()}')

import transformers
print(f'✅ transformers: {transformers.__version__}')

import scanpy as sc
print('✅ scanpy imported successfully')

try:
    import geneformer
    print('✅ geneformer imported successfully')
except ImportError as e:
    print(f'⚠️  geneformer import issue: {e}')
    print('This may be normal - will check specific modules in demo script')

print('🎉 Geneformer environment setup COMPLETE!')
"

# Get an optimized version with flash attention
# To try and speed things up, let's 
# Get pytorch version
python -c "
import torch
print(f'✅ PyTorch: {torch.__version__}, CUDA: {torch.cuda.is_available()}')
"
# Upgrade pytorch to cuda 12.8
pip3 install torch torchvision torchaudio --index-url https://download.pytorch.org/whl/cu128 --force-reinstall
# Install nvcc and gcc and g++
mamba install -c conda-forge -c nvidia cuda-runtime==12.8* cuda-nvcc==12.8*  mamba install cuda-toolkit==12.8 -y-y
# Install gcc and g++
mamba install -c conda-forge gcc gxx libstdcxx-ng -y
# Install cuda runtime
# Export gcc and g++ variables
export CC=$CONDA_PREFIX/bin/gcc
export CXX=$CONDA_PREFIX/bin/g++
# Export the include path
export C_INCLUDE_PATH=$CONDA_PREFIX/include
export CPLUS_INCLUDE_PATH=$CONDA_PREFIX/include
# Install the rest of the requirements
cd $FULL_ENV_PATH
cd Geneformer
module load cuda
mamba install -c conda-forge ninja -y
git clone https://github.com/Dao-AILab/flash-attention
cd flash-attention
cd hopper
export MAX_JOBS=$SLURM_CPUS_PER_TASK
python setup.py install

echo "✅ Geneformer environment ready at: $FULL_ENV_PATH"
echo "To activate: mamba activate $FULL_ENV_PATH"
echo "Model repository cloned to: $FULL_ENV_PATH/Geneformer"
