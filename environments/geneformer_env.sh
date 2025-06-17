mamba create --prefix /gpfs/scratch/nk4167/miniconda/envs/geneformer_env python=3.10 -y
conda activate /gpfs/scratch/nk4167/miniconda/envs/geneformer_env
# Install PyTorch with CUDA 12.1 support
pip3 install torch torchvision torchaudio --index-url https://download.pytorch.org/whl/cu121
# Install git through mamba
mamba install git -y
# Install git-lfs for large file support
mamba install git-lfs -y
# Initialize git-lfs
git lfs install
# Clone the Geneformer repository and install it
git clone https://huggingface.co/ctheodoris/Geneformer
cd Geneformer
pip install .