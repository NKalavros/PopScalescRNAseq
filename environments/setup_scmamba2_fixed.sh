ENV_PREFIX="/gpfs/scratch/nk4167/miniconda/envs"
ENV_NAME="scmamba2_env"
FULL_ENV_PATH="$ENV_PREFIX/$ENV_NAME"
rm -rf $FULL_ENV_PATH
conda create --prefix "$FULL_ENV_PATH" python=3.10 -y
conda activate "$FULL_ENV_PATH"
conda install -c conda-forge mamba -y
# Install cuda 12.8 from conda-forge
mamba install -c conda-forge -c nvidia cuda-toolkit=12.8 cuda-nvcc=12.8 -y
mamba install -c conda-forge gcc gxx -y
pip3 install torch torchvision torchaudio --index-url https://download.pytorch.org/whl/cu128
cd $FULL_ENV_PATH
git clone https://github.com/xtalpi-xic/SC-MAMBA2
cd SC-MAMBA2
# Install the requirements
# Relax the faiss-gpu requirement
#sed -i 's/faiss-gpu==1.7.2/faiss-gpu/g' requirements.txt
# Relax netron
sed -i 's/netron==7.8.0/netron/g' requirements.txt
# Relax tiledbsoma
sed -i 's/tiledbsoma==1.11.3/tiledbsoma/g' requirements.txt
sed -i 's/cellxgene-census==1.14.0/cellxgene-census/g' requirements.txt
# Relax transformers
sed -i 's/transformers==4.45.0.dev0/transformers/g' requirements.txt
# Relax numpy
sed -i 's/numpy==1.23.2/numpy/g' requirements.txt
# Relax fsspec
sed -i 's/fsspec==2024.5.0/fsspec/g' requirements.txt
sed -i 's/s3fs==2024.5.0/s3fs/g' requirements.txt
# Install from different packaging
#pip install faiss-gpu-cu12
mamba install -c conda-forge ninja -y #For fast flash-attn compilation
export MAX_JOBS=$SLURM_CPUS_PER_TASK
export CC=$CONDA_PREFIX/bin/gcc
export CXX=$CONDA_PREFIX/bin/g++
pip install -r requirements.txt
pip install flash_attn
pip install causal-conv1d
pip install mamba-ssm
# remove lines for requirements.txt for mamba-ssm and causal-conv1d
sed -i '/flash-attn/d' requirements.txt
sed -i '/mamba-ssm/d' requirements.txt
sed -i '/causal-conv1d/d' requirements.txt

# Install the rest of the requirements
pip install -r requirements.txt