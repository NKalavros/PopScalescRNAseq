mamba create --prefix /gpfs/scratch/nk4167/miniconda/envs/scgpt_env python=3.10 -y
conda activate /gpfs/scratch/nk4167/miniconda/envs/scgpt_env
# Install PyTorch with CUDA 12.8 support
# Note: Ensure that the CUDA version matches your system's installed version with module load cuda/12.8
pip3 install torch torchvision torchaudio --index-url https://download.pytorch.org/whl/cu128
# Install nvcc
mamba install -c conda-forge cudatoolkit-dev -y
# Install flash-attn with a version constraint
# Note: flash-attn<1.0.5 is required for compatibility with scgpt
pip install scgpt "flash-attn<1.0.5"
# Install gdown
pip install gdown
# Download the model from https://drive.google.com/drive/folders/1oWh_-ZRdhtoGQ2Fw24HP41FgLoomVo-y
gdown --folder https://drive.google.com/drive/folders/1oWh_-ZRdhtoGQ2Fw24HP41FgLoomVo-yg -O /gpfs/scratch/nk4167/miniconda/envs/scgpt_env/
# The model is now available in /gpfs/scratch/nk4167/miniconda/envs/scgpt_env/
# You can now use scgpt in this environment