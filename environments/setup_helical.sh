ENV_PREFIX="/gpfs/scratch/nk4167/miniconda/envs"
ENV_NAME="helical-package"
FULL_ENV_PATH="$ENV_PREFIX/$ENV_NAME"
conda create --prefix "$FULL_ENV_PATH" python=3.11.8 -y
conda activate "$FULL_ENV_PATH"
conda install -c conda-forge mamba -y
# Install torch and cuda
pip install torch==2.1.1+cu118 torchvision==0.16.1+cu118 torchaudio==2.1.1+cu118 --index-url https://download.pytorch.org/whl/cu118
pip install helical
# Ensure pytorch is available
python -c "
import torch
print(f'✅ PyTorch: {torch.__version__}, CUDA: {torch.cuda.is_available()}')
"