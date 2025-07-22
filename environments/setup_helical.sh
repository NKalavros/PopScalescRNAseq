ENV_PREFIX="/gpfs/scratch/nk4167/miniconda/envs"
ENV_NAME="helical-package"
FULL_ENV_PATH="$ENV_PREFIX/$ENV_NAME"
mamba create --prefix "$FULL_ENV_PATH" python=3.11.8 -y
conda activate "$FULL_ENV_PATH"
mamba install -c conda-forge mamba -y
# Install cuda 12.8 from conda-forge
mamba install -c conda-forge -c nvidia cuda-toolkit=12.6 -y
# Install torch and cuda
pip install helical
pip install helical[mamba-ssm] #Second installatoin to get causal conv
mamba install -c conda-forge -c nvidia cuda-toolkit=12.6 gcc gxx -y
pip install --no-build-isolation "git+https://github.com/state-spaces/mamba@v2.2.4"
pip install helical[mamba-ssm] #Second installatoin to get causal conv
pip install ipython
# Ensure pytorch is available
python -c "
import torch
print(f'✅ PyTorch: {torch.__version__}, CUDA: {torch.cuda.is_available()}')
"