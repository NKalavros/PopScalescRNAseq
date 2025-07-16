# macOS-compatible scGPT environment setup
mamba create --prefix ~/miniconda/envs/scgpt_env python=3.10 -y
conda activate ~/miniconda/envs/scgpt_env

# Install CPU-only PyTorch for macOS (no CUDA support)
mamba install -c pytorch -c conda-forge pytorch=2.3.0 torchvision torchaudio cpuonly -y

# Install flash-attn with a version constraint (CPU version)
pip install "flash-attn<=2" --no-build-isolation

# Install scgpt
pip install scgpt

# Install additional dependencies
pip install gdown ipython

# Install scvi-tools
pip install --upgrade scvi-tools

# Install wandb
pip install wandb

# Remove torchtext if it causes compatibility issues
pip uninstall torchtext -y

# Create local model directory
mkdir -p ~/miniconda/envs/scgpt_env/models

# Download the model from Google Drive
gdown --folder https://drive.google.com/drive/folders/1oWh_-ZRdhtoGQ2Fw24HP41FgLoomVo-y -O ~/miniconda/envs/scgpt_env/models/

# Download additional model files
gdown --id 1x1SfmFdI-zcocmqWAd7ZTC9CTEAVfKZq -O ~/miniconda/envs/scgpt_env/models/scgpt_model.ckpt
gdown --id 1jfT_T5n8WNbO9QZcLWObLdRG8lYFKH-Q -O ~/miniconda/envs/scgpt_env/models/vocab.json
gdown --id 15TEZmd2cZCrHwgfE424fgQkGUZCXiYrR -O ~/miniconda/envs/scgpt_env/models/args.json

echo "scGPT environment setup complete for macOS!"
echo "Activate with: conda activate ~/miniconda/envs/scgpt_env"
