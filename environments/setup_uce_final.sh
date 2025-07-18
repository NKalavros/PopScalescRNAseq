ENV_PREFIX="/gpfs/scratch/nk4167/miniconda/envs"
ENV_NAME="uce_env"
FULL_ENV_PATH="$ENV_PREFIX/$ENV_NAME"

echo "🚀 Setting up scGPT environment..."

# Remove existing environment if it exists
if [ -d "$FULL_ENV_PATH" ]; then
    echo "Removing existing environment..."
    mamba env remove --prefix "$FULL_ENV_PATH" -y || true
fi

# Create environment from YAML with correct prefix
mamba create -c conda-forge --prefix "$FULL_ENV_PATH" python=3.10 -y
mamba install --prefix "$FULL_ENV_PATH" -c conda-forge pip -y

# Activate environment
echo "Activating environment..."
eval "$(conda shell.bash hook)"
conda activate "$FULL_ENV_PATH"

# Clone and enter UCE directory
git clone https://github.com/snap-stanford/UCE/
cd UCE
# Install the CUDA packages (CUDA 11)
echo "Installing PyTorch with CUDA 11.8..."
pip install torch==2.1.1+cu118 torchvision==0.16.1+cu118 torchaudio==2.1.1+cu118 --index-url https://download.pytorch.org/whl/cu118
pip install numpy==1.26.4
# Check that cuda is available
python -c "
import torch
print(f'✅ PyTorch: {torch.__version__}, CUDA: {torch.cuda.is_available()}')
"
# Install from the requirements file
echo "Installing requirements..."
pip install -r requirements.txt
# Run the example (downloads requisite files)
python eval_single_anndata.py
# Download larger model from https://figshare.com/ndownloader/articles/24320806/versions/5 with a user agent for firefox
echo "Downloading larger model..."
# Use curl to download the file with a user agent
# This is necessary because figshare requires a user agent to download files
curl -L -o 33l_8ep_1024t_1280.torch https://figshare.com/ndownloader/articles/24320806/versions/5 -H "User-Agent: Mozilla/5.0 (Windows NT 10.0; Win64; x64; rv:109.0) Gecko/20100101 Firefox/115.0"