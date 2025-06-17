mamba create --prefix /gpfs/scratch/nk4167/miniconda/envs/CellFM python=3.9 -y
conda activate /gpfs/scratch/nk4167/miniconda/envs/CellFM
pip install mindspore==2.2.10 scanpy==1.10 scib==1.1.5
# Older pytorch here since it's 3.9
pip3 install torch torchvision torchaudio --index-url https://download.pytorch.org/whl/cu121
