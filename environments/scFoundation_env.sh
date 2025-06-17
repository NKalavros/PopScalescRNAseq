mamba create --prefix /gpfs/scratch/nk4167/miniconda/envs/scFoundation_env python=3.10 -y
conda activate /gpfs/scratch/nk4167/miniconda/envs/scFoundation_env
pip install numpy
pip install pandas
pip install scipy
pip install torch==2.5.0+cu121 --index-url https://download.pytorch.org/whl/cu121
pip install einops
pip install scanpy
pip install local_attention
#Download model from this URL: https://hopebio2020.sharepoint.com/sites/PublicSharedfiles/Shared%20Documents/Forms/AllItems.aspx?id=%2Fsites%2FPublicSharedfiles%2FShared%20Documents%2FPublic%20Shared%20files&p=true&ga=1
mkdir -p /gpfs/scratch/nk4167/miniconda/envs/scFoundation_env/models
wget -O /gpfs/scratch/nk4167/miniconda/envs/scFoundation_env/models/scFoundation_v1.0.tar.gz "https://hopebio2020.sharepoint.com/sites/PublicSharedfiles/Shared%20Documents/Public%20Shared%20files/scFoundation_v1.0.tar.gz?download=1"
tar -xzf /gpfs/scratch/nk4167/miniconda/envs/scFoundation_env/models/scFoundation_v1.0.tar.gz -C /gpfs/scratch/nk4167/miniconda/envs/scFoundation_env/models/