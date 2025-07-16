mamba create --prefix /gpfs/scratch/nk4167/miniconda/envs/scimilarity_env python=3.10 -y
conda activate /gpfs/scratch/nk4167/miniconda/envs/scimilarity_env
pip install scimilarity
curl -L 'https://zenodo.org/records/10685499/files/model_v1.1.tar.gz?download=1' -o /gpfs/scratch/nk4167/miniconda/envs/scimilarity_env/model_v1.1.tar.gz
tar -xzf /gpfs/scratch/nk4167/miniconda/envs/scimilarity_env/model_v1.1.tar.gz -C /gpfs/scratch/nk4167/miniconda/envs/scimilarity_env/
rm /gpfs/scratch/nk4167/miniconda/envs/scimilarity_env/model_v1.1.tar.gz
# The model is now available at /gpfs/scratch/nk4167/miniconda/envs/scimilarity_env/model_v1.1
