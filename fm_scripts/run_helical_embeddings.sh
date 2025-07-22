

import logging
import warnings
from datasets import load_dataset

logging.getLogger().setLevel(logging.ERROR)

warnings.filterwarnings("ignore")

dataset = load_dataset("helical-ai/yolksac_human", split="train[:10%]", trust_remote_code=True, download_mode="reuse_cache_if_exists")
labels = dataset["LVL1"]



from helical.utils import get_anndata_from_hf_dataset

ann_data = get_anndata_from_hf_dataset(dataset)
ann_data




from helical import models
import pkgutil
for model in pkgutil.iter_modules(models.__path__):
    if model.name != "base_models" and model.name != "fine_tune":
        print("Model :",model.name)




from helical.models.geneformer import Geneformer, GeneformerConfig
import torch

device = "cuda" if torch.cuda.is_available() else "cpu"
model_config = GeneformerConfig(batch_size=10,device=device)
geneformer = Geneformer(configurer=model_config)



#from helical.models.uce.model import UCE,UCEConfig

#model_config=UCEConfig(batch_size=10)
#uce = UCE(configurer=model_config)

dataset = geneformer.process_data(ann_data[:100], gene_names="gene_name")



import umap
import seaborn as sns
import pandas as pd
import matplotlib.pyplot as plt



reducer = umap.UMAP(min_dist=0.1, n_components=2, n_epochs=None,n_neighbors=4)
mapper = reducer.fit(embeddings)

labels = labels[:100]

plot_df = pd.DataFrame(mapper.embedding_,columns=['px','py'])
plot_df['Cell Type'] = labels

plt.figure(figsize=(10,10))
ax = plt.axes()
sns.set_style('dark')
plt.style.use("dark_background")

sns.scatterplot(data = plot_df,x='px',y='py',hue='Cell Type',sizes=(50,200),ax=ax,palette="pastel")
plt.title('UMAP of Reference Data')
plt.show()

