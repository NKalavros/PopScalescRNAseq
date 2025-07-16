#!/usr/bin/env python
"""
Create realistic test dataset with actual gene symbols for foundation model testing
"""

import numpy as np
import pandas as pd
from scipy.sparse import csr_matrix
import anndata as ad

def create_realistic_test_data():
    """Create a test .h5ad file with realistic gene symbols"""
    print("Creating realistic test dataset...")
    
    # Common human gene symbols (a subset that's likely to overlap with foundation models)
    common_genes = [
        'ACTB', 'GAPDH', 'RPL23', 'RPS18', 'EEF1A1', 'TPT1', 'FTL', 'FTH1',
        'B2M', 'CALM1', 'CALM2', 'CALM3', 'CD74', 'HLA-A', 'HLA-B', 'HLA-C',
        'MALAT1', 'NEAT1', 'XIST', 'MT-CO1', 'MT-CO2', 'MT-CO3', 'MT-ND1',
        'MT-ND2', 'MT-ND4', 'MT-ATP6', 'MT-ATP8', 'MT-CYB', 'MT-ND4L',
        'CD3D', 'CD3E', 'CD3G', 'CD4', 'CD8A', 'CD8B', 'CD19', 'CD20',
        'IL2RA', 'FOXP3', 'GZMB', 'PRF1', 'IFNG', 'TNF', 'IL1B', 'IL6',
        'CCL5', 'CXCL10', 'CSF1R', 'CD68', 'CD163', 'ARG1', 'NOS2',
        'ALB', 'AFP', 'CRP', 'APOE', 'LDLR', 'PCSK9', 'SREBF1', 'SREBF2',
        'INS', 'IGF1', 'INSR', 'IRS1', 'PIK3CA', 'AKT1', 'MTOR', 'TSC1',
        'TP53', 'MDM2', 'CDKN1A', 'RB1', 'E2F1', 'MYC', 'BCL2', 'BAX',
        'BRCA1', 'BRCA2', 'ATM', 'CHEK1', 'CHEK2', 'DNA-PK', 'PARP1',
        'EGFR', 'ERBB2', 'PDGFRA', 'KIT', 'FLT3', 'VEGFA', 'VEGFR1',
        'TGFB1', 'SMAD2', 'SMAD3', 'SMAD4', 'WNT1', 'CTNNB1', 'APC',
        'NOTCH1', 'HES1', 'SOX2', 'OCT4', 'NANOG', 'KLF4', 'MYC',
        'VIM', 'CDH1', 'CDH2', 'SNAI1', 'TWIST1', 'ZEB1', 'ZEB2',
        # Add more common genes to reach a good number
        'ACTA2', 'COL1A1', 'COL3A1', 'FN1', 'PECAM1', 'VWF', 'KDR',
        'PTPRC', 'CD45', 'LYZ', 'S100A8', 'S100A9', 'ELANE', 'MPO',
        'FCGR3A', 'NCAM1', 'KLRB1', 'GNLY', 'NKG7', 'GZMA', 'GZMH',
        'CCR7', 'SELL', 'LEF1', 'TCF7', 'IL7R', 'CD27', 'CD28',
        'ICOS', 'CTLA4', 'PDCD1', 'LAG3', 'HAVCR2', 'TIGIT', 'TOX'
    ]
    
    # Extend to get more genes (repeat with slight variations)
    extended_genes = common_genes.copy()
    for i in range(len(common_genes), 2000):
        base_gene = common_genes[i % len(common_genes)]
        extended_genes.append(f"{base_gene}P{i//len(common_genes)}")
    
    gene_names = extended_genes[:2000]  # Use first 2000
    
    # Create synthetic data
    n_obs = 1000  # cells
    n_vars = len(gene_names)  # genes
    
    # Generate more realistic expression data (log-normal distribution)
    np.random.seed(42)
    
    # Create expression matrix with realistic patterns
    # Some genes highly expressed (housekeeping), others lowly expressed
    base_expression = np.random.lognormal(mean=1.0, sigma=1.5, size=(n_obs, n_vars))
    
    # Add some structure - make some genes consistently high/low
    housekeeping_indices = [i for i, gene in enumerate(gene_names) 
                           if gene in ['ACTB', 'GAPDH', 'RPL23', 'RPS18', 'EEF1A1']]
    for idx in housekeeping_indices:
        base_expression[:, idx] *= 5  # Make housekeeping genes higher
    
    # Add some noise and convert to counts-like data
    expression_data = (base_expression * 10).astype(np.int32)
    X = csr_matrix(expression_data.astype(np.float32))
    
    # Create cell barcodes
    cell_barcodes = [f"Cell_{i:04d}" for i in range(n_obs)]
    
    # Create AnnData object
    adata = ad.AnnData(
        X=X,
        obs=pd.DataFrame(index=cell_barcodes),
        var=pd.DataFrame(index=gene_names)
    )
    
    # Add gene names (scimilarity expects gene symbols in var.index)
    adata.var["gene_name"] = gene_names
    adata.var_names = gene_names  # Ensure var.index has gene symbols
    
    print(f"Realistic test data shape: {adata.shape}")
    print(f"Sample gene names: {adata.var_names[:10].tolist()}")
    print(f"Expression range: {adata.X.min():.2f} - {adata.X.max():.2f}")
    
    # Save test data
    output_path = "test_data_realistic.h5ad"
    adata.write_h5ad(output_path)
    print(f"✅ Realistic test data saved to {output_path}")
    
    return output_path

if __name__ == "__main__":
    create_realistic_test_data()