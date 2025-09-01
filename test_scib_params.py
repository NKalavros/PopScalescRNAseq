#!/usr/bin/env python3
"""
Test script to validate scIB parameter fixes.
"""
import sys
import os
import numpy as np
import pandas as pd
import scanpy as sc

# Add the viz directory to the path
viz_dir = os.path.join(os.path.dirname(__file__), 'viz')
sys.path.insert(0, viz_dir)

def test_scib_parameters():
    """Test that we can call scIB metrics with correct parameter names."""
    try:
        import scib
        from scib.metrics import metrics as scib_metrics
        print("✅ scIB imported successfully")
    except ImportError as e:
        print(f"❌ scIB not available: {e}")
        return False
        
    # Create some minimal test data
    np.random.seed(42)
    n_cells = 100
    n_genes = 50
    
    # Original data
    X = np.random.negative_binomial(5, 0.3, (n_cells, n_genes))
    obs = pd.DataFrame({
        'batch': np.random.choice(['batch1', 'batch2'], n_cells),
        'cell_type': np.random.choice(['TypeA', 'TypeB', 'TypeC'], n_cells)
    })
    var = pd.DataFrame(index=[f'gene_{i}' for i in range(n_genes)])
    adata = sc.AnnData(X=X, obs=obs, var=var)
    
    # "Integrated" data (just some random embedding)
    emb = np.random.normal(0, 1, (n_cells, 10))
    adata_int = sc.AnnData(X=emb, obs=obs.copy())
    
    print("✅ Test data created")
    
    # Test the corrected parameter names
    try:
        results = scib_metrics(
            adata,
            adata_int,
            batch_key='batch',
            label_key='cell_type',
            ari_=True,  # Corrected parameter name
            nmi_=True,  # Corrected parameter name
            silhouette_=True,  # Corrected parameter name
            graph_conn_=False,  # Disable slow ones for testing
            pcr_=False,
            ilisi_=False,
            clisi_=False,
            kBET_=False,
            isolated_labels_f1_=False,
            isolated_labels_asw_=False,
            hvg_score_=False,
            cell_cycle_=False,
            trajectory_=False,
        )
        print("✅ scIB metrics call succeeded with corrected parameter names")
        print(f"✅ Results shape: {results.shape}")
        print(f"✅ Columns: {list(results.columns)}")
        return True
    except Exception as e:
        print(f"❌ scIB metrics call failed: {e}")
        return False

if __name__ == "__main__":
    print("🧪 Testing scIB parameter fixes...")
    
    if test_scib_parameters():
        print("\n🎉 scIB parameter test passed!")
        print("The parameter names are now correct.")
    else:
        print("\n❌ scIB parameter test failed")
        sys.exit(1)
