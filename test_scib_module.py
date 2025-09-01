#!/usr/bin/env python3
"""
Test script to validate our scIB integration module.
"""
import sys
import os

# Add the viz directory to the path
viz_dir = os.path.join(os.path.dirname(__file__), 'viz')
sys.path.insert(0, viz_dir)

def test_scib_import():
    """Test that we can import the scIB module."""
    try:
        import importlib.util
        spec = importlib.util.spec_from_file_location(
            "interstudy_scib_mod", 
            os.path.join(viz_dir, "250811_interstudy_scib.py")
        )
        if spec is None or spec.loader is None:
            raise ImportError("Could not create module spec")
        scib_mod = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(scib_mod)
        print("✅ Successfully imported scIB module")
        
        # Check if we have the function
        if hasattr(scib_mod, 'run_scib_benchmark'):
            print("✅ run_scib_benchmark function found")
        else:
            print("❌ run_scib_benchmark function not found")
            return False
            
        return True
    except Exception as e:
        print(f"❌ Failed to import scIB module: {e}")
        return False

def test_metric_flags():
    """Test that the metric flags are correctly formatted."""
    try:
        import importlib.util
        spec = importlib.util.spec_from_file_location(
            "interstudy_scib_mod", 
            os.path.join(viz_dir, "250811_interstudy_scib.py")
        )
        scib_mod = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(scib_mod)
        
        # We can't directly access METRIC_FLAGS_BASE as it's inside the function,
        # but we can examine the source code to verify our changes
        with open(os.path.join(viz_dir, "250811_interstudy_scib.py"), 'r') as f:
            content = f.read()
            
        # Check for the corrected parameter names with underscores
        corrected_params = [
            'ari_=True',
            'nmi_=True', 
            'silhouette_=True',
            'graph_conn_=True',
            'pcr_=True',
            'ilisi_=',
            'clisi_=',
            'kBET_=',
            'isolated_labels_f1_=True',
            'isolated_labels_asw_=True',
            'hvg_score_=False',
            'cell_cycle_=False',
            'trajectory_=False'
        ]
        
        missing_params = []
        for param in corrected_params:
            if param not in content:
                missing_params.append(param)
                
        if missing_params:
            print(f"❌ Missing corrected parameters: {missing_params}")
            return False
        else:
            print("✅ All parameter names have been corrected with underscores")
            return True
            
    except Exception as e:
        print(f"❌ Failed to check metric flags: {e}")
        return False

if __name__ == "__main__":
    print("🧪 Testing scIB integration module...")
    
    success = True
    
    if not test_scib_import():
        success = False
        
    if not test_metric_flags():
        success = False
        
    if success:
        print("\n🎉 All scIB module tests passed!")
        print("The parameter names have been corrected to use underscores.")
    else:
        print("\n❌ scIB module tests failed")
        sys.exit(1)
