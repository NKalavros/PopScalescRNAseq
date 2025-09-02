#!/usr/bin/env python3
"""
Test script to debug the scPoli axis index error with our fixes.
"""
import sys
import os

# Add the viz directory to the path
viz_dir = os.path.join(os.path.dirname(__file__), 'viz')
sys.path.insert(0, viz_dir)

def test_scpoli_with_fixes():
    """Test scPoli with the debugging fixes we added."""
    print("🧪 Testing scPoli with index error fixes...")
    
    try:
        # Import the main evaluation module
        import importlib.util
        spec = importlib.util.spec_from_file_location(
            "interstudy_eval", 
            os.path.join(viz_dir, "250811_evaluate_inter_studies_kidney.py")
        )
        if spec is None or spec.loader is None:
            raise ImportError("Could not create module spec")
        interstudy_module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(interstudy_module)
        print("✅ Successfully imported inter-study evaluation module")
        
        # Verify we're in debug mode
        if not interstudy_module.DEBUG_MODE:
            print("❌ DEBUG_MODE is not enabled. Please enable it in the main script.")
            return False
            
        print(f"✅ DEBUG_MODE enabled - using max {interstudy_module.DEBUG_N_CELLS} cells per study")
        
        # Load the data (just like the main script does)
        print("\n📊 Loading study data in debug mode...")
        adatas = interstudy_module.load_studies(interstudy_module.STUDIES)
        if not adatas:
            print("❌ No studies loaded")
            return False
            
        print(f"✅ Loaded {len(adatas)} studies")
        
        # Concatenate the data
        print("\n🔗 Concatenating studies...")
        combined = interstudy_module.concatenate(adatas)
        print(f"✅ Combined data shape: {combined.shape}")
        
        # Clean cell type labels
        interstudy_module.clean_cell_type_labels(combined, interstudy_module.CELL_TYPE_KEY)
        print("✅ Cell type labels cleaned")
        
        # Sanitize obs data
        interstudy_module.sanitize_obs(combined)
        print("✅ Observation data sanitized")
        
        # Now test scPoli specifically
        print("\n🔬 Testing scPoli with our fixes...")
        
        # Import scPoli module
        scpoli_spec = importlib.util.spec_from_file_location(
            'interstudy_scpoli_mod', 
            os.path.join(viz_dir, '250811_interstudy_scpoli.py')
        )
        if scpoli_spec is None or scpoli_spec.loader is None:
            raise ImportError("Could not create scPoli module spec")
        scpoli_mod = importlib.util.module_from_spec(scpoli_spec)
        scpoli_spec.loader.exec_module(scpoli_mod)
        print("✅ Successfully imported scPoli module with fixes")
        
        # Check if scPoli is available
        if scpoli_mod.scPoli is None:
            print("❌ scPoli not available in environment")
            print(f"Import error: {scpoli_mod._SCPOLI_IMPORT_ERROR}")
            return False
            
        print("✅ scPoli library is available")
        
        # Test with debug parameters (smaller, faster)
        print("\n🚀 Running scPoli with debug parameters...")
        latent_dim = 10  # Very small for testing
        pretraining_epochs = 5  # Very small for testing  
        total_epochs = 10  # Very small for testing
        
        try:
            scpoli_metrics_df, scpoli_ref = scpoli_mod.run_scpoli_interstudy(
                combined,
                study_key=interstudy_module.STUDY_KEY,
                batch_key=interstudy_module.LIBRARY_KEY,
                cell_type_key=interstudy_module.CELL_TYPE_KEY,
                embedding_key='X_scpoli',
                latent_dim=latent_dim,
                pretraining_epochs=pretraining_epochs,
                total_epochs=total_epochs,
                max_queries=2,  # Limit to 2 queries for testing
            )
            
            print(f"🎉 scPoli completed successfully!")
            print(f"📊 Metrics shape: {scpoli_metrics_df.shape}")
            print(f"📍 Reference study: {scpoli_ref}")
            print(f"✅ No axis index errors occurred!")
            
            # Check if embedding was added
            if 'X_scpoli' in combined.obsm:
                print(f"✅ scPoli embedding added to adata.obsm: {combined.obsm['X_scpoli'].shape}")
            else:
                print("⚠️ scPoli embedding not found in adata.obsm")
            
            return True
            
        except Exception as e:
            print(f"❌ scPoli failed with error: {e}")
            import traceback
            print("\n📋 Full error traceback:")
            traceback.print_exc()
            return False
            
    except Exception as e:
        print(f"❌ Test setup failed: {e}")
        import traceback
        traceback.print_exc()
        return False

if __name__ == "__main__":
    print("🧪 Testing scPoli axis index error fixes...")
    
    success = test_scpoli_with_fixes()
    
    if success:
        print("\n🎉 All tests passed! The scPoli fixes appear to work.")
        print("💡 The axis index error should now be resolved.")
    else:
        print("\n❌ Tests failed. The error may still exist or there are other issues.")
        sys.exit(1)
