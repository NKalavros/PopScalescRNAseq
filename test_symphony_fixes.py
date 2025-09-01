#!/usr/bin/env python3
"""
Test script to validate our Symphony fixes.
"""
import sys
import os

# Add the viz directory to the path
viz_dir = os.path.join(os.path.dirname(__file__), 'viz')
sys.path.insert(0, viz_dir)

def test_symphony_import():
    """Test that we can import the Symphony module."""
    try:
        import importlib.util
        spec = importlib.util.spec_from_file_location(
            "interstudy_symphony_mod", 
            os.path.join(viz_dir, "250811_interstudy_symphony.py")
        )
        if spec is None or spec.loader is None:
            raise ImportError("Could not create module spec")
        symphony_mod = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(symphony_mod)
        print("✅ Successfully imported Symphony module")
        
        # Check if we have the function
        if hasattr(symphony_mod, 'run_symphony_interstudy'):
            print("✅ run_symphony_interstudy function found")
        else:
            print("❌ run_symphony_interstudy function not found")
            return False
            
        return True
    except Exception as e:
        print(f"❌ Failed to import Symphony module: {e}")
        return False

def test_symphony_requirements():
    """Test Symphony requirements."""
    try:
        import symphonypy as sp
        print("✅ symphonypy is available")
    except ImportError as e:
        print(f"❌ symphonypy not available: {e}")
        return False
        
    try:
        from sklearn.metrics import f1_score, accuracy_score
        print("✅ scikit-learn metrics available")
    except ImportError as e:
        print(f"❌ scikit-learn metrics not available: {e}")
        return False
        
    return True

if __name__ == "__main__":
    print("🧪 Testing Symphony fixes...")
    
    if not test_symphony_requirements():
        print("\n❌ Symphony requirements not met")
        sys.exit(1)
        
    if not test_symphony_import():
        print("\n❌ Symphony import failed")
        sys.exit(1)
        
    print("\n🎉 All Symphony tests passed!")
    print("Ready to test with actual data.")
