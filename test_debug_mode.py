#!/usr/bin/env python3
"""
Quick test script to verify debug mode functionality in the inter-study evaluation script.
This will test the debug mode configuration and output filename generation.
"""
import sys
import os

# Add the viz directory to the path
viz_dir = os.path.join(os.path.dirname(__file__), 'viz')
sys.path.insert(0, viz_dir)

# Import the module (this will also test for syntax errors)
try:
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
except Exception as e:
    print(f"❌ Failed to import module: {e}")
    sys.exit(1)

# Test debug mode configuration
print(f"\n📊 Debug mode configuration:")
print(f"  DEBUG_MODE: {interstudy_module.DEBUG_MODE}")
print(f"  DEBUG_N_CELLS: {interstudy_module.DEBUG_N_CELLS}")

# Test output filename generation
print(f"\n📁 Output filenames:")
print(f"  Combined data: {interstudy_module.OUTPUT_COMBINED_FILENAME}")
print(f"  Label mapping: {interstudy_module.LABEL_MAPPING_JSON}")
print(f"  Symphony metrics: {interstudy_module.SYMPHONY_METRICS_CSV}")
print(f"  scPoli metrics: {interstudy_module.SCPOLI_METRICS_CSV}")
print(f"  scArches metrics: {interstudy_module.SCARCHES_METRICS_CSV}")
print(f"  scIB metrics: {interstudy_module.SCIB_METRICS_CSV}")

# Test that all debug filenames have the debug suffix
expected_debug_files = [
    'interstudy_combined_debug.h5ad',
    'interstudy_celltype_mapping_debug.json',
    'symphony_interstudy_metrics_debug.csv',
    'scpoli_interstudy_metrics_debug.csv',
    'scarches_interstudy_metrics_debug.csv',
    'scib_metrics_long_debug.csv'
]

actual_files = [
    interstudy_module.OUTPUT_COMBINED_FILENAME,
    interstudy_module.LABEL_MAPPING_JSON,
    interstudy_module.SYMPHONY_METRICS_CSV,
    interstudy_module.SCPOLI_METRICS_CSV,
    interstudy_module.SCARCHES_METRICS_CSV,
    interstudy_module.SCIB_METRICS_CSV
]

print(f"\n🧪 Testing debug filename logic:")
all_correct = True
for expected, actual in zip(expected_debug_files, actual_files):
    if expected == actual:
        print(f"  ✅ {actual}")
    else:
        print(f"  ❌ Expected: {expected}, Got: {actual}")
        all_correct = False

if all_correct:
    print(f"\n🎉 All tests passed! Debug mode is correctly configured.")
else:
    print(f"\n❌ Some tests failed. Check the configuration.")
    sys.exit(1)

# Test that we can access the load_studies function
if hasattr(interstudy_module, 'load_studies'):
    print(f"\n✅ load_studies function is accessible")
else:
    print(f"\n❌ load_studies function not found")
    sys.exit(1)

print(f"\n📋 Studies to load: {interstudy_module.STUDIES}")
print(f"📂 Base directory: {interstudy_module.BASE_DIR}")
print(f"\n🚀 Ready to run inter-study evaluation in debug mode!")
