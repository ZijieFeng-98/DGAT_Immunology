#!/usr/bin/env python3
"""
Installation Test Script for GBM Single-Cell Pipeline
======================================================

Run this script to verify that all required packages are installed
and functioning correctly.

Usage:
    python test_installation.py
"""

import sys
from typing import Tuple, List

def test_package(package_name: str, import_name: str = None) -> Tuple[bool, str]:
    """Test if a package can be imported."""
    if import_name is None:
        import_name = package_name
    
    try:
        __import__(import_name)
        return True, "✓ Installed"
    except ImportError as e:
        return False, f"✗ Not installed: {str(e)}"

def main():
    """Run all installation tests."""
    print("="*70)
    print("GBM Single-Cell Pipeline - Installation Test")
    print("="*70)
    print()
    
    # Required packages
    required_packages = [
        ('scanpy', 'scanpy'),
        ('anndata', 'anndata'),
        ('numpy', 'numpy'),
        ('pandas', 'pandas'),
        ('scipy', 'scipy'),
        ('matplotlib', 'matplotlib'),
        ('seaborn', 'seaborn'),
    ]
    
    # Optional but recommended packages
    optional_packages = [
        ('scrublet', 'scrublet'),
        ('harmonypy', 'harmonypy'),
        ('infercnvpy', 'infercnvpy'),
        ('liana-py', 'liana'),
        ('lifelines', 'lifelines'),
        ('scikit-learn', 'sklearn'),
    ]
    
    print("REQUIRED PACKAGES:")
    print("-" * 70)
    
    required_results = []
    for pkg_name, import_name in required_packages:
        success, message = test_package(pkg_name, import_name)
        required_results.append(success)
        status = "✓" if success else "✗"
        print(f"  {status} {pkg_name:<20} {message}")
    
    print()
    print("OPTIONAL PACKAGES:")
    print("-" * 70)
    
    optional_results = []
    for pkg_name, import_name in optional_packages:
        success, message = test_package(pkg_name, import_name)
        optional_results.append(success)
        status = "✓" if success else "✗"
        print(f"  {status} {pkg_name:<20} {message}")
    
    print()
    print("="*70)
    
    # Summary
    req_passed = sum(required_results)
    req_total = len(required_results)
    opt_passed = sum(optional_results)
    opt_total = len(optional_results)
    
    print("SUMMARY:")
    print(f"  Required: {req_passed}/{req_total} packages installed")
    print(f"  Optional: {opt_passed}/{opt_total} packages installed")
    print()
    
    if req_passed == req_total:
        print("✓ All required packages installed! You're ready to run the pipeline.")
        
        if opt_passed < opt_total:
            print()
            print("Note: Some optional packages are missing. You can still run the")
            print("basic pipeline, but some advanced features may be unavailable.")
            print()
            print("To install missing optional packages:")
            missing = [pkg for pkg, _ in optional_packages if not optional_results[optional_packages.index((pkg, _))]]
            for pkg in missing:
                print(f"  pip install {pkg}")
        
        return 0
    else:
        print("✗ Some required packages are missing. Please install them:")
        print()
        print("  pip install -r requirements.txt")
        print()
        print("Or install individually:")
        missing = [pkg for pkg, _ in required_packages if not required_results[required_packages.index((pkg, _))]]
        for pkg in missing:
            print(f"  pip install {pkg}")
        
        return 1
    
def test_scanpy_functionality():
    """Test basic scanpy functionality."""
    print()
    print("="*70)
    print("TESTING SCANPY FUNCTIONALITY")
    print("="*70)
    
    try:
        import scanpy as sc
        import numpy as np
        
        # Create a small test dataset
        print("Creating test dataset...")
        n_obs = 100
        n_vars = 50
        X = np.random.negative_binomial(5, 0.3, (n_obs, n_vars))
        
        adata = sc.AnnData(X)
        adata.obs['sample'] = ['sample1'] * 50 + ['sample2'] * 50
        adata.var_names = [f'Gene{i}' for i in range(n_vars)]
        
        # Test basic operations
        print("Testing normalization...")
        sc.pp.normalize_total(adata, target_sum=1e4)
        
        print("Testing log transformation...")
        sc.pp.log1p(adata)
        
        print("Testing HVG selection...")
        sc.pp.highly_variable_genes(adata, n_top_genes=20)
        
        print("Testing PCA...")
        sc.tl.pca(adata, n_comps=10)
        
        print("Testing neighbors...")
        sc.pp.neighbors(adata, n_neighbors=10, n_pcs=10)
        
        print("Testing UMAP...")
        sc.tl.umap(adata)
        
        print("Testing clustering...")
        sc.tl.leiden(adata, resolution=0.5)
        
        print()
        print("✓ All Scanpy functionality tests passed!")
        print()
        print(f"Test dataset: {adata.n_obs} cells x {adata.n_vars} genes")
        print(f"Clusters found: {len(adata.obs['leiden'].unique())}")
        
        return True
        
    except Exception as e:
        print()
        print(f"✗ Scanpy functionality test failed: {str(e)}")
        print()
        return False

def test_visualization():
    """Test visualization capabilities."""
    print()
    print("="*70)
    print("TESTING VISUALIZATION")
    print("="*70)
    
    try:
        import matplotlib
        import matplotlib.pyplot as plt
        import seaborn as sns
        import numpy as np
        
        # Check matplotlib backend
        print(f"Matplotlib backend: {matplotlib.get_backend()}")
        
        # Create simple test plot
        print("Creating test plot...")
        fig, ax = plt.subplots(figsize=(6, 4))
        x = np.linspace(0, 10, 100)
        y = np.sin(x)
        ax.plot(x, y)
        ax.set_title("Test Plot")
        plt.close(fig)
        
        print("✓ Visualization test passed!")
        
        return True
        
    except Exception as e:
        print(f"✗ Visualization test failed: {str(e)}")
        return False

if __name__ == '__main__':
    # Run installation tests
    exit_code = main()
    
    if exit_code == 0:
        # Run functionality tests
        scanpy_ok = test_scanpy_functionality()
        viz_ok = test_visualization()
        
        print()
        print("="*70)
        print("FINAL RESULT")
        print("="*70)
        
        if scanpy_ok and viz_ok:
            print("✓ All tests passed! Your installation is ready.")
            print()
            print("Next steps:")
            print("  1. Prepare your 10x data")
            print("  2. Run: python sc_glioblastoma_dgat1_pipeline.py --help")
            print("  3. Check QUICKSTART.md for examples")
        else:
            print("⚠ Installation is incomplete or some tests failed.")
            print("  The pipeline may still work, but with limited functionality.")
        
    sys.exit(exit_code)
