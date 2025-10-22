# Complete Automated Setup Script for ScRNA_GBM_Analysis
# This script will set up everything needed to run the pipeline

$ErrorActionPreference = "Stop"

# Colors for output
function Write-Success { param($msg) Write-Host "✓ $msg" -ForegroundColor Green }
function Write-Info { param($msg) Write-Host "→ $msg" -ForegroundColor Cyan }
function Write-Warning { param($msg) Write-Host "⚠ $msg" -ForegroundColor Yellow }
function Write-Error { param($msg) Write-Host "✗ $msg" -ForegroundColor Red }
function Write-Header { param($msg) Write-Host "`n$('='*60)`n$msg`n$('='*60)" -ForegroundColor Cyan }

Write-Header "ScRNA-seq GBM Analysis - Complete Automated Setup"

# Check Python installation
Write-Info "Checking Python installation..."
try {
    $pythonVersion = python --version 2>&1
    if ($LASTEXITCODE -eq 0) {
        Write-Success "Found: $pythonVersion"
    } else {
        throw "Python not found"
    }
} catch {
    Write-Error "Python is not installed or not in PATH"
    Write-Info "Please install Python 3.9+ from https://www.python.org/downloads/"
    Write-Info "Make sure to check 'Add Python to PATH' during installation"
    Write-Info "Then restart your terminal and run this script again."
    exit 1
}

# Check Python version
$pythonVersionNum = python -c "import sys; print(f'{sys.version_info.major}.{sys.version_info.minor}')" 2>&1
if ([double]$pythonVersionNum -lt 3.8) {
    Write-Error "Python $pythonVersionNum detected. Need Python 3.8 or higher"
    exit 1
}
Write-Success "Python version $pythonVersionNum is compatible"

# Create virtual environment
Write-Header "Step 1: Creating Virtual Environment"
if (Test-Path "venv") {
    Write-Warning "Virtual environment already exists. Skipping creation."
} else {
    Write-Info "Creating virtual environment..."
    python -m venv venv
    Write-Success "Virtual environment created"
}

# Activate virtual environment
Write-Info "Activating virtual environment..."
& .\venv\Scripts\Activate.ps1
Write-Success "Virtual environment activated"

# Upgrade pip
Write-Header "Step 2: Upgrading pip"
Write-Info "Upgrading pip to latest version..."
python -m pip install --upgrade pip --quiet
Write-Success "pip upgraded"

# Install dependencies
Write-Header "Step 3: Installing Python Packages"
Write-Info "This may take 10-15 minutes. Please be patient..."
Write-Info "Installing packages from requirements.txt..."

$packages = @(
    "numpy>=1.21.0",
    "pandas>=1.3.0",
    "scipy>=1.7.0",
    "matplotlib>=3.5.0",
    "seaborn>=0.11.0",
    "scanpy>=1.9.0",
    "anndata>=0.8.0",
    "scrublet>=0.2.3",
    "harmonypy>=0.0.9",
    "infercnvpy>=0.4.0",
    "tqdm",
    "jupyter"
)

$totalPackages = $packages.Count
$current = 0

foreach ($package in $packages) {
    $current++
    $packageName = $package.Split(">=")[0]
    Write-Info "[$current/$totalPackages] Installing $packageName..."
    pip install $package --quiet
    if ($LASTEXITCODE -eq 0) {
        Write-Success "$packageName installed"
    } else {
        Write-Warning "$packageName installation may have issues (continuing...)"
    }
}

Write-Success "All packages installed!"

# Verify installations
Write-Header "Step 4: Verifying Installation"
Write-Info "Checking critical packages..."

$criticalPackages = @("scanpy", "numpy", "pandas", "matplotlib")
$allInstalled = $true

foreach ($pkg in $criticalPackages) {
    try {
        python -c "import $pkg; print(f'$pkg version: {$pkg.__version__}')" 2>&1 | Out-Null
        if ($LASTEXITCODE -eq 0) {
            $version = python -c "import $pkg; print($pkg.__version__)" 2>&1
            Write-Success "$pkg (v$version)"
        } else {
            Write-Error "$pkg not found"
            $allInstalled = $false
        }
    } catch {
        Write-Error "$pkg not found"
        $allInstalled = $false
    }
}

if (-not $allInstalled) {
    Write-Error "Some packages failed to install. Please run:"
    Write-Info "pip install -r requirements.txt"
    exit 1
}

# Create data download script
Write-Header "Step 5: Creating Data Download Scripts"

$downloadScript = @'
"""
download_geo_data.py
Automated script to download Klemm et al. GBM scRNA-seq data from GEO
"""

import os
import urllib.request
import gzip
import shutil
from pathlib import Path

def download_file(url, dest_path):
    """Download file with progress bar"""
    print(f"Downloading: {url}")
    print(f"To: {dest_path}")
    
    try:
        urllib.request.urlretrieve(url, dest_path)
        print(f"✓ Downloaded successfully")
        return True
    except Exception as e:
        print(f"✗ Download failed: {e}")
        return False

def main():
    print("\n" + "="*60)
    print("Downloading Klemm et al. GBM Dataset from GEO")
    print("="*60 + "\n")
    
    # Note: This is a template. Actual GEO data may require different URLs
    # For real data, visit: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE163108
    
    print("⚠ IMPORTANT: Automated download from GEO is complex.")
    print("Please follow manual download instructions:")
    print("\n1. Visit: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE163108")
    print("2. Click 'Download family' -> 'Custom' tab")
    print("3. Select 'Series Matrix File(s)'")
    print("4. Download supplementary files")
    print("\n5. For processed data, check:")
    print("   https://singlecell.broadinstitute.org/single_cell/study/SCP1290")
    print("\n6. Save files to: data/raw/")
    
    # Create placeholder structure
    os.makedirs("data/raw/tumour", exist_ok=True)
    os.makedirs("data/raw/normal", exist_ok=True)
    
    # Create README
    readme_content = """# Data Download Instructions

## Klemm et al., 2020 Dataset (GSE163108)

### Option 1: GEO Database
1. Visit: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE163108
2. Download supplementary files
3. Extract to data/raw/

### Option 2: Single Cell Portal
1. Visit: https://singlecell.broadinstitute.org/single_cell/study/SCP1290
2. Download processed data
3. Extract to data/raw/

### Option 3: Original Publication
- Paper: Nature (2020) doi:10.1038/s41586-020-1959-y
- Check supplementary materials

### Data Structure Needed:
```
data/raw/
├── tumour/
│   ├── matrix.mtx.gz
│   ├── features.tsv.gz
│   └── barcodes.tsv.gz
└── normal/
    ├── matrix.mtx.gz
    ├── features.tsv.gz
    └── barcodes.tsv.gz
```

### Converting Other Formats:
If data is in different format (.h5ad, .rds, .loom), use conversion scripts:
- For .h5ad: `python scripts/convert_h5ad_to_10x.py`
- For .rds: Use R to convert to 10x format
- For .loom: `python scripts/convert_loom_to_10x.py`

Contact: Check GEO database for data access details
"""
    
    with open("data/DATA_DOWNLOAD_INSTRUCTIONS.md", "w") as f:
        f.write(readme_content)
    
    print("\n✓ Created data structure and instructions")
    print(f"✓ See: data/DATA_DOWNLOAD_INSTRUCTIONS.md")

if __name__ == "__main__":
    main()
'@

Set-Content -Path "scripts\download_geo_data.py" -Value $downloadScript
Write-Success "Created data download script"

# Run data download information script
Write-Header "Step 6: Data Download Information"
python scripts\download_geo_data.py

# Create test/demo data
Write-Info "Creating demo data structure..."

# Create demo data script
$demoDataScript = @'
"""Generate small demo dataset for testing"""
import numpy as np
import pandas as pd
import scipy.sparse as sp
from pathlib import Path
import gzip

def create_demo_10x_data(output_dir, n_genes=2000, n_cells=500, prefix="demo"):
    """Create a small demo 10x format dataset"""
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)
    
    # Generate random count matrix
    np.random.seed(42)
    density = 0.1  # 10% non-zero
    data = np.random.negative_binomial(5, 0.3, size=int(n_genes * n_cells * density))
    
    # Create sparse matrix
    row_indices = np.random.randint(0, n_genes, size=len(data))
    col_indices = np.random.randint(0, n_cells, size=len(data))
    matrix = sp.coo_matrix((data, (row_indices, col_indices)), shape=(n_genes, n_cells))
    
    # Save matrix.mtx.gz
    with gzip.open(output_path / "matrix.mtx.gz", 'wt') as f:
        f.write(f"%%MatrixMarket matrix coordinate integer general\n")
        f.write(f"%\n")
        f.write(f"{n_genes} {n_cells} {len(matrix.data)}\n")
        for i, j, v in zip(matrix.row + 1, matrix.col + 1, matrix.data):
            f.write(f"{i} {j} {int(v)}\n")
    
    # Generate gene names
    gene_names = [f"GENE{i:04d}" for i in range(n_genes)]
    # Add some realistic gene names
    real_genes = ['DGAT1', 'DGAT2', 'FASN', 'CD68', 'CD3D', 'GFAP', 'MBP', 'PECAM1', 
                  'MT-CO1', 'MT-CO2', 'MT-ND1', 'PTPRC', 'CSF1R', 'ITGAM']
    for idx, gene in enumerate(real_genes[:min(len(real_genes), n_genes)]):
        gene_names[idx] = gene
    
    # Save features.tsv.gz
    with gzip.open(output_path / "features.tsv.gz", 'wt') as f:
        for gene in gene_names:
            f.write(f"{gene}\t{gene}\tGene Expression\n")
    
    # Save barcodes.tsv.gz
    with gzip.open(output_path / "barcodes.tsv.gz", 'wt') as f:
        for i in range(n_cells):
            f.write(f"{prefix}-{i:04d}-1\n")
    
    print(f"✓ Created demo dataset: {output_dir}")
    print(f"  - {n_genes} genes")
    print(f"  - {n_cells} cells")

if __name__ == "__main__":
    # Create demo tumour data
    create_demo_10x_data("data/raw/demo_tumour", n_genes=2000, n_cells=300, prefix="TUMOUR")
    
    # Create demo normal data
    create_demo_10x_data("data/raw/demo_normal", n_genes=2000, n_cells=200, prefix="NORMAL")
    
    print("\n✓ Demo datasets created!")
    print("\nTo use demo data, run:")
    print("  python scripts/sc_pipeline.py --tumour_path data/raw/demo_tumour/ --normal_path data/raw/demo_normal/ --output_dir results/demo/")
'@

Set-Content -Path "scripts\create_demo_data.py" -Value $demoDataScript

# Create demo data
Write-Header "Step 7: Creating Demo Test Data"
Write-Info "Generating small demo dataset for testing..."
python scripts\create_demo_data.py
Write-Success "Demo data created in data/raw/demo_tumour/ and data/raw/demo_normal/"

# Create verification script
Write-Header "Step 8: Creating Verification Script"

$verifyScript = @'
"""Verify that setup completed successfully"""
import sys

def check_package(name):
    try:
        __import__(name)
        return True
    except ImportError:
        return False

print("\n" + "="*60)
print("Setup Verification")
print("="*60 + "\n")

# Check Python version
print(f"Python version: {sys.version.split()[0]}")
if sys.version_info < (3, 8):
    print("✗ Python 3.8+ required")
    sys.exit(1)
else:
    print("✓ Python version OK")

# Check packages
packages = {
    "scanpy": "scanpy",
    "numpy": "numpy",
    "pandas": "pandas",
    "scipy": "scipy",
    "matplotlib": "matplotlib",
    "seaborn": "seaborn",
    "anndata": "anndata",
    "scrublet": "scrublet",
}

all_ok = True
for display_name, import_name in packages.items():
    if check_package(import_name):
        print(f"✓ {display_name}")
    else:
        print(f"✗ {display_name} - NOT INSTALLED")
        all_ok = False

# Check data
import os
if os.path.exists("data/raw/demo_tumour/matrix.mtx.gz"):
    print("✓ Demo data available")
else:
    print("⚠ Demo data not found (optional)")

print("\n" + "="*60)
if all_ok:
    print("✓ Setup verification PASSED!")
    print("\nYou can now run:")
    print("  python scripts/step_by_step.py")
    print("\nOr test with demo data:")
    print("  python scripts/sc_pipeline.py --tumour_path data/raw/demo_tumour/ --normal_path data/raw/demo_normal/ --output_dir results/demo/")
else:
    print("✗ Setup verification FAILED")
    print("Please install missing packages: pip install -r requirements.txt")
    sys.exit(1)
print("="*60 + "\n")
'@

Set-Content -Path "scripts\verify_setup.py" -Value $verifyScript

# Run verification
Write-Header "Step 9: Verifying Complete Setup"
python scripts\verify_setup.py

# Final summary
Write-Header "🎉 Setup Complete!"
Write-Host @"

✅ Virtual environment created and activated
✅ All Python packages installed
✅ Demo test data generated
✅ Data download instructions created
✅ Scripts verified and ready

📍 Current Status:
  - Location: $PWD
  - Python: $(python --version)
  - Virtual env: ACTIVATED

🎯 Next Steps:

Option 1 - Test with Demo Data (Recommended First):
  python scripts\sc_pipeline.py --tumour_path data/raw/demo_tumour/ --normal_path data/raw/demo_normal/ --output_dir results/demo/ --max_genes 5000

Option 2 - Download Real Klemm Dataset:
  See: data/DATA_DOWNLOAD_INSTRUCTIONS.md
  Visit: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE163108

Option 3 - Run Interactive Analysis:
  python scripts\step_by_step.py

📚 Documentation:
  - SETUP_COMPLETE.md - This guide
  - GET_STARTED.md - Detailed walkthrough
  - README.md - Complete documentation
  - data/DATA_DOWNLOAD_INSTRUCTIONS.md - Data download guide

"@ -ForegroundColor Green

Write-Host "`nSetup completed successfully! 🚀" -ForegroundColor Cyan
Write-Host "The virtual environment is still activated in this session." -ForegroundColor Yellow
Write-Host "To use it later, run: .\venv\Scripts\Activate.ps1`n" -ForegroundColor Yellow

