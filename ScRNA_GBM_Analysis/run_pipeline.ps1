# PowerShell script to run the single-cell pipeline
# Usage: .\run_pipeline.ps1

Write-Host "=" -NoNewline -ForegroundColor Cyan
Write-Host "="*60 -ForegroundColor Cyan
Write-Host "Single-Cell GBM Analysis Pipeline" -ForegroundColor Green
Write-Host "=" -NoNewline -ForegroundColor Cyan
Write-Host "="*60 -ForegroundColor Cyan

# Check if Python is installed
Write-Host "`nChecking Python installation..." -ForegroundColor Yellow
try {
    $pythonVersion = python --version 2>&1
    Write-Host "✓ $pythonVersion" -ForegroundColor Green
} catch {
    Write-Host "✗ Python not found. Please install Python 3.8+ first." -ForegroundColor Red
    exit 1
}

# Check if virtual environment exists
if (!(Test-Path "venv")) {
    Write-Host "`nCreating virtual environment..." -ForegroundColor Yellow
    python -m venv venv
    Write-Host "✓ Virtual environment created" -ForegroundColor Green
}

# Activate virtual environment
Write-Host "`nActivating virtual environment..." -ForegroundColor Yellow
& .\venv\Scripts\Activate.ps1
Write-Host "✓ Virtual environment activated" -ForegroundColor Green

# Check if dependencies are installed
Write-Host "`nChecking dependencies..." -ForegroundColor Yellow
$scanpyInstalled = pip list 2>&1 | Select-String "scanpy"
if (!$scanpyInstalled) {
    Write-Host "Installing dependencies (this may take a few minutes)..." -ForegroundColor Yellow
    pip install -r requirements.txt
    Write-Host "✓ Dependencies installed" -ForegroundColor Green
} else {
    Write-Host "✓ Dependencies already installed" -ForegroundColor Green
}

# Check if data exists
Write-Host "`nChecking for data..." -ForegroundColor Yellow
if (!(Test-Path "data\raw\tumour") -or !(Test-Path "data\raw\normal")) {
    Write-Host "✗ Data not found in data/raw/" -ForegroundColor Red
    Write-Host "`nPlease organize your 10x data as:" -ForegroundColor Yellow
    Write-Host "  data/raw/tumour/  (matrix.mtx.gz, features.tsv.gz, barcodes.tsv.gz)"
    Write-Host "  data/raw/normal/  (matrix.mtx.gz, features.tsv.gz, barcodes.tsv.gz)"
    Write-Host "`nWould you like to continue anyway? (y/n): " -NoNewline
    $response = Read-Host
    if ($response -ne "y") {
        Write-Host "Exiting..." -ForegroundColor Yellow
        exit 0
    }
} else {
    Write-Host "✓ Data directories found" -ForegroundColor Green
}

# Ask user which mode to run
Write-Host "`n" + "="*60 -ForegroundColor Cyan
Write-Host "How would you like to run the pipeline?" -ForegroundColor Green
Write-Host "="*60 -ForegroundColor Cyan
Write-Host "  1. Interactive step-by-step mode (recommended)" -ForegroundColor White
Write-Host "  2. Automatic pipeline with default parameters" -ForegroundColor White
Write-Host "  3. Custom parameters (manual)" -ForegroundColor White
Write-Host "  4. Open Jupyter notebook" -ForegroundColor White
Write-Host "  q. Quit" -ForegroundColor White
Write-Host ""
$choice = Read-Host "Enter your choice (1-4 or q)"

switch ($choice) {
    "1" {
        Write-Host "`nStarting interactive step-by-step mode..." -ForegroundColor Green
        python scripts\step_by_step.py
    }
    "2" {
        Write-Host "`nStarting automatic pipeline with default parameters..." -ForegroundColor Green
        python scripts\sc_pipeline.py `
            --tumour_path data\raw\tumour `
            --normal_path data\raw\normal `
            --output_dir results `
            --min_genes 200 `
            --max_genes 2500 `
            --max_mito 10 `
            --resolution 0.5
    }
    "3" {
        Write-Host "`nEnter custom parameters:" -ForegroundColor Green
        $tumourPath = Read-Host "Tumour data path [data\raw\tumour]"
        if ([string]::IsNullOrWhiteSpace($tumourPath)) { $tumourPath = "data\raw\tumour" }
        
        $normalPath = Read-Host "Normal data path [data\raw\normal]"
        if ([string]::IsNullOrWhiteSpace($normalPath)) { $normalPath = "data\raw\normal" }
        
        $outputDir = Read-Host "Output directory [results]"
        if ([string]::IsNullOrWhiteSpace($outputDir)) { $outputDir = "results" }
        
        $minGenes = Read-Host "Minimum genes per cell [200]"
        if ([string]::IsNullOrWhiteSpace($minGenes)) { $minGenes = "200" }
        
        $maxGenes = Read-Host "Maximum genes per cell [2500]"
        if ([string]::IsNullOrWhiteSpace($maxGenes)) { $maxGenes = "2500" }
        
        $maxMito = Read-Host "Maximum mitochondrial % [10]"
        if ([string]::IsNullOrWhiteSpace($maxMito)) { $maxMito = "10" }
        
        $resolution = Read-Host "Clustering resolution [0.5]"
        if ([string]::IsNullOrWhiteSpace($resolution)) { $resolution = "0.5" }
        
        Write-Host "`nStarting pipeline with custom parameters..." -ForegroundColor Green
        python scripts\sc_pipeline.py `
            --tumour_path $tumourPath `
            --normal_path $normalPath `
            --output_dir $outputDir `
            --min_genes $minGenes `
            --max_genes $maxGenes `
            --max_mito $maxMito `
            --resolution $resolution
    }
    "4" {
        Write-Host "`nStarting Jupyter Notebook..." -ForegroundColor Green
        jupyter notebook notebooks\
    }
    "q" {
        Write-Host "`nExiting..." -ForegroundColor Yellow
        exit 0
    }
    default {
        Write-Host "`n✗ Invalid choice. Exiting..." -ForegroundColor Red
        exit 1
    }
}

Write-Host "`n" + "="*60 -ForegroundColor Cyan
Write-Host "Pipeline Complete!" -ForegroundColor Green
Write-Host "="*60 -ForegroundColor Cyan
Write-Host "`nCheck the results/ directory for outputs." -ForegroundColor White

