# Progress Monitor for DGAT1-Inclusive Analysis
# Monitors the re-run with 10,000 HVGs

Write-Host "`n" -NoNewline
Write-Host "="*70 -ForegroundColor Green
Write-Host "DGAT1-INCLUSIVE ANALYSIS - PROGRESS MONITOR" -ForegroundColor Green
Write-Host "="*70 -ForegroundColor Green
Write-Host ""

$outputDir = "results\with_dgat1"
$expectedFile = "all_samples_processed.h5ad"

Write-Host "[TARGET] Getting DGAT1 by keeping 10,000 HVGs (vs 3,000)`n" -ForegroundColor Yellow

if (-not (Test-Path $outputDir)) {
    Write-Host "[INFO] Output directory not created yet" -ForegroundColor Yellow
    Write-Host "Analysis is starting...`n" -ForegroundColor Yellow
    exit
}

Write-Host "[CHECKING] Output: $outputDir`n" -ForegroundColor Yellow

# Check for main file
if (Test-Path (Join-Path $outputDir $expectedFile)) {
    $mainFile = Get-Item (Join-Path $outputDir $expectedFile)
    $sizeMB = [math]::Round($mainFile.Length / 1MB, 2)
    $sizeGB = [math]::Round($mainFile.Length / 1GB, 2)
    
    Write-Host "="*70 -ForegroundColor Green
    Write-Host "[SUCCESS] ANALYSIS COMPLETE!" -ForegroundColor Green
    Write-Host "="*70 -ForegroundColor Green
    
    Write-Host "`nDataset ready:" -ForegroundColor Green
    Write-Host "  File: $expectedFile" -ForegroundColor White
    Write-Host "  Size: $sizeGB GB ($sizeMB MB)" -ForegroundColor White
    Write-Host "  Location: $outputDir\" -ForegroundColor White
    Write-Host "  Modified: $($mainFile.LastWriteTime)" -ForegroundColor White
    
    Write-Host "`n[DGAT1] Checking if included..." -ForegroundColor Cyan
    Write-Host "  Running quick check..." -ForegroundColor Yellow
    
    # Quick Python check for DGAT1
    $checkCmd = @"
import scanpy as sc
adata = sc.read_h5ad('$outputDir/$expectedFile')
dgat1_present = 'DGAT1' in adata.var_names
print(f'DGAT1 in dataset: {dgat1_present}')
if dgat1_present:
    print(f'Cells: {adata.n_obs:,}')
    print(f'Genes: {adata.n_vars:,}')
    print('SUCCESS: DGAT1 is available for analysis!')
else:
    print('WARNING: DGAT1 still not in dataset')
"@
    
    py -c $checkCmd
    
    Write-Host "`nNext steps:" -ForegroundColor Cyan
    Write-Host "  1. Generate individual figures with DGAT1:" -ForegroundColor White
    Write-Host "     py scripts\generate_individual_figures.py --input $outputDir\$expectedFile --output results\final_figures" -ForegroundColor Gray
    Write-Host "`n  2. Analyze DGAT1 expression:" -ForegroundColor White
    Write-Host "     py scripts\analyze_dgat1.py --input $outputDir\$expectedFile" -ForegroundColor Gray
    
} else {
    Write-Host "[RUNNING] Analysis in progress...`n" -ForegroundColor Yellow
    
    # Check for any output files
    $files = Get-ChildItem -Path $outputDir -File -ErrorAction SilentlyContinue
    
    if ($files.Count -gt 0) {
        Write-Host "[FOUND] Intermediate files:" -ForegroundColor Green
        foreach ($file in $files) {
            $sizeMB = [math]::Round($file.Length / 1MB, 2)
            $age = (Get-Date) - $file.LastWriteTime
            $ageMin = [math]::Round($age.TotalMinutes, 1)
            Write-Host "  - $($file.Name) ($sizeMB MB) [$ageMin min ago]" -ForegroundColor White
        }
    } else {
        Write-Host "[INFO] No files yet - analysis just started" -ForegroundColor Yellow
    }
    
    Write-Host "`nEstimated stages (total ~40-50 min):" -ForegroundColor Cyan
    Write-Host "  1. Loading 18 samples      (~5 min)" -ForegroundColor Gray
    Write-Host "  2. Quality control         (~3 min)" -ForegroundColor Gray
    Write-Host "  3. Normalization           (~5 min)" -ForegroundColor Gray
    Write-Host "  4. Find 10,000 HVGs        (~8 min) [DGAT1 included!]" -ForegroundColor Yellow
    Write-Host "  5. Batch correction        (~15 min) [Longest step]" -ForegroundColor Gray
    Write-Host "  6. Clustering & UMAP       (~5 min)" -ForegroundColor Gray
    Write-Host "  7. Cell type annotation    (~3 min)" -ForegroundColor Gray
    
    Write-Host "`nRun this script again to check progress" -ForegroundColor Cyan
}

Write-Host ""
$now = Get-Date -Format "HH:mm:ss"
Write-Host "Last checked: $now`n" -ForegroundColor DarkGray

