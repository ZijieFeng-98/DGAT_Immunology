# Progress Monitor for All 21 Samples Analysis
# Monitors the comprehensive GSE222520 analysis

Write-Host "`n" -NoNewline
Write-Host "="*70 -ForegroundColor Cyan
Write-Host "ALL 21 SAMPLES ANALYSIS - PROGRESS MONITOR" -ForegroundColor Cyan
Write-Host "="*70 -ForegroundColor Cyan
Write-Host ""

$outputDir = "results\all_samples"
$expectedFile = "all_samples_processed.h5ad"
$summaryFile = "sample_summary.csv"

# Check if output directory exists
if (-not (Test-Path $outputDir)) {
    Write-Host "[ERROR] Output directory not found: $outputDir" -ForegroundColor Red
    Write-Host "Analysis may not have started yet.`n" -ForegroundColor Yellow
    exit
}

Write-Host "[CHECKING] Output directory: $outputDir`n" -ForegroundColor Yellow

# Get all files in output directory
$files = Get-ChildItem -Path $outputDir -File -ErrorAction SilentlyContinue | Sort-Object LastWriteTime -Descending

if ($files.Count -eq 0) {
    Write-Host "[STATUS] No output files yet..." -ForegroundColor Yellow
    Write-Host "Analysis is likely in early stages (loading samples)`n" -ForegroundColor Yellow
} else {
    Write-Host "[FOUND] Output files ($($files.Count) total):`n" -ForegroundColor Green
    
    foreach ($file in $files) {
        $sizeMB = [math]::Round($file.Length / 1MB, 2)
        $time = $file.LastWriteTime.ToString("HH:mm:ss")
        $age = (Get-Date) - $file.LastWriteTime
        $ageMin = [math]::Round($age.TotalMinutes, 1)
        
        Write-Host "  - $($file.Name)" -ForegroundColor White -NoNewline
        Write-Host " ($sizeMB MB)" -ForegroundColor Gray -NoNewline
        Write-Host " [Modified: $time, $ageMin min ago]" -ForegroundColor DarkGray
    }
    Write-Host ""
}

# Check for expected outputs
Write-Host "[EXPECTED] Key output files:`n" -ForegroundColor Cyan

$status = @{}

# Check main output file
if (Test-Path (Join-Path $outputDir $expectedFile)) {
    $mainFile = Get-Item (Join-Path $outputDir $expectedFile)
    $sizeMB = [math]::Round($mainFile.Length / 1MB, 2)
    Write-Host "  [OK] " -ForegroundColor Green -NoNewline
    Write-Host "$expectedFile ($sizeMB MB)" -ForegroundColor White
    $status['main'] = $true
} else {
    Write-Host "  [ ] " -ForegroundColor Yellow -NoNewline
    Write-Host "$expectedFile (not found - still processing)" -ForegroundColor Gray
    $status['main'] = $false
}

# Check summary file
if (Test-Path (Join-Path $outputDir $summaryFile)) {
    Write-Host "  [OK] " -ForegroundColor Green -NoNewline
    Write-Host "$summaryFile" -ForegroundColor White
    $status['summary'] = $true
    
    # Read and display sample summary
    $csv = Import-Csv (Join-Path $outputDir $summaryFile)
    Write-Host "`n[SAMPLES] Loaded ($($csv.Count) samples):" -ForegroundColor Cyan
    foreach ($row in $csv) {
        $cells = [int]$row.n_cells
        Write-Host "  - $($row.sample) ($($row.group)): $($cells.ToString('N0')) cells" -ForegroundColor White
    }
} else {
    Write-Host "  [ ] " -ForegroundColor Yellow -NoNewline
    Write-Host "$summaryFile (not found)" -ForegroundColor Gray
    $status['summary'] = $false
}

Write-Host ""

# Overall status
Write-Host "="*70 -ForegroundColor Cyan

if ($status['main'] -and $status['summary']) {
    Write-Host "[SUCCESS] ANALYSIS COMPLETE!" -ForegroundColor Green
    Write-Host "="*70 -ForegroundColor Green
    
    # Get file info
    $mainFile = Get-Item (Join-Path $outputDir $expectedFile)
    $sizeMB = [math]::Round($mainFile.Length / 1MB, 2)
    
    Write-Host "`nResults ready:" -ForegroundColor Green
    Write-Host "  - File: $expectedFile" -ForegroundColor White
    Write-Host "  - Size: $sizeMB MB" -ForegroundColor White
    Write-Host "  - Location: $outputDir\" -ForegroundColor White
    
    Write-Host "`nNext steps:" -ForegroundColor Cyan
    Write-Host "  1. Generate Gupta-style figures:" -ForegroundColor White
    Write-Host "     py scripts\gupta_style_simple.py --input $outputDir\$expectedFile --output results\all_samples_figures" -ForegroundColor Gray
    Write-Host "`n  2. Analyze DGAT1 across all groups" -ForegroundColor White
    Write-Host "     py scripts\analyze_dgat1_all_samples.py --input $outputDir\$expectedFile" -ForegroundColor Gray
    Write-Host "`n  3. Compare Primary vs Recurrent" -ForegroundColor White
    Write-Host "     py scripts\compare_groups.py --input $outputDir\$expectedFile" -ForegroundColor Gray
    
} elseif ($status['summary']) {
    Write-Host "[RUNNING] Analysis in progress - samples loaded" -ForegroundColor Yellow
    Write-Host "="*70 -ForegroundColor Yellow
    
    Write-Host "`nCurrent stage: Processing cells (QC, normalization, clustering)" -ForegroundColor Yellow
    Write-Host "Estimated remaining time: 20-40 minutes" -ForegroundColor Yellow
    
    Write-Host "`nStages:" -ForegroundColor Cyan
    Write-Host "  [OK] Loading samples" -ForegroundColor Green
    Write-Host "  [>>] Quality control" -ForegroundColor Yellow
    Write-Host "  [ ] Normalization" -ForegroundColor Gray
    Write-Host "  [ ] Batch correction (Harmony)" -ForegroundColor Gray
    Write-Host "  [ ] Clustering & UMAP" -ForegroundColor Gray
    Write-Host "  [ ] Cell type annotation" -ForegroundColor Gray
    
} else {
    Write-Host "[RUNNING] Analysis in progress - early stage" -ForegroundColor Yellow
    Write-Host "="*70 -ForegroundColor Yellow
    
    Write-Host "`nCurrent stage: Loading and combining samples" -ForegroundColor Yellow
    Write-Host "Estimated remaining time: 30-60 minutes" -ForegroundColor Yellow
    
    Write-Host "`nStages:" -ForegroundColor Cyan
    Write-Host "  [>>] Loading samples (21 total)" -ForegroundColor Yellow
    Write-Host "  [ ] Quality control" -ForegroundColor Gray
    Write-Host "  [ ] Normalization" -ForegroundColor Gray
    Write-Host "  [ ] Batch correction (Harmony)" -ForegroundColor Gray
    Write-Host "  [ ] Clustering & UMAP" -ForegroundColor Gray
    Write-Host "  [ ] Cell type annotation" -ForegroundColor Gray
}

Write-Host ""

# Show timestamp
$now = Get-Date -Format "yyyy-MM-dd HH:mm:ss"
Write-Host "Last checked: $now" -ForegroundColor DarkGray
Write-Host "`nRun this script again to check progress:`n  .\check_all_samples_progress.ps1`n" -ForegroundColor Cyan

