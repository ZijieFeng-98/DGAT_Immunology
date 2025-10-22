# Check progress of real data analysis
# Run this anytime to see what's been generated

Write-Host "`n" -NoNewline
Write-Host "="*70 -ForegroundColor Cyan
Write-Host "Real GBM Data Analysis - Progress Check" -ForegroundColor Green
Write-Host "="*70 -ForegroundColor Cyan

$resultsDir = "results\real_pilot"

if (Test-Path $resultsDir) {
    Write-Host "`n[DATA] Analysis Directory: $resultsDir" -ForegroundColor Yellow
    
    # Check H5AD files (checkpoints)
    $h5adFiles = Get-ChildItem $resultsDir -Filter "*.h5ad" -ErrorAction SilentlyContinue
    if ($h5adFiles) {
        Write-Host "`n[CHECKPOINT] H5AD Files Created:" -ForegroundColor Green
        foreach ($file in $h5adFiles) {
            $sizeMB = [math]::Round($file.Length / 1MB, 2)
            Write-Host "  - $($file.Name) ($sizeMB MB) - $($file.LastWriteTime)" -ForegroundColor White
        }
    }
    
    # Check figures
    $figures = Get-ChildItem "$resultsDir\figures" -Filter "*.png" -ErrorAction SilentlyContinue 2>$null
    if ($figures) {
        Write-Host "`n[PLOT] Figures Generated:" -ForegroundColor Green
        foreach ($fig in $figures) {
            Write-Host "  - $($fig.Name)" -ForegroundColor White
        }
    }
    
    # Check CSV files
    $csvFiles = Get-ChildItem $resultsDir -Filter "*.csv" -ErrorAction SilentlyContinue
    if ($csvFiles) {
        Write-Host "`n[FILE] CSV Files:" -ForegroundColor Green
        foreach ($csv in $csvFiles) {
            Write-Host "  - $($csv.Name)" -ForegroundColor White
        }
    }
    
    # Estimate progress based on files
    $totalExpected = 8  # 8 checkpoint files expected
    $currentCheckpoints = $h5adFiles.Count
    $progress = [math]::Round(($currentCheckpoints / $totalExpected) * 100, 0)
    
    Write-Host "`n[PROGRESS] Estimated: $progress% complete ($currentCheckpoints/$totalExpected steps)" -ForegroundColor Cyan
    
    if ($currentCheckpoints -ge 8) {
        Write-Host "`n[OK] ANALYSIS COMPLETE!" -ForegroundColor Green
        Write-Host "Run validation: py scripts\validate_and_summarize.py" -ForegroundColor Yellow
    } elseif ($currentCheckpoints -gt 0) {
        Write-Host "`n[INFO] Analysis in progress..." -ForegroundColor Yellow
        Write-Host "Current step: ~Step $currentCheckpoints of 8" -ForegroundColor White
    } else {
        Write-Host "`n[INFO] Analysis starting or not yet begun" -ForegroundColor Yellow
    }
    
} else {
    Write-Host "`n[INFO] Results directory not created yet" -ForegroundColor Yellow
    Write-Host "Analysis may be starting or data loading in progress..." -ForegroundColor White
}

Write-Host "`n" -NoNewline
Write-Host "="*70 -ForegroundColor Cyan
Write-Host "Tip: Run this script again to check updated progress" -ForegroundColor Gray
Write-Host "="*70 -ForegroundColor Cyan
Write-Host ""

