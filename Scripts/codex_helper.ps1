# Codex Helper for DGAT Immunology Project
# Launch Codex AI agent for code assistance

param(
    [Parameter(ValueFromRemainingArguments=$true)]
    [string[]]$Prompt
)

Write-Host "========================================" -ForegroundColor Cyan
Write-Host "Codex Agent for DGAT Immunology Project" -ForegroundColor Cyan
Write-Host "========================================" -ForegroundColor Cyan
Write-Host ""

# Ensure we're in the project directory
$projectDir = "D:\DGAT_Immunology"
Set-Location $projectDir

if ($Prompt) {
    $promptText = $Prompt -join " "
    Write-Host "Running Codex with prompt: $promptText" -ForegroundColor Yellow
    wsl -e bash -lic "cd /mnt/d/DGAT_Immunology && codex '$promptText'"
} else {
    Write-Host "Starting interactive Codex session..." -ForegroundColor Green
    Write-Host "Working directory: $projectDir" -ForegroundColor Gray
    Write-Host ""
    wsl -e bash -lic "cd /mnt/d/DGAT_Immunology && codex"
}

