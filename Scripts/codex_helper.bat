@echo off
:: Codex Helper - Launch Codex for DGAT Immunology Project
:: Usage: codex_helper.bat [optional prompt]

echo ========================================
echo Codex Agent for DGAT Immunology Project
echo ========================================
echo.

if "%~1"=="" (
    echo Starting interactive Codex session...
    wsl -e bash -lic "cd /mnt/d/DGAT_Immunology && codex"
) else (
    echo Running Codex with prompt: %*
    wsl -e bash -lic "cd /mnt/d/DGAT_Immunology && codex '%*'"
)

