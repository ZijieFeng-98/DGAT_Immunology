# Codex Agent - Quick Reference

## Launch Codex

### Method 1: Project Helper (Recommended)
```powershell
cd D:\DGAT_Immunology
.\Scripts\codex_helper.ps1
```

### Method 2: Direct Command
```powershell
codex
```

### Method 3: From Anywhere
```powershell
D:\codex
```

## Common Tasks

### 📊 Data Analysis
```bash
codex "analyze the TCGA_GBM survival data"
codex "help me understand the batch correction results"
codex "create QC plots for CPTAC proteomics"
```

### 🔬 R Script Help
```bash
codex "review Scripts/Transcriptomics/03_immune_marker_correlation.R"
codex "optimize the deconvolution analysis script"
codex "add error handling to data processing pipeline"
```

### 📝 Documentation
```bash
codex "document the master pipeline workflow"
codex "create a README for Processed_Data/TCGA_GBM/"
codex "explain what project_config.R does"
```

### 🔍 Code Review
```bash
codex "review all transcriptomics scripts for consistency"
codex "check for deprecated R functions in Scripts/"
codex "suggest improvements for proteomics pipeline"
```

### 🛠️ Troubleshooting
```bash
codex "why is my PCA plot not showing batch effects?"
codex "help debug the survival analysis script"
codex "fix errors in the lipid marker correlation"
```

## In-Codex Commands

While Codex is running:

- `/help` - Show available commands
- `/approvals` - Change approval mode
- `/exit` or Ctrl+C - Exit Codex
- `/clear` - Clear conversation history
- `/undo` - Undo last action

## Approval Modes

- **Auto**: ✅ Read/edit in working dir | ❓ Ask for network/outside access
- **Read Only**: ✅ Inspect/explain | ❌ No changes
- **Full Access**: ✅ Everything without approval | ⚠️ Use carefully

## Best Practices

1. **Be specific**: "Review X script for Y issue" vs "help me"
2. **One task at a time**: Break complex requests into steps
3. **Review changes**: Always check what Codex modified
4. **Use context**: Reference specific files/directories
5. **Iterate**: Refine requests based on Codex responses

## File Locations

- **Codex Launcher**: `D:\codex.bat`
- **Helper Scripts**: `D:\DGAT_Immunology\Scripts\codex_helper.*`
- **Config Guide**: `D:\DGAT_Immunology\.codex_config.md`
- **This File**: `D:\DGAT_Immunology\Scripts\CODEX_QUICK_REFERENCE.md`

## Example Workflow

```powershell
# 1. Navigate to project
cd D:\DGAT_Immunology

# 2. Launch Codex
codex

# 3. In Codex prompt:
> Review Scripts/Transcriptomics/ for code quality

# 4. Follow up:
> Optimize the most inefficient script you found

# 5. Exit when done:
> /exit
```

## Authentication

First time running Codex:
1. You'll be prompted to authenticate
2. **Option A**: Sign in with ChatGPT account (Plus/Pro/Team/Enterprise)
3. **Option B**: Use OpenAI API key

## Getting Help

- Codex docs: https://developers.openai.com/codex/cli/
- This project: See `.codex_config.md`
- Issues: Check Codex output for error messages

