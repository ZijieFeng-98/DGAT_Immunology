# 📝 Automated Daily Work Logging System

A simple, automated system to track your daily research activities for the DGAT Immunology project.

---

## 🚀 Quick Start

### Option 1: Log from Command Line (Fastest)
```bash
# Log a quick note
./log_today.sh "Analyzed TCGA immune correlations"

# View today's log
./log_today.sh --view

# View this week's summary
./log_today.sh --weekly

# Open today's log in editor
./log_today.sh
```

### Option 2: Log from R Scripts (Automatic)
Add this to the top of any R analysis script:

```r
# At the start of your script
source("Scripts/log_activity.R")
log_script_start("01_dgat_immune_analysis.R", "Running TCGA immune correlation analysis")

# During your analysis
log_data("TCGA_GBM", "Loaded expression data", "286 samples")
log_result("Found 12 significant immune correlations (FDR < 0.05)")

# At the end
log_script_end("01_dgat_immune_analysis.R", "Analysis complete, saved to Results/")
```

### Option 3: Auto-logging (Session-wide)
Add this to your R session startup:

```r
source("Scripts/log_activity.R")
auto_track()  # Automatically logs data loading operations
```

---

## 📊 Log Categories

The system uses different categories with emojis for easy visual scanning:

| Category | Emoji | Usage |
|----------|-------|-------|
| `script` | 🔧 | Running R/Python scripts |
| `analysis` | 📊 | Statistical analyses |
| `data` | 💾 | Data loading/processing |
| `result` | ✅ | Key findings/outputs |
| `error` | ❌ | Errors or issues |
| `meeting` | 👥 | Lab meetings, discussions |
| `reading` | 📖 | Papers, documentation |
| `writing` | ✍️ | Manuscript, reports |
| `planning` | 📋 | Planning next steps |

---

## 📁 Log Storage

All logs are saved to: `Notes/Daily_Logs/`

Format: `WorkLog_YYYY-MM-DD.md`

Example: `WorkLog_2025-10-20.md`

---

## 💡 Example Workflows

### Daily Research Session

```bash
# Morning: Start your day
./log_today.sh "Starting work on DGAT-immune correlations"

# Run your analyses (they auto-log if you source log_activity.R)
Rscript Scripts/Transcriptomics/03_analysis/01_dgat_immune_analysis.R

# Quick notes throughout the day
./log_today.sh "Meeting with Dr. Guo - discussed SREBP1 findings"
./log_today.sh "Read Smith et al. 2024 paper on lipid droplets in GBM"

# Evening: View your day
./log_today.sh --view
```

### Weekly Review

```bash
# See what you accomplished this week
./log_today.sh --weekly
```

---

## 🔧 Advanced Features

### Manual Logging Functions

```r
source("Scripts/log_activity.R")

# Log any activity
log_activity("Fixed batch correction pipeline", category = "analysis")

# Log data operations
log_data("CGGA_GBM", "Preprocessed expression data", "694 samples")

# Log results
log_result("DGAT1 high group shows increased M2 macrophages (p < 0.001)")

# Log errors/issues
log_error("GSVA failed on immune_checkpoints geneset - too few genes")

# Add end-of-day summary
add_daily_summary(
  accomplishments = "
  - Completed TCGA immune correlation analysis
  - Found significant M2 macrophage association
  - Generated all publication figures
  ",
  next_steps = "
  - Validate findings in CGGA cohort
  - Add partial correlation analysis adjusting for purity
  - Compare with proteomics data
  ",
  notes = "
  - SREBP1 correlation is very strong (r=0.68)
  - May need to check for outliers in purity < 0.3 samples
  "
)
```

### View Past Logs

```r
# View specific date
view_log(as.Date("2025-10-15"))

# Weekly summary
weekly_summary()
```

---

## ⚙️ Automatic Integration

### Add to Your .Rprofile (Optional)

Add this to automatically load logging in every R session:

```r
# Add to ~/.Rprofile or project .Rprofile
if (file.exists("Scripts/log_activity.R")) {
  source("Scripts/log_activity.R")
  log_activity("Started R session", category = "general")
}
```

### Shell Aliases (Optional)

Add to your `~/.zshrc` or `~/.bashrc`:

```bash
# Quick logging aliases
alias logwork='cd "/Users/zijiefeng/Desktop/Guo'\''s lab/My_Research/DGAT_Immunology" && ./log_today.sh'
alias viewlog='cd "/Users/zijiefeng/Desktop/Guo'\''s lab/My_Research/DGAT_Immunology" && ./log_today.sh --view'
alias weeklog='cd "/Users/zijiefeng/Desktop/Guo'\''s lab/My_Research/DGAT_Immunology" && ./log_today.sh --weekly'
```

Then you can use from anywhere:
```bash
logwork "Analyzed proteomics data"
viewlog
weeklog
```

---

## 📋 Example Log Output

```markdown
# Research Work Log - 2025-10-20

**Project:** DGAT Immunology Analysis  
**Researcher:** Zijie Feng  
**Date:** October 20, 2025  
**Day:** Monday

---

## Today's Activities

**[09:15:30]** • Started work session

**[09:20:15]** 🔧 Started running 01_dgat_immune_analysis.R - Running TCGA immune correlation analysis  
  _Script: `Scripts/Transcriptomics/03_analysis/01_dgat_immune_analysis.R`_

**[09:21:03]** 💾 TCGA_GBM: Loaded expression data (286 samples, 20K genes)

**[09:45:22]** ✅ Found 12 significant immune correlations (FDR < 0.05)

**[09:46:10]** ✅ Completed 01_dgat_immune_analysis.R - Analysis complete

**[11:00:00]** 👥 Lab meeting with Dr. Guo - discussed SREBP1 findings

**[14:30:00]** 📖 Read Smith et al. 2024 paper on lipid droplets in GBM

**[16:00:00]** ✍️ Updated manuscript methods section with batch correction details

---

## Daily Summary

### ✅ Accomplishments
- Completed TCGA immune correlation analysis
- Found significant M2 macrophage association with DGAT1
- Generated all publication figures

### 📋 Next Steps
- Validate findings in CGGA cohort
- Add partial correlation analysis adjusting for purity
- Compare with proteomics data

### 📝 Notes
- SREBP1 correlation is very strong (r=0.68)
- May need to check for outliers in purity < 0.3 samples

---
*Log generated: 2025-10-20 17:30:45*
```

---

## 🎯 Benefits

1. **📅 Track Progress**: See exactly what you did each day
2. **🔍 Reproducibility**: Know which scripts you ran and when
3. **📊 Lab Updates**: Quickly generate weekly summaries for meetings
4. **✍️ Manuscript Writing**: Easy to recall what analyses were done
5. **🧠 Memory Aid**: Never forget important findings or next steps

---

## 🆘 Troubleshooting

**Problem**: Can't run `./log_today.sh`  
**Solution**: Make it executable: `chmod +x log_today.sh`

**Problem**: Logs not saving  
**Solution**: Check that `Notes/Daily_Logs/` directory exists (it should auto-create)

**Problem**: Can't find log functions in R  
**Solution**: Make sure to `source("Scripts/log_activity.R")` first

---

## 📧 Questions?

This system is simple and customizable. Feel free to modify the logging functions or categories to fit your workflow!

**Location of main script**: `Scripts/log_activity.R`  
**Location of logs**: `Notes/Daily_Logs/`  
**Shell wrapper**: `log_today.sh`

