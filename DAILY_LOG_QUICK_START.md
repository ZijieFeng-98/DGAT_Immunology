# 📝 Daily Logging - Quick Start

## ⚡ Super Quick Commands

```bash
# Log something you just did
./log_today.sh "Analyzed TCGA immune correlations"

# View today's log
./log_today.sh --view

# View this week
./log_today.sh --weekly

# Open log in editor
./log_today.sh
```

---

## 🎯 Most Common Use Cases

### 1. Log from Terminal (Easiest!)
```bash
./log_today.sh "Completed survival analysis for DGAT1"
./log_today.sh "Meeting with Dr. Guo about SREBP1 findings"
./log_today.sh "Fixed batch correction bug in preprocessing"
```

### 2. Log from R Script
Add to your R analysis scripts:

```r
# At the top
source("Scripts/log_activity.R")

# Log progress
log_data("TCGA_GBM", "Loaded 286 samples")
log_result("Found 12 significant immune correlations")
log_activity("Updated figures for manuscript")
```

### 3. Categories
| Type | Category |
|------|----------|
| Running analyses | `analysis` |
| Loading/processing data | `data` |
| Key findings | `result` |
| Meetings | `meeting` |
| Reading papers | `reading` |
| Writing | `writing` |

```r
log_activity("Read Smith 2024 lipid droplet paper", category = "reading")
log_activity("Updated Methods section", category = "writing")
```

---

## 📂 Where Are Logs Saved?

**Location**: `Notes/Daily_Logs/WorkLog_YYYY-MM-DD.md`

**Example**: `Notes/Daily_Logs/WorkLog_2025-10-20.md`

---

## 💡 Pro Tips

1. **Log as you go** - Don't wait until end of day
2. **Be specific** - "Ran TCGA analysis" → "Found DGAT1 correlates with M2 macrophages (r=0.45)"
3. **Use categories** - Makes logs easier to scan later
4. **Review weekly** - Use `./log_today.sh --weekly` on Fridays

---

## 🆘 Troubleshooting

**Can't run `./log_today.sh`?**
```bash
chmod +x log_today.sh
```

**Want to view past logs?**
```bash
ls Notes/Daily_Logs/
cat Notes/Daily_Logs/WorkLog_2025-10-15.md
```

---

## 📖 Full Documentation

See: `Notes/DAILY_LOGGING_SYSTEM.md` for complete guide

---

**That's it! Start logging your daily work now:**
```bash
./log_today.sh "Started using the daily logging system!"
```

