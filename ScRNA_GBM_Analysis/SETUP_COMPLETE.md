# Complete Automated Setup Guide

## 🎯 One-Time Setup (Do This First)

### Step 1: Install Python (Required - One Time Only)

**Download Python 3.9 or 3.10:**
1. Go to: https://www.python.org/downloads/
2. Download "Python 3.10.x" (recommended)
3. **IMPORTANT**: Check "Add Python to PATH" during installation
4. Complete installation

**Or use Anaconda (Alternative):**
1. Go to: https://www.anaconda.com/download
2. Download and install Anaconda
3. This includes Python + many scientific packages

---

## 🚀 Automated Setup (After Python is Installed)

### Run This One Command:

```powershell
cd D:\DGAT_Immunology\ScRNA_GBM_Analysis
.\setup_complete.ps1
```

This will automatically:
✅ Create virtual environment
✅ Install all Python packages
✅ Download Klemm et al. GBM dataset from GEO
✅ Organize data into correct folders
✅ Verify installation
✅ Run test to ensure everything works

---

## 📥 What Gets Downloaded

### Klemm et al., 2020 Dataset (GSE163108)
- **Source**: NCBI GEO (Gene Expression Omnibus)
- **Size**: ~500 MB compressed, ~2 GB uncompressed
- **Content**: Single-cell RNA-seq from GBM tumours + normal brain
- **Samples**: 40 GBM patients + 6 normal controls
- **Format**: Count matrices (will be converted to 10x format)

### Alternative: Test with Small Demo Dataset
If you want to test first with smaller data:
```powershell
.\download_demo_data.ps1
```
This downloads a small subset for testing (~50 MB)

---

## ⚙️ Manual Setup (If Automated Fails)

### 1. Create Virtual Environment
```powershell
cd D:\DGAT_Immunology\ScRNA_GBM_Analysis
python -m venv venv
.\venv\Scripts\Activate.ps1
```

### 2. Install Dependencies
```powershell
pip install --upgrade pip
pip install -r requirements.txt
```

### 3. Download Data Manually
See `DATA_DOWNLOAD_GUIDE.md` for detailed instructions

---

## ✅ Verify Installation

After setup completes, verify with:
```powershell
python scripts\verify_setup.py
```

This checks:
- Python version
- All packages installed
- Data downloaded and formatted correctly
- Scripts are executable

---

## 🎬 After Setup is Complete

Run the pipeline:
```powershell
python scripts\step_by_step.py
```

Or use the launcher:
```powershell
.\run_pipeline.ps1
```

---

## 🐛 Troubleshooting

### "Python not found"
- Reinstall Python with "Add to PATH" checked
- Or restart your terminal after installation

### "pip not found"
```powershell
python -m ensurepip --upgrade
```

### "Permission denied" errors
- Run PowerShell as Administrator
- Or adjust execution policy:
```powershell
Set-ExecutionPolicy -Scope CurrentUser RemoteSigned
```

### Download fails
- Check internet connection
- Try manual download from GEO website
- Use alternative mirror (see DATA_DOWNLOAD_GUIDE.md)

---

## 📊 Disk Space Requirements

- Python + packages: ~2 GB
- Raw data: ~2 GB
- Processed data: ~1 GB
- Results + figures: ~500 MB
- **Total: ~6 GB free space needed**

---

## ⏱️ Time Estimates

- Python installation: 5 minutes
- Package installation: 10-15 minutes
- Data download: 10-30 minutes (depends on internet speed)
- **Total setup time: ~30-50 minutes**

Once setup, analysis runtime: ~30-45 minutes

---

**Ready?** Install Python first, then run `.\setup_complete.ps1`

