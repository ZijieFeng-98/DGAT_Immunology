# Fix Python PATH Issue

## ⚠️ Problem: Python is installed but not in PATH

You have Python installed, but it's not accessible from the command line.

---

## 🔧 Quick Fixes (Choose One)

### **Option 1: Find Your Python Installation**

Run this in PowerShell to find Python:

```powershell
# Search for Python executable
Get-ChildItem -Path C:\ -Filter python.exe -Recurse -ErrorAction SilentlyContinue | Select-Object FullName
```

Once you find it (e.g., `C:\Users\YourName\AppData\Local\Programs\Python\Python310\python.exe`), use the full path:

```powershell
# Test with full path
C:\Users\YourName\AppData\Local\Programs\Python\Python310\python.exe --version

# Then run setup with full path
cd D:\DGAT_Immunology\ScRNA_GBM_Analysis
C:\Users\YourName\AppData\Local\Programs\Python\Python310\python.exe -m venv venv
.\venv\Scripts\Activate
pip install -r requirements.txt
```

---

### **Option 2: Add Python to PATH Manually**

1. Press `Win + X` → Choose "System"
2. Click "Advanced system settings"
3. Click "Environment Variables"
4. Under "User variables", select "Path" → Click "Edit"
5. Click "New" and add your Python path (e.g., `C:\Python310\`)
6. Click "New" again and add Scripts folder (e.g., `C:\Python310\Scripts\`)
7. Click "OK" on all windows
8. **Restart PowerShell**
9. Test: `python --version`

---

### **Option 3: Reinstall Python (Recommended)**

**This is the easiest solution:**

1. Download Python 3.10 from: https://www.python.org/downloads/
2. Run installer
3. ✅ **CHECK "Add Python to PATH"** (most important!)
4. Choose "Install Now"
5. Wait for installation
6. **Restart PowerShell**
7. Test: `python --version`

---

### **Option 4: Use Microsoft Store Python (Alternative)**

1. Open Microsoft Store
2. Search "Python 3.10"
3. Click "Get"
4. Wait for installation
5. **Restart PowerShell**
6. Test: `python --version`

---

### **Option 5: Use Anaconda (Alternative)**

If you prefer Anaconda:

1. Download from: https://www.anaconda.com/download
2. Install Anaconda
3. Open "Anaconda Prompt" (not regular PowerShell)
4. Navigate to project:
   ```bash
   cd D:\DGAT_Immunology\ScRNA_GBM_Analysis
   ```
5. Create environment:
   ```bash
   conda create -n scRNA python=3.10
   conda activate scRNA
   pip install -r requirements.txt
   ```

---

## ✅ After Python is Working

Once `python --version` works, run:

```powershell
cd D:\DGAT_Immunology\ScRNA_GBM_Analysis
.\setup_complete.ps1
```

This will automatically set up everything!

---

## 🆘 Still Having Issues?

### Manual Setup (if automated script fails):

```powershell
# Navigate to project
cd D:\DGAT_Immunology\ScRNA_GBM_Analysis

# Create virtual environment
python -m venv venv

# Activate it
.\venv\Scripts\Activate.ps1

# Upgrade pip
python -m pip install --upgrade pip

# Install packages
pip install scanpy numpy pandas matplotlib seaborn scipy
pip install scrublet harmonypy infercnvpy anndata
pip install tqdm jupyter

# Create demo data
python scripts\create_demo_data.py

# Verify setup
python scripts\verify_setup.py

# Test pipeline
python scripts\sc_pipeline.py --tumour_path data\raw\demo_tumour\ --normal_path data\raw\demo_normal\ --output_dir results\demo\ --max_genes 5000
```

---

## 📞 Quick Diagnosis

Run this to diagnose:

```powershell
# Check Python
python --version
python3 --version
py --version

# Check conda
conda --version

# Search for Python
Get-Command python -ErrorAction SilentlyContinue
where python
```

---

**Recommended**: Just reinstall Python with "Add to PATH" checked. Fastest solution! ✨

