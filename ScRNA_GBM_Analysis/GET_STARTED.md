# 🚀 GET STARTED - ScRNA_GBM_Analysis

## ✅ Project Setup Complete!

Your single-cell GBM analysis project is ready to go!

---

## 📋 What We've Built

### ✓ Complete Pipeline Implementation
- **Main pipeline** (`scripts/sc_pipeline.py`) - 443 lines, fully documented
- **Interactive runner** (`scripts/step_by_step.py`) - 8 modular steps with checkpoints
- **Configuration** (`config.yaml`) - Customizable parameters
- **Documentation** - README, Quick Start, Protocol Guide

### ✓ Analysis Workflow
1. **Data Loading** - 10x Genomics format
2. **Quality Control** - Filter low-quality cells
3. **Doublet Detection** - Scrublet algorithm
4. **Normalization** - Log-transform, HVG selection
5. **Batch Correction** - Harmony integration
6. **Clustering** - Leiden algorithm, UMAP
7. **CNV Inference** - Identify malignant cells
8. **Cell Annotation** - Marker-based typing
9. **DGAT1 Analysis** - Lipid metabolism profiling

### ✓ Project Structure
```
ScRNA_GBM_Analysis/
├── 📄 Documentation (README, QUICKSTART, etc.)
├── 📜 Scripts (pipeline + interactive)
├── 📁 data/ (raw & processed)
├── 📊 results/ (outputs)
├── 📈 figures/ (plots)
├── 📓 notebooks/ (Jupyter)
└── 📝 logs/
```

---

## 🎯 Quick Start (3 Options)

### Option 1: Windows Launcher (Easiest)
```powershell
cd D:\DGAT_Immunology\ScRNA_GBM_Analysis
.\run_pipeline.ps1
```
This will guide you through setup and execution.

### Option 2: Interactive Python
```powershell
cd D:\DGAT_Immunology\ScRNA_GBM_Analysis
python scripts\step_by_step.py
```
Run each step manually, inspect results, adjust parameters.

### Option 3: Direct Pipeline
```powershell
cd D:\DGAT_Immunology\ScRNA_GBM_Analysis
python scripts\sc_pipeline.py --tumour_path data\raw\tumour\ --normal_path data\raw\normal\ --output_dir results\
```
Fully automated with default parameters.

---

## 📥 Before You Start

### 1. Install Python Dependencies
```powershell
pip install -r requirements.txt
```

**Required packages:**
- scanpy (single-cell analysis)
- scrublet (doublet detection)
- harmonypy (batch correction)
- infercnvpy (CNV inference)
- matplotlib, seaborn (visualization)
- numpy, pandas, scipy

### 2. Prepare Your Data

Your data should be in **10x Genomics format**:

```
data/raw/
├── tumour/
│   ├── matrix.mtx.gz
│   ├── features.tsv.gz
│   └── barcodes.tsv.gz
└── normal/
    ├── matrix.mtx.gz
    ├── features.tsv.gz
    └── barcodes.tsv.gz
```

#### Don't have the Klemm data yet?

**Option A: Download from GEO**
- Search for **GSE163108** (Klemm et al., 2020) on NCBI GEO
- Download processed count matrices
- Convert to 10x format if needed

**Option B: Use your own GBM scRNA-seq data**
- Must be in 10x format
- Must have tumour + normal samples
- Ensure gene symbols are annotated

**Option C: Test with existing data**
- Use the Darmanis dataset already in your project:
  ```
  D:\DGAT_Immunology\Raw_Data\scRNA\Glioma\Darmanis_GSE84465_scRNA.rds
  ```
  (You'll need to convert from RDS to 10x format first)

---

## 🎬 Step-by-Step First Run

### Step 1: Navigate to Project
```powershell
cd D:\DGAT_Immunology\ScRNA_GBM_Analysis
```

### Step 2: Install Dependencies
```powershell
# Create virtual environment (recommended)
python -m venv venv
.\venv\Scripts\Activate

# Install packages
pip install -r requirements.txt
```

### Step 3: Add Your Data
Place 10x data in:
- `data/raw/tumour/`
- `data/raw/normal/`

### Step 4: Run Interactive Mode
```powershell
python scripts\step_by_step.py
```

### Step 5: Follow the Menu
```
SINGLE-CELL GBM ANALYSIS - INTERACTIVE MENU
============================================
Current step: 0/8

Available steps:
  1. Load data
  2. Quality control
  3. Doublet detection
  4. Normalization & integration
  5. Clustering
  6. CNV inference
  7. Cell-type annotation
  8. DGAT1 analysis

Other options:
  9. Run all remaining steps
  s. Show current statistics
  l. Load checkpoint
  q. Quit
```

### Step 6: Check Results
After each step, check:
- `results/checkpoints/` - Saved progress
- `results/qc_plots/` - Quality control plots
- Console output - Statistics and metrics

---

## 📊 Expected Timeline

| Step | Duration | Memory | Notes |
|------|----------|---------|-------|
| 1. Load data | 1-2 min | ~2GB | Fast |
| 2. QC filtering | 2-3 min | ~2GB | Fast |
| 3. Doublets | 3-5 min | ~3GB | Moderate |
| 4. Normalize & integrate | 5-10 min | ~4GB | Moderate |
| 5. Clustering | 2-3 min | ~3GB | Fast |
| 6. CNV inference | 10-15 min | ~6GB | **SLOW** ⚠️ |
| 7. Annotation | 1-2 min | ~2GB | Fast |
| 8. DGAT1 analysis | 2-3 min | ~2GB | Fast |
| **TOTAL** | **~30-45 min** | **Peak ~6GB** | |

*Times are for ~10,000 cells. Scale accordingly.*

---

## 🎨 What You'll Get

### Output Files
- `processed_adata.h5ad` - Final annotated dataset
- `DGAT1_expression_summary.csv` - Expression stats
- `cell_counts_summary.csv` - Cell type composition

### Figures
- UMAP colored by sample, clusters, cell types
- QC violin plots
- DGAT1 expression plots
- CNV scores
- Marker gene dot plots

### Checkpoints (for resuming)
- `01_loaded.h5ad`
- `02_qc_filtered.h5ad`
- `03_doublets_removed.h5ad`
- `04_normalized_integrated.h5ad`
- `05_clustered.h5ad`
- `06_cnv_inferred.h5ad`
- `07_annotated.h5ad`
- `08_final.h5ad`

---

## 🔧 Customization

### Adjust QC Thresholds
Edit `config.yaml`:
```yaml
qc:
  min_genes: 200    # Lower = keep more cells
  max_genes: 2500   # Higher = keep doublets
  max_mito: 10.0    # Lower = stricter QC
```

### Change Clustering Resolution
```yaml
clustering:
  resolution: 0.5   # 0.1-2.0 (higher = more clusters)
```

### Add Custom Cell-Type Markers
```yaml
markers:
  My_Custom_Cell_Type:
    - MARKER1
    - MARKER2
    - MARKER3
```

---

## 💡 Pro Tips

✅ **Start with interactive mode** - Best for learning and troubleshooting

✅ **Check QC plots** - Look for bimodal distributions, outliers

✅ **Save checkpoints** - Don't re-run everything if step 7 fails

✅ **Start with defaults** - Optimize parameters later

✅ **Monitor memory** - Close other programs during CNV inference

---

## 🐛 Troubleshooting

### "Module not found"
```powershell
pip install -r requirements.txt
```

### "Data not found"
- Check paths in `config.yaml`
- Ensure files are named correctly (matrix.mtx.gz, etc.)

### "Out of memory"
- Close other programs
- Reduce `n_top_genes` to 2000
- Process samples separately

### "Too many/few clusters"
- Adjust `resolution` parameter (0.1-2.0)
- Check in `config.yaml` or pass to script

### Pipeline is slow
- Normal! CNV inference takes 10-15 min
- Use checkpoints to avoid re-running
- Consider skipping CNV if not needed

---

## 📚 Documentation

- **QUICKSTART.md** - Fast getting started guide
- **README.md** - Complete documentation
- **config.yaml** - All parameters explained
- **PROJECT_STATUS.md** - Current project state
- **../Protocols/scRNA_Klemm_Pipeline_Guide.md** - Detailed methodology

---

## 🎓 Next Steps After First Run

1. **Explore results in Jupyter**
   ```powershell
   jupyter notebook notebooks/
   ```

2. **Customize analysis**
   - Add more marker genes
   - Try different clustering resolutions
   - Subset cell types of interest

3. **Advanced analyses**
   - Cell-cell communication (CellChat)
   - Trajectory inference (Monocle 3)
   - RNA velocity (scVelo)
   - Metabolic flux (METAFlux)

4. **Integrate with bulk data**
   - Compare with TCGA GBM in main project
   - Validate findings across platforms

5. **Publication**
   - Follow FAIR principles
   - Share code and data
   - Document all parameters

---

## 🆘 Need Help?

1. Check the error message carefully
2. Review documentation files
3. Inspect QC plots
4. Try with default parameters first
5. Check Scanpy documentation: https://scanpy.readthedocs.io/
6. Review sc-best-practices: https://www.sc-best-practices.org/

---

## ✨ Ready to Go!

Your single-cell GBM analysis pipeline is fully set up and ready to run!

**Start now:**
```powershell
cd D:\DGAT_Immunology\ScRNA_GBM_Analysis
.\run_pipeline.ps1
```

**Or go step-by-step:**
```powershell
python scripts\step_by_step.py
```

Good luck with your analysis! 🚀

---

**Created**: October 21, 2025  
**Location**: `D:\DGAT_Immunology\ScRNA_GBM_Analysis`  
**Status**: ✅ Ready to Run  
**Parent Project**: DGAT_Immunology

