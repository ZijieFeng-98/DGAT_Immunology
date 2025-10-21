# Quick Start Guide

**Last Updated:** 2025-10-10  
**Status:** ✅ Clean & Organized Structure

---

## 📂 Directory Structure

```
Scripts/
├── Transcriptomics/          RNA-seq Analysis
├── Proteomics/              Mass Spec Analysis
├── _Archive/                Old/obsolete scripts
├── GDCdata/                 Raw TCGA downloads
└── archive/                 Original archive
```

---

## 🚀 Quick Start

### For Transcriptomics (RNA-seq)

```bash
# 1. Read the overview
cat Scripts/Transcriptomics/README.md

# 2. Setup environment (first time only)
Rscript Scripts/Transcriptomics/00_setup/setup_environment.R

# 3. Download data (or use existing)
Rscript Scripts/Transcriptomics/01_download/01_download_tcga.R
Rscript Scripts/Transcriptomics/01_download/02_download_cgga.R

# 4. Preprocess
Rscript Scripts/Transcriptomics/01_preprocessing/01_tcga_fpkm_preprocessing.R
Rscript Scripts/Transcriptomics/01_preprocessing/02_cgga_fpkm_preprocessing.R

# 5. Batch correction
Rscript Scripts/Transcriptomics/02_batch_correction/01_tcga_batch_correction.R
Rscript Scripts/Transcriptomics/02_batch_correction/02_cgga_batch_correction.R

# 6. Analysis
Rscript Scripts/Transcriptomics/03_analysis/01_dgat_immune_analysis.R
Rscript Scripts/Transcriptomics/03_analysis/02_survival_analysis.R
```

### For Proteomics (Mass Spec)

```bash
# 1. Read the overview
cat Scripts/Proteomics/README.md

# 2. Process CPTAC data
Rscript Scripts/Proteomics/02_preprocessing/01_cptac_preprocessing.R

# 3. If needed, troubleshoot
Rscript Scripts/Proteomics/02_preprocessing/02_cptac_troubleshooting.R

# 4. Analysis
Rscript Scripts/Proteomics/03_analysis/01_protein_survival.R
```

---

## 📖 Documentation

| Location | Content |
|----------|---------|
| `Scripts/README.md` | Master overview of all scripts |
| `Scripts/REORGANIZATION_PLAN.md` | Details of reorganization |
| `Transcriptomics/README.md` | Complete RNA-seq pipeline guide |
| `Proteomics/README.md` | Complete proteomics pipeline guide |
| Each subdirectory | Specific section documentation |

---

## 🔍 Finding Scripts

### By Purpose:

| Need to... | Go to... |
|------------|----------|
| Setup R environment | `Transcriptomics/00_setup/` |
| Download data | `Transcriptomics/01_download/` |
| Preprocess RNA-seq | `Transcriptomics/01_preprocessing/` |
| Batch correction | `Transcriptomics/02_batch_correction/` |
| DGAT-immune analysis | `Transcriptomics/03_analysis/01_dgat_immune_analysis.R` |
| Survival analysis | `Transcriptomics/03_analysis/02_survival_analysis.R` |
| Process proteomics | `Proteomics/02_preprocessing/` |
| Protein survival | `Proteomics/03_analysis/` |

### By Data Type:

| Data Type | Directory |
|-----------|-----------|
| **TCGA RNA-seq** | `Transcriptomics/` (all scripts with `tcga` in name) |
| **CGGA RNA-seq** | `Transcriptomics/` (all scripts with `cgga` in name) |
| **CPTAC Proteomics** | `Proteomics/` (all scripts with `cptac` in name) |
| **HPA Protein Data** | `Proteomics/01_download/download_hpa.R` |

---

## 🎯 Common Tasks

### Validate Processed Data
```bash
Rscript Scripts/Transcriptomics/01_preprocessing/03_validate_clean_data.R
```

### Check DGAT1 in Proteomics
```r
protein <- readRDS("Processed_Data/CPTAC_GBM_Proteomics/protein_matrix_cleaned.rds")
"DGAT1" %in% rownames(protein)  # Should return TRUE
```

### Troubleshoot Missing Genes
```bash
Rscript Scripts/Proteomics/02_preprocessing/02_cptac_troubleshooting.R
```

---

## 📊 Data Flow

### Transcriptomics
```
Raw_Data/TCGA_GBM/ or CGGA_GBM/
    ↓ [01_download]
    ↓ [01_preprocessing]
Processed_Data/.../Master_Pipeline/
    ↓ [02_batch_correction]
Processed_Data/.../Batch_Corrected/
    ↓ [03_analysis]
Results/Bulk/
```

### Proteomics
```
Raw_Data/HPA_Protein/CPTAC3_...
    ↓ [02_preprocessing]
Processed_Data/CPTAC_GBM_Proteomics/
    ↓ [03_analysis]
Results/Proteomics/
```

---

## 🛠️ Troubleshooting

### Script Not Found?
- Check you're in the correct omics folder (Transcriptomics vs Proteomics)
- Scripts are numbered - follow the order: 00 → 01 → 02 → 03...

### Data Not Found?
- Make sure you've run previous steps in order
- Check `Processed_Data/` for intermediate outputs
- Verify input file paths in script headers

### DGAT1 Not Detected?
- For **RNA-seq:** Check gene symbol mapping
- For **Proteomics:** Run troubleshooting script

### Old Scripts?
- All old/duplicate scripts are in `_Archive/`
- Use NEW organized scripts in `Transcriptomics/` and `Proteomics/`

---

## 📝 Script Naming Convention

- **Numbered prefixes:** Show execution order (01, 02, 03...)
- **Dataset names:** `tcga`, `cgga`, `cptac`, `hpa`
- **Purpose:** `preprocessing`, `batch_correction`, `analysis`, etc.

**Example:** `01_tcga_fpkm_preprocessing.R`
- `01` = First in workflow
- `tcga` = TCGA dataset
- `fpkm` = Data type
- `preprocessing` = Purpose

---

## ✅ Current Status

| Component | Scripts | READMEs | Status |
|-----------|---------|---------|--------|
| **Transcriptomics** | 14 | 5 | ✅ Complete |
| **Proteomics** | 5 | 3 | ✅ Complete |
| **Documentation** | - | 12 | ✅ Complete |
| **Duplicates** | 0 | - | ✅ Removed |

---

## 🎓 Learning Path

**New to the project?** Follow this order:

1. Read `Scripts/README.md` - Overview
2. Read `Transcriptomics/README.md` - RNA-seq details
3. Read `Proteomics/README.md` - Proteomics details
4. Start with section-specific READMEs as needed

**Ready to analyze?** Follow the workflow:

1. **Transcriptomics:** 00_setup → 01_download → 01_preprocessing → 02_batch_correction → 03_analysis
2. **Proteomics:** 01_download → 02_preprocessing → 03_analysis

---

## 💡 Tips

✅ **Always read the README first** before running scripts  
✅ **Check input/output paths** in script headers  
✅ **Run scripts in order** (follow numbered directories)  
✅ **Validate outputs** using validation scripts  
✅ **Check QC plots** after processing  

---

**Questions?** Check the relevant README file or master `Scripts/README.md`

**Last Updated:** 2025-10-10  
**Structure:** Clean & Duplicate-Free ✅


