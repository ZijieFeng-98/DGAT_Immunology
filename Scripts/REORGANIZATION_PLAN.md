# Scripts Reorganization Plan

**Date:** 2025-10-10  
**Purpose:** Separate Transcriptomics and Proteomics with clear workflow structure

---

## Current Structure Analysis

### ✅ Keep (Already organized)
- `Transcriptomics/` - Newly created, being populated
- `Proteomics/` - Newly created, being populated

### 📦 Archive (Obsolete/Documented)
- `00_organize/` → Archive (obsolete, analysis documented)
- `Analysis_Framework/` → Archive (old framework, superseded)
- `archive/` → Keep as is

### 🔄 Reorganize into Transcriptomics

| Current Location | Purpose | New Location |
|-----------------|---------|--------------|
| `01_setup/00_setup_env_simple.R` | R environment setup | `Transcriptomics/00_setup/setup_environment.R` |
| `02_data_download/01_get_tcga.R` | TCGA download | `Transcriptomics/01_download/01_download_tcga.R` |
| `02_data_download/02_get_cgga.R` | CGGA download | `Transcriptomics/01_download/02_download_cgga.R` |
| `02_data_download/03_get_gtex.R` | GTEx download | `Transcriptomics/01_download/03_download_gtex.R` |
| `03_data_processing/gbm_fpkm_simple.R` | TCGA preprocessing | `Transcriptomics/02_preprocessing/01_tcga_fpkm_preprocessing.R` ✅ |
| `03_data_processing/cgga_fpkm_simple.R` | CGGA preprocessing | `Transcriptomics/02_preprocessing/02_cgga_fpkm_preprocessing.R` ✅ |
| `03_data_processing/01_validate_clean_data.R` | Validation | `Transcriptomics/02_preprocessing/03_validate_clean_data.R` ✅ |
| `03_data_processing/batch_correction_customized.R` | TCGA batch correction | `Transcriptomics/03_batch_correction/01_tcga_batch_correction.R` ✅ |
| `03_data_processing/cgga_batch_correction_customized.R` | CGGA batch correction | `Transcriptomics/03_batch_correction/02_cgga_batch_correction.R` ✅ |
| `04_analysis/02_bulk_immune_analysis.R` | DGAT-immune analysis | `Transcriptomics/04_analysis/01_dgat_immune_analysis.R` ✅ |
| `04_analysis/01_bulk_survival_bestcut_auto.R` | Survival analysis | `Transcriptomics/04_analysis/02_survival_analysis.R` ✅ |
| `05_visualization/theme_heatmap.R` | Plotting theme | `Transcriptomics/05_visualization/theme_heatmap.R` |
| `06_utils/*.R` | Helper functions | `Transcriptomics/utils/` |

### 🔄 Reorganize into Proteomics

| Current Location | Purpose | New Location |
|-----------------|---------|--------------|
| `05_protein/01_get_hpa_dgat.R` | HPA query | `Proteomics/01_download/download_hpa.R` |
| `05_protein/01_cptac_gbm_proteomics_processing.R` | CPTAC processing | `Proteomics/02_preprocessing/01_cptac_preprocessing.R` ✅ |
| `05_protein/02_cptac_troubleshooting.R` | Troubleshooting | `Proteomics/02_preprocessing/02_cptac_troubleshooting.R` ✅ |
| `05_protein/03_protein_survival_analysis.R` | Protein survival | `Proteomics/03_analysis/01_protein_survival.R` |

### 🗑️ Delete (Redundant/Obsolete)

| File | Reason | Status |
|------|--------|--------|
| `03_data_processing/batch_correction_pipeline.R` | Superseded by customized version | ❌ Delete |
| `03_data_processing/batch_correction_simple.R` | Analysis only, not correction | ❌ Delete |
| `03_data_processing/gbm_fpkm_prep_master.R` | Obsolete, superseded by simple | ❌ Delete |
| `03_data_processing/02_dgat_basic_analysis.R` | Superseded by immune analysis | ❌ Delete |
| `02_data_download/00_run_all_downloads.R` | Wrapper script, can recreate | ❌ Delete |
| `02_data_download/00_run_remaining_downloads.R` | Wrapper script | ❌ Delete |
| `02_data_download/04_download_missing_datasets.R` | Ad-hoc script | ❌ Delete |

### 📂 Reference Only (Don't Move)
- `GDCdata/` - Contains raw downloaded TCGA data files

---

## New Structure

```
Scripts/
├── Transcriptomics/
│   ├── 00_setup/
│   │   ├── setup_environment.R
│   │   └── README.md
│   ├── 01_download/
│   │   ├── 01_download_tcga.R
│   │   ├── 02_download_cgga.R
│   │   ├── 03_download_gtex.R
│   │   └── README.md
│   ├── 02_preprocessing/
│   │   ├── 01_tcga_fpkm_preprocessing.R ✅
│   │   ├── 02_cgga_fpkm_preprocessing.R ✅
│   │   ├── 03_validate_clean_data.R ✅
│   │   └── README.md
│   ├── 03_batch_correction/
│   │   ├── 01_tcga_batch_correction.R ✅
│   │   ├── 02_cgga_batch_correction.R ✅
│   │   └── README.md
│   ├── 04_analysis/
│   │   ├── 01_dgat_immune_analysis.R ✅
│   │   ├── 02_survival_analysis.R ✅
│   │   └── README.md
│   ├── 05_visualization/
│   │   ├── theme_heatmap.R
│   │   └── README.md
│   ├── utils/
│   │   ├── utils_io.R
│   │   ├── utils_signatures.R
│   │   └── README.md
│   └── README.md ✅
│
├── Proteomics/
│   ├── 01_download/
│   │   ├── download_hpa.R
│   │   └── README.md
│   ├── 02_preprocessing/
│   │   ├── 01_cptac_preprocessing.R ✅
│   │   ├── 02_cptac_troubleshooting.R ✅
│   │   └── README.md
│   ├── 03_analysis/
│   │   ├── 01_protein_survival.R
│   │   └── README.md
│   └── README.md ✅
│
├── _Archive/
│   ├── 00_organize/ (moved)
│   ├── Analysis_Framework/ (moved)
│   └── old_scripts/ (moved)
│
└── REORGANIZATION_PLAN.md (this file)
```

---

## Implementation Steps

1. ✅ Create new directory structure
2. ⏳ Copy/move scripts to new locations
3. ⏳ Update README files for each section
4. ⏳ Create wrapper scripts for common workflows
5. ⏳ Archive obsolete directories
6. ⏳ Delete redundant scripts
7. ⏳ Update main project README

---

## Benefits

✅ **Clear separation** - Transcriptomics vs Proteomics  
✅ **Logical workflow** - Setup → Download → Process → Analyze  
✅ **Easy navigation** - Numbered directories show order  
✅ **Reduced clutter** - Obsolete scripts archived  
✅ **Better documentation** - README in each section  
✅ **Maintainable** - Clear structure for future additions

---

**Status:** In Progress  
**Next:** Execute reorganization


