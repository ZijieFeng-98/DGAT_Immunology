# Transcriptomics Analysis Pipeline

**Data Type:** RNA-seq / Gene Expression (FPKM, TPM, Raw Counts)  
**Datasets:** TCGA-GBM, CGGA-GBM  
**Purpose:** Process and analyze gene expression data for DGAT-immune relationships

---

## Directory Structure

```
Transcriptomics/
├── 01_preprocessing/          # Data loading, QC, and preprocessing
├── 02_batch_correction/       # Batch effect removal
└── 03_analysis/              # Downstream analysis (immune, survival)
```

---

## 01_preprocessing/

**Purpose:** Load, clean, and prepare RNA-seq data for analysis

### Scripts:

1. **`01_tcga_fpkm_preprocessing.R`**
   - **Input:** `Raw_Data/TCGA_GBM/TCGA_GBM_Expression_FPKM_with_symbols.rds`
   - **Output:** `Processed_Data/TCGA_GBM_Master_Pipeline/`
   - **Functions:**
     - Load FPKM expression data
     - Filter primary tumors only
     - Deduplicate patients
     - Filter genes (ribosomal/mitochondrial, low detection, low variance)
     - Generate PCA plots
   - **Runtime:** ~2-3 minutes

2. **`02_cgga_fpkm_preprocessing.R`**
   - **Input:** `Raw_Data/CGGA_GBM/CGGA_693_Expression_Cleaned.tsv`
   - **Output:** `Processed_Data/CGGA_GBM_Master_Pipeline/`
   - **Functions:**
     - Same as TCGA preprocessing, adapted for CGGA data
     - Handles CGGA-specific metadata structure
   - **Runtime:** ~2-3 minutes

3. **`03_validate_clean_data.R`**
   - **Input:** Processed data from above scripts
   - **Output:** Console validation report
   - **Functions:**
     - Validate data dimensions
     - Check DGAT1/DGAT2 presence
     - Verify metadata completeness
     - QC summary statistics

---

## 02_batch_correction/

**Purpose:** Remove batch effects while preserving biological signal

### Scripts:

1. **`01_tcga_batch_correction.R`**
   - **Input:** `Processed_Data/TCGA_GBM_Master_Pipeline/`
   - **Output:** `Processed_Data/TCGA_GBM_Batch_Corrected/`
   - **Method:** ComBat batch correction
   - **Covariates preserved:** Age, purity
   - **Features:**
     - Handles singleton batches
     - Imputes missing age/purity values
     - Pre/post-correction comparison plots
     - Correlation heatmaps
   - **Runtime:** ~3-5 minutes

2. **`02_cgga_batch_correction.R`**
   - **Input:** `Processed_Data/CGGA_GBM_Master_Pipeline/`
   - **Output:** `Processed_Data/CGGA_GBM_Batch_Corrected/`
   - **Method:** ComBat batch correction
   - **Covariates preserved:** Age only (CGGA has limited metadata)
   - **Runtime:** ~3-5 minutes

---

## 03_analysis/

**Purpose:** Downstream biological analysis

### Scripts:

1. **`01_dgat_immune_analysis.R`**
   - **Input:** Batch-corrected expression data
   - **Output:** `Results/Bulk/Immune_Analysis/`
   - **Analyses:**
     - DGAT1/DGAT2 expression calculation
     - Immune deconvolution (GSVA with 28 curated gene sets)
     - Differential immune infiltration (DGAT-high vs DGAT-low)
     - Spearman correlations (DGAT vs immune signatures)
     - Volcano plots, heatmaps, boxplots
   - **Gene Sets:** T cells, NK cells, Myeloid, B cells, checkpoints, IFN, complement
   - **Runtime:** ~5-10 minutes per dataset

2. **`02_survival_analysis.R`**
   - **Input:** Batch-corrected expression + clinical data
   - **Output:** `Results/Bulk/Survival/`
   - **Analyses:**
     - Auto-detection of survival columns
     - Optimal cutpoint determination (`surv_cutpoint`)
     - Kaplan-Meier curves
     - Cox proportional hazards regression
     - Multivariate analysis (age, IDH, MGMT adjustment)
   - **Runtime:** ~3-5 minutes per dataset

---

## Typical Workflow

```r
# Step 1: Preprocess TCGA data
source("Scripts/Transcriptomics/01_preprocessing/01_tcga_fpkm_preprocessing.R")

# Step 2: Batch correction
source("Scripts/Transcriptomics/02_batch_correction/01_tcga_batch_correction.R")

# Step 3: Validate processed data
source("Scripts/Transcriptomics/01_preprocessing/03_validate_clean_data.R")

# Step 4: DGAT-immune analysis
source("Scripts/Transcriptomics/03_analysis/01_dgat_immune_analysis.R")

# Step 5: Survival analysis
source("Scripts/Transcriptomics/03_analysis/02_survival_analysis.R")

# Repeat for CGGA dataset
```

---

## Data Flow

```
Raw_Data/
  ├── TCGA_GBM/
  └── CGGA_GBM/
       ↓
01_preprocessing/
       ↓
Processed_Data/.../Master_Pipeline/
       ↓
02_batch_correction/
       ↓
Processed_Data/.../Batch_Corrected/
       ↓
03_analysis/
       ↓
Results/Bulk/
```

---

## Key Features

- **Modular design:** Each script is standalone and well-documented
- **Robust QC:** Multiple quality control checkpoints
- **Batch effect handling:** ComBat with covariate preservation
- **Reproducible:** All parameters logged and saved
- **Flexible:** Easy to adapt for new datasets

---

## Dependencies

**R Packages Required:**
- `data.table`, `dplyr` - Data manipulation
- `edgeR`, `limma` - Gene filtering
- `sva` - Batch correction (ComBat)
- `GSVA` - Immune deconvolution
- `survival`, `survminer` - Survival analysis
- `ggplot2`, `pheatmap` - Visualization

---

## Notes

- All scripts use absolute paths to ensure reproducibility
- TCGA and CGGA pipelines are parallel but dataset-specific
- Batch-corrected data is recommended for all downstream analyses
- CGGA has limited clinical metadata compared to TCGA

---

**Last Updated:** 2025-10-10  
**Status:** Production-ready


