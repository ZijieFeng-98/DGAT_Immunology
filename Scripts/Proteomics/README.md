# Proteomics Analysis Pipeline

**Data Type:** Protein Expression (Mass Spectrometry)  
**Datasets:** CPTAC3 GBM Discovery Cohort (TMT11-plex)  
**Purpose:** Validate DGAT-immune relationships at the protein level

---

## Directory Structure

```
Proteomics/
├── 01_cptac_preprocessing.R        # Main processing pipeline
├── 02_cptac_troubleshooting.R      # Diagnostic tool
└── 03_hpa_dgat_query.R             # HPA database query
```

---

## Scripts

### 1. `01_cptac_preprocessing.R` ✅ PRODUCTION-READY

**Purpose:** Process and QC CPTAC3 GBM proteomics data

**Input:**
- `Raw_Data/HPA_Protein/CPTAC3_Glioblastoma_Multiforme_Proteome.tmt11.tsv`
- `Raw_Data/HPA_Protein/S048_CPTAC_GBM_Discovery_Cohort_TMT11_CaseID_SampleID_AliquotID_Map_Dec2019_r1.xlsx`
- `Raw_Data/HPA_Protein/S048_CPTAC_GBM_Discovery_Cohort_Clinical_Data_Dec2019_r1.xlsx`

**Output:**
- `Processed_Data/CPTAC_GBM_Proteomics/protein_matrix_cleaned.rds`
- `Processed_Data/CPTAC_GBM_Proteomics/clinical_data_matched.rds`
- `Processed_Data/CPTAC_GBM_Proteomics/QC_Plots/` (5 plots)
- `Processed_Data/CPTAC_GBM_Proteomics/QC_Summary.csv`

**Processing Steps:**
1. Load TMT11-plex proteomics data (10,980 proteins)
2. Map sample IDs (Aliquot → Case ID)
3. Filter samples (remove QC/pooled samples)
4. Deduplicate patients (one sample per patient)
5. Filter proteins (>50% missing → removed)
6. Match with clinical data (110/110 patients)
7. Generate comprehensive QC plots

**Final Output:**
- **Proteins:** 10,804 (after QC)
- **Samples:** 110 unique GBM patients
- **DGAT1:** Detected (18.2% missing)
- **Immune markers:** 4/5 detected (CD274, HAVCR2, MRC1, ARG1)
- **Clinical match:** 100%

**QC Plots Generated:**
1. Missing value distribution
2. Sample-sample correlations
3. PCA analysis
4. DGAT1 expression distribution
5. DGAT-immune protein correlation heatmap

**Runtime:** ~2-3 minutes

**Status:** ✅ Validated and production-ready  
**Bugs:** All fixed (see troubleshooting script)

---

### 2. `02_cptac_troubleshooting.R` 🔍 DIAGNOSTIC TOOL

**Purpose:** Diagnose missing genes and data loading issues

**Use Cases:**
- Gene not found in protein matrix
- Clinical data matching problems
- Excel file structure issues
- Gene symbol verification

**Features:**
- **6 search strategies for missing genes:**
  1. Exact match
  2. Case-insensitive search
  3. Partial match
  4. Alternative gene names
  5. Fuzzy matching (approximate grep)
  6. Missing value check (filtered genes)

- **Clinical data diagnosis:**
  - Tests multiple Excel loading strategies
  - Identifies correct sheet names
  - Finds proper Case ID columns

- **Gene discovery:**
  - Searches for DGAT1/DGAT2
  - Searches for 8 immune checkpoint markers
  - Searches for 10 lipid metabolism markers

**Output:**
- `Processed_Data/CPTAC_GBM_Proteomics/Troubleshooting/troubleshooting_summary.csv`
- `Processed_Data/CPTAC_GBM_Proteomics/Troubleshooting/troubleshooting_findings.rds`
- `Processed_Data/CPTAC_GBM_Proteomics/Troubleshooting/immune_markers_found.csv`

**When to Run:**
- After preprocessing if key genes are missing
- Before updating gene symbol lists
- When clinical data doesn't match
- For any data loading diagnostics

**Runtime:** ~1 minute

---

### 3. `03_hpa_dgat_query.R` 🔎 DATABASE QUERY

**Purpose:** Query Human Protein Atlas (HPA) for DGAT protein expression

**Data Source:** HPA database (www.proteinatlas.org)

**Features:**
- Download DGAT1/DGAT2 protein expression data
- Tissue-specific expression profiles
- Immunohistochemistry images
- Subcellular localization

**Output:** HPA query results and protein expression profiles

**Note:** For exploratory analysis and literature review

---

## Typical Workflow

```r
# Step 1: Process CPTAC proteomics data
source("Scripts/Proteomics/01_cptac_preprocessing.R")

# Step 2: (If needed) Troubleshoot missing genes
source("Scripts/Proteomics/02_cptac_troubleshooting.R")

# Step 3: Validate processed data
protein <- readRDS("Processed_Data/CPTAC_GBM_Proteomics/protein_matrix_cleaned.rds")
clinical <- readRDS("Processed_Data/CPTAC_GBM_Proteomics/clinical_data_matched.rds")

# Check DGAT1 detection
"DGAT1" %in% rownames(protein)  # Should be TRUE

# Step 4: Analyze DGAT-immune protein correlations
# (Use protein data for correlation analysis, survival analysis, etc.)
```

---

## Data Flow

```
Raw_Data/HPA_Protein/
  ├── CPTAC3_Glioblastoma_Multiforme_Proteome.tmt11.tsv
  ├── S048_..._TMT11_CaseID_SampleID_AliquotID_Map_...xlsx
  └── S048_..._Clinical_Data_...xlsx
       ↓
01_cptac_preprocessing.R
       ↓
Processed_Data/CPTAC_GBM_Proteomics/
  ├── protein_matrix_cleaned.rds      (10,804 × 110)
  ├── clinical_data_matched.rds       (110 × 43)
  ├── QC_Plots/ (5 plots)
  └── QC_Summary.csv
       ↓
Downstream Analysis:
  • DGAT-immune protein correlations
  • Protein-based survival analysis
  • mRNA-protein concordance
  • Multi-omic integration
```

---

## Key Features

- **Robust QC:** 5 comprehensive quality control plots
- **Complete clinical matching:** 110/110 patients (100%)
- **Gene detection:** DGAT1 + 4 immune markers validated
- **Troubleshooting tools:** Built-in diagnostics for missing genes
- **Production-ready:** All critical bugs fixed and validated

---

## Data Characteristics

**Technology:** TMT11-plex Mass Spectrometry  
**Normalization:** Log2 ratio to pooled reference  
**Coverage:** ~11,000 proteins quantified  
**Batch structure:** 12 TMT batches  
**Missing values:** Controlled (<50% per protein threshold)

---

## Validated Markers

**DGAT Proteins:**
- ✅ DGAT1 (median = 0.124, coverage = 81.8%)
- ❌ DGAT2 (not detected - below MS detection limit)

**Immune Checkpoint Markers:**
- ✅ CD274 (PD-L1) - 100% coverage
- ✅ HAVCR2 (TIM-3) - 81.8% coverage
- ✅ MRC1 (CD206, M2 macrophage) - 100% coverage
- ✅ ARG1 (immunosuppressive) - 100% coverage
- ❌ CTLA4 (not detected)
- ❌ LAG3 (not detected)

**Lipid Metabolism Markers:**
- ✅ All 8 markers detected: DGAT1, PLIN2, PLIN3, CPT1A, FASN, ACACA, SCD, SREBF1

---

## Known Limitations

1. **DGAT2 not detected** - Likely below MS detection limit or not expressed in GBM
2. **Some immune markers missing** - CTLA4, LAG3 not in CPTAC panel
3. **Missing values** - Proteomics inherently has more missing data than RNA-seq
4. **Smaller sample size** - 110 patients vs ~150-600 in TCGA/CGGA RNA-seq

---

## Integration with Transcriptomics

CPTAC samples are a **subset** of TCGA-GBM cohort, enabling:

- mRNA-protein correlation analysis
- Validation of transcriptomic findings at protein level
- Post-transcriptional regulation studies
- Multi-omic integration

**Note:** Sample matching requires CPTAC-TCGA mapping file (not yet obtained)

---

## Dependencies

**R Packages Required:**
- `data.table`, `dplyr` - Data manipulation
- `readxl` - Excel file reading
- `ggplot2`, `pheatmap`, `RColorBrewer` - Visualization
- `survival`, `survminer` - Survival analysis (for downstream)

---

## Troubleshooting

**Problem:** DGAT1 not found in protein matrix  
**Solution:** Run `02_cptac_troubleshooting.R` to diagnose  
**Likely cause:** Gene symbol mismatch or data loading bug  
**Fix:** Check troubleshooting output for correct gene name

**Problem:** Clinical data not matching  
**Solution:** Run troubleshooting to identify correct Excel sheet and column names  
**Current fix:** Use sheet="Clinical_Attributes", column="case_id"

**Problem:** Immune markers not detected  
**Solution:** Run troubleshooting to find actual gene symbols in dataset  
**Note:** Some markers may not be in proteomics panel

---

## Citation

**Original Publication:**
```
Wang LB, et al. (2021). Proteogenomic and metabolomic characterization 
of human glioblastoma. Cancer Cell, 39(4):509-528.e20.
DOI: 10.1016/j.ccell.2021.01.006
```

**Data Source:**
```
Clinical Proteomic Tumor Analysis Consortium (CPTAC)
CPTAC3 Glioblastoma Discovery Study
National Cancer Institute, December 2019
```

---

**Last Updated:** 2025-10-10  
**Status:** ✅ Production-ready & Validated  
**Validation Certificate:** See `Processed_Data/CPTAC_GBM_Proteomics/DATA_VALIDATION_STAMP.txt`


