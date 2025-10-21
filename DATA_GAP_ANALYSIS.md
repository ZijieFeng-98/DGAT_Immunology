# 📊 Data Gap Analysis - Pipeline Requirements Assessment

**Date:** October 20, 2025  
**Analysis:** Comparison of existing processed data vs standards requirements  
**Reference:** `PROJECT_STANDARDS_CHECKLIST.md`

---

## ✅ EXECUTIVE SUMMARY

**Overall Status:** 🟡 **70% Ready** - Good foundation but critical covariates missing

**What's Working:**
- ✅ Batch correction complete and validated
- ✅ Expression matrices clean and ready
- ✅ Tumor purity estimates available (TCGA)
- ✅ Basic survival data available
- ✅ Sample size adequate (TCGA: 285, CGGA: 693, CPTAC: 110)

**Critical Gaps (MUST FIX for Week 1):**
- ❌ **TCGA: Missing IDH, MGMT, WHO Grade in batch-corrected metadata**
- ❌ **CGGA: Missing Gender (Sex) in processed clinical file**
- ❌ **CGGA: Missing tumor purity estimates**
- ❌ **Metadata files incomplete** - need to merge with full clinical data

**Impact:** Cannot run multivariable Cox models (Week 1 Priority #1) until these gaps are fixed.

---

## 📋 DETAILED ASSESSMENT BY COHORT

### 1. TCGA-GBM (N=285)

#### ✅ What You Have (GOOD):

| Component | Status | Location |
|-----------|--------|----------|
| **Expression data** | ✅ Complete | `Processed_Data/TCGA_GBM_Batch_Corrected/expression_batch_corrected.rds` |
| **Batch correction** | ✅ Validated | PC1-batch: 0.376 → 0.144 ✓ |
| **Sample size** | ✅ Adequate | 285 samples (excellent) |
| **Survival data** | ✅ Available | OS_days, OS_status in clinical file |
| **Age** | ✅ Available | In batch-corrected metadata (needs conversion from days) |
| **Gender** | ✅ Available | In batch-corrected metadata |
| **Tumor purity** | ✅ Available | ESTIMATE scores in metadata (0.53-0.89) |
| **Batch variable** | ✅ Available | For including in models if needed |

#### ❌ Critical Gaps (MUST FIX):

| Required Variable | Current Status | Location if Available | Action Needed |
|-------------------|----------------|------------------------|---------------|
| **IDH mutation status** | ❌ NOT in batch metadata | May be in GDC clinical | **PRIORITY 1: Extract from GDC or raw clinical** |
| **MGMT methylation** | ❌ NOT in batch metadata | May be in GDC clinical | **PRIORITY 1: Extract from GDC or raw clinical** |
| **WHO Grade** | ❌ NOT in batch metadata | Should be "GBM = Grade IV" mostly | **PRIORITY 1: Extract from GDC or raw clinical** |
| **Age (years)** | 🟡 Available but wrong format | Currently in DAYS (e.g., 20932) | **Convert: days/365.25 = years** |

**Current batch-corrected metadata columns:**
```
patient_id, tss, batch, age (in DAYS!), gender, purity
```

**Needed metadata columns for multivariable Cox:**
```
patient_id, age_years, gender, batch, purity, IDH_status, MGMT_status, WHO_grade, OS_days, OS_years, OS_status
```

---

### 2. CGGA-GBM (N=693)

#### ✅ What You Have (GOOD):

| Component | Status | Location |
|-----------|--------|----------|
| **Expression data** | ✅ Complete | `Processed_Data/CGGA_GBM_Batch_Corrected/expression_batch_corrected.rds` |
| **Batch correction** | ✅ Done | Batch effects corrected |
| **Sample size** | ✅ Excellent | 693 samples (very large cohort) |
| **Survival data** | ✅ Available | OS_time_years, OS_status in raw clinical |
| **Age** | ✅ Available | In both raw and processed |
| **IDH status** | ✅ AVAILABLE | In raw clinical: `Raw_Data/CGGA_GBM/CGGA_693_Clinical_Cleaned.tsv` |
| **MGMT status** | ✅ AVAILABLE | In raw clinical file (MGMTp_methylation_status) |
| **WHO Grade** | ✅ AVAILABLE | In raw clinical file (Grade: WHO II/III/IV) |

#### ❌ Critical Gaps (MUST FIX):

| Required Variable | Current Status | Location if Available | Action Needed |
|-------------------|----------------|------------------------|---------------|
| **Gender (Sex)** | ❌ Missing from processed | ✅ Available in RAW clinical | **PRIORITY 1: Merge from raw clinical** |
| **Tumor purity** | ❌ ALL NA | Not in CGGA data | **PRIORITY 2: Calculate ESTIMATE scores** |
| **Complete metadata** | ❌ Incomplete | Raw has all clinical data | **PRIORITY 1: Re-merge clinical data** |

**Current batch-corrected metadata columns:**
```
patient_id, age, gender (empty!), batch, purity (all NA)
```

**Available in raw clinical but NOT in processed:**
```
Grade, IDH_mutation_status, MGMTp_methylation_status, Gender, OS, Censor, Radio_status, Chemo_status, 1p19q_codeletion_status
```

**This is a DATA MERGE ISSUE** - all the data exists but wasn't properly merged into batch-corrected metadata!

---

### 3. CPTAC-GBM Proteomics (N=110)

#### ✅ What You Have (GOOD):

| Component | Status | Location |
|-----------|--------|----------|
| **Protein matrix** | ✅ Complete | 10,804 proteins × 110 samples |
| **Sample size** | ✅ Adequate | 110 patients (good for proteomics) |
| **DGAT1 protein** | ✅ Detected | 90/110 samples with data |
| **Quality control** | ✅ Validated | QC plots and validation done |
| **Processing** | ✅ Proper | TMT log2 ratios, cleaned |

#### 🟡 Gaps (Less Critical - CPTAC is exploratory):

| Variable | Status | Notes |
|----------|--------|-------|
| **Clinical data** | 🟡 Limited | README notes "0 matched in current run" |
| **Survival data** | ❓ Unknown | Not clear if available |
| **Molecular covariates** | ❓ Unknown | IDH/MGMT/Grade availability unclear |

**Note:** Per standards document, CPTAC is **optional enhancement** (Section: "Enhancement 3"). Prioritize TCGA/CGGA fixes first.

---

## 🚨 IMMEDIATE ACTION REQUIRED (BEFORE Week 1 Priorities)

You **CANNOT proceed with multivariable Cox models** (Week 1 Priority #1) until these data gaps are fixed.

### Action Plan: Data Integration (Est. 1-2 hours)

#### Task 1: Fix TCGA Metadata (30-45 min)

**Goal:** Create complete metadata file with all multivariable Cox covariates

**Steps:**
1. Load batch-corrected metadata (current)
2. Load full TCGA clinical data from GDC (has IDH, MGMT, Grade)
3. Merge the two by `patient_id`
4. Convert age from days to years: `age_years = age/365.25`
5. Clean IDH status (Mutant/Wildtype)
6. Clean MGMT status (Methylated/Unmethylated)
7. Add WHO grade (most will be "Grade IV" for GBM)
8. Add survival data (OS_days, OS_years, OS_status)
9. Save as: `Processed_Data/TCGA_GBM_Batch_Corrected/metadata_complete.csv`

**Expected output columns:**
```
patient_id, age_years, gender, batch, purity, IDH_status, MGMT_status, WHO_grade, OS_days, OS_years, OS_status
```

---

#### Task 2: Fix CGGA Metadata (30-45 min)

**Goal:** Properly merge raw clinical data into batch-corrected metadata

**Steps:**
1. Load batch-corrected metadata (current - incomplete)
2. Load raw clinical: `Raw_Data/CGGA_GBM/CGGA_693_Clinical_Cleaned.tsv`
3. Merge by `CGGA_ID` / `patient_id`
4. Extract:
   - Gender (from raw: "Gender" column)
   - IDH_mutation_status → standardize to "Mutant"/"Wildtype"
   - MGMTp_methylation_status → standardize to "Methylated"/"Unmethylated"  
   - Grade → keep "WHO II", "WHO III", "WHO IV"
   - OS → rename to OS_days (already in days)
   - Censor → convert to OS_status (0→0, 1→1)
5. Calculate tumor purity using ESTIMATE (if needed)
6. Save as: `Processed_Data/CGGA_GBM_Batch_Corrected/metadata_complete.csv`

**Expected output columns:**
```
patient_id, age_years, gender, batch, purity, IDH_status, MGMT_status, WHO_grade, OS_days, OS_years, OS_status
```

---

#### Task 3: Validate Merged Data (10-15 min)

**Quality checks to run:**
```r
# Check completeness
metadata <- read.csv("metadata_complete.csv")

# Required checks:
1. N samples matches expression matrix
2. All patient IDs in expression matrix have metadata
3. Key variables have <10% missing:
   - age_years: should be 0% missing
   - gender: should be 0% missing
   - IDH_status: <10% missing acceptable
   - MGMT_status: <20% missing acceptable (known issue in TCGA/CGGA)
   - OS_status: should be 0% missing
4. Survival data makes sense:
   - Event rate: 60-90% for GBM
   - Median OS: 0.8-2.5 years
5. Age distribution: 18-90 years (not 6000-25000 days!)
```

---

## 📊 STANDARDS REQUIREMENTS CHECKLIST

### Phase 1: Data Preprocessing ✅→🟡 (Need metadata fix)

| Requirement | TCGA | CGGA | CPTAC | Status |
|-------------|------|------|-------|--------|
| Sample QC | ✅ | ✅ | ✅ | **DONE** |
| Gene filtering | ✅ | ✅ | ✅ | **DONE** |
| Batch correction | ✅ | ✅ | N/A | **DONE** |
| Batch validation | 🟡 | 🟡 | N/A | **Need metrics table (Week 1)** |
| **Complete metadata** | ❌ | ❌ | 🟡 | **FIX FIRST** |
| Tumor purity | ✅ | ❌ | ❓ | **TCGA done, CGGA needs ESTIMATE** |

### Phase 2: Core Analyses (Depends on complete metadata)

| Requirement | TCGA | CGGA | CPTAC | Dependency |
|-------------|------|------|-------|------------|
| Univariate Cox | ✅ | ❓ | ❓ | **Needs survival data** ✓ |
| **Multivariable Cox** | ❌ | ❌ | ❌ | **Needs complete metadata** ← **BLOCKING** |
| GSVA | ✅ | ❓ | N/A | Can proceed |
| CIBERSORTx | ❌ | ❌ | N/A | **Week 1 priority** |
| Purity assessment | 🟡 | ❌ | ❓ | **Needs complete metadata** ← **BLOCKING** |

---

## 💡 RECOMMENDATION: DATA FIX BEFORE WEEK 1

### Revised Timeline:

**TODAY (Day 0): Data Integration** ← **ADD THIS STEP**
- ✅ Fix TCGA metadata (Task 1): 30-45 min
- ✅ Fix CGGA metadata (Task 2): 30-45 min  
- ✅ Validate merged data (Task 3): 10-15 min
- ✅ Calculate CGGA purity (ESTIMATE): 30 min
- **Total time: 2-3 hours**

**THEN Week 1 Priorities (as planned):**
- Day 1: Multivariable Cox models
- Day 2: Batch validation metrics + Purity assessment
- Day 3: CIBERSORTx validation

---

## 🔧 SCRIPTS NEEDED

I can generate these for you:

### 1. `Scripts/00_data_integration/01_merge_tcga_clinical.R`
- Merges batch-corrected data with full GDC clinical
- Extracts IDH, MGMT, Grade
- Converts age to years
- Adds survival data
- Validates completeness

### 2. `Scripts/00_data_integration/02_merge_cgga_clinical.R`
- Merges batch-corrected data with raw clinical
- Standardizes variable names and formats
- Calculates OS in years
- Validates completeness

### 3. `Scripts/00_data_integration/03_calculate_purity.R`
- Calculates ESTIMATE scores for CGGA
- Validates purity estimates
- Adds to metadata

### 4. `Scripts/00_data_integration/04_validate_metadata.R`
- Comprehensive QC checks
- Generates summary reports
- Confirms analysis-ready status

---

## ✅ NEXT STEP

**Would you like me to generate the data integration scripts?**

Once we fix the metadata (2-3 hours), you'll be 95% ready for Week 1 priorities.

**Command to start:**
```
"Generate the TCGA clinical data merge script"
```

This will unblock your entire Week 1 plan! 🚀

---

**Status after data fix:** 70% → 90% ready for Week 1 priorities

