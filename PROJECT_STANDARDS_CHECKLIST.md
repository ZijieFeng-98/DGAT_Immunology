# DGAT1 Project Standards & Implementation Checklist

**Last Updated:** October 20, 2025  
**Status:** Active workflow baseline for all analyses  
**Reference:** Cancer Immunology Bioinformatics Standards v1.1

---

## 📋 CURRENT PROJECT STATUS

### ✅ Completed (Don't Redo)
- [x] Sample QC and filtering (TCGA, CGGA, CPTAC)
- [x] Batch correction (ComBat-seq)
- [x] GSVA immune signatures
- [x] Univariate Cox survival analysis
- [x] Basic correlation analyses
- [x] Data preprocessing pipelines
- [x] Daily logging system installed

### ⭐⭐⭐ WEEK 1 PRIORITIES (Critical - Do First)

#### Priority 1: Multivariable Cox Models (Day 1)
**Status:** NOT STARTED  
**Time:** 1 day  
**Required in:** 70% of papers  

**Implementation Checklist:**
- [ ] Define all clinical covariates table:
  - [ ] Age (continuous, years)
  - [ ] Sex (categorical, M/F)
  - [ ] WHO Grade (categorical, III/IV)
  - [ ] IDH Status (categorical, mut/WT)
  - [ ] MGMT Methylation (categorical, meth/unmeth)
  - [ ] Tumor Purity (continuous, ESTIMATE scores)
- [ ] Build multivariable Cox: `coxph(Surv ~ DGAT1 + age + sex + grade + IDH + MGMT + purity)`
- [ ] Test proportional hazards assumption (Schoenfeld residuals)
- [ ] Generate forest plot (HR with 95% CI)
- [ ] Compare univariate vs multivariable results
- [ ] Document any violations
- [ ] Repeat for both TCGA and CGGA

**Output Files:**
- `Results/Survival/TCGA_GBM/multivariable_cox_results.csv`
- `Results/Survival/TCGA_GBM/forest_plot_multivariable.pdf`
- `Results/Survival/CGGA_GBM/multivariable_cox_results.csv`

---

#### Priority 2: Batch Correction Validation Metrics (Day 1)
**Status:** PARTIALLY DONE (need formal table)  
**Time:** 2 hours  
**Required in:** 80-90% of papers  

**Implementation Checklist:**
- [ ] Calculate PC1-PC10 vs batch correlations (before/after)
- [ ] Kruskal-Wallis tests (PC ~ batch)
- [ ] Variance explained by PC1-PC10 (before/after)
- [ ] Document ComBat-seq parameters used
- [ ] Create publication-ready table
- [ ] Confirm biology preserved (IDH/Grade clustering)

**Output Files:**
- `Processed_Data/TCGA_GBM_Batch_Corrected/batch_validation_metrics.csv`
- `Processed_Data/CGGA_GBM_Batch_Corrected/batch_validation_metrics.csv`
- `Results/QC_Reports/batch_correction_validation_table.csv`

**Table Format:**
```
Metric                                  Before      After
─────────────────────────────────────────────────────────
PC1-batch correlation (Spearman)       0.376***    0.144*
PC2-batch correlation                  0.089       0.067
PC1 variance explained                 26.7%       22.4%
Kruskal-Wallis p-value (PC1~batch)    <0.001      0.089
Known biology preserved                -           IDH clustering ✓
```

---

#### Priority 3: Tumor Purity Assessment (Day 2)
**Status:** NOT STARTED  
**Time:** 2 hours  
**Required in:** 60-70% of papers  

**Implementation Checklist:**
- [ ] Calculate ESTIMATE purity scores (if not already done)
- [ ] Test: `cor.test(DGAT1_expr, purity, method="spearman")`
- [ ] Test: `cor.test(DGAT1_expr, immune_infiltration, method="spearman")`
- [ ] Document results:
  - [ ] If r < 0.3 AND p > 0.05: "No evidence of purity confounding"
  - [ ] If r > 0.3 OR p < 0.05: Plan sensitivity analysis (Week 2)
- [ ] Include purity as covariate in ALL Cox models
- [ ] Report in Results section

**Output Files:**
- `Results/Purity_Analysis/DGAT1_vs_purity_correlation.csv`
- `Results/Purity_Analysis/purity_scatter_plots.pdf`

**Decision Logic:**
```R
if (abs(cor_DGAT1_purity) < 0.3 & p_value > 0.05) {
  action <- "Include purity as covariate only"
  report <- "No evidence of purity confounding"
} else {
  action <- "Include purity as covariate + sensitivity analysis"
  report <- "DGAT1 shows correlation with purity; adjusted in models"
}
```

---

#### Priority 4: CIBERSORTx Validation (Day 3)
**Status:** NOT STARTED  
**Time:** 1 day  
**Required in:** 40-60% of papers  

**Implementation Checklist:**
- [ ] Download LM22 signature matrix
- [ ] Run CIBERSORTx (absolute mode, 1000 permutations)
- [ ] Filter samples: p < 0.05
- [ ] Extract key cell types:
  - [ ] M2 macrophages
  - [ ] TAMs
  - [ ] CD8 T cells
  - [ ] MDSCs
- [ ] Correlate DGAT1 with cell proportions
- [ ] Cross-validate with GSVA:
  - [ ] GSVA TAM vs CIBERSORTx M2: expect r > 0.5
  - [ ] GSVA CD8 vs CIBERSORTx CD8: expect r > 0.6
- [ ] Create concordance table
- [ ] Report both methods in manuscript

**Output Files:**
- `Results/Immune_Deconvolution/CIBERSORTx_results.csv`
- `Results/Immune_Deconvolution/GSVA_vs_CIBERSORTx_concordance.csv`
- `Results/Immune_Deconvolution/method_comparison_plots.pdf`

---

### ⭐⭐ WEEK 2 PRIORITIES (High Priority)

#### Priority 5: Patient Characteristics Table (Day 1)
**Status:** NOT STARTED  
**Time:** 1 day  
**Required in:** >95% of papers (Table 1)  

**Implementation Checklist:**
- [ ] Create Table 1 for each cohort (TCGA, CGGA, CPTAC)
- [ ] Include for each:
  - [ ] Sample size (N)
  - [ ] Age: median (IQR)
  - [ ] Sex: N (%) male/female
  - [ ] WHO Grade: N (%) III/IV
  - [ ] IDH mutation: N (%) mut/WT
  - [ ] MGMT methylation: N (%) meth/unmeth
  - [ ] Median follow-up (months)
  - [ ] Overall survival: median (95% CI)
  - [ ] Event rate: N (%) deaths
  - [ ] DGAT1 expression: median (IQR)
- [ ] Statistical comparison across cohorts (chi-square, Kruskal-Wallis)
- [ ] Format as publication-ready table

**Output Files:**
- `Results/Tables/Table1_Patient_Characteristics.csv`
- `Results/Tables/Table1_Patient_Characteristics.docx`

---

#### Priority 6: CGGA Validation Cohort (Days 2-4)
**Status:** 30% DONE (preprocessing)  
**Time:** 2-3 days  
**Required in:** 50-70% of papers  

**Implementation Checklist:**
- [ ] Apply TCGA-derived median cutpoint to CGGA
- [ ] Multivariable Cox with same covariates
- [ ] Compare HR direction and magnitude
- [ ] Meta-analysis:
  - [ ] Random effects model (DerSimonian-Laird)
  - [ ] Calculate pooled HR with 95% CI
  - [ ] Assess heterogeneity (I² statistic)
  - [ ] Forest plot (TCGA + CGGA + pooled)
- [ ] Interpret results:
  - [ ] Both p<0.05, same direction → Strong evidence ✓
  - [ ] One p<0.05, same direction → Supportive
  - [ ] Opposite direction → Investigate heterogeneity

**Output Files:**
- `Results/Validation/CGGA_survival_validation.csv`
- `Results/Validation/meta_analysis_results.csv`
- `Results/Validation/forest_plot_meta_analysis.pdf`

---

#### Priority 7: Sensitivity Analyses (Days 5-7)
**Status:** NOT STARTED  
**Time:** 2-3 days  
**Required in:** 40-60% of papers  

**Implementation Checklist:**

**A. Cutpoint Selection:**
- [ ] Median (primary analysis) ✓
- [ ] Tertiles
- [ ] Quartiles
- [ ] Optimal (survminer::surv_cutpoint)
- [ ] Bootstrap stability (1000 iterations)

**B. IDH Stratification:**
- [ ] Separate analyses: IDH-mutant vs IDH-wildtype
- [ ] Repeat survival and correlation analyses
- [ ] Document associations within strata

**C. Purity Stratification (if r > 0.3):**
- [ ] Stratify by purity tertiles (low/medium/high)
- [ ] Repeat survival and correlation within strata
- [ ] Report: "Associations remained significant"

**D. Method Concordance:**
- [ ] Compare GSVA vs CIBERSORTx
- [ ] Correlation statistics
- [ ] Document consistency

**Output Files:**
- `Results/Sensitivity/cutpoint_comparison.csv`
- `Results/Sensitivity/IDH_stratified_analyses.csv`
- `Results/Sensitivity/purity_stratified_analyses.csv`
- `Results/Sensitivity/method_concordance.csv`

---

### ⭐ WEEK 2-3: FIGURES & DOCUMENTATION

#### Priority 8: Main Figures (Days 1-3)

**Figure 1: Study Design & QC**
- [ ] Panel A: CONSORT flow diagram
- [ ] Panel B: PCA before batch correction (colored by plate)
- [ ] Panel C: PCA after batch correction
- [ ] Panel D: PCA colored by IDH (biology preserved)
- [ ] Panel E: Batch validation metrics

**Figure 2: Survival Analysis**
- [ ] Panel A: Kaplan-Meier (TCGA)
- [ ] Panel B: Forest plot (multivariable Cox, TCGA)
- [ ] Panel C: Kaplan-Meier (CGGA validation)
- [ ] Panel D: Meta-analysis forest plot
- [ ] Panel E: DGAT1 distribution by survival group

**Figure 3: Immune Landscape**
- [ ] Panel A: Heatmap (GSVA scores)
- [ ] Panel B: Violin plots (M2 TAM, CD8 by DGAT1)
- [ ] Panel C: CIBERSORTx proportions
- [ ] Panel D: Correlation network
- [ ] Panel E: DGAT1 vs DGAT2 comparison

**Supplementary Figures:**
- [ ] S1: Patient characteristics comparison
- [ ] S2: Known gene validation
- [ ] S3: Sensitivity analyses
- [ ] S4: Purity correlations
- [ ] S5: Full correlation matrix
- [ ] S6: IDH-stratified analyses

---

#### Priority 9: Methods Writing (Days 4-5)

**Required Sections:**
- [ ] Data Acquisition and QC
- [ ] Batch Correction
- [ ] Immune Profiling (GSVA + CIBERSORTx)
- [ ] Survival Analysis (Univariate + Multivariable)
- [ ] Confounder Assessment
- [ ] Statistical Analysis
- [ ] Sensitivity Analyses (NEW)
- [ ] Limitations (NEW)
- [ ] Data Availability

**Use Templates From:** Standards document Section 7

---

#### Priority 10: Code Repository (Days 6-7)

**GitHub Structure:**
```
DGAT1-GBM-Immunology/
├── README.md
├── environment.yml
├── sessionInfo.txt
├── data/
│   └── README.md (links only, no raw data)
├── scripts/
│   ├── 01_preprocessing/
│   ├── 02_analysis/
│   ├── 03_validation/
│   └── 04_figures/
├── results/
│   ├── tables/
│   ├── figures/
│   └── supplementary/
└── docs/
    └── methods_detailed.md
```

**Checklist:**
- [ ] Clean and document all code
- [ ] Add README with usage instructions
- [ ] Include sessionInfo.txt
- [ ] Create environment.yml
- [ ] Add .gitignore for large files
- [ ] Commit with descriptive messages
- [ ] Tag release: v1.0-submission

---

## 🔄 OPTIONAL ENHANCEMENTS (If Time Permits or Reviewers Request)

### Enhancement 1: Multi-Algorithm Deconvolution
**Status:** OPTIONAL (reviewer may request)  
**Time:** 3-5 days  
**Found in:** 30-40% of papers  

- [ ] Add TIMER2
- [ ] Add EPIC
- [ ] Add xCell
- [ ] Calculate ensemble scores
- [ ] Report cross-method concordance

**Decision:** Skip for initial submission. Plan for revision if requested.

---

### Enhancement 2: Expression Pre-Adjustment for Confounders
**Status:** OPTIONAL (<10% of papers do this)  
**Time:** 2-3 days  

```R
# Limma residuals approach (if requested)
design <- model.matrix(~ tumor_purity + IDH + MGMT, data = metadata)
fit <- lmFit(expr_matrix, design)
expr_adjusted <- residuals(fit, expr_matrix)
```

**Decision:** Skip. We include covariates in multivariable models (standard practice).

---

### Enhancement 3: TPM Re-Normalization
**Status:** OPTIONAL (FPKM acceptable when documented)  
**Time:** 1-2 days  

**Decision:** Skip. FPKM is GDC standard format, widely accepted.

---

## 📊 MANDATORY ANALYSIS STANDARDS

### For ALL Survival Analyses:
✅ **Always include:**
- Univariate Cox (screening)
- **Multivariable Cox** with clinical covariates (primary claim)
- Kaplan-Meier curves
- Log-rank test
- Proportional hazards assumption testing
- HR with 95% CI

### For ALL Correlation Analyses:
✅ **Always include:**
- Spearman correlation (primary)
- FDR correction (Benjamini-Hochberg)
- Sample size (N)
- Report r and p-value

### For ALL Batch Corrections:
✅ **Always include:**
- Method description (ComBat-seq parameters)
- PC-batch correlation (before/after)
- PCA plots (before/after)
- Biology preservation check
- Validation metrics table

### For ALL Immune Analyses:
✅ **Always include:**
- Clear method description (GSVA/CIBERSORTx)
- Gene sets/signatures used
- Validation approach
- Confounder assessment
- FDR correction

---

## 📝 METHODS TEXT TEMPLATES

### Template: Multivariable Cox Model

```
Multivariable Cox proportional hazards regression was performed to 
evaluate the independent prognostic value of DGAT1 expression while 
adjusting for established clinical covariates: age (continuous, years), 
sex (male/female), WHO grade (III/IV), IDH mutation status (mutant/
wild-type), MGMT promoter methylation status (methylated/unmethylated), 
and tumor purity (continuous, ESTIMATE scores). Hazard ratios (HR) with 
95% confidence intervals (CI) were calculated. The proportional hazards 
assumption was tested using Schoenfeld residuals. All survival analyses 
were performed using the survival (v3.5-5) and survminer (v0.4.9) R 
packages.
```

### Template: Batch Correction Validation

```
Technical batch effects from sequencing plates were assessed by principal 
component analysis and correlation testing between the first 10 principal 
components and plate assignments. Significant batch effects were identified 
(PC1-plate Spearman correlation = 0.376, p < 0.001). Batch correction was 
performed using ComBat-seq from the sva package (v3.48.0) on log2-
transformed FPKM values. Correction efficacy was validated by: (1) reduced 
PC-batch correlation (0.376 → 0.144), (2) visual inspection of PCA plots, 
and (3) preservation of known biological relationships (IDH mutation status 
clustering). [Report specific metrics from validation table]
```

### Template: Tumor Purity Assessment

```
To ensure observed associations between DGAT1 expression and immune 
infiltration were not driven by tumor microenvironment composition, we 
tested correlations between DGAT1 and tumor purity (ESTIMATE algorithm). 
DGAT1 expression showed [no significant/weak/moderate] correlation with 
tumor purity (Spearman r = [X.XX], p = [X.XX]), suggesting immune 
associations were [largely] independent of tumor cell fraction. All 
multivariable Cox models included tumor purity as a continuous covariate 
to adjust for potential confounding effects on survival.
```

---

## 🚨 QUALITY CONTROL GATES

### Before Proceeding to Next Phase:

**Gate 1: After Week 1 Priorities**
- [ ] Multivariable Cox models complete for TCGA
- [ ] Batch validation table generated
- [ ] Tumor purity assessment done
- [ ] CIBERSORTx validation complete
- [ ] All outputs documented and saved

**Gate 2: After Week 2 Priorities**
- [ ] CGGA validation complete
- [ ] Patient characteristics table done
- [ ] Sensitivity analyses run
- [ ] All main figures generated
- [ ] Methods section written

**Gate 3: Before Submission**
- [ ] All figures finalized (PDF + PNG)
- [ ] All tables in publication format
- [ ] Methods section complete
- [ ] Limitations section written
- [ ] Code repository on GitHub
- [ ] Data availability statement ready
- [ ] All citations formatted
- [ ] Supplementary materials organized

---

## 📈 PUBLICATION READINESS TRACKER

| Component | Required % | Status | Priority | Timeline |
|-----------|-----------|--------|----------|----------|
| Data preprocessing | Universal | ✅ 95% | - | Complete |
| Batch correction | 80-90% | ✅ 95% | ⭐ | Add table (2h) |
| **Multivariable Cox** | **70%** | **❌ 0%** | **⭐⭐⭐** | **Week 1 (1d)** |
| GSVA signatures | 60-80% | ✅ 100% | - | Complete |
| **CIBERSORTx** | **40-60%** | **❌ 0%** | **⭐⭐** | **Week 1 (1d)** |
| **Purity assessment** | **60-70%** | **❌ 0%** | **⭐⭐** | **Week 1 (2h)** |
| Patient table | >95% | ❌ 0% | ⭐ | Week 2 (1d) |
| CGGA validation | 50-70% | 🟡 30% | ⭐⭐ | Week 2 (2-3d) |
| Sensitivity analyses | 40-60% | ❌ 0% | ⭐ | Week 2 (2-3d) |
| Main figures | Universal | 🟡 40% | ⭐⭐ | Week 2 (2-3d) |
| Methods text | Universal | 🟡 20% | ⭐ | Week 2 (1-2d) |
| Code repository | 70-80% | 🟡 10% | ⭐ | Week 2 (1d) |

**Overall Status:**
- **Current:** 60% publication-ready
- **After Week 1:** 85% publication-ready
- **After Week 2:** 95% publication-ready

---

## 🎯 CRITICAL SUCCESS FACTORS

### For Nature/Cell-Level Publication:

1. **Multivariable adjustment is NON-NEGOTIABLE** for primary survival claims
2. **Independent validation** (CGGA) massively strengthens conclusions
3. **Method validation** (GSVA + CIBERSORTx) demonstrates rigor
4. **Transparent documentation** of all methodological choices
5. **Confounder assessment** shows you understand potential biases
6. **Reproducibility** through code/data sharing

### What Reviewers WILL Ask For:

1. "Did you adjust for clinical covariates in your survival model?" → **YES (Week 1)**
2. "How did you validate your batch correction?" → **Metrics table (Week 1)**
3. "Could tumor purity confound these associations?" → **Tested (Week 1)**
4. "Can you validate with another deconvolution method?" → **CIBERSORTx (Week 1)**
5. "Do findings replicate in independent cohort?" → **CGGA (Week 2)**
6. "Are results robust to methodological choices?" → **Sensitivity (Week 2)**

---

## 📚 KEY CITATIONS (Copy-Paste Ready)

**Batch Correction:**
- Zhang Y et al. ComBat-seq: batch effect adjustment for RNA-seq count data. NAR Genomics Bioinform. 2020

**Immune Deconvolution:**
- Newman AM et al. Robust enumeration of cell subsets from tissue expression profiles. Nat Methods. 2015
- Hänzelmann S et al. GSVA: gene set variation analysis. BMC Bioinformatics. 2013

**Tumor Purity:**
- Yoshihara K et al. Inferring tumour purity and stromal and immune cell admixture. Nat Commun. 2013

**Survival Analysis:**
- Therneau TM. A Package for Survival Analysis in R. 2023

**FDR Correction:**
- Benjamini Y, Hochberg Y. Controlling the false discovery rate. J Royal Stat Soc B. 1995

---

## 🔄 VERSION CONTROL

**This Document:**
- Version: 1.0
- Created: October 20, 2025
- Last Updated: October 20, 2025
- Status: Active baseline for all analyses
- Next Review: After Week 1 priorities completion

**Project Status:**
- Phase: Implementation (Week 1 priorities)
- Target Submission: [Target date after Week 2 completion]
- Journal Target: Nature Communications / Cancer Cell / Cell Metabolism

---

**All future analyses and scripts MUST follow these standards.**

**End of Standards Checklist**

