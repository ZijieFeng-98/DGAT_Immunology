# CPTAC3 Glioblastoma Multiforme Proteomics Data

**Dataset:** CPTAC3 GBM Discovery Cohort  
**Source:** Clinical Proteomic Tumor Analysis Consortium (CPTAC)  
**Release:** December 2019 (Revision 1)  
**Technology:** TMT11-plex Mass Spectrometry

---

## Overview

This directory contains **protein-level quantification data** for glioblastoma multiforme (GBM) from the CPTAC3 Discovery Cohort. The data provides complementary proteomic information to the RNA-seq data (TCGA-GBM) in this project, enabling multi-omic analysis of DGAT expression and immune relationships at both the transcript and protein levels.

---

## Dataset Contents

### 📊 Proteomics Quantification Data

#### 1. **`CPTAC3_Glioblastoma_Multiforme_Proteome.tmt11.tsv`**
- **Size:** 10,981 rows × 227 columns
- **Content:** Protein abundance data (log2 ratios)
- **Format:** Tab-separated values

**Structure:**
- **Row 1 (Header):** Gene symbols followed by sample IDs with suffixes
  - Example columns: `Gene`, `CPT0204410003 Log Ratio`, `CPT0204410003 Unshared Log Ratio`, ...
- **Rows 2-10,981:** Protein quantification data
  - One protein (gene) per row
  - Log2 ratios normalized to pooled reference samples
  - Both "Log Ratio" and "Unshared Log Ratio" for each sample
    - **Log Ratio:** Abundance relative to pooled reference (all peptides)
    - **Unshared Log Ratio:** Abundance using unique peptides only (more specific)

**Summary Statistics (Row 2-4):**
```
Mean, Median, StdDev (calculated across all samples)
```

**Proteins:**
- 10,978 quantified proteins (genes)
- Examples: A1BG, A2M, AAAS, AARS, DGAT1, DGAT2, etc.
- Last columns: NCBIGeneID, Authority, Description, Organism, Chromosome, Locus

**Sample Coverage:**
- 100 GBM patient samples
- Each sample measured in duplicates (Log Ratio + Unshared Log Ratio)
- Total: 227 columns (1 Gene + 2 × 100 samples + 13 metadata columns)

---

#### 2. **`CPTAC3_Glioblastoma_Multiforme_Proteome.summary.tsv`**
- **Size:** 10,978 rows
- **Content:** Comprehensive protein metadata and quantification summary
- **Format:** Tab-separated values

**Contains:**
- Protein IDs (NCBI Gene ID, Gene Symbol)
- Protein descriptions
- Genomic locations (Chromosome, Locus)
- Summary statistics across all samples

---

### 🧬 Sample Metadata

#### 3. **`CPTAC3_Glioblastoma_Multiforme_Proteome.sample.txt`**
- **Size:** 13 rows (12 TMT batches + 1 header)
- **Content:** TMT11-plex experimental design

**Structure:**
```
FileNameRegEx | AnalyticalSample | 126C | 127N | 127C | ... | 131C | LabelReagent | Ratios
```

**Columns:**
- **FileNameRegEx:** Pattern to match raw MS files (e.g., `^01CPTAC_GBM_W_PNNL_20190123_B1S1`)
- **AnalyticalSample:** Batch identifier
- **126C - 131C:** Sample IDs assigned to each TMT channel (10 channels + 1 pooled reference)
  - **126C (Channel 1):** POOL (pooled reference for normalization)
  - **127N - 131C (Channels 2-11):** Individual patient samples (e.g., CPT0204410003, CPT0206670004, ...)
- **LabelReagent:** TMT reagent lot number (SG252258_SH258846)
- **Ratios:** Ratio calculation formula (all channels normalized to 126C/POOL)

**TMT Batches:**
- **12 batches total** (B1S1-B1S4, B2S1-B2S4, B3S1-B3S3)
  - Batch 1: 4 runs (samples processed 2019-01-23)
  - Batch 2: 4 runs (samples processed 2019-03-06)
  - Batch 3: 3 runs (samples processed 2019-05-01)
- **~110 unique patient samples** (10 samples per run × 11 runs)
- **Pooled reference** used in every run for batch effect correction

---

#### 4. **`S048_CPTAC_GBM_Discovery_Cohort_Clinical_Data_Dec2019_r1.xlsx`**
- **Content:** Clinical and demographic information
- **Format:** Excel spreadsheet

**Expected Variables:**
- Patient demographics (Age, Gender)
- Tumor characteristics (Grade, Histology, Location)
- Treatment history (Surgery, Chemotherapy, Radiation)
- Survival outcomes (Overall Survival, Progression-Free Survival)
- Molecular markers (IDH1/2 mutation, MGMT methylation, 1p/19q codeletion)

---

#### 5. **`S048_CPTAC_GBM_Discovery_Cohort_TMT11_CaseID_SampleID_AliquotID_Map_Dec2019_r1.xlsx`**
- **Content:** Complete sample mapping and provenance
- **Format:** Excel spreadsheet

**Purpose:**
- Links **Case IDs** (patients) to **Sample IDs** (tissue aliquots) to **Aliquot IDs** (proteomics runs)
- Maps TMT channel assignments to patient identifiers
- Tracks batch and plate information
- Provides quality control metrics

---

## Data Characteristics

### Technology Details

**TMT11-plex Mass Spectrometry:**
- **Labeling:** Tandem Mass Tag (TMT) isobaric labeling
- **Channels:** 11-plex (1 pooled reference + 10 samples per run)
- **Platform:** Likely Orbitrap Fusion or similar high-resolution MS
- **Quantification:** Reporter ion intensities normalized to pooled reference

**Advantages:**
- High-throughput (10 samples per run)
- Reduced missing values (vs. label-free)
- Controlled batch effects (pooled reference normalization)
- Protein-level isoform resolution

### Sample Characteristics

**Cohort Size:**
- ~100-110 GBM patient samples
- Discovery cohort (exploratory analysis)
- Matched with TCGA-GBM genomics data

**Tissue Type:**
- Primary glioblastoma tumors
- Frozen tissue samples
- CPTAC standardized sample processing

**Coverage:**
- ~11,000 proteins quantified
- Includes low-abundance proteins (e.g., transcription factors, signaling proteins)
- Suitable for targeted protein analysis (DGAT1, DGAT2, immune markers)

---

## Data Quality & Normalization

### Normalization Strategy

1. **Within-run normalization:**
   - All samples normalized to pooled reference (Channel 126C)
   - Log2 transformation applied
   
2. **Across-run normalization:**
   - Pooled reference ensures consistency across 12 TMT batches
   - Allows direct comparison of samples from different batches

3. **Peptide-level filtering:**
   - **Log Ratio:** Uses all peptides mapping to a protein
   - **Unshared Log Ratio:** Uses unique peptides only (removes shared peptides that could cause ambiguity)

### Quality Metrics

**Expected characteristics:**
- **Mean log ratio ≈ 0** for most proteins (centered on pooled reference)
- **Median similar to mean** (symmetric distribution)
- **Standard deviation** reflects biological + technical variability
- **Missing values:** Some proteins may not be quantified in all samples

---

## Relevance to DGAT Immunology Project

### Primary Uses

1. **DGAT Protein Expression:**
   - Validate RNA-seq findings at protein level
   - Check if DGAT1/DGAT2 mRNA correlates with protein abundance
   - Identify post-transcriptional regulation

2. **Multi-Omic Integration:**
   - Combine with TCGA RNA-seq data
   - Identify transcript-protein discordances
   - Validate immune signatures at protein level

3. **Immune Protein Markers:**
   - Quantify immune checkpoint proteins (PD-L1, CTLA-4, LAG3, etc.)
   - Measure cytokines and chemokines
   - Assess infiltrating immune cell markers

4. **Biomarker Discovery:**
   - Identify protein biomarkers associated with DGAT expression
   - Protein-level survival predictors
   - Therapeutic target validation

---

## Data Access & Citation

### Source
- **Database:** [CPTAC Data Portal](https://proteomics.cancer.gov/programs/cptac)
- **Study:** CPTAC3 GBM Discovery Cohort
- **Release Date:** December 2019
- **Revision:** r1 (first release)

### Related Publications

**CPTAC GBM Study:**
- Wang, L.B., et al. (2021). "Proteogenomic and metabolomic characterization of human glioblastoma." *Cancer Cell*, 39(4), 509-528.e20.
  - DOI: [10.1016/j.ccell.2021.01.006](https://doi.org/10.1016/j.ccell.2021.01.006)
  - PMID: [33577785](https://pubmed.ncbi.nlm.nih.gov/33577785/)

**CPTAC Program:**
- Gillette, M.A., et al. (2020). "Proteogenomic Characterization Reveals Therapeutic Vulnerabilities in Lung Adenocarcinoma." *Cell*, 182(1), 200-225.e35.
  - DOI: [10.1016/j.cell.2020.06.013](https://doi.org/10.1016/j.cell.2020.06.013)

### Citation

If using this data, please cite:

```
Wang LB, et al. (2021). Proteogenomic and metabolomic characterization of 
human glioblastoma. Cancer Cell, 39(4):509-528.e20.

Clinical Proteomic Tumor Analysis Consortium (CPTAC). CPTAC3 Glioblastoma 
Discovery Study. National Cancer Institute, December 2019.
```

---

## Processing Recommendations

### Loading the Data (R)

```r
# Load proteomics data
library(data.table)
library(dplyr)

# Read protein abundance data
proteome <- fread("CPTAC3_Glioblastoma_Multiforme_Proteome.tmt11.tsv")

# Read sample mapping
sample_map <- read.delim("CPTAC3_Glioblastoma_Multiforme_Proteome.sample.txt")

# Extract protein matrix (log2 ratios only, exclude unshared)
# Columns with "Log Ratio" but not "Unshared"
log_ratio_cols <- grep("Log Ratio$", colnames(proteome), value = TRUE)
log_ratio_cols <- log_ratio_cols[!grepl("Unshared", log_ratio_cols)]

# Create clean matrix
protein_matrix <- proteome[4:nrow(proteome), c("Gene", log_ratio_cols), with = FALSE]
rownames(protein_matrix) <- protein_matrix$Gene
protein_matrix$Gene <- NULL

# Convert to numeric
protein_matrix <- as.data.frame(protein_matrix)
for (col in colnames(protein_matrix)) {
  protein_matrix[[col]] <- as.numeric(protein_matrix[[col]])
}

# Clean column names (extract sample IDs)
colnames(protein_matrix) <- gsub(" Log Ratio$", "", colnames(protein_matrix))

# Now you have a proteins × samples matrix ready for analysis
```

### Filtering & Quality Control

```r
# Remove proteins with >50% missing values
missing_prop <- rowSums(is.na(protein_matrix)) / ncol(protein_matrix)
protein_matrix_filtered <- protein_matrix[missing_prop < 0.5, ]

# Check for DGAT genes
dgat_proteins <- protein_matrix_filtered[c("DGAT1", "DGAT2"), ]
print(dgat_proteins)

# Summary statistics
summary_stats <- data.frame(
  Mean = rowMeans(protein_matrix_filtered, na.rm = TRUE),
  Median = apply(protein_matrix_filtered, 1, median, na.rm = TRUE),
  SD = apply(protein_matrix_filtered, 1, sd, na.rm = TRUE),
  Missing = rowSums(is.na(protein_matrix_filtered))
)
```

### Batch Effect Correction

```r
# Extract batch information from sample names
# Sample names follow pattern: CPT0204410003
# Map to batches using sample_map file

# Use ComBat (sva package) for batch correction
library(sva)

# Create batch vector (requires manual mapping from sample_map)
# batch <- ... # Extract from sample_map

# Run ComBat
# proteome_corrected <- ComBat(
#   dat = as.matrix(protein_matrix_filtered),
#   batch = batch,
#   mod = model.matrix(~1, data = data.frame(row.names = colnames(protein_matrix_filtered)))
# )
```

---

## Integration with Other Data

### Matching with TCGA-GBM RNA-seq

**Sample overlap:**
- CPTAC samples are a **subset** of TCGA-GBM cohort
- Sample IDs (CPT...) can be matched to TCGA barcodes using CPTAC metadata
- Typically ~80-100 samples have both RNA-seq and proteomics data

**Integration strategy:**
```r
# 1. Match sample IDs between CPTAC and TCGA
# 2. Extract common samples from both datasets
# 3. Correlation analysis: mRNA vs. protein for each gene
# 4. Identify transcript-protein discordances
```

### Complementary Datasets in This Project

| Dataset | Location | Type | Samples | Purpose |
|---------|----------|------|---------|---------|
| TCGA-GBM RNA-seq | `Raw_Data/TCGA_GBM/` | Transcriptomics | ~170 | Primary analysis |
| CGGA-GBM RNA-seq | `Raw_Data/CGGA_GBM/` | Transcriptomics | ~100 | Validation cohort |
| CPTAC3-GBM Proteomics | `Raw_Data/HPA_Protein/` | Proteomics | ~100 | Protein validation |
| GTEx Brain | `Raw_Data/GTEx_Brain/` | Transcriptomics | Normal | Baseline comparison |

---

## Known Issues & Limitations

### Missing Values
- Not all proteins quantified in all samples
- Low-abundance proteins may have >50% missing values
- Consider imputation for downstream analysis

### Batch Effects
- 12 TMT batches may introduce technical variation
- Pooled reference mitigates but doesn't eliminate batch effects
- Recommend ComBat or similar batch correction

### Sample Size
- Discovery cohort (~100 samples) is smaller than TCGA RNA-seq cohort (~170 samples)
- May limit statistical power for rare proteins or small effect sizes

### Isoform Ambiguity
- Protein-level quantification (not isoform-specific)
- Multiple isoforms of same gene may be aggregated
- "Unshared Log Ratio" reduces ambiguity but may lose data

---

## Suggested Analyses

### 1. DGAT Protein Expression Analysis
```r
# Extract DGAT1 and DGAT2 protein abundance
dgat1_protein <- protein_matrix_filtered["DGAT1", ]
dgat2_protein <- protein_matrix_filtered["DGAT2", ]

# Compare with RNA-seq data (if samples matched)
# cor(dgat1_rna, dgat1_protein, use = "complete.obs")
```

### 2. Immune Protein Markers
```r
# Define immune checkpoint proteins
immune_checkpoints <- c("CD274", "PDCD1", "CTLA4", "LAG3", "HAVCR2", "TIGIT")

# Extract their protein abundance
immune_protein <- protein_matrix_filtered[immune_checkpoints, ]
immune_protein <- immune_protein[complete.cases(immune_protein), ]

# Correlation with DGAT proteins
# cor(t(dgat_proteins), t(immune_protein), use = "pairwise.complete.obs")
```

### 3. Multi-Omic Integration
```r
# For matched samples:
# 1. mRNA-protein correlation for all genes
# 2. Identify genes with high mRNA but low protein (translational repression)
# 3. Identify genes with low mRNA but high protein (protein stability)
# 4. Focus on DGAT-correlated immune genes
```

---

## File Provenance

**Download Date:** Unknown (present in project as of 2025-10-10)  
**Source:** CPTAC Data Portal  
**Version:** December 2019, Revision 1  
**Processing:** Unprocessed raw data (as downloaded)

---

## Contact & Support

For questions about:
- **CPTAC data:** [CPTAC Support](https://proteomics.cancer.gov/contact)
- **This project:** See main project README

---

**Last Updated:** 2025-10-10  
**Created By:** DGAT Immunology Project Analysis Pipeline

