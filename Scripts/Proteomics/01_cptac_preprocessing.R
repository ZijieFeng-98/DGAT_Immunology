#!/usr/bin/env Rscript
################################################################################
# CPTAC GBM Proteomics Data Cleanup Pipeline
# Purpose: Process CPTAC3 GBM Discovery Cohort proteomics data for DGAT 
#          cancer immunology analysis
# Reference: Wang et al. (2021) Cancer Cell 39(4):509-528
# Date: 2025-10-10
################################################################################

# Load required libraries
suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(readxl)
  library(ggplot2)
  library(pheatmap)
  library(RColorBrewer)
})

# Set working directory and paths
BASE_DIR <- "/Users/zijiefeng/Desktop/Guo's lab/My_Research/DGAT_Immunology"
RAW_DATA_DIR <- file.path(BASE_DIR, "Raw_Data/HPA_Protein")
OUTPUT_DIR <- file.path(BASE_DIR, "Processed_Data/CPTAC_GBM_Proteomics")

# Create output directories
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(OUTPUT_DIR, "QC_Plots"), recursive = TRUE, showWarnings = FALSE)

################################################################################
# STEP 1: LOAD RAW DATA
################################################################################

cat("========================================================================\n")
cat("CPTAC GBM PROTEOMICS DATA PROCESSING PIPELINE\n")
cat("========================================================================\n\n")

cat("Step 1: Loading CPTAC data files...\n")

# Load protein abundance data (TMT11-plex)
proteome_file <- file.path(RAW_DATA_DIR, "CPTAC3_Glioblastoma_Multiforme_Proteome.tmt11.tsv")
proteome <- fread(proteome_file)
cat("  - Proteome data loaded:", nrow(proteome), "rows x", ncol(proteome), "columns\n")

# Load sample mapping (CRITICAL for ID matching)
mapping_file <- file.path(RAW_DATA_DIR, "S048_CPTAC_GBM_Discovery_Cohort_TMT11_CaseID_SampleID_AliquotID_Map_Dec2019_r1.xlsx")
mapping <- read_excel(mapping_file, skip = 6)
cat("  - Sample mapping loaded:", nrow(mapping), "samples\n")

# Load clinical data (use Clinical_Attributes sheet)
clinical_file <- file.path(RAW_DATA_DIR, "S048_CPTAC_GBM_Discovery_Cohort_Clinical_Data_Dec2019_r1.xlsx")
clinical <- read_excel(clinical_file, sheet = "Clinical_Attributes")
cat("  - Clinical data loaded:", nrow(clinical), "patients\n\n")

################################################################################
# STEP 2: EXTRACT AND CLEAN PROTEIN MATRIX
################################################################################

cat("Step 2: Extracting protein matrix...\n")

# Identify Log Ratio columns (NOT Unshared Log Ratio)
log_ratio_cols <- grep("Log Ratio$", colnames(proteome), value = TRUE)
log_ratio_cols <- log_ratio_cols[!grepl("Unshared", log_ratio_cols)]
cat("  - Found", length(log_ratio_cols), "protein abundance columns\n")

# Remove summary statistics rows (rows 1-3: Mean, Median, StdDev)
proteome_data <- proteome[4:nrow(proteome), ]
cat("  - Removed summary rows, kept", nrow(proteome_data), "proteins\n")

# Create protein matrix
protein_matrix <- proteome_data[, c("Gene", log_ratio_cols), with = FALSE]

# Convert to data.frame FIRST (before setting rownames)
protein_matrix <- as.data.frame(protein_matrix)

# Then set rownames
rownames(protein_matrix) <- protein_matrix$Gene
protein_matrix$Gene <- NULL

# Convert to numeric matrix
for(col in colnames(protein_matrix)) {
  protein_matrix[[col]] <- as.numeric(protein_matrix[[col]])
}

cat("  - Protein matrix dimensions:", nrow(protein_matrix), "proteins x", 
    ncol(protein_matrix), "samples\n\n")

################################################################################
# STEP 3: SAMPLE ID MAPPING (Aliquot ID → Case ID)
################################################################################

cat("Step 3: Mapping sample IDs...\n")

# Extract aliquot IDs from column names
aliquot_ids <- gsub(" Log Ratio$", "", colnames(protein_matrix))

# Create mapping lookup (using correct column names)
aliquot_to_case <- setNames(mapping$`Case ID (Participant ID)`, mapping$`Aliquot ID`)

# Map to case IDs
case_ids <- aliquot_to_case[aliquot_ids]

# Check for unmapped samples
n_unmapped <- sum(is.na(case_ids))
if(n_unmapped > 0) {
  cat("  WARNING:", n_unmapped, "samples could not be mapped to Case IDs\n")
}

cat("  - Successfully mapped", sum(!is.na(case_ids)), "samples\n\n")

################################################################################
# STEP 4: SAMPLE FILTERING
################################################################################

cat("Step 4: Filtering samples...\n")

# Identify samples to keep
valid_samples <- !is.na(case_ids) & 
                 !grepl("POOL|QC|Ref|Reference", aliquot_ids, ignore.case = TRUE)

protein_matrix_filtered <- protein_matrix[, valid_samples]
case_ids_filtered <- case_ids[valid_samples]
aliquot_ids_filtered <- aliquot_ids[valid_samples]

cat("  - Removed", sum(!valid_samples), "QC/pooled/unmapped samples\n")
cat("  - Retained", ncol(protein_matrix_filtered), "samples\n")

# Deduplicate patients (keep one sample per patient)
# Strategy: keep first occurrence
unique_patients <- !duplicated(case_ids_filtered)
protein_matrix_unique <- protein_matrix_filtered[, unique_patients]
case_ids_unique <- case_ids_filtered[unique_patients]

cat("  - After deduplication:", ncol(protein_matrix_unique), "unique patients\n\n")

################################################################################
# STEP 5: PROTEIN-LEVEL FILTERING
################################################################################

cat("Step 5: Filtering proteins...\n")

# Calculate missing value proportion per protein
missing_prop <- rowSums(is.na(protein_matrix_unique)) / ncol(protein_matrix_unique)

# Filter proteins with >50% missing values (Wang et al. 2021 criterion)
protein_filtered <- protein_matrix_unique[missing_prop <= 0.5, ]

cat("  - Removed", sum(missing_prop > 0.5), "proteins with >50% missing values\n")
cat("  - Retained", nrow(protein_filtered), "proteins\n\n")

################################################################################
# STEP 6: CHECK KEY PROTEINS
################################################################################

cat("Step 6: Checking key proteins for DGAT immunology analysis...\n")

# DGAT proteins
dgat_genes <- c("DGAT1", "DGAT2")
dgat_present <- dgat_genes[dgat_genes %in% rownames(protein_filtered)]
cat("  - DGAT proteins detected:", paste(dgat_present, collapse = ", "), "\n")

# Immune checkpoint proteins (based on troubleshooting findings)
immune_markers <- c("CD274",   # PD-L1
                    "HAVCR2",  # TIM-3
                    "CD8A",    # CD8
                    "MRC1",    # CD206
                    "ARG1")    # Arginase-1
immune_present <- immune_markers[immune_markers %in% rownames(protein_filtered)]
cat("  - Immune markers detected:", length(immune_present), "of", 
    length(immune_markers), "\n")
cat("    ", paste(immune_present, collapse = ", "), "\n")

# Lipid metabolism proteins (based on troubleshooting findings)
lipid_markers <- c("DGAT1", "PLIN2", "PLIN3", "CPT1A", 
                   "FASN", "ACACA", "SCD", "SREBF1")
lipid_present <- lipid_markers[lipid_markers %in% rownames(protein_filtered)]
cat("  - Lipid markers detected:", length(lipid_present), "of", 
    length(lipid_markers), "\n")
cat("    ", paste(lipid_present, collapse = ", "), "\n\n")

################################################################################
# STEP 7: QUALITY CONTROL PLOTS
################################################################################

cat("Step 7: Generating QC plots...\n")

# 7.1 Missing value distribution
pdf(file.path(OUTPUT_DIR, "QC_Plots/01_missing_values_distribution.pdf"), width = 10, height = 8)
missing_matrix <- is.na(protein_filtered)
protein_missing <- rowSums(missing_matrix) / ncol(missing_matrix) * 100
sample_missing <- colSums(missing_matrix) / nrow(missing_matrix) * 100

par(mfrow = c(1, 2))
hist(protein_missing, breaks = 50, 
     main = "Missing Values per Protein", 
     xlab = "% Missing", col = "steelblue")
hist(sample_missing, breaks = 30,
     main = "Missing Values per Sample",
     xlab = "% Missing", col = "coral")
dev.off()
cat("  - Saved: 01_missing_values_distribution.pdf\n")

# 7.2 Sample-sample correlation
pdf(file.path(OUTPUT_DIR, "QC_Plots/02_sample_correlations.pdf"), width = 8, height = 6)
cor_matrix <- cor(protein_filtered, use = "pairwise.complete.obs")
hist(cor_matrix[lower.tri(cor_matrix)], breaks = 50,
     main = "Sample-Sample Correlations",
     xlab = "Pearson r", col = "darkgreen", xlim = c(0, 1))
abline(v = median(cor_matrix[lower.tri(cor_matrix)]), col = "red", lwd = 2, lty = 2)
legend("topleft", legend = paste("Median r =", 
       round(median(cor_matrix[lower.tri(cor_matrix)]), 3)),
       col = "red", lty = 2, lwd = 2)
dev.off()
cat("  - Saved: 02_sample_correlations.pdf\n")

# 7.3 PCA plot
pdf(file.path(OUTPUT_DIR, "QC_Plots/03_pca_plot.pdf"), width = 10, height = 8)
pca_data <- prcomp(t(na.omit(protein_filtered)), scale. = TRUE)
pca_var <- round(summary(pca_data)$importance[2, 1:2] * 100, 1)

plot(pca_data$x[, 1], pca_data$x[, 2], 
     pch = 19, col = "steelblue", cex = 2,
     xlab = paste0("PC1 (", pca_var[1], "% variance)"),
     ylab = paste0("PC2 (", pca_var[2], "% variance)"),
     main = "PCA: CPTAC GBM Proteomics Data")
text(pca_data$x[, 1], pca_data$x[, 2], 
     labels = 1:nrow(pca_data$x), pos = 3, cex = 0.6)
dev.off()
cat("  - Saved: 03_pca_plot.pdf\n")

# 7.4 DGAT1 distribution
if("DGAT1" %in% rownames(protein_filtered)) {
  pdf(file.path(OUTPUT_DIR, "QC_Plots/04_dgat1_distribution.pdf"), width = 8, height = 6)
  dgat1_values <- as.numeric(protein_filtered["DGAT1", ])
  hist(dgat1_values, breaks = 30, 
       main = "DGAT1 Protein Expression Distribution",
       xlab = "Log2 Ratio (relative to pooled reference)",
       col = "purple", xlim = c(min(dgat1_values, na.rm = TRUE) - 0.5,
                                 max(dgat1_values, na.rm = TRUE) + 0.5))
  abline(v = median(dgat1_values, na.rm = TRUE), col = "red", lwd = 2, lty = 2)
  legend("topright", legend = paste("Median =", 
         round(median(dgat1_values, na.rm = TRUE), 3)),
         col = "red", lty = 2, lwd = 2)
  dev.off()
  cat("  - Saved: 04_dgat1_distribution.pdf\n")
}

# 7.5 Correlation heatmap (DGAT + key immune/lipid markers)
key_proteins <- c("DGAT1", "CD274", "HAVCR2", "CD8A", "MRC1", 
                  "ARG1", "PLIN2", "PLIN3", "CPT1A", "FASN", "ACACA")
key_proteins_present <- key_proteins[key_proteins %in% rownames(protein_filtered)]

if(length(key_proteins_present) > 2) {
  pdf(file.path(OUTPUT_DIR, "QC_Plots/05_dgat_immune_correlation.pdf"), width = 8, height = 7)
  key_matrix <- protein_filtered[key_proteins_present, ]
  cor_key <- cor(t(key_matrix), use = "pairwise.complete.obs")
  pheatmap(cor_key, 
           main = "DGAT-Immune Protein Correlations",
           color = colorRampPalette(rev(brewer.pal(11, "RdBu")))(100),
           breaks = seq(-1, 1, length.out = 101),
           display_numbers = TRUE,
           number_format = "%.2f",
           fontsize_number = 8)
  dev.off()
  cat("  - Saved: 05_dgat_immune_correlation.pdf\n")
}

cat("\n")

################################################################################
# STEP 8: MATCH WITH CLINICAL DATA
################################################################################

cat("Step 8: Matching with clinical data...\n")

# Rename columns to Case IDs
colnames(protein_filtered) <- case_ids_unique

# Match clinical data
# Note: Clinical data uses "case_id" column (from Clinical_Attributes sheet)
case_col <- "case_id"
if(case_col %in% colnames(clinical)) {
  clinical_matched <- clinical %>%
    filter(!!sym(case_col) %in% case_ids_unique) %>%
    arrange(match(!!sym(case_col), case_ids_unique))
} else {
  # Fallback: try to find column with case IDs
  case_col <- colnames(clinical)[1]
  clinical_matched <- clinical %>%
    filter(!!sym(case_col) %in% case_ids_unique) %>%
    arrange(match(!!sym(case_col), case_ids_unique))
}

# Verify match
if(nrow(clinical_matched) != ncol(protein_filtered)) {
  cat("  WARNING: Clinical data doesn't match all samples!\n")
  cat("  Protein samples:", ncol(protein_filtered), "\n")
  cat("  Clinical samples:", nrow(clinical_matched), "\n")
} else {
  cat("  - Successfully matched", nrow(clinical_matched), "samples\n")
}

################################################################################
# STEP 9: SAVE CLEANED DATA
################################################################################

cat("\nStep 9: Saving cleaned data...\n")

# Save protein matrix
saveRDS(protein_filtered, file.path(OUTPUT_DIR, "protein_matrix_cleaned.rds"))
fwrite(as.data.table(protein_filtered, keep.rownames = "Gene"), 
       file.path(OUTPUT_DIR, "protein_matrix_cleaned.csv"))
cat("  - Saved: protein_matrix_cleaned.rds/.csv\n")

# Save clinical data
saveRDS(clinical_matched, file.path(OUTPUT_DIR, "clinical_data_matched.rds"))
fwrite(clinical_matched, file.path(OUTPUT_DIR, "clinical_data_matched.csv"))
cat("  - Saved: clinical_data_matched.rds/.csv\n")

# Save sample metadata
sample_metadata <- data.frame(
  Aliquot_ID = aliquot_ids_filtered[unique_patients],
  Case_ID = case_ids_unique,
  In_Clinical = case_ids_unique %in% clinical_matched[[case_col]]
)
fwrite(sample_metadata, file.path(OUTPUT_DIR, "sample_metadata.csv"))
cat("  - Saved: sample_metadata.csv\n")

# Save QC summary
qc_summary <- data.frame(
  Metric = c(
    "Raw proteins",
    "Raw samples", 
    "Proteins after missing filter",
    "Unique patients",
    "DGAT1 detected",
    "DGAT2 detected",
    "Immune markers detected",
    "Median sample correlation"
  ),
  Value = c(
    nrow(protein_matrix),
    ncol(protein_matrix),
    nrow(protein_filtered),
    ncol(protein_filtered),
    "DGAT1" %in% rownames(protein_filtered),
    "DGAT2" %in% rownames(protein_filtered),
    length(immune_present),
    round(median(cor_matrix[lower.tri(cor_matrix)]), 3)
  )
)
fwrite(qc_summary, file.path(OUTPUT_DIR, "QC_Summary.csv"))
cat("  - Saved: QC_Summary.csv\n")

# Save processing log
log_content <- c(
  "CPTAC GBM Proteomics Processing Log",
  paste("Date:", Sys.time()),
  paste("R version:", R.version.string),
  "",
  "Input Files:",
  paste("  - Proteome:", proteome_file),
  paste("  - Mapping:", mapping_file),
  paste("  - Clinical:", clinical_file),
  "",
  "Processing Steps:",
  paste("  1. Loaded", nrow(proteome), "x", ncol(proteome), "raw proteome data"),
  paste("  2. Extracted", length(log_ratio_cols), "Log Ratio columns"),
  paste("  3. Filtered to", ncol(protein_matrix_filtered), "valid samples"),
  paste("  4. Deduplicated to", ncol(protein_matrix_unique), "unique patients"),
  paste("  5. Filtered proteins (>50% missing):", nrow(protein_filtered), "retained"),
  "",
  "Output Files:",
  paste("  - Protein matrix:", nrow(protein_filtered), "proteins x", ncol(protein_filtered), "samples"),
  paste("  - Clinical data:", nrow(clinical_matched), "samples"),
  paste("  - QC plots: 5 files"),
  "",
  "Key Proteins Detected:",
  paste("  - DGAT1:", "DGAT1" %in% rownames(protein_filtered)),
  paste("  - DGAT2:", "DGAT2" %in% rownames(protein_filtered)),
  paste("  - Immune markers:", length(immune_present), "of", length(immune_markers)),
  paste("  - Lipid markers:", length(lipid_present), "of", length(lipid_markers))
)
writeLines(log_content, file.path(OUTPUT_DIR, "Processing_Log.txt"))
cat("  - Saved: Processing_Log.txt\n")

################################################################################
# FINAL SUMMARY
################################################################################

cat("\n")
cat("========================================================================\n")
cat("CPTAC GBM PROTEOMICS DATA PROCESSING COMPLETE\n")
cat("========================================================================\n")
cat("Final dataset dimensions:\n")
cat("  - Proteins:", nrow(protein_filtered), "\n")
cat("  - Patients:", ncol(protein_filtered), "\n")
cat("  - DGAT1 present:", "DGAT1" %in% rownames(protein_filtered), "\n")
cat("  - DGAT2 present:", "DGAT2" %in% rownames(protein_filtered), "\n")
cat("  - Clinical data matched:", nrow(clinical_matched), "patients\n")
cat("\nOutput directory:", OUTPUT_DIR, "\n")
cat("  - Processed data: protein_matrix_cleaned.rds/.csv\n")
cat("  - Clinical data: clinical_data_matched.rds/.csv\n")
cat("  - QC plots: QC_Plots/ (5 plots)\n")
cat("  - QC summary: QC_Summary.csv\n")
cat("\nNext steps:\n")
cat("  1. Review QC plots for any issues\n")
cat("  2. Run survival analysis (if clinical data includes OS/PFS)\n")
cat("  3. Correlate DGAT1/2 with immune markers\n")
cat("  4. Compare with TCGA RNA-seq data\n")
cat("========================================================================\n")

