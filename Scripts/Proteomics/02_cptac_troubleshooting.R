#!/usr/bin/env Rscript
################################################################################
# CPTAC GBM Proteomics - Troubleshooting & Gene Discovery Script
# Purpose: Diagnose missing DGAT1/2 and immune markers in proteomics data
# Date: 2025-10-10
################################################################################

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(readxl)
})

# Set working directory and paths
BASE_DIR <- "/Users/zijiefeng/Desktop/Guo's lab/My_Research/DGAT_Immunology"
RAW_DATA_DIR <- file.path(BASE_DIR, "Raw_Data/HPA_Protein")
OUTPUT_DIR <- file.path(BASE_DIR, "Processed_Data/CPTAC_GBM_Proteomics")

# Create output directory for troubleshooting
dir.create(file.path(OUTPUT_DIR, "Troubleshooting"), showWarnings = FALSE)

################################################################################
# PART 1: SEARCH FOR DGAT GENES IN RAW DATA
################################################################################

cat(paste(rep("=", 78), collapse=""), "\n")
cat("PART 1: SEARCHING FOR DGAT GENES\n")
cat(paste(rep("=", 78), collapse=""), "\n\n")

# Load raw proteomics data
proteome_raw <- fread(file.path(RAW_DATA_DIR, "CPTAC3_Glioblastoma_Multiforme_Proteome.tmt11.tsv"))

cat("Raw data dimensions:", nrow(proteome_raw), "x", ncol(proteome_raw), "\n\n")

# Search for DGAT with multiple strategies
cat("Strategy 1: Exact match search\n")
dgat_exact <- proteome_raw$Gene[grep("^DGAT[12]$", proteome_raw$Gene, ignore.case = FALSE)]
cat("  Exact matches (DGAT1, DGAT2):", paste(dgat_exact, collapse = ", "), 
    ifelse(length(dgat_exact) == 0, "(NONE FOUND)", ""), "\n\n")

cat("Strategy 2: Case-insensitive search\n")
dgat_case <- proteome_raw$Gene[grep("dgat", proteome_raw$Gene, ignore.case = TRUE)]
cat("  Case-insensitive matches:", paste(dgat_case, collapse = ", "), 
    ifelse(length(dgat_case) == 0, "(NONE FOUND)", ""), "\n\n")

cat("Strategy 3: Partial match search (contains 'DGAT')\n")
dgat_partial <- proteome_raw$Gene[grep("DGAT", proteome_raw$Gene, ignore.case = TRUE)]
cat("  Partial matches:", paste(dgat_partial, collapse = ", "), 
    ifelse(length(dgat_partial) == 0, "(NONE FOUND)", ""), "\n\n")

cat("Strategy 4: Alternative names search\n")
dgat_alt_names <- c("DGAT1", "DGAT2", "DGAT-1", "DGAT-2", 
                    "ARAT", "DGAT1L", "DGAT2L",
                    "diacylglycerol acyltransferase")
dgat_alt <- proteome_raw$Gene[proteome_raw$Gene %in% dgat_alt_names]
cat("  Alternative name matches:", paste(dgat_alt, collapse = ", "), 
    ifelse(length(dgat_alt) == 0, "(NONE FOUND)", ""), "\n\n")

cat("Strategy 5: Fuzzy matching (approximate grep)\n")
all_genes <- proteome_raw$Gene[4:nrow(proteome_raw)]  # Skip summary rows
dgat1_fuzzy <- agrep("DGAT1", all_genes, max.distance = 2, value = TRUE)
dgat2_fuzzy <- agrep("DGAT2", all_genes, max.distance = 2, value = TRUE)
cat("  DGAT1 fuzzy matches:", paste(dgat1_fuzzy, collapse = ", "), 
    ifelse(length(dgat1_fuzzy) == 0, "(NONE FOUND)", ""), "\n")
cat("  DGAT2 fuzzy matches:", paste(dgat2_fuzzy, collapse = ", "), 
    ifelse(length(dgat2_fuzzy) == 0, "(NONE FOUND)", ""), "\n\n")

# Check if DGAT genes have high missing values
cat("Strategy 6: Check if DGAT genes were filtered due to missing values\n")
if(length(dgat_exact) > 0) {
  for(gene in dgat_exact) {
    gene_row <- proteome_raw[proteome_raw$Gene == gene, ]
    log_ratio_cols <- grep("Log Ratio$", colnames(gene_row), value = TRUE)
    log_ratio_cols <- log_ratio_cols[!grepl("Unshared", log_ratio_cols)]
    
    values <- as.numeric(gene_row[, ..log_ratio_cols])
    missing_pct <- sum(is.na(values)) / length(values) * 100
    
    cat("  ", gene, "- Missing values:", round(missing_pct, 1), "%\n")
    if(missing_pct > 50) {
      cat("    ⚠️ WARNING: >50% missing - would be filtered!\n")
    }
  }
} else {
  cat("  Cannot check - no DGAT genes found\n")
}

cat("\n")

################################################################################
# PART 2: SEARCH FOR IMMUNE MARKERS
################################################################################

cat(paste(rep("=", 78), collapse=""), "\n")
cat("PART 2: SEARCHING FOR IMMUNE MARKERS\n")
cat(paste(rep("=", 78), collapse=""), "\n\n")

# Define immune markers with alternative names
immune_markers_list <- list(
  PD_L1 = c("CD274", "PD-L1", "PDL1", "PDCD1L1", "PDCD1LG1"),
  PD_1 = c("PDCD1", "PD-1", "PD1", "CD279"),
  CTLA4 = c("CTLA4", "CTLA-4", "CD152"),
  LAG3 = c("LAG3", "LAG-3", "CD223"),
  TIM3 = c("HAVCR2", "TIM3", "TIM-3", "TIMD3"),
  CD8 = c("CD8A", "CD8", "CD8a"),
  CD206 = c("MRC1", "CD206", "MRC1L1"),
  ARG1 = c("ARG1", "arginase-1", "arginase 1")
)

# Search for each marker
immune_found <- list()
for(marker_name in names(immune_markers_list)) {
  aliases <- immune_markers_list[[marker_name]]
  pattern <- paste(aliases, collapse = "|")
  
  found <- proteome_raw$Gene[grep(pattern, proteome_raw$Gene, ignore.case = TRUE)]
  
  cat(marker_name, "(", paste(aliases[1:2], collapse = "/"), "):\n")
  if(length(found) > 0) {
    cat("  ✓ FOUND:", paste(found, collapse = ", "), "\n")
    immune_found[[marker_name]] <- found
  } else {
    cat("  ✗ NOT FOUND\n")
  }
}

cat("\nSummary: Found", length(immune_found), "out of", 
    length(immune_markers_list), "immune markers\n\n")

################################################################################
# PART 3: SEARCH FOR LIPID METABOLISM MARKERS
################################################################################

cat(paste(rep("=", 78), collapse=""), "\n")
cat("PART 3: SEARCHING FOR LIPID METABOLISM MARKERS\n")
cat(paste(rep("=", 78), collapse=""), "\n\n")

lipid_markers <- c("PLIN2", "PLIN3", "TIP47", "CPT1A", "CPT1B", 
                   "FASN", "ACACA", "SCD", "SREBF1", "PPARG")

lipid_found <- c()
for(marker in lipid_markers) {
  found <- proteome_raw$Gene[grep(marker, proteome_raw$Gene, ignore.case = TRUE)]
  
  if(length(found) > 0) {
    cat("  ✓", marker, "- Found as:", paste(found, collapse = ", "), "\n")
    lipid_found <- c(lipid_found, found)
  } else {
    cat("  ✗", marker, "- NOT FOUND\n")
  }
}

cat("\nSummary: Found", length(lipid_found), "lipid metabolism markers\n\n")

################################################################################
# PART 4: CHECK CLINICAL DATA STRUCTURE
################################################################################

cat(paste(rep("=", 78), collapse=""), "\n")
cat("PART 4: DIAGNOSING CLINICAL DATA MATCHING ISSUE\n")
cat(paste(rep("=", 78), collapse=""), "\n\n")

# Try loading clinical file with different parameters
clinical_file <- file.path(RAW_DATA_DIR, "S048_CPTAC_GBM_Discovery_Cohort_Clinical_Data_Dec2019_r1.xlsx")

cat("Attempt 1: Load without skip\n")
tryCatch({
  clinical1 <- read_excel(clinical_file)
  cat("  Dimensions:", nrow(clinical1), "rows x", ncol(clinical1), "columns\n")
  cat("  Column names:", paste(colnames(clinical1)[1:min(5, ncol(clinical1))], 
                                collapse = ", "), "...\n")
  cat("  First column name:", colnames(clinical1)[1], "\n")
  cat("  Sample first 3 IDs:", paste(head(clinical1[[1]], 3), collapse = ", "), "\n\n")
}, error = function(e) {
  cat("  ERROR:", e$message, "\n\n")
})

cat("Attempt 2: Load with skip=1\n")
tryCatch({
  clinical2 <- read_excel(clinical_file, skip = 1)
  cat("  Dimensions:", nrow(clinical2), "rows x", ncol(clinical2), "columns\n")
  cat("  Column names:", paste(colnames(clinical2)[1:min(5, ncol(clinical2))], 
                                collapse = ", "), "...\n")
  cat("  First column name:", colnames(clinical2)[1], "\n")
  cat("  Sample first 3 IDs:", paste(head(clinical2[[1]], 3), collapse = ", "), "\n\n")
}, error = function(e) {
  cat("  ERROR:", e$message, "\n\n")
})

cat("Attempt 3: Check Excel sheet names\n")
tryCatch({
  sheets <- excel_sheets(clinical_file)
  cat("  Available sheets:", paste(sheets, collapse = ", "), "\n\n")
  
  # Try loading each sheet
  for(sheet in sheets) {
    cat("  Sheet:", sheet, "\n")
    temp <- read_excel(clinical_file, sheet = sheet)
    cat("    Dimensions:", nrow(temp), "x", ncol(temp), "\n")
    cat("    First column:", colnames(temp)[1], "\n")
  }
}, error = function(e) {
  cat("  ERROR:", e$message, "\n\n")
})

# Load sample metadata to check Case IDs
cat("\nProtein sample Case IDs:\n")
sample_meta <- fread(file.path(OUTPUT_DIR, "sample_metadata.csv"))
cat("  Total samples:", nrow(sample_meta), "\n")
cat("  Example Case IDs:", paste(head(sample_meta$Case_ID, 5), collapse = ", "), "\n")
cat("  Case ID format:", class(sample_meta$Case_ID), "\n\n")

################################################################################
# PART 5: GENERATE SOLUTIONS & RECOMMENDATIONS
################################################################################

cat(paste(rep("=", 78), collapse=""), "\n")
cat("PART 5: RECOMMENDATIONS & NEXT STEPS\n")
cat(paste(rep("=", 78), collapse=""), "\n\n")

cat("📋 ISSUE 1: DGAT Genes Not Found\n")
if(length(dgat_exact) > 0) {
  cat("  ✓ Solution: DGAT genes ARE present in raw data\n")
  cat("  → Check if filtered due to missing values (see Strategy 6 above)\n")
  cat("  → If >50% missing, consider lowering threshold to 60-70%\n\n")
  
  cat("  Code to re-run with relaxed filter:\n")
  cat("  # In cleanup script, change line:\n")
  cat("  # protein_filtered <- protein_matrix_unique[missing_prop <= 0.5, ]\n")
  cat("  # TO:\n")
  cat("  # protein_filtered <- protein_matrix_unique[missing_prop <= 0.7, ]\n\n")
  
} else {
  cat("  ✗ Problem: DGAT genes NOT in proteomics dataset\n")
  cat("  → These proteins may be below detection limit for mass spec\n")
  cat("  → Consider using RNA-seq data (TCGA) as primary analysis\n")
  cat("  → CPTAC data can validate other immune/lipid markers\n\n")
}

cat("📋 ISSUE 2: Immune Markers\n")
if(length(immune_found) > 0) {
  cat("  ✓ Solution: Some immune markers found!\n")
  cat("  → Update cleanup script with correct gene names\n")
  cat("  → Use aliases found in Part 2 above\n\n")
  
  cat("  Code to update immune marker list:\n")
  cat("  immune_markers <- c(\n")
  for(marker_name in names(immune_found)) {
    cat("    '", immune_found[[marker_name]][1], "',  # ", marker_name, "\n", sep = "")
  }
  cat("  )\n\n")
  
} else {
  cat("  ⚠️ Warning: No immune markers found with common aliases\n")
  cat("  → This is unusual - double-check raw data\n")
  cat("  → May need to examine actual gene symbols in proteome_raw$Gene\n\n")
}

cat("📋 ISSUE 3: Clinical Data Matching\n")
cat("  → Check results from Part 4 above\n")
cat("  → Likely needs skip parameter or specific sheet name\n")
cat("  → May need to rename Case ID column\n\n")

cat("  Code to fix (based on Part 4 findings):\n")
cat("  # In cleanup script, replace:\n")
cat("  # clinical <- read_excel('clinical_file.xlsx')\n")
cat("  # WITH:\n")
cat("  # clinical <- read_excel('clinical_file.xlsx', skip = 1)  # or skip = 0\n")
cat("  # OR:\n")
cat("  # clinical <- read_excel('clinical_file.xlsx', sheet = 'sheet_name')\n\n")

################################################################################
# PART 6: SAVE FINDINGS
################################################################################

cat(paste(rep("=", 78), collapse=""), "\n")
cat("PART 6: SAVING TROUBLESHOOTING RESULTS\n")
cat(paste(rep("=", 78), collapse=""), "\n\n")

# Create comprehensive findings report
findings <- list(
  dgat_genes_found = if(length(dgat_exact) > 0) dgat_exact else "NONE",
  immune_markers_found = immune_found,
  lipid_markers_found = lipid_found,
  total_proteins_in_raw = nrow(proteome_raw) - 3,  # Exclude summary rows
  dgat_search_strategies = list(
    exact = dgat_exact,
    case_insensitive = dgat_case,
    partial = dgat_partial,
    fuzzy_dgat1 = dgat1_fuzzy,
    fuzzy_dgat2 = dgat2_fuzzy
  )
)

saveRDS(findings, file.path(OUTPUT_DIR, "Troubleshooting/troubleshooting_findings.rds"))
cat("Saved detailed findings to: Troubleshooting/troubleshooting_findings.rds\n\n")

# Create summary table
summary_table <- data.frame(
  Category = c("DGAT Genes", "Immune Markers", "Lipid Markers"),
  Expected = c(2, length(immune_markers_list), length(lipid_markers)),
  Found = c(
    length(dgat_exact),
    length(immune_found),
    length(lipid_found)
  ),
  Status = c(
    ifelse(length(dgat_exact) >= 1, "✓ FOUND", "✗ MISSING"),
    ifelse(length(immune_found) >= 3, "✓ PARTIAL", "✗ LOW"),
    ifelse(length(lipid_found) >= 5, "✓ GOOD", "✓ PARTIAL")
  )
)

fwrite(summary_table, file.path(OUTPUT_DIR, "Troubleshooting/troubleshooting_summary.csv"))
cat("Saved summary table to: Troubleshooting/troubleshooting_summary.csv\n\n")

# Save gene lists for easy reference
if(length(immune_found) > 0) {
  immune_genes_df <- data.frame(
    Marker = rep(names(immune_found), sapply(immune_found, length)),
    Gene_Symbol = unlist(immune_found)
  )
  fwrite(immune_genes_df, file.path(OUTPUT_DIR, "Troubleshooting/immune_markers_found.csv"))
  cat("Saved immune marker list to: Troubleshooting/immune_markers_found.csv\n\n")
}

cat(paste(rep("=", 78), collapse=""), "\n")
cat("TROUBLESHOOTING COMPLETE\n")
cat(paste(rep("=", 78), collapse=""), "\n")
cat("\nNext steps:\n")
cat("1. Review findings above carefully\n")
cat("2. Update cleanup script with correct gene names\n")
cat("3. Re-run cleanup with adjusted parameters if needed\n")
cat("4. If DGAT not found, focus analysis on RNA-seq data\n")
cat(paste(rep("=", 78), collapse=""), "\n")

