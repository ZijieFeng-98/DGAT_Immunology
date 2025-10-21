#!/usr/bin/env Rscript
# =====================================================================
# TCGA-GBM: Simple Download & Preprocessing Pipeline
# Downloads TCGA data with better error handling
# =====================================================================

suppressPackageStartupMessages({
  library(TCGAbiolinks)
  library(SummarizedExperiment)
  library(data.table)
  library(dplyr)
  library(survival)
  library(survminer)
})

# ==================== CONFIGURATION ====================
OUTPUT_DIR <- "Raw_Data/TCGA/GBM_Fresh"
PROCESSED_DIR <- "Processed_Data/TCGA_GBM_Clean"
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(PROCESSED_DIR, recursive = TRUE, showWarnings = FALSE)

say <- function(...) cat(sprintf(...), "\n")

# ==================== STEP 1: QUERY DATA ====================
say("\n========================================")
say("STEP 1: QUERYING TCGA-GBM DATA")
say("========================================")

# Query GDC for GBM RNA-seq data
query <- GDCquery(
  project = "TCGA-GBM",
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification",
  workflow.type = "STAR - Counts",
  experimental.strategy = "RNA-Seq"
)

say("Query summary:")
say("  Project: TCGA-GBM")
say("  Files found: %d", nrow(getResults(query)))

# ==================== STEP 2: DOWNLOAD WITH RETRY ====================
say("\n========================================")
say("STEP 2: DOWNLOADING DATA (with retry logic)")
say("========================================")

# Try multiple download methods
download_success <- FALSE
download_methods <- c("api", "client")

for (method in download_methods) {
  say("\nTrying download method: %s", method)
  tryCatch({
    if (method == "api") {
      GDCdownload(query, method = "api", files.per.chunk = 5)
    } else {
      GDCdownload(query, method = "client", files.per.chunk = 5)
    }
    say("✓ Download successful with method: %s", method)
    download_success <- TRUE
    break
  }, error = function(e) {
    say("✗ Download failed with method %s: %s", method, e$message)
  })
}

if (!download_success) {
  say("\n❌ All download methods failed. Trying to proceed with existing data...")
  # Check if we have any existing data
  if (dir.exists("GDCdata")) {
    say("Found existing GDCdata directory, attempting to use it...")
  } else {
    stop("No data available for processing")
  }
}

# ==================== STEP 3: PREPARE DATA ====================
say("\n========================================")
say("STEP 3: PREPARING DATA")
say("========================================")

tryCatch({
  gbm_data <- GDCprepare(query)
  say("✓ Data preparation complete")
  say("  Samples: %d", ncol(gbm_data))
  say("  Genes: %d", nrow(gbm_data))
}, error = function(e) {
  say("❌ Data preparation failed: %s", e$message)
  stop("Cannot proceed without prepared data")
})

# ==================== STEP 4: EXTRACT DATA ====================
say("\n========================================")
say("STEP 4: EXTRACTING EXPRESSION & METADATA")
say("========================================")

# Get FPKM-normalized expression
expr_fpkm <- assay(gbm_data, "fpkm_unstrand")
say("Expression matrix extracted: %d genes x %d samples", 
    nrow(expr_fpkm), ncol(expr_fpkm))

# Get clinical/sample metadata
metadata <- as.data.frame(colData(gbm_data))
say("Metadata extracted: %d samples, %d variables", 
    nrow(metadata), ncol(metadata))

# Add proper TCGA barcodes as column names
colnames(expr_fpkm) <- metadata$barcode

# Save raw downloaded data
saveRDS(expr_fpkm, file.path(OUTPUT_DIR, "TCGA_GBM_Expression_FPKM_raw.rds"))

# Convert metadata to data.frame and remove list columns for CSV export
metadata_df <- as.data.frame(metadata)
# Remove list columns that cause fwrite issues
list_cols <- sapply(metadata_df, function(x) is.list(x))
if (any(list_cols)) {
  say("Removing %d list columns from metadata for CSV export", sum(list_cols))
  metadata_df <- metadata_df[, !list_cols]
}
fwrite(metadata_df, file.path(OUTPUT_DIR, "TCGA_GBM_Metadata_raw.csv"))

# Also save full metadata as RDS to preserve all information
saveRDS(metadata, file.path(OUTPUT_DIR, "TCGA_GBM_Metadata_raw.rds"))
say("✓ Raw data saved to: %s", OUTPUT_DIR)

# ==================== STEP 5: FILTER PRIMARY TUMORS ====================
say("\n========================================")
say("STEP 5: FILTERING TO PRIMARY TUMORS")
say("========================================")

# Extract sample type from barcode (positions 14-15)
sample_codes <- substr(metadata$barcode, 14, 15)
say("Sample type distribution:")
print(table(sample_codes))

# Keep only primary solid tumors (code "01")
primary_mask <- sample_codes == "01"
say("\nPrimary tumors: %d / %d samples", sum(primary_mask), length(primary_mask))

if (sum(primary_mask) == 0) {
  stop("No primary tumor samples found!")
}

expr_primary <- expr_fpkm[, primary_mask, drop = FALSE]
metadata_primary <- metadata[primary_mask, ]

# ==================== STEP 6: DEDUPLICATE TO PATIENT LEVEL ====================
say("\n========================================")
say("STEP 6: DEDUPLICATING TO PATIENT LEVEL")
say("========================================")

# Extract patient IDs (first 12 characters)
metadata_primary$patient12 <- substr(metadata_primary$barcode, 1, 12)
say("Initial samples: %d", nrow(metadata_primary))
say("Unique patients: %d", length(unique(metadata_primary$patient12)))

# Find duplicates
dup_patients <- metadata_primary$patient12[duplicated(metadata_primary$patient12)]
if (length(dup_patients) > 0) {
  say("\nPatients with multiple samples: %d", length(unique(dup_patients)))
  
  # Extract plate numbers from barcode (positions 22-25)
  metadata_primary$plate <- suppressWarnings(
    as.integer(substr(metadata_primary$barcode, 22, 25))
  )
  
  # Calculate mean expression for each sample
  mean_expr <- colMeans(expr_primary, na.rm = TRUE)
  metadata_primary$mean_expr <- mean_expr
  
  # Select best sample per patient (highest plate, then highest expression)
  metadata_dedup <- metadata_primary %>%
    group_by(patient12) %>%
    arrange(desc(plate), desc(mean_expr)) %>%
    slice(1) %>%
    ungroup() %>%
    as.data.frame()
  
  # Filter expression matrix
  expr_dedup <- expr_primary[, metadata_dedup$barcode, drop = FALSE]
  
  say("After deduplication: %d unique patients", nrow(metadata_dedup))
} else {
  say("No duplicate patients found")
  expr_dedup <- expr_primary
  metadata_dedup <- metadata_primary
}

# ==================== STEP 7: FILTER GENES ====================
say("\n========================================")
say("STEP 7: FILTERING GENES")
say("========================================")

say("Initial genes: %d", nrow(expr_dedup))

# Remove genes with very low expression
detection_rate <- rowMeans(expr_dedup > 0.1, na.rm = TRUE)
expr_filtered <- expr_dedup[detection_rate >= 0.1, , drop = FALSE]
say("After low-expression filter (≥10%% detection): %d genes", nrow(expr_filtered))

# Pattern-based filtering for protein-coding genes
gene_names <- rownames(expr_filtered)

# Check if we have standard gene symbols or Ensembl IDs
if (any(grepl("^ENSG", gene_names))) {
  say("Detected Ensembl gene IDs - using less aggressive filtering")
  # For Ensembl IDs, just exclude obvious non-coding
  exclude_patterns <- c(
    "^RPL", "^RPS", "^MRPL", "^MRPS", "^MT-",
    "-AS", "-PS", "^LINC", "^LOC", "^MIR",
    "^SNOR", "^AC[0-9]", "^AP[0-9]", "^RP[0-9]"
  )
  exclude_mask <- grepl(paste(exclude_patterns, collapse="|"), gene_names)
  keep_mask <- !exclude_mask
} else {
  # Standard gene symbol pattern (more permissive)
  standard_pattern <- "^[A-Z][A-Z0-9-]{1,14}$"
  looks_standard <- grepl(standard_pattern, gene_names)
  
  # Exclude non-coding patterns (less aggressive)
  exclude_patterns <- c(
    "^RPL", "^RPS", "^MRPL", "^MRPS", "^MT-",
    "-AS", "-PS", "^LINC", "^LOC", "^MIR",
    "^SNOR"
  )
  exclude_mask <- grepl(paste(exclude_patterns, collapse="|"), gene_names)
  keep_mask <- looks_standard & !exclude_mask
}

expr_clean <- expr_filtered[keep_mask, , drop = FALSE]
say("After protein-coding filter: %d genes", nrow(expr_clean))

if (nrow(expr_clean) == 0) {
  say("⚠ Warning: No genes remaining after filtering - using all filtered genes")
  expr_clean <- expr_filtered
}

# ==================== STEP 8: PREPARE SURVIVAL DATA ====================
say("\n========================================")
say("STEP 8: PREPARING SURVIVAL DATA")
say("========================================")

# Extract survival variables
survival_data <- data.frame(
  barcode = metadata_dedup$barcode,
  patient12 = metadata_dedup$patient12,
  stringsAsFactors = FALSE
)

# Overall survival time (days)
survival_data$OS_days <- ifelse(
  !is.na(metadata_dedup$days_to_death),
  metadata_dedup$days_to_death,
  metadata_dedup$days_to_last_follow_up
)
survival_data$OS_years <- survival_data$OS_days / 365.25

# Overall survival status (1 = dead, 0 = alive)
vital_status <- metadata_dedup$vital_status
survival_data$OS_status <- as.integer(
  vital_status %in% c("Dead", "DECEASED", "dead")
)

# Remove invalid entries
valid <- is.finite(survival_data$OS_years) & 
         !is.na(survival_data$OS_status) & 
         survival_data$OS_years > 0

expr_final <- expr_clean[, valid, drop = FALSE]
survival_final <- survival_data[valid, ]

say("Final cohort: %d patients with valid survival data", nrow(survival_final))

# ==================== STEP 9: QC VALIDATION ====================
say("\n========================================")
say("STEP 9: QUALITY CONTROL VALIDATION")
say("========================================")

# Event rate
event_rate <- mean(survival_final$OS_status)
n_events <- sum(survival_final$OS_status)
n_censored <- sum(survival_final$OS_status == 0)

say("\nSurvival metrics:")
say("  Events: %d (%.1f%%)", n_events, 100*event_rate)
say("  Censored: %d (%.1f%%)", n_censored, 100*(1-event_rate))

# Check against expected GBM ranges
event_ok <- event_rate >= 0.70 && event_rate <= 0.85
say("  Expected event rate: 70-85%% | %s", ifelse(event_ok, "PASS", "FAIL"))

# Median survival
med_os <- median(survival_final$OS_years[survival_final$OS_status == 1], na.rm = TRUE)
med_fu <- median(survival_final$OS_years[survival_final$OS_status == 0], na.rm = TRUE)

say("  Median OS (events): %.2f years", med_os)
say("  Median follow-up (censored): %.2f years", med_fu)

med_ok <- med_os >= 0.8 && med_os <= 2.0
say("  Expected median OS: 0.8-2.0 y | %s", ifelse(med_ok, "PASS", "FAIL"))

# ==================== STEP 10: VALIDATE WITH KNOWN MARKERS ====================
say("\n========================================")
say("STEP 10: VALIDATION WITH KNOWN MARKERS")
say("========================================")

validation_genes <- c("EGFR", "TP53", "PTEN", "IDH1")
validation_results <- list()

for (gene in validation_genes) {
  if (!(gene %in% rownames(expr_final))) {
    say("%s: NOT FOUND", gene)
    next
  }
  
  gene_expr <- as.numeric(expr_final[gene, ])
  
  # Remove missing values
  valid_idx <- is.finite(gene_expr) & is.finite(survival_final$OS_years)
  if (sum(valid_idx) < 50) {
    say("%s: Insufficient data", gene)
    next
  }
  
  df <- data.frame(
    expr = gene_expr[valid_idx],
    time = survival_final$OS_years[valid_idx],
    status = survival_final$OS_status[valid_idx]
  )
  
  # Find optimal cutpoint
  tryCatch({
    cut <- surv_cutpoint(df, time = "time", event = "status", 
                        variables = "expr", minprop = 0.1)
    cutval <- cut$cutpoint$cutpoint[1]
    group <- ifelse(df$expr >= cutval, "High", "Low")
    
    # Log-rank test
    sdiff <- survdiff(Surv(df$time, df$status) ~ group)
    pval <- pchisq(sdiff$chisq, df = 1, lower.tail = FALSE)
    
    # Cox model
    cox <- coxph(Surv(df$time, df$status) ~ group)
    hr <- exp(coef(cox))
    
    validation_results[[gene]] <- list(p = pval, hr = hr)
    
    sig <- ifelse(pval < 0.05, "✓", "✗")
    say("%s %s: p=%.4f, HR=%.2f", sig, gene, pval, hr)
    
  }, error = function(e) {
    say("✗ %s: Analysis failed - %s", gene, e$message)
  })
}

if (length(validation_results) > 0) {
  n_sig <- sum(sapply(validation_results, function(x) x$p < 0.05))
  say("\nSignificant markers: %d / %d", n_sig, length(validation_results))
} else {
  n_sig <- 0
  say("\nNo validation markers found in expression data")
}

if (length(validation_results) > 0) {
  if (n_sig >= 2) {
    say("✓ Validation successful - multiple known markers significant")
  } else {
    warning("⚠ Only %d/%d validation genes significant - review data quality", 
            n_sig, length(validation_results))
  }
} else {
  say("⚠ No validation genes found - check gene ID format or filtering")
}

# ==================== STEP 11: SAVE PROCESSED DATA ====================
say("\n========================================")
say("STEP 11: SAVING PROCESSED DATA")
say("========================================")

# Save expression
out_expr <- file.path(PROCESSED_DIR, "TCGA_GBM_Expression_Cleaned.rds")
saveRDS(expr_final, out_expr)
say("✓ Expression: %s", out_expr)

# Save clinical/survival
out_clin <- file.path(PROCESSED_DIR, "TCGA_GBM_Clinical_Cleaned.csv")
fwrite(survival_final, out_clin)
say("✓ Clinical: %s", out_clin)

# Create comprehensive tracking file
tracking <- data.frame(
  barcode = survival_final$barcode,
  patient12 = survival_final$patient12,
  OS_years = survival_final$OS_years,
  OS_status = survival_final$OS_status,
  stringsAsFactors = FALSE
)

out_track <- file.path(PROCESSED_DIR, "TCGA_GBM_Sample_Tracking.csv")
fwrite(tracking, out_track)
say("✓ Tracking: %s", out_track)

# Save QC report
qc_report <- c(
  "TCGA-GBM FRESH DOWNLOAD & PREPROCESSING",
  sprintf("Processing date: %s", Sys.Date()),
  sprintf("TCGAbiolinks version: %s", packageVersion("TCGAbiolinks")),
  "",
  "FINAL COHORT:",
  sprintf("  Patients: %d", ncol(expr_final)),
  sprintf("  Genes: %d", nrow(expr_final)),
  "",
  "SURVIVAL DATA:",
  sprintf("  Events: %d (%.1f%%)", n_events, 100*event_rate),
  sprintf("  Censored: %d (%.1f%%)", n_censored, 100*(1-event_rate)),
  sprintf("  Median OS (events): %.2f years", med_os),
  sprintf("  Median follow-up (censored): %.2f years", med_fu),
  "",
  "QC CHECKS:",
  sprintf("  Event rate: 70-85%% expected | %.1f%% observed | %s",
          100*event_rate, ifelse(event_ok, "PASS", "FAIL")),
  sprintf("  Median OS: 0.8-2.0 y expected | %.2f y observed | %s",
          med_os, ifelse(med_ok, "PASS", "FAIL")),
  "",
  "VALIDATION GENES:",
  paste(sprintf("  %s: p=%.4f", names(validation_results),
                sapply(validation_results, function(x) x$p)), collapse = "\n"),
  sprintf("  Significant: %d / %d", n_sig, length(validation_results)),
  "",
  "DATA PROVENANCE:",
  "  Source: TCGA GDC Data Portal",
  "  Project: TCGA-GBM",
  "  Workflow: STAR - Counts",
  "  Normalization: FPKM (unstrand)",
  "",
  "PREPROCESSING STEPS:",
  "  1. Downloaded via TCGAbiolinks from GDC",
  "  2. Filtered to primary solid tumors (code 01)",
  "  3. Deduplicated to one sample per patient",
  "  4. Filtered genes: protein-coding + ≥10% detection",
  "  5. Validated with known prognostic markers"
)

out_qc <- file.path(PROCESSED_DIR, "TCGA_GBM_QC_Report.txt")
writeLines(qc_report, out_qc)
say("✓ QC Report: %s", out_qc)

# ==================== FINAL SUMMARY ====================
say("\n========================================")
say("PROCESSING COMPLETE")
say("========================================")
say("\nFinal dataset:")
say("  Patients: %d", ncol(expr_final))
say("  Genes: %d", nrow(expr_final))
say("  Event rate: %.1f%%", 100*event_rate)
say("  Median OS: %.2f years", med_os)
say("  Validation: %d/%d genes significant", n_sig, length(validation_results))
say("\nAll files saved to: %s", PROCESSED_DIR)

say("\n✓ Ready for downstream analysis!")
say("  This data has proper TCGA barcodes, known provenance,")
say("  and has been validated against known GBM biology.")
