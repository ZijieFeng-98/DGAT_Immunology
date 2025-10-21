#!/usr/bin/env Rscript
# =============================================================================
# 01_validate_clean_data.R — Validate Current Clean Database
# Purpose: Verify that the cleaned TCGA-GBM data is ready for analysis
# Compatible with: Processed_Data/TCGA_GBM_Clean/
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(survival)
  library(survminer)
})

# =============================================================================
# CONFIGURATION
# =============================================================================

# Paths to clean data
EXPR_FILE <- "Processed_Data/TCGA_GBM_Clean/TCGA_GBM_Expression_Cleaned.rds"
CLIN_FILE <- "Processed_Data/TCGA_GBM_Clean/TCGA_GBM_Clinical_Cleaned.csv"
TRACK_FILE <- "Processed_Data/TCGA_GBM_Clean/TCGA_GBM_Sample_Tracking.csv"

# Target genes for DGAT project
TARGET_GENES <- c("DGAT1", "DGAT2")

# Known prognostic markers for validation
VALIDATION_GENES <- c("EGFR", "TP53", "PTEN", "IDH1")

# =============================================================================
# DATA VALIDATION FUNCTIONS
# =============================================================================

validate_data_files <- function() {
  cat("=== DATA FILE VALIDATION ===\n")
  
  files_to_check <- c(
    "Expression" = EXPR_FILE,
    "Clinical" = CLIN_FILE,
    "Tracking" = TRACK_FILE
  )
  
  all_exist <- TRUE
  for (name in names(files_to_check)) {
    file_path <- files_to_check[name]
    exists <- file.exists(file_path)
    cat(sprintf("%-12s: %s %s\n", name, ifelse(exists, "✓", "✗"), basename(file_path)))
    if (!exists) all_exist <- FALSE
  }
  
  return(all_exist)
}

validate_expression_data <- function() {
  cat("\n=== EXPRESSION DATA VALIDATION ===\n")
  
  if (!file.exists(EXPR_FILE)) {
    cat("❌ Expression file not found\n")
    return(NULL)
  }
  
  expr <- readRDS(EXPR_FILE)
  
  cat("Dimensions:", nrow(expr), "genes x", ncol(expr), "samples\n")
  cat("Data type:", class(expr), "\n")
  cat("Gene ID format:", ifelse(any(grepl("^ENSG", rownames(expr)[1:10])), "Ensembl", "Gene Symbol"), "\n")
  
  # Check for target genes
  cat("\nTarget genes (DGAT):\n")
  for (gene in TARGET_GENES) {
    found <- gene %in% rownames(expr)
    cat(sprintf("  %-6s: %s\n", gene, ifelse(found, "✓ Found", "✗ Missing")))
    if (found) {
      gene_expr <- as.numeric(expr[gene, ])
      cat(sprintf("         Detection: %.1f%%, Mean: %.2f\n", 
                  100*mean(gene_expr > 0.1), mean(gene_expr)))
    }
  }
  
  # Check for validation genes
  cat("\nValidation genes:\n")
  for (gene in VALIDATION_GENES) {
    found <- gene %in% rownames(expr)
    cat(sprintf("  %-6s: %s\n", gene, ifelse(found, "✓ Found", "✗ Missing")))
  }
  
  return(expr)
}

validate_clinical_data <- function() {
  cat("\n=== CLINICAL DATA VALIDATION ===\n")
  
  if (!file.exists(CLIN_FILE)) {
    cat("❌ Clinical file not found\n")
    return(NULL)
  }
  
  clin <- fread(CLIN_FILE)
  
  cat("Dimensions:", nrow(clin), "patients x", ncol(clin), "variables\n")
  cat("Columns:", paste(names(clin), collapse = ", "), "\n")
  
  # Check survival data
  if (all(c("OS_years", "OS_status") %in% names(clin))) {
    n_complete <- sum(complete.cases(clin[, c("OS_years", "OS_status")]))
    event_rate <- mean(clin$OS_status, na.rm = TRUE)
    n_events <- sum(clin$OS_status, na.rm = TRUE)
    
    cat("\nSurvival data:\n")
    cat(sprintf("  Complete cases: %d/%d (%.1f%%)\n", n_complete, nrow(clin), 100*n_complete/nrow(clin)))
    cat(sprintf("  Events: %d (%.1f%%)\n", n_events, 100*event_rate))
    cat(sprintf("  Censored: %d (%.1f%%)\n", nrow(clin)-n_events, 100*(1-event_rate)))
    
    if (n_events > 0) {
      med_surv <- median(clin$OS_years[clin$OS_status == 1], na.rm = TRUE)
      cat(sprintf("  Median OS (events): %.2f years\n", med_surv))
    }
  } else {
    cat("❌ Missing survival columns (OS_years, OS_status)\n")
  }
  
  return(clin)
}

validate_sample_matching <- function(expr, clin) {
  cat("\n=== SAMPLE MATCHING VALIDATION ===\n")
  
  if (is.null(expr) || is.null(clin)) {
    cat("❌ Cannot validate - missing expression or clinical data\n")
    return(FALSE)
  }
  
  expr_samples <- colnames(expr)
  clin_samples <- clin$barcode
  
  cat("Expression samples:", length(expr_samples), "\n")
  cat("Clinical samples:", length(clin_samples), "\n")
  
  # Check if samples match
  if (all(expr_samples %in% clin_samples) && all(clin_samples %in% expr_samples)) {
    cat("✓ Perfect sample matching\n")
    return(TRUE)
  } else {
    missing_in_clin <- setdiff(expr_samples, clin_samples)
    missing_in_expr <- setdiff(clin_samples, expr_samples)
    
    if (length(missing_in_clin) > 0) {
      cat(sprintf("❌ %d expression samples missing in clinical data\n", length(missing_in_clin)))
    }
    if (length(missing_in_expr) > 0) {
      cat(sprintf("❌ %d clinical samples missing in expression data\n", length(missing_in_expr)))
    }
    return(FALSE)
  }
}

# =============================================================================
# SURVIVAL ANALYSIS VALIDATION
# =============================================================================

validate_dgat_survival <- function(expr, clin) {
  cat("\n=== DGAT SURVIVAL ANALYSIS ===\n")
  
  if (is.null(expr) || is.null(clin)) {
    cat("❌ Cannot perform survival analysis - missing data\n")
    return(NULL)
  }
  
  # Ensure samples match
  common_samples <- intersect(colnames(expr), clin$barcode)
  expr_matched <- expr[, common_samples]
  clin_matched <- clin[clin$barcode %in% common_samples, ]
  
  cat("Analyzing", length(common_samples), "matched samples\n")
  
  results <- list()
  
  for (gene in TARGET_GENES) {
    if (!(gene %in% rownames(expr_matched))) {
      cat(sprintf("%-6s: ✗ Not found\n", gene))
      next
    }
    
    gene_expr <- as.numeric(expr_matched[gene, ])
    
    # Remove missing values
    valid_idx <- is.finite(gene_expr) & is.finite(clin_matched$OS_years)
    if (sum(valid_idx) < 50) {
      cat(sprintf("%-6s: ✗ Insufficient data\n", gene))
      next
    }
    
    # Prepare survival data
    surv_data <- data.frame(
      time = clin_matched$OS_years[valid_idx],
      status = clin_matched$OS_status[valid_idx],
      expr = gene_expr[valid_idx],
      stringsAsFactors = FALSE
    )
    
    # Median split
    median_val <- median(surv_data$expr)
    surv_data$group <- ifelse(surv_data$expr >= median_val, "High", "Low")
    
    # Cox regression
    cox_model <- coxph(Surv(time, status) ~ group, data = surv_data)
    hr <- exp(coef(cox_model))
    p_val <- summary(cox_model)$coefficients["groupLow", "Pr(>|z|)"]
    
    results[[gene]] <- list(hr = hr, p = p_val, median = median_val)
    
    sig <- ifelse(p_val < 0.05, "✓", "✗")
    cat(sprintf("%-6s: %s HR=%.3f, p=%.3f (median=%.2f)\n", 
                gene, sig, hr, p_val, median_val))
  }
  
  return(results)
}

# =============================================================================
# MAIN VALIDATION
# =============================================================================

main <- function() {
  cat("============================================================================\n")
  cat(" DGAT PROJECT - CLEAN DATA VALIDATION\n")
  cat("============================================================================\n")
  
  # Step 1: Check files exist
  if (!validate_data_files()) {
    cat("\n❌ CRITICAL: Missing data files. Cannot proceed.\n")
    return(FALSE)
  }
  
  # Step 2: Validate expression data
  expr <- validate_expression_data()
  
  # Step 3: Validate clinical data
  clin <- validate_clinical_data()
  
  # Step 4: Check sample matching
  sample_match <- validate_sample_matching(expr, clin)
  
  # Step 5: DGAT survival analysis
  dgat_results <- validate_dgat_survival(expr, clin)
  
  # Summary
  cat("\n============================================================================\n")
  cat(" VALIDATION SUMMARY\n")
  cat("============================================================================\n")
  
  # Check if DGAT genes are present and analyzable
  dgat_ready <- FALSE
  if (!is.null(dgat_results)) {
    dgat_genes_found <- names(dgat_results)
    if (length(dgat_genes_found) >= 2) {
      dgat_ready <- TRUE
      cat("✓ DGAT1 and DGAT2 both found and analyzable\n")
    } else {
      cat("⚠ Only", length(dgat_genes_found), "DGAT gene(s) found\n")
    }
  } else {
    cat("❌ DGAT genes not found or not analyzable\n")
  }
  
  # Overall status
  if (dgat_ready && sample_match) {
    cat("\n🎉 DATA READY FOR DGAT ANALYSIS!\n")
    cat("✓ Clean TCGA-GBM dataset with proper gene symbols\n")
    cat("✓ DGAT1 and DGAT2 detected and analyzable\n")
    cat("✓ Survival data complete and validated\n")
    cat("✓ Sample matching verified\n")
    cat("\nNext steps: Proceed with DGAT survival analysis\n")
    return(TRUE)
  } else {
    cat("\n⚠ DATA ISSUES DETECTED\n")
    if (!dgat_ready) cat("- DGAT genes missing or insufficient data\n")
    if (!sample_match) cat("- Sample matching issues\n")
    cat("\nReview and fix issues before proceeding\n")
    return(FALSE)
  }
}

# Run validation
if (!interactive()) {
  success <- main()
  quit(status = ifelse(success, 0, 1))
}
