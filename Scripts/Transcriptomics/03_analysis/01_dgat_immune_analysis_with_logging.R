#!/usr/bin/env Rscript
# =============================================================================
# 01_dgat_immune_analysis_with_logging.R — EXAMPLE with Automated Logging
# =============================================================================
#
# This is an EXAMPLE showing how to integrate the logging system
# into your existing analysis scripts.
#
# Simply add log_activity() calls at key points in your workflow!
# =============================================================================

# ============================================================================
# INITIALIZE LOGGING (ADD THIS TO ANY SCRIPT)
# ============================================================================
source("Scripts/log_activity.R")
log_script_start("01_dgat_immune_analysis.R", "Comprehensive DGAT-Immune Analysis")

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(ggplot2)
  library(ggrepel)
  library(pheatmap)
  library(GSVA)
})

log_activity("Loaded required R packages", category = "script")

# =============================================================================
# YOUR ORIGINAL CODE WITH LOGGING ADDED AT KEY POINTS
# =============================================================================

CONFIG <- list(
  tcga_expr = "Processed_Data/TCGA_GBM_Batch_Corrected/expression_batch_corrected.rds",
  cgga_expr = "Processed_Data/CGGA_GBM_Batch_Corrected/expression_batch_corrected.rds",
  output_dir = "Results/Bulk/Immune_Analysis",
  dgat_gene = "DGAT1",
  min_geneset_size = 3,
  max_geneset_size = 500,
  correlation_method = "spearman",
  fdr_threshold = 0.05
)

dir.create(CONFIG$output_dir, recursive = TRUE, showWarnings = FALSE)

log_activity("Initialized analysis configuration", category = "planning")

# ... [Rest of your original functions here - I'll show key logging points]

# Example: In your load_expression_data function
load_expression_data <- function(expr_path, cohort_name) {
  cat("\n📖 Loading", cohort_name, "expression data...\n")
  
  if (!file.exists(expr_path)) {
    log_error(paste("File not found:", expr_path))
    stop("File not found: ", expr_path)
  }
  
  mat <- readRDS(expr_path)
  if (nrow(mat) < ncol(mat)) mat <- t(mat)
  mat <- as.matrix(mat)
  storage.mode(mat) <- "numeric"
  
  # LOG DATA LOADING
  log_data(cohort_name, 
           sprintf("Loaded expression data: %d genes × %d samples", 
                   nrow(mat), ncol(mat)))
  
  if (any(duplicated(rownames(mat)))) {
    n_dup <- sum(duplicated(rownames(mat)))
    DT <- as.data.table(mat, keep.rownames = "gene")
    mat <- as.matrix(DT[, lapply(.SD, mean, na.rm = TRUE), by = gene] %>% 
                     tibble::column_to_rownames("gene"))
    log_activity(sprintf("Aggregated %d duplicated genes in %s", n_dup, cohort_name),
                category = "data")
  }
  
  return(mat)
}

# Example: In your run_gsva_analysis function
run_gsva_analysis <- function(expr_mat, genesets, cohort_name) {
  cat("\n🔬 Running GSVA analysis for", cohort_name, "...\n")
  
  # ... filtering code ...
  
  log_activity(sprintf("Running GSVA on %d immune gene sets for %s", 
                      length(filtered_genesets), cohort_name),
              category = "analysis")
  
  gsva_scores <- tryCatch({
    gsva(expr_mat, filtered_genesets, method = "gsva", kcdf = "Gaussian")
  }, error = function(e) {
    log_error(sprintf("GSVA failed for %s: %s", cohort_name, e$message))
    NULL
  })
  
  if (!is.null(gsva_scores)) {
    log_result(sprintf("GSVA complete for %s: %d gene sets scored", 
                      cohort_name, nrow(gsva_scores)))
  }
  
  return(gsva_scores)
}

# Example: In your perform_differential_analysis function
perform_differential_analysis <- function(gsva_scores, groups) {
  cat("\n📊 Performing differential analysis...\n")
  
  log_activity("Running differential analysis (Wilcoxon test)", category = "analysis")
  
  # ... your differential analysis code ...
  
  n_sig <- sum(results$FDR < CONFIG$fdr_threshold)
  log_result(sprintf("Differential analysis: %d significant gene sets (FDR < %.2f)", 
                    n_sig, CONFIG$fdr_threshold))
  
  if (n_sig > 0) {
    top_result <- results[1, ]
    log_result(sprintf("Top hit: %s (Δ=%.3f, FDR=%.2e)", 
                      top_result$GeneSet, top_result$Delta, top_result$FDR))
  }
  
  return(results)
}

# Example: In your perform_correlation_analysis function
perform_correlation_analysis <- function(gsva_scores, dgat_expr) {
  cat("\n📈 Performing correlation analysis...\n")
  
  log_activity(sprintf("Calculating %s correlations with DGAT1", 
                      CONFIG$correlation_method),
              category = "analysis")
  
  # ... your correlation code ...
  
  n_sig_cor <- sum(cor_results$FDR < CONFIG$fdr_threshold)
  log_result(sprintf("Correlation analysis: %d significant correlations (FDR < %.2f)",
                    n_sig_cor, CONFIG$fdr_threshold))
  
  if (n_sig_cor > 0) {
    strongest <- cor_results[which.max(abs(cor_results$rho)), ]
    log_result(sprintf("Strongest correlation: %s (ρ=%.3f, FDR=%.2e)",
                      strongest$GeneSet, strongest$rho, strongest$FDR))
  }
  
  return(cor_results)
}

# Example: Main workflow with logging
run_immune_analysis <- function(expr_path, cohort_name) {
  cat("\n", strrep("=", 70), "\n")
  cat("🔬 Running immune analysis for", cohort_name, "\n")
  cat(strrep("=", 70), "\n")
  
  log_activity(sprintf("Starting immune analysis workflow for %s", cohort_name),
              category = "analysis")
  
  # ... your workflow code ...
  
  # At the end
  log_activity(sprintf("Completed immune analysis for %s", cohort_name),
              category = "result")
  log_data(cohort_name, 
           sprintf("Saved results to %s", 
                   file.path(CONFIG$output_dir, cohort_name)))
  
  return(list(...))
}

# =============================================================================
# EXECUTION
# =============================================================================

log_activity("Starting TCGA analysis", category = "analysis")
if (file.exists(CONFIG$tcga_expr)) {
  tcga_results <- run_immune_analysis(CONFIG$tcga_expr, "TCGA_GBM")
} else {
  log_error("TCGA data not found")
}

log_activity("Starting CGGA analysis", category = "analysis")
if (file.exists(CONFIG$cgga_expr)) {
  cgga_results <- run_immune_analysis(CONFIG$cgga_expr, "CGGA_GBM")
} else {
  log_error("CGGA data not found")
}

# =============================================================================
# END OF SCRIPT - LOG COMPLETION
# =============================================================================
log_script_end("01_dgat_immune_analysis.R", 
               "Completed DGAT-immune analysis for TCGA and CGGA")

cat("\n✅ Check your daily log: Notes/Daily_Logs/WorkLog_", Sys.Date(), ".md\n\n", sep="")

