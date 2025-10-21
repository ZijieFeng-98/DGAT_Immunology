#!/usr/bin/env Rscript
# =============================================================================
# Script 01: Data Preparation & Cleaning
# Purpose: Load, clean, and prepare TCGA-GBM and CGGA-GBM data for analysis
# Features: GDC sample filtering, batch effect detection, purity estimation
# =============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(data.table)
  library(SummarizedExperiment)
  library(ggplot2)
  library(gridExtra)
  library(estimate)
  library(org.Hs.eg.db)
  library(AnnotationDbi)
})

# =============================================================================
# CONFIGURATION
# =============================================================================

# Data paths
TCGA_EXPR_PATH <- "Processed_Data/TCGA_GBM_Clean/TCGA_GBM_Expression_Cleaned.rds"
TCGA_CLIN_PATH <- "Processed_Data/TCGA_GBM_Clean/TCGA_GBM_Clinical_Cleaned.csv"
CGGA_EXPR_PATH <- "Processed_Data/Clean_CGGA_Cohort/mRNAseq_693_GBM/CGGA_mRNAseq_693_GBM_expr_clean.rds"
CGGA_CLIN_PATH <- "Processed_Data/Clean_CGGA_Cohort/mRNAseq_693_GBM/CGGA_mRNAseq_693_GBM_clin_clean.csv"

# Output directories
OUTPUT_DIR <- "Results/Analysis_Framework"
TCGA_OUT_DIR <- file.path(OUTPUT_DIR, "TCGA_GBM/01_cleaned_data")
CGGA_OUT_DIR <- file.path(OUTPUT_DIR, "CGGA_GBM/01_cleaned_data")
COMBINED_OUT_DIR <- file.path(OUTPUT_DIR, "Combined/01_cleaned_data")

# Create output directories
dir.create(TCGA_OUT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(CGGA_OUT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(COMBINED_OUT_DIR, recursive = TRUE, showWarnings = FALSE)

# Target genes
TARGET_GENES <- c("DGAT1", "DGAT2")

# Quality control thresholds
QC_THRESHOLDS <- list(
  min_library_size = 1e6,
  min_genes_detected = 5000,
  max_mt_fraction = 0.20,
  min_detection_rate = 0.10
)

say <- function(...) cat(sprintf(...), "\n")

# =============================================================================
# HELPER FUNCTIONS
# =============================================================================

# Batch effect detection using PCA
pca_batch_check <- function(expr, metadata, title = "PCA Batch Check") {
  # Remove genes with zero variance
  gene_vars <- apply(expr, 1, var, na.rm = TRUE)
  expr_filtered <- expr[gene_vars > 0, ]
  
  # Perform PCA
  pca <- prcomp(t(expr_filtered), scale = TRUE)
  pca_df <- data.frame(
    PC1 = pca$x[, 1],
    PC2 = pca$x[, 2],
    PC3 = pca$x[, 3],
    Sample = colnames(expr)
  )
  
  # Add metadata if available
  if (!is.null(metadata)) {
    # Extract plate/TSS information from sample IDs if available
    if (any(grepl("TCGA-", colnames(expr)))) {
      # TCGA barcode parsing
      pca_df$Plate <- substr(colnames(expr), 22, 25)
      pca_df$TSS <- substr(colnames(expr), 6, 7)
    }
    
    # Add other metadata
    for (col in c("age", "sex", "grade")) {
      if (col %in% names(metadata)) {
        pca_df[[col]] <- metadata[[col]][match(pca_df$Sample, metadata$sample_id)]
      }
    }
  }
  
  # Create plots
  plots <- list()
  
  if ("Plate" %in% names(pca_df)) {
    plots$plate <- ggplot(pca_df, aes(PC1, PC2, color = Plate)) +
      geom_point(size = 2, alpha = 0.7) +
      labs(title = paste(title, "- Plate Effect"),
           x = paste0("PC1 (", round(summary(pca)$importance[2, 1] * 100, 1), "%)"),
           y = paste0("PC2 (", round(summary(pca)$importance[2, 2] * 100, 1), "%)")) +
      theme_minimal()
  }
  
  if ("TSS" %in% names(pca_df)) {
    plots$tss <- ggplot(pca_df, aes(PC1, PC2, color = TSS)) +
      geom_point(size = 2, alpha = 0.7) +
      labs(title = paste(title, "- TSS Effect"),
           x = paste0("PC1 (", round(summary(pca)$importance[2, 1] * 100, 1), "%)"),
           y = paste0("PC2 (", round(summary(pca)$importance[2, 2] * 100, 1), "%)")) +
      theme_minimal()
  }
  
  # PC1 vs PC3 for additional insight
  plots$pc13 <- ggplot(pca_df, aes(PC1, PC3)) +
    geom_point(size = 2, alpha = 0.7) +
    labs(title = paste(title, "- PC1 vs PC3"),
         x = paste0("PC1 (", round(summary(pca)$importance[2, 1] * 100, 1), "%)"),
         y = paste0("PC3 (", round(summary(pca)$importance[2, 3] * 100, 1), "%)")) +
    theme_minimal()
  
  return(list(
    plots = plots,
    pca_data = pca_df,
    variance_explained = summary(pca)$importance,
    pca_object = pca
  ))
}

# Calculate ESTIMATE purity scores
calculate_estimate_purity <- function(expr_matrix) {
  say("Calculating ESTIMATE purity scores...")
  
  # ESTIMATE requires gene symbols and specific format
  estimate_result <- tryCatch({
    # Convert to ESTIMATE format
    expr_df <- data.frame(
      GeneSymbol = rownames(expr_matrix),
      expr_matrix,
      stringsAsFactors = FALSE
    )
    
    # Calculate ESTIMATE scores
    estimate_scores <- ESTIMATE::estimateScore(
      expr_df,
      platform = "illumina"  # TCGA uses Illumina HiSeq
    )
    
    return(estimate_scores)
  }, error = function(e) {
    say("ESTIMATE calculation failed: %s", e$message)
    return(NULL)
  })
  
  return(estimate_result)
}

# Quality control metrics
calculate_qc_metrics <- function(expr_matrix, dataset_name) {
  say("Calculating QC metrics for %s...", dataset_name)
  
  # Library sizes
  library_sizes <- colSums(expr_matrix, na.rm = TRUE)
  
  # Genes detected per sample
  genes_detected <- colSums(expr_matrix > 0, na.rm = TRUE)
  
  # Mitochondrial gene fraction
  mt_genes <- grepl("^MT-", rownames(expr_matrix))
  mt_fraction <- if (any(mt_genes)) {
    colSums(expr_matrix[mt_genes, ], na.rm = TRUE) / library_sizes
  } else {
    rep(0, ncol(expr_matrix))
  }
  
  # Detection rate per gene
  detection_rate <- rowMeans(expr_matrix > 0, na.rm = TRUE)
  
  # QC summary
  qc_summary <- data.frame(
    Sample = colnames(expr_matrix),
    Library_Size = library_sizes,
    Genes_Detected = genes_detected,
    MT_Fraction = mt_fraction,
    Dataset = dataset_name
  )
  
  # Flag low-quality samples
  qc_summary$Low_Library <- library_sizes < QC_THRESHOLDS$min_library_size
  qc_summary$Low_Genes <- genes_detected < QC_THRESHOLDS$min_genes_detected
  qc_summary$High_MT <- mt_fraction > QC_THRESHOLDS$max_mt_fraction
  qc_summary$QC_Pass <- !(qc_summary$Low_Library | qc_summary$Low_Genes | qc_summary$High_MT)
  
  # Gene-level QC
  gene_qc <- data.frame(
    Gene = rownames(expr_matrix),
    Detection_Rate = detection_rate,
    Mean_Expression = rowMeans(expr_matrix, na.rm = TRUE),
    Median_Expression = apply(expr_matrix, 1, median, na.rm = TRUE)
  )
  
  return(list(
    sample_qc = qc_summary,
    gene_qc = gene_qc,
    summary_stats = list(
      n_samples = ncol(expr_matrix),
      n_genes = nrow(expr_matrix),
      n_passed_qc = sum(qc_summary$QC_Pass),
      median_library_size = median(library_sizes),
      median_genes_detected = median(genes_detected),
      median_mt_fraction = median(mt_fraction)
    )
  ))
}

# =============================================================================
# TCGA-GBM DATA PREPARATION
# =============================================================================

say("\n=== TCGA-GBM Data Preparation ===")

# Load TCGA data
tcga_expr <- readRDS(TCGA_EXPR_PATH)
tcga_clin <- fread(TCGA_CLIN_PATH)

say("TCGA Expression: %d genes × %d samples", nrow(tcga_expr), ncol(tcga_expr))
say("TCGA Clinical: %d patients", nrow(tcga_clin))

# Log2 transform expression (add pseudocount)
tcga_expr_log2 <- log2(tcga_expr + 1)

# Calculate QC metrics
tcga_qc <- calculate_qc_metrics(tcga_expr_log2, "TCGA-GBM")

# Check target gene presence
tcga_target_genes <- TARGET_GENES[TARGET_GENES %in% rownames(tcga_expr_log2)]
say("TCGA Target genes present: %s", paste(tcga_target_genes, collapse = ", "))

# Calculate ESTIMATE purity
tcga_estimate <- calculate_estimate_purity(tcga_expr_log2)

# Batch effect analysis
tcga_batch <- pca_batch_check(tcga_expr_log2, tcga_clin, "TCGA-GBM")

# Save TCGA results
saveRDS(tcga_expr_log2, file.path(TCGA_OUT_DIR, "tcga_gbm_expression_log2.rds"))
fwrite(tcga_clin, file.path(TCGA_OUT_DIR, "tcga_gbm_clinical.csv"))
fwrite(tcga_qc$sample_qc, file.path(TCGA_OUT_DIR, "tcga_gbm_sample_qc.csv"))
fwrite(tcga_qc$gene_qc, file.path(TCGA_OUT_DIR, "tcga_gbm_gene_qc.csv"))

if (!is.null(tcga_estimate)) {
  fwrite(tcga_estimate, file.path(TCGA_OUT_DIR, "tcga_gbm_estimate_purity.csv"))
}

# Save PCA plots
if (length(tcga_batch$plots) > 0) {
  png(file.path(TCGA_OUT_DIR, "tcga_gbm_pca_batch_effects.png"), 
      width = 12, height = 4, units = "in", res = 300)
  do.call(grid.arrange, c(tcga_batch$plots, ncol = length(tcga_batch$plots)))
  dev.off()
}

# Save QC summary
tcga_summary <- c(
  "TCGA-GBM Data Preparation Summary",
  sprintf("Date: %s", Sys.time()),
  "",
  "Dataset Statistics:",
  sprintf("  Samples: %d", tcga_qc$summary_stats$n_samples),
  sprintf("  Genes: %d", tcga_qc$summary_stats$n_genes),
  sprintf("  QC Passed: %d (%.1f%%)", tcga_qc$summary_stats$n_passed_qc, 
          100 * tcga_qc$summary_stats$n_passed_qc / tcga_qc$summary_stats$n_samples),
  "",
  "Quality Metrics:",
  sprintf("  Median Library Size: %.0f", tcga_qc$summary_stats$median_library_size),
  sprintf("  Median Genes Detected: %.0f", tcga_qc$summary_stats$median_genes_detected),
  sprintf("  Median MT Fraction: %.3f", tcga_qc$summary_stats$median_mt_fraction),
  "",
  "Target Genes:",
  sprintf("  Present: %s", paste(tcga_target_genes, collapse = ", ")),
  sprintf("  Missing: %s", paste(setdiff(TARGET_GENES, tcga_target_genes), collapse = ", ")),
  "",
  "PCA Variance Explained:",
  sprintf("  PC1: %.1f%%", tcga_batch$variance_explained[2, 1] * 100),
  sprintf("  PC2: %.1f%%", tcga_batch$variance_explained[2, 2] * 100),
  sprintf("  PC3: %.1f%%", tcga_batch$variance_explained[2, 3] * 100)
)

writeLines(tcga_summary, file.path(TCGA_OUT_DIR, "tcga_gbm_preparation_summary.txt"))

# =============================================================================
# CGGA-GBM DATA PREPARATION
# =============================================================================

say("\n=== CGGA-GBM Data Preparation ===")

# Load CGGA data
cgga_expr <- readRDS(CGGA_EXPR_PATH)
cgga_clin <- fread(CGGA_CLIN_PATH)

say("CGGA Expression: %d genes × %d samples", nrow(cgga_expr), ncol(cgga_expr))
say("CGGA Clinical: %d patients", nrow(cgga_clin))

# Log2 transform expression (add pseudocount)
cgga_expr_log2 <- log2(cgga_expr + 1)

# Calculate QC metrics
cgga_qc <- calculate_qc_metrics(cgga_expr_log2, "CGGA-GBM")

# Check target gene presence
cgga_target_genes <- TARGET_GENES[TARGET_GENES %in% rownames(cgga_expr_log2)]
say("CGGA Target genes present: %s", paste(cgga_target_genes, collapse = ", "))

# Calculate ESTIMATE purity
cgga_estimate <- calculate_estimate_purity(cgga_expr_log2)

# Batch effect analysis
cgga_batch <- pca_batch_check(cgga_expr_log2, cgga_clin, "CGGA-GBM")

# Save CGGA results
saveRDS(cgga_expr_log2, file.path(CGGA_OUT_DIR, "cgga_gbm_expression_log2.rds"))
fwrite(cgga_clin, file.path(CGGA_OUT_DIR, "cgga_gbm_clinical.csv"))
fwrite(cgga_qc$sample_qc, file.path(CGGA_OUT_DIR, "cgga_gbm_sample_qc.csv"))
fwrite(cgga_qc$gene_qc, file.path(CGGA_OUT_DIR, "cgga_gbm_gene_qc.csv"))

if (!is.null(cgga_estimate)) {
  fwrite(cgga_estimate, file.path(CGGA_OUT_DIR, "cgga_gbm_estimate_purity.csv"))
}

# Save PCA plots
if (length(cgga_batch$plots) > 0) {
  png(file.path(CGGA_OUT_DIR, "cgga_gbm_pca_batch_effects.png"), 
      width = 12, height = 4, units = "in", res = 300)
  do.call(grid.arrange, c(cgga_batch$plots, ncol = length(cgga_batch$plots)))
  dev.off()
}

# Save QC summary
cgga_summary <- c(
  "CGGA-GBM Data Preparation Summary",
  sprintf("Date: %s", Sys.time()),
  "",
  "Dataset Statistics:",
  sprintf("  Samples: %d", cgga_qc$summary_stats$n_samples),
  sprintf("  Genes: %d", cgga_qc$summary_stats$n_genes),
  sprintf("  QC Passed: %d (%.1f%%)", cgga_qc$summary_stats$n_passed_qc, 
          100 * cgga_qc$summary_stats$n_passed_qc / cgga_qc$summary_stats$n_samples),
  "",
  "Quality Metrics:",
  sprintf("  Median Library Size: %.0f", cgga_qc$summary_stats$median_library_size),
  sprintf("  Median Genes Detected: %.0f", cgga_qc$summary_stats$median_genes_detected),
  sprintf("  Median MT Fraction: %.3f", cgga_qc$summary_stats$median_mt_fraction),
  "",
  "Target Genes:",
  sprintf("  Present: %s", paste(cgga_target_genes, collapse = ", ")),
  sprintf("  Missing: %s", paste(setdiff(TARGET_GENES, cgga_target_genes), collapse = ", ")),
  "",
  "PCA Variance Explained:",
  sprintf("  PC1: %.1f%%", cgga_batch$variance_explained[2, 1] * 100),
  sprintf("  PC2: %.1f%%", cgga_batch$variance_explained[2, 2] * 100),
  sprintf("  PC3: %.1f%%", cgga_batch$variance_explained[2, 3] * 100)
)

writeLines(cgga_summary, file.path(CGGA_OUT_DIR, "cgga_gbm_preparation_summary.txt"))

# =============================================================================
# COMBINED ANALYSIS SETUP
# =============================================================================

say("\n=== Combined Analysis Setup ===")

# Find common genes between datasets
common_genes <- intersect(rownames(tcga_expr_log2), rownames(cgga_expr_log2))
say("Common genes between TCGA and CGGA: %d", length(common_genes))

# Create combined expression matrix (genes as rows, samples as columns)
# Note: This is for cross-dataset analysis, not merging samples
combined_info <- data.frame(
  Gene = common_genes,
  TCGA_Present = TRUE,
  CGGA_Present = TRUE,
  stringsAsFactors = FALSE
)

# Add target gene information
combined_info$Is_Target <- combined_info$Gene %in% TARGET_GENES

# Save combined information
fwrite(combined_info, file.path(COMBINED_OUT_DIR, "common_genes_info.csv"))

# Create overall summary
overall_summary <- c(
  "DGAT Immunology Analysis - Data Preparation Complete",
  sprintf("Date: %s", Sys.time()),
  "",
  "Dataset Comparison:",
  sprintf("  TCGA-GBM: %d samples, %d genes", 
          ncol(tcga_expr_log2), nrow(tcga_expr_log2)),
  sprintf("  CGGA-GBM: %d samples, %d genes", 
          ncol(cgga_expr_log2), nrow(cgga_expr_log2)),
  sprintf("  Common genes: %d", length(common_genes)),
  "",
  "Target Gene Status:",
  sprintf("  TCGA: %s", paste(tcga_target_genes, collapse = ", ")),
  sprintf("  CGGA: %s", paste(cgga_target_genes, collapse = ", ")),
  "",
  "Next Steps:",
  "  1. Run 02_deconvolution.R for immune cell estimation",
  "  2. Run 03_gsva.R for pathway scoring",
  "  3. Run 04_correlations.R for purity-adjusted correlations",
  "  4. Run 05_models.R for survival analysis",
  "  5. Run 06_figures.R for final visualization"
)

writeLines(overall_summary, file.path(COMBINED_OUT_DIR, "data_preparation_summary.txt"))

say("\n=== Data Preparation Complete ===")
say("TCGA-GBM outputs saved to: %s", TCGA_OUT_DIR)
say("CGGA-GBM outputs saved to: %s", CGGA_OUT_DIR)
say("Combined analysis info saved to: %s", COMBINED_OUT_DIR)
say("\nReady for Script 02: Immune Deconvolution")
