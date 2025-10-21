#!/usr/bin/env Rscript
# =============================================================================
# cgga_fpkm_simple.R
# Simplified pipeline for CGGA FPKM data processing
# Mirrors the TCGA pipeline for consistency
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(ggplot2)
  library(matrixStats)
})

# =========================
# CONFIG
# =========================
CONFIG <- list(
  expr_file     = "Raw_Data/CGGA_GBM/CGGA_693_Expression_Cleaned.tsv",
  clinical_file = "Raw_Data/CGGA_GBM/CGGA_693_Clinical_Cleaned.tsv",
  outdir        = "Processed_Data/CGGA_GBM_Master_Pipeline"
)

dir.create(CONFIG$outdir, recursive = TRUE, showWarnings = FALSE)

# =========================
# Load data
# =========================
cat("Loading CGGA expression data...\n")
expr_data <- fread(CONFIG$expr_file)
clinical <- fread(CONFIG$clinical_file)

# Extract gene symbols and convert to matrix
gene_symbols <- expr_data$GeneSymbol
expr_matrix <- as.matrix(expr_data[, -1])
rownames(expr_matrix) <- gene_symbols

cat("Expression data dimensions:", nrow(expr_matrix), "genes x", ncol(expr_matrix), "samples\n")
cat("Clinical data dimensions:", nrow(clinical), "samples x", ncol(clinical), "variables\n")

# =========================
# Sample filtering
# =========================
cat("\n=== SAMPLE FILTERING ===\n")

# CGGA data is already GBM-only, no need for tumor type filtering
# Check for duplicate patients using CGGA_ID
if ("CGGA_ID" %in% colnames(clinical)) {
  patient_ids <- clinical$CGGA_ID
  unique_patients <- unique(patient_ids)
  
  if (length(unique_patients) < length(patient_ids)) {
    cat("Found duplicates, keeping one per patient...\n")
    # Keep first occurrence of each patient
    keep_idx <- !duplicated(patient_ids)
    clinical <- clinical[keep_idx, ]
    expr_matrix <- expr_matrix[, patient_ids[keep_idx], drop = FALSE]
  }
}

# Align expression and clinical data
common_samples <- intersect(colnames(expr_matrix), clinical$CGGA_ID)
expr_matrix <- expr_matrix[, common_samples, drop = FALSE]
clinical <- clinical[match(common_samples, clinical$CGGA_ID), ]

cat("After alignment:", ncol(expr_matrix), "samples\n")

# =========================
# Gene filtering
# =========================
cat("\n=== GENE FILTERING ===\n")

# Remove ribosomal and mitochondrial genes
ribo_mito <- grepl("^(RPL|RPS|MRPL|MRPS|MT-)", rownames(expr_matrix))
expr_clean <- expr_matrix[!ribo_mito, , drop = FALSE]
cat("Removed ribosomal/mitochondrial genes:", sum(ribo_mito), "\n")

# Detection filter (genes expressed in at least 10% of samples)
# CGGA data appears to be raw counts/FPKM, threshold at 0.1
detection_thresh <- 0.1
detected_prop <- rowMeans(expr_clean > detection_thresh, na.rm = TRUE)
detected_genes <- detected_prop >= 0.1
expr_detected <- expr_clean[detected_genes, , drop = FALSE]

cat("After detection filter:", nrow(expr_detected), "genes\n")

# Log2 transform (adding 1 to avoid log(0))
expr_log <- log2(expr_detected + 1)

# Variance filter
gene_vars <- rowVars(as.matrix(expr_log))
var_thresh <- quantile(gene_vars, 0.1, na.rm = TRUE)
var_genes <- gene_vars >= var_thresh
expr_final <- expr_log[var_genes, , drop = FALSE]

cat("After variance filter:", nrow(expr_final), "genes\n")

# =========================
# Quality control
# =========================
cat("\n=== QUALITY CONTROL ===\n")

# Create metadata for QC
qc_meta <- data.frame(
  sample_id = colnames(expr_final),
  patient_id = colnames(expr_final),
  stringsAsFactors = FALSE
)

# Add clinical variables if available
if ("Age" %in% colnames(clinical)) {
  qc_meta$age <- clinical$Age[match(colnames(expr_final), clinical$CGGA_ID)]
}

if ("Gender" %in% colnames(clinical)) {
  qc_meta$gender <- clinical$Gender[match(colnames(expr_final), clinical$CGGA_ID)]
}

# Create batch variable from patient ID patterns (if applicable)
# CGGA IDs typically have patterns like CGGA_1002, CGGA_P100, etc.
# Extract prefix as potential batch
qc_meta$batch <- gsub("^(CGGA_[A-Z]?).*", "\\1", qc_meta$patient_id)
# If all same, use a single batch
if (length(unique(qc_meta$batch)) == 1) {
  qc_meta$batch <- "batch1"
}

# Add purity if available (CGGA doesn't typically have this)
qc_meta$purity <- NA_real_

cat("QC metadata created for", nrow(qc_meta), "samples\n")

# =========================
# PCA analysis
# =========================
cat("\n=== PCA ANALYSIS ===\n")

# PCA analysis
pr <- prcomp(t(expr_final), scale. = TRUE)
pcs <- as.data.frame(pr$x[, 1:10])

# Add metadata to PCA results
pcs$sample_id <- rownames(pcs)
pcs$batch <- qc_meta$batch
if ("age" %in% colnames(qc_meta)) {
  pcs$age <- qc_meta$age
}
if ("gender" %in% colnames(qc_meta)) {
  pcs$gender <- qc_meta$gender
}

# Calculate variance explained
var_explained <- pr$sdev^2 / sum(pr$sdev^2)

# Create PCA plots
p1 <- ggplot(pcs, aes(PC1, PC2, color = batch)) +
  geom_point(size = 2) +
  labs(title = "PCA by Batch",
       subtitle = sprintf("PC1: %.1f%% | PC2: %.1f%%", var_explained[1]*100, var_explained[2]*100)) +
  theme_minimal()

# PCA by gender if available
if ("gender" %in% colnames(pcs)) {
  p2 <- ggplot(pcs, aes(PC1, PC2, color = gender)) +
    geom_point(size = 2) +
    labs(title = "PCA by Gender",
         subtitle = sprintf("PC1: %.1f%% | PC2: %.1f%%", var_explained[1]*100, var_explained[2]*100)) +
    theme_minimal()
  ggsave(file.path(CONFIG$outdir, "PCA_by_Gender.png"), p2, width = 8, height = 6, dpi = 300)
}

# PCA by age if available
if ("age" %in% colnames(pcs)) {
  pcs$age_numeric <- suppressWarnings(as.numeric(pcs$age))
  if (sum(!is.na(pcs$age_numeric)) > 0) {
    p3 <- ggplot(pcs[!is.na(pcs$age_numeric), ], aes(PC1, PC2, color = age_numeric)) +
      geom_point(size = 2) +
      scale_color_gradient(low = "blue", high = "red") +
      labs(title = "PCA by Age",
           subtitle = sprintf("PC1: %.1f%% | PC2: %.1f%%", var_explained[1]*100, var_explained[2]*100),
           color = "Age") +
      theme_minimal()
    ggsave(file.path(CONFIG$outdir, "PCA_by_Age.png"), p3, width = 8, height = 6, dpi = 300)
  }
}

# Save batch plot
ggsave(file.path(CONFIG$outdir, "PCA_by_Batch.png"), p1, width = 8, height = 6, dpi = 300)

cat("PCA plots saved\n")

# =========================
# Save results
# =========================
cat("\n=== SAVING RESULTS ===\n")

# Save processed expression data
saveRDS(expr_final, file.path(CONFIG$outdir, "expression_processed.rds"))

# Save metadata
fwrite(qc_meta, file.path(CONFIG$outdir, "metadata_processed.csv"))

# Save PCA results
fwrite(pcs, file.path(CONFIG$outdir, "pca_results.csv"))

# Create summary report
summary_report <- data.frame(
  Metric = c("Original samples", "After filtering", "Final samples", 
             "Original genes", "After detection filter", "After variance filter"),
  Count = c(ncol(expr_matrix), ncol(expr_clean), ncol(expr_final),
            nrow(expr_matrix), nrow(expr_detected), nrow(expr_final)),
  stringsAsFactors = FALSE
)

fwrite(summary_report, file.path(CONFIG$outdir, "processing_summary.csv"))

# Print summary
cat("\n=== PROCESSING SUMMARY ===\n")
print(summary_report)

# Check for DGAT genes
dgat_genes <- c("DGAT1", "DGAT2")
dgat_present <- dgat_genes %in% rownames(expr_final)
cat("\nDGAT gene presence:\n")
for (i in seq_along(dgat_genes)) {
  cat(sprintf("  %s: %s\n", dgat_genes[i], ifelse(dgat_present[i], "✓ Present", "✗ Missing")))
}

# Check expression levels of DGAT genes
if (all(dgat_present)) {
  cat("\nDGAT gene expression levels (log2 scale):\n")
  for (gene in dgat_genes) {
    expr_vals <- expr_final[gene, ]
    cat(sprintf("  %s: Mean=%.3f, Median=%.3f, Range=[%.3f, %.3f]\n", 
                gene, mean(expr_vals), median(expr_vals), min(expr_vals), max(expr_vals)))
  }
}

cat("\n=== PIPELINE COMPLETE ===\n")
cat("Results saved to:", CONFIG$outdir, "\n")
cat("Processed expression data:", nrow(expr_final), "genes x", ncol(expr_final), "samples\n")


