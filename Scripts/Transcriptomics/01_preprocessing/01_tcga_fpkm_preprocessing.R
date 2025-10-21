#!/usr/bin/env Rscript
# =============================================================================
# gbm_fpkm_simple.R
# Simplified pipeline for FPKM data processing
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
  fpkm_file     = "Raw_Data/TCGA_GBM/TCGA_GBM_Expression_FPKM_with_symbols.rds",
  metadata_file = "Raw_Data/TCGA_GBM/TCGA_GBM_Metadata_with_symbols.csv",
  outdir        = "Processed_Data/TCGA_GBM_Master_Pipeline"
)

dir.create(CONFIG$outdir, recursive = TRUE, showWarnings = FALSE)

# =========================
# Load data
# =========================
cat("Loading FPKM data...\n")
fpkm_data <- readRDS(CONFIG$fpkm_file)
meta <- fread(CONFIG$metadata_file)

# Convert to matrix if needed
if (inherits(fpkm_data, "data.frame")) {
  fpkm_data <- as.matrix(fpkm_data)
}

cat("FPKM data dimensions:", nrow(fpkm_data), "genes x", ncol(fpkm_data), "samples\n")
cat("Metadata dimensions:", nrow(meta), "samples x", ncol(meta), "variables\n")

# =========================
# Sample filtering
# =========================
cat("\n=== SAMPLE FILTERING ===\n")

# Filter for primary tumors (sample type "01")
primary_samples <- meta$sample_type_id == "01" | grepl("-01[A-Z]-", meta$barcode)
meta_filtered <- meta[primary_samples, ]
fpkm_filtered <- fpkm_data[, meta_filtered$barcode, drop = FALSE]

cat("Primary tumor samples:", ncol(fpkm_filtered), "out of", ncol(fpkm_data), "\n")

# Deduplicate patients (keep one sample per patient)
patient_ids <- substr(colnames(fpkm_filtered), 1, 12)
unique_patients <- unique(patient_ids)
keep_samples <- character(0)

for (pat in unique_patients) {
  pat_samples <- colnames(fpkm_filtered)[patient_ids == pat]
  if (length(pat_samples) == 1) {
    keep_samples <- c(keep_samples, pat_samples)
  } else {
    # Prefer -01A- samples, then largest library size
    pref_order <- order(!grepl("-01A-", pat_samples))
    lib_sizes <- colSums(fpkm_filtered[, pat_samples, drop = FALSE])
    size_order <- order(-lib_sizes)
    final_order <- order(pref_order, size_order)
    keep_samples <- c(keep_samples, pat_samples[final_order[1]])
  }
}

fpkm_dedup <- fpkm_filtered[, keep_samples, drop = FALSE]
meta_dedup <- meta_filtered[meta_filtered$barcode %in% keep_samples, ]

cat("After deduplication:", ncol(fpkm_dedup), "samples\n")

# =========================
# Gene filtering
# =========================
cat("\n=== GENE FILTERING ===\n")

# Remove ribosomal and mitochondrial genes
ribo_mito <- grepl("^(RPL|RPS|MRPL|MRPS|MT-)", rownames(fpkm_dedup))
fpkm_clean <- fpkm_dedup[!ribo_mito, , drop = FALSE]
cat("Removed ribosomal/mitochondrial genes:", sum(ribo_mito), "\n")

# Detection filter (genes expressed in at least 10% of samples)
detection_thresh <- 0.1  # FPKM > 0.1
detected_prop <- rowMeans(fpkm_clean > detection_thresh, na.rm = TRUE)
detected_genes <- detected_prop >= 0.1
fpkm_detected <- fpkm_clean[detected_genes, , drop = FALSE]

cat("After detection filter:", nrow(fpkm_detected), "genes\n")

# Variance filter
fpkm_log <- log2(fpkm_detected + 1)
gene_vars <- rowVars(as.matrix(fpkm_log))
var_thresh <- quantile(gene_vars, 0.1, na.rm = TRUE)
var_genes <- gene_vars >= var_thresh
fpkm_final <- fpkm_log[var_genes, , drop = FALSE]

cat("After variance filter:", nrow(fpkm_final), "genes\n")

# =========================
# Quality control
# =========================
cat("\n=== QUALITY CONTROL ===\n")

# Create metadata for QC
qc_meta <- data.frame(
  sample_id = colnames(fpkm_final),
  patient_id = substr(colnames(fpkm_final), 1, 12),
  tss = substr(colnames(fpkm_final), 6, 7),
  batch = substr(colnames(fpkm_final), 6, 7),  # Use TSS as batch
  age = meta_dedup$age_at_diagnosis[match(colnames(fpkm_final), meta_dedup$barcode)],
  gender = meta_dedup$gender[match(colnames(fpkm_final), meta_dedup$barcode)],
  stringsAsFactors = FALSE
)

# Add purity estimates if available
if ("paper_ABSOLUTE.purity" %in% colnames(meta_dedup)) {
  qc_meta$purity <- meta_dedup$paper_ABSOLUTE.purity[match(colnames(fpkm_final), meta_dedup$barcode)]
  cat("Purity data available for", sum(!is.na(qc_meta$purity)), "samples\n")
} else {
  qc_meta$purity <- NA_real_
  cat("No purity data available\n")
}

# =========================
# PCA analysis
# =========================
cat("\n=== PCA ANALYSIS ===\n")

# PCA analysis
pr <- prcomp(t(fpkm_final), scale. = TRUE)
pcs <- as.data.frame(pr$x[, 1:10])

# Add metadata to PCA results
pcs$sample_id <- rownames(pcs)
pcs$tss <- qc_meta$tss
pcs$batch <- qc_meta$batch
pcs$age <- qc_meta$age
pcs$gender <- qc_meta$gender

# Create PCA plots
p1 <- ggplot(pcs, aes(PC1, PC2, color = batch)) +
  geom_point(size = 2) +
  ggtitle("PCA by Batch (TSS)") +
  theme_minimal()

p2 <- ggplot(pcs, aes(PC1, PC2, color = gender)) +
  geom_point(size = 2) +
  ggtitle("PCA by Gender") +
  theme_minimal()

p3 <- ggplot(pcs, aes(PC1, PC2, color = age)) +
  geom_point(size = 2) +
  scale_color_gradient(low = "blue", high = "red") +
  ggtitle("PCA by Age") +
  theme_minimal()

# Save PCA plots
ggsave(file.path(CONFIG$outdir, "PCA_by_Batch.png"), p1, width = 8, height = 6, dpi = 300)
ggsave(file.path(CONFIG$outdir, "PCA_by_Gender.png"), p2, width = 8, height = 6, dpi = 300)
ggsave(file.path(CONFIG$outdir, "PCA_by_Age.png"), p3, width = 8, height = 6, dpi = 300)

cat("PCA plots saved\n")

# =========================
# Save results
# =========================
cat("\n=== SAVING RESULTS ===\n")

# Save processed expression data
saveRDS(fpkm_final, file.path(CONFIG$outdir, "expression_processed.rds"))

# Save metadata
fwrite(qc_meta, file.path(CONFIG$outdir, "metadata_processed.csv"))

# Save PCA results
fwrite(pcs, file.path(CONFIG$outdir, "pca_results.csv"))

# Create summary report
summary_report <- data.frame(
  Metric = c("Original samples", "Primary tumor samples", "After deduplication", 
             "Final samples", "Original genes", "After filtering", "Final genes"),
  Count = c(ncol(fpkm_data), ncol(fpkm_filtered), ncol(fpkm_dedup),
            ncol(fpkm_final), nrow(fpkm_data), nrow(fpkm_detected), nrow(fpkm_final)),
  stringsAsFactors = FALSE
)

fwrite(summary_report, file.path(CONFIG$outdir, "processing_summary.csv"))

# Print summary
cat("\n=== PROCESSING SUMMARY ===\n")
print(summary_report)

# Check for DGAT genes
dgat_genes <- c("DGAT1", "DGAT2")
dgat_present <- dgat_genes %in% rownames(fpkm_final)
cat("\nDGAT gene presence:\n")
for (i in seq_along(dgat_genes)) {
  cat(sprintf("  %s: %s\n", dgat_genes[i], ifelse(dgat_present[i], "✓ Present", "✗ Missing")))
}

# Check expression levels of DGAT genes
if (all(dgat_present)) {
  cat("\nDGAT gene expression levels:\n")
  for (gene in dgat_genes) {
    expr_vals <- fpkm_final[gene, ]
    cat(sprintf("  %s: Mean=%.3f, Median=%.3f, Range=[%.3f, %.3f]\n", 
                gene, mean(expr_vals), median(expr_vals), min(expr_vals), max(expr_vals)))
  }
}

cat("\n=== PIPELINE COMPLETE ===\n")
cat("Results saved to:", CONFIG$outdir, "\n")
cat("Processed expression data:", nrow(fpkm_final), "genes x", ncol(fpkm_final), "samples\n")
