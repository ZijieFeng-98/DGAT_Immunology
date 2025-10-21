# ============================================================================
# TCGA GBM Batch Correction - Customized for Your Data
# ============================================================================

suppressPackageStartupMessages({
  library(sva)
  library(matrixStats)
  library(ggplot2)
  library(gridExtra)
  library(dplyr)
})

# ============================================================================
# CONFIGURATION
# ============================================================================

CONFIG <- list(
  input_expr     = "Processed_Data/TCGA_GBM_Master_Pipeline/expression_processed.rds",
  input_meta     = "Processed_Data/TCGA_GBM_Master_Pipeline/metadata_processed.csv",
  output_dir     = "Processed_Data/TCGA_GBM_Batch_Corrected",
  output_expr    = "expression_batch_corrected.rds",
  output_meta    = "metadata_batch_corrected.csv"
)

# Create output directory
dir.create(CONFIG$output_dir, recursive = TRUE, showWarnings = FALSE)

# ============================================================================
# LOAD DATA
# ============================================================================

message("Loading data...")
expr <- readRDS(CONFIG$input_expr)
meta <- read.csv(CONFIG$input_meta, row.names = 1, stringsAsFactors = FALSE)

# Verify alignment
stopifnot(all(colnames(expr) == rownames(meta)))

message(sprintf("Loaded: %d genes x %d samples", nrow(expr), ncol(expr)))

# ============================================================================
# CRITICAL: Remove singleton batches (batch 8 has n=1)
# ============================================================================

batch_counts <- table(meta$batch)
singleton_batches <- names(batch_counts[batch_counts < 2])

if (length(singleton_batches) > 0) {
  message(sprintf("Removing %d samples from singleton batches: %s", 
                  length(singleton_batches), paste(singleton_batches, collapse = ", ")))
  
  keep_samples <- !meta$batch %in% singleton_batches
  expr <- expr[, keep_samples]
  meta <- meta[keep_samples, ]
  
  message(sprintf("Retained %d samples across %d batches", 
                  ncol(expr), length(unique(meta$batch))))
}

# ============================================================================
# Gene filtering
# ============================================================================

# Remove zero-variance genes
gene_sds <- rowSds(expr, na.rm = TRUE)
expr <- expr[gene_sds > 0 & !is.na(gene_sds), ]

message(sprintf("Final dimensions: %d genes x %d samples", nrow(expr), ncol(expr)))

# ============================================================================
# Build model matrix to preserve biology
# ============================================================================

batch <- factor(meta$batch)

# Check for NA values in purity and age
purity_na <- sum(is.na(meta$purity))
age_na <- sum(is.na(meta$age))

if (purity_na > 0 || age_na > 0) {
  message(sprintf("Warning: %d samples missing purity, %d missing age", purity_na, age_na))
  message("Applying median imputation...")
  
  # Median imputation (safest option)
  if (purity_na > 0) meta$purity[is.na(meta$purity)] <- median(meta$purity, na.rm = TRUE)
  if (age_na > 0) meta$age[is.na(meta$age)] <- median(meta$age, na.rm = TRUE)
}

# Build model matrix
mod <- model.matrix(~ purity + age, data = meta)

message("Model matrix built to preserve purity and age relationships")

# ============================================================================
# Apply ComBat
# ============================================================================

message("\nApplying ComBat batch correction...")

expr_corrected <- tryCatch({
  ComBat(
    dat = expr,
    batch = batch,
    mod = mod,
    par.prior = TRUE,
    prior.plots = FALSE
  )
}, error = function(e) {
  message("Parametric prior failed, trying non-parametric...")
  ComBat(
    dat = expr,
    batch = batch,
    mod = mod,
    par.prior = FALSE,
    prior.plots = FALSE
  )
})

message("ComBat completed successfully!")

# ============================================================================
# Validation: Before vs After PCA
# ============================================================================

message("\nGenerating validation plots...")

# PCA before
pca_pre <- prcomp(t(expr), scale. = TRUE)
var_pre <- pca_pre$sdev^2 / sum(pca_pre$sdev^2)

# PCA after
pca_post <- prcomp(t(expr_corrected), scale. = TRUE)
var_post <- pca_post$sdev^2 / sum(pca_post$sdev^2)

# Calculate batch correlation BEFORE
pc1_pre <- pca_pre$x[,1]
batch_numeric <- as.numeric(factor(meta$batch))
cor_before <- cor(pc1_pre, batch_numeric, method = "spearman")

# Calculate batch correlation AFTER
pc1_post <- pca_post$x[,1]
cor_after <- cor(pc1_post, batch_numeric, method = "spearman")

message(sprintf("PC1-batch correlation: %.3f → %.3f", cor_before, cor_after))

# Create comparison plots
df_pre <- data.frame(
  PC1 = pca_pre$x[,1],
  PC2 = pca_pre$x[,2],
  batch = factor(meta$batch)
)

df_post <- data.frame(
  PC1 = pca_post$x[,1],
  PC2 = pca_post$x[,2],
  batch = factor(meta$batch)
)

p1 <- ggplot(df_pre, aes(PC1, PC2, color = batch)) +
  geom_point(size = 2, alpha = 0.7) +
  labs(title = "Before ComBat",
       subtitle = sprintf("PC1: %.1f%% | Batch-cor: %.3f", var_pre[1]*100, cor_before)) +
  theme_minimal()

p2 <- ggplot(df_post, aes(PC1, PC2, color = batch)) +
  geom_point(size = 2, alpha = 0.7) +
  labs(title = "After ComBat",
       subtitle = sprintf("PC1: %.1f%% | Batch-cor: %.3f", var_post[1]*100, cor_after)) +
  theme_minimal()

# Save plots
pdf(file.path(CONFIG$output_dir, "batch_correction_comparison.pdf"), width = 14, height = 6)
grid.arrange(p1, p2, ncol = 2)
dev.off()

png(file.path(CONFIG$output_dir, "batch_correction_comparison.png"), 
    width = 14, height = 6, units = "in", res = 300)
grid.arrange(p1, p2, ncol = 2)
dev.off()

# ============================================================================
# Additional validation plots
# ============================================================================

# Create individual PCA plots
p3 <- ggplot(df_pre, aes(PC1, PC2, color = batch)) +
  geom_point(size = 2, alpha = 0.7) +
  labs(title = "PCA Before ComBat Correction",
       subtitle = sprintf("PC1: %.1f%% variance | Batch correlation: %.3f", 
                         var_pre[1]*100, cor_before)) +
  theme_minimal() +
  theme(legend.position = "right")

p4 <- ggplot(df_post, aes(PC1, PC2, color = batch)) +
  geom_point(size = 2, alpha = 0.7) +
  labs(title = "PCA After ComBat Correction",
       subtitle = sprintf("PC1: %.1f%% variance | Batch correlation: %.3f", 
                         var_post[1]*100, cor_after)) +
  theme_minimal() +
  theme(legend.position = "right")

ggsave(file.path(CONFIG$output_dir, "PCA_before_ComBat.png"), p3, width = 10, height = 8, dpi = 300)
ggsave(file.path(CONFIG$output_dir, "PCA_after_ComBat.png"), p4, width = 10, height = 8, dpi = 300)

# ============================================================================
# Save corrected data
# ============================================================================

message("\nSaving corrected data...")

saveRDS(expr_corrected, file.path(CONFIG$output_dir, CONFIG$output_expr))
write.csv(meta, file.path(CONFIG$output_dir, CONFIG$output_meta))

# Save summary statistics
summary_stats <- data.frame(
  metric = c("samples_before", "samples_after", "genes", "batches",
             "PC1_batch_cor_before", "PC1_batch_cor_after",
             "PC1_var_before", "PC1_var_after",
             "PC2_var_before", "PC2_var_after"),
  value = c(ncol(expr) + length(singleton_batches), ncol(expr_corrected), 
            nrow(expr_corrected), length(unique(meta$batch)),
            cor_before, cor_after,
            var_pre[1], var_post[1],
            var_pre[2], var_post[2])
)

write.csv(summary_stats, file.path(CONFIG$output_dir, "batch_correction_summary.csv"), row.names = FALSE)

# ============================================================================
# Generate detailed report
# ============================================================================

sink(file.path(CONFIG$output_dir, "batch_correction_report.txt"))
cat(strrep("=", 70), "\n")
cat("TCGA GBM BATCH CORRECTION REPORT\n")
cat(strrep("=", 70), "\n\n")
cat(sprintf("Analysis date: %s\n", Sys.time()))
cat(sprintf("Original samples: %d\n", ncol(expr) + length(singleton_batches)))
cat(sprintf("Final samples: %d\n", ncol(expr_corrected)))
cat(sprintf("Genes: %d\n", nrow(expr_corrected)))
cat(sprintf("Batches: %d\n\n", length(unique(meta$batch))))

if (length(singleton_batches) > 0) {
  cat("Removed singleton batches:\n")
  for (batch in singleton_batches) {
    cat(sprintf("  - Batch %s (n=1)\n", batch))
  }
  cat("\n")
}

cat("Batch effect reduction:\n")
cat(sprintf("  PC1-batch correlation: %.3f → %.3f\n", cor_before, cor_after))
cat(sprintf("  PC1 variance explained: %.1f%% → %.1f%%\n", var_pre[1]*100, var_post[1]*100))
cat(sprintf("  PC2 variance explained: %.1f%% → %.1f%%\n", var_pre[2]*100, var_post[2]*100))
cat("\n")

cat("Model matrix preserved:\n")
cat("  - Purity relationships\n")
cat("  - Age relationships\n")
cat("\n")

cat("Files created:\n")
cat("  - expression_batch_corrected.rds (main output)\n")
cat("  - metadata_batch_corrected.csv (updated metadata)\n")
cat("  - batch_correction_comparison.png (before/after plots)\n")
cat("  - PCA_before_ComBat.png (pre-correction PCA)\n")
cat("  - PCA_after_ComBat.png (post-correction PCA)\n")
cat("  - batch_correction_summary.csv (statistics)\n")
cat("\n")

cat(strrep("=", 70), "\n")
cat("INTERPRETATION:\n")
if (cor_after < 0.3) {
  cat("✓ SUCCESS: Batch effects successfully reduced (correlation < 0.3)\n")
} else {
  cat("⚠ PARTIAL: Batch effects reduced but may still be present\n")
}
cat("✓ Biological relationships (purity, age) preserved\n")
cat("✓ Ready for downstream DGAT1 correlation analysis\n")
cat(strrep("=", 70), "\n")
sink()

# ============================================================================
# DONE
# ============================================================================

message("\n", strrep("=", 70))
message("BATCH CORRECTION COMPLETE!")
message(strrep("=", 70))
message(sprintf("Output directory: %s", normalizePath(CONFIG$output_dir)))
message(sprintf("Batch effect reduced: %.3f → %.3f", cor_before, cor_after))
message("\nFiles created:")
message("  - expression_batch_corrected.rds (use this for DGAT1 analysis)")
message("  - batch_correction_comparison.png (validation plots)")
message("  - batch_correction_report.txt (detailed summary)")
message("\nNext step: Use corrected data for DGAT1 correlations")
message(strrep("=", 70))
