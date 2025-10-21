#!/usr/bin/env Rscript
# =============================================================================
# 00_validate_immune_pipeline.R — Validate Script 02 immune signatures
# - Test known immune checkpoint genes vs immune cell signatures
# - Expected: POSITIVE correlations (high checkpoint → high immune infiltration)
# - If negative correlations appear, investigate signature validity
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
})

# ----------------------------- CONFIG ----------------------------------------
BASE_DIR    <- "/Users/zijiefeng/Desktop/Guo's lab/My_Research/DGAT_Immunology"
RESULTS_DIR <- file.path(BASE_DIR, "Results")
DECONV_DIR  <- file.path(RESULTS_DIR, "Deconvolution", "GBM")
EXPR_FILE   <- file.path(RESULTS_DIR, "Analysis_Framework", "TCGA_GBM", "01_cleaned_data", "tcga_gbm_expression_log2.rds")
SIG_FILE    <- file.path(DECONV_DIR, "immune_signature_scores.csv")

# Known immune checkpoint/activation genes that SHOULD correlate with immune cells
TEST_GENES <- c(
  "CD274",    # PD-L1 (tumor/immune)
  "PDCD1",    # PD-1 (T cells)
  "CTLA4",    # CTLA4 (Tregs/T cells)
  "HAVCR2",   # TIM-3 (exhaustion)
  "IDO1",     # IDO1 (immunosuppression)
  "LAG3",     # LAG3 (exhaustion)
  "CD8A",     # CD8 (T cells)
  "CD4",      # CD4 (T cells)
  "FOXP3",    # Tregs
  "CD68",     # Macrophages
  "CD163"     # M2 Macrophages
)

msg <- function(...) cat(sprintf(...), "\n")

# ----------------------------- LOAD DATA -------------------------------------
msg("\n=== Loading Data ===")
expr <- readRDS(EXPR_FILE)
msg("Expression: %d genes × %d samples", nrow(expr), ncol(expr))

sig_scores <- fread(SIG_FILE, data.table = FALSE)
msg("Signatures: %d samples × %d cell types", nrow(sig_scores), ncol(sig_scores) - 1)

# Check which test genes are present
present_genes <- TEST_GENES[TEST_GENES %in% rownames(expr)]
missing_genes <- TEST_GENES[!TEST_GENES %in% rownames(expr)]

msg("\nTest genes present: %d/%d", length(present_genes), length(TEST_GENES))
if (length(missing_genes) > 0) {
  msg("Missing genes: %s", paste(missing_genes, collapse = ", "))
}

# ----------------------------- CORRELATIONS ----------------------------------
msg("\n=== Computing Correlations ===\n")

# Align samples
common_samples <- intersect(colnames(expr), sig_scores$sample)
expr_aligned <- expr[, common_samples, drop = FALSE]
sig_aligned <- sig_scores[match(common_samples, sig_scores$sample), ]

# Store results
results <- list()

for (gene in present_genes) {
  gene_expr <- as.numeric(expr_aligned[gene, ])
  
  # Get signature columns (exclude 'sample' column)
  sig_cols <- setdiff(names(sig_aligned), "sample")
  
  gene_results <- data.frame(
    gene = character(),
    signature = character(),
    rho = numeric(),
    p_value = numeric(),
    expected_sign = character(),
    validation = character(),
    stringsAsFactors = FALSE
  )
  
  for (sig_name in sig_cols) {
    sig_score <- sig_aligned[[sig_name]]
    
    # Spearman correlation
    cor_test <- cor.test(gene_expr, sig_score, method = "spearman", exact = FALSE)
    rho <- cor_test$estimate
    p_val <- cor_test$p.value
    
    # Expected correlation sign based on gene biology
    expected_sign <- case_when(
      gene %in% c("CD274", "IDO1") ~ "positive",  # Expressed on tumor/immune interface
      gene %in% c("PDCD1", "CTLA4", "HAVCR2", "LAG3") ~ "positive",  # T cell markers
      gene %in% c("CD8A", "CD4", "FOXP3") ~ "positive",  # T cell subset markers
      gene %in% c("CD68", "CD163") ~ "positive",  # Myeloid markers
      TRUE ~ "positive"  # Default: immune genes correlate with immune cells
    )
    
    # Validation status
    validation <- if (expected_sign == "positive" && rho > 0) {
      "✓ PASS"
    } else if (expected_sign == "positive" && rho < 0) {
      "✗ FAIL"
    } else {
      "? UNKNOWN"
    }
    
    gene_results <- rbind(gene_results, data.frame(
      gene = gene,
      signature = sig_name,
      rho = rho,
      p_value = p_val,
      expected_sign = expected_sign,
      validation = validation,
      stringsAsFactors = FALSE
    ))
  }
  
  results[[gene]] <- gene_results
  
  # Print summary for this gene
  sig_count <- sum(gene_results$p_value < 0.05)
  pos_count <- sum(gene_results$rho > 0 & gene_results$p_value < 0.05)
  neg_count <- sum(gene_results$rho < 0 & gene_results$p_value < 0.05)
  fail_count <- sum(gene_results$validation == "✗ FAIL")
  
  msg("%-8s: %2d sig correlations (%2d positive, %2d negative) | %s",
      gene, sig_count, pos_count, neg_count,
      if (fail_count == 0) "✓ ALL PASS" else sprintf("✗ %d FAIL", fail_count))
}

# ----------------------------- COMBINE & SAVE --------------------------------
all_results <- do.call(rbind, results)
rownames(all_results) <- NULL

# Add FDR correction per gene
all_results <- all_results %>%
  group_by(gene) %>%
  mutate(fdr = p.adjust(p_value, method = "fdr")) %>%
  ungroup() %>%
  arrange(gene, p_value)

# Save full results
out_file <- file.path(DECONV_DIR, "validation_checkpoint_correlations.csv")
fwrite(all_results, out_file)
msg("\nSaved full results: %s", out_file)

# ----------------------------- SUMMARY TABLE ---------------------------------
msg("\n=== Validation Summary ===\n")

# Per-gene summary
gene_summary <- all_results %>%
  group_by(gene) %>%
  summarise(
    n_sig = sum(p_value < 0.05),
    n_pos = sum(rho > 0 & p_value < 0.05),
    n_neg = sum(rho < 0 & p_value < 0.05),
    n_fail = sum(validation == "✗ FAIL"),
    max_rho = max(rho),
    min_p = min(p_value),
    status = if_else(n_fail == 0, "✓ PASS", "✗ FAIL")
  ) %>%
  arrange(desc(n_pos))

print(as.data.frame(gene_summary))

# Overall validation
total_tests <- nrow(all_results)
total_fails <- sum(all_results$validation == "✗ FAIL")
pass_rate <- (total_tests - total_fails) / total_tests * 100

msg("\n=== Overall Validation ===")
msg("Total tests: %d", total_tests)
msg("Failures: %d", total_fails)
msg("Pass rate: %.1f%%", pass_rate)

if (pass_rate >= 80) {
  msg("\n✓✓✓ PIPELINE VALIDATED - Immune signatures behave as expected ✓✓✓")
} else if (pass_rate >= 60) {
  msg("\n⚠ WARNING - Some unexpected correlations detected, review required")
} else {
  msg("\n✗✗✗ PIPELINE FAILED - Major issues with immune signatures ✗✗✗")
}

# ----------------------------- TOP HITS --------------------------------------
msg("\n=== Top Positive Correlations (Expected Pattern) ===\n")
top_positive <- all_results %>%
  filter(rho > 0, p_value < 0.05) %>%
  arrange(desc(abs(rho))) %>%
  head(20)

print(as.data.frame(select(top_positive, gene, signature, rho, p_value, fdr, validation)))

msg("\n=== Unexpected Negative Correlations (Review Required) ===\n")
unexpected_neg <- all_results %>%
  filter(validation == "✗ FAIL") %>%
  arrange(p_value)

if (nrow(unexpected_neg) > 0) {
  print(as.data.frame(select(unexpected_neg, gene, signature, rho, p_value, fdr)))
  msg("\n⚠ %d unexpected negative correlations found - investigate these!", nrow(unexpected_neg))
} else {
  msg("None found - all correlations match expected patterns! ✓")
}

# ----------------------------- VISUALIZATION ---------------------------------
# Create a heatmap of gene-signature correlations
msg("\n=== Creating Validation Heatmap ===")

# Reshape to wide format for heatmap
cor_matrix <- all_results %>%
  select(gene, signature, rho) %>%
  tidyr::pivot_wider(names_from = signature, values_from = rho) %>%
  tibble::column_to_rownames("gene") %>%
  as.matrix()

# Simple heatmap using base R
png(file.path(DECONV_DIR, "validation_checkpoint_heatmap.png"), 
    width = 1200, height = 800, res = 120)
par(mar = c(10, 10, 4, 2))

# Color palette: blue (negative) to white (0) to red (positive)
col_breaks <- seq(-1, 1, length.out = 101)
col_palette <- colorRampPalette(c("blue", "white", "red"))(100)

image(t(cor_matrix), 
      col = col_palette,
      axes = FALSE,
      main = "Validation: Checkpoint Genes vs Immune Signatures\n(Expected: Positive Correlations)",
      cex.main = 1.2)

# Add axes
axis(1, at = seq(0, 1, length.out = ncol(cor_matrix)), 
     labels = colnames(cor_matrix), las = 2, cex.axis = 0.8)
axis(2, at = seq(0, 1, length.out = nrow(cor_matrix)), 
     labels = rownames(cor_matrix), las = 2, cex.axis = 0.8)

# Add color scale legend
legend("topright", 
       legend = c("Strong Neg", "Weak", "Strong Pos"), 
       fill = c("blue", "white", "red"),
       title = "Correlation",
       cex = 0.8)

dev.off()
msg("Saved heatmap: %s", file.path(DECONV_DIR, "validation_checkpoint_heatmap.png"))

msg("\n=== Validation Complete ===\n")

