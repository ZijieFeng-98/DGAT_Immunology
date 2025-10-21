#!/usr/bin/env Rscript
# =============================================================================
# 03_purity_adjustment_test.R — Test if immune correlations are purity artifacts
# Critical analysis: Do DGAT1-immune associations disappear after purity control?
# =============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(ppcor)
  library(ggplot2)
  library(readr)
})

# ====== CONFIG =================================================================
BASE_DIR    <- "/Users/zijiefeng/Desktop/Guo's lab/My_Research/DGAT_Immunology"
RESULTS_DIR <- file.path(BASE_DIR, "Results")
DECONV_DIR  <- file.path(RESULTS_DIR, "Deconvolution")
OUT_DIR     <- file.path(RESULTS_DIR, "Purity_Analysis")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

CANCER_TYPE <- "GBM"
DGAT_GENE   <- "DGAT1"

# =============================================================================
# ANALYSIS FUNCTIONS
# =============================================================================

msg <- function(...) cat(paste0(..., "\n"))

# Compute ESTIMATE purity
compute_estimate_purity <- function(Xlog) {
  msg("Computing ESTIMATE purity scores...")
  
  if (!requireNamespace("estimate", quietly = TRUE)) {
    stop(
      "estimate package required but not installed.\n",
      "Install from: BiocManager::install('estimate')\n"
    )
  }
  suppressPackageStartupMessages(library(estimate))
  
  # Remove duplicate gene names first
  dup_genes <- duplicated(rownames(Xlog))
  if (any(dup_genes)) {
    msg("  Removing ", sum(dup_genes), " duplicate genes...")
    Xlog <- Xlog[!dup_genes, , drop = FALSE]
  }
  
  # ESTIMATE expects data.frame with Description column
  # Reverse log2(x+1) transform to get approximate counts/FPKM
  approx_expr <- pmax(2^Xlog - 1, 0)
  
  # Convert to data.frame format ESTIMATE expects
  estimate_input <- data.frame(
    Description = rownames(approx_expr),
    approx_expr,
    row.names = rownames(approx_expr),
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  
  # Use temp files for ESTIMATE workflow
  tmp_in  <- tempfile(fileext = ".txt")
  tmp_out <- tempfile(fileext = ".gct")
  tmp_fil <- tempfile(fileext = ".gct")
  
  # Write input file
  write.table(estimate_input, tmp_in, sep="\t", quote=FALSE, row.names=TRUE, col.names=NA)
  
  # Filter to common genes and calculate scores
  estimate::filterCommonGenes(input.f = tmp_in, output.f = tmp_fil, id = "GeneSymbol")
  estimate::estimateScore(tmp_fil, tmp_out, platform = "illumina")
  
  # Read GCT output file manually
  # GCT format: #1.2 header, then N M dimensions, then header row, then data
  lines <- readLines(tmp_out)
  # Skip first 2 lines (#1.2 and dimensions)
  header_idx <- 3
  data_start <- 4
  
  # Parse header
  header <- unlist(strsplit(lines[header_idx], "\t"))
  sample_names <- header[3:length(header)]  # Skip NAME and Description
  
  # Parse data rows
  data_lines <- lines[data_start:length(lines)]
  score_list <- lapply(data_lines, function(line) {
    parts <- unlist(strsplit(line, "\t"))
    row_name <- parts[1]
    values <- as.numeric(parts[3:length(parts)])
    list(name = row_name, values = values)
  })
  
  # Build matrix
  score_matrix <- do.call(rbind, lapply(score_list, function(x) x$values))
  rownames(score_matrix) <- sapply(score_list, function(x) x$name)
  colnames(score_matrix) <- sample_names
  
  msg("  ESTIMATE scores calculated:")
  msg("    Row names: ", paste(rownames(score_matrix), collapse=", "))
  
  # Extract purity - ESTIMATE uses different row names depending on version
  # Try multiple possible names
  purity_row <- NULL
  purity_candidates <- c("TumorPurity", "Purity", "ESTIMATE.Purity", "tumor_purity")
  for (cand in purity_candidates) {
    if (cand %in% rownames(score_matrix)) {
      purity_row <- cand
      break
    }
  }
  
  # If no purity row, calculate from ESTIMATE score
  # Purity formula: cos(0.6049872018 + 0.0001467884*ESTIMATEScore)
  if (is.null(purity_row)) {
    if ("ESTIMATEScore" %in% rownames(score_matrix)) {
      msg("  Calculating purity from ESTIMATEScore...")
      estimate_scores <- as.numeric(score_matrix["ESTIMATEScore", ])
      purity <- cos(0.6049872018 + 0.0001467884 * estimate_scores)
      names(purity) <- colnames(score_matrix)
    } else {
      warning("Cannot calculate purity - no ESTIMATEScore row found.")
      return(NULL)
    }
  } else {
    purity <- as.numeric(score_matrix[purity_row, ])
    names(purity) <- colnames(score_matrix)
  }
  
  purity
}

# Load all data
load_data <- function(ct) {
  msg("Loading data for ", ct, "...")
  
  # Expression data
  expr_file <- file.path(RESULTS_DIR, "Analysis_Framework", "TCGA_GBM", 
                         "01_cleaned_data", "tcga_gbm_expression_log2.rds")
  if (!file.exists(expr_file)) {
    stop("Expression file not found: ", expr_file)
  }
  expr <- readRDS(expr_file)
  
  # Purity estimates - compute on the fly using ESTIMATE
  msg("Purity estimates not found - computing using ESTIMATE...")
  purity_vec <- compute_estimate_purity(expr)
  
  purity_df <- data.frame(
    sample = names(purity_vec),
    purity_estimate = as.numeric(purity_vec),
    stringsAsFactors = FALSE
  )
  
  # Remove any non-sample rows (like "Description" header artifacts)
  purity_df <- purity_df[!grepl("^Description", purity_df$sample), , drop = FALSE]
  
  # Standardize sample names: dots to dashes for TCGA format
  purity_df$sample <- gsub("\\.", "-", purity_df$sample)
  
  # Save for future use
  purity_out <- file.path(RESULTS_DIR, "Analysis_Framework", "TCGA_GBM",
                          "01_cleaned_data", "tcga_gbm_purity_estimates.csv")
  write_csv(purity_df, purity_out)
  msg("Saved purity estimates to: ", purity_out)
  
  # Immune signature scores
  immune_file <- file.path(DECONV_DIR, ct, "immune_signature_scores.csv")
  if (!file.exists(immune_file)) {
    stop("Immune scores not found. Run Script 02 first.")
  }
  immune_df <- read_csv(immune_file, show_col_types = FALSE)
  
  list(expr = expr, purity = purity_df, immune = immune_df)
}

# Prepare merged dataset
prepare_dataset <- function(data_list) {
  msg("Preparing merged dataset...")
  
  expr <- data_list$expr
  purity <- data_list$purity
  immune <- data_list$immune
  
  # Extract DGAT1 expression
  if (!(DGAT_GENE %in% rownames(expr))) {
    stop(DGAT_GENE, " not found in expression matrix")
  }
  
  dgat1_vec <- expr[DGAT_GENE, ]
  dgat1_df <- data.frame(
    sample = names(dgat1_vec),
    dgat1 = as.numeric(dgat1_vec)
  )
  
  # Merge all data
  merged <- immune %>%
    inner_join(dgat1_df, by = "sample") %>%
    inner_join(purity %>% dplyr::select(sample, purity_estimate), by = "sample") %>%
    filter(!is.na(purity_estimate))
  
  msg("  Samples with all data: ", nrow(merged))
  msg("  Immune signatures: ", ncol(immune) - 1)
  
  return(merged)
}

# Raw correlations
compute_raw_correlations <- function(df, signatures) {
  msg("\nComputing raw correlations (DGAT1 vs immune)...")
  
  results <- lapply(signatures, function(sig) {
    if (!(sig %in% names(df))) return(NULL)
    
    # Remove NAs
    valid <- !is.na(df[[sig]]) & !is.na(df$dgat1)
    x <- df[[sig]][valid]
    y <- df$dgat1[valid]
    
    if (length(x) < 10) return(NULL)
    
    # Spearman correlation
    test <- cor.test(x, y, method = "spearman")
    
    data.frame(
      Signature = sig,
      rho_raw = test$estimate,
      p_raw = test$p.value,
      n = length(x)
    )
  })
  
  do.call(rbind, results[!sapply(results, is.null)])
}

# Partial correlations (controlling for purity)
compute_partial_correlations <- function(df, signatures) {
  msg("\nComputing partial correlations (controlling for purity)...")
  
  results <- lapply(signatures, function(sig) {
    if (!(sig %in% names(df))) return(NULL)
    
    # Remove NAs
    valid <- !is.na(df[[sig]]) & !is.na(df$dgat1) & !is.na(df$purity_estimate)
    
    if (sum(valid) < 10) return(NULL)
    
    x <- df[[sig]][valid]
    y <- df$dgat1[valid]
    z <- df$purity_estimate[valid]
    
    # Partial Spearman
    test <- try(
      pcor.test(x, y, z, method = "spearman"),
      silent = TRUE
    )
    
    if (inherits(test, "try-error")) return(NULL)
    
    data.frame(
      Signature = sig,
      rho_partial = test$estimate,
      p_partial = test$p.value,
      n = sum(valid)
    )
  })
  
  do.call(rbind, results[!sapply(results, is.null)])
}

# Test purity-DGAT1 correlation
test_purity_dgat_correlation <- function(df) {
  msg("\n=== CRITICAL TEST: Does DGAT1 correlate with purity? ===")
  
  valid <- !is.na(df$dgat1) & !is.na(df$purity_estimate)
  
  test <- cor.test(
    df$dgat1[valid], 
    df$purity_estimate[valid], 
    method = "spearman"
  )
  
  msg("  Spearman rho: ", round(test$estimate, 3))
  msg("  P-value: ", format(test$p.value, scientific = TRUE))
  
  if (abs(test$estimate) > 0.3 & test$p.value < 0.05) {
    msg("\n  ⚠️ WARNING: Strong DGAT1-purity correlation detected!")
    msg("  This suggests purity confounding is likely.")
  } else if (abs(test$estimate) > 0.15 & test$p.value < 0.05) {
    msg("\n  ⚠️ CAUTION: Moderate DGAT1-purity correlation.")
    msg("  Partial correlations are recommended.")
  } else {
    msg("\n  ✓ DGAT1-purity correlation is weak.")
    msg("  Purity confounding is unlikely.")
  }
  
  # Plot
  p <- ggplot(df, aes(purity_estimate, dgat1)) +
    geom_point(alpha = 0.5, size = 2) +
    geom_smooth(method = "lm", color = "red") +
    labs(
      title = "DGAT1 vs Tumor Purity",
      subtitle = sprintf("Spearman ρ = %.3f, p = %.2e", 
                        test$estimate, test$p.value),
      x = "Tumor Purity (ESTIMATE)",
      y = paste(DGAT_GENE, "Expression (log2)")
    ) +
    theme_minimal(base_size = 12)
  
  ggsave(
    file.path(OUT_DIR, "dgat1_vs_purity.png"),
    p, width = 7, height = 6, dpi = 300
  )
  
  return(test)
}

# Compare raw vs partial correlations
compare_correlations <- function(raw_cor, partial_cor) {
  msg("\n=== Comparing raw vs purity-adjusted correlations ===")
  
  # Merge
  comparison <- raw_cor %>%
    inner_join(partial_cor, by = "Signature") %>%
    mutate(
      fdr_raw = p.adjust(p_raw, method = "BH"),
      fdr_partial = p.adjust(p_partial, method = "BH"),
      delta_rho = rho_partial - rho_raw,
      sign_change = sign(rho_raw) != sign(rho_partial),
      sig_raw = fdr_raw < 0.05,
      sig_partial = fdr_partial < 0.05,
      sig_status = case_when(
        sig_raw & sig_partial ~ "Both significant",
        sig_raw & !sig_partial ~ "Lost after adjustment",
        !sig_raw & sig_partial ~ "Gained after adjustment",
        TRUE ~ "Neither significant"
      )
    ) %>%
    arrange(fdr_raw)
  
  # Summary
  msg("\n  Total signatures tested: ", nrow(comparison))
  msg("  Significant (raw): ", sum(comparison$sig_raw))
  msg("  Significant (purity-adjusted): ", sum(comparison$sig_partial))
  msg("  Lost significance: ", sum(comparison$sig_status == "Lost after adjustment"))
  msg("  Gained significance: ", sum(comparison$sig_status == "Gained after adjustment"))
  msg("  Sign changes: ", sum(comparison$sign_change))
  
  # Key findings
  msg("\n  Top 5 correlations that LOST significance:")
  lost <- comparison %>%
    filter(sig_status == "Lost after adjustment") %>%
    head(5)
  
  if (nrow(lost) > 0) {
    for (i in 1:nrow(lost)) {
      msg(sprintf("    %s: ρ %.3f → %.3f (FDR %.3f → %.3f)",
                  lost$Signature[i],
                  lost$rho_raw[i],
                  lost$rho_partial[i],
                  lost$fdr_raw[i],
                  lost$fdr_partial[i]))
    }
  } else {
    msg("    (None)")
  }
  
  msg("\n  Correlations that REMAIN significant:")
  remain <- comparison %>%
    filter(sig_status == "Both significant")
  
  if (nrow(remain) > 0) {
    for (i in 1:min(5, nrow(remain))) {
      msg(sprintf("    %s: ρ %.3f → %.3f (both FDR < 0.05)",
                  remain$Signature[i],
                  remain$rho_raw[i],
                  remain$rho_partial[i]))
    }
  } else {
    msg("    (None - all associations are purity artifacts!)")
  }
  
  return(comparison)
}

# Visualization: raw vs partial correlations
plot_comparison <- function(comparison) {
  msg("\nGenerating comparison plots...")
  
  # Scatter plot: raw vs partial
  p1 <- ggplot(comparison, aes(rho_raw, rho_partial, color = sig_status)) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
    geom_point(size = 3, alpha = 0.7) +
    geom_text(
      data = comparison %>% filter(abs(rho_raw) > 0.25 | abs(rho_partial) > 0.25),
      aes(label = Signature),
      vjust = -0.5, hjust = 0.5, size = 3, check_overlap = TRUE
    ) +
    scale_color_manual(
      values = c(
        "Both significant" = "#D73027",
        "Lost after adjustment" = "#FDB863",
        "Neither significant" = "grey60",
        "Gained after adjustment" = "#4575B4"
      )
    ) +
    labs(
      title = "Raw vs Purity-Adjusted Correlations",
      subtitle = "DGAT1 vs Immune Signatures",
      x = "Raw Spearman ρ",
      y = "Purity-adjusted Spearman ρ",
      color = "Significance Status"
    ) +
    theme_minimal(base_size = 12) +
    theme(legend.position = "bottom")
  
  ggsave(
    file.path(OUT_DIR, "raw_vs_partial_scatter.png"),
    p1, width = 9, height = 8, dpi = 300
  )
  
  # Heatmap: before and after
  top_sigs <- comparison %>%
    filter(sig_raw | sig_partial) %>%
    pull(Signature)
  
  if (length(top_sigs) > 0) {
    heatmap_data <- comparison %>%
      filter(Signature %in% top_sigs) %>%
      dplyr::select(Signature, rho_raw, rho_partial) %>%
      pivot_longer(-Signature, names_to = "Type", values_to = "Correlation") %>%
      mutate(Type = recode(Type, 
                          rho_raw = "Raw",
                          rho_partial = "Purity-adjusted"))
    
    p2 <- ggplot(heatmap_data, aes(Type, Signature, fill = Correlation)) +
      geom_tile(color = "white") +
      geom_text(aes(label = sprintf("%.2f", Correlation)), size = 3) +
      scale_fill_gradient2(
        low = "#4575B4", mid = "white", high = "#D73027",
        midpoint = 0, limits = c(-0.5, 0.5)
      ) +
      labs(
        title = "Effect of Purity Adjustment",
        subtitle = "Significant associations only"
      ) +
      theme_minimal(base_size = 11) +
      theme(axis.text.x = element_text(angle = 0, hjust = 0.5))
    
    ggsave(
      file.path(OUT_DIR, "purity_adjustment_heatmap.png"),
      p2, width = 6, height = max(4, length(top_sigs) * 0.3), dpi = 300
    )
  }
}

# Generate final interpretation
interpret_results <- function(purity_test, comparison) {
  msg("\n" , paste(rep("=", 70), collapse = ""))
  msg("FINAL INTERPRETATION")
  msg(paste(rep("=", 70), collapse = ""))
  
  purity_rho <- abs(purity_test$estimate)
  lost_pct <- mean(comparison$sig_status == "Lost after adjustment") * 100
  remain_n <- sum(comparison$sig_status == "Both significant")
  
  if (purity_rho > 0.3 & lost_pct > 50) {
    msg("\n🔴 CONCLUSION: Strong purity confounding detected")
    msg("\n  Evidence:")
    msg("  - DGAT1 correlates strongly with purity (|ρ| = ", round(purity_rho, 2), ")")
    msg("  - ", round(lost_pct, 1), "% of associations lost after adjustment")
    msg("  - Only ", remain_n, " associations remain significant")
    msg("\n  Interpretation:")
    msg("  Your negative DGAT1-immune correlations are largely artifacts of")
    msg("  tumor purity. DGAT1-high tumors have higher purity (fewer stromal/")
    msg("  immune cells), creating spurious negative correlations.")
    msg("\n  Recommendation:")
    msg("  - Report purity-adjusted results as primary findings")
    msg("  - Discuss tumor purity as major confounder in methods/discussion")
    msg("  - Consider stratifying by purity tertiles for sensitivity analysis")
    
  } else if (purity_rho > 0.15 & lost_pct > 30) {
    msg("\n🟡 CONCLUSION: Moderate purity confounding")
    msg("\n  Evidence:")
    msg("  - DGAT1 moderately correlates with purity (ρ = ", round(purity_rho, 2), ")")
    msg("  - ", round(lost_pct, 1), "% of associations weakened after adjustment")
    msg("  - ", remain_n, " associations remain robust")
    msg("\n  Interpretation:")
    msg("  Purity plays a role but doesn't fully explain DGAT1-immune associations.")
    msg("  Some effects are real biology, others are confounded.")
    msg("\n  Recommendation:")
    msg("  - Report both raw and adjusted results")
    msg("  - Focus on associations that survive adjustment")
    msg("  - Acknowledge purity as partial confounder")
    
  } else {
    msg("\n🟢 CONCLUSION: Minimal purity confounding")
    msg("\n  Evidence:")
    msg("  - DGAT1-purity correlation is weak (ρ = ", round(purity_rho, 2), ")")
    msg("  - Only ", round(lost_pct, 1), "% of associations changed")
    msg("  - ", remain_n, " associations remain significant")
    msg("\n  Interpretation:")
    msg("  Your DGAT1-immune correlations reflect real biology, not purity artifacts.")
    msg("  DGAT1-high tumors genuinely have different immune landscapes.")
    msg("\n  Recommendation:")
    msg("  - This is a real biological finding!")
    msg("  - Investigate WHY DGAT1-high tumors are immune-excluded")
    msg("  - Consider molecular subtype analysis (mesenchymal vs proneural)")
  }
  
  msg("\n" , paste(rep("=", 70), collapse = ""))
}

# =============================================================================
# MAIN
# =============================================================================

main <- function() {
  msg("\n#####################################################################")
  msg("  PURITY CONFOUNDING ANALYSIS")
  msg("  Testing if DGAT1-immune associations are purity artifacts")
  msg("#####################################################################\n")
  
  # Load data
  data <- load_data(CANCER_TYPE)
  df <- prepare_dataset(data)
  
  # Get signature names
  signatures <- setdiff(names(df), c("sample", "dgat1", "purity_estimate"))
  
  # Test 1: DGAT1 vs purity correlation
  purity_test <- test_purity_dgat_correlation(df)
  
  # Test 2: Raw correlations
  raw_cor <- compute_raw_correlations(df, signatures)
  
  # Test 3: Partial correlations
  partial_cor <- compute_partial_correlations(df, signatures)
  
  # Test 4: Compare
  comparison <- compare_correlations(raw_cor, partial_cor)
  
  # Save results
  write_csv(comparison, file.path(OUT_DIR, "raw_vs_partial_comparison.csv"))
  write_csv(
    data.frame(
      metric = c("dgat1_purity_rho", "dgat1_purity_p"),
      value = c(purity_test$estimate, purity_test$p.value)
    ),
    file.path(OUT_DIR, "purity_test_results.csv")
  )
  
  # Visualize
  plot_comparison(comparison)
  
  # Final interpretation
  interpret_results(purity_test, comparison)
  
  msg("\n✓ Analysis complete. Results saved to: ", OUT_DIR)
  msg("\nKey files:")
  msg("  - raw_vs_partial_comparison.csv")
  msg("  - dgat1_vs_purity.png")
  msg("  - raw_vs_partial_scatter.png")
  msg("  - purity_adjustment_heatmap.png")
}

# Run
main()

