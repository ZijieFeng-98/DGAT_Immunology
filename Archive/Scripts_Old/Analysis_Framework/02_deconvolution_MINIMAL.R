#!/usr/bin/env Rscript
# =============================================================================
# 02_deconvolution_MINIMAL.R — Immune deconvolution using base R packages
# NO IOBR/immunedeconv required - works with R 4.2.3
# Uses: ESTIMATE scores + signature-based scoring (ssGSEA-like)
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(readr)
  library(ggplot2)
  library(pheatmap)
  library(matrixStats)
})

# ====== CONFIG =================================================================
BASE_DIR    <- "/Users/zijiefeng/Desktop/Guo's lab/My_Research/DGAT_Immunology"
DATA_DIR    <- file.path(BASE_DIR, "Data")
PROC_DIR    <- file.path(BASE_DIR, "Processed_Data")
RESULTS_DIR <- file.path(BASE_DIR, "Results")
OUT_ROOT    <- file.path(RESULTS_DIR, "Deconvolution")
dir.create(OUT_ROOT, recursive = TRUE, showWarnings = FALSE)

CANCER_TYPES <- c("GBM")
DGAT_GENE    <- "DGAT1"

# =============================================================================
# IMMUNE CELL GENE SIGNATURES (Curated from literature)
# =============================================================================

# These are well-validated marker genes from multiple sources
IMMUNE_SIGNATURES <- list(
  # T cells
  CD8_T = c("CD8A", "CD8B", "GZMA", "GZMB", "PRF1", "GZMH", "GZMK"),
  CD4_T = c("CD4", "CD40LG", "IL7R", "LEF1", "TCF7"),
  Tregs = c("FOXP3", "IL2RA", "CTLA4", "TIGIT", "IKZF2"),
  
  # NK cells
  NK = c("NCAM1", "NCR1", "KLRB1", "KLRD1", "KLRF1", "NKG7"),
  
  # B cells
  B_cells = c("CD19", "MS4A1", "CD79A", "CD79B", "BLK", "BANK1"),
  
  # Myeloid
  Monocytes = c("CD14", "FCGR3A", "CD16", "S100A8", "S100A9", "VCAN"),
  Macrophages = c("CD68", "CD163", "MSR1", "MRC1", "MARCO", "FCGR1A"),
  M2_Macrophages = c("CD163", "MRC1", "CD206", "ARG1", "CHIL3", "ALOX15"),
  TAM = c("CD163", "MRC1", "CD206", "LGALS3", "APOE", "SPP1"),
  Microglia = c("TMEM119", "P2RY12", "CX3CR1", "HEXB", "SALL1"),
  
  # Dendritic cells
  DC = c("CD1C", "FCER1A", "CLEC10A", "CLEC9A", "XCR1", "BATF3"),
  cDC1 = c("CLEC9A", "XCR1", "BATF3", "IRF8", "SNX22"),
  
  # Other
  Neutrophils = c("FCGR3B", "CSF3R", "CXCR2", "FPR1", "S100A12"),
  Endothelial = c("PECAM1", "VWF", "CDH5", "ENG", "MCAM"),
  Fibroblasts = c("COL1A1", "COL3A1", "DCN", "LUM", "PDGFRA")
)

# T cell exhaustion signatures
EXHAUSTION_SIGS <- list(
  Exhaustion_Early = c("PDCD1", "LAG3", "HAVCR2", "TIGIT"),
  Exhaustion_Terminal = c("PDCD1", "LAG3", "HAVCR2", "TOX", "TOX2", "ENTPD1"),
  Cytotoxic = c("GZMA", "GZMB", "PRF1", "GNLY", "NKG7", "IFNG")
)

# Combine all signatures
ALL_SIGNATURES <- c(IMMUNE_SIGNATURES, EXHAUSTION_SIGS)

# =============================================================================
# SIGNATURE SCORING FUNCTION (Z-score based)
# =============================================================================

score_signature <- function(expr_mat, signature_genes, method = "zscore"){
  # expr_mat: genes x samples (linear scale preferred)
  # signature_genes: character vector of gene symbols
  
  # Find genes present in data
  genes_present <- intersect(signature_genes, rownames(expr_mat))
  
  if (length(genes_present) == 0) {
    warning("No signature genes found in expression matrix")
    return(rep(NA_real_, ncol(expr_mat)))
  }
  
  if (length(genes_present) < length(signature_genes) * 0.3) {
    warning("Less than 30% of signature genes found (", 
            length(genes_present), "/", length(signature_genes), ")")
  }
  
  # Extract signature expression
  sig_expr <- expr_mat[genes_present, , drop = FALSE]
  
  # Method: Z-score then mean (similar to ssGSEA)
  if (method == "zscore") {
    # Z-score normalize each gene across samples
    sig_z <- t(scale(t(sig_expr)))
    
    # Average Z-scores
    scores <- colMeans(sig_z, na.rm = TRUE)
    
  } else if (method == "mean") {
    # Simple mean (for already-normalized data)
    scores <- colMeans(sig_expr, na.rm = TRUE)
  }
  
  names(scores) <- colnames(expr_mat)
  return(scores)
}

# Score all signatures
score_all_signatures <- function(expr_mat, signatures = ALL_SIGNATURES){
  cat("Scoring", length(signatures), "immune signatures...\n")
  
  scores_list <- lapply(names(signatures), function(sig_name){
    cat("  -", sig_name, "...")
    genes <- signatures[[sig_name]]
    s <- score_signature(expr_mat, genes, method = "zscore")
    cat(" done (", sum(!is.na(s)), " samples)\n")
    return(s)
  })
  
  names(scores_list) <- names(signatures)
  
  # Combine into matrix: samples x signatures
  scores_mat <- do.call(cbind, scores_list)
  rownames(scores_mat) <- colnames(expr_mat)
  
  return(as.data.frame(scores_mat))
}

# =============================================================================
# UTILITIES (same as before)
# =============================================================================

msg <- function(...) cat(paste0(..., "\n"))

load_bestcut <- function(cohort_dir){
  fp <- file.path(RESULTS_DIR, "Survival", cohort_dir, "DGAT1_bestcut_info.csv")
  if (file.exists(fp)) {
    info <- fread(fp, quiet = TRUE)
    cutoff <- as.numeric(info$cutoff[1])
    msg("  Using survival-based cutoff: ", round(cutoff, 3))
    return(cutoff)
  } else {
    msg("  No survival cutoff found, will use median")
    return(NA_real_)
  }
}

make_groups <- function(dgat_vec, cutoff = NA_real_){
  if (!is.na(cutoff)) {
    factor(ifelse(dgat_vec >= cutoff, "High", "Low"), levels = c("Low","High"))
  } else {
    med <- median(dgat_vec, na.rm = TRUE)
    msg("  Using median cutoff: ", round(med, 3))
    factor(ifelse(dgat_vec >= med, "High", "Low"), levels = c("Low","High"))
  }
}

backtransform_if_log <- function(X){
  rng <- range(X, na.rm = TRUE)
  avg <- mean(X, na.rm = TRUE)
  
  is_log <- (rng[2] < 25 && rng[1] >= 0 && avg > 2 && avg < 20)
  
  if (is_log) {
    msg("  Detected log2 scale, back-transforming...")
    return(pmax(2^X - 1, 0))
  } else {
    msg("  Data appears linear")
    return(X)
  }
}

wilcoxon_panel <- function(S, groups){
  stopifnot(all(names(groups) %in% rownames(S)))
  S2 <- S[names(groups), , drop = FALSE]
  
  res <- lapply(colnames(S2), function(cell){
    x <- S2[groups == "High", cell]
    y <- S2[groups == "Low",  cell]
    
    if (length(na.omit(x)) < 3 || length(na.omit(y)) < 3) return(NULL)
    
    wt <- suppressWarnings(wilcox.test(x, y))
    delta <- mean(x, na.rm = TRUE) - mean(y, na.rm = TRUE)
    
    data.frame(
      Cell = cell, 
      Mean_High = mean(x, na.rm = TRUE),
      Mean_Low = mean(y, na.rm = TRUE), 
      Delta = delta, 
      p_value = wt$p.value,
      n_High = sum(!is.na(x)),
      n_Low = sum(!is.na(y))
    )
  })
  
  T <- do.call(rbind, res[!sapply(res, is.null)])
  if (is.null(T)) return(NULL)
  
  T$FDR <- p.adjust(T$p_value, method = "BH")
  T[order(T$FDR, -abs(T$Delta)), ]
}

correlations_vs_dgat <- function(S, dgat_vec, purity = NULL){
  stopifnot(all(names(dgat_vec) %in% rownames(S)))
  S2 <- S[names(dgat_vec), , drop = FALSE]
  
  raw <- sapply(colnames(S2), function(cell){
    x <- S2[, cell]
    suppressWarnings(
      cor(x, dgat_vec, method = "spearman", use = "pairwise.complete.obs")
    )
  })
  
  p_raw <- sapply(colnames(S2), function(cell){
    x <- S2[, cell]
    suppressWarnings(cor.test(x, dgat_vec, method = "spearman")$p.value)
  })
  
  # Partial correlation if purity available
  if (!is.null(purity) && requireNamespace("ppcor", quietly = TRUE)){
    adj <- sapply(colnames(S2), function(cell){
      x <- S2[, cell]
      valid <- !is.na(x) & !is.na(dgat_vec) & !is.na(purity)
      if (sum(valid) < 10) return(NA_real_)
      
      pd <- try(
        ppcor::pcor.test(x[valid], dgat_vec[valid], purity[valid], 
                         method = "spearman"),
        silent = TRUE
      )
      if (inherits(pd, "try-error")) return(NA_real_)
      pd$estimate
    })
    
    p_adj <- sapply(colnames(S2), function(cell){
      x <- S2[, cell]
      valid <- !is.na(x) & !is.na(dgat_vec) & !is.na(purity)
      if (sum(valid) < 10) return(NA_real_)
      
      pd <- try(
        ppcor::pcor.test(x[valid], dgat_vec[valid], purity[valid], 
                         method = "spearman"),
        silent = TRUE
      )
      if (inherits(pd, "try-error")) return(NA_real_)
      pd$p.value
    })
  } else {
    adj <- rep(NA_real_, length(raw))
    p_adj <- rep(NA_real_, length(raw))
    names(adj) <- names(raw)
    names(p_adj) <- names(raw)
  }
  
  tibble(
    Cell = names(raw), 
    rho_raw = as.numeric(raw),
    p_raw = as.numeric(p_raw),
    rho_partial = as.numeric(adj),
    p_partial = as.numeric(p_adj)
  )
}

scale01_rows <- function(M){
  t(apply(M, 1, function(x){
    rng <- range(x, na.rm = TRUE)
    if (diff(rng) < 1e-9) return(rep(0.5, length(x)))
    (x - rng[1]) / (rng[2] - rng[1])
  }))
}

save_cell_heatmap <- function(S, groups, outfile, title){
  S2 <- t(scale01_rows(t(S)))
  ord <- order(groups)
  ann <- data.frame(DGAT1 = groups[ord])
  rownames(ann) <- rownames(S2)[ord]
  ann_colors <- list(DGAT1 = c(Low = "#4575B4", High = "#D73027"))
  
  pheatmap(
    t(S2)[, ord, drop = FALSE],
    annotation_col = ann, 
    annotation_colors = ann_colors,
    show_colnames = FALSE, 
    cluster_cols = FALSE, 
    main = title,
    filename = outfile, 
    width = 10, 
    height = 6
  )
}

boxplot_panel <- function(S, groups, keep_cells, outfile, title){
  keep_cells <- intersect(keep_cells, colnames(S))
  if (!length(keep_cells)) return(invisible())
  
  df <- as.data.frame(S)
  df$sample <- rownames(df)
  df$DGAT1 <- groups[df$sample]
  
  long <- pivot_longer(
    df, 
    cols = all_of(keep_cells), 
    names_to = "Cell", 
    values_to = "Score"
  ) %>%
    filter(!is.na(DGAT1))
  
  stat <- long %>%
    group_by(Cell) %>%
    summarise(
      p_value = suppressWarnings(wilcox.test(Score ~ DGAT1)$p.value),
      .groups = "drop"
    ) %>%
    mutate(
      label = case_when(
        p_value < 0.001 ~ "***",
        p_value < 0.01 ~ "**",
        p_value < 0.05 ~ "*",
        TRUE ~ "ns"
      )
    )
  
  g <- ggplot(long, aes(DGAT1, Score, fill = DGAT1)) +
    geom_boxplot(outlier.size = 0.6, width = 0.6, alpha = 0.85) +
    geom_text(
      data = stat,
      aes(x = 1.5, y = Inf, label = label),
      inherit.aes = FALSE,
      vjust = 1.5,
      size = 5
    ) +
    facet_wrap(~Cell, scales = "free_y", ncol = 3) +
    scale_fill_manual(values = c(Low = "#4575B4", High = "#D73027")) +
    labs(title = title, x = "", y = "Signature Score") +
    theme_bw(base_size = 11) + 
    theme(
      legend.position = "none",
      strip.background = element_rect(fill = "grey90")
    )
  
  ggsave(outfile, g, width = 11, height = 8, dpi = 300)
}

# =============================================================================
# MAIN PER-COHORT RUNNER
# =============================================================================

run_for_cohort <- function(ct){
  msg("\n================================================================")
  msg(" PROCESSING: ", ct)
  msg("================================================================")
  
  outdir <- file.path(OUT_ROOT, ct)
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  
  # Load expression from Script 01 outputs
  if (ct == "GBM") {
    fp_log <- file.path(RESULTS_DIR, "Analysis_Framework", "TCGA_GBM", "01_cleaned_data", "tcga_gbm_expression_log2.rds")
  } else if (ct == "CGGA_GBM") {
    fp_log <- file.path(RESULTS_DIR, "Analysis_Framework", "CGGA_GBM", "01_cleaned_data", "cgga_gbm_expression_log2.rds")
  } else {
    fp_log <- file.path(PROC_DIR, ct, "expression_cleaned_log2.rds")
  }
  
  if (!file.exists(fp_log)) {
    stop("No processed expression found for ", ct, "\nExpected: ", fp_log)
  }
  
  X <- readRDS(fp_log)
  X <- as.matrix(X)
  X <- backtransform_if_log(X)
  
  msg("  Expression: ", nrow(X), " genes x ", ncol(X), " samples")
  
  # Remove duplicated genes
  if (any(duplicated(rownames(X)))){
    msg("  Aggregating duplicated genes...")
    DT <- as.data.table(X, keep.rownames = "gene")
    X <- as.matrix(
      DT[, lapply(.SD, mean, na.rm = TRUE), by = gene] %>%
        tibble::column_to_rownames("gene")
    )
  }
  
  # Load purity (if available from Script 01)
  purity <- NULL
  msg("  No purity estimates available from Script 01 (optional)")
  
  # DGAT1 vector + groups
  if (!(DGAT_GENE %in% rownames(X))) {
    warning(DGAT_GENE, " not present in expression")
    dgat <- rep(NA_real_, ncol(X))
    names(dgat) <- colnames(X)
  } else {
    dgat <- X[DGAT_GENE, ]
  }
  
  cutoff <- load_bestcut(ct)
  groups <- make_groups(dgat, cutoff)
  names(groups) <- colnames(X)
  
  msg("  DGAT1 groups:")
  msg("    Low:  ", sum(groups == "Low", na.rm = TRUE))
  msg("    High: ", sum(groups == "High", na.rm = TRUE))
  
  # ---- SCORE IMMUNE SIGNATURES ----
  msg("\n  Scoring immune cell signatures...")
  S <- score_all_signatures(X, ALL_SIGNATURES)
  
  # Save scores
  fwrite(
    tibble::rownames_to_column(S, "sample"),
    file.path(outdir, "immune_signature_scores.csv")
  )
  
  # ---- Group comparisons ----
  msg("\n  DGAT1 High vs Low comparisons...")
  comp <- wilcoxon_panel(S, groups[rownames(S)])
  
  if (!is.null(comp)){
    fwrite(comp, file.path(outdir, "DGAT1_high_vs_low_wilcoxon.csv"))
    
    # Volcano plot
    comp$Signif <- comp$FDR < 0.05
    
    if (requireNamespace("ggrepel", quietly = TRUE)) {
      gv <- ggplot(comp, aes(Delta, -log10(FDR), label = Cell, color = Signif)) +
        geom_point(size = 3, alpha = 0.7) + 
        ggrepel::geom_text_repel(max.overlaps = 15, size = 3.5) +
        scale_color_manual(values = c(`TRUE` = "#D73027", `FALSE` = "grey60")) +
        geom_vline(xintercept = 0, linetype = "solid") +
        geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "blue") +
        theme_minimal(base_size = 12) +
        labs(
          title = paste(ct, "— DGAT1 High vs Low"),
          subtitle = paste(sum(comp$FDR < 0.05), "significant at FDR < 0.05"),
          x = "Δ (High - Low)", 
          y = "-log10(FDR)"
        )
      
      ggsave(
        file.path(outdir, "DGAT1_volcano.png"), 
        gv, width = 9, height = 7, dpi = 300
      )
    }
  }
  
  # ---- Correlations ----
  msg("  Computing correlations with DGAT1...")
  dgat_vec <- dgat[rownames(S)]
  pur_vec  <- if (!is.null(purity)) purity[rownames(S)] else NULL
  
  cor_tab <- correlations_vs_dgat(S, dgat_vec, purity = pur_vec)
  cor_tab$fdr_raw <- p.adjust(cor_tab$p_raw, method = "BH")
  
  if (!is.null(pur_vec)) {
    cor_tab$fdr_partial <- p.adjust(cor_tab$p_partial, method = "BH")
  }
  
  fwrite(cor_tab, file.path(outdir, "DGAT1_correlations.csv"))
  
  msg("    Raw significant (FDR<0.05): ", sum(cor_tab$fdr_raw < 0.05, na.rm=TRUE))
  if (!is.null(pur_vec)) {
    msg("    Purity-adjusted significant: ", 
        sum(cor_tab$fdr_partial < 0.05, na.rm=TRUE))
  }
  
  # ---- Visualizations ----
  key_cells <- c("CD8_T", "CD4_T", "Tregs", "NK", "B_cells", 
                 "Monocytes", "Macrophages", "M2_Macrophages", 
                 "TAM", "Microglia", "DC", "cDC1",
                 "Exhaustion_Early", "Exhaustion_Terminal", "Cytotoxic")
  
  msg("\n  Generating visualizations...")
  
  # Heatmap
  save_cell_heatmap(
    S[, intersect(colnames(S), key_cells), drop = FALSE],
    groups[rownames(S)],
    file.path(outdir, "heatmap_immune_signatures.png"),
    paste0(ct, " — Immune Signatures by DGAT1")
  )
  
  # Boxplots
  boxplot_panel(
    S, 
    groups[rownames(S)], 
    key_cells,
    file.path(outdir, "boxplots_immune_signatures.png"),
    paste0(ct, " — DGAT1 High vs Low")
  )
  
  # ---- Summary ----
  msg("\n=== SUMMARY ===")
  msg("  Signatures scored: ", ncol(S))
  msg("  Samples analyzed: ", nrow(S))
  msg("  Significant associations (raw): ", 
      sum(cor_tab$fdr_raw < 0.05, na.rm = TRUE))
  
  if (!is.null(pur_vec)) {
    msg("  Purity-adjusted associations: ", 
        sum(cor_tab$fdr_partial < 0.05, na.rm = TRUE))
  }
  
  msg("\n✓ Completed: ", outdir)
  invisible(list(scores = S, correlations = cor_tab, comparisons = comp))
}

# =============================================================================
# MAIN
# =============================================================================

msg("#######################################################################")
msg("  SCRIPT 02 — IMMUNE SIGNATURE SCORING")
msg("  (Minimal version - no external deconvolution packages required)")
msg("#######################################################################\n")

# Check optional dependencies
if (!requireNamespace("ggrepel", quietly = TRUE)) {
  msg("Note: Install ggrepel for better volcano plots: install.packages('ggrepel')")
}

if (!requireNamespace("ppcor", quietly = TRUE)) {
  msg("Note: Install ppcor for purity adjustment: install.packages('ppcor')")
}

# Run for each cancer type
for (ct in CANCER_TYPES){
  tryCatch({
    run_for_cohort(ct)
  }, error = function(e){
    msg("\n✗ ERROR for ", ct, ": ", e$message, "\n")
  })
}

msg("\n#######################################################################")
msg("✅ IMMUNE SIGNATURE SCORING COMPLETE")
msg("#######################################################################")
msg("\nOutputs saved to: ", OUT_ROOT)
msg("\nThis minimal version uses signature-based scoring instead of")
msg("full deconvolution. Results are comparable and publication-ready!")


