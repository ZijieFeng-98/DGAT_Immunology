#!/usr/bin/env Rscript
# =============================================================================
# 02_deconvolution.R — Immune deconvolution (xCell, MCPcounter, optional GBMdeconvoluteR)
# FIXED VERSION - Corrected file name mismatches from Script 01
# =============================================================================

suppressPackageStartupMessages({
  library(data.table); library(dplyr); library(tidyr); library(readr)
  library(ggplot2); library(pheatmap); library(stringr)
  library(matrixStats); library(ppcor)
})

# ====== CONFIG =================================================================
BASE_DIR    <- "/Users/zijiefeng/Desktop/Guo's lab/My_Research/DGAT_Immunology"
DATA_DIR    <- file.path(BASE_DIR, "Data")
PROC_DIR    <- file.path(BASE_DIR, "Processed_Data")  # Updated to match your structure
RESULTS_DIR <- file.path(BASE_DIR, "Results")
OUT_ROOT    <- file.path(RESULTS_DIR, "Deconvolution")
dir.create(OUT_ROOT, recursive = TRUE, showWarnings = FALSE)

CANCER_TYPES <- c("GBM")  # Add others: "BRCA", "OV", "PAAD", "CGGA_GBM"
DGAT_GENE    <- "DGAT1"

# Optional: path to GBMdeconvoluteR markers
GBMDEC_MARKER_FILE <- NULL

# =============================================================================
# UTILITIES
# =============================================================================

msg <- function(...) cat(paste0(..., "\n"))

# Best-cut loader (from survival analysis if present)
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

# IMPROVED: More robust back-transform detection
backtransform_if_log <- function(X){
  rng <- range(X, na.rm = TRUE)
  avg <- mean(X, na.rm = TRUE)
  
  # Multiple indicators for log2 scale:
  # 1. Max value < 25 (typical for log2(FPKM+1))
  # 2. Min >= 0 (log2(x+1) is always positive)
  # 3. Mean in typical log range (3-15)
  
  is_log <- (rng[2] < 25 && rng[1] >= 0 && avg > 2 && avg < 20)
  
  if (is_log) {
    msg("  Detected log2(x+1) scale, back-transforming...")
    msg("    Range before: [", round(rng[1], 2), ", ", round(rng[2], 2), "]")
    X_linear <- pmax(2^X - 1, 0)
    new_rng <- range(X_linear, na.rm = TRUE)
    msg("    Range after:  [", round(new_rng[1], 2), ", ", round(new_rng[2], 2), "]")
    return(X_linear)
  } else {
    msg("  Data appears linear (range: [", round(rng[1], 2), ", ", round(rng[2], 2), "])")
    return(X)
  }
}

# Harmonize labels across methods
standardize_cell_labels <- function(df){
  labs <- colnames(df)
  map <- c(
    "CD8_T"      = "CD8_T|CD8\\+?\\s*T.*|T.*CD8|Tcytotoxic|Tc\\b",
    "CD4_T"      = "CD4\\+?\\s*T.*|T.*CD4",
    "Tregs"      = "Tregs?|Regulatory.*T|T.*regulatory",
    "NK"         = "^NK|Natural killer",
    "B_cells"    = "^B-?cells?|B lineage|Bcell",
    "Monocytes"  = "Monocytes?",
    "Macrophages"= "Macrophages?$|^Macro|Mono.*Macro",
    "M2_Macrophages" = "M2.*Macro|Macro.*M2",
    "TAM"        = "TAMs?|Tumor.*macro",
    "DC"         = "Dendritic|^DC\\b|Myeloid.*dendritic",
    "Neutrophils"= "Neutrophils?",
    "Endothelial"= "Endothelial",
    "Fibroblasts"= "Fibroblasts?|CAFs?",
    "Microglia"  = "Microglia"
  )
  
  newlabs <- labs
  for (std in names(map)){
    pat <- map[[std]]
    hits <- grepl(pat, labs, ignore.case = TRUE)
    if (any(hits)) newlabs[hits] <- std
  }
  
  colnames(df) <- newlabs
  
  # Collapse duplicated standardized columns by mean
  df <- as.data.frame(df)
  df <- as.data.frame(sapply(split.default(df, names(df)), function(z){
    Z <- as.data.frame(z)
    if (ncol(Z) == 1) return(Z[[1]])
    rowMeans(Z, na.rm = TRUE)
  }))
  
  as.data.frame(df)
}

# Simple 0–1 row scaling
scale01_rows <- function(M){
  t(apply(M, 1, function(x){
    rng <- range(x, na.rm = TRUE)
    if (diff(rng) < 1e-9) return(rep(0.5, length(x)))
    (x - rng[1]) / (rng[2] - rng[1])
  }))
}

# Method agreement: mean Spearman across common cell types
method_agreement <- function(score_list){
  if (length(score_list) < 2) {
    return(matrix(1, 1, 1, dimnames = list(names(score_list), names(score_list))))
  }
  
  meths <- names(score_list)
  M <- matrix(NA_real_, length(meths), length(meths), 
              dimnames = list(meths, meths))
  
  for (i in seq_along(meths)){
    for (j in seq_along(meths)){
      A <- score_list[[i]]; B <- score_list[[j]]
      common_cells <- intersect(colnames(A), colnames(B))
      
      if (length(common_cells) < 2) { 
        M[i,j] <- NA
        next 
      }
      
      cors <- sapply(common_cells, function(cn) {
        suppressWarnings(
          cor(A[, cn], B[, cn], method = "spearman", use = "pairwise.complete.obs")
        )
      })
      M[i,j] <- mean(cors, na.rm = TRUE)
    }
  }
  M
}

# Wilcoxon High vs Low per cell-type
wilcoxon_panel <- function(S, groups){
  stopifnot(all(names(groups) %in% rownames(S)))
  S2 <- S[names(groups), , drop = FALSE]
  
  res <- lapply(colnames(S2), function(cell){
    x <- S2[groups == "High", cell]
    y <- S2[groups == "Low",  cell]
    
    # Need at least 3 samples per group
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

# FIXED: Correlations with proper purity adjustment
correlations_vs_dgat <- function(S, dgat_vec, purity = NULL){
  stopifnot(all(names(dgat_vec) %in% rownames(S)))
  S2 <- S[names(dgat_vec), , drop = FALSE]
  
  # Raw Spearman
  raw <- sapply(colnames(S2), function(cell){
    x <- S2[, cell]
    suppressWarnings(
      cor(x, dgat_vec, method = "spearman", use = "pairwise.complete.obs")
    )
  })
  
  # Raw p-values
  p_raw <- sapply(colnames(S2), function(cell){
    x <- S2[, cell]
    suppressWarnings(
      cor.test(x, dgat_vec, method = "spearman")$p.value
    )
  })
  
  # Purity-adjusted if available
  if (!is.null(purity)){
    adj <- sapply(colnames(S2), function(cell){
      x <- S2[, cell]
      # Remove NAs pairwise
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

# Save heatmap
save_cell_heatmap <- function(S, groups, outfile, title){
  # Scale columns (cell types) to 0-1
  S2 <- t(scale01_rows(t(S)))
  
  # Order by DGAT1 group
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

# Enhanced boxplots with statistics
boxplot_panel <- function(S, groups, keep_cells, outfile, title){
  keep_cells <- intersect(keep_cells, colnames(S))
  if (!length(keep_cells)) { 
    msg("  (boxplot) No matching cells")
    return(invisible())
  }
  
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
  
  # Calculate stats
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
    labs(title = title, x = "", y = "Score / abundance") +
    theme_bw(base_size = 11) + 
    theme(
      legend.position = "none",
      strip.background = element_rect(fill = "grey90")
    )
  
  ggsave(outfile, g, width = 11, height = 8, dpi = 300)
}

# =============================================================================
# DECONVOLUTION BACKENDS
# =============================================================================

run_xcell_iobr <- function(expr_mat){
  if (!requireNamespace("IOBR", quietly = TRUE)) {
    stop("IOBR package not installed. Install: remotes::install_github('IOBR/IOBR')")
  }
  suppressPackageStartupMessages(library(IOBR))
  
  res <- IOBR::xCell.matrix(expr = expr_mat)
  res <- as.data.frame(t(res))  # samples x cells
  
  # Ensure sample order matches
  if (!all(rownames(res) == colnames(expr_mat))) {
    msg("    Reordering xCell samples to match expression...")
    res <- res[colnames(expr_mat), , drop=FALSE]
  }
  
  res
}

run_mcp_iobr <- function(expr_mat){
  if (!requireNamespace("IOBR", quietly = TRUE)) {
    stop("IOBR package not installed.")
  }
  suppressPackageStartupMessages(library(IOBR))
  
  res <- t(IOBR::MCPcounter.estimate(expr_mat))
  res <- as.data.frame(res)
  
  if (!all(rownames(res) == colnames(expr_mat))) {
    res <- res[colnames(expr_mat), , drop=FALSE]
  }
  
  res
}

run_quantiseq_iobr <- function(expr_mat){
  if (!requireNamespace("IOBR", quietly = TRUE)) {
    stop("IOBR package not installed.")
  }
  suppressPackageStartupMessages(library(IOBR))
  
  res <- IOBR::quanTIseq_deconvolving(expr_mat, tumor = TRUE)
  
  # Returns cell x sample; transpose
  if (nrow(res) < ncol(res)) res <- t(res)
  res <- as.data.frame(res)
  
  if (!all(rownames(res) == colnames(expr_mat))) {
    res <- res[colnames(expr_mat), , drop=FALSE]
  }
  
  res
}

run_immunedeconv_safe <- function(expr_mat, method = c("xcell","mcp_counter","quantiseq")){
  if (!requireNamespace("immunedeconv", quietly = TRUE)) return(NULL)
  suppressPackageStartupMessages(library(immunedeconv))
  
  df <- as.data.frame(expr_mat) %>% tibble::rownames_to_column("gene")
  method <- match.arg(method)
  
  r <- try(
    immunedeconv::deconvolute(df, method = method, tumor = TRUE), 
    silent = TRUE
  )
  
  if (inherits(r, "try-error")) return(NULL)
  
  r <- as.data.frame(r)
  rownames(r) <- r$cell_type
  r$cell_type <- NULL
  
  t(r)  # samples x cells
}

run_gbmdeconvoluter <- function(expr_mat, marker_csv){
  if (is.null(marker_csv) || !file.exists(marker_csv)) return(NULL)
  
  M <- suppressWarnings(readr::read_csv(marker_csv, show_col_types = FALSE))
  
  gene_col <- names(M)[which(tolower(names(M)) == "gene")]
  M[[gene_col]] <- toupper(M[[gene_col]])
  
  expr <- expr_mat
  rownames(expr) <- toupper(rownames(expr))
  
  cells <- setdiff(names(M), gene_col)
  
  scores <- sapply(cells, function(cn){
    markers_df <- data.frame(
      gene = M[[gene_col]], 
      weight = M[[cn]]
    ) %>%
      filter(!is.na(weight) & weight != 0)
    
    genes_present <- intersect(markers_df$gene, rownames(expr))
    
    if (length(genes_present) == 0) {
      warning("No markers found for: ", cn)
      return(rep(NA_real_, ncol(expr)))
    }
    
    # Weighted or simple mean
    if (all(markers_df$weight == 1)) {
      colMeans(expr[genes_present, , drop=FALSE], na.rm=TRUE)
    } else {
      weights_vec <- markers_df$weight[match(genes_present, markers_df$gene)]
      apply(expr[genes_present, , drop=FALSE], 2, function(x) {
        weighted.mean(x, weights_vec, na.rm=TRUE)
      })
    }
  })
  
  as.data.frame(scores)
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
  
  # ---- FIXED: Load expression from Script 01 outputs ----
  # Script 01 creates files in Results/Analysis_Framework/<CANCER>/01_cleaned_data/
  if (ct == "GBM") {
    fp_log <- file.path(RESULTS_DIR, "Analysis_Framework", "TCGA_GBM", "01_cleaned_data", "tcga_gbm_expression_log2.rds")
  } else if (ct == "CGGA_GBM") {
    fp_log <- file.path(RESULTS_DIR, "Analysis_Framework", "CGGA_GBM", "01_cleaned_data", "cgga_gbm_expression_log2.rds")
  } else {
    # Fallback to original path
    fp_log <- file.path(PROC_DIR, ct, "expression_cleaned_log2.rds")
  }
  
  if (!file.exists(fp_log)) {
    stop("No processed expression found for ", ct, 
         "\nExpected: ", fp_log)
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
  
  # ---- FIXED: Load purity from Script 01 outputs ----
  purity <- NULL
  # Script 01 doesn't create purity files yet, so we'll skip this for now
  msg("  No purity estimates available from Script 01 (optional)")
  
  # ---- DGAT1 vector + groups ----
  if (!(DGAT_GENE %in% rownames(X))) {
    warning(DGAT_GENE, " not present in expression for ", ct)
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
  
  # ---- Run deconvolution methods ----
  scores <- list()
  
  # xCell via IOBR
  if (requireNamespace("IOBR", quietly = TRUE)) {
    msg("\n  Running xCell via IOBR...")
    xres <- try(run_xcell_iobr(X), silent = TRUE)
    if (!inherits(xres, "try-error")) {
      scores[["xCell"]] <- xres
      msg("    ✓ xCell: ", ncol(xres), " cell types")
    } else {
      msg("    ✗ xCell failed: ", as.character(xres))
    }
  }
  
  # MCPcounter via IOBR
  if (requireNamespace("IOBR", quietly = TRUE)) {
    msg("  Running MCPcounter via IOBR...")
    mres <- try(run_mcp_iobr(X), silent = TRUE)
    if (!inherits(mres, "try-error")) {
      scores[["MCPcounter"]] <- mres
      msg("    ✓ MCPcounter: ", ncol(mres), " cell types")
    } else {
      msg("    ✗ MCPcounter failed")
    }
  }
  
  # quanTIseq via IOBR (optional)
  if (requireNamespace("IOBR", quietly = TRUE)) {
    msg("  Running quanTIseq via IOBR (optional)...")
    qres <- try(run_quantiseq_iobr(X), silent = TRUE)
    if (!inherits(qres, "try-error")) {
      scores[["quanTIseq"]] <- qres
      msg("    ✓ quanTIseq: ", ncol(qres), " cell types")
    } else {
      msg("    ✗ quanTIseq failed (optional)")
    }
  }
  
  # immunedeconv fallback
  if (length(scores) == 0 && requireNamespace("immunedeconv", quietly = TRUE)) {
    msg("  Running xCell via immunedeconv (fallback)...")
    x2 <- run_immunedeconv_safe(X, "xcell")
    if (!is.null(x2)) {
      scores[["xCell"]] <- x2
      msg("    ✓ xCell (immunedeconv): ", ncol(x2), " cell types")
    }
  }
  
  # GBMdeconvoluteR (optional)
  if (!is.null(GBMDEC_MARKER_FILE)) {
    msg("  Running GBMdeconvoluteR (marker scoring)...")
    gsc <- try(run_gbmdeconvoluter(X, GBMDEC_MARKER_FILE), silent = TRUE)
    if (!inherits(gsc, "try-error") && !is.null(gsc)) {
      scores[["GBMdeconv"]] <- gsc
      msg("    ✓ GBMdeconv: ", ncol(gsc), " cell types")
    }
  }
  
  if (length(scores) == 0) {
    stop("No deconvolution method succeeded. Install IOBR: remotes::install_github('IOBR/IOBR')")
  }
  
  # ---- Save raw + standardized tables ----
  msg("\n  Saving deconvolution results...")
  std_scores <- list()
  
  for (nm in names(scores)){
    tab <- scores[[nm]]
    
    # Align samples
    tab <- tab[intersect(rownames(tab), colnames(X)), , drop = FALSE]
    
    # Save raw
    fwrite(
      tibble::rownames_to_column(as.data.frame(tab), "sample"),
      file.path(outdir, paste0("deconv_", nm, "_raw.csv"))
    )
    
    # Standardize and save
    std <- standardize_cell_labels(tab)
    fwrite(
      tibble::rownames_to_column(as.data.frame(std), "sample"),
      file.path(outdir, paste0("deconv_", nm, "_standardized.csv"))
    )
    
    std_scores[[nm]] <- std
  }
  
  # ---- Method agreement ----
  if (length(std_scores) >= 2){
    msg("  Computing method agreement...")
    M <- method_agreement(std_scores)
    fwrite(
      as.data.frame(M), 
      file.path(outdir, "method_agreement_mean_spearman.csv")
    )
    
    pheatmap(
      M, 
      display_numbers = TRUE, 
      main = paste(ct, "Method Agreement (mean ρ)"),
      filename = file.path(outdir, "method_agreement_heatmap.png"),
      width = 6, 
      height = 5
    )
  }
  
  # ---- Pick primary table (xCell preferred) ----
  pick <- if ("xCell" %in% names(std_scores)) "xCell" else names(std_scores)[1]
  msg("  Using ", pick, " as primary method for visualizations")
  
  S <- std_scores[[pick]]
  S <- S[intersect(rownames(S), colnames(X)), , drop = FALSE]
  S <- S[colnames(X), , drop = FALSE]
  
  # ---- Group comparisons ----
  msg("\n  DGAT1 High vs Low comparisons...")
  comp <- wilcoxon_panel(S, groups[colnames(X)])
  
  if (!is.null(comp)){
    fwrite(
      comp, 
      file.path(outdir, paste0("DGAT1_high_vs_low_", pick, "_wilcoxon.csv"))
    )
    
    # Volcano plot
    comp$Signif <- comp$FDR < 0.05
    
    gv <- ggplot(comp, aes(Delta, -log10(FDR), label = Cell, color = Signif)) +
      geom_point(size = 3, alpha = 0.7) + 
      ggrepel::geom_text_repel(max.overlaps = 15, size = 3.5) +
      scale_color_manual(values = c(`TRUE` = "#D73027", `FALSE` = "grey60")) +
      geom_vline(xintercept = c(-0.1, 0.1), linetype = "dashed", color = "grey40") +
      geom_vline(xintercept = 0, linetype = "solid", color = "black") +
      geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "blue") +
      annotate("text", x = Inf, y = -log10(0.05), label = "FDR = 0.05", 
               hjust = 1.1, vjust = -0.5, size = 3, color = "blue") +
      theme_minimal(base_size = 12) +
      labs(
        title = paste(ct, pick, "— DGAT1 High vs Low"),
        subtitle = paste(sum(comp$FDR < 0.05), "significant at FDR < 0.05"),
        x = "Δ (High - Low)", 
        y = "-log10(FDR)"
      )
    
    ggsave(
      file.path(outdir, paste0("DGAT1_volcano_", pick, ".png")), 
      gv, width = 9, height = 7, dpi = 300
    )
  }
  
  # ---- Correlations with DGAT1 ----
  msg("  Computing correlations with DGAT1...")
  dgat_vec <- dgat[colnames(S)]
  pur_vec  <- if (!is.null(purity)) purity[colnames(S)] else NULL
  
  cor_tab <- correlations_vs_dgat(S, dgat_vec, purity = pur_vec)
  
  # Add FDR
  cor_tab$fdr_raw <- p.adjust(cor_tab$p_raw, method = "BH")
  
  if (!is.null(pur_vec)) {
    cor_tab$fdr_partial <- p.adjust(cor_tab$p_partial, method = "BH")
  } else {
    cor_tab$fdr_partial <- NA_real_
  }
  
  fwrite(
    cor_tab, 
    file.path(outdir, paste0("DGAT1_correlations_", pick, ".csv"))
  )
  
  msg("    Raw significant (FDR<0.05): ", sum(cor_tab$fdr_raw < 0.05, na.rm=TRUE))
  if (!is.null(pur_vec)) {
    msg("    Purity-adjusted significant: ", sum(cor_tab$fdr_partial < 0.05, na.rm=TRUE))
  }
  
  # ---- Visualizations ----
  key_cells <- c("CD8_T", "Tregs", "NK", "Monocytes", "Macrophages", 
                 "M2_Macrophages", "TAM", "DC", "B_cells", "Microglia")
  
  # Heatmap
  msg("\n  Generating heatmap...")
  save_cell_heatmap(
    S[, intersect(colnames(S), key_cells), drop = FALSE],
    groups[colnames(S)],
    file.path(outdir, paste0("heatmap_keycells_", pick, ".png")),
    paste0(ct, " ", pick, " — Key Immune Populations by DGAT1")
  )
  
  # Boxplots
  msg("  Generating boxplots...")
  boxplot_panel(
    S, 
    groups[colnames(S)], 
    key_cells,
    file.path(outdir, paste0("boxplots_keycells_", pick, ".png")),
    paste0(ct, " ", pick, " — DGAT1 High vs Low")
  )
  
  # ---- Validation Summary ----
  msg("\n=== VALIDATION ===")
  msg("  Methods run: ", paste(names(scores), collapse = ", "))
  msg("  Samples in analysis: ", nrow(S))
  msg("  Cell types detected: ", ncol(S))
  msg("  Significant DGAT1 associations (raw): ", 
      sum(cor_tab$fdr_raw < 0.05, na.rm = TRUE))
  
  if (!is.null(pur_vec)) {
    msg("  Purity-adjusted associations: ", 
        sum(cor_tab$fdr_partial < 0.05, na.rm = TRUE))
  }
  
  msg("\n✓ Completed: ", outdir)
  invisible(list(scores = std_scores, correlations = cor_tab, comparisons = comp))
}

# =============================================================================
# MAIN
# =============================================================================

msg("#######################################################################")
msg("  SCRIPT 02 — IMMUNE DECONVOLUTION")
msg("#######################################################################\n")

# Check dependencies
if (!requireNamespace("ggrepel", quietly = TRUE)) {
  stop("Please install ggrepel: install.packages('ggrepel')")
}

if (!requireNamespace("ppcor", quietly = TRUE)) {
  stop("Please install ppcor: install.packages('ppcor')")
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
msg("✅ DECONVOLUTION COMPLETE")
msg("#######################################################################")
msg("\nOutputs saved to: ", OUT_ROOT)
msg("\nNext step: Run Script 03 for GSVA pathway analysis")
