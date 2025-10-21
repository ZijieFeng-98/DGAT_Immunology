#!/usr/bin/env Rscript
# =============================================================================
# 02_deconvolution.R — Immune deconvolution (xCell, MCPcounter, optional GBMdeconvoluteR)
# Inputs: Data/Processed/<CANCER>/expression_cleaned.rds (or *_log2.rds), clinical_standardized.csv
# Outputs: Results/Deconvolution/<CANCER>/*
# - Per-method cell scores
# - Method agreement matrix
# - DGAT1 High/Low comparisons (Wilcoxon + FDR)
# - Raw + purity-adjusted (partial) correlations vs DGAT1
# - Heatmaps/boxplots (publication-ready defaults)
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

CANCER_TYPES <- c("GBM", "CGGA_GBM")       # You can add BRCA/OV/PAAD later if desired
DGAT_GENE    <- "DGAT1"

# Optional: path to GBMdeconvoluteR markers (CSV with genes x cell-types, values 0/1 or weights)
GBMDEC_MARKER_FILE <- NULL  # e.g., "Resources/GBMdeconvoluter_markers.csv" if available

# =============================================================================
# UTILITIES
# =============================================================================

# Log helpers
msg <- function(...) cat(paste0(..., "\n"))

# Best-cut loader (from Step 1 survival output if present)
load_bestcut <- function(cohort_dir){
  fp <- file.path(RESULTS_DIR, "Bulk", "Survival", cohort_dir, "DGAT1_bestcut_info.csv")
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

# Back-transform from log2(x+1) if needed
backtransform_if_log <- function(X){
  # Heuristic: if many values < 25 and some small decimals, likely log2+1
  rng <- range(X, na.rm = TRUE)
  if (rng[2] < 25 && rng[1] >= 0) {
    return(pmax(2^X - 1, 0))
  }
  X
}

# Harmonize labels across methods
standardize_cell_labels <- function(df){
  # df: samples x cell-types; colnames are raw labels
  labs <- colnames(df)
  map <- c(
    "CD8_T"      = "CD8 T cells|CD8\\+?\\s*T.*|Tcytotoxic|Tc",
    "CD4_T"      = "CD4\\+?\\s*T.*",
    "Tregs"      = "Tregs?|Regulatory.*",
    "NK"         = "^NK|Natural killer",
    "B_cells"    = "^B-?cells?|B lineage|Bcell",
    "Monocytes"  = "Monocytes?",
    "Macrophages"= "Macrophages?$|TAMs?|Mono/Macro",
    "DC"         = "Dendritic|Myeloid dendritic",
    "Neutrophils"= "Neutrophils?",
    "Endothelial"= "Endothelial",
    "Fibroblasts"= "Fibroblasts?|CAFs?",
    "Microglia"  = "Microglia",
    "TAM"        = "TAMs?"
  )
  newlabs <- labs
  for (std in names(map)){
    pat <- map[[std]]
    hits <- grepl(pat, labs, ignore.case = TRUE)
    newlabs[hits] <- std
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

# Simple 0–1 row scaling for heatmaps
scale01_rows <- function(M){
  t(apply(M, 1, function(x){
    rng <- range(x, na.rm = TRUE)
    if (diff(rng) < 1e-9) return(rep(0.5, length(x)))
    (x - rng[1]) / (rng[2] - rng[1])
  }))
}

# Method agreement: mean Spearman rho across common cell types
method_agreement <- function(score_list){
  if (length(score_list) < 2) return(matrix(1, 1, 1, dimnames = list(names(score_list), names(score_list))))
  meths <- names(score_list)
  M <- matrix(NA_real_, length(meths), length(meths), dimnames = list(meths, meths))
  for (i in seq_along(meths)){
    for (j in seq_along(meths)){
      A <- score_list[[i]]; B <- score_list[[j]]
      common_cells <- intersect(colnames(A), colnames(B))
      if (length(common_cells) < 2) { M[i,j] <- NA; next }
      cors <- sapply(common_cells, function(cn)
        suppressWarnings(cor(A[, cn], B[, cn], method = "spearman", use = "pairwise.complete.obs")))
      M[i,j] <- mean(cors, na.rm = TRUE)
    }
  }
  M
}

# Wilcoxon compare High vs Low per cell-type
wilcoxon_panel <- function(S, groups){
  stopifnot(all(names(groups) %in% rownames(S)))
  S2 <- S[names(groups), , drop = FALSE]
  res <- lapply(colnames(S2), function(cell){
    x <- S2[groups == "High", cell]
    y <- S2[groups == "Low",  cell]
    if (length(na.omit(x)) < 3 || length(na.omit(y)) < 3) return(NULL)
    wt <- suppressWarnings(wilcox.test(x, y))
    delta <- mean(x, na.rm = TRUE) - mean(y, na.rm = TRUE)
    data.frame(Cell = cell, Mean_High = mean(x, na.rm = TRUE),
               Mean_Low = mean(y, na.rm = TRUE), Delta = delta, p_value = wt$p.value)
  })
  T <- do.call(rbind, res[!sapply(res, is.null)])
  if (is.null(T)) return(NULL)
  T$FDR <- p.adjust(T$p_value, method = "BH")
  T[order(T$FDR, -abs(T$Delta)), ]
}

# Correlations vs DGAT1 (raw and purity-adjusted)
correlations_vs_dgat <- function(S, dgat_vec, purity = NULL){
  # S: samples x cell-types
  stopifnot(all(names(dgat_vec) %in% rownames(S)))
  S2 <- S[names(dgat_vec), , drop = FALSE]
  raw <- sapply(colnames(S2), function(cell){
    x <- S2[, cell]; suppressWarnings(cor(x, dgat_vec, method = "spearman", use = "pairwise.complete.obs"))
  })
  if (!is.null(purity)){
    adj <- sapply(colnames(S2), function(cell){
      x <- S2[, cell]
      pd <- try(ppcor::pcor.test(x, dgat_vec, purity, method = "spearman"), silent = TRUE)
      if (inherits(pd, "try-error")) return(NA_real_)
      pd$estimate
    })
  } else {
    adj <- rep(NA_real_, length(raw)); names(adj) <- names(raw)
  }
  tibble(Cell = names(raw), rho_raw = as.numeric(raw), rho_partial = as.numeric(adj))
}

# Save heatmap of selected cells
save_cell_heatmap <- function(S, groups, outfile, title){
  S2 <- t(scale01_rows(t(S)))  # scale per cell-type (columns)
  ord <- order(groups)
  ann <- data.frame(DGAT1 = groups[ord]); rownames(ann) <- rownames(S2)[ord]
  ann_colors <- list(DGAT1 = c(Low = "#4575B4", High = "#D73027"))
  pheatmap(t(S2)[, ord, drop = FALSE],
           annotation_col = ann, annotation_colors = ann_colors,
           show_colnames = FALSE, cluster_cols = FALSE, main = title,
           filename = outfile, width = 10, height = 6)
}

# Boxplots for key panels
boxplot_panel <- function(S, groups, keep_cells, outfile, title){
  keep_cells <- intersect(keep_cells, colnames(S))
  if (!length(keep_cells)) { msg("  (boxplot) No matching cells"); return(invisible())}
  df <- as.data.frame(S); df$sample <- rownames(df); df$DGAT1 <- groups[df$sample]
  long <- pivot_longer(df, cols = all_of(keep_cells), names_to = "Cell", values_to = "Score") %>%
    filter(!is.na(DGAT1))
  stat <- long %>%
    group_by(Cell) %>%
    summarise(p_value = suppressWarnings(wilcox.test(Score ~ DGAT1)$p.value),
              .groups = "drop") %>%
    mutate(label = case_when(p_value < 0.001 ~ "***", p_value < 0.01 ~ "**",
                             p_value < 0.05 ~ "*", TRUE ~ ""))

  g <- ggplot(long, aes(DGAT1, Score, fill = DGAT1)) +
    geom_boxplot(outlier.size = 0.6, width = 0.7, alpha = 0.85) +
    facet_wrap(~Cell, scales = "free_y") +
    scale_fill_manual(values = c(Low = "#4575B4", High = "#D73027")) +
    labs(title = title, x = "", y = "Score / abundance") +
    theme_bw(base_size = 11) + theme(legend.position = "none")
  ggsave(outfile, g, width = 11, height = 7, dpi = 300)
}

# =============================================================================
# BACKENDS
# =============================================================================
run_xcell_iobr <- function(expr_mat){
  if (!requireNamespace("IOBR", quietly = TRUE)) stop("IOBR package not installed.")
  suppressPackageStartupMessages(library(IOBR))
  # IOBR::xCell.matrix expects genes as rownames; TPM-like scale is OK
  res <- IOBR::xCell.matrix(expr = expr_mat)
  # xCell returns cell x sample; make samples x cells
  as.data.frame(t(res))
}

run_mcp_iobr <- function(expr_mat){
  if (!requireNamespace("IOBR", quietly = TRUE)) stop("IOBR package not installed.")
  suppressPackageStartupMessages(library(IOBR))
  res <- t(IOBR::MCPcounter.estimate(expr_mat))
  as.data.frame(res)
}

run_quantiseq_iobr <- function(expr_mat){
  if (!requireNamespace("IOBR", quietly = TRUE)) stop("IOBR package not installed.")
  suppressPackageStartupMessages(library(IOBR))
  res <- IOBR::quanTIseq_deconvolving(expr_mat, tumor = TRUE)
  # returns cell x sample; transpose to samples x cells if needed
  if (nrow(res) < ncol(res)) res <- t(res)
  as.data.frame(res)
}

run_immunedeconv_safe <- function(expr_mat, method = c("xcell","mcp_counter","quantiseq")){
  if (!requireNamespace("immunedeconv", quietly = TRUE)) return(NULL)
  suppressPackageStartupMessages(library(immunedeconv))
  df <- as.data.frame(expr_mat) %>% tibble::rownames_to_column("gene")
  method <- match.arg(method)
  r <- try(immunedeconv::deconvolute(df, method = method, tumor = TRUE), silent = TRUE)
  if (inherits(r, "try-error")) return(NULL)
  r <- as.data.frame(r); rownames(r) <- r$cell_type; r$cell_type <- NULL
  t(r)  # samples x cells
}

run_gbmdeconvoluter <- function(expr_mat, marker_csv){
  # Simple marker-based scoring: mean of marker genes per cell-type
  if (is.null(marker_csv) || !file.exists(marker_csv)) return(NULL)
  M <- suppressWarnings(readr::read_csv(marker_csv, show_col_types = FALSE))
  # Expect columns: gene, <cell1>, <cell2>, ... with 0/1 or weights
  stopifnot("gene" %in% tolower(names(M)))
  gene_col <- names(M)[which(tolower(names(M)) == "gene")]
  M[[gene_col]] <- toupper(M[[gene_col]])
  expr <- expr_mat
  rownames(expr) <- toupper(rownames(expr))
  cells <- setdiff(names(M), gene_col)
  scores <- sapply(cells, function(cn){
    weights <- M[[cn]]
    g <- M[[gene_col]][!is.na(weights) & weights != 0]
    gg <- intersect(g, rownames(expr))
    if (!length(gg)) return(rep(NA_real_, ncol(expr)))
    colMeans(expr[gg, , drop = FALSE], na.rm = TRUE)
  })
  as.data.frame(t(scores)) %>% t() %>% as.data.frame() # samples x cells
}

# =============================================================================
# MAIN PER-COHORT RUNNER
# =============================================================================
run_for_cohort <- function(ct){
  msg("\n=== ", ct, " ===")
  outdir <- file.path(OUT_ROOT, ct)
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

  # ---- Load cleaned expression (from Script 01) ----
  # Updated paths to match your existing structure
  if (ct == "GBM") {
    fp1 <- file.path(PROC_DIR, "TCGA_GBM_Clean", "TCGA_GBM_Expression_Cleaned.rds")
  } else if (ct == "CGGA_GBM") {
    fp1 <- file.path(PROC_DIR, "Clean_CGGA_Cohort", "mRNAseq_693_GBM", "CGGA_mRNAseq_693_GBM_expr_clean.rds")
  } else {
    # Fallback to original paths
    fp1 <- file.path(PROC_DIR, ct, "expression_cleaned.rds")
    fp2 <- file.path(PROC_DIR, ct, "expression_cleaned_log2.rds")
  }
  
  if (file.exists(fp1)) {
    X <- readRDS(fp1)
    X <- as.matrix(X)
    # if linear scale already, keep; many methods prefer non-log TPM/FPKM
  } else if (exists("fp2") && file.exists(fp2)) {
    X <- readRDS(fp2)
    X <- as.matrix(X)
    X <- backtransform_if_log(X)
  } else {
    stop("No processed expression found for ", ct, " (expected *_cleaned.rds or *_log2.rds)")
  }
  msg("  Expression: ", nrow(X), " genes x ", ncol(X), " samples")

  # Remove duplicated genes by mean
  if (any(duplicated(rownames(X)))){
    msg("  Aggregating duplicated genes...")
    DT <- as.data.table(X, keep.rownames = "gene")
    X <- as.matrix(DT[, lapply(.SD, mean, na.rm = TRUE), by = gene] %>%
                     tibble::column_to_rownames("gene"))
  }

  # ---- Load clinical + purity (from Script 01) ----
  # Updated paths to match your existing structure
  if (ct == "GBM") {
    clin_fp <- file.path(PROC_DIR, "TCGA_GBM_Clean", "TCGA_GBM_Clinical_Cleaned.csv")
    est_fp <- file.path(RESULTS_DIR, "Analysis_Framework", "TCGA_GBM", "01_cleaned_data", "tcga_gbm_purity_estimate.csv")
  } else if (ct == "CGGA_GBM") {
    clin_fp <- file.path(PROC_DIR, "Clean_CGGA_Cohort", "mRNAseq_693_GBM", "CGGA_mRNAseq_693_GBM_clin_clean.csv")
    est_fp <- file.path(RESULTS_DIR, "Analysis_Framework", "CGGA_GBM", "01_cleaned_data", "cgga_gbm_purity_estimate.csv")
  } else {
    clin_fp <- file.path(PROC_DIR, ct, "clinical_standardized.csv")
    est_fp <- file.path(PROC_DIR, ct, "estimate_purity.csv")
  }
  
  clin <- if (file.exists(clin_fp)) readr::read_csv(clin_fp, show_col_types = FALSE) else NULL

  purity <- NULL
  if (file.exists(est_fp)) {
    est <- suppressWarnings(readr::read_csv(est_fp, show_col_types = FALSE))
    if (all(c("sample","purity") %in% names(est))) {
      purity <- est$purity; names(purity) <- est$sample
    }
  }

  # ---- DGAT1 vector + groups ----
  if (!(DGAT_GENE %in% rownames(X))) {
    warning("DGAT1 not present in expression for ", ct)
    dgat <- rep(NA_real_, ncol(X)); names(dgat) <- colnames(X)
  } else {
    dgat <- X[DGAT_GENE, ]
  }

  cutoff <- load_bestcut(ct)
  groups <- make_groups(dgat, cutoff)
  names(groups) <- colnames(X)
  msg("  Groups: ", sum(groups == "Low", na.rm = TRUE), " Low; ",
      sum(groups == "High", na.rm = TRUE), " High")

  # ---- Run deconvolution methods ----
  scores <- list()

  # Prefer IOBR backends (robust, no license) — xCell
  xcell_ok <- FALSE
  if (requireNamespace("IOBR", quietly = TRUE)) {
    msg("  Running xCell via IOBR ...")
    xres <- try(run_xcell_iobr(X), silent = TRUE)
    if (!inherits(xres, "try-error")) { scores[["xCell"]] <- xres; xcell_ok <- TRUE; msg("    ✓ xCell ok") }
  }
  # MCPcounter via IOBR
  if (requireNamespace("IOBR", quietly = TRUE)) {
    msg("  Running MCPcounter via IOBR ...")
    mres <- try(run_mcp_iobr(X), silent = TRUE)
    if (!inherits(mres, "try-error")) { scores[["MCPcounter"]] <- mres; msg("    ✓ MCPcounter ok") }
  }
  # quanTIseq via IOBR (optional; slower)
  if (requireNamespace("IOBR", quietly = TRUE)) {
    msg("  Running quanTIseq via IOBR (optional) ...")
    qres <- try(run_quantiseq_iobr(X), silent = TRUE)
    if (!inherits(qres, "try-error")) { scores[["quanTIseq"]] <- qres; msg("    ✓ quanTIseq ok") }
  }

  # immunedeconv fallback if IOBR not available
  if (!xcell_ok && requireNamespace("immunedeconv", quietly = TRUE)) {
    msg("  Running xCell via immunedeconv (fallback) ...")
    x2 <- run_immunedeconv_safe(X, "xcell")
    if (!is.null(x2)) { scores[["xCell"]] <- x2; msg("    ✓ xCell ok (immunedeconv)") }
  }

  # Optional GBMdeconvoluteR (marker mean scoring)
  if (!is.null(GBMDEC_MARKER_FILE)) {
    msg("  Running GBMdeconvoluteR (marker scoring) ...")
    gsc <- try(run_gbmdeconvoluter(X, GBMDEC_MARKER_FILE), silent = TRUE)
    if (!inherits(gsc, "try-error") && !is.null(gsc)) {
      scores[["GBMdeconv"]] <- gsc; msg("    ✓ GBMdeconvoluter ok")
    }
  }

  if (length(scores) == 0) stop("No deconvolution method succeeded. Install IOBR or immunedeconv.")

  # ---- Save raw per-method tables (standardized labels too) ----
  std_scores <- list()
  for (nm in names(scores)){
    tab <- scores[[nm]]
    # Ensure samples align with expression
    tab <- tab[intersect(rownames(tab), colnames(X)), , drop = FALSE]
    fwrite(tibble::rownames_to_column(as.data.frame(tab), "sample"),
           file.path(outdir, paste0("deconv_", nm, "_raw.csv")))
    # Standardize labels & save
    std <- standardize_cell_labels(tab)
    fwrite(tibble::rownames_to_column(as.data.frame(std), "sample"),
           file.path(outdir, paste0("deconv_", nm, "_standardized.csv")))
    std_scores[[nm]] <- std
  }

  # ---- Method agreement ----
  if (length(std_scores) >= 2){
    M <- method_agreement(std_scores)
    fwrite(as.data.frame(M), file.path(outdir, "method_agreement_mean_spearman.csv"))
    # quickheat
    pheatmap(M, display_numbers = TRUE, main = paste(ct, "Method agreement (ρ)"),
             filename = file.path(outdir, "method_agreement_heatmap.png"),
             width = 6, height = 5)
  }

  # ---- Pick a primary table for visualization (xCell preferred) ----
  pick <- if ("xCell" %in% names(std_scores)) "xCell" else names(std_scores)[1]
  S <- std_scores[[pick]]
  # Align to expression sample order
  S <- S[intersect(rownames(S), colnames(X)), , drop = FALSE]
  S <- S[colnames(X), , drop = FALSE]

  # ---- Group comparisons ----
  msg("  Group comparisons on ", pick, " scores ...")
  comp <- wilcoxon_panel(S, groups[colnames(X)])
  if (!is.null(comp)){
    fwrite(comp, file.path(outdir, paste0("DGAT1_high_vs_low_", pick, "_wilcoxon.csv")))
    # Volcano-like plot
    comp$Signif <- comp$FDR < 0.05
    gv <- ggplot(comp, aes(Delta, -log10(FDR), label = Cell, color = Signif)) +
      geom_point() + ggrepel::geom_text_repel(max.overlaps = 12, size = 3) +
      scale_color_manual(values = c(`TRUE` = "#D73027", `FALSE` = "grey60")) +
      geom_vline(xintercept = 0, linetype = "dashed") +
      geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
      theme_minimal(base_size = 12) +
      labs(title = paste(ct, pick, "DGAT1 High vs Low"), x = "Δ (High - Low)", y = "-log10(FDR)")
    ggsave(file.path(outdir, paste0("DGAT1_volcano_", pick, ".png")), gv, width = 8, height = 6, dpi = 300)
  }

  # ---- Correlations with DGAT1 (raw and partial by purity) ----
  msg("  Correlations vs DGAT1 (raw + purity-adjusted) ...")
  dgat_vec <- dgat[colnames(S)]
  pur_vec  <- if (!is.null(purity)) purity[colnames(S)] else NULL
  cor_tab <- correlations_vs_dgat(S, dgat_vec, purity = pur_vec)
  # BH over multiple cells
  cor_tab$p_raw <- sapply(colnames(S), function(cell){
    suppressWarnings(cor.test(S[,cell], dgat_vec, method = "spearman")$p.value)
  })
  cor_tab$fdr_raw <- p.adjust(cor_tab$p_raw, method = "BH")
  if (!is.null(pur_vec)){
    cor_tab$p_partial <- sapply(colnames(S), function(cell){
      pd <- try(ppcor::pcor.test(S[,cell], dgat_vec, pur_vec, method = "spearman"), silent = TRUE)
      if (inherits(pd, "try-error")) return(NA_real_)
      pd$p.value
    })
    cor_tab$fdr_partial <- p.adjust(cor_tab$p_partial, method = "BH")
  } else {
    cor_tab$p_partial <- NA_real_; cor_tab$fdr_partial <- NA_real_
  }
  fwrite(cor_tab, file.path(outdir, paste0("DGAT1_correlations_", pick, ".csv")))

  # ---- Visuals: heatmap + boxplots for key cells ----
  key_cells <- c("CD8_T","Tregs","NK","Monocytes","Macrophages","DC","B_cells","Microglia","TAM")
  # Heatmap
  save_cell_heatmap(t(S[, intersect(colnames(S), key_cells), drop = FALSE]),
                    groups[colnames(S)],
                    file.path(outdir, paste0("heatmap_keycells_", pick, ".png")),
                    paste0(ct, " ", pick, " — key populations by DGAT1 group"))
  # Boxplots
  boxplot_panel(S, groups[colnames(S)], key_cells,
                file.path(outdir, paste0("boxplots_keycells_", pick, ".png")),
                paste0(ct, " ", pick, " — key populations vs DGAT1 group"))

  msg("  ✓ Completed: ", outdir)
}
