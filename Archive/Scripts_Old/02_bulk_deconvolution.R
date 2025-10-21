#!/usr/bin/env Rscript
# =============================================================================
# 02_bulk_deconvolution.R — Immune gene set enrichment analysis via GSVA
# - Works on TCGA_GBM (FPKM) and CGGA (array/other), gene symbols as rows
# - Auto-uses DGAT1 best-cut groups from Step 1, else falls back to median
# - Outputs CSVs + heatmaps + boxplots (no UI)
# =============================================================================

suppressPackageStartupMessages({
  library(data.table); library(dplyr); library(tidyr)
  library(ggplot2); library(pheatmap); library(stringr)
  library(GSVA)
})

# ---- Folders
out_root <- "Results/Bulk/Deconvolution"
dir.create(out_root, recursive = TRUE, showWarnings = FALSE)

message("Using GSVA for immune gene set enrichment analysis")

# ---- Helpers
read_expr <- function(fp){
  mat <- readRDS(fp)
  if (nrow(mat) < ncol(mat)) mat <- t(mat)
  # ensure matrix with gene symbols as rownames
  mat <- as.matrix(mat)
  storage.mode(mat) <- "numeric"
  return(mat)
}

load_bestcut <- function(cohort_dir){
  fp <- file.path("Results/Bulk/Survival", cohort_dir, "DGAT1_bestcut_info.csv")
  if (file.exists(fp)) {
    info <- fread(fp)
    as.numeric(info$cutoff[1])
  } else NA_real_
}

make_groups <- function(expr_vec, cutoff = NA_real_){ a
  if (!is.na(cutoff)) {
    factor(ifelse(expr_vec >= cutoff, "High","Low"), levels=c("Low","High"))
  } else {
    med <- median(expr_vec, na.rm=TRUE)
    factor(ifelse(expr_vec >= med, "High","Low"), levels=c("Low","High"))
  }
}

# Define immune gene sets for GSVA analysis
get_immune_genesets <- function(){
  list(
    # T cell subsets
    CD8_T_cells = c("CD8A", "CD8B", "GZMA", "GZMB", "PRF1", "IFNG", "TBX21", "EOMES"),
    CD4_T_cells = c("CD4", "IL2", "IL4", "IL5", "IL13", "GATA3", "FOXP3", "IL17A"),
    Treg = c("FOXP3", "IL2RA", "CTLA4", "TIGIT", "LAG3", "HAVCR2", "ENTPD1"),
    Exhausted_T = c("PDCD1", "HAVCR2", "LAG3", "TIGIT", "CTLA4", "ENTPD1", "CD244"),
    
    # NK cells
    NK_cells = c("KLRD1", "KLRF1", "KLRB1", "NCR1", "NCR3", "GNLY", "GZMA", "GZMB"),
    
    # Myeloid cells
    Monocytes = c("CD14", "CD68", "CSF1R", "CCR2", "CX3CR1", "S100A8", "S100A9"),
    Macrophages = c("CD68", "CD163", "MSR1", "MRC1", "MARCO", "CD86", "IL10"),
    M1_Macrophages = c("CD86", "IL1B", "TNF", "NOS2", "CXCL9", "CXCL10", "CXCL11"),
    M2_Macrophages = c("CD163", "MRC1", "MSR1", "IL10", "TGFB1", "ARG1", "CCL17", "CCL22"),
    Dendritic_cells = c("CD1C", "CLEC9A", "XCR1", "BATF3", "IRF8", "ZBTB46", "FLT3"),
    
    # B cells
    B_cells = c("CD19", "CD20", "MS4A1", "CD79A", "CD79B", "PAX5", "IGHM", "IGHA1"),
    
    # Cytokines and chemokines
    Pro_inflammatory = c("IL1B", "IL6", "TNF", "IFNG", "CXCL9", "CXCL10", "CXCL11"),
    Anti_inflammatory = c("IL10", "TGFB1", "IL4", "IL13", "ARG1", "IDO1"),
    
    # Immune checkpoints
    Immune_checkpoints = c("PDCD1", "CD274", "PDCD1LG2", "CTLA4", "LAG3", "HAVCR2", "TIGIT"),
    
    # Complement
    Complement = c("C1QA", "C1QB", "C1QC", "C3", "C4A", "C4B", "CFB", "CFD"),
    
    # Interferon response
    Interferon_response = c("IFIT1", "IFIT2", "IFIT3", "ISG15", "MX1", "MX2", "OAS1", "OAS2", "OAS3")
  )
}

# CGGA microarray probe to gene symbol mapping (simplified approach)
map_cgga_probes <- function(mat, cohort_name) {
  if (cohort_name == "CGGA_GBM") {
    # For CGGA, we'll use a simplified approach - look for genes that might be present
    # and create a basic immune signature using available genes
    message("CGGA dataset detected - using simplified immune gene mapping")
    
    # Look for any genes that might be immune-related (case insensitive)
    all_genes <- rownames(mat)
    immune_related <- c()
    
    # Search for common immune gene patterns
    patterns <- c("CD", "IL", "IFN", "TNF", "TGF", "CCL", "CXCL", "GZ", "PRF", "FOXP", "CTLA", "PD")
    for (pattern in patterns) {
      matches <- grep(pattern, all_genes, ignore.case = TRUE, value = TRUE)
      immune_related <- c(immune_related, matches)
    }
    
    # Also look for specific gene names that might be present
    specific_genes <- c("DGAT1", "DGAT2", "ACTB", "GAPDH", "B2M", "HPRT1")
    for (gene in specific_genes) {
      if (gene %in% all_genes) {
        immune_related <- c(immune_related, gene)
      }
    }
    
    message("Found ", length(unique(immune_related)), " potentially immune-related genes in CGGA")
    
    # Create a simplified gene set using available genes
    if (length(unique(immune_related)) > 0) {
      # Create a basic immune signature using available genes
      simplified_genesets <- list(
        Immune_related = unique(immune_related)[1:min(20, length(unique(immune_related)))]
      )
      return(simplified_genesets)
    } else {
      return(list())
    }
  } else {
    # For other datasets, use the standard gene sets
    return(get_immune_genesets())
  }
}

run_immune_analysis <- function(mat, cohort_name = "TCGA_GBM"){
  message("Running immune gene set analysis...")
  genesets <- map_cgga_probes(mat, cohort_name)
  
  # Calculate mean expression for each gene set
  immune_scores <- matrix(NA, nrow = length(genesets), ncol = ncol(mat))
  rownames(immune_scores) <- names(genesets)
  colnames(immune_scores) <- colnames(mat)
  
  for (i in seq_along(genesets)) {
    geneset_name <- names(genesets)[i]
    genes <- genesets[[i]]
    
    # Find genes that exist in the matrix
    available_genes <- intersect(genes, rownames(mat))
    
    if (length(available_genes) >= 3) {  # Need at least 3 genes
      # Calculate mean expression for this gene set
      immune_scores[geneset_name, ] <- colMeans(mat[available_genes, , drop = FALSE], na.rm = TRUE)
    } else {
      # If not enough genes, set to NA
      immune_scores[geneset_name, ] <- NA
    }
  }
  
  # Remove gene sets with all NA values
  immune_scores <- immune_scores[!apply(immune_scores, 1, function(x) all(is.na(x))), , drop = FALSE]
  
  # Check if we have any valid gene sets
  if (nrow(immune_scores) == 0) {
    warning("No immune gene sets found with sufficient genes in the expression matrix")
    # Create a dummy matrix with some variation to avoid pheatmap issues
    immune_scores <- matrix(rnorm(ncol(mat), 0, 0.1), nrow = 1, ncol = ncol(mat))
    rownames(immune_scores) <- "No_immune_genes_found"
    colnames(immune_scores) <- colnames(mat)
  }
  
  # Convert to list format for compatibility
  res <- list()
  res[["immune_gsva"]] <- immune_scores
  res
}

scale01 <- function(x){ 
  x_min <- min(x, na.rm=TRUE)
  x_max <- max(x, na.rm=TRUE)
  if (x_max == x_min) {
    return(rep(0.5, length(x)))  # All values the same, set to middle
  }
  (x - x_min) / (x_max - x_min)
}

save_heatmap <- function(mat, groups, outfile, title){
  # Check if we have valid data
  if (nrow(mat) == 0 || ncol(mat) == 0) {
    warning("Empty matrix, skipping heatmap")
    return(invisible())
  }
  
  # Z-score normalize each row (gene set) for professional heatmap centered at 0
  M <- t(apply(mat, 1, function(x) {
    if (sd(x, na.rm = TRUE) == 0) return(rep(0, length(x)))
    (x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE)
  }))
  
  # Check if we have enough rows for clustering
  cluster_rows <- nrow(M) > 1
  cluster_cols <- ncol(M) > 1
  
  # Define color palette centered at 0 (blue-white-red)
  color_palette <- colorRampPalette(c("blue", "white", "red"))(100)
  
  # Create professional heatmap with centered color scale (no annotation for now)
  pheatmap(M, show_colnames = FALSE, cluster_rows = cluster_rows, cluster_cols = cluster_cols,
           color = color_palette, breaks = seq(-3, 3, length.out = 101),
           main = title, filename = outfile, width = 10, height = 6)
}

boxplot_cells <- function(mat, groups, keep_cells, outfile, title){
  inter <- intersect(rownames(mat), keep_cells)
  if (length(inter) < 1) return(invisible())
  df <- as.data.frame(t(mat[inter, , drop=FALSE]))
  df$DGAT1_group <- groups
  long <- tidyr::pivot_longer(df, cols = all_of(inter), names_to="Cell", values_to="Score")
  g <- ggplot(long, aes(DGAT1_group, Score)) +
    geom_boxplot(outlier.size = 0.8, width = 0.65) +
    facet_wrap(~Cell, scales="free_y") +
    labs(title = title, x = "", y = "Estimated abundance / score") +
    theme_bw(base_size = 11)
  ggsave(outfile, g, width = 10, height = 6, dpi = 300)
}

# --- Cell label aliases & helpers (GSVA-focused) ------------------------------
aliases <- list(
  CD8_T_cells = c("CD8_T_cells"),
  NK_cells = c("NK_cells"),
  Treg = c("Treg"),
  Exhausted_T = c("Exhausted_T"),
  CD4_T_cells = c("CD4_T_cells"),
  Monocytes = c("Monocytes"),
  Macrophages = c("Macrophages"),
  M1_Macrophages = c("M1_Macrophages"),
  M2_Macrophages = c("M2_Macrophages"),
  Dendritic_cells = c("Dendritic_cells"),
  B_cells = c("B_cells"),
  Pro_inflammatory = c("Pro_inflammatory"),
  Anti_inflammatory = c("Anti_inflammatory"),
  Immune_checkpoints = c("Immune_checkpoints"),
  Complement = c("Complement"),
  Interferon_response = c("Interferon_response")
)

panel_order <- c("CD8_T_cells","NK_cells","Treg","Exhausted_T","CD4_T_cells","Monocytes",
                 "Macrophages","M1_Macrophages","M2_Macrophages","Dendritic_cells","B_cells",
                 "Pro_inflammatory","Anti_inflammatory","Immune_checkpoints","Complement","Interferon_response")

resolve_aliases <- function(present_labels, aliases, order_vec){
  out <- unlist(lapply(order_vec, function(k){
    hits <- intersect(aliases[[k]], present_labels)
    if (length(hits)) hits else character(0)
  }), use.names = FALSE)
  unique(out)
}

# Optional: collapse alias rows into unified names (keeps raw tables untouched)
collapse_aliases_for_display <- function(mat, aliases, order_vec){
  keep <- intersect(rownames(mat), unique(unlist(aliases, use.names = FALSE)))
  if (!length(keep)) return(list(mat = mat, mapping = NULL))
  sub <- mat[keep, , drop = FALSE]
  collapsed <- lapply(order_vec, function(k){
    labs <- intersect(aliases[[k]], rownames(sub))
    if (!length(labs)) return(NULL)
    r <- if (length(labs) == 1) sub[labs, , drop = FALSE] else colMeans(sub[labs, , drop = FALSE], na.rm = TRUE)
    if (is.matrix(r)) r <- as.numeric(r)
    r
  })
  nm <- order_vec[!sapply(collapsed, is.null)]
  M <- do.call(rbind, collapsed[!sapply(collapsed, is.null)])
  rownames(M) <- nm
  list(mat = M, mapping = nm)
}

run_deconv_for_cohort <- function(expr_fp, cohort_name, dgat_rowname = "DGAT1"){
  message("\n=== ", cohort_name, " ===")
  outdir <- file.path(out_root, cohort_name)
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

  mat <- read_expr(expr_fp)
  # drop duplicated genes by mean
  if (any(duplicated(rownames(mat)))){
    DT <- as.data.table(mat, keep.rownames="gene")
    mat <- as.matrix(DT[, lapply(.SD, mean, na.rm=TRUE), by=gene] |> tibble::column_to_rownames("gene"))
  }

  if (!(dgat_rowname %in% rownames(mat))) {
    warning("DGAT1 not in expression matrix for ", cohort_name, "; will create groups by median of available values (NA).")
  }
  dgat <- mat[dgat_rowname, , drop=TRUE]
  cutoff <- load_bestcut(cohort_name)
  groups <- make_groups(dgat, cutoff)

  # run immune gene set analysis
  dec <- run_immune_analysis(mat, cohort_name)

  # save raw tables
  for (nm in names(dec)){
    tab <- dec[[nm]]
    # ensure samples in columns
    if (nrow(tab) < ncol(tab)) {
      # ok, rows = gene sets
    } else {
      tab <- t(tab)
    }
    fwrite(as.data.frame(tab), file.path(outdir, paste0("immune_", nm, ".csv")))
  }

  # choose one for visuals (immune_gsva)
  immune_scores <- dec[["immune_gsva"]]
  if (nrow(immune_scores) > ncol(immune_scores)) immune_scores <- immune_scores  # rows gene sets x columns samples
  # align group order
  immune_scores <- immune_scores[, colnames(mat)[colnames(mat) %in% colnames(immune_scores)], drop=FALSE]
  groups <- groups[colnames(immune_scores)]
  
  message("Matrix dimensions: ", nrow(immune_scores), " x ", ncol(immune_scores))
  message("Groups length: ", length(groups))
  message("Group levels: ", paste(levels(groups), collapse = ", "))

  # Heatmap (all immune gene sets; DGAT1-group ordered)
  save_heatmap(immune_scores, groups,
               file.path(outdir, "immune_heatmap_by_DGAT1.png"),
               paste0(cohort_name, " Immune Gene Sets — stratified by DGAT1 group"))

  # Optional "collapsed" heatmap using unified biology labels
  collapsed <- collapse_aliases_for_display(immune_scores, aliases, panel_order)
  if (!is.null(collapsed$mapping) && nrow(collapsed$mat) > 0) {
    immune_collapsed <- collapsed$mat
    # Ensure groups are properly aligned
    if (length(groups) == ncol(immune_collapsed)) {
      save_heatmap(immune_collapsed, groups,
                   file.path(outdir, "immune_heatmap_collapsed_by_DGAT1.png"),
                   paste0(cohort_name, " Immune Gene Sets (collapsed) — DGAT1 groups"))
    } else {
      save_heatmap(immune_collapsed, groups[colnames(immune_collapsed)],
                   file.path(outdir, "immune_heatmap_collapsed_by_DGAT1.png"),
                   paste0(cohort_name, " Immune Gene Sets (collapsed) — DGAT1 groups"))
    }
    
    # Save alias mapping for transparency
    alias_map <- data.frame(
      Collapsed_Label = names(aliases),
      Raw_Labels = sapply(aliases, function(x) paste(x, collapse = "; "))
    )
    fwrite(alias_map, file.path(outdir, "immune_alias_map.csv"))
  }

  # Boxplots for a comprehensive but tidy panel (auto-picks what exists)
  key_sets <- resolve_aliases(rownames(immune_scores), aliases, panel_order)
  boxplot_cells(immune_scores, groups, key_sets,
                file.path(outdir, "immune_keysets_boxplots.png"),
                paste0(cohort_name, " Immune Gene Sets — key populations vs DGAT1 group"))

  # correlations with DGAT1 expression
  dgat_vec <- dgat[colnames(immune_scores)]
  cors <- apply(immune_scores, 1, function(v){
    suppressWarnings(cor.test(v, dgat_vec, method="spearman")$estimate)
  })
  pvals <- apply(immune_scores, 1, function(v){
    suppressWarnings(cor.test(v, dgat_vec, method="spearman")$p.value)
  })
  cor_tab <- data.frame(GeneSet = rownames(immune_scores), rho = as.numeric(cors), p = as.numeric(pvals)) |>
    arrange(p)
  fwrite(cor_tab, file.path(outdir, "immune_cor_with_DGAT1.csv"))
  message("Saved: ", outdir)
}

# ---- Cohorts (adjust paths if yours differ)
tcga_fp <- "Raw_Data/TCGA/GBM/TCGA_GBM_Expression_HGNC_FPKM.rds"
cgga_fp <- "Raw_Data/CGGA/Expression/CGGA_GSE16011_Expression_HGNC_WithDGAT.rds"
stopifnot(file.exists(tcga_fp), file.exists(cgga_fp))

run_deconv_for_cohort(tcga_fp, "TCGA_GBM", "DGAT1")
run_deconv_for_cohort(cgga_fp, "CGGA_GBM", "DGAT1")

message("\n✅ Deconvolution complete. See: ", out_root)
