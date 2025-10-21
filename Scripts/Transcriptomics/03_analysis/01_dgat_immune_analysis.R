#!/usr/bin/env Rscript
# =============================================================================
# 02_bulk_immune_analysis.R — Comprehensive DGAT-Immune Analysis
# =============================================================================
#
# Analyzes DGAT gene expression and immune cell infiltration relationships.
# Combines GSVA immune deconvolution, differential analysis, and visualization.
#
# Features:
# - 28 curated immune gene sets (Neftel 2019, Wang 2021)
# - Integrates optimal survival cutoffs from Script 01
# - Publication-quality visualizations (volcano, heatmap, boxplots)
# - Differential analysis (Wilcoxon + effect size)
# - Correlation analysis (Spearman)
# - Robust error handling (standalone, no external dependencies)
#
# Usage: Rscript Scripts/04_analysis/02_bulk_immune_analysis.R
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(ggplot2)
  library(ggrepel)
  library(pheatmap)
  library(GSVA)
})

cat("🚀 Starting DGAT-Immune Analysis\n")
cat("=====================================\n\n")

# =============================================================================
# CONFIGURATION
# =============================================================================

CONFIG <- list(
  # Input paths (adjust for your batch-corrected data)
  tcga_expr = "Processed_Data/TCGA_GBM_Batch_Corrected/expression_batch_corrected.rds",
  cgga_expr = "Processed_Data/CGGA_GBM_Batch_Corrected/expression_batch_corrected.rds",
  
  # Output directory
  output_dir = "Results/Bulk/Immune_Analysis",
  
  # Analysis parameters
  dgat_gene = "DGAT1",
  min_geneset_size = 3,
  max_geneset_size = 500,
  correlation_method = "spearman",
  fdr_threshold = 0.05
)

dir.create(CONFIG$output_dir, recursive = TRUE, showWarnings = FALSE)

# =============================================================================
# 1. CURATED IMMUNE GENE SETS
# =============================================================================

get_immune_genesets <- function() {
  cat("Loading curated immune gene sets...\n")
  
  list(
    # T cell subsets (Neftel 2019, Cell)
    CD8_T_cells = c("CD8A", "CD8B", "GZMA", "GZMB", "GZMH", "GZMK", "PRF1", "GNLY"),
    CD4_T_cells = c("CD4", "IL7R", "CCR7", "TCF7", "SELL", "LEF1"),
    Tregs = c("FOXP3", "IL2RA", "CTLA4", "IKZF2", "IL10", "TGFB1"),
    Th1 = c("IFNG", "TBX21", "STAT1", "STAT4", "IL12RB2"),
    Th2 = c("IL4", "IL5", "IL13", "GATA3", "STAT6"),
    Exhausted_T = c("PDCD1", "CTLA4", "LAG3", "HAVCR2", "TIGIT", "TOX", "ENTPD1"),
    
    # NK cells
    NK_cells = c("NCAM1", "KLRD1", "KLRF1", "KLRC1", "GNLY", "NKG7", "GZMA", "GZMB"),
    NK_activated = c("IFNG", "TNF", "GZMB", "PRF1", "FASLG"),
    
    # Myeloid cells (Wang 2021, Cell)
    Monocytes = c("CD14", "FCGR3A", "S100A8", "S100A9", "LYZ", "VCAN"),
    Macrophages = c("CD68", "CD163", "MSR1", "MRC1", "MARCO", "MERTK"),
    M1_Macrophages = c("CD80", "CD86", "IL1B", "IL6", "TNF", "NOS2", "CXCL9", "CXCL10"),
    M2_Macrophages = c("CD163", "MRC1", "MSR1", "IL10", "TGFB1", "ARG1", "CHI3L1"),
    TAMs = c("APOE", "C1QA", "C1QB", "C1QC", "TREM2", "SPP1"),
    MDSCs = c("CD14", "CD33", "S100A8", "S100A9", "ARG1", "IDO1"),
    Dendritic_cells = c("ITGAX", "ITGAE", "CD1C", "CLEC9A", "XCR1", "BATF3"),
    pDCs = c("CLEC4C", "IL3RA", "NRP1", "GZMB", "IRF7", "IRF8"),
    
    # B cells
    B_cells = c("CD19", "MS4A1", "CD79A", "CD79B", "PAX5", "CD22"),
    Plasma_cells = c("SDC1", "CD38", "IGHG1", "IGKC", "MZB1", "XBP1"),
    
    # Functional programs
    Cytotoxicity = c("GZMA", "GZMB", "GZMH", "GZMK", "PRF1", "GNLY", "FASLG"),
    Antigen_presentation = c("HLA-A", "HLA-B", "HLA-C", "B2M", "TAP1", "TAP2", "PSMB9"),
    Pro_inflammatory = c("IL1B", "IL6", "TNF", "IFNG", "IL12A", "IL12B", "IL23A"),
    Anti_inflammatory = c("IL10", "TGFB1", "IL4", "IL13", "ARG1", "IDO1", "VEGFA"),
    
    # Immune checkpoints
    Immune_checkpoints = c("PDCD1", "CD274", "PDCD1LG2", "CTLA4", "LAG3", "HAVCR2", "TIGIT", "VSIR"),
    Costimulatory = c("CD28", "ICOS", "CD27", "CD40LG", "TNFRSF4", "TNFRSF9"),
    
    # Interferon response
    IFN_alpha = c("IFIT1", "IFIT2", "IFIT3", "ISG15", "MX1", "OAS1", "OAS2", "OAS3"),
    IFN_gamma = c("CXCL9", "CXCL10", "CXCL11", "IDO1", "STAT1", "IRF1", "GBP1"),
    
    # Complement
    Complement = c("C1QA", "C1QB", "C1QC", "C3", "C3AR1", "C5AR1", "CFB", "CFD")
  )
}

# =============================================================================
# 2. DATA LOADING & PREPROCESSING
# =============================================================================

load_expression_data <- function(expr_path, cohort_name) {
  cat("\n📖 Loading", cohort_name, "expression data...\n")
  
  if (!file.exists(expr_path)) {
    stop("File not found: ", expr_path)
  }
  
  # Load and ensure proper format
  mat <- readRDS(expr_path)
  if (nrow(mat) < ncol(mat)) mat <- t(mat)
  mat <- as.matrix(mat)
  storage.mode(mat) <- "numeric"
  
  # Handle duplicated genes
  if (any(duplicated(rownames(mat)))) {
    n_dup <- sum(duplicated(rownames(mat)))
    cat("  Aggregating", n_dup, "duplicated genes...\n")
    DT <- as.data.table(mat, keep.rownames = "gene")
    mat <- as.matrix(DT[, lapply(.SD, mean, na.rm = TRUE), by = gene] %>% 
                     tibble::column_to_rownames("gene"))
  }
  
  cat("  ✅ Loaded:", nrow(mat), "genes x", ncol(mat), "samples\n")
  return(mat)
}

load_survival_cutoff <- function(cohort_name) {
  fp <- file.path("Results/Bulk/Survival", cohort_name, "DGAT1_bestcut_info.csv")
  
  if (file.exists(fp)) {
    info <- fread(fp)
    cutoff <- as.numeric(info$cutoff[1])
    cat("  ✅ Using optimal survival cutoff:", round(cutoff, 3), "\n")
    return(cutoff)
  } else {
    cat("  ⚠️  No survival cutoff found, will use median\n")
    return(NA_real_)
  }
}

create_dgat_groups <- function(expr_vec, cutoff = NA_real_) {
  n_na <- sum(is.na(expr_vec))
  if (n_na > 0) cat("  Warning:", n_na, "samples with NA DGAT expression\n")
  
  if (!is.na(cutoff)) {
    groups <- factor(ifelse(expr_vec >= cutoff, "High", "Low"), 
                    levels = c("Low", "High"))
    cat("  Using optimal cutoff:", round(cutoff, 3), "\n")
  } else {
    med <- median(expr_vec, na.rm = TRUE)
    groups <- factor(ifelse(expr_vec >= med, "High", "Low"), 
                    levels = c("Low", "High"))
    cat("  Using median cutoff:", round(med, 3), "\n")
  }
  
  cat("  Groups: Low =", sum(groups == "Low", na.rm = TRUE), 
      "| High =", sum(groups == "High", na.rm = TRUE), "\n")
  
  return(groups)
}

# =============================================================================
# 3. GSVA IMMUNE DECONVOLUTION
# =============================================================================

run_gsva_analysis <- function(expr_mat, genesets, cohort_name) {
  cat("\n🔬 Running GSVA analysis for", cohort_name, "...\n")
  
  # Filter gene sets to those with sufficient genes
  filtered_genesets <- list()
  for (gs_name in names(genesets)) {
    genes_present <- intersect(genesets[[gs_name]], rownames(expr_mat))
    
    if (length(genes_present) >= CONFIG$min_geneset_size) {
      filtered_genesets[[gs_name]] <- genes_present
      cat("  ✓", gs_name, ":", length(genes_present), "genes\n")
    } else {
      cat("  ✗", gs_name, ": skipped (only", length(genes_present), "genes)\n")
    }
  }
  
  if (length(filtered_genesets) == 0) {
    stop("No gene sets with sufficient genes!")
  }
  
  cat("\n  Running GSVA on", length(filtered_genesets), "gene sets...\n")
  
  # Run GSVA with error handling
  gsva_scores <- tryCatch({
    gsva(expr_mat, filtered_genesets,
         method = "gsva",
         kcdf = "Gaussian",
         min.sz = CONFIG$min_geneset_size,
         max.sz = CONFIG$max_geneset_size,
         parallel.sz = 1)
  }, error = function(e) {
    cat("  ⚠️  GSVA failed, using mean expression approach:", e$message, "\n")
    
    # Fallback: mean expression
    scores <- matrix(NA, nrow = length(filtered_genesets), ncol = ncol(expr_mat))
    rownames(scores) <- names(filtered_genesets)
    colnames(scores) <- colnames(expr_mat)
    
    for (i in seq_along(filtered_genesets)) {
      genes <- filtered_genesets[[i]]
      scores[i, ] <- colMeans(expr_mat[genes, , drop = FALSE], na.rm = TRUE)
    }
    
    scores
  })
  
  cat("  ✅ GSVA complete:", nrow(gsva_scores), "gene sets scored\n")
  return(gsva_scores)
}

# =============================================================================
# 4. DIFFERENTIAL ANALYSIS
# =============================================================================

perform_differential_analysis <- function(gsva_scores, groups) {
  cat("\n📊 Performing differential analysis...\n")
  
  stopifnot(length(groups) == ncol(gsva_scores))
  
  results <- lapply(rownames(gsva_scores), function(gs) {
    vals <- gsva_scores[gs, ]
    g_high <- vals[groups == "High"]
    g_low <- vals[groups == "Low"]
    
    if (length(g_high) < 3 || length(g_low) < 3) return(NULL)
    
    # Wilcoxon test
    wt <- suppressWarnings(wilcox.test(g_high, g_low))
    
    # Effect size (Cohen's d)
    pooled_sd <- sqrt((var(g_high, na.rm = TRUE) + var(g_low, na.rm = TRUE)) / 2)
    effect_size <- (mean(g_high, na.rm = TRUE) - mean(g_low, na.rm = TRUE)) / pooled_sd
    
    data.frame(
      GeneSet = gs,
      Mean_High = mean(g_high, na.rm = TRUE),
      Mean_Low = mean(g_low, na.rm = TRUE),
      Delta = mean(g_high, na.rm = TRUE) - mean(g_low, na.rm = TRUE),
      Effect_Size = effect_size,
      p_value = wt$p.value,
      stringsAsFactors = FALSE
    )
  })
  
  results <- do.call(rbind, results[!sapply(results, is.null)])
  results$FDR <- p.adjust(results$p_value, method = "BH")
  results <- results[order(results$FDR, -abs(results$Delta)), ]
  
  # Print top results
  cat("\n  Top enriched in DGAT High (FDR <", CONFIG$fdr_threshold, "):\n")
  top_high <- results %>% filter(FDR < CONFIG$fdr_threshold, Delta > 0) %>% head(5)
  if (nrow(top_high) > 0) {
    for (i in 1:nrow(top_high)) {
      cat("    ", top_high$GeneSet[i], ": Δ =", round(top_high$Delta[i], 3), 
          ", FDR =", format(top_high$FDR[i], digits = 3), "\n")
    }
  } else {
    cat("    (none significant)\n")
  }
  
  cat("\n  Top enriched in DGAT Low (FDR <", CONFIG$fdr_threshold, "):\n")
  top_low <- results %>% filter(FDR < CONFIG$fdr_threshold, Delta < 0) %>% head(5)
  if (nrow(top_low) > 0) {
    for (i in 1:nrow(top_low)) {
      cat("    ", top_low$GeneSet[i], ": Δ =", round(top_low$Delta[i], 3), 
          ", FDR =", format(top_low$FDR[i], digits = 3), "\n")
    }
  } else {
    cat("    (none significant)\n")
  }
  
  return(results)
}

# =============================================================================
# 5. CORRELATION ANALYSIS
# =============================================================================

perform_correlation_analysis <- function(gsva_scores, dgat_expr) {
  cat("\n📈 Performing correlation analysis...\n")
  
  # Align samples
  common_samples <- intersect(colnames(gsva_scores), names(dgat_expr))
  gsva_scores <- gsva_scores[, common_samples, drop = FALSE]
  dgat_expr <- dgat_expr[common_samples]
  
  # Calculate correlations
  cor_results <- lapply(rownames(gsva_scores), function(gs) {
    vals <- gsva_scores[gs, ]
    
    if (sd(vals, na.rm = TRUE) == 0) return(NULL)
    
    ct <- cor.test(vals, dgat_expr, method = CONFIG$correlation_method)
    
    data.frame(
      GeneSet = gs,
      rho = as.numeric(ct$estimate),
      p_value = ct$p.value,
      stringsAsFactors = FALSE
    )
  })
  
  cor_results <- do.call(rbind, cor_results[!sapply(cor_results, is.null)])
  cor_results$FDR <- p.adjust(cor_results$p_value, method = "BH")
  cor_results <- cor_results[order(cor_results$FDR, -abs(cor_results$rho)), ]
  
  cat("  ✅ Calculated correlations for", nrow(cor_results), "gene sets\n")
  
  return(cor_results)
}

# =============================================================================
# 6. VISUALIZATION FUNCTIONS
# =============================================================================

create_volcano_plot <- function(diff_results, outfile, title) {
  cat("  Creating volcano plot...\n")
  
  if (is.null(diff_results) || nrow(diff_results) == 0) return(invisible())
  
  # Add significance labels
  diff_results$Significance <- dplyr::case_when(
    diff_results$FDR < CONFIG$fdr_threshold & diff_results$Delta > 0 ~ "Up in High",
    diff_results$FDR < CONFIG$fdr_threshold & diff_results$Delta < 0 ~ "Up in Low",
    TRUE ~ "NS"
  )
  
  # Label top genes
  diff_results$Label <- ""
  top_up <- diff_results %>% filter(Significance == "Up in High") %>% head(5)
  top_down <- diff_results %>% filter(Significance == "Up in Low") %>% head(5)
  diff_results$Label[diff_results$GeneSet %in% c(top_up$GeneSet, top_down$GeneSet)] <-
    diff_results$GeneSet[diff_results$GeneSet %in% c(top_up$GeneSet, top_down$GeneSet)]
  
  # Create plot
  p <- ggplot(diff_results, aes(x = Delta, y = -log10(FDR))) +
    geom_point(aes(color = Significance), size = 2, alpha = 0.7) +
    scale_color_manual(values = c("Up in High" = "#D73027",
                                  "Up in Low" = "#4575B4",
                                  "NS" = "grey70")) +
    geom_hline(yintercept = -log10(CONFIG$fdr_threshold), 
               linetype = "dashed", alpha = 0.5) +
    geom_vline(xintercept = 0, linetype = "dashed", alpha = 0.5) +
    geom_text_repel(aes(label = Label), size = 3, max.overlaps = 15) +
    labs(title = title,
         x = "GSVA Score Difference (High - Low DGAT1)",
         y = "-log10(FDR)",
         color = "Enrichment") +
    theme_minimal(base_size = 12) +
    theme(legend.position = c(0.85, 0.85),
          plot.title = element_text(face = "bold", hjust = 0.5))
  
  ggsave(outfile, p, width = 8, height = 6, dpi = 300)
}

create_heatmap <- function(gsva_scores, groups, outfile, title) {
  cat("  Creating heatmap...\n")
  
  # Z-score normalization
  mat_z <- t(scale(t(gsva_scores)))
  
  # Order by groups
  ord <- order(groups)
  mat_z <- mat_z[, ord]
  groups_ordered <- groups[ord]
  
  # Annotation
  ann <- data.frame(DGAT1 = groups_ordered, row.names = colnames(mat_z))
  ann_colors <- list(DGAT1 = c(Low = "#4575B4", High = "#D73027"))
  
  # Cap extreme values
  mat_z[mat_z > 2.5] <- 2.5
  mat_z[mat_z < -2.5] <- -2.5
  
  # Create heatmap
  pheatmap(mat_z,
           annotation_col = ann,
           annotation_colors = ann_colors,
           show_colnames = FALSE,
           cluster_cols = FALSE,
           clustering_distance_rows = "correlation",
           clustering_method = "ward.D2",
           color = colorRampPalette(c("#4575B4", "white", "#D73027"))(100),
           breaks = seq(-2.5, 2.5, length.out = 101),
           main = title,
           filename = outfile,
           width = 12,
           height = 8)
}

create_boxplots <- function(gsva_scores, groups, outfile, title) {
  cat("  Creating boxplots for key signatures...\n")
  
  # Key signatures to plot
  key_sets <- c("CD8_T_cells", "Tregs", "Exhausted_T", "M1_Macrophages",
                "M2_Macrophages", "TAMs", "MDSCs", "Immune_checkpoints")
  
  available_sets <- intersect(key_sets, rownames(gsva_scores))
  
  if (length(available_sets) == 0) {
    cat("    No key signatures found, skipping boxplots\n")
    return(invisible())
  }
  
  # Prepare data
  plot_data <- data.frame(t(gsva_scores[available_sets, , drop = FALSE]), 
                         DGAT1 = groups) %>%
    pivot_longer(cols = -DGAT1, names_to = "GeneSet", values_to = "Score")
  
  # Calculate statistics
  stat_data <- plot_data %>% 
    group_by(GeneSet) %>%
    summarise(
      p_value = wilcox.test(Score ~ DGAT1)$p.value,
      max_score = max(Score),
      .groups = 'drop'
    ) %>%
    mutate(
      label = case_when(
        p_value < 0.001 ~ "***",
        p_value < 0.01 ~ "**",
        p_value < 0.05 ~ "*",
        TRUE ~ ""
      ),
      y_pos = max_score * 1.1
    )
  
  # Create plot
  p <- ggplot(plot_data, aes(x = DGAT1, y = Score, fill = DGAT1)) +
    geom_boxplot(alpha = 0.7, outlier.size = 1) +
    geom_jitter(width = 0.2, alpha = 0.3, size = 0.8) +
    geom_text(data = stat_data, aes(x = 1.5, y = y_pos, label = label), 
              inherit.aes = FALSE, size = 4) +
    facet_wrap(~GeneSet, scales = "free_y", ncol = 4) +
    scale_fill_manual(values = c(Low = "#4575B4", High = "#D73027")) +
    labs(title = title, x = "", y = "GSVA Score") +
    theme_minimal(base_size = 11) +
    theme(legend.position = "none",
          strip.text = element_text(face = "bold"),
          plot.title = element_text(face = "bold", hjust = 0.5))
  
  ggsave(outfile, p, width = 12, height = 8, dpi = 300)
}

# =============================================================================
# 7. MAIN ANALYSIS WORKFLOW
# =============================================================================

run_immune_analysis <- function(expr_path, cohort_name) {
  cat("\n", strrep("=", 70), "\n")
  cat("🔬 Running immune analysis for", cohort_name, "\n")
  cat(strrep("=", 70), "\n")
  
  # Create output directory
  outdir <- file.path(CONFIG$output_dir, cohort_name)
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  
  # 1. Load data
  expr_mat <- load_expression_data(expr_path, cohort_name)
  
  # 2. Check for DGAT gene
  if (!(CONFIG$dgat_gene %in% rownames(expr_mat))) {
    warning(CONFIG$dgat_gene, " not found in ", cohort_name)
    return(NULL)
  }
  
  dgat_expr <- expr_mat[CONFIG$dgat_gene, , drop = TRUE]
  
  # 3. Load survival cutoff and create groups
  cutoff <- load_survival_cutoff(cohort_name)
  groups <- create_dgat_groups(dgat_expr, cutoff)
  
  # 4. Get immune gene sets and run GSVA
  immune_genesets <- get_immune_genesets()
  gsva_scores <- run_gsva_analysis(expr_mat, immune_genesets, cohort_name)
  
  if (is.null(gsva_scores)) {
    warning("GSVA failed for ", cohort_name)
    return(NULL)
  }
  
  # 5. Align samples
  common_samples <- intersect(colnames(expr_mat), colnames(gsva_scores))
  gsva_scores <- gsva_scores[, common_samples, drop = FALSE]
  groups <- groups[common_samples]
  dgat_expr <- dgat_expr[common_samples]
  
  # 6. Differential analysis
  diff_results <- perform_differential_analysis(gsva_scores, groups)
  
  # 7. Correlation analysis
  cor_results <- perform_correlation_analysis(gsva_scores, dgat_expr)
  
  # 8. Save results
  cat("\n💾 Saving results...\n")
  fwrite(as.data.frame(gsva_scores) %>% tibble::rownames_to_column("GeneSet"),
         file.path(outdir, "gsva_scores.csv"))
  fwrite(diff_results, file.path(outdir, "differential_analysis.csv"))
  fwrite(cor_results, file.path(outdir, "correlation_analysis.csv"))
  
  # 9. Create visualizations
  cat("\n📊 Creating visualizations...\n")
  
  create_volcano_plot(diff_results,
                     file.path(outdir, "volcano_plot.png"),
                     paste0(cohort_name, " - Immune Programs (DGAT1 High vs Low)"))
  
  create_heatmap(gsva_scores, groups,
                file.path(outdir, "heatmap_gsva.png"),
                paste0(cohort_name, " - GSVA Immune Signatures"))
  
  create_boxplots(gsva_scores, groups,
                 file.path(outdir, "boxplots_key_signatures.png"),
                 paste0(cohort_name, " - Key Immune Programs"))
  
  cat("\n✅ Analysis complete for", cohort_name, "\n")
  cat("   Output directory:", outdir, "\n")
  
  return(list(
    gsva_scores = gsva_scores,
    differential = diff_results,
    correlations = cor_results,
    groups = groups
  ))
}

# =============================================================================
# 8. EXECUTE ANALYSIS
# =============================================================================

cat("\n", strrep("=", 70), "\n")
cat("DGAT-IMMUNE ANALYSIS PIPELINE\n")
cat(strrep("=", 70), "\n\n")

# Run for TCGA
if (file.exists(CONFIG$tcga_expr)) {
  tcga_results <- run_immune_analysis(CONFIG$tcga_expr, "TCGA_GBM")
} else {
  cat("⚠️  TCGA data not found:", CONFIG$tcga_expr, "\n")
}

# Run for CGGA
if (file.exists(CONFIG$cgga_expr)) {
  cgga_results <- run_immune_analysis(CONFIG$cgga_expr, "CGGA_GBM")
} else {
  cat("⚠️  CGGA data not found:", CONFIG$cgga_expr, "\n")
}

cat("\n", strrep("=", 70), "\n")
cat("🎉 ANALYSIS COMPLETE!\n")
cat(strrep("=", 70), "\n")
cat("Results saved to:", CONFIG$output_dir, "\n\n")
cat("Output files per cohort:\n")
cat("  - gsva_scores.csv (all immune signature scores)\n")
cat("  - differential_analysis.csv (DGAT High vs Low)\n")
cat("  - correlation_analysis.csv (DGAT expression correlations)\n")
cat("  - volcano_plot.png (differential enrichment)\n")
cat("  - heatmap_gsva.png (all signatures)\n")
cat("  - boxplots_key_signatures.png (key programs)\n")
cat(strrep("=", 70), "\n")


