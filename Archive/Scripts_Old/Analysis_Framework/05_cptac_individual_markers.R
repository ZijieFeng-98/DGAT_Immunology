#!/usr/bin/env Rscript
# =============================================================================
# CPTAC GBM: Individual Immune Marker Analysis - PRODUCTION VERSION
# Robust handling of proteomics data with comprehensive QC
# =============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(ggrepel)
})

# Optional packages
if (requireNamespace("ComplexHeatmap", quietly = TRUE)) {
  library(ComplexHeatmap)
  library(circlize)
  COMPLEX_HEATMAP_AVAILABLE <- TRUE
} else {
  COMPLEX_HEATMAP_AVAILABLE <- FALSE
  cat("ComplexHeatmap not available - will skip advanced heatmap\n")
}

# ====== CONFIG =================================================================
BASE_DIR <- "/Users/zijiefeng/Desktop/Guo's lab/My_Research/DGAT_Immunology"
DATA_DIR <- file.path(BASE_DIR, "Raw_Data", "Proteome", "CPTAC")
OUT_DIR  <- file.path(BASE_DIR, "Results", "CPTAC_Individual_Markers")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

PROTEIN_FILE <- file.path(DATA_DIR, "CPTAC3_Glioblastoma_Multiforme_Proteome.tmt11.tsv")

# QC thresholds
MIN_DETECTION_RATE <- 0.30  # Require 30% non-NA values
WINSORIZE_PERCENTILE <- 0.01  # Trim 1% tails for outliers
USE_PARTIAL_CORRELATION <- FALSE  # Set TRUE if purity proxy available

# Individual immune markers with cell-type annotations
MARKER_DEFINITIONS <- data.frame(
  Protein = c(
    # Macrophages
    "CD68", "CD163", "MRC1", "MSR1", "MARCO", "LGALS3", "APOE",
    # T cells
    "CD8A", "CD8B", "CD4", "CD3E", "CD3D", "CD3G",
    # Tregs
    "FOXP3", "IL2RA", "CTLA4", "TIGIT",
    # TAM-specific
    "SPP1", "ARG1", "IL10",
    # NK cells  
    "NCAM1", "KLRD1", "KLRF1", "NCR1",
    # Dendritic cells
    "CD1C", "FCER1A", "CLEC9A", "BATF3",
    # B cells
    "CD19", "MS4A1", "CD79A", "CD79B",
    # Neutrophils
    "FCGR3B", "CSF3R", "S100A12",
    # Checkpoint
    "PDCD1", "CD274", "LAG3", "HAVCR2"
  ),
  CellType = c(
    rep("Macrophage", 7), rep("T cell", 6), rep("Treg", 4),
    rep("TAM", 3), rep("NK cell", 4), rep("Dendritic cell", 4),
    rep("B cell", 4), rep("Neutrophil", 3), rep("Checkpoint", 4)
  ),
  Description = c(
    "Pan-macrophage marker", "M2 macrophage marker", "Mannose receptor (CD206)", 
    "Scavenger receptor", "Scavenger receptor", "Galectin-3", "Apolipoprotein E",
    "CD8+ T cell marker", "CD8+ T cell marker", "CD4+ T cell marker",
    "T cell receptor component", "T cell receptor component", "T cell receptor component",
    "Treg master transcription factor", "IL-2 receptor alpha (CD25)", 
    "Immune checkpoint", "Inhibitory receptor",
    "Osteopontin (TAM marker)", "Arginase-1 (immunosuppressive)", "IL-10 (immunosuppressive)",
    "Neural cell adhesion molecule (CD56)", "Killer cell lectin receptor",
    "Killer cell lectin receptor", "Natural cytotoxicity receptor",
    "DC marker", "DC marker", "Cross-presenting DC marker", "DC transcription factor",
    "B cell marker", "B cell marker (CD20)", "B cell receptor component", 
    "B cell receptor component",
    "Neutrophil Fc receptor (CD16b)", "G-CSF receptor", "Neutrophil marker",
    "PD-1", "PD-L1", "LAG3", "TIM-3"
  ),
  stringsAsFactors = FALSE
)

# =============================================================================
# UTILITY FUNCTIONS
# =============================================================================

# FIX 1: Robust file loading with duplicate handling
load_protein_data <- function(filepath) {
  cat("Loading proteomics data...\n")
  prot <- read_tsv(filepath, show_col_types = FALSE)
  
  # Identify gene column
  gene_cols <- c("Gene", "gene", "Gene_Symbol")
  gene_col <- gene_cols[gene_cols %in% names(prot)][1]
  
  if (is.na(gene_col)) {
    stop("Could not identify gene symbol column")
  }
  
  # Identify sample columns (numeric data only)
  meta_cols <- c(gene_col, "Gene_ID", "Protein_ID", "Description")
  sample_cols <- setdiff(names(prot), meta_cols)
  
  # Coerce to numeric, warn about non-numeric
  numeric_cols <- sample_cols[sapply(prot[sample_cols], function(x) {
    is.numeric(x) || all(is.na(suppressWarnings(as.numeric(as.character(x)))))
  })]
  
  cat("  Found", length(numeric_cols), "numeric sample columns\n")
  
  # Extract matrix and convert to numeric
  prot_mat <- as.matrix(prot[, numeric_cols])
  rownames(prot_mat) <- prot[[gene_col]]
  
  # Convert to numeric, handling any non-numeric values
  prot_mat <- apply(prot_mat, 2, function(x) {
    as.numeric(as.character(x))
  })
  rownames(prot_mat) <- prot[[gene_col]]
  
  # Remove non-gene rows (like "Mean", "Median", "StdDev")
  non_gene_rows <- c("Mean", "Median", "StdDev", "Gene")
  prot_mat <- prot_mat[!rownames(prot_mat) %in% non_gene_rows, , drop = FALSE]
  
  # FIX 1: Handle duplicates by taking MEDIAN across duplicated genes
  if (any(duplicated(rownames(prot_mat)))) {
    dup_genes <- unique(rownames(prot_mat)[duplicated(rownames(prot_mat))])
    cat("  WARNING:", length(dup_genes), "duplicated gene symbols detected\n")
    cat("  Collapsing duplicates using MEDIAN...\n")
    
    prot_df <- as.data.frame(prot_mat)
    prot_df$Gene <- rownames(prot_mat)
    
    prot_collapsed <- prot_df %>%
      group_by(Gene) %>%
      summarise(across(everything(), ~median(.x, na.rm = TRUE)), .groups = "drop")
    
    prot_mat <- as.matrix(prot_collapsed[, -1])
    rownames(prot_mat) <- prot_collapsed$Gene
  }
  
  cat("  Final matrix:", nrow(prot_mat), "proteins x", ncol(prot_mat), "samples\n\n")
  
  prot_mat
}

# FIX 3: DGAT1 distribution check with optional winsorization
check_dgat1_distribution <- function(dgat_vec, winsorize = TRUE, percentile = 0.01) {
  cat("\n=== DGAT1 Distribution Check ===\n")
  cat("Samples with data:", sum(!is.na(dgat_vec)), "/", length(dgat_vec), "\n")
  cat("Range: [", round(min(dgat_vec, na.rm = TRUE), 3), ",", 
      round(max(dgat_vec, na.rm = TRUE), 3), "]\n")
  cat("Mean ± SD:", round(mean(dgat_vec, na.rm = TRUE), 3), "±", 
      round(sd(dgat_vec, na.rm = TRUE), 3), "\n")
  cat("Median [IQR]:", round(median(dgat_vec, na.rm = TRUE), 3), "[",
      round(quantile(dgat_vec, 0.25, na.rm = TRUE), 3), "-",
      round(quantile(dgat_vec, 0.75, na.rm = TRUE), 3), "]\n")
  
  if (winsorize && percentile > 0) {
    lower_bound <- quantile(dgat_vec, percentile, na.rm = TRUE)
    upper_bound <- quantile(dgat_vec, 1 - percentile, na.rm = TRUE)
    
    n_winsorized <- sum(dgat_vec < lower_bound | dgat_vec > upper_bound, na.rm = TRUE)
    
    if (n_winsorized > 0) {
      cat("\nWinsorizing", n_winsorized, "extreme values to [", 
          round(lower_bound, 3), ",", round(upper_bound, 3), "]\n")
      
      dgat_vec[dgat_vec < lower_bound] <- lower_bound
      dgat_vec[dgat_vec > upper_bound] <- upper_bound
    }
  }
  
  cat("================================\n\n")
  
  dgat_vec
}

# FIX 2: Detection filtering for sparse proteomics data
filter_markers_by_detection <- function(prot_mat, markers_df, min_rate = 0.30) {
  cat("Filtering markers by detection rate (>", min_rate * 100, "%)...\n")
  
  present <- markers_df$Protein[markers_df$Protein %in% rownames(prot_mat)]
  missing <- markers_df$Protein[!markers_df$Protein %in% rownames(prot_mat)]
  
  if (length(present) == 0) {
    stop("No immune markers found in proteomics data")
  }
  
  # Calculate detection rates
  detection_rates <- sapply(present, function(p) {
    sum(!is.na(prot_mat[p, ])) / ncol(prot_mat)
  })
  
  passed <- present[detection_rates >= min_rate]
  failed <- present[detection_rates < min_rate]
  
  cat("  Total markers in definition:", nrow(markers_df), "\n")
  cat("  Found in data:", length(present), "\n")
  cat("  Passed detection filter:", length(passed), "\n")
  cat("  Failed detection filter:", length(failed), "\n")
  cat("  Missing from data:", length(missing), "\n\n")
  
  # FIX 4: Save coverage reports
  present_df <- markers_df %>%
    filter(Protein %in% passed) %>%
    mutate(DetectionRate = detection_rates[Protein])
  
  failed_df <- markers_df %>%
    filter(Protein %in% failed) %>%
    mutate(DetectionRate = detection_rates[Protein])
  
  missing_df <- markers_df %>%
    filter(Protein %in% missing) %>%
    mutate(DetectionRate = NA)
  
  list(
    passed = passed,
    present_df = present_df,
    failed_df = failed_df,
    missing_df = missing_df
  )
}

# Individual marker analysis
analyze_individual_markers <- function(prot_mat, dgat_vec, passed_markers, 
                                       markers_df, use_partial = FALSE, 
                                       purity_proxy = NULL) {
  
  cat("Analyzing", length(passed_markers), "markers...\n")
  
  results <- lapply(passed_markers, function(protein) {
    vec <- prot_mat[protein, ]
    valid <- !is.na(vec) & !is.na(dgat_vec)
    
    if (sum(valid) < 10) return(NULL)
    
    # Standard Spearman correlation
    test <- cor.test(vec[valid], dgat_vec[valid], method = "spearman")
    
    # FIX 7: Optional partial correlation (if purity proxy provided)
    partial_rho <- NA
    partial_p <- NA
    
    if (use_partial && !is.null(purity_proxy)) {
      if (requireNamespace("ppcor", quietly = TRUE)) {
        purity_aligned <- purity_proxy[valid]
        if (sum(!is.na(purity_aligned)) >= 10) {
          pcor_data <- data.frame(
            dgat = dgat_vec[valid],
            marker = vec[valid],
            purity = purity_aligned
          )
          pcor_data <- na.omit(pcor_data)
          
          if (nrow(pcor_data) >= 10) {
            pcor_result <- ppcor::pcor.test(
              pcor_data$dgat, 
              pcor_data$marker, 
              pcor_data$purity
            )
            partial_rho <- pcor_result$estimate
            partial_p <- pcor_result$p.value
          }
        }
      }
    }
    
    marker_info <- markers_df[markers_df$Protein == protein, ]
    
    data.frame(
      Protein = protein,
      CellType = marker_info$CellType,
      Description = marker_info$Description,
      Rho = unname(test$estimate),
      P = test$p.value,
      N = sum(valid),
      Partial_Rho = partial_rho,
      Partial_P = partial_p,
      stringsAsFactors = FALSE
    )
  })
  
  tab <- do.call(rbind, results[!sapply(results, is.null)])
  tab$FDR <- p.adjust(tab$P, method = "BH")
  tab$Significant <- tab$FDR < 0.05
  
  if (use_partial && any(!is.na(tab$Partial_P))) {
    tab$Partial_FDR <- p.adjust(tab$Partial_P, method = "BH")
  }
  
  tab[order(tab$FDR), ]
}

# =============================================================================
# PLOTTING FUNCTIONS
# =============================================================================

# FIX 5: Volcano with -log10(P) on y-axis, FDR for coloring
plot_volcano_colored <- function(cor_results, outfile) {
  
  cell_colors <- c(
    "Macrophage" = "#E41A1C", "TAM" = "#984EA3", 
    "T cell" = "#4DAF4A", "Treg" = "#377EB8",
    "NK cell" = "#FF7F00", "Dendritic cell" = "#FFFF33",
    "B cell" = "#A65628", "Neutrophil" = "#F781BF",
    "Checkpoint" = "#999999"
  )
  
  cor_results$CellType <- factor(cor_results$CellType, levels = names(cell_colors))
  
  # FIX 5: Label by FDR < 0.1 OR extreme correlation
  cor_results$Label <- ifelse(
    cor_results$FDR < 0.1 | abs(cor_results$Rho) > 0.15,
    cor_results$Protein,
    ""
  )
  
  # Color by FDR significance
  cor_results$Sig_Level <- case_when(
    cor_results$FDR < 0.001 ~ "FDR < 0.001",
    cor_results$FDR < 0.01 ~ "FDR < 0.01",
    cor_results$FDR < 0.05 ~ "FDR < 0.05",
    cor_results$P < 0.05 ~ "P < 0.05",
    TRUE ~ "NS"
  )
  cor_results$Sig_Level <- factor(cor_results$Sig_Level, 
                                   levels = c("FDR < 0.001", "FDR < 0.01", 
                                             "FDR < 0.05", "P < 0.05", "NS"))
  
  p <- ggplot(cor_results, aes(x = Rho, y = -log10(P))) +
    geom_point(aes(color = CellType, size = Sig_Level, alpha = Sig_Level)) +
    geom_text_repel(
      aes(label = Label),
      size = 3,
      max.overlaps = 20,
      box.padding = 0.5,
      segment.color = "grey50"
    ) +
    scale_color_manual(values = cell_colors, name = "Cell Type") +
    scale_size_manual(
      values = c("FDR < 0.001" = 4, "FDR < 0.01" = 3.5, "FDR < 0.05" = 3, 
                 "P < 0.05" = 2.5, "NS" = 2),
      name = "Significance"
    ) +
    scale_alpha_manual(
      values = c("FDR < 0.001" = 1, "FDR < 0.01" = 0.9, "FDR < 0.05" = 0.8, 
                 "P < 0.05" = 0.6, "NS" = 0.4),
      name = "Significance"
    ) +
    geom_vline(xintercept = 0, linetype = "solid", color = "black") +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "blue", alpha = 0.5) +
    annotate("text", x = Inf, y = -log10(0.05), 
             label = "P = 0.05", hjust = 1.1, vjust = -0.5, 
             size = 3, color = "blue") +
    labs(
      title = "DGAT1 Protein vs Individual Immune Markers",
      subtitle = "Y-axis: nominal P-value | Point size/color: FDR-adjusted significance",
      x = "Spearman ρ (DGAT1 correlation)",
      y = "-log10(P-value)"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      legend.position = "right",
      plot.title = element_text(face = "bold", size = 14)
    )
  
  ggsave(outfile, p, width = 14, height = 8, dpi = 300)
  
  return(p)
}

# Dot plot (unchanged, already good)
plot_dotplot_annotated <- function(cor_results, outfile) {
  
  cor_results <- cor_results %>%
    arrange(Rho) %>%
    mutate(Protein = factor(Protein, levels = Protein))
  
  cor_results$SigLabel <- case_when(
    cor_results$FDR < 0.001 ~ "***",
    cor_results$FDR < 0.01 ~ "**",
    cor_results$FDR < 0.05 ~ "*",
    cor_results$P < 0.05 ~ "†",
    TRUE ~ ""
  )
  
  p <- ggplot(cor_results, aes(y = Protein, x = Rho)) +
    geom_segment(aes(yend = Protein, xend = 0, color = CellType), 
                 linewidth = 0.8) +
    geom_point(aes(color = CellType, size = -log10(P))) +
    geom_text(aes(label = SigLabel), 
              hjust = -0.5, size = 5, fontface = "bold") +
    geom_vline(xintercept = 0, linetype = "solid", color = "black") +
    scale_color_manual(
      values = c(
        "Macrophage" = "#E41A1C", "TAM" = "#984EA3", 
        "T cell" = "#4DAF4A", "Treg" = "#377EB8",
        "NK cell" = "#FF7F00", "Dendritic cell" = "#FFFF33",
        "B cell" = "#A65628", "Neutrophil" = "#F781BF",
        "Checkpoint" = "#999999"
      )
    ) +
    scale_size_continuous(range = c(2, 6), name = "-log10(P)") +
    labs(
      title = "DGAT1 Correlation with Individual Immune Markers",
      subtitle = "† p<0.05, * FDR<0.05, ** FDR<0.01, *** FDR<0.001",
      x = "Spearman ρ",
      y = "",
      color = "Cell Type"
    ) +
    theme_minimal(base_size = 11) +
    theme(
      axis.text.y = element_text(size = 9),
      legend.position = "right",
      panel.grid.major.y = element_line(color = "grey90")
    )
  
  ggsave(outfile, p, width = 10, height = max(6, nrow(cor_results) * 0.25), 
         dpi = 300, limitsize = FALSE)
  
  return(p)
}

# FIX 6: Safer heatmap with alignment checks
plot_clustered_heatmap <- function(prot_mat, dgat_vec, cor_results, outfile) {
  
  if (!COMPLEX_HEATMAP_AVAILABLE) {
    cat("  ⚠ ComplexHeatmap not available, skipping heatmap\n")
    return(NULL)
  }
  
  present <- cor_results$Protein[cor_results$Protein %in% rownames(prot_mat)]
  
  if (length(present) < 2) {
    cat("  ⚠ Too few markers for heatmap (need ≥2)\n")
    return(NULL)
  }
  
  # FIX 6: Scale rows, remove all-NA samples
  marker_mat <- prot_mat[present, , drop = FALSE]
  marker_z <- t(scale(t(marker_mat)))
  
  valid_samples <- colSums(!is.na(marker_z)) > 0
  marker_z <- marker_z[, valid_samples, drop = FALSE]
  dgat_aligned <- dgat_vec[colnames(marker_z)]
  
  # FIX 6: Ensure row annotations match exactly
  row_anno <- cor_results %>%
    filter(Protein %in% rownames(marker_z)) %>%
    arrange(match(Protein, rownames(marker_z)))
  
  # Verify alignment
  if (!all(row_anno$Protein == rownames(marker_z))) {
    warning("Row annotation mismatch, reordering...")
    row_anno <- row_anno[match(rownames(marker_z), row_anno$Protein), ]
  }
  
  ha_row <- rowAnnotation(
    CellType = row_anno$CellType,
    Correlation = anno_barplot(
      row_anno$Rho,
      border = FALSE,
      gp = gpar(fill = ifelse(row_anno$Rho > 0, "#D73027", "#4575B4"))
    ),
    col = list(
      CellType = c(
        "Macrophage" = "#E41A1C", "TAM" = "#984EA3",
        "T cell" = "#4DAF4A", "Treg" = "#377EB8",
        "NK cell" = "#FF7F00", "Checkpoint" = "#999999"
      )
    )
  )
  
  dgat_groups <- ifelse(dgat_aligned >= median(dgat_aligned, na.rm = TRUE), 
                        "High", "Low")
  
  ha_col <- HeatmapAnnotation(
    DGAT1 = dgat_groups,
    col = list(DGAT1 = c("High" = "#D73027", "Low" = "#4575B4"))
  )
  
  ht <- Heatmap(
    marker_z,
    name = "Z-score",
    col = colorRamp2(c(-2, 0, 2), c("#4575B4", "white", "#D73027")),
    cluster_rows = TRUE,
    cluster_columns = TRUE,
    left_annotation = ha_row,
    top_annotation = ha_col,
    show_row_names = TRUE,
    show_column_names = FALSE,
    row_names_gp = gpar(fontsize = 9),
    column_title = "Individual Immune Marker Expression",
    column_title_gp = gpar(fontsize = 14, fontface = "bold")
  )
  
  png(outfile, width = 12, height = max(6, nrow(marker_z) * 0.3), 
      units = "in", res = 300)
  draw(ht)
  dev.off()
  
  return(ht)
}

# =============================================================================
# MAIN
# =============================================================================

main <- function() {
  cat("\n", paste(rep("=", 70), collapse = ""), "\n")
  cat("CPTAC GBM: Individual Immune Marker Analysis (PRODUCTION)\n")
  cat(paste(rep("=", 70), collapse = ""), "\n\n")
  
  # Load data with robust handling
  prot_mat <- load_protein_data(PROTEIN_FILE)
  
  if (!"DGAT1" %in% rownames(prot_mat)) {
    stop("DGAT1 not found in proteomics data")
  }
  
  # FIX 3: Check DGAT1 distribution and winsorize
  dgat <- prot_mat["DGAT1", ]
  dgat <- check_dgat1_distribution(dgat, winsorize = TRUE, 
                                     percentile = WINSORIZE_PERCENTILE)
  
  # FIX 2 & 4: Filter markers by detection + save coverage reports
  filter_results <- filter_markers_by_detection(
    prot_mat, MARKER_DEFINITIONS, MIN_DETECTION_RATE
  )
  
  write_csv(filter_results$present_df, 
            file.path(OUT_DIR, "markers_present_passed.csv"))
  write_csv(filter_results$failed_df, 
            file.path(OUT_DIR, "markers_present_failed_detection.csv"))
  write_csv(filter_results$missing_df, 
            file.path(OUT_DIR, "markers_missing_from_data.csv"))
  
  # Analyze correlations
  cat("Running correlation analysis...\n")
  cor_results <- analyze_individual_markers(
    prot_mat, dgat, filter_results$passed, MARKER_DEFINITIONS,
    use_partial = USE_PARTIAL_CORRELATION
  )
  
  write_csv(cor_results, file.path(OUT_DIR, "individual_marker_correlations.csv"))
  
  # Summary
  cat("\n=== RESULTS SUMMARY ===\n")
  cat("Markers analyzed:", nrow(cor_results), "\n")
  cat("Significant (FDR < 0.05):", sum(cor_results$FDR < 0.05), "\n")
  cat("Nominal (P < 0.05):", sum(cor_results$P < 0.05), "\n\n")
  
  summary_table <- cor_results %>%
    group_by(CellType) %>%
    summarise(
      N_markers = n(),
      Mean_rho = mean(Rho),
      N_positive = sum(Rho > 0),
      N_negative = sum(Rho < 0),
      N_sig_FDR = sum(FDR < 0.05),
      .groups = "drop"
    )
  print(summary_table)
  
  cat("\nTop 10 correlations:\n")
  print(head(cor_results %>% select(Protein, CellType, Rho, P, FDR), 10))
  
  # Generate plots
  cat("\nGenerating visualizations...\n")
  
  plot_volcano_colored(cor_results,
                       file.path(OUT_DIR, "volcano_individual_markers.png"))
  cat("  ✓ Volcano plot\n")
  
  plot_dotplot_annotated(cor_results,
                         file.path(OUT_DIR, "dotplot_individual_markers.png"))
  cat("  ✓ Dot plot\n")
  
  plot_clustered_heatmap(prot_mat, dgat, cor_results,
                         file.path(OUT_DIR, "heatmap_clustered_markers.png"))
  cat("  ✓ Heatmap\n")
  
  cat("\n", paste(rep("=", 70), collapse = ""), "\n")
  cat("Analysis complete. Results in:", OUT_DIR, "\n")
  cat(paste(rep("=", 70), collapse = ""), "\n\n")
  
  invisible(cor_results)
}

# Run
if (!interactive()) {
  results <- main()
}
