#!/usr/bin/env Rscript
# =============================================================================
# TCGA GBM: IMMUNE MARKERS Analysis Only
# Dual plotting: Cell Type vs. Immune Function
# RNA-seq expression data (TPM/FPKM)
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
DATA_DIR <- file.path(BASE_DIR, "Processed_Data", "TCGA_GBM_Clean")
OUT_DIR  <- file.path(BASE_DIR, "Results", "TCGA_Immune_Markers")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

EXPRESSION_FILE <- file.path(DATA_DIR, "TCGA_GBM_Expression_Cleaned.rds")

MIN_DETECTION_RATE <- 0.30
WINSORIZE_PERCENTILE <- 0.01

# =============================================================================
# IMMUNE MARKERS ONLY (No Lipid Metabolism)
# =============================================================================

IMMUNE_MARKERS <- tribble(
  ~Gene,        ~CellType,            ~ImmuneFunction,        ~Description,
  
  # Macrophages & Microglia
  "CD68",       "Macrophage",         "Pan-Myeloid",          "Pan-macrophage",
  "CD163",      "Macrophage",         "Immunosuppressive",    "M2 macrophage",
  "MRC1",       "Macrophage",         "Immunosuppressive",    "CD206 (M2)",
  "APOE",       "Macrophage",         "Immunosuppressive",    "Lipid-loaded TAMs",
  "ARG1",       "TAM",                "Immunosuppressive",    "Arginase-1",
  "IL10",       "TAM",                "Immunosuppressive",    "IL-10",
  "TGFB1",      "TAM",                "Immunosuppressive",    "TGF-beta",
  "SPP1",       "TAM",                "Immunosuppressive",    "Osteopontin",
  "VEGFA",      "TAM",                "Immunosuppressive",    "VEGF-A",
  "NOS2",       "M1-Macrophage",      "Pro-inflammatory",     "iNOS (M1)",
  "IL12A",      "M1-Macrophage",      "Anti-tumor",           "IL-12p35",
  "CXCL10",     "M1-Macrophage",      "Anti-tumor",           "IP-10",
  "P2RY12",     "Microglia",          "Pan-Myeloid",          "Microglia-specific",
  "TMEM119",    "Microglia",          "Pan-Myeloid",          "Microglia-specific",
  
  # MDSCs
  "S100A8",     "MDSC",               "Immunosuppressive",    "M-MDSC",
  "S100A9",     "MDSC",               "Immunosuppressive",    "M-MDSC",
  "CD14",       "MDSC",               "Immunosuppressive",    "Monocytic MDSC",
  "ARG2",       "MDSC",               "Immunosuppressive",    "Arginase-2",
  "IDO1",       "MDSC",               "Immunosuppressive",    "IDO1",
  
  # T Cells - Effector
  "CD8A",       "CD8-T",              "Anti-tumor",           "CD8+ T cell",
  "CD8B",       "CD8-T",              "Anti-tumor",           "CD8+ T cell",
  "CD4",        "CD4-T",              "Anti-tumor",           "CD4+ T helper",
  "GZMA",       "Cytotoxic",          "Anti-tumor",           "Granzyme A",
  "GZMB",       "Cytotoxic",          "Anti-tumor",           "Granzyme B",
  "PRF1",       "Cytotoxic",          "Anti-tumor",           "Perforin",
  "IFNG",       "Th1",                "Anti-tumor",           "IFN-gamma",
  "CD3E",       "Pan-T",              "Anti-tumor",           "TCR component",
  
  # T Cells - Exhausted
  "PDCD1",      "Exhausted-T",        "Dysfunctional",        "PD-1",
  "LAG3",       "Exhausted-T",        "Dysfunctional",        "LAG3",
  "HAVCR2",     "Exhausted-T",        "Dysfunctional",        "TIM-3",
  "TIGIT",      "Exhausted-T",        "Dysfunctional",        "TIGIT",
  "CTLA4",      "Exhausted-T",        "Dysfunctional",        "CTLA-4",
  "TOX",        "Exhausted-T",        "Dysfunctional",        "TOX TF",
  
  # Regulatory T Cells
  "FOXP3",      "Treg",               "Immunosuppressive",    "Treg master TF",
  "IL2RA",      "Treg",               "Immunosuppressive",    "CD25",
  
  # NK Cells
  "NCAM1",      "NK-cell",            "Anti-tumor",           "CD56",
  "KLRD1",      "NK-cell",            "Anti-tumor",           "CD94",
  "KLRK1",      "NK-cell",            "Anti-tumor",           "NKG2D",
  "NCR1",       "NK-cell",            "Anti-tumor",           "NKp46",
  
  # Dendritic Cells
  "CLEC9A",     "DC",                 "Anti-tumor",           "cDC1",
  "BATF3",      "DC",                 "Anti-tumor",           "cDC1 TF",
  "CD1C",       "DC",                 "Anti-tumor",           "cDC2",
  "FCER1A",     "DC",                 "Anti-tumor",           "cDC2",
  
  # B Cells
  "CD19",       "B-cell",             "Adaptive-immune",      "B cell",
  "MS4A1",      "B-cell",             "Adaptive-immune",      "CD20",
  
  # Neutrophils
  "FCGR3B",     "Neutrophil",         "Context-dependent",    "CD16b",
  "CSF3R",      "Neutrophil",         "Context-dependent",    "G-CSF receptor",
  
  # Checkpoints
  "CD274",      "Checkpoint",         "Immunosuppressive",    "PD-L1",
  "PDCD1LG2",   "Checkpoint",         "Immunosuppressive",    "PD-L2",
  "CD276",      "Checkpoint",         "Immunosuppressive",    "B7-H3",
  
  # Cytokines & Chemokines
  "IL6",        "Cytokine",           "Pro-inflammatory",     "IL-6",
  "IL12B",      "Cytokine",           "Anti-tumor",           "IL-12p40",
  "IL15",       "Cytokine",           "Anti-tumor",           "IL-15",
  "TNF",        "Cytokine",           "Pro-inflammatory",     "TNF-alpha",
  "IL4",        "Cytokine",           "Immunosuppressive",    "IL-4",
  "CCL5",       "Chemokine",          "Anti-tumor",           "RANTES",
  "CCL2",       "Chemokine",          "Immunosuppressive",    "MCP-1",
  
  # Antigen Presentation
  "HLA-A",      "MHC-I",              "Anti-tumor",           "MHC class I",
  "HLA-DRA",    "MHC-II",             "Anti-tumor",           "MHC class II",
  "B2M",        "MHC-I",              "Anti-tumor",           "Beta-2-microglobulin",
  "TAP1",       "Antigen-Process",    "Anti-tumor",           "Peptide transporter"
)

# =============================================================================
# COLOR PALETTES
# =============================================================================

CELLTYPE_COLORS <- c(
  "Macrophage" = "#E41A1C", "M1-Macrophage" = "#FF7F00", "TAM" = "#984EA3",
  "Microglia" = "#A65628", "MDSC" = "#8B4513", "CD8-T" = "#4DAF4A",
  "CD4-T" = "#377EB8", "Pan-T" = "#4DAF4A", "Cytotoxic" = "#006400",
  "Th1" = "#228B22", "Exhausted-T" = "#B0B0B0", "Treg" = "#9370DB",
  "NK-cell" = "#FF6347", "DC" = "#FFD700", "B-cell" = "#00CED1",
  "Neutrophil" = "#F781BF", "Checkpoint" = "#696969", "Cytokine" = "#FF69B4",
  "Chemokine" = "#DB7093", "MHC-I" = "#4682B4", "MHC-II" = "#5F9EA0", 
  "Antigen-Process" = "#6495ED"
)

FUNCTION_COLORS <- c(
  "Anti-tumor" = "#228B22", "Immunosuppressive" = "#B22222",
  "Dysfunctional" = "#808080", "Pro-inflammatory" = "#FF8C00",
  "Pan-Myeloid" = "#D3D3D3", "Adaptive-immune" = "#4169E1",
  "Context-dependent" = "#DAA520"
)

# =============================================================================
# FUNCTIONS
# =============================================================================

load_expression_matrix <- function(filepath) {
  cat("Loading RNA-seq data...\n")
  if (grepl("\\.rds$", filepath, ignore.case = TRUE)) {
    expr <- readRDS(filepath)
  } else if (grepl("\\.csv$", filepath, ignore.case = TRUE)) {
    expr <- read_csv(filepath, show_col_types = FALSE)
  } else {
    expr <- read_tsv(filepath, show_col_types = FALSE)
  }
  
  if (is.data.frame(expr)) {
    gene_cols <- c("Gene", "gene", "Gene_Symbol", "SYMBOL")
    gene_col <- intersect(names(expr), gene_cols)[1]
    if (is.na(gene_col)) stop("Could not identify gene column")
    sample_cols <- setdiff(names(expr), c(gene_col, "Gene_ID", "Ensembl_ID"))
    expr_mat <- as.matrix(expr[, sample_cols])
    rownames(expr_mat) <- expr[[gene_col]]
    expr <- expr_mat
  }
  
  if (any(duplicated(rownames(expr)))) {
    cat("  Collapsing", sum(duplicated(rownames(expr))), "duplicates by median\n")
    expr_df <- as.data.frame(expr) %>%
      mutate(Gene = rownames(expr)) %>%
      group_by(Gene) %>%
      summarise(across(everything(), ~median(.x, na.rm = TRUE)), .groups = "drop")
    expr <- as.matrix(expr_df[, -1])
    rownames(expr) <- expr_df$Gene
  }
  
  cat("  Matrix:", nrow(expr), "genes x", ncol(expr), "samples\n\n")
  return(expr)
}

check_dgat1_distribution <- function(dgat_vec, winsorize = TRUE, percentile = 0.01) {
  cat("\n=== DGAT1 Distribution ===\n")
  cat("Samples:", sum(!is.na(dgat_vec)), "/", length(dgat_vec), "\n")
  cat("Range: [", round(min(dgat_vec, na.rm = TRUE), 3), ",", 
      round(max(dgat_vec, na.rm = TRUE), 3), "]\n")
  cat("Mean ± SD:", round(mean(dgat_vec, na.rm = TRUE), 3), "±", 
      round(sd(dgat_vec, na.rm = TRUE), 3), "\n")
  cat("Median [IQR]:", round(median(dgat_vec, na.rm = TRUE), 3), "[",
      round(quantile(dgat_vec, 0.25, na.rm = TRUE), 3), "-",
      round(quantile(dgat_vec, 0.75, na.rm = TRUE), 3), "]\n")
  
  if (winsorize && percentile > 0) {
    lower <- quantile(dgat_vec, percentile, na.rm = TRUE)
    upper <- quantile(dgat_vec, 1 - percentile, na.rm = TRUE)
    n_winsorized <- sum(dgat_vec < lower | dgat_vec > upper, na.rm = TRUE)
    if (n_winsorized > 0) {
      cat("\nWinsorizing", n_winsorized, "outliers\n")
      dgat_vec[dgat_vec < lower] <- lower
      dgat_vec[dgat_vec > upper] <- upper
    }
  }
  cat("================================\n\n")
  return(dgat_vec)
}

filter_markers_by_detection <- function(expr_mat, markers_df, min_rate = 0.30) {
  cat("\n=== Filtering Markers ===\n")
  present <- markers_df$Gene[markers_df$Gene %in% rownames(expr_mat)]
  missing <- markers_df$Gene[!markers_df$Gene %in% rownames(expr_mat)]
  if (length(present) == 0) stop("No markers found")
  
  detection_rates <- sapply(present, function(g) {
    sum(!is.na(expr_mat[g, ])) / ncol(expr_mat)
  })
  
  passed <- present[detection_rates >= min_rate]
  failed <- present[detection_rates < min_rate]
  
  cat("  Total:", nrow(markers_df), "| Found:", length(present), 
      "| Passed:", length(passed), "| Failed:", length(failed), 
      "| Missing:", length(missing), "\n\n")
  
  list(
    passed = passed,
    present_df = markers_df %>% filter(Gene %in% passed) %>% 
                  mutate(DetectionRate = detection_rates[Gene]),
    failed_df = markers_df %>% filter(Gene %in% failed) %>% 
                  mutate(DetectionRate = detection_rates[Gene]),
    missing_df = markers_df %>% filter(Gene %in% missing) %>% 
                  mutate(DetectionRate = NA)
  )
}

analyze_markers <- function(expr_mat, dgat_vec, passed, markers_df) {
  cat("Analyzing", length(passed), "markers...\n")
  results <- lapply(passed, function(gene) {
    vec <- expr_mat[gene, ]
    valid <- !is.na(vec) & !is.na(dgat_vec)
    if (sum(valid) < 10) return(NULL)
    
    test <- cor.test(vec[valid], dgat_vec[valid], method = "spearman")
    marker_info <- markers_df[markers_df$Gene == gene, ]
    
    data.frame(
      Gene = gene,
      CellType = marker_info$CellType,
      ImmuneFunction = marker_info$ImmuneFunction,
      Description = marker_info$Description,
      Rho = unname(test$estimate),
      P = test$p.value,
      N = sum(valid),
      stringsAsFactors = FALSE
    )
  })
  
  tab <- do.call(rbind, results[!sapply(results, is.null)])
  tab$FDR <- p.adjust(tab$P, method = "BH")
  tab$Significant <- tab$FDR < 0.05
  tab[order(tab$FDR), ]
}

plot_volcano_celltype <- function(cor_results, outfile) {
  if (is.null(cor_results) || nrow(cor_results) == 0) return(NULL)
  cor_results$CellType <- factor(cor_results$CellType, levels = names(CELLTYPE_COLORS))
  cor_results$Label <- ifelse(cor_results$FDR < 0.1 | abs(cor_results$Rho) > 0.2, 
                               cor_results$Gene, "")
  
  p <- ggplot(cor_results, aes(x = Rho, y = -log10(P))) +
    geom_point(aes(color = CellType, size = -log10(P), 
                   alpha = ifelse(FDR < 0.05, 1, 0.5))) +
    geom_text_repel(aes(label = Label), size = 3, max.overlaps = 25) +
    scale_color_manual(values = CELLTYPE_COLORS, name = "Cell Type") +
    scale_size_continuous(range = c(2, 6), guide = "none") +
    scale_alpha_identity() +
    geom_vline(xintercept = 0, color = "black") +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "blue") +
    labs(title = "DGAT1 vs Immune Markers (TCGA GBM)",
         subtitle = "Colored by Cell Type",
         x = "Spearman ρ", y = "-log10(P-value)") +
    theme_minimal(base_size = 12) +
    theme(legend.position = "right", plot.title = element_text(face = "bold"))
  
  ggsave(outfile, p, width = 14, height = 8, dpi = 300)
  cat("  ✓ Cell type volcano\n")
  return(p)
}

plot_volcano_function <- function(cor_results, outfile) {
  if (is.null(cor_results) || nrow(cor_results) == 0) return(NULL)
  cor_results$ImmuneFunction <- factor(cor_results$ImmuneFunction, levels = names(FUNCTION_COLORS))
  cor_results$Label <- ifelse(cor_results$FDR < 0.1 | abs(cor_results$Rho) > 0.2, 
                               cor_results$Gene, "")
  
  p <- ggplot(cor_results, aes(x = Rho, y = -log10(P))) +
    geom_point(aes(color = ImmuneFunction, size = -log10(P), 
                   alpha = ifelse(FDR < 0.05, 1, 0.5))) +
    geom_text_repel(aes(label = Label), size = 3, max.overlaps = 25) +
    scale_color_manual(values = FUNCTION_COLORS, name = "Immune Function") +
    scale_size_continuous(range = c(2, 6), guide = "none") +
    scale_alpha_identity() +
    geom_vline(xintercept = 0, color = "black") +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "blue") +
    labs(title = "DGAT1 vs Immune Markers (TCGA GBM)",
         subtitle = "Colored by Immune Function (Pro- vs Anti-tumor)",
         x = "Spearman ρ", y = "-log10(P-value)") +
    theme_minimal(base_size = 12) +
    theme(legend.position = "right", plot.title = element_text(face = "bold"))
  
  ggsave(outfile, p, width = 14, height = 8, dpi = 300)
  cat("  ✓ Function volcano\n")
  return(p)
}

plot_dotplot_celltype <- function(cor_results, outfile) {
  cor_results <- cor_results %>% arrange(Rho) %>% mutate(Gene = factor(Gene, levels = Gene))
  cor_results$SigLabel <- case_when(
    cor_results$FDR < 0.001 ~ "***", cor_results$FDR < 0.01 ~ "**",
    cor_results$FDR < 0.05 ~ "*", cor_results$P < 0.05 ~ "†", TRUE ~ ""
  )
  
  p <- ggplot(cor_results, aes(y = Gene, x = Rho)) +
    geom_segment(aes(yend = Gene, xend = 0, color = CellType), linewidth = 0.8) +
    geom_point(aes(color = CellType, size = -log10(P))) +
    geom_text(aes(label = SigLabel), hjust = -0.5, size = 5, fontface = "bold") +
    geom_vline(xintercept = 0, color = "black") +
    scale_color_manual(values = CELLTYPE_COLORS) +
    scale_size_continuous(range = c(2, 6), name = "-log10(P)") +
    labs(title = "DGAT1 Correlation (TCGA)", subtitle = "By Cell Type",
         x = "Spearman ρ", y = "") +
    theme_minimal(base_size = 11) +
    theme(axis.text.y = element_text(size = 9), legend.position = "right")
  
  ggsave(outfile, p, width = 10, height = max(6, nrow(cor_results) * 0.25), 
         dpi = 300, limitsize = FALSE)
  cat("  ✓ Cell type dot plot\n")
  return(p)
}

plot_dotplot_function <- function(cor_results, outfile) {
  cor_results <- cor_results %>% arrange(Rho) %>% mutate(Gene = factor(Gene, levels = Gene))
  cor_results$SigLabel <- case_when(
    cor_results$FDR < 0.001 ~ "***", cor_results$FDR < 0.01 ~ "**",
    cor_results$FDR < 0.05 ~ "*", cor_results$P < 0.05 ~ "†", TRUE ~ ""
  )
  
  p <- ggplot(cor_results, aes(y = Gene, x = Rho)) +
    geom_segment(aes(yend = Gene, xend = 0, color = ImmuneFunction), linewidth = 0.8) +
    geom_point(aes(color = ImmuneFunction, size = -log10(P))) +
    geom_text(aes(label = SigLabel), hjust = -0.5, size = 5, fontface = "bold") +
    geom_vline(xintercept = 0, color = "black") +
    scale_color_manual(values = FUNCTION_COLORS) +
    scale_size_continuous(range = c(2, 6), name = "-log10(P)") +
    labs(title = "DGAT1 Correlation (TCGA)", subtitle = "By Function",
         x = "Spearman ρ", y = "") +
    theme_minimal(base_size = 11) +
    theme(axis.text.y = element_text(size = 9), legend.position = "right")
  
  ggsave(outfile, p, width = 10, height = max(6, nrow(cor_results) * 0.25), 
         dpi = 300, limitsize = FALSE)
  cat("  ✓ Function dot plot\n")
  return(p)
}

plot_heatmap <- function(expr_mat, dgat_vec, cor_results, outfile) {
  if (!COMPLEX_HEATMAP_AVAILABLE) {
    cat("  ⚠ ComplexHeatmap not available\n")
    return(NULL)
  }
  
  present <- cor_results$Gene[cor_results$Gene %in% rownames(expr_mat)]
  if (length(present) < 2) return(NULL)
  
  marker_mat <- expr_mat[present, ]
  if (max(marker_mat, na.rm = TRUE) > 100) marker_mat <- log2(marker_mat + 1)
  marker_z <- t(scale(t(marker_mat)))
  valid_samples <- colSums(!is.na(marker_z)) > 0
  marker_z <- marker_z[, valid_samples]
  dgat_aligned <- dgat_vec[colnames(marker_z)]
  
  row_anno <- cor_results %>% filter(Gene %in% rownames(marker_z)) %>%
    arrange(match(Gene, rownames(marker_z)))
  
  ha_row <- rowAnnotation(
    CellType = row_anno$CellType,
    Function = row_anno$ImmuneFunction,
    Correlation = anno_barplot(row_anno$Rho, border = FALSE,
      gp = gpar(fill = ifelse(row_anno$Rho > 0, "#D73027", "#4575B4"))),
    col = list(CellType = CELLTYPE_COLORS, Function = FUNCTION_COLORS)
  )
  
  dgat_groups <- ifelse(dgat_aligned >= median(dgat_aligned, na.rm = TRUE), "High", "Low")
  ha_col <- HeatmapAnnotation(
    DGAT1 = dgat_groups,
    col = list(DGAT1 = c("High" = "#D73027", "Low" = "#4575B4"))
  )
  
  ht <- Heatmap(
    marker_z, name = "Z-score",
    col = colorRamp2(c(-2, 0, 2), c("#4575B4", "white", "#D73027")),
    cluster_rows = TRUE, cluster_columns = TRUE,
    left_annotation = ha_row, top_annotation = ha_col,
    show_row_names = TRUE, show_column_names = FALSE,
    row_names_gp = gpar(fontsize = 9),
    column_title = "Immune Marker Expression (TCGA GBM)"
  )
  
  png(outfile, width = 14, height = max(8, nrow(marker_z) * 0.3), units = "in", res = 300)
  draw(ht)
  dev.off()
  cat("  ✓ Heatmap\n")
  return(ht)
}

# =============================================================================
# MAIN
# =============================================================================

main <- function() {
  cat("\n", paste(rep("=", 70), collapse = ""), "\n")
  cat("TCGA GBM: Immune Marker Analysis\n")
  cat(paste(rep("=", 70), collapse = ""), "\n\n")
  
  expr_mat <- load_expression_matrix(EXPRESSION_FILE)
  if (!"DGAT1" %in% rownames(expr_mat)) stop("DGAT1 not found")
  
  dgat <- expr_mat["DGAT1", ]
  dgat <- check_dgat1_distribution(dgat, winsorize = TRUE, percentile = WINSORIZE_PERCENTILE)
  
  filter_res <- filter_markers_by_detection(expr_mat, IMMUNE_MARKERS, MIN_DETECTION_RATE)
  write_csv(filter_res$present_df, file.path(OUT_DIR, "markers_present.csv"))
  write_csv(filter_res$failed_df, file.path(OUT_DIR, "markers_failed.csv"))
  write_csv(filter_res$missing_df, file.path(OUT_DIR, "markers_missing.csv"))
  
  cor_results <- analyze_markers(expr_mat, dgat, filter_res$passed, IMMUNE_MARKERS)
  write_csv(cor_results, file.path(OUT_DIR, "immune_marker_correlations.csv"))
  
  cat("\nMarkers:", nrow(cor_results), "| Sig (FDR<0.05):", sum(cor_results$FDR < 0.05), 
      "| Nominal (P<0.05):", sum(cor_results$P < 0.05), "\n\n")
  
  cat("By Cell Type:\n")
  print(cor_results %>% group_by(CellType) %>%
    summarise(N = n(), Mean_rho = mean(Rho), N_sig = sum(FDR < 0.05), .groups = "drop"))
  
  cat("\nBy Immune Function:\n")
  print(cor_results %>% group_by(ImmuneFunction) %>%
    summarise(N = n(), Mean_rho = mean(Rho), N_sig = sum(FDR < 0.05), .groups = "drop"))
  
  cat("\nTop 10:\n")
  print(head(cor_results %>% select(Gene, CellType, ImmuneFunction, Rho, P, FDR), 10))
  
  cat("\n=== Generating Plots ===\n")
  plot_volcano_celltype(cor_results, file.path(OUT_DIR, "volcano_by_celltype.png"))
  plot_volcano_function(cor_results, file.path(OUT_DIR, "volcano_by_function.png"))
  plot_dotplot_celltype(cor_results, file.path(OUT_DIR, "dotplot_by_celltype.png"))
  plot_dotplot_function(cor_results, file.path(OUT_DIR, "dotplot_by_function.png"))
  plot_heatmap(expr_mat, dgat, cor_results, file.path(OUT_DIR, "heatmap_markers.png"))
  
  cat("\n", paste(rep("=", 70), collapse = ""), "\n")
  cat("Complete! Results in:", OUT_DIR, "\n")
  cat(paste(rep("=", 70), collapse = ""), "\n\n")
  
  invisible(cor_results)
}

if (!interactive()) results <- main()
