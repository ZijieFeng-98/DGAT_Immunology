#!/usr/bin/env Rscript
# =============================================================================
# CPTAC GBM: Streamlined Individual Immune & Metabolic Marker Analysis
# Dual plotting strategy: Cell Type vs. Immune Function
# Publication-ready for DGAT1 cancer immunology
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
OUT_DIR  <- file.path(BASE_DIR, "Results", "CPTAC_Streamlined_Dual_Analysis")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

PROTEIN_FILE <- file.path(DATA_DIR, "CPTAC3_Glioblastoma_Multiforme_Proteome.tmt11.tsv")

# QC thresholds
MIN_DETECTION_RATE <- 0.30  # Require 30% non-NA values
WINSORIZE_PERCENTILE <- 0.01  # Trim 1% tails for outliers

# =============================================================================
# STREAMLINED MARKER DEFINITIONS (74 markers total)
# Same as RNA analysis for cross-platform consistency
# =============================================================================

CURATED_MARKERS <- tribble(
  ~Gene,        ~CellType,            ~ImmuneFunction,        ~Description,
  
  # =========================================================================
  # MACROPHAGES & MICROGLIA
  # =========================================================================
  "CD68",       "Macrophage",         "Pan-Myeloid",          "Pan-macrophage/microglia marker",
  "CD163",      "Macrophage",         "Immunosuppressive",    "M2 macrophage, scavenger receptor",
  "MRC1",       "Macrophage",         "Immunosuppressive",    "CD206, mannose receptor (M2)",
  "APOE",       "Macrophage",         "Immunosuppressive",    "Apolipoprotein E, lipid-loaded TAMs",
  
  # M2/TAM markers
  "ARG1",       "TAM",                "Immunosuppressive",    "Arginase-1, depletes arginine",
  "IL10",       "TAM",                "Immunosuppressive",    "Anti-inflammatory cytokine",
  "TGFB1",      "TAM",                "Immunosuppressive",    "TGF-beta, suppresses T cells",
  "SPP1",       "TAM",                "Immunosuppressive",    "Osteopontin, TAM marker",
  "VEGFA",      "TAM",                "Immunosuppressive",    "Angiogenic factor",
  
  # M1 markers
  "NOS2",       "M1-Macrophage",      "Pro-inflammatory",     "iNOS, M1 marker",
  "IL12A",      "M1-Macrophage",      "Anti-tumor",           "IL-12p35",
  "CXCL10",     "M1-Macrophage",      "Anti-tumor",           "IP-10, T cell recruitment",
  
  # Microglia-specific
  "P2RY12",     "Microglia",          "Pan-Myeloid",          "Microglia-specific",
  "TMEM119",    "Microglia",          "Pan-Myeloid",          "Microglia-specific",
  
  # =========================================================================
  # MDSCs
  # =========================================================================
  "S100A8",     "MDSC",               "Immunosuppressive",    "M-MDSC marker",
  "S100A9",     "MDSC",               "Immunosuppressive",    "M-MDSC marker",
  "CD14",       "MDSC",               "Immunosuppressive",    "Monocytic MDSC",
  "ARG2",       "MDSC",               "Immunosuppressive",    "Arginase-2",
  "IDO1",       "MDSC",               "Immunosuppressive",    "Tryptophan catabolism",
  
  # =========================================================================
  # T CELLS - EFFECTOR
  # =========================================================================
  "CD8A",       "CD8-T",              "Anti-tumor",           "CD8+ T cell",
  "CD8B",       "CD8-T",              "Anti-tumor",           "CD8+ T cell",
  "CD4",        "CD4-T",              "Anti-tumor",           "CD4+ T helper",
  "GZMA",       "Cytotoxic",          "Anti-tumor",           "Granzyme A",
  "GZMB",       "Cytotoxic",          "Anti-tumor",           "Granzyme B",
  "PRF1",       "Cytotoxic",          "Anti-tumor",           "Perforin",
  "IFNG",       "Th1",                "Anti-tumor",           "IFN-gamma",
  "CD3E",       "Pan-T",              "Anti-tumor",           "T cell receptor",
  
  # =========================================================================
  # T CELLS - EXHAUSTED
  # =========================================================================
  "PDCD1",      "Exhausted-T",        "Dysfunctional",        "PD-1",
  "LAG3",       "Exhausted-T",        "Dysfunctional",        "LAG3",
  "HAVCR2",     "Exhausted-T",        "Dysfunctional",        "TIM-3",
  "TIGIT",      "Exhausted-T",        "Dysfunctional",        "TIGIT",
  "CTLA4",      "Exhausted-T",        "Dysfunctional",        "CTLA-4",
  "TOX",        "Exhausted-T",        "Dysfunctional",        "TOX TF",
  
  # =========================================================================
  # REGULATORY T CELLS
  # =========================================================================
  "FOXP3",      "Treg",               "Immunosuppressive",    "Treg master TF",
  "IL2RA",      "Treg",               "Immunosuppressive",    "CD25",
  
  # =========================================================================
  # NK CELLS
  # =========================================================================
  "NCAM1",      "NK-cell",            "Anti-tumor",           "CD56",
  "KLRD1",      "NK-cell",            "Anti-tumor",           "CD94",
  "KLRK1",      "NK-cell",            "Anti-tumor",           "NKG2D",
  "NCR1",       "NK-cell",            "Anti-tumor",           "NKp46",
  
  # =========================================================================
  # DENDRITIC CELLS
  # =========================================================================
  "CLEC9A",     "DC",                 "Anti-tumor",           "cDC1",
  "BATF3",      "DC",                 "Anti-tumor",           "cDC1 TF",
  "CD1C",       "DC",                 "Anti-tumor",           "cDC2",
  "FCER1A",     "DC",                 "Anti-tumor",           "cDC2",
  
  # =========================================================================
  # B CELLS
  # =========================================================================
  "CD19",       "B-cell",             "Adaptive-immune",      "B cell marker",
  "MS4A1",      "B-cell",             "Adaptive-immune",      "CD20",
  
  # =========================================================================
  # NEUTROPHILS
  # =========================================================================
  "FCGR3B",     "Neutrophil",         "Context-dependent",    "CD16b",
  "CSF3R",      "Neutrophil",         "Context-dependent",    "G-CSF receptor",
  
  # =========================================================================
  # IMMUNE CHECKPOINTS
  # =========================================================================
  "CD274",      "Checkpoint",         "Immunosuppressive",    "PD-L1",
  "PDCD1LG2",   "Checkpoint",         "Immunosuppressive",    "PD-L2",
  "CD276",      "Checkpoint",         "Immunosuppressive",    "B7-H3",
  
  # =========================================================================
  # CYTOKINES & CHEMOKINES
  # =========================================================================
  "IL6",        "Cytokine",           "Pro-inflammatory",     "IL-6",
  "IL12B",      "Cytokine",           "Anti-tumor",           "IL-12p40",
  "IL15",       "Cytokine",           "Anti-tumor",           "IL-15",
  "TNF",        "Cytokine",           "Pro-inflammatory",     "TNF-alpha",
  "IL4",        "Cytokine",           "Immunosuppressive",    "IL-4",
  "CCL5",       "Chemokine",          "Anti-tumor",           "RANTES",
  "CCL2",       "Chemokine",          "Immunosuppressive",    "MCP-1",
  
  # =========================================================================
  # LIPID METABOLISM (7 KEY GENES)
  # =========================================================================
  "DGAT1",      "Lipid-LD",           "Metabolic",            "Lipid droplet formation",
  "SREBF1",     "Lipid-Synthesis",    "Metabolic",            "SREBP-1",
  "SCAP",       "Lipid-Regulation",   "Metabolic",            "SREBP escort",
  "ACACA",      "Lipid-Synthesis",    "Metabolic",            "ACC",
  "FASN",       "Lipid-Synthesis",    "Metabolic",            "Fatty acid synthase",
  "SCD",        "Lipid-Synthesis",    "Metabolic",            "SCD1",
  "HMGCR",      "Lipid-Synthesis",    "Metabolic",            "HMG-CoA reductase",
  
  # =========================================================================
  # HYPOXIA
  # =========================================================================
  "HIF1A",      "Hypoxia",            "Metabolic",            "HIF-1alpha",
  
  # =========================================================================
  # ANTIGEN PRESENTATION
  # =========================================================================
  "HLA-A",      "MHC-I",              "Anti-tumor",           "MHC class I",
  "HLA-DRA",    "MHC-II",             "Anti-tumor",           "MHC class II",
  "B2M",        "MHC-I",              "Anti-tumor",           "Beta-2-microglobulin",
  "TAP1",       "Antigen-Process",    "Anti-tumor",           "Peptide transporter"
)

# =============================================================================
# COLOR PALETTES
# =============================================================================

# Cell Type Colors
CELLTYPE_COLORS <- c(
  "Macrophage" = "#E41A1C",
  "M1-Macrophage" = "#FF7F00",
  "TAM" = "#984EA3",
  "Microglia" = "#A65628",
  "MDSC" = "#8B4513",
  "CD8-T" = "#4DAF4A",
  "CD4-T" = "#377EB8",
  "Pan-T" = "#4DAF4A",
  "Cytotoxic" = "#006400",
  "Th1" = "#228B22",
  "Exhausted-T" = "#B0B0B0",
  "Treg" = "#9370DB",
  "NK-cell" = "#FF6347",
  "DC" = "#FFD700",
  "B-cell" = "#00CED1",
  "Neutrophil" = "#F781BF",
  "Checkpoint" = "#696969",
  "Cytokine" = "#FF69B4",
  "Chemokine" = "#DB7093",
  "Lipid-LD" = "#8B0000",
  "Lipid-Synthesis" = "#DC143C",
  "Lipid-Regulation" = "#CD5C5C",
  "Hypoxia" = "#2F4F4F",
  "MHC-I" = "#4682B4",
  "MHC-II" = "#5F9EA0",
  "Antigen-Process" = "#6495ED"
)

# Immune Function Colors
FUNCTION_COLORS <- c(
  "Anti-tumor" = "#228B22",
  "Immunosuppressive" = "#B22222",
  "Dysfunctional" = "#808080",
  "Pro-inflammatory" = "#FF8C00",
  "Pan-Myeloid" = "#D3D3D3",
  "Adaptive-immune" = "#4169E1",
  "Context-dependent" = "#DAA520",
  "Metabolic" = "#8B4513"
)

# =============================================================================
# FUNCTIONS
# =============================================================================

load_protein_data <- function(filepath) {
  cat("Loading proteomics data...\n")
  prot <- read_tsv(filepath, show_col_types = FALSE)
  
  gene_cols <- c("Gene", "gene", "Gene_Symbol")
  gene_col <- gene_cols[gene_cols %in% names(prot)][1]
  
  sample_cols <- setdiff(names(prot), c(gene_col, "Gene_ID", "Protein_ID"))
  
  prot_mat <- as.matrix(prot[, sample_cols])
  rownames(prot_mat) <- prot[[gene_col]]
  
  # Convert to numeric, handling any non-numeric values
  prot_mat <- apply(prot_mat, 2, function(x) {
    as.numeric(as.character(x))
  })
  rownames(prot_mat) <- prot[[gene_col]]
  
  # Remove non-gene rows (like "Mean", "Median", "StdDev")
  non_gene_rows <- c("Mean", "Median", "StdDev", "Gene")
  prot_mat <- prot_mat[!rownames(prot_mat) %in% non_gene_rows, , drop = FALSE]
  
  # Handle duplicates by median
  if (any(duplicated(rownames(prot_mat)))) {
    cat("  WARNING: Found", sum(duplicated(rownames(prot_mat))), "duplicated genes\n")
    cat("  Collapsing by median...\n")
    
    prot_df <- as.data.frame(prot_mat) %>%
      mutate(Gene = rownames(prot_mat)) %>%
      group_by(Gene) %>%
      summarise(across(everything(), ~median(.x, na.rm = TRUE)), .groups = "drop")
    
    prot_mat <- as.matrix(prot_df[, -1])
    rownames(prot_mat) <- prot_df$Gene
  }
  
  cat("  Matrix:", nrow(prot_mat), "proteins x", ncol(prot_mat), "samples\n\n")
  
  return(prot_mat)
}

check_dgat1_distribution <- function(dgat_vec, winsorize = TRUE, percentile = 0.01) {
  cat("\n=== DGAT1 Distribution ===\n")
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
  
  return(dgat_vec)
}

filter_markers_by_detection <- function(prot_mat, markers_df, min_rate = 0.30) {
  cat("\n=== Filtering Markers by Detection ===\n")
  
  present <- markers_df$Gene[markers_df$Gene %in% rownames(prot_mat)]
  missing <- markers_df$Gene[!markers_df$Gene %in% rownames(prot_mat)]
  
  if (length(present) == 0) {
    stop("No immune markers found in proteomics data")
  }
  
  detection_rates <- sapply(present, function(p) {
    sum(!is.na(prot_mat[p, ])) / ncol(prot_mat)
  })
  
  passed <- present[detection_rates >= min_rate]
  failed <- present[detection_rates < min_rate]
  
  cat("  Total markers:", nrow(markers_df), "\n")
  cat("  Found in data:", length(present), "\n")
  cat("  Passed detection (>", min_rate * 100, "%):", length(passed), "\n")
  cat("  Failed detection:", length(failed), "\n")
  cat("  Missing from data:", length(missing), "\n\n")
  
  present_df <- markers_df %>%
    filter(Gene %in% passed) %>%
    mutate(DetectionRate = detection_rates[Gene])
  
  failed_df <- markers_df %>%
    filter(Gene %in% failed) %>%
    mutate(DetectionRate = detection_rates[Gene])
  
  missing_df <- markers_df %>%
    filter(Gene %in% missing) %>%
    mutate(DetectionRate = NA)
  
  list(
    passed = passed,
    present_df = present_df,
    failed_df = failed_df,
    missing_df = missing_df
  )
}

analyze_individual_markers <- function(prot_mat, dgat_vec, passed_markers, markers_df) {
  
  cat("Analyzing", length(passed_markers), "markers...\n")
  
  results <- lapply(passed_markers, function(protein) {
    vec <- prot_mat[protein, ]
    valid <- !is.na(vec) & !is.na(dgat_vec)
    
    if (sum(valid) < 10) return(NULL)
    
    test <- cor.test(vec[valid], dgat_vec[valid], method = "spearman")
    
    marker_info <- markers_df[markers_df$Gene == protein, ]
    
    data.frame(
      Gene = protein,
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

# =============================================================================
# DUAL PLOTTING FUNCTIONS
# =============================================================================

# PLOT 1: Color by Cell Type
plot_volcano_by_celltype <- function(cor_results, outfile) {
  
  if (is.null(cor_results) || nrow(cor_results) == 0) {
    cat("  Skipping cell type volcano (no data)\n")
    return(NULL)
  }
  
  cor_results$CellType <- factor(cor_results$CellType, 
                                  levels = names(CELLTYPE_COLORS))
  
  cor_results$Label <- ifelse(
    cor_results$FDR < 0.1 | abs(cor_results$Rho) > 0.2,
    cor_results$Gene,
    ""
  )
  
  cor_results$Sig_Level <- case_when(
    cor_results$FDR < 0.001 ~ "FDR < 0.001",
    cor_results$FDR < 0.01 ~ "FDR < 0.01",
    cor_results$FDR < 0.05 ~ "FDR < 0.05",
    cor_results$P < 0.05 ~ "P < 0.05",
    TRUE ~ "NS"
  )
  
  p <- ggplot(cor_results, aes(x = Rho, y = -log10(P))) +
    geom_point(aes(color = CellType, size = Sig_Level, alpha = Sig_Level)) +
    geom_text_repel(
      aes(label = Label),
      size = 3,
      max.overlaps = 25,
      box.padding = 0.5,
      segment.color = "grey50"
    ) +
    scale_color_manual(values = CELLTYPE_COLORS, name = "Cell Type") +
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
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", 
               color = "blue", alpha = 0.5) +
    labs(
      title = "DGAT1 Protein vs Individual Markers",
      subtitle = "Colored by Cell Type (Biological Classification)",
      x = "Spearman ρ (DGAT1 correlation)",
      y = "-log10(P-value)"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      legend.position = "right",
      plot.title = element_text(face = "bold", size = 14)
    )
  
  ggsave(outfile, p, width = 14, height = 8, dpi = 300)
  cat("  ✓ Cell type volcano saved\n")
  
  return(p)
}

# PLOT 2: Color by Immune Function
plot_volcano_by_function <- function(cor_results, outfile) {
  
  if (is.null(cor_results) || nrow(cor_results) == 0) {
    cat("  Skipping function volcano (no data)\n")
    return(NULL)
  }
  
  cor_results$ImmuneFunction <- factor(cor_results$ImmuneFunction, 
                                        levels = names(FUNCTION_COLORS))
  
  cor_results$Label <- ifelse(
    cor_results$FDR < 0.1 | abs(cor_results$Rho) > 0.2,
    cor_results$Gene,
    ""
  )
  
  cor_results$Sig_Level <- case_when(
    cor_results$FDR < 0.001 ~ "FDR < 0.001",
    cor_results$FDR < 0.01 ~ "FDR < 0.01",
    cor_results$FDR < 0.05 ~ "FDR < 0.05",
    cor_results$P < 0.05 ~ "P < 0.05",
    TRUE ~ "NS"
  )
  
  p <- ggplot(cor_results, aes(x = Rho, y = -log10(P))) +
    geom_point(aes(color = ImmuneFunction, size = Sig_Level, alpha = Sig_Level)) +
    geom_text_repel(
      aes(label = Label),
      size = 3,
      max.overlaps = 25,
      box.padding = 0.5,
      segment.color = "grey50"
    ) +
    scale_color_manual(values = FUNCTION_COLORS, name = "Immune Function") +
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
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", 
               color = "blue", alpha = 0.5) +
    labs(
      title = "DGAT1 Protein vs Individual Markers",
      subtitle = "Colored by Immune Function (Pro- vs Anti-tumor)",
      x = "Spearman ρ (DGAT1 correlation)",
      y = "-log10(P-value)"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      legend.position = "right",
      plot.title = element_text(face = "bold", size = 14)
    )
  
  ggsave(outfile, p, width = 14, height = 8, dpi = 300)
  cat("  ✓ Immune function volcano saved\n")
  
  return(p)
}

# Dot plot by cell type
plot_dotplot_by_celltype <- function(cor_results, outfile) {
  
  cor_results <- cor_results %>%
    arrange(Rho) %>%
    mutate(Gene = factor(Gene, levels = Gene))
  
  cor_results$SigLabel <- case_when(
    cor_results$FDR < 0.001 ~ "***",
    cor_results$FDR < 0.01 ~ "**",
    cor_results$FDR < 0.05 ~ "*",
    cor_results$P < 0.05 ~ "†",
    TRUE ~ ""
  )
  
  p <- ggplot(cor_results, aes(y = Gene, x = Rho)) +
    geom_segment(aes(yend = Gene, xend = 0, color = CellType), 
                 linewidth = 0.8) +
    geom_point(aes(color = CellType, size = -log10(P))) +
    geom_text(aes(label = SigLabel), 
              hjust = -0.5, size = 5, fontface = "bold") +
    geom_vline(xintercept = 0, linetype = "solid", color = "black") +
    scale_color_manual(values = CELLTYPE_COLORS) +
    scale_size_continuous(range = c(2, 6), name = "-log10(P)") +
    labs(
      title = "DGAT1 Correlation with Individual Markers",
      subtitle = "Colored by Cell Type | † p<0.05, * FDR<0.05, ** FDR<0.01, *** FDR<0.001",
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
  
  cat("  ✓ Dot plot saved\n")
  
  return(p)
}

# Dot plot by immune function
plot_dotplot_by_function <- function(cor_results, outfile) {
  
  cor_results <- cor_results %>%
    arrange(Rho) %>%
    mutate(Gene = factor(Gene, levels = Gene))
  
  cor_results$SigLabel <- case_when(
    cor_results$FDR < 0.001 ~ "***",
    cor_results$FDR < 0.01 ~ "**",
    cor_results$FDR < 0.05 ~ "*",
    cor_results$P < 0.05 ~ "†",
    TRUE ~ ""
  )
  
  p <- ggplot(cor_results, aes(y = Gene, x = Rho)) +
    geom_segment(aes(yend = Gene, xend = 0, color = ImmuneFunction), 
                 linewidth = 0.8) +
    geom_point(aes(color = ImmuneFunction, size = -log10(P))) +
    geom_text(aes(label = SigLabel), 
              hjust = -0.5, size = 5, fontface = "bold") +
    geom_vline(xintercept = 0, linetype = "solid", color = "black") +
    scale_color_manual(values = FUNCTION_COLORS) +
    scale_size_continuous(range = c(2, 6), name = "-log10(P)") +
    labs(
      title = "DGAT1 Correlation with Individual Markers",
      subtitle = "Colored by Immune Function | † p<0.05, * FDR<0.05, ** FDR<0.01, *** FDR<0.001",
      x = "Spearman ρ",
      y = "",
      color = "Immune Function"
    ) +
    theme_minimal(base_size = 11) +
    theme(
      axis.text.y = element_text(size = 9),
      legend.position = "right",
      panel.grid.major.y = element_line(color = "grey90")
    )
  
  ggsave(outfile, p, width = 10, height = max(6, nrow(cor_results) * 0.25), 
         dpi = 300, limitsize = FALSE)
  
  cat("  ✓ Function dot plot saved\n")
  
  return(p)
}

# Heatmap
plot_clustered_heatmap <- function(prot_mat, dgat_vec, cor_results, outfile) {
  
  if (!COMPLEX_HEATMAP_AVAILABLE) {
    cat("  ⚠ ComplexHeatmap not available, skipping\n")
    return(NULL)
  }
  
  present <- cor_results$Gene[cor_results$Gene %in% rownames(prot_mat)]
  
  if (length(present) < 2) {
    cat("  ⚠ Too few markers for heatmap\n")
    return(NULL)
  }
  
  marker_mat <- prot_mat[present, ]
  marker_z <- t(scale(t(marker_mat)))
  
  valid_samples <- colSums(!is.na(marker_z)) > 0
  marker_z <- marker_z[, valid_samples]
  dgat_aligned <- dgat_vec[colnames(marker_z)]
  
  row_anno <- cor_results %>%
    filter(Gene %in% rownames(marker_z)) %>%
    arrange(match(Gene, rownames(marker_z)))
  
  if (!all(row_anno$Gene == rownames(marker_z))) {
    row_anno <- row_anno[match(rownames(marker_z), row_anno$Gene), ]
  }
  
  ha_row <- rowAnnotation(
    CellType = row_anno$CellType,
    Function = row_anno$ImmuneFunction,
    Correlation = anno_barplot(
      row_anno$Rho,
      border = FALSE,
      gp = gpar(fill = ifelse(row_anno$Rho > 0, "#D73027", "#4575B4"))
    ),
    col = list(
      CellType = CELLTYPE_COLORS,
      Function = FUNCTION_COLORS
    ),
    show_annotation_name = TRUE
  )
  
  dgat_groups <- ifelse(dgat_aligned >= median(dgat_aligned, na.rm = TRUE), 
                        "High", "Low")
  
  ha_col <- HeatmapAnnotation(
    DGAT1 = dgat_groups,
    col = list(DGAT1 = c("High" = "#D73027", "Low" = "#4575B4")),
    show_annotation_name = TRUE
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
    column_title = "Individual Marker Expression (Protein Level)",
    column_title_gp = gpar(fontsize = 14, fontface = "bold")
  )
  
  png(outfile, width = 14, height = max(8, nrow(marker_z) * 0.3), 
      units = "in", res = 300)
  draw(ht)
  dev.off()
  
  cat("  ✓ Heatmap saved\n")
  
  return(ht)
}

# =============================================================================
# MAIN
# =============================================================================

main <- function() {
  cat("\n", paste(rep("=", 70), collapse = ""), "\n")
  cat("CPTAC GBM: Streamlined Individual Marker Analysis\n")
  cat(paste(rep("=", 70), collapse = ""), "\n\n")
  
  # Load data
  prot_mat <- load_protein_data(PROTEIN_FILE)
  
  if (!"DGAT1" %in% rownames(prot_mat)) {
    stop("DGAT1 not found")
  }
  
  dgat <- prot_mat["DGAT1", ]
  dgat <- check_dgat1_distribution(dgat, winsorize = TRUE, 
                                     percentile = WINSORIZE_PERCENTILE)
  
  # Filter markers
  filter_res <- filter_markers_by_detection(prot_mat, CURATED_MARKERS, 
                                              MIN_DETECTION_RATE)
  
  # Save coverage reports
  write_csv(filter_res$present_df, 
            file.path(OUT_DIR, "markers_present.csv"))
  write_csv(filter_res$failed_df, 
            file.path(OUT_DIR, "markers_failed_detection.csv"))
  write_csv(filter_res$missing_df, 
            file.path(OUT_DIR, "markers_missing.csv"))
  
  # Analyze correlations
  cat("Analyzing individual markers...\n")
  cor_results <- analyze_individual_markers(prot_mat, dgat, 
                                             filter_res$passed, CURATED_MARKERS)
  
  write_csv(cor_results, file.path(OUT_DIR, "individual_marker_correlations.csv"))
  
  cat("\nMarkers detected:", nrow(cor_results), "\n")
  cat("Significant (FDR < 0.05):", sum(cor_results$FDR < 0.05), "\n")
  cat("Nominal (P < 0.05):", sum(cor_results$P < 0.05), "\n\n")
  
  # Summary by cell type
  cat("Summary by Cell Type:\n")
  summary_celltype <- cor_results %>%
    group_by(CellType) %>%
    summarise(
      N_markers = n(),
      Mean_rho = mean(Rho),
      N_positive = sum(Rho > 0),
      N_negative = sum(Rho < 0),
      N_sig = sum(FDR < 0.05),
      .groups = "drop"
    )
  print(summary_celltype)
  
  cat("\n\nSummary by Immune Function:\n")
  summary_function <- cor_results %>%
    group_by(ImmuneFunction) %>%
    summarise(
      N_markers = n(),
      Mean_rho = mean(Rho),
      N_positive = sum(Rho > 0),
      N_negative = sum(Rho < 0),
      N_sig = sum(FDR < 0.05),
      .groups = "drop"
    )
  print(summary_function)
  
  cat("\nTop 10 correlations:\n")
  print(head(cor_results %>% select(Gene, CellType, ImmuneFunction, Rho, P, FDR), 10))
  
  # Generate dual plots
  cat("\n=== Generating Dual Visualizations ===\n")
  
  # Volcano plots
  plot_volcano_by_celltype(
    cor_results,
    file.path(OUT_DIR, "volcano_by_celltype.png")
  )
  
  plot_volcano_by_function(
    cor_results,
    file.path(OUT_DIR, "volcano_by_function.png")
  )
  
  # Dot plots
  plot_dotplot_by_celltype(
    cor_results,
    file.path(OUT_DIR, "dotplot_by_celltype.png")
  )
  
  plot_dotplot_by_function(
    cor_results,
    file.path(OUT_DIR, "dotplot_by_function.png")
  )
  
  # Heatmap
  plot_clustered_heatmap(
    prot_mat, dgat, cor_results,
    file.path(OUT_DIR, "heatmap_markers.png")
  )
  
  cat("\n", paste(rep("=", 70), collapse = ""), "\n")
  cat("Analysis complete. Results in:", OUT_DIR, "\n")
  cat(paste(rep("=", 70), collapse = ""), "\n\n")
  
  invisible(cor_results)
}

# Run
if (!interactive()) {
  results <- main()
}
