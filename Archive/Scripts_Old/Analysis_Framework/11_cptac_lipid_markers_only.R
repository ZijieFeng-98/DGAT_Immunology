#!/usr/bin/env Rscript
# =============================================================================
# CPTAC GBM: LIPID METABOLISM MARKERS Analysis Only
# Proteomics data
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
OUT_DIR  <- file.path(BASE_DIR, "Results", "CPTAC_Lipid_Markers")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

PROTEIN_FILE <- file.path(DATA_DIR, "CPTAC3_Glioblastoma_Multiforme_Proteome.tmt11.tsv")

MIN_DETECTION_RATE <- 0.30
WINSORIZE_PERCENTILE <- 0.01

# =============================================================================
# LIPID METABOLISM MARKERS ONLY
# =============================================================================

LIPID_MARKERS <- tribble(
  ~Gene,        ~Pathway,             ~Category,              ~Description,
  
  # Core lipid droplet formation
  "DGAT1",      "Lipid Droplet",      "LD-Formation",         "Lipid droplet formation (target gene)",
  
  # SREBP pathway
  "SREBF1",     "SREBP-Pathway",      "Master-Regulator",     "SREBP-1 (master lipogenesis regulator)",
  "SCAP",       "SREBP-Pathway",      "Escort-Protein",       "SREBP cleavage-activating protein",
  
  # Fatty acid synthesis
  "ACACA",      "FA-Synthesis",       "Rate-Limiting",        "Acetyl-CoA carboxylase (ACC)",
  "FASN",       "FA-Synthesis",       "De-Novo-Synthesis",    "Fatty acid synthase",
  "SCD",        "FA-Synthesis",       "Desaturation",         "Stearoyl-CoA desaturase 1 (SCD1)",
  
  # Cholesterol synthesis
  "HMGCR",      "Cholesterol",        "Rate-Limiting",        "HMG-CoA reductase",
  
  # Hypoxia
  "HIF1A",      "Hypoxia",            "Metabolic-Rewiring",   "Hypoxia-inducible factor 1-alpha"
)

# =============================================================================
# COLOR PALETTES
# =============================================================================

PATHWAY_COLORS <- c(
  "Lipid Droplet" = "#8B0000",
  "SREBP-Pathway" = "#DC143C",
  "FA-Synthesis" = "#CD5C5C",
  "Cholesterol" = "#F08080",
  "Hypoxia" = "#2F4F4F"
)

CATEGORY_COLORS <- c(
  "LD-Formation" = "#8B0000",
  "Master-Regulator" = "#DC143C",
  "Escort-Protein" = "#CD5C5C",
  "Rate-Limiting" = "#B22222",
  "De-Novo-Synthesis" = "#F08080",
  "Desaturation" = "#FA8072",
  "Metabolic-Rewiring" = "#2F4F4F"
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
  
  if (any(duplicated(rownames(prot_mat)))) {
    cat("  Collapsing", sum(duplicated(rownames(prot_mat))), "duplicates\n")
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
  cat("Samples:", sum(!is.na(dgat_vec)), "/", length(dgat_vec), "\n")
  cat("Range: [", round(min(dgat_vec, na.rm = TRUE), 3), ",", 
      round(max(dgat_vec, na.rm = TRUE), 3), "]\n")
  cat("Mean ± SD:", round(mean(dgat_vec, na.rm = TRUE), 3), "±", 
      round(sd(dgat_vec, na.rm = TRUE), 3), "\n")
  
  if (winsorize && percentile > 0) {
    lower <- quantile(dgat_vec, percentile, na.rm = TRUE)
    upper <- quantile(dgat_vec, 1 - percentile, na.rm = TRUE)
    n_winsorized <- sum(dgat_vec < lower | dgat_vec > upper, na.rm = TRUE)
    if (n_winsorized > 0) {
      cat("Winsorizing", n_winsorized, "outliers\n")
      dgat_vec[dgat_vec < lower] <- lower
      dgat_vec[dgat_vec > upper] <- upper
    }
  }
  cat("================================\n\n")
  return(dgat_vec)
}

filter_markers <- function(prot_mat, markers_df, min_rate = 0.30) {
  cat("\n=== Filtering Markers ===\n")
  present <- markers_df$Gene[markers_df$Gene %in% rownames(prot_mat)]
  missing <- markers_df$Gene[!markers_df$Gene %in% rownames(prot_mat)]
  
  detection_rates <- sapply(present, function(p) {
    sum(!is.na(prot_mat[p, ])) / ncol(prot_mat)
  })
  
  passed <- present[detection_rates >= min_rate]
  failed <- present[detection_rates < min_rate]
  
  cat("  Total:", nrow(markers_df), "| Found:", length(present), 
      "| Passed:", length(passed), "| Missing:", length(missing), "\n\n")
  
  list(
    passed = passed,
    present_df = markers_df %>% filter(Gene %in% passed) %>% 
                  mutate(DetectionRate = detection_rates[Gene])
  )
}

analyze_markers <- function(prot_mat, dgat_vec, passed, markers_df) {
  cat("Analyzing", length(passed), "lipid markers...\n")
  results <- lapply(passed, function(protein) {
    vec <- prot_mat[protein, ]
    valid <- !is.na(vec) & !is.na(dgat_vec)
    if (sum(valid) < 10) return(NULL)
    
    test <- cor.test(vec[valid], dgat_vec[valid], method = "spearman")
    marker_info <- markers_df[markers_df$Gene == protein, ]
    
    data.frame(
      Gene = protein,
      Pathway = marker_info$Pathway,
      Category = marker_info$Category,
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

plot_volcano_pathway <- function(cor_results, outfile) {
  if (is.null(cor_results) || nrow(cor_results) == 0) return(NULL)
  cor_results$Pathway <- factor(cor_results$Pathway, levels = names(PATHWAY_COLORS))
  cor_results$Label <- cor_results$Gene
  
  p <- ggplot(cor_results, aes(x = Rho, y = -log10(P))) +
    geom_point(aes(color = Pathway, size = -log10(P), 
                   alpha = ifelse(FDR < 0.05, 1, 0.5))) +
    geom_text_repel(aes(label = Label), size = 4, max.overlaps = 20, fontface = "bold") +
    scale_color_manual(values = PATHWAY_COLORS, name = "Pathway") +
    scale_size_continuous(range = c(3, 8), guide = "none") +
    scale_alpha_identity() +
    geom_vline(xintercept = 0, color = "black", linewidth = 1) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "blue") +
    labs(title = "DGAT1 Protein vs Lipid Metabolism Genes (CPTAC GBM)",
         subtitle = "Colored by Metabolic Pathway",
         x = "Spearman ρ", y = "-log10(P-value)") +
    theme_minimal(base_size = 14) +
    theme(legend.position = "right", plot.title = element_text(face = "bold", size = 16))
  
  ggsave(outfile, p, width = 12, height = 8, dpi = 300)
  cat("  ✓ Pathway volcano\n")
  return(p)
}

plot_volcano_category <- function(cor_results, outfile) {
  if (is.null(cor_results) || nrow(cor_results) == 0) return(NULL)
  cor_results$Category <- factor(cor_results$Category, levels = names(CATEGORY_COLORS))
  cor_results$Label <- cor_results$Gene
  
  p <- ggplot(cor_results, aes(x = Rho, y = -log10(P))) +
    geom_point(aes(color = Category, size = -log10(P), 
                   alpha = ifelse(FDR < 0.05, 1, 0.5))) +
    geom_text_repel(aes(label = Label), size = 4, max.overlaps = 20, fontface = "bold") +
    scale_color_manual(values = CATEGORY_COLORS, name = "Function") +
    scale_size_continuous(range = c(3, 8), guide = "none") +
    scale_alpha_identity() +
    geom_vline(xintercept = 0, color = "black", linewidth = 1) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "blue") +
    labs(title = "DGAT1 Protein vs Lipid Metabolism Genes (CPTAC GBM)",
         subtitle = "Colored by Functional Category",
         x = "Spearman ρ", y = "-log10(P-value)") +
    theme_minimal(base_size = 14) +
    theme(legend.position = "right", plot.title = element_text(face = "bold", size = 16))
  
  ggsave(outfile, p, width = 12, height = 8, dpi = 300)
  cat("  ✓ Category volcano\n")
  return(p)
}

plot_barplot <- function(cor_results, outfile) {
  cor_results <- cor_results %>% arrange(Rho) %>% 
    mutate(Gene = factor(Gene, levels = Gene))
  
  cor_results$SigLabel <- case_when(
    cor_results$FDR < 0.001 ~ "***", cor_results$FDR < 0.01 ~ "**",
    cor_results$FDR < 0.05 ~ "*", cor_results$P < 0.05 ~ "†", TRUE ~ ""
  )
  
  p <- ggplot(cor_results, aes(y = Gene, x = Rho)) +
    geom_col(aes(fill = Pathway), alpha = 0.8) +
    geom_text(aes(label = SigLabel), hjust = ifelse(cor_results$Rho > 0, -0.2, 1.2), 
              size = 6, fontface = "bold") +
    geom_vline(xintercept = 0, color = "black", linewidth = 1) +
    scale_fill_manual(values = PATHWAY_COLORS) +
    labs(title = "DGAT1 Protein Correlation with Lipid Genes (CPTAC)",
         subtitle = "† p<0.05, * FDR<0.05, ** FDR<0.01, *** FDR<0.001",
         x = "Spearman ρ", y = "") +
    theme_minimal(base_size = 13) +
    theme(axis.text.y = element_text(size = 12, face = "bold"),
          legend.position = "right",
          panel.grid.major.y = element_blank())
  
  ggsave(outfile, p, width = 10, height = 6, dpi = 300)
  cat("  ✓ Barplot\n")
  return(p)
}

plot_heatmap <- function(prot_mat, dgat_vec, cor_results, outfile) {
  if (!COMPLEX_HEATMAP_AVAILABLE) {
    cat("  ⚠ ComplexHeatmap not available\n")
    return(NULL)
  }
  
  present <- cor_results$Gene[cor_results$Gene %in% rownames(prot_mat)]
  if (length(present) < 2) return(NULL)
  
  marker_mat <- prot_mat[present, ]
  marker_z <- t(scale(t(marker_mat)))
  valid_samples <- colSums(!is.na(marker_z)) > 0
  marker_z <- marker_z[, valid_samples]
  dgat_aligned <- dgat_vec[colnames(marker_z)]
  
  row_anno <- cor_results %>% filter(Gene %in% rownames(marker_z)) %>%
    arrange(match(Gene, rownames(marker_z)))
  
  ha_row <- rowAnnotation(
    Pathway = row_anno$Pathway,
    Correlation = anno_barplot(row_anno$Rho, border = FALSE,
      gp = gpar(fill = ifelse(row_anno$Rho > 0, "#D73027", "#4575B4"))),
    col = list(Pathway = PATHWAY_COLORS)
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
    row_names_gp = gpar(fontsize = 11, fontface = "bold"),
    column_title = "Lipid Metabolism Protein Expression (CPTAC GBM)",
    column_title_gp = gpar(fontsize = 14, fontface = "bold")
  )
  
  png(outfile, width = 12, height = 8, units = "in", res = 300)
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
  cat("CPTAC GBM: Lipid Metabolism Marker Analysis\n")
  cat(paste(rep("=", 70), collapse = ""), "\n\n")
  
  prot_mat <- load_protein_data(PROTEIN_FILE)
  if (!"DGAT1" %in% rownames(prot_mat)) stop("DGAT1 not found")
  
  dgat <- prot_mat["DGAT1", ]
  dgat <- check_dgat1_distribution(dgat, winsorize = TRUE, percentile = WINSORIZE_PERCENTILE)
  
  filter_res <- filter_markers(prot_mat, LIPID_MARKERS, MIN_DETECTION_RATE)
  write_csv(filter_res$present_df, file.path(OUT_DIR, "lipid_markers_present.csv"))
  
  cor_results <- analyze_markers(prot_mat, dgat, filter_res$passed, LIPID_MARKERS)
  write_csv(cor_results, file.path(OUT_DIR, "lipid_marker_correlations.csv"))
  
  cat("\nMarkers:", nrow(cor_results), "| Sig (FDR<0.05):", sum(cor_results$FDR < 0.05), "\n\n")
  
  cat("By Pathway:\n")
  print(cor_results %>% group_by(Pathway) %>%
    summarise(N = n(), Mean_rho = mean(Rho), N_sig = sum(FDR < 0.05), .groups = "drop"))
  
  cat("\nAll Results:\n")
  print(cor_results %>% select(Gene, Pathway, Category, Rho, P, FDR))
  
  cat("\n=== Generating Plots ===\n")
  plot_volcano_pathway(cor_results, file.path(OUT_DIR, "volcano_by_pathway.png"))
  plot_volcano_category(cor_results, file.path(OUT_DIR, "volcano_by_category.png"))
  plot_barplot(cor_results, file.path(OUT_DIR, "barplot_correlations.png"))
  plot_heatmap(prot_mat, dgat, cor_results, file.path(OUT_DIR, "heatmap_lipid_proteins.png"))
  
  cat("\n", paste(rep("=", 70), collapse = ""), "\n")
  cat("Complete! Results in:", OUT_DIR, "\n")
  cat(paste(rep("=", 70), collapse = ""), "\n\n")
  
  invisible(cor_results)
}

if (!interactive()) results <- main()
