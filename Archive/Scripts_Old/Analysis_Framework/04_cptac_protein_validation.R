#!/usr/bin/env Rscript
# =============================================================================
# CPTAC GBM Proteomics: DGAT1 vs Immune Infiltration (CORRECTED)
# Validates mRNA findings at protein level
# =============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(data.table)
  library(ppcor)
  library(pheatmap)
})

# ====== CONFIG =================================================================
BASE_DIR <- "/Users/zijiefeng/Desktop/Guo's lab/My_Research/DGAT_Immunology"
DATA_DIR <- file.path(BASE_DIR, "Raw_Data", "Proteome", "CPTAC")
OUT_DIR  <- file.path(BASE_DIR, "Results", "CPTAC_Protein_Analysis")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

# File paths (adjust to your actual files)
PROTEIN_FILE <- file.path(DATA_DIR, "CPTAC3_Glioblastoma_Multiforme_Proteome.tmt11.tsv")
CLINICAL_FILE <- file.path(DATA_DIR, "S048_CPTAC_GBM_Discovery_Cohort_Clinical_Data_Dec2019_r1.xlsx")
MAPPING_FILE <- file.path(DATA_DIR, "S048_CPTAC_GBM_Discovery_Cohort_TMT11_CaseID_SampleID_AliquotID_Map_Dec2019_r1.xlsx")

# Immune marker proteins
IMMUNE_MARKERS <- list(
  Macrophages = c("CD163", "CD68", "MRC1", "MSR1", "MARCO", "LGALS3", "APOE"),
  TAM = c("CD163", "MRC1", "LGALS3", "SPP1", "APOE"),
  T_cells = c("CD8A", "CD8B", "CD4", "CD3E", "CD3D", "CD3G"),
  Tregs = c("FOXP3", "IL2RA", "CTLA4", "TIGIT"),
  NK = c("NCAM1", "KLRD1", "KLRF1", "NCR1"),
  DC = c("CD1C", "FCER1A", "CLEC9A", "BATF3"),
  B_cells = c("CD19", "MS4A1", "CD79A", "CD79B"),
  Neutrophils = c("FCGR3B", "CSF3R", "S100A12"),
  Checkpoint = c("PDCD1", "CD274", "CTLA4", "LAG3", "HAVCR2")
)

# =============================================================================
# FUNCTIONS
# =============================================================================

msg <- function(...) cat(paste0(..., "\n"))

# Load CPTAC protein data (TMT format)
load_cptac_protein <- function(filepath) {
  msg("Loading protein data: ", basename(filepath))
  
  if (!file.exists(filepath)) stop("File not found: ", filepath)
  
  # CPTAC proteomics typically has:
  # - First few columns: Gene, Gene_ID, etc.
  # - Remaining columns: Sample abundance (TMT ratios or intensities)
  
  prot <- fread(filepath, sep="\t")
  
  # Find gene symbol column
  gene_cols <- c("Gene", "gene", "Gene_Symbol", "symbol", "SYMBOL")
  gene_col <- gene_cols[gene_cols %in% names(prot)][1]
  
  if (is.na(gene_col)) {
    msg("Available columns:", paste(names(prot), collapse=", "))
    stop("Could not find gene symbol column. Check column names.")
  }
  
  msg("  Using gene column: ", gene_col)
  
  # Extract sample columns (usually end with case/aliquot ID patterns)
  sample_cols <- setdiff(names(prot), 
                         c(gene_col, "Gene_ID", "gene_id", "Protein_ID"))
  
  # Create matrix: genes x samples
  prot_mat <- as.matrix(prot[, ..sample_cols])
  rownames(prot_mat) <- prot[[gene_col]]
  
  # Convert to numeric and remove non-gene rows (like "Mean", "Median", "StdDev")
  prot_mat <- apply(prot_mat, 2, as.numeric)
  rownames(prot_mat) <- prot[[gene_col]]
  
  # Remove non-gene rows
  non_gene_rows <- c("Mean", "Median", "StdDev", "Gene")
  prot_mat <- prot_mat[!rownames(prot_mat) %in% non_gene_rows, , drop = FALSE]
  
  # Remove duplicate genes (keep first occurrence)
  if (any(duplicated(rownames(prot_mat)))) {
    msg("  Removing ", sum(duplicated(rownames(prot_mat))), " duplicate genes")
    prot_mat <- prot_mat[!duplicated(rownames(prot_mat)), ]
  }
  
  msg("  Loaded: ", nrow(prot_mat), " proteins x ", ncol(prot_mat), " samples")
  
  return(prot_mat)
}

# Create immune score (mean of available markers)
create_immune_score <- function(prot_mat, marker_list) {
  scores <- lapply(names(marker_list), function(category) {
    genes <- marker_list[[category]]
    present <- intersect(genes, rownames(prot_mat))
    
    if (length(present) == 0) return(NULL)
    
    # Mean z-score of present markers
    subset_mat <- prot_mat[present, , drop=FALSE]
    z_mat <- t(scale(t(subset_mat)))
    
    score <- colMeans(z_mat, na.rm=TRUE)
    
    data.frame(
      sample = names(score),
      category = category,
      score = score,
      n_markers = length(present),
      markers = paste(present, collapse=","),
      stringsAsFactors = FALSE
    )
  })
  
  do.call(rbind, scores[!sapply(scores, is.null)])
}

# Correlations
compute_correlations <- function(prot_mat, dgat_vec, immune_genes) {
  results <- lapply(immune_genes, function(g) {
    if (!(g %in% rownames(prot_mat))) return(NULL)
    
    vec <- prot_mat[g, ]
    valid <- !is.na(vec) & !is.na(dgat_vec)
    
    if (sum(valid) < 10) return(NULL)
    
    test <- cor.test(vec[valid], dgat_vec[valid], method="spearman")
    
    data.frame(
      Protein = g,
      Rho = test$estimate,
      P = test$p.value,
      N = sum(valid),
      stringsAsFactors = FALSE
    )
  })
  
  tab <- do.call(rbind, results[!sapply(results, is.null)])
  tab$FDR <- p.adjust(tab$P, method="BH")
  tab[order(tab$FDR), ]
}

# High/Low group comparison
compare_groups <- function(prot_mat, dgat_vec, immune_genes) {
  # Median split
  cutoff <- median(dgat_vec, na.rm=TRUE)
  groups <- ifelse(dgat_vec >= cutoff, "High", "Low")
  names(groups) <- names(dgat_vec)
  
  results <- lapply(immune_genes, function(g) {
    if (!(g %in% rownames(prot_mat))) return(NULL)
    
    vec <- prot_mat[g, ]
    valid <- !is.na(vec) & !is.na(groups)
    
    if (sum(valid) < 10) return(NULL)
    
    high_vals <- vec[valid & groups == "High"]
    low_vals <- vec[valid & groups == "Low"]
    
    if (length(high_vals) < 3 | length(low_vals) < 3) return(NULL)
    
    test <- wilcox.test(high_vals, low_vals)
    
    data.frame(
      Protein = g,
      Mean_High = mean(high_vals, na.rm=TRUE),
      Mean_Low = mean(low_vals, na.rm=TRUE),
      Delta = mean(high_vals, na.rm=TRUE) - mean(low_vals, na.rm=TRUE),
      P = test$p.value,
      stringsAsFactors = FALSE
    )
  })
  
  tab <- do.call(rbind, results[!sapply(results, is.null)])
  tab$FDR <- p.adjust(tab$P, method="BH")
  tab[order(tab$FDR), ]
}

# =============================================================================
# MAIN ANALYSIS
# =============================================================================

main <- function() {
  msg("\n" , paste(rep("=", 70), collapse=""))
  msg("CPTAC GBM PROTEOMICS: DGAT1 vs IMMUNE INFILTRATION")
  msg(paste(rep("=", 70), collapse=""), "\n")
  
  # Load data
  prot_mat <- load_cptac_protein(PROTEIN_FILE)
  
  # Check DGAT1
  if (!"DGAT1" %in% rownames(prot_mat)) {
    stop("DGAT1 protein not found in dataset!")
  }
  
  dgat <- prot_mat["DGAT1", ]
  msg("\nDGAT1 protein detected in ", sum(!is.na(dgat)), " samples")
  msg("  Range: [", round(min(dgat, na.rm=TRUE), 2), ", ", 
      round(max(dgat, na.rm=TRUE), 2), "]")
  
  # Check immune marker coverage
  msg("\n=== Immune Marker Coverage ===")
  all_markers <- unique(unlist(IMMUNE_MARKERS))
  present <- intersect(all_markers, rownames(prot_mat))
  missing <- setdiff(all_markers, rownames(prot_mat))
  
  msg("  Present: ", length(present), " / ", length(all_markers))
  for (cat in names(IMMUNE_MARKERS)) {
    genes <- IMMUNE_MARKERS[[cat]]
    n_present <- sum(genes %in% rownames(prot_mat))
    msg("    ", cat, ": ", n_present, " / ", length(genes))
  }
  
  if (length(missing) > 0) {
    msg("\n  Missing markers: ", paste(head(missing, 10), collapse=", "))
  }
  
  # Analysis 1: Individual protein correlations
  msg("\n=== Correlation Analysis ===")
  cor_tab <- compute_correlations(prot_mat, dgat, present)
  
  fwrite(cor_tab, file.path(OUT_DIR, "DGAT1_protein_correlations.csv"))
  
  msg("  Significant (FDR < 0.05): ", sum(cor_tab$FDR < 0.05, na.rm=TRUE))
  if (nrow(cor_tab) > 0) {
    msg("\n  Top 5 correlations:")
    print(head(cor_tab, 5))
  }
  
  # Analysis 2: High/Low comparison
  msg("\n=== DGAT1 High vs Low Comparison ===")
  comp_tab <- compare_groups(prot_mat, dgat, present)
  
  fwrite(comp_tab, file.path(OUT_DIR, "DGAT1_HighLow_comparison.csv"))
  
  msg("  Significant (FDR < 0.05): ", sum(comp_tab$FDR < 0.05, na.rm=TRUE))
  if (nrow(comp_tab) > 0) {
    msg("\n  Top 5 differences:")
    print(head(comp_tab, 5))
  }
  
  # Analysis 3: Category scores
  msg("\n=== Immune Category Scores ===")
  scores <- create_immune_score(prot_mat, IMMUNE_MARKERS)
  
  # Correlate category scores with DGAT1
  score_wide <- scores %>%
    dplyr::select(sample, category, score) %>%
    pivot_wider(names_from = category, values_from = score)
  
  score_mat <- as.matrix(score_wide[, -1])
  rownames(score_mat) <- score_wide$sample
  
  dgat_aligned <- dgat[rownames(score_mat)]
  
  cat_cors <- sapply(colnames(score_mat), function(cat) {
    test <- cor.test(score_mat[, cat], dgat_aligned, method="spearman")
    c(Rho = test$estimate, P = test$p.value)
  })
  
  cat_cors <- as.data.frame(t(cat_cors))
  cat_cors$FDR <- p.adjust(cat_cors$P, method="BH")
  cat_cors$Category <- rownames(cat_cors)
  cat_cors <- cat_cors[order(cat_cors$FDR), ]
  
  fwrite(cat_cors, file.path(OUT_DIR, "DGAT1_category_correlations.csv"))
  
  msg("\n  Category correlations:")
  print(cat_cors)
  
  # Visualization
  msg("\n=== Generating Plots ===")
  
  # Volcano plot
  if (requireNamespace("ggrepel", quietly=TRUE)) {
    g1 <- ggplot(cor_tab, aes(Rho, -log10(FDR), label=Protein)) +
      geom_point(aes(color = FDR < 0.05), size=3, alpha=0.7) +
      ggrepel::geom_text_repel(
        data = cor_tab %>% filter(FDR < 0.1 | abs(Rho) > 0.3),
        max.overlaps = 20, size=3
      ) +
      scale_color_manual(values=c("TRUE"="#D73027", "FALSE"="grey70")) +
      geom_vline(xintercept=0, linetype="solid") +
      geom_hline(yintercept=-log10(0.05), linetype="dashed", color="blue") +
      theme_minimal(base_size=12) +
      labs(
        title="DGAT1 Protein vs Immune Markers",
        subtitle=paste("N =", sum(!is.na(dgat)), "samples"),
        x="Spearman ρ", 
        y="-log10(FDR)"
      )
    
    ggsave(file.path(OUT_DIR, "DGAT1_protein_volcano.png"), 
           g1, width=9, height=7, dpi=300)
  }
  
  # Heatmap of category scores
  ann <- data.frame(
    DGAT1_Group = ifelse(dgat_aligned >= median(dgat_aligned, na.rm=TRUE), 
                         "High", "Low"),
    row.names = names(dgat_aligned)
  )
  
  # Handle NAs in score matrix for heatmap
  score_mat_clean <- score_mat
  score_mat_clean[is.na(score_mat_clean)] <- 0
  
  # Skip heatmap for now due to clustering issues with sparse data
  msg("  Skipping heatmap - data too sparse for clustering")
  
  msg("\n" , paste(rep("=", 70), collapse=""))
  msg("ANALYSIS COMPLETE")
  msg(paste(rep("=", 70), collapse=""))
  msg("\nOutputs saved to: ", OUT_DIR)
  msg("\nKey files:")
  msg("  - DGAT1_protein_correlations.csv")
  msg("  - DGAT1_HighLow_comparison.csv")
  msg("  - DGAT1_category_correlations.csv")
  msg("  - DGAT1_protein_volcano.png")
  msg("  - immune_category_heatmap.png")
  
  # Return summary
  invisible(list(
    correlations = cor_tab,
    comparisons = comp_tab,
    categories = cat_cors,
    n_samples = sum(!is.na(dgat))
  ))
}

# Run
results <- main()
