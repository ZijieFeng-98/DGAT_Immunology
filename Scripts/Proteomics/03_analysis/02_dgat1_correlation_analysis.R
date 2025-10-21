#!/usr/bin/env Rscript
# ============================================================================
# DGAT1 in Cancer Immunology - CPTAC GBM Proteomics Analysis
# ============================================================================
# Purpose: Analyze DGAT1 correlations with immune and lipid metabolism markers
# Data: CPTAC3 GBM Discovery Cohort (10,804 proteins × 110 samples)
# Date: 2025-10-10
# Status: VALIDATED & ANALYSIS-READY
# ============================================================================

# Load required libraries
suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(pheatmap)
  library(RColorBrewer)
  library(ggrepel)  # For volcano plot labels
})

# Set working directory
BASE_DIR <- "/Users/zijiefeng/Desktop/Guo's lab/My_Research/DGAT_Immunology"
DATA_DIR <- file.path(BASE_DIR, "Processed_Data/CPTAC_GBM_Proteomics")
OUTPUT_DIR <- file.path(BASE_DIR, "Results/Proteomics/DGAT1_Correlations")

# Create output directory
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)
setwd(OUTPUT_DIR)

# Set seed for reproducibility
set.seed(1)

# Harden IO & environment - check files exist
protein_file <- file.path(DATA_DIR, "protein_matrix_cleaned.rds")
metadata_file <- file.path(DATA_DIR, "sample_metadata.csv")

stopifnot(file.exists(protein_file), file.exists(metadata_file))

# ============================================================================
# PART 1: DATA LOADING AND QUALITY CHECK
# ============================================================================

cat("\n=== LOADING CPTAC GBM PROTEOMICS DATA ===\n")

# Load protein matrix
protein_matrix <- readRDS(protein_file)

# Coerce to matrix/numeric to ensure proper data type
protein_matrix <- as.matrix(protein_matrix)
storage.mode(protein_matrix) <- "numeric"

# Collapse duplicate gene symbols by median (if any)
if (any(duplicated(rownames(protein_matrix)))) {
  message("Collapsing duplicated gene symbols by median...")
  protein_dt <- as.data.frame(protein_matrix)
  protein_dt$Gene <- rownames(protein_matrix)
  protein_dt <- aggregate(. ~ Gene, data = protein_dt, FUN = function(x) median(x, na.rm=TRUE))
  rownames(protein_dt) <- protein_dt$Gene
  protein_dt$Gene <- NULL
  protein_matrix <- as.matrix(protein_dt)
}

# Uppercase gene symbols for consistency
rownames(protein_matrix) <- toupper(rownames(protein_matrix))

cat("Protein matrix dimensions:", dim(protein_matrix), "\n")

# Load sample metadata
sample_meta <- fread(metadata_file)

# Align samples between metadata and matrix
id_col <- intersect(colnames(sample_meta), c("Case_ID", "SampleID","sample_id","SpecimenID","ProteomicsID"))
if (length(id_col) >= 1) {
  sid <- sample_meta[[id_col[1]]]
  overlap <- intersect(colnames(protein_matrix), sid)
  if (length(overlap) < 10) {
    warning("Very few overlapping samples between metadata and matrix: ", length(overlap))
  } else {
    cat("Sample overlap:", length(overlap), "samples matched\n")
  }
  # Reorder protein_matrix to metadata order where possible
  keep <- sid[sid %in% colnames(protein_matrix)]
  protein_matrix <- protein_matrix[, keep, drop=FALSE]
} else {
  warning("No sample ID column recognized in metadata; proceeding without alignment.")
}

cat("Sample metadata:", nrow(sample_meta), "samples\n")

# Check DGAT1 availability - CRITICAL
if (!("DGAT1" %in% rownames(protein_matrix))) {
  stop("DGAT1 not found in proteomics rownames. Check gene symbol mapping.")
}

dgat1_values <- as.numeric(protein_matrix["DGAT1", ])

if (all(is.na(dgat1_values))) {
  stop("DGAT1 present but all values are NA. Aborting correlation analysis.")
}

# DGAT1 summary
cat("DGAT1 detected: TRUE\n")
cat("DGAT1 summary:\n")
print(summary(dgat1_values))
cat("Non-NA values:", sum(!is.na(dgat1_values)), "/", length(dgat1_values), "\n")

# Minimum threshold check
if (sum(!is.na(dgat1_values)) < 30) {
  warning("DGAT1 has <30 non-NA samples. Correlation power may be limited.")
}

# ============================================================================
# PART 2: DEFINE MARKER GENE SETS
# ============================================================================

cat("\n=== DEFINING MARKER GENE SETS ===\n")

# -----------------------------
# IMMUNE MARKERS (Based on project knowledge)
# -----------------------------

immune_markers <- list(
  # TAM (Tumor-Associated Macrophages) - M2 markers
  TAM_M2 = c("CD163", "MRC1", "ARG1", "CHI3L1", "IL10", "TGFB1", 
             "PTGS2", "SEPP1", "CLEC7A"),
  
  # Pan-macrophage markers
  TAM_Pan = c("CD68", "CD14", "ITGAM", "AIF1"),
  
  # M1 macrophage markers (should be low in GBM)
  TAM_M1 = c("NOS2", "TNF", "IL1B", "IL12B", "FCGR3A"),
  
  # Immune checkpoints - Ligands
  Checkpoint_Ligand = c("CD274", "PDCD1LG2", "LGALS9", "NECTIN2", "VSIR"),
  
  # Immune checkpoints - Receptors (NOTE: ICOS not CD278)
  Checkpoint_Receptor = c("PDCD1", "CTLA4", "HAVCR2", "LAG3", "TIGIT", "ICOS"),
  
  # CD8+ T cell markers (effector)
  CD8_Effector = c("CD8A", "CD8B", "GZMB", "PRF1", "IFNG", "TNF"),
  
  # T cell exhaustion markers
  T_Exhaustion = c("PDCD1", "HAVCR2", "LAG3", "TIGIT", "TOX"),
  
  # Immunosuppressive cytokines/chemokines
  Immunosuppressive = c("IL10", "TGFB1", "CCL2", "CSF1", "CXCL12", "IDO1"),
  
  # MDSCs markers
  MDSC = c("S100A8", "S100A9", "ARG1", "NOS2")
)

# Uppercase all markers for consistency
immune_markers <- lapply(immune_markers, toupper)

# Flatten all immune markers
all_immune_markers <- unique(unlist(immune_markers))
cat("Total unique immune markers defined:", length(all_immune_markers), "\n")

# -----------------------------
# LIPID METABOLISM MARKERS
# -----------------------------

lipid_markers <- list(
  # Lipid droplet proteins (Deliang Guo's lab focus)
  Lipid_Droplet = c("PLIN1", "PLIN2", "PLIN3", "PLIN4", "PLIN5"),
  
  # Fatty acid oxidation (FAO)
  FAO = c("CPT1A", "CPT1B", "CPT1C", "CPT2", "ACADVL", "HADHA", "HADHB"),
  
  # Fatty acid synthesis (Lipogenesis)
  Lipogenesis = c("FASN", "ACACA", "ACACB", "SCD", "SCD5", "ACLY", "ELOVL6"),
  
  # Fatty acid uptake and activation
  FA_Uptake = c("CD36", "FABP3", "FABP4", "FABP5", "FABP7", 
                "ACSL1", "ACSL3", "ACSL4", "ACSL5"),
  
  # ER stress markers
  ER_Stress = c("HSPA5", "DDIT3", "ATF4", "ATF6", "XBP1", "ERN1"),
  
  # SREBP pathway (DGAT2 related)
  SREBP = c("SREBF1", "SREBF2", "SCAP", "INSIG1", "INSIG2"),
  
  # Cholesterol metabolism
  Cholesterol = c("HMGCR", "HMGCS1", "LDLR", "SCARB1", "NPC1", "SOAT1"),
  
  # Triglyceride synthesis (DGAT family)
  TG_Synthesis = c("DGAT1", "DGAT2", "MOGAT2", "MOGAT3", "GPAM", "AGPAT2"),
  
  # Lipolysis
  Lipolysis = c("PNPLA2", "LIPE", "MGLL", "ABHD5"),
  
  # Mitochondrial function
  Mitochondrial = c("PPARGC1A", "PPARGC1B", "TFAM", "NRF1", "ESRRA")
)

# Uppercase all markers for consistency
lipid_markers <- lapply(lipid_markers, toupper)

# Flatten all lipid markers
all_lipid_markers <- unique(unlist(lipid_markers))
cat("Total unique lipid markers defined:", length(all_lipid_markers), "\n")

# ============================================================================
# PART 3: CHECK MARKER AVAILABILITY IN PROTEOMICS DATA
# ============================================================================

cat("\n=== CHECKING MARKER AVAILABILITY ===\n")
cat("NOTE: Cytokines (IL10, TGFB1, IFNG, TNF) and enzymes (NOS2, PTGS2)\n")
cat("      are often undetected in discovery proteomics - this is expected.\n\n")

# Function to check marker availability
check_markers <- function(marker_list, protein_names) {
  results <- list()
  for (category in names(marker_list)) {
    markers <- marker_list[[category]]
    available <- markers[markers %in% protein_names]
    missing <- markers[!markers %in% protein_names]
    
    results[[category]] <- list(
      available = available,
      missing = missing,
      n_available = length(available),
      n_total = length(markers),
      pct_available = round(100 * length(available) / length(markers), 1)
    )
    
    cat(sprintf("\n%s: %d/%d (%.1f%%) available\n", 
                category, length(available), length(markers),
                100 * length(available) / length(markers)))
    if (length(missing) > 0 && length(missing) <= 5) {
      cat("  Missing:", paste(missing, collapse=", "), "\n")
    }
  }
  return(results)
}

# Check immune markers
cat("\n--- IMMUNE MARKERS ---")
immune_available <- check_markers(immune_markers, rownames(protein_matrix))

# Check lipid markers
cat("\n--- LIPID MARKERS ---")
lipid_available <- check_markers(lipid_markers, rownames(protein_matrix))

# Get all available markers for correlation analysis
immune_markers_detected <- unique(unlist(lapply(immune_available, function(x) x$available)))
lipid_markers_detected <- unique(unlist(lapply(lipid_available, function(x) x$available)))

cat("\n=== SUMMARY ===\n")
cat("Total immune markers detected:", length(immune_markers_detected), "\n")
cat("Total lipid markers detected:", length(lipid_markers_detected), "\n")

# Save list of tested genes for provenance
fwrite(data.frame(Gene = immune_markers_detected, Category = "Immune"),
       "DGAT1_Analysis_Immune_Genes_Tested.csv")
fwrite(data.frame(Gene = lipid_markers_detected, Category = "Lipid"),
       "DGAT1_Analysis_Lipid_Genes_Tested.csv")
cat("Gene lists saved for reproducibility.\n")

# ============================================================================
# PART 4: DGAT1 CORRELATION WITH IMMUNE MARKERS
# ============================================================================

cat("\n=== ANALYZING DGAT1-IMMUNE CORRELATIONS ===\n")

# Function to calculate correlations
calc_correlations <- function(target_gene, marker_genes, protein_mat) {
  target_expr <- as.numeric(protein_mat[target_gene, ])
  
  results <- data.frame(
    Gene = character(),
    Correlation = numeric(),
    P_value = numeric(),
    N_samples = integer(),
    stringsAsFactors = FALSE
  )
  
  for (gene in marker_genes) {
    if (gene %in% rownames(protein_mat)) {
      marker_expr <- as.numeric(protein_mat[gene, ])
      
      # Remove NA pairs
      valid_idx <- !is.na(target_expr) & !is.na(marker_expr)
      
      if (sum(valid_idx) >= 10) {  # Minimum 10 samples
        cor_test <- cor.test(target_expr[valid_idx], 
                            marker_expr[valid_idx], 
                            method = "spearman")
        
        results <- rbind(results, data.frame(
          Gene = gene,
          Correlation = cor_test$estimate,
          P_value = cor_test$p.value,
          N_samples = sum(valid_idx),
          stringsAsFactors = FALSE
        ))
      }
    }
  }
  
  # Adjust p-values for multiple testing
  if (nrow(results) > 0) {
    results$P_adjusted <- p.adjust(results$P_value, method = "BH")
    results$Significant <- results$P_adjusted < 0.05
    results <- results[order(abs(results$Correlation), decreasing = TRUE), ]
  }
  
  return(results)
}

# Calculate DGAT1-immune correlations
dgat1_immune_cor <- calc_correlations("DGAT1", immune_markers_detected, protein_matrix)

cat("\nDGAT1-Immune Marker Correlations (Top 20):\n")
print(head(dgat1_immune_cor, 20))

# Save results
fwrite(dgat1_immune_cor, "DGAT1_Immune_Correlations.csv")

# ============================================================================
# PART 5: DGAT1 CORRELATION WITH LIPID MARKERS
# ============================================================================

cat("\n=== ANALYZING DGAT1-LIPID CORRELATIONS ===\n")

# Calculate DGAT1-lipid correlations
dgat1_lipid_cor <- calc_correlations("DGAT1", lipid_markers_detected, protein_matrix)

cat("\nDGAT1-Lipid Marker Correlations (Top 20):\n")
print(head(dgat1_lipid_cor, 20))

# Save results
fwrite(dgat1_lipid_cor, "DGAT1_Lipid_Correlations.csv")

# ============================================================================
# PART 6: VISUALIZATIONS
# ============================================================================

cat("\n=== GENERATING VISUALIZATIONS ===\n")

# Create output directory for plots
plot_dir <- "DGAT1_Analysis_Plots"
dir.create(plot_dir, showWarnings = FALSE)

# ---------------------
# 6.1: Correlation Bar Plots
# ---------------------

# Function to create correlation bar plot
plot_correlation_bars <- function(cor_data, title, top_n = 30) {
  # Get significant correlations only
  cor_data_sig <- cor_data[cor_data$Significant, ]
  
  if (nrow(cor_data_sig) == 0) {
    warning("No significant correlations to plot for: ", title)
    return(NULL)
  }
  
  # Balance top positive and negative
  pos <- cor_data_sig[order(cor_data_sig$Correlation, decreasing = TRUE), ]
  neg <- cor_data_sig[order(cor_data_sig$Correlation, decreasing = FALSE), ]
  need_each <- max(1, floor(top_n/2))
  pos_top <- head(pos, min(need_each, nrow(pos)))
  neg_top <- head(neg, min(need_each, nrow(neg)))
  cor_data_plot <- unique(rbind(pos_top, neg_top))
  
  cor_data_plot$Gene <- factor(cor_data_plot$Gene, 
                               levels = cor_data_plot$Gene[order(cor_data_plot$Correlation)])
  cor_data_plot$Direction <- ifelse(cor_data_plot$Correlation > 0, "Positive", "Negative")
  
  p <- ggplot(cor_data_plot, aes(x = Gene, y = Correlation, fill = Direction)) +
    geom_bar(stat = "identity") +
    scale_fill_manual(values = c("Positive" = "#D73027", "Negative" = "#4575B4")) +
    coord_flip() +
    labs(title = title,
         subtitle = sprintf("Spearman correlation (FDR < 0.05, n=%d)", nrow(cor_data_plot)),
         x = "", y = "Spearman Correlation") +
    theme_minimal() +
    theme(
      legend.position = "bottom",
      plot.title = element_text(face = "bold", size = 14),
      axis.text.y = element_text(size = 10)
    ) +
    geom_hline(yintercept = 0, linetype = "dashed")
  
  return(p)
}

# Plot immune correlations
p1 <- plot_correlation_bars(dgat1_immune_cor, 
                            "DGAT1 vs Immune Markers", 
                            top_n = 30)
if (!is.null(p1)) {
  ggsave(file.path(plot_dir, "01_DGAT1_Immune_Correlations_Barplot.pdf"), 
         p1, width = 10, height = 12)
}

# Plot lipid correlations
p2 <- plot_correlation_bars(dgat1_lipid_cor, 
                            "DGAT1 vs Lipid Metabolism Markers", 
                            top_n = 30)
if (!is.null(p2)) {
  ggsave(file.path(plot_dir, "02_DGAT1_Lipid_Correlations_Barplot.pdf"), 
         p2, width = 10, height = 12)
}

# ---------------------
# 6.2: Volcano Plots
# ---------------------

# Function to create volcano plot
plot_volcano <- function(cor_data, title) {
  # Add p-value floor to prevent -log10(0) = Inf
  eps <- 1e-300
  cor_data$P_adjusted[is.na(cor_data$P_adjusted)] <- 1
  cor_data$NegLog10P <- -log10(pmax(cor_data$P_adjusted, eps))
  
  cor_data$Label <- ifelse(cor_data$Significant & abs(cor_data$Correlation) > 0.3,
                           cor_data$Gene, "")
  
  p <- ggplot(cor_data, aes(x = Correlation, y = NegLog10P)) +
    geom_point(aes(color = Significant), alpha = 0.6, size = 2) +
    scale_color_manual(values = c("FALSE" = "gray70", "TRUE" = "#D73027")) +
    geom_vline(xintercept = c(-0.3, 0.3), linetype = "dashed", color = "gray40") +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "gray40") +
    ggrepel::geom_text_repel(aes(label = Label), size = 3, max.overlaps = 20) +
    labs(title = title,
         x = "Spearman Correlation with DGAT1",
         y = "-log10(Adjusted P-value)") +
    theme_minimal() +
    theme(
      legend.position = "right",
      plot.title = element_text(face = "bold", size = 14)
    )
  
  return(p)
}

# Volcano plot for immune markers
p3 <- plot_volcano(dgat1_immune_cor, "DGAT1 vs Immune Markers - Volcano Plot")
ggsave(file.path(plot_dir, "03_DGAT1_Immune_Volcano.pdf"), 
       p3, width = 12, height = 8)

# Volcano plot for lipid markers
p4 <- plot_volcano(dgat1_lipid_cor, "DGAT1 vs Lipid Markers - Volcano Plot")
ggsave(file.path(plot_dir, "04_DGAT1_Lipid_Volcano.pdf"), 
       p4, width = 12, height = 8)

# ---------------------
# 6.3: Correlation Heatmaps
# ---------------------

# Function to create correlation heatmap
create_correlation_heatmap <- function(target_gene, marker_genes, protein_mat, title) {
  # Filter available markers (excluding target gene)
  available_markers <- marker_genes[marker_genes %in% rownames(protein_mat)]
  available_markers <- setdiff(available_markers, target_gene)
  
  if (length(available_markers) < 2) {
    message("Not enough markers available for heatmap: ", title)
    return(NULL)
  }
  
  # Get expression matrix
  genes_for_heatmap <- c(target_gene, available_markers)
  expr_subset <- protein_mat[genes_for_heatmap, , drop=FALSE]
  
  # Calculate correlation matrix
  cor_matrix <- suppressWarnings(
    cor(t(expr_subset), use = "pairwise.complete.obs", method = "spearman")
  )
  
  # Plot heatmap
  pheatmap(cor_matrix,
           color = colorRampPalette(c("#4575B4", "white", "#D73027"))(100),
           breaks = seq(-1, 1, length.out = 101),
           main = title,
           fontsize = 8,
           fontsize_row = 7,
           fontsize_col = 7,
           filename = file.path(plot_dir, paste0(gsub(" ", "_", title), ".pdf")),
           width = 12,
           height = 11)
  
  return(cor_matrix)
}

# Heatmap for top immune markers
top_immune <- head(dgat1_immune_cor[dgat1_immune_cor$Significant, ]$Gene, 40)
if (length(top_immune) > 0) {
  create_correlation_heatmap("DGAT1", top_immune, protein_matrix,
                            "05_DGAT1_Top_Immune_Markers_Heatmap")
}

# Heatmap for top lipid markers
top_lipid <- head(dgat1_lipid_cor[dgat1_lipid_cor$Significant, ]$Gene, 40)
if (length(top_lipid) > 0) {
  create_correlation_heatmap("DGAT1", top_lipid, protein_matrix,
                            "06_DGAT1_Top_Lipid_Markers_Heatmap")
}

# ---------------------
# 6.4: Category-wise Box Plots
# ---------------------

# Function to create category-wise correlation summary
plot_category_summary <- function(cor_data, marker_list, title) {
  # Assign categories
  cor_data$Category <- NA
  for (cat_name in names(marker_list)) {
    cor_data$Category[cor_data$Gene %in% marker_list[[cat_name]]] <- cat_name
  }
  
  cor_data_cat <- cor_data[!is.na(cor_data$Category), ]
  
  if (nrow(cor_data_cat) == 0) return(NULL)
  
  p <- ggplot(cor_data_cat, aes(x = Category, y = Correlation, fill = Category)) +
    geom_boxplot(alpha = 0.7, outlier.shape = NA) +
    geom_jitter(width = 0.2, alpha = 0.5, size = 2) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    scale_fill_brewer(palette = "Set3") +
    labs(title = title,
         subtitle = "Distribution of DGAT1 correlations by category",
         x = "", y = "Spearman Correlation") +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
      legend.position = "none",
      plot.title = element_text(face = "bold", size = 14)
    )
  
  return(p)
}

# Category summary for immune markers
p5 <- plot_category_summary(dgat1_immune_cor, immune_markers, 
                            "DGAT1 Correlations by Immune Category")
if (!is.null(p5)) {
  ggsave(file.path(plot_dir, "07_DGAT1_Immune_Category_Summary.pdf"), 
         p5, width = 12, height = 8)
}

# Category summary for lipid markers
p6 <- plot_category_summary(dgat1_lipid_cor, lipid_markers, 
                            "DGAT1 Correlations by Lipid Metabolism Category")
if (!is.null(p6)) {
  ggsave(file.path(plot_dir, "08_DGAT1_Lipid_Category_Summary.pdf"), 
         p6, width = 12, height = 8)
}

# ============================================================================
# PART 7: KEY FINDINGS SUMMARY
# ============================================================================

sep70 <- paste(rep("=", 70), collapse = "")

cat("\n", sep70, "\n")
cat("=== KEY FINDINGS SUMMARY ===\n")
cat(sep70, "\n\n")

# Immune findings
cat("IMMUNE MARKERS:\n")
cat("---------------\n")
dgat1_immune_sig <- dgat1_immune_cor[dgat1_immune_cor$Significant, ]
cat("Significant correlations:", nrow(dgat1_immune_sig), "\n")
cat("Positive correlations:", sum(dgat1_immune_sig$Correlation > 0), "\n")
cat("Negative correlations:", sum(dgat1_immune_sig$Correlation < 0), "\n\n")

cat("Top 5 POSITIVE immune correlations:\n")
print(head(dgat1_immune_sig[dgat1_immune_sig$Correlation > 0, c("Gene", "Correlation", "P_adjusted")], 5))

cat("\nTop 5 NEGATIVE immune correlations:\n")
print(head(dgat1_immune_sig[dgat1_immune_sig$Correlation < 0, c("Gene", "Correlation", "P_adjusted")], 5))

# Lipid findings
cat("\n\nLIPID MARKERS:\n")
cat("---------------\n")
dgat1_lipid_sig <- dgat1_lipid_cor[dgat1_lipid_cor$Significant, ]
cat("Significant correlations:", nrow(dgat1_lipid_sig), "\n")
cat("Positive correlations:", sum(dgat1_lipid_sig$Correlation > 0), "\n")
cat("Negative correlations:", sum(dgat1_lipid_sig$Correlation < 0), "\n\n")

cat("Top 5 POSITIVE lipid correlations:\n")
print(head(dgat1_lipid_sig[dgat1_lipid_sig$Correlation > 0, c("Gene", "Correlation", "P_adjusted")], 5))

cat("\nTop 5 NEGATIVE lipid correlations:\n")
print(head(dgat1_lipid_sig[dgat1_lipid_sig$Correlation < 0, c("Gene", "Correlation", "P_adjusted")], 5))

# ============================================================================
# PART 8: GENERATE FINAL SUMMARY REPORT
# ============================================================================

# Create summary report
summary_report <- list(
  Analysis_Date = Sys.Date(),
  Dataset = "CPTAC3 GBM Discovery Cohort",
  N_Samples = ncol(protein_matrix),
  N_Proteins = nrow(protein_matrix),
  DGAT1_Detected = "DGAT1" %in% rownames(protein_matrix),
  DGAT1_Non_NA = sum(!is.na(dgat1_values)),
  
  # Log tested genes
  Immune_Genes_Tested = immune_markers_detected,
  Lipid_Genes_Tested = lipid_markers_detected,
  
  Immune_Analysis = list(
    N_Markers_Tested = length(immune_markers_detected),
    N_Significant = nrow(dgat1_immune_sig),
    N_Positive_Sig = sum(dgat1_immune_sig$Correlation > 0),
    N_Negative_Sig = sum(dgat1_immune_sig$Correlation < 0),
    Top_Positive = head(dgat1_immune_sig[order(-dgat1_immune_sig$Correlation), "Gene"], 5),
    Top_Negative = head(dgat1_immune_sig[order(dgat1_immune_sig$Correlation), "Gene"], 5)
  ),
  
  Lipid_Analysis = list(
    N_Markers_Tested = length(lipid_markers_detected),
    N_Significant = nrow(dgat1_lipid_sig),
    N_Positive_Sig = sum(dgat1_lipid_sig$Correlation > 0),
    N_Negative_Sig = sum(dgat1_lipid_sig$Correlation < 0),
    Top_Positive = head(dgat1_lipid_sig[order(-dgat1_lipid_sig$Correlation), "Gene"], 5),
    Top_Negative = head(dgat1_lipid_sig[order(dgat1_lipid_sig$Correlation), "Gene"], 5)
  )
)

# Save summary
saveRDS(summary_report, "DGAT1_Analysis_Summary.rds")

# Save session info for reproducibility
writeLines(capture.output(sessionInfo()), con = "DGAT1_Analysis_sessionInfo.txt")

cat("\n✅ Analysis complete! All results saved.\n")
cat("\nOutput directory:", OUTPUT_DIR, "\n")
cat("\nOutput files:\n")
cat("  - DGAT1_Immune_Correlations.csv\n")
cat("  - DGAT1_Lipid_Correlations.csv\n")
cat("  - DGAT1_Analysis_Summary.rds\n")
cat("  - DGAT1_Analysis_sessionInfo.txt\n")
cat("  - DGAT1_Analysis_Immune_Genes_Tested.csv\n")
cat("  - DGAT1_Analysis_Lipid_Genes_Tested.csv\n")
cat("  - DGAT1_Analysis_Plots/ (8 visualizations)\n")

cat("\n", sep70, "\n")
cat("END OF ANALYSIS\n")
cat(sep70, "\n")


