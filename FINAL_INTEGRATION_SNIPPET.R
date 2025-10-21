# Final Integration Snippet for DGAT Immunology Analysis
# Combines results from multiple analysis components

library(dplyr)
library(ggplot2)
library(readr)

# Load all analysis results
load_integration_results <- function() {
  
  # Load survival analysis results
  survival_results <- list()
  cancer_types <- c("BRCA", "GBM", "OV", "PAAD")
  
  for (cancer in cancer_types) {
    survival_file <- paste0("Results/Bulk/Survival/TCGA_", cancer, "/survival_summary.csv")
    if (file.exists(survival_file)) {
      survival_results[[cancer]] <- read_csv(survival_file)
    }
  }
  
  # Load correlation analysis results
  correlation_results <- list()
  for (cancer in cancer_types) {
    cor_file <- paste0("Results/Bulk/", cancer, "/Top20_DGAT_Immune_Correlations.csv")
    if (file.exists(cor_file)) {
      correlation_results[[cancer]] <- read_csv(cor_file)
    }
  }
  
  # Load protein analysis results
  protein_results <- NULL
  protein_file <- "Results/Proteome/Reports/DGAT_protein_analysis_report.csv"
  if (file.exists(protein_file)) {
    protein_results <- read_csv(protein_file)
  }
  
  return(list(
    survival = survival_results,
    correlation = correlation_results,
    protein = protein_results
  ))
}

# Create integrated summary
create_integrated_summary <- function(results) {
  
  # Combine survival results across cancer types
  all_survival <- do.call(rbind, results$survival)
  
  # Summary statistics
  cat("=== INTEGRATED ANALYSIS SUMMARY ===\n")
  cat("Total survival analyses:", nrow(all_survival), "\n")
  cat("Cancer types analyzed:", length(unique(all_survival$Cancer_Type)), "\n")
  cat("Genes analyzed:", length(unique(all_survival$Gene)), "\n")
  
  # Significant associations
  significant <- all_survival[all_survival$P_Value < 0.05, ]
  cat("Significant associations (p < 0.05):", nrow(significant), "\n")
  
  # Top correlations across cancer types
  all_correlations <- do.call(rbind, results$correlation)
  top_correlations <- all_correlations %>%
    arrange(desc(abs(correlation))) %>%
    head(10)
  
  cat("\nTop 10 DGAT-Immune Correlations:\n")
  print(top_correlations)
  
  return(list(
    survival_summary = all_survival,
    top_correlations = top_correlations,
    significant_associations = significant
  ))
}

# Main integration function
run_integration <- function() {
  cat("Loading integration results...\n")
  results <- load_integration_results()
  
  cat("Creating integrated summary...\n")
  summary <- create_integrated_summary(results)
  
  # Save integrated results
  write_csv(summary$survival_summary, "Results/integrated_survival_summary.csv")
  write_csv(summary$top_correlations, "Results/integrated_top_correlations.csv")
  
  cat("Integration complete! Results saved to Results/ directory.\n")
  
  return(summary)
}

# Run integration if script is executed directly
if (!interactive()) {
  integration_results <- run_integration()
}