#!/usr/bin/env Rscript
# -----------------------------------------------------------------------------
# Fixed TCGA Analysis Script for DGAT Immunology
# Updated to work with actual clean TCGA-GBM data structure
# -----------------------------------------------------------------------------

# Load required libraries
suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(survival)
  library(survminer)
  library(readr)
  library(data.table)
  library(pheatmap)
  library(RColorBrewer)
})

# Set working directory
setwd(".")

# Configuration - Updated for actual data structure
config <- list(
  data_dir = "Processed_Data/TCGA_GBM_Clean",
  results_dir = "Results/Legacy_Analysis",
  cancer_types = c("GBM"),  # Only GBM available
  dgat_genes = c("DGAT1", "DGAT2")
)

# Create results directory if it doesn't exist
if (!dir.exists(config$results_dir)) {
  dir.create(config$results_dir, recursive = TRUE)
}

# Function to load TCGA data - Fixed for actual file structure
load_tcga_data <- function(cancer_type) {
  cat("Loading TCGA data for", cancer_type, "...\n")
  
  # Load expression data - Updated file path and format
  expression_file <- file.path(config$data_dir, "TCGA_GBM_Expression_Cleaned.rds")
  
  if (!file.exists(expression_file)) {
    cat("Expression file not found:", expression_file, "\n")
    return(NULL)
  }
  
  expression_data <- readRDS(expression_file)
  cat("Expression data loaded:", nrow(expression_data), "genes x", ncol(expression_data), "samples\n")
  
  # Load clinical data - Updated file path and format
  clinical_file <- file.path(config$data_dir, "TCGA_GBM_Clinical_Cleaned.csv")
  
  if (!file.exists(clinical_file)) {
    cat("Clinical file not found:", clinical_file, "\n")
    return(NULL)
  }
  
  clinical_data <- fread(clinical_file)
  cat("Clinical data loaded:", nrow(clinical_data), "patients\n")
  
  return(list(
    expression = expression_data,
    clinical = clinical_data
  ))
}

# Function to extract DGAT expression - Fixed for matrix format
extract_dgat_expression <- function(expression_data) {
  cat("Extracting DGAT expression...\n")
  
  # Check if DGAT genes exist
  dgat_genes_found <- config$dgat_genes[config$dgat_genes %in% rownames(expression_data)]
  
  if (length(dgat_genes_found) == 0) {
    cat("No DGAT genes found in expression data\n")
    cat("Available genes:", length(rownames(expression_data)), "total\n")
    cat("First 10 genes:", paste(head(rownames(expression_data), 10), collapse = ", "), "\n")
    return(NULL)
  }
  
  cat("Found DGAT genes:", paste(dgat_genes_found, collapse = ", "), "\n")
  
  # Extract DGAT expression data
  dgat_expr <- expression_data[dgat_genes_found, , drop = FALSE]
  
  # Convert to long format
  dgat_long <- data.frame(
    sample_id = rep(colnames(dgat_expr), each = nrow(dgat_expr)),
    gene_name = rep(rownames(dgat_expr), ncol(dgat_expr)),
    expression = as.vector(dgat_expr),
    stringsAsFactors = FALSE
  )
  
  cat("DGAT expression extracted for", nrow(dgat_long), "gene-sample combinations\n")
  
  return(dgat_long)
}

# Function to perform survival analysis - Fixed column names and logic
perform_survival_analysis <- function(dgat_data, clinical_data, cancer_type) {
  cat("Performing survival analysis for", cancer_type, "...\n")
  
  # Merge DGAT expression with clinical data
  survival_data <- clinical_data %>%
    left_join(dgat_data, by = c("barcode" = "sample_id"))
  
  # Remove samples with missing survival data - Fixed column names
  survival_data <- survival_data %>%
    filter(!is.na(OS_years) & !is.na(OS_status))
  
  if (nrow(survival_data) == 0) {
    cat("No samples with survival data found\n")
    return(NULL)
  }
  
  cat("Survival data prepared:", nrow(survival_data), "samples\n")
  
  results <- list()
  
  # Analyze each DGAT gene
  for (gene in config$dgat_genes) {
    cat("Analyzing", gene, "...\n")
    
    # Filter data for current gene
    gene_data <- survival_data %>%
      filter(gene_name == gene) %>%
      filter(!is.na(expression))
    
    if (nrow(gene_data) == 0) {
      cat("No expression data for", gene, "\n")
      next
    }
    
    cat("  Samples for", gene, ":", nrow(gene_data), "\n")
    
    # Create survival object - Fixed formula
    surv_obj <- Surv(time = gene_data$OS_years, event = gene_data$OS_status)
    
    # Cox proportional hazards model - Fixed formula
    cox_formula <- Surv(OS_years, OS_status) ~ expression
    cox_model <- coxph(cox_formula, data = gene_data)
    
    # Kaplan-Meier analysis
    gene_data$expression_group <- ifelse(gene_data$expression > median(gene_data$expression, na.rm = TRUE),
                                       "High", "Low")
    
    km_fit <- survfit(Surv(OS_years, OS_status) ~ expression_group, data = gene_data)
    
    # Store results
    results[[gene]] <- list(
      cox_model = cox_model,
      km_fit = km_fit,
      data = gene_data,
      n_samples = nrow(gene_data),
      median_expr = median(gene_data$expression, na.rm = TRUE),
      hr = exp(coef(cox_model)),
      p_value = summary(cox_model)$coefficients["expression", "Pr(>|z|)"]
    )
    
    cat("  ", gene, ": HR =", round(results[[gene]]$hr, 3), 
        ", p =", format(results[[gene]]$p_value, digits = 3), "\n")
  }
  
  return(results)
}

# Function to create survival plots - Fixed for new data structure
create_survival_plots <- function(survival_results, cancer_type) {
  cat("Creating survival plots for", cancer_type, "...\n")
  
  plots <- list()
  
  for (gene in names(survival_results)) {
    if (is.null(survival_results[[gene]])) next
    
    # Create Kaplan-Meier plot
    km_plot <- ggsurvplot(
      survival_results[[gene]]$km_fit,
      data = survival_results[[gene]]$data,
      pval = TRUE,
      conf.int = TRUE,
      title = paste(gene, "Expression and Survival -", cancer_type),
      subtitle = paste("HR =", round(survival_results[[gene]]$hr, 3),
                      ", p =", format(survival_results[[gene]]$p_value, digits = 3)),
      xlab = "Time (years)",
      ylab = "Overall Survival Probability",
      legend.labs = c("High Expression", "Low Expression"),
      palette = c("#E7B800", "#2E9FDF"),
      risk.table = TRUE,
      risk.table.height = 0.3
    )
    
    plots[[gene]] <- km_plot
    
    # Save plot
    plot_file <- file.path(config$results_dir, 
                          paste0(gene, "_survival_", cancer_type, ".png"))
    ggsave(plot_file, km_plot$plot, width = 10, height = 8, dpi = 300)
    cat("  Saved plot for", gene, "\n")
  }
  
  return(plots)
}

# Function to perform correlation analysis - Simplified without immune signatures
perform_correlation_analysis <- function(expression_data, cancer_type) {
  cat("Performing correlation analysis for", cancer_type, "...\n")
  
  # Check if DGAT genes exist
  dgat_genes_found <- config$dgat_genes[config$dgat_genes %in% rownames(expression_data)]
  
  if (length(dgat_genes_found) == 0) {
    cat("No DGAT genes found for correlation analysis\n")
    return(NULL)
  }
  
  # Extract DGAT expression
  dgat_expr <- expression_data[dgat_genes_found, , drop = FALSE]
  
  # Calculate correlation between DGAT genes
  if (nrow(dgat_expr) >= 2) {
    dgat_cor <- cor(t(dgat_expr))
    cat("DGAT gene correlation matrix:\n")
    print(round(dgat_cor, 3))
    
    # Create correlation heatmap
    if (nrow(dgat_cor) > 1) {
      cor_plot <- pheatmap(
        dgat_cor,
        cluster_rows = FALSE,
        cluster_cols = FALSE,
        display_numbers = TRUE,
        number_format = "%.3f",
        main = "DGAT Gene Expression Correlation",
        color = colorRampPalette(c("#2E9FDF", "white", "#E7B800"))(100)
      )
      
      # Save correlation plot
      cor_file <- file.path(config$results_dir, "DGAT_correlation_heatmap.png")
      png(cor_file, width = 600, height = 500, res = 300)
      print(cor_plot)
      dev.off()
      cat("  Saved correlation heatmap\n")
    }
    
    return(dgat_cor)
  } else {
    cat("Need at least 2 DGAT genes for correlation analysis\n")
    return(NULL)
  }
}

# Function to create summary report
create_summary_report <- function(survival_results, cancer_type) {
  cat("Creating summary report for", cancer_type, "...\n")
  
  # Create summary data frame
  summary_data <- data.frame(
    Gene = character(),
    N_Samples = integer(),
    Median_Expression = numeric(),
    Hazard_Ratio = numeric(),
    P_Value = numeric(),
    Significant = logical(),
    stringsAsFactors = FALSE
  )
  
  for (gene in names(survival_results)) {
    if (!is.null(survival_results[[gene]])) {
      summary_data <- rbind(summary_data, data.frame(
        Gene = gene,
        N_Samples = survival_results[[gene]]$n_samples,
        Median_Expression = survival_results[[gene]]$median_expr,
        Hazard_Ratio = survival_results[[gene]]$hr,
        P_Value = survival_results[[gene]]$p_value,
        Significant = survival_results[[gene]]$p_value < 0.05,
        stringsAsFactors = FALSE
      ))
    }
  }
  
  # Save summary
  summary_file <- file.path(config$results_dir, paste0("DGAT_summary_", cancer_type, ".csv"))
  fwrite(summary_data, summary_file)
  cat("  Saved summary report:", summary_file, "\n")
  
  # Print summary
  cat("\n=== DGAT SURVIVAL ANALYSIS SUMMARY ===\n")
  print(summary_data)
  
  return(summary_data)
}

# Main analysis function
main_analysis <- function() {
  cat("============================================================================\n")
  cat(" DGAT TCGA ANALYSIS - FIXED VERSION\n")
  cat("============================================================================\n")
  
  # Load data
  tcga_data <- load_tcga_data("GBM")
  if (is.null(tcga_data)) {
    cat("❌ Failed to load TCGA data\n")
    return(FALSE)
  }
  
  # Extract DGAT expression
  dgat_data <- extract_dgat_expression(tcga_data$expression)
  if (is.null(dgat_data)) {
    cat("❌ Failed to extract DGAT expression\n")
    return(FALSE)
  }
  
  # Perform survival analysis
  survival_results <- perform_survival_analysis(dgat_data, tcga_data$clinical, "GBM")
  if (is.null(survival_results) || length(survival_results) == 0) {
    cat("❌ Failed to perform survival analysis\n")
    return(FALSE)
  }
  
  # Create survival plots
  plots <- create_survival_plots(survival_results, "GBM")
  
  # Perform correlation analysis
  cor_results <- perform_correlation_analysis(tcga_data$expression, "GBM")
  
  # Create summary report
  summary_data <- create_summary_report(survival_results, "GBM")
  
  cat("\n============================================================================\n")
  cat(" ANALYSIS COMPLETE\n")
  cat("============================================================================\n")
  cat("Results saved to:", config$results_dir, "\n")
  
  # Final summary
  n_significant <- sum(summary_data$Significant, na.rm = TRUE)
  cat(sprintf("Significant associations: %d/%d DGAT genes\n", n_significant, nrow(summary_data)))
  
  if (n_significant > 0) {
    cat("✓ DGAT genes show prognostic significance in GBM\n")
  } else {
    cat("⚠ No significant prognostic associations found\n")
  }
  
  return(TRUE)
}

# Run analysis
if (!interactive()) {
  success <- main_analysis()
  quit(status = ifelse(success, 0, 1))
}
