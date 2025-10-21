#!/usr/bin/env Rscript
# =============================================================================
# TCGA GBM: DGAT1 Survival Analysis
# Overall survival with optimal cutoff selection
# =============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(survival)
  library(survminer)
  library(ggrepel)
})

# Optional packages
if (requireNamespace("cowplot", quietly = TRUE)) {
  library(cowplot)
  COWPLOT_AVAILABLE <- TRUE
} else {
  COWPLOT_AVAILABLE <- FALSE
  cat("cowplot not available - will use simple plot layout\n")
}

# Optional packages
if (requireNamespace("maxstat", quietly = TRUE)) {
  library(maxstat)
  MAXSTAT_AVAILABLE <- TRUE
} else {
  MAXSTAT_AVAILABLE <- FALSE
  cat("maxstat package not available - will use median cutoff\n")
}

# ====== CONFIG =================================================================
BASE_DIR <- "/Users/zijiefeng/Desktop/Guo's lab/My_Research/DGAT_Immunology"
DATA_DIR <- file.path(BASE_DIR, "Processed_Data", "TCGA_GBM_Clean")
OUT_DIR  <- file.path(BASE_DIR, "Results", "TCGA_Survival")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

# Input files
EXPRESSION_FILE <- file.path(DATA_DIR, "TCGA_GBM_Expression_Cleaned.rds")
CLINICAL_FILE <- file.path(DATA_DIR, "TCGA_GBM_Clinical_Cleaned.csv")

# =============================================================================
# FUNCTIONS
# =============================================================================

load_expression_data <- function(filepath) {
  cat("Loading expression data...\n")
  if (grepl("\\.rds$", filepath)) {
    expr <- readRDS(filepath)
  } else if (grepl("\\.csv$", filepath)) {
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
    expr_df <- as.data.frame(expr) %>%
      mutate(Gene = rownames(expr)) %>%
      group_by(Gene) %>%
      summarise(across(everything(), ~median(.x, na.rm = TRUE)), .groups = "drop")
    expr <- as.matrix(expr_df[, -1])
    rownames(expr) <- expr_df$Gene
  }
  
  cat("  Expression matrix:", nrow(expr), "genes x", ncol(expr), "samples\n")
  return(expr)
}

load_clinical_data <- function(filepath) {
  cat("\nLoading clinical data...\n")
  clinical <- read_csv(filepath, show_col_types = FALSE)
  
  # Standardize column names
  names(clinical) <- tolower(names(clinical))
  
  # Try to identify key columns
  sample_col <- names(clinical)[grep("sample|patient|case|barcode", names(clinical), ignore.case = TRUE)][1]
  os_time_col <- names(clinical)[grep("os\\.time|survival|days|time", names(clinical), ignore.case = TRUE)][1]
  os_status_col <- names(clinical)[grep("os\\.status|vital|status|death", names(clinical), ignore.case = TRUE)][1]
  
  if (is.na(sample_col) | is.na(os_time_col) | is.na(os_status_col)) {
    cat("\n  Available columns:", paste(names(clinical), collapse = ", "), "\n")
    stop("Could not identify required clinical columns. Please check column names.")
  }
  
  # Rename for standardization
  clinical <- clinical %>%
    rename(
      sample_id = !!sym(sample_col),
      os_time = !!sym(os_time_col),
      os_status = !!sym(os_status_col)
    )
  
  # Convert OS time to months if in days
  if (max(clinical$os_time, na.rm = TRUE) > 365) {
    clinical$os_time <- clinical$os_time / 30.44  # Convert days to months
    cat("  Converted OS time from days to months\n")
  }
  
  # Standardize OS status (1 = death, 0 = censored)
  if (is.character(clinical$os_status)) {
    clinical$os_status <- ifelse(
      tolower(clinical$os_status) %in% c("dead", "deceased", "1", "yes"), 1, 0
    )
  }
  
  # Extract covariates if available
  age_col <- names(clinical)[grep("age", names(clinical), ignore.case = TRUE)][1]
  sex_col <- names(clinical)[grep("sex|gender", names(clinical), ignore.case = TRUE)][1]
  
  if (!is.na(age_col)) clinical$age <- clinical[[age_col]]
  if (!is.na(sex_col)) clinical$sex <- clinical[[sex_col]]
  
  cat("  Clinical data:", nrow(clinical), "samples\n")
  cat("  Events (deaths):", sum(clinical$os_status == 1, na.rm = TRUE), "\n")
  cat("  Median follow-up:", round(median(clinical$os_time, na.rm = TRUE), 1), "months\n")
  
  return(clinical)
}

prepare_survival_data <- function(expr_mat, clinical, gene = "DGAT1") {
  cat("\nPreparing survival data for", gene, "...\n")
  
  if (!gene %in% rownames(expr_mat)) {
    stop(paste(gene, "not found in expression data"))
  }
  
  dgat1_expr <- data.frame(
    sample_id = colnames(expr_mat),
    dgat1 = as.numeric(expr_mat[gene, ]),
    stringsAsFactors = FALSE
  )
  
  # Merge with clinical
  surv_data <- clinical %>%
    inner_join(dgat1_expr, by = "sample_id") %>%
    filter(!is.na(os_time) & !is.na(os_status) & !is.na(dgat1))
  
  cat("  Merged data:", nrow(surv_data), "samples with complete survival + expression\n")
  
  return(surv_data)
}

find_optimal_cutoff <- function(surv_data) {
  cat("\nFinding optimal cutoff...\n")
  
  if (MAXSTAT_AVAILABLE) {
    # Use maximally selected rank statistics
    tryCatch({
      maxstat_test <- maxstat.test(
        Surv(os_time, os_status) ~ dgat1,
        data = surv_data,
        smethod = "LogRank",
        pmethod = "condMC",
        B = 9999
      )
      
      optimal_cutoff <- maxstat_test$estimate
      maxstat_p <- maxstat_test$p.value
      
      cat("  Optimal cutoff (maxstat):", round(optimal_cutoff, 3), "\n")
      cat("  Maxstat p-value:", format.pval(maxstat_p, digits = 3), "\n")
      
      return(list(cutoff = optimal_cutoff, p_value = maxstat_p, method = "maxstat"))
    }, error = function(e) {
      cat("  Maxstat failed, using median cutoff\n")
      median_cutoff <- median(surv_data$dgat1, na.rm = TRUE)
      return(list(cutoff = median_cutoff, p_value = NA, method = "median"))
    })
  } else {
    cat("  Using median cutoff (maxstat not available)\n")
    median_cutoff <- median(surv_data$dgat1, na.rm = TRUE)
    return(list(cutoff = median_cutoff, p_value = NA, method = "median"))
  }
}

perform_survival_analysis <- function(surv_data, cutoff_result) {
  cat("\nPerforming survival analysis...\n")
  
  # Create high/low groups
  surv_data <- surv_data %>%
    mutate(
      dgat1_group = ifelse(dgat1 >= cutoff_result$cutoff, "High", "Low"),
      dgat1_group = factor(dgat1_group, levels = c("Low", "High"))
    )
  
  cat("  Low DGAT1:", sum(surv_data$dgat1_group == "Low"), "samples\n")
  cat("  High DGAT1:", sum(surv_data$dgat1_group == "High"), "samples\n")
  
  # Kaplan-Meier analysis
  surv_obj <- Surv(surv_data$os_time, surv_data$os_status)
  fit <- survfit(surv_obj ~ dgat1_group, data = surv_data)
  
  # Log-rank test
  logrank <- survdiff(surv_obj ~ dgat1_group, data = surv_data)
  logrank_p <- 1 - pchisq(logrank$chisq, df = 1)
  
  # Median survival
  median_surv <- summary(fit)$table[, "median"]
  
  cat("\n  Median survival:\n")
  cat("    Low DGAT1:", round(median_surv[1], 1), "months\n")
  cat("    High DGAT1:", round(median_surv[2], 1), "months\n")
  cat("  Log-rank p-value:", format.pval(logrank_p, digits = 3), "\n")
  
  # Cox proportional hazards - univariate
  cox_uni <- coxph(surv_obj ~ dgat1_group, data = surv_data)
  cox_uni_summary <- summary(cox_uni)
  
  hr_uni <- cox_uni_summary$conf.int[1, 1]
  hr_lower_uni <- cox_uni_summary$conf.int[1, 3]
  hr_upper_uni <- cox_uni_summary$conf.int[1, 4]
  cox_p_uni <- cox_uni_summary$coefficients[1, 5]
  
  cat("\n  Univariate Cox:\n")
  cat("    HR:", round(hr_uni, 2), "[", round(hr_lower_uni, 2), "-", 
      round(hr_upper_uni, 2), "]\n")
  cat("    p-value:", format.pval(cox_p_uni, digits = 3), "\n")
  
  # Multivariate Cox (if covariates available)
  cox_multi <- NULL
  if ("age" %in% names(surv_data) && "sex" %in% names(surv_data)) {
    surv_data_complete <- surv_data %>%
      filter(!is.na(age) & !is.na(sex))
    
    if (nrow(surv_data_complete) > 50) {
      cox_multi <- coxph(
        Surv(os_time, os_status) ~ dgat1_group + age + sex,
        data = surv_data_complete
      )
      cox_multi_summary <- summary(cox_multi)
      
      hr_multi <- cox_multi_summary$conf.int[1, 1]
      hr_lower_multi <- cox_multi_summary$conf.int[1, 3]
      hr_upper_multi <- cox_multi_summary$conf.int[1, 4]
      cox_p_multi <- cox_multi_summary$coefficients[1, 5]
      
      cat("\n  Multivariate Cox (adjusted for age, sex):\n")
      cat("    HR:", round(hr_multi, 2), "[", round(hr_lower_multi, 2), "-", 
          round(hr_upper_multi, 2), "]\n")
      cat("    p-value:", format.pval(cox_p_multi, digits = 3), "\n")
    }
  }
  
  return(list(
    surv_data = surv_data,
    fit = fit,
    logrank_p = logrank_p,
    cox_uni = cox_uni,
    cox_multi = cox_multi,
    median_surv = median_surv
  ))
}

plot_km_curve <- function(results, cutoff_result, outfile) {
  cat("\nGenerating publication-ready Kaplan-Meier plot...\n")
  
  # Calculate additional statistics
  n_low <- sum(results$surv_data$dgat1_group == "Low")
  n_high <- sum(results$surv_data$dgat1_group == "High")
  events_low <- sum(results$surv_data$dgat1_group == "Low" & results$surv_data$os_status == 1)
  events_high <- sum(results$surv_data$dgat1_group == "High" & results$surv_data$os_status == 1)
  
  # Extract survival data from fit object
  surv_summary <- summary(results$fit)
  
  # Create data frame for plotting
  surv_plot_data <- data.frame(
    time = surv_summary$time,
    surv = surv_summary$surv,
    lower = surv_summary$lower,
    upper = surv_summary$upper,
    group = rep(c("Low", "High"), times = c(sum(surv_summary$strata == "dgat1_group=Low"), 
                                           sum(surv_summary$strata == "dgat1_group=High")))
  )
  
  # Extract HR info
  cox_summary <- summary(results$cox_uni)
  hr <- cox_summary$conf.int[1, 1]
  hr_lower <- cox_summary$conf.int[1, 3]
  hr_upper <- cox_summary$conf.int[1, 4]
  hr_p <- cox_summary$coefficients[1, 5]
  
  # Create main survival plot
  p <- ggplot(surv_plot_data, aes(x = time, y = surv, color = group)) +
    geom_step(linewidth = 2) +
    geom_ribbon(aes(ymin = lower, ymax = upper, fill = group), alpha = 0.2, color = NA) +
    scale_color_manual(
      values = c("Low" = "#377EB8", "High" = "#E41A1C"),
      labels = c(
        paste0("Low (n=", n_low, ", events=", events_low, ")"),
        paste0("High (n=", n_high, ", events=", events_high, ")")
      )
    ) +
    scale_fill_manual(
      values = c("Low" = "#377EB8", "High" = "#E41A1C"),
      guide = "none"
    ) +
    labs(
      title = "TCGA-GBM: DGAT1 Expression and Overall Survival",
      x = "Time (months)",
      y = "Overall Survival Probability",
      color = "DGAT1 Expression"
    ) +
    theme_minimal(base_size = 14) +
    theme(
      plot.title = element_text(hjust = 0.5, size = 18, face = "bold"),
      axis.title.x = element_text(size = 16, face = "bold", margin = margin(t = 15)),
      axis.title.y = element_text(size = 16, face = "bold", margin = margin(r = 15)),
      axis.text = element_text(size = 14, color = "black"),
      legend.position = c(0.85, 0.35),
      legend.background = element_rect(fill = "white", color = "black", linewidth = 0.5),
      legend.key = element_rect(fill = "white"),
      panel.grid.major = element_line(color = "grey90", linewidth = 0.5),
      panel.grid.minor = element_blank(),
      panel.background = element_rect(fill = "white"),
      plot.margin = margin(20, 20, 20, 20)
    ) +
    ylim(0, 1) +
    scale_x_continuous(breaks = seq(0, max(surv_plot_data$time), 12))
  
  # Create annotation text box
  annotation_text <- paste0(
    "Hazard Ratio: ", sprintf("%.2f", hr),
    " (95% CI: ", sprintf("%.2f", hr_lower), "-", sprintf("%.2f", hr_upper), ")\n",
    "Log-rank p = ", format.pval(results$logrank_p, digits = 3), "\n",
    "Cox p = ", format.pval(hr_p, digits = 3), "\n",
    "\n",
    "Median survival:\n",
    "  Low: ", sprintf("%.1f", results$median_surv[1]), " months\n",
    "  High: ", sprintf("%.1f", results$median_surv[2]), " months"
  )
  
  # Add enhanced annotations
  p <- p + 
    # Text box with statistics
    annotate("rect", 
             xmin = max(surv_plot_data$time) * 0.55, 
             xmax = max(surv_plot_data$time) * 0.98,
             ymin = 0.65, ymax = 0.98,
             fill = "white", color = "black", linewidth = 0.8, alpha = 0.9) +
    annotate("text", 
             x = max(surv_plot_data$time) * 0.765, 
             y = 0.815,
             label = annotation_text,
             size = 4.5, hjust = 0.5, vjust = 0.5,
             family = "sans", fontface = "plain") +
    # Cutoff value annotation
    annotate("text",
             x = max(surv_plot_data$time) * 0.02,
             y = 0.05,
             label = paste0("Optimal cutoff: ", sprintf("%.2f", cutoff_result$cutoff)),
             size = 4, hjust = 0, family = "sans", fontface = "italic") +
    # Sample size annotation
    annotate("text",
             x = max(surv_plot_data$time) * 0.02,
             y = 0.12,
             label = paste0("Total n = ", nrow(results$surv_data)),
             size = 4, hjust = 0, family = "sans", fontface = "italic")
  
  # Create risk table
  risk_times <- seq(0, max(surv_plot_data$time), 12)
  risk_data <- data.frame(
    time = risk_times,
    low_risk = sapply(risk_times, function(t) sum(results$surv_data$dgat1_group == "Low" & results$surv_data$os_time >= t)),
    high_risk = sapply(risk_times, function(t) sum(results$surv_data$dgat1_group == "High" & results$surv_data$os_time >= t))
  )
  
  risk_table <- ggplot(risk_data, aes(x = time)) +
    geom_text(aes(y = 0.5, label = low_risk), color = "#377EB8", size = 4, fontface = "bold") +
    geom_text(aes(y = 0, label = high_risk), color = "#E41A1C", size = 4, fontface = "bold") +
    scale_x_continuous(breaks = seq(0, max(surv_plot_data$time), 12)) +
    scale_y_continuous(breaks = c(0, 0.5), labels = c("High", "Low")) +
    labs(title = "Number at risk", x = "Time (months)", y = "") +
    theme_minimal(base_size = 12) +
    theme(
      plot.title = element_text(size = 14, face = "bold"),
      axis.text = element_text(size = 12, color = "black", face = "bold"),
      axis.title = element_text(size = 12, face = "bold"),
      panel.grid = element_blank(),
      plot.margin = margin(5, 20, 5, 20)
    )
  
  # Combine plots
  if (COWPLOT_AVAILABLE) {
    combined_plot <- cowplot::plot_grid(
      p, risk_table,
      ncol = 1, rel_heights = c(3, 1),
      align = "v"
    )
    # Save high-resolution figure
    ggsave(outfile, combined_plot, width = 12, height = 10, dpi = 300)
    cat("  Saved:", outfile, "\n")
    return(combined_plot)
  } else {
    # Save main plot only
    ggsave(outfile, p, width = 12, height = 8, dpi = 300)
    cat("  Saved:", outfile, "\n")
    return(p)
  }
}

create_forest_plot <- function(results, outfile) {
  cat("\nGenerating forest plot...\n")
  
  # Extract HR data
  cox_uni_summary <- summary(results$cox_uni)
  
  hr_data <- data.frame(
    model = "Univariate",
    hr = cox_uni_summary$conf.int[1, 1],
    lower = cox_uni_summary$conf.int[1, 3],
    upper = cox_uni_summary$conf.int[1, 4],
    p = cox_uni_summary$coefficients[1, 5]
  )
  
  if (!is.null(results$cox_multi)) {
    cox_multi_summary <- summary(results$cox_multi)
    hr_data <- rbind(
      hr_data,
      data.frame(
        model = "Multivariate (age, sex)",
        hr = cox_multi_summary$conf.int[1, 1],
        lower = cox_multi_summary$conf.int[1, 3],
        upper = cox_multi_summary$conf.int[1, 4],
        p = cox_multi_summary$coefficients[1, 5]
      )
    )
  }
  
  # Create forest plot
  p <- ggplot(hr_data, aes(y = model, x = hr)) +
    geom_point(size = 4, color = "#E41A1C") +
    geom_errorbarh(aes(xmin = lower, xmax = upper), height = 0.2, linewidth = 1) +
    geom_vline(xintercept = 1, linetype = "dashed", color = "black") +
    geom_text(aes(label = paste0("HR=", round(hr, 2), 
                                  " [", round(lower, 2), "-", round(upper, 2), "]")),
              hjust = -0.1, size = 4) +
    geom_text(aes(label = paste0("p=", format.pval(p, digits = 3))),
              hjust = -0.1, vjust = -1.5, size = 4) +
    scale_x_continuous(trans = "log10", breaks = c(0.5, 1, 2, 4)) +
    labs(
      title = "DGAT1 High vs Low: Hazard Ratios for Overall Survival",
      subtitle = "TCGA GBM Cohort",
      x = "Hazard Ratio (log scale)",
      y = ""
    ) +
    theme_minimal(base_size = 14) +
    theme(
      plot.title = element_text(face = "bold", size = 16),
      axis.text.y = element_text(size = 12)
    )
  
  ggsave(outfile, p, width = 10, height = 4, dpi = 300)
  cat("  Saved:", outfile, "\n")
  
  return(p)
}

generate_summary_report <- function(results, cutoff_result, outfile) {
  cat("\nGenerating summary report...\n")
  
  report <- c(
    "# TCGA GBM: DGAT1 Survival Analysis Summary",
    "",
    "## Sample Characteristics",
    paste0("- Total samples: ", nrow(results$surv_data)),
    paste0("- Events (deaths): ", sum(results$surv_data$os_status == 1)),
    paste0("- Censored: ", sum(results$surv_data$os_status == 0)),
    paste0("- Median follow-up: ", round(median(results$surv_data$os_time), 1), " months"),
    "",
    "## DGAT1 Expression",
    paste0("- Optimal cutoff: ", round(cutoff_result$cutoff, 3)),
    paste0("- Low DGAT1: ", sum(results$surv_data$dgat1_group == "Low"), " samples"),
    paste0("- High DGAT1: ", sum(results$surv_data$dgat1_group == "High"), " samples"),
    "",
    "## Survival Outcomes",
    paste0("- Median survival (Low DGAT1): ", round(results$median_surv[1], 1), " months"),
    paste0("- Median survival (High DGAT1): ", round(results$median_surv[2], 1), " months"),
    paste0("- Log-rank p-value: ", format.pval(results$logrank_p, digits = 3)),
    "",
    "## Cox Regression",
    "### Univariate Model",
    paste0("- HR (High vs Low): ", round(summary(results$cox_uni)$conf.int[1, 1], 2),
           " [", round(summary(results$cox_uni)$conf.int[1, 3], 2), "-",
           round(summary(results$cox_uni)$conf.int[1, 4], 2), "]"),
    paste0("- p-value: ", format.pval(summary(results$cox_uni)$coefficients[1, 5], digits = 3))
  )
  
  if (!is.null(results$cox_multi)) {
    report <- c(
      report,
      "",
      "### Multivariate Model (adjusted for age, sex)",
      paste0("- HR (High vs Low): ", round(summary(results$cox_multi)$conf.int[1, 1], 2),
             " [", round(summary(results$cox_multi)$conf.int[1, 3], 2), "-",
             round(summary(results$cox_multi)$conf.int[1, 4], 2), "]"),
      paste0("- p-value: ", format.pval(summary(results$cox_multi)$coefficients[1, 5], digits = 3))
    )
  }
  
  writeLines(report, outfile)
  cat("  Saved:", outfile, "\n")
}

# =============================================================================
# MAIN
# =============================================================================

main <- function() {
  cat("\n", paste(rep("=", 70), collapse = ""), "\n")
  cat("TCGA GBM: DGAT1 Survival Analysis\n")
  cat(paste(rep("=", 70), collapse = ""), "\n\n")
  
  # Load data
  expr_mat <- load_expression_data(EXPRESSION_FILE)
  clinical <- load_clinical_data(CLINICAL_FILE)
  
  # Prepare survival data
  surv_data <- prepare_survival_data(expr_mat, clinical)
  
  # Find optimal cutoff
  cutoff_result <- find_optimal_cutoff(surv_data)
  
  # Perform survival analysis
  results <- perform_survival_analysis(surv_data, cutoff_result)
  
  # Save results
  write_csv(results$surv_data, file.path(OUT_DIR, "survival_data_with_groups.csv"))
  
  # Generate plots
  plot_km_curve(results, cutoff_result, 
                file.path(OUT_DIR, "kaplan_meier_curve.png"))
  
  create_forest_plot(results, 
                     file.path(OUT_DIR, "forest_plot_hazard_ratios.png"))
  
  # Generate summary report
  generate_summary_report(results, cutoff_result, 
                          file.path(OUT_DIR, "survival_analysis_summary.txt"))
  
  cat("\n", paste(rep("=", 70), collapse = ""), "\n")
  cat("Analysis complete! Results in:", OUT_DIR, "\n")
  cat(paste(rep("=", 70), collapse = ""), "\n\n")
  
  invisible(results)
}

if (!interactive()) {
  results <- main()
}
