# Complete Survival Analysis Workflow for DGAT Immunology

This protocol provides a comprehensive workflow for conducting survival analysis using both TCGA (mRNA) and CPTAC (protein) data for DGAT immunology research.

## 📋 Overview

This workflow covers:
1. **mRNA-based survival analysis** using TCGA data
2. **Protein-based survival analysis** using CPTAC data
3. **Multi-omics integration** and comparison
4. **Statistical validation** and visualization

## 🧬 Part I: mRNA-Based Survival Analysis (TCGA)

### Step 1: Data Preparation
```r
# Load required libraries
library(survival)
library(survminer)
library(dplyr)
library(ggplot2)
library(readr)

# Load TCGA data
tcga_expression <- read_csv("Raw_Data/TCGA/TCGA_mRNA_expression.csv")
tcga_clinical <- read_csv("Raw_Data/TCGA/TCGA_clinical.csv")

# Extract DGAT1 and DGAT2 expression
dgat_expression <- tcga_expression %>%
  filter(gene_name %in% c("DGAT1", "DGAT2")) %>%
  pivot_longer(cols = -gene_name, names_to = "sample_id", values_to = "expression")

# Merge with clinical data
survival_data <- tcga_clinical %>%
  left_join(dgat_expression, by = "sample_id")
```

### Step 2: Survival Analysis Setup
```r
# Create survival object
survival_obj <- Surv(time = survival_data$OS_time, 
                     event = survival_data$OS_event)

# Define DGAT expression groups (median cut)
survival_data$DGAT1_group <- ifelse(survival_data$DGAT1_expression > median(survival_data$DGAT1_expression, na.rm = TRUE), 
                                   "High", "Low")
survival_data$DGAT2_group <- ifelse(survival_data$DGAT2_expression > median(survival_data$DGAT2_expression, na.rm = TRUE), 
                                   "High", "Low")
```

### Step 3: Cox Proportional Hazards Analysis
```r
# DGAT1 Cox regression
dgat1_cox <- coxph(survival_obj ~ DGAT1_expression + age + gender + stage, 
                   data = survival_data)

# DGAT2 Cox regression
dgat2_cox <- coxph(survival_obj ~ DGAT2_expression + age + gender + stage, 
                   data = survival_data)

# Summary of results
summary(dgat1_cox)
summary(dgat2_cox)
```

### Step 4: Kaplan-Meier Analysis
```r
# DGAT1 Kaplan-Meier
dgat1_fit <- survfit(survival_obj ~ DGAT1_group, data = survival_data)

# DGAT2 Kaplan-Meier
dgat2_fit <- survfit(survival_obj ~ DGAT2_group, data = survival_data)

# Plot Kaplan-Meier curves
p1 <- ggsurvplot(dgat1_fit, 
                 data = survival_data,
                 pval = TRUE,
                 conf.int = TRUE,
                 title = "DGAT1 Expression and Overall Survival",
                 xlab = "Time (months)",
                 ylab = "Survival Probability")

p2 <- ggsurvplot(dgat2_fit, 
                 data = survival_data,
                 pval = TRUE,
                 conf.int = TRUE,
                 title = "DGAT2 Expression and Overall Survival",
                 xlab = "Time (months)",
                 ylab = "Survival Probability")
```

## 🧬 Part II: Protein-Based Survival Analysis (CPTAC)

### Step 1: CPTAC Data Preparation
```r
# Load CPTAC proteomics data
cptac_protein <- read_tsv("Raw_Data/Proteome/CPTAC/CPTAC3_Glioblastoma_Multiforme_Proteome.tmt11.tsv")
cptac_clinical <- read_excel("Raw_Data/Proteome/CPTAC/S048_CPTAC_GBM_Discovery_Cohort_Clinical_Data_Dec2019_r1.xlsx")

# Extract DGAT1 and DGAT2 protein levels
dgat_protein <- cptac_protein %>%
  filter(Gene %in% c("DGAT1", "DGAT2")) %>%
  select(Gene, starts_with("TCGA"))

# Transpose for analysis
dgat_protein_long <- dgat_protein %>%
  pivot_longer(cols = -Gene, names_to = "sample_id", values_to = "protein_level") %>%
  pivot_wider(names_from = Gene, values_from = protein_level)

# Merge with clinical data
protein_survival_data <- cptac_clinical %>%
  left_join(dgat_protein_long, by = "sample_id")
```

### Step 2: Protein Survival Analysis
```r
# Create survival object for protein data
protein_survival_obj <- Surv(time = protein_survival_data$OS_time, 
                             event = protein_survival_data$OS_event)

# DGAT1 protein Cox regression
dgat1_protein_cox <- coxph(protein_survival_obj ~ DGAT1 + age + gender + stage, 
                           data = protein_survival_data)

# DGAT2 protein Cox regression
dgat2_protein_cox <- coxph(protein_survival_obj ~ DGAT2 + age + gender + stage, 
                           data = protein_survival_data)

# Summary
summary(dgat1_protein_cox)
summary(dgat2_protein_cox)
```

### Step 3: Protein Kaplan-Meier Analysis
```r
# Define protein expression groups
protein_survival_data$DGAT1_protein_group <- ifelse(protein_survival_data$DGAT1 > median(protein_survival_data$DGAT1, na.rm = TRUE), 
                                                    "High", "Low")
protein_survival_data$DGAT2_protein_group <- ifelse(protein_survival_data$DGAT2 > median(protein_survival_data$DGAT2, na.rm = TRUE), 
                                                    "High", "Low")

# Protein Kaplan-Meier curves
dgat1_protein_fit <- survfit(protein_survival_obj ~ DGAT1_protein_group, data = protein_survival_data)
dgat2_protein_fit <- survfit(protein_survival_obj ~ DGAT2_protein_group, data = protein_survival_data)

# Plot protein survival curves
p3 <- ggsurvplot(dgat1_protein_fit, 
                 data = protein_survival_data,
                 pval = TRUE,
                 conf.int = TRUE,
                 title = "DGAT1 Protein and Overall Survival",
                 xlab = "Time (months)",
                 ylab = "Survival Probability")

p4 <- ggsurvplot(dgat2_protein_fit, 
                 data = protein_survival_data,
                 pval = TRUE,
                 conf.int = TRUE,
                 title = "DGAT2 Protein and Overall Survival",
                 xlab = "Time (months)",
                 ylab = "Survival Probability")
```

## 🔬 Part III: Multi-omics Integration

### Step 1: RNA-Protein Correlation
```r
# Calculate correlation between mRNA and protein levels
# (This requires matched samples between TCGA and CPTAC)

# For demonstration, create correlation analysis
correlation_data <- data.frame(
  mRNA_DGAT1 = survival_data$DGAT1_expression,
  mRNA_DGAT2 = survival_data$DGAT2_expression,
  # Note: Protein data would need to be matched by sample ID
  protein_DGAT1 = NA,  # Would be matched protein levels
  protein_DGAT2 = NA   # Would be matched protein levels
)

# Calculate correlations (when data is available)
# cor_mrna_protein_dgat1 <- cor(correlation_data$mRNA_DGAT1, correlation_data$protein_DGAT1, use = "complete.obs")
# cor_mrna_protein_dgat2 <- cor(correlation_data$mRNA_DGAT2, correlation_data$protein_DGAT2, use = "complete.obs")
```

### Step 2: Combined Analysis
```r
# Create forest plot comparing mRNA vs protein results
library(forestplot)

# Prepare data for forest plot
forest_data <- data.frame(
  Analysis = c("DGAT1 mRNA", "DGAT2 mRNA", "DGAT1 Protein", "DGAT2 Protein"),
  HR = c(exp(coef(dgat1_cox)["DGAT1_expression"]),
         exp(coef(dgat2_cox)["DGAT2_expression"]),
         exp(coef(dgat1_protein_cox)["DGAT1"]),
         exp(coef(dgat2_protein_cox)["DGAT2"])),
  CI_lower = c(exp(confint(dgat1_cox)["DGAT1_expression", 1]),
               exp(confint(dgat2_cox)["DGAT2_expression", 1]),
               exp(confint(dgat1_protein_cox)["DGAT1", 1]),
               exp(confint(dgat2_protein_cox)["DGAT2", 1])),
  CI_upper = c(exp(confint(dgat1_cox)["DGAT1_expression", 2]),
               exp(confint(dgat2_cox)["DGAT2_expression", 2]),
               exp(confint(dgat1_protein_cox)["DGAT1", 2]),
               exp(confint(dgat2_protein_cox)["DGAT2", 2])),
  P_value = c(summary(dgat1_cox)$coefficients["DGAT1_expression", "Pr(>|z|)"],
              summary(dgat2_cox)$coefficients["DGAT2_expression", "Pr(>|z|)"],
              summary(dgat1_protein_cox)$coefficients["DGAT1", "Pr(>|z|)"],
              summary(dgat2_protein_cox)$coefficients["DGAT2", "Pr(>|z|)"])
)

# Create forest plot
forestplot(forest_data$Analysis,
           forest_data$HR,
           forest_data$CI_lower,
           forest_data$CI_upper,
           title = "DGAT Expression and Survival: mRNA vs Protein",
           xlab = "Hazard Ratio (95% CI)",
           zero = 1)
```

## 📊 Part IV: Statistical Validation

### Step 1: Model Assumptions
```r
# Check proportional hazards assumption
cox.zph(dgat1_cox)
cox.zph(dgat2_cox)
cox.zph(dgat1_protein_cox)
cox.zph(dgat2_protein_cox)

# Check model fit
plot(survfit(dgat1_cox), main = "DGAT1 Cox Model Fit")
plot(survfit(dgat2_cox), main = "DGAT2 Cox Model Fit")
```

### Step 2: Sensitivity Analysis
```r
# Different cutpoints for expression groups
# Quartile-based grouping
survival_data$DGAT1_quartile <- ntile(survival_data$DGAT1_expression, 4)
survival_data$DGAT2_quartile <- ntile(survival_data$DGAT2_expression, 4)

# Survival by quartiles
dgat1_quartile_fit <- survfit(survival_obj ~ DGAT1_quartile, data = survival_data)
dgat2_quartile_fit <- survfit(survival_obj ~ DGAT2_quartile, data = survival_data)

# Plot quartile survival curves
p5 <- ggsurvplot(dgat1_quartile_fit, 
                 data = survival_data,
                 pval = TRUE,
                 title = "DGAT1 Expression Quartiles and Survival")
```

## 📈 Part V: Results Summary and Visualization

### Step 1: Create Summary Tables
```r
# Create summary table of all analyses
summary_table <- data.frame(
  Gene = c("DGAT1", "DGAT2", "DGAT1", "DGAT2"),
  Data_Type = c("mRNA", "mRNA", "Protein", "Protein"),
  HR = c(exp(coef(dgat1_cox)["DGAT1_expression"]),
         exp(coef(dgat2_cox)["DGAT2_expression"]),
         exp(coef(dgat1_protein_cox)["DGAT1"]),
         exp(coef(dgat2_protein_cox)["DGAT2"])),
  CI_95 = c(paste0("(", round(exp(confint(dgat1_cox)["DGAT1_expression", 1]), 2), 
                   "-", round(exp(confint(dgat1_cox)["DGAT1_expression", 2]), 2), ")"),
            paste0("(", round(exp(confint(dgat2_cox)["DGAT2_expression", 1]), 2), 
                   "-", round(exp(confint(dgat2_cox)["DGAT2_expression", 2]), 2), ")"),
            paste0("(", round(exp(confint(dgat1_protein_cox)["DGAT1", 1]), 2), 
                   "-", round(exp(confint(dgat1_protein_cox)["DGAT1", 2]), 2), ")"),
            paste0("(", round(exp(confint(dgat2_protein_cox)["DGAT2", 1]), 2), 
                   "-", round(exp(confint(dgat2_protein_cox)["DGAT2", 2]), 2), ")")),
  P_value = c(summary(dgat1_cox)$coefficients["DGAT1_expression", "Pr(>|z|)"],
              summary(dgat2_cox)$coefficients["DGAT2_expression", "Pr(>|z|)"],
              summary(dgat1_protein_cox)$coefficients["DGAT1", "Pr(>|z|)"],
              summary(dgat2_protein_cox)$coefficients["DGAT2", "Pr(>|z|)"])
)

# Save summary table
write_csv(summary_table, "Results/Survival/survival_analysis_summary.csv")
```

### Step 2: Combine All Plots
```r
# Arrange all survival plots
library(gridExtra)

# Combine mRNA plots
mrna_plots <- grid.arrange(p1$plot, p2$plot, ncol = 2)

# Combine protein plots
protein_plots <- grid.arrange(p3$plot, p4$plot, ncol = 2)

# Save plots
ggsave("Results/Survival/DGAT_mRNA_survival_curves.png", mrna_plots, width = 12, height = 6)
ggsave("Results/Survival/DGAT_protein_survival_curves.png", protein_plots, width = 12, height = 6)
```

## 🔍 Part VI: Advanced Analysis

### Step 1: Time-Dependent Analysis
```r
# Time-dependent Cox regression
library(timeROC)

# Calculate time-dependent AUC
time_roc_dgat1 <- timeROC(T = survival_data$OS_time,
                          delta = survival_data$OS_event,
                          marker = survival_data$DGAT1_expression,
                          cause = 1,
                          times = c(12, 24, 36, 48, 60))

# Plot time-dependent ROC
plot(time_roc_dgat1, time = 24, col = "red", title = "DGAT1 Time-dependent ROC (24 months)")
```

### Step 2: Subgroup Analysis
```r
# Subgroup analysis by clinical characteristics
# Age subgroups
survival_data$age_group <- ifelse(survival_data$age > 65, ">65", "≤65")

# Survival analysis by age group
dgat1_age_fit <- survfit(survival_obj ~ DGAT1_group + age_group, data = survival_data)

# Plot stratified survival curves
ggsurvplot(dgat1_age_fit, 
           data = survival_data,
           pval = TRUE,
           title = "DGAT1 Survival by Age Group",
           facet.by = "age_group")
```

## 📝 Part VII: Reporting and Documentation

### Step 1: Generate Report
```r
# Create comprehensive survival analysis report
report_text <- paste0(
  "DGAT Survival Analysis Report\n",
  "============================\n\n",
  "Analysis Date: ", Sys.Date(), "\n",
  "Total Samples: ", nrow(survival_data), "\n\n",
  
  "mRNA Analysis Results:\n",
  "- DGAT1 HR: ", round(exp(coef(dgat1_cox)["DGAT1_expression"]), 3), 
  " (95% CI: ", round(exp(confint(dgat1_cox)["DGAT1_expression", 1]), 3), 
  "-", round(exp(confint(dgat1_cox)["DGAT1_expression", 2]), 3), 
  "), p = ", round(summary(dgat1_cox)$coefficients["DGAT1_expression", "Pr(>|z|)"], 3), "\n",
  
  "- DGAT2 HR: ", round(exp(coef(dgat2_cox)["DGAT2_expression"]), 3), 
  " (95% CI: ", round(exp(confint(dgat2_cox)["DGAT2_expression", 1]), 3), 
  "-", round(exp(confint(dgat2_cox)["DGAT2_expression", 2]), 3), 
  "), p = ", round(summary(dgat2_cox)$coefficients["DGAT2_expression", "Pr(>|z|)"], 3), "\n\n",
  
  "Protein Analysis Results:\n",
  "- DGAT1 HR: ", round(exp(coef(dgat1_protein_cox)["DGAT1"]), 3), 
  " (95% CI: ", round(exp(confint(dgat1_protein_cox)["DGAT1", 1]), 3), 
  "-", round(exp(confint(dgat1_protein_cox)["DGAT1", 2]), 3), 
  "), p = ", round(summary(dgat1_protein_cox)$coefficients["DGAT1", "Pr(>|z|)"], 3), "\n",
  
  "- DGAT2 HR: ", round(exp(coef(dgat2_protein_cox)["DGAT2"]), 3), 
  " (95% CI: ", round(exp(confint(dgat2_protein_cox)["DGAT2", 1]), 3), 
  "-", round(exp(confint(dgat2_protein_cox)["DGAT2", 2]), 3), 
  "), p = ", round(summary(dgat2_protein_cox)$coefficients["DGAT2", "Pr(>|z|)"], 3), "\n"
)

# Save report
writeLines(report_text, "Results/Survival/survival_analysis_report.txt")
```

## ✅ Quality Control Checklist

- [ ] **Data Quality**: Verify all survival data is complete and accurate
- [ ] **Model Assumptions**: Check proportional hazards assumption
- [ ] **Statistical Power**: Ensure adequate sample size for analysis
- [ ] **Multiple Testing**: Apply appropriate correction for multiple comparisons
- [ ] **Reproducibility**: Document all analysis parameters and versions
- [ ] **Validation**: Cross-validate findings with independent datasets

## 📚 References

1. Therneau, T. M., & Grambsch, P. M. (2000). Modeling survival data: extending the Cox model.
2. Kassambara, A., Kosinski, M., & Biecek, P. (2021). survminer: Drawing Survival Curves using 'ggplot2'.
3. Blighe, K., Rana, S., & Lewis, M. (2021). EnhancedVolcano: Publication-ready volcano plots.

---

**This protocol provides a comprehensive framework for conducting robust survival analysis in DGAT immunology research using both mRNA and protein data.**
