# =============================================================================
# project_config.R — Project Configuration Settings
# =============================================================================
#
# Central configuration file for the DGAT Immunology analysis pipeline.
# Modify these settings to customize the analysis for different datasets
# or analysis parameters.
#
# =============================================================================

# Project metadata
PROJECT_NAME <- "DGAT_Immunology"
PROJECT_VERSION <- "1.0.0"
PROJECT_DATE <- "2024-09-29"

# Directory paths (relative to project root)
DIRS <- list(
  raw_data = "Raw_Data",
  processed_data = "Processed_Data", 
  results = "Results",
  scripts = "Scripts",
  docs = "Docs",
  figures = "Figures"
)

# Cancer types to analyze
CANCER_TYPES <- c("BRCA", "GBM", "OV", "PAAD")

# DGAT genes of interest
DGAT_GENES <- c(
  "DGAT1", "DGAT2",           # Main DGAT enzymes
  "ACAT1", "ACAT2",           # Acetyl-CoA acetyltransferases
  "SCD", "FASN", "ACC1",      # Lipid metabolism genes
  "LIPE", "ATGL",             # Lipolysis genes
  "CPT1A", "CPT1B", "CPT1C", # Carnitine palmitoyltransferases
  "SREBF1", "SREBF2",         # Sterol regulatory element binding factors
  "PPARA", "PPARG"            # Peroxisome proliferator-activated receptors
)

# Immune signature sources
IMMUNE_SIGNATURES <- list(
  iobr = list(
    t_cells = "sig_tcell",
    b_cells = "sig_bcell", 
    nk_cells = "sig_nk",
    macrophages = "sig_macrophage",
    dendritic_cells = "sig_dc",
    neutrophils = "sig_neutrophil",
    mast_cells = "sig_mast",
    eosinophils = "sig_eosinophil",
    monocytes = "sig_monocyte"
  ),
  custom = list(
    # Add custom signatures here
  )
)

# Analysis parameters
ANALYSIS_PARAMS <- list(
  # Expression filtering
  min_gene_expression = 1,      # Minimum CPM/TPM for gene filtering
  min_sample_genes = 5000,      # Minimum genes per sample
  
  # Correlation analysis
  correlation_method = "spearman",
  min_sample_size = 50,
  
  # Survival analysis
  survival_time_col = "OS.time",
  survival_status_col = "OS", 
  survival_cutoff = "median",
  
  # GSVA parameters
  gsva_method = "ssgsea",
  min_geneset_size = 10,
  max_geneset_size = 500,
  
  # Plotting parameters
  plot_width = 10,
  plot_height = 8,
  plot_dpi = 300,
  plot_theme = "minimal"
)

# Database download settings
DOWNLOAD_PARAMS <- list(
  # TCGA
  tcga = list(
    data_category = "Transcriptome Profiling",
    data_type = "Gene Expression Quantification",
    workflow_type = "HTSeq - Counts",
    sample_type = c("Primary Tumor", "Solid Tissue Normal")
  ),
  
  # CGGA
  cgga = list(
    base_url = "http://www.cgga.org.cn/",
    datasets = c("CGGA_325", "CGGA_693", "CGGA_301")
  ),
  
  # GTEx
  gtex = list(
    version = "v8",
    tissue_types = c("Brain", "Liver", "Pancreas", "Ovary", "Breast")
  )
)

# File naming conventions
FILE_PATTERNS <- list(
  expression_matrix = "{dataset}_{cancer_type}_expression_matrix_{normalization}.csv",
  clinical_data = "{dataset}_{cancer_type}_clinical_data.csv",
  correlation_results = "{cancer_type}_DGAT_immune_correlations.csv",
  survival_results = "{cancer_type}_DGAT_survival_results.csv",
  gsva_scores = "{cancer_type}_immune_gsva_scores.csv",
  heatmap_plot = "{cancer_type}_DGAT_immune_heatmap.png",
  survival_plot = "{cancer_type}_DGAT_survival_km.png"
)

# Quality control thresholds
QC_THRESHOLDS <- list(
  # Gene filtering
  min_gene_count = 1000,
  max_gene_count = 50000,
  
  # Sample filtering  
  min_lib_size = 1000000,
  max_lib_size = 100000000,
  
  # Mitochondrial genes
  max_mito_percent = 20,
  
  # Ribosomal genes
  max_ribo_percent = 50
)

# Statistical parameters
STATS_PARAMS <- list(
  # Multiple testing correction
  fdr_method = "BH",
  alpha_level = 0.05,
  
  # Bootstrap parameters
  n_bootstrap = 1000,
  bootstrap_seed = 42,
  
  # Cross-validation
  cv_folds = 5,
  cv_seed = 42
)

# Visualization settings
PLOT_SETTINGS <- list(
  # Color palettes
  colors = list(
    cancer_types = c("BRCA" = "#E31A1C", "GBM" = "#1F78B4", 
                     "OV" = "#33A02C", "PAAD" = "#FF7F00"),
    immune_cells = c("T_cells" = "#E31A1C", "B_cells" = "#1F78B4",
                     "NK_cells" = "#33A02C", "Macrophages" = "#FF7F00",
                     "Dendritic_cells" = "#6A3D9A", "Neutrophils" = "#B15928"),
    survival_groups = c("High" = "#E31A1C", "Low" = "#1F78B4")
  ),
  
  # Theme settings
  theme_base_size = 12,
  theme_font_family = "Arial",
  
  # Heatmap settings
  heatmap_cluster_rows = TRUE,
  heatmap_cluster_cols = TRUE,
  heatmap_scale = "row"
)

# Print configuration summary
cat("📋 DGAT Immunology Project Configuration Loaded\n")
cat("   Project:", PROJECT_NAME, "v", PROJECT_VERSION, "\n")
cat("   Cancer types:", paste(CANCER_TYPES, collapse = ", "), "\n")
cat("   DGAT genes:", length(DGAT_GENES), "genes\n")
cat("   Analysis date:", PROJECT_DATE, "\n")
cat("✅ Configuration ready\n\n")
