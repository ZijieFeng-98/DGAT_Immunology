#!/usr/bin/env Rscript
# =============================================================================
# STREAMLINED IMMUNE & METABOLIC MARKER SET
# Publication-ready with dual classification:
#   1) Cell Type (for biological context)
#   2) Immune Function (Suppressive vs. Anti-tumor)
# =============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
})

# =============================================================================
# CURATED MARKER SET - QUALITY OVER QUANTITY
# =============================================================================

CURATED_MARKERS <- tribble(
  ~Gene,        ~CellType,            ~ImmuneFunction,        ~Description,
  
  # =========================================================================
  # MACROPHAGES & MICROGLIA
  # =========================================================================
  # Pan-macrophage markers
  "CD68",       "Macrophage",         "Pan-Myeloid",          "Pan-macrophage/microglia marker",
  "CD163",      "Macrophage",         "Immunosuppressive",    "M2 macrophage, scavenger receptor",
  "MRC1",       "Macrophage",         "Immunosuppressive",    "CD206, mannose receptor (M2)",
  "APOE",       "Macrophage",         "Immunosuppressive",    "Apolipoprotein E, lipid-loaded TAMs",
  
  # M2/TAM markers (immunosuppressive)
  "ARG1",       "TAM",                "Immunosuppressive",    "Arginase-1, depletes arginine, suppresses T cells",
  "IL10",       "TAM",                "Immunosuppressive",    "Anti-inflammatory cytokine",
  "TGFB1",      "TAM",                "Immunosuppressive",    "TGF-beta, suppresses effector T cells",
  "SPP1",       "TAM",                "Immunosuppressive",    "Osteopontin, TAM marker, pro-tumor",
  "VEGFA",      "TAM",                "Immunosuppressive",    "Angiogenic factor, immune evasion",
  
  # M1 markers (pro-inflammatory, anti-tumor potential)
  "NOS2",       "M1-Macrophage",      "Pro-inflammatory",     "iNOS, produces NO, M1 marker",
  "IL12A",      "M1-Macrophage",      "Anti-tumor",           "IL-12p35, Th1 polarization",
  "CXCL10",     "M1-Macrophage",      "Anti-tumor",           "IP-10, T cell recruitment",
  
  # Microglia-specific
  "P2RY12",     "Microglia",          "Pan-Myeloid",          "Microglia-specific marker",
  "TMEM119",    "Microglia",          "Pan-Myeloid",          "Microglia-specific marker",
  
  # =========================================================================
  # MDSCs (MYELOID-DERIVED SUPPRESSOR CELLS)
  # =========================================================================
  "S100A8",     "MDSC",               "Immunosuppressive",    "Alarmin, M-MDSC marker",
  "S100A9",     "MDSC",               "Immunosuppressive",    "Alarmin, M-MDSC marker",
  "CD14",       "MDSC",               "Immunosuppressive",    "Monocytic MDSC marker",
  "ARG2",       "MDSC",               "Immunosuppressive",    "Arginase-2, T cell suppression",
  "IDO1",       "MDSC",               "Immunosuppressive",    "Tryptophan catabolism, T cell suppression",
  
  # =========================================================================
  # T CELLS - EFFECTOR (ANTI-TUMOR)
  # =========================================================================
  # CD8+ cytotoxic T cells
  "CD8A",       "CD8-T",              "Anti-tumor",           "CD8+ T cell marker",
  "CD8B",       "CD8-T",              "Anti-tumor",           "CD8+ T cell marker",
  
  # CD4+ helper T cells
  "CD4",        "CD4-T",              "Anti-tumor",           "CD4+ T helper marker",
  
  # Cytotoxic effectors
  "GZMA",       "Cytotoxic",          "Anti-tumor",           "Granzyme A, kills tumor cells",
  "GZMB",       "Cytotoxic",          "Anti-tumor",           "Granzyme B, kills tumor cells",
  "PRF1",       "Cytotoxic",          "Anti-tumor",           "Perforin, pore-forming cytotoxic",
  "IFNG",       "Th1",                "Anti-tumor",           "IFN-gamma, Th1 effector cytokine",
  
  # T cell activation
  "CD3E",       "Pan-T",              "Anti-tumor",           "T cell receptor component",
  
  # =========================================================================
  # T CELLS - EXHAUSTED (DYSFUNCTIONAL)
  # =========================================================================
  "PDCD1",      "Exhausted-T",        "Dysfunctional",        "PD-1, exhaustion marker",
  "LAG3",       "Exhausted-T",        "Dysfunctional",        "LAG3, exhaustion marker",
  "HAVCR2",     "Exhausted-T",        "Dysfunctional",        "TIM-3, exhaustion marker",
  "TIGIT",      "Exhausted-T",        "Dysfunctional",        "TIGIT, exhaustion/suppression",
  "CTLA4",      "Exhausted-T",        "Dysfunctional",        "CTLA-4, immune checkpoint",
  "TOX",        "Exhausted-T",        "Dysfunctional",        "TOX, exhaustion transcription factor",
  
  # =========================================================================
  # REGULATORY T CELLS (IMMUNOSUPPRESSIVE)
  # =========================================================================
  "FOXP3",      "Treg",               "Immunosuppressive",    "Treg master transcription factor",
  "IL2RA",      "Treg",               "Immunosuppressive",    "CD25, IL-2 receptor alpha",
  
  # =========================================================================
  # NK CELLS (ANTI-TUMOR)
  # =========================================================================
  "NCAM1",      "NK-cell",            "Anti-tumor",           "CD56, NK cell marker",
  "KLRD1",      "NK-cell",            "Anti-tumor",           "CD94, NK receptor",
  "KLRK1",      "NK-cell",            "Anti-tumor",           "NKG2D, activating NK receptor",
  "NCR1",       "NK-cell",            "Anti-tumor",           "NKp46, NK activating receptor",
  
  # =========================================================================
  # DENDRITIC CELLS (ANTIGEN PRESENTATION)
  # =========================================================================
  "CLEC9A",     "DC",                 "Anti-tumor",           "cDC1, cross-presenting DC",
  "BATF3",      "DC",                 "Anti-tumor",           "cDC1 transcription factor",
  "CD1C",       "DC",                 "Anti-tumor",           "cDC2 marker",
  "FCER1A",     "DC",                 "Anti-tumor",           "cDC2 marker",
  
  # =========================================================================
  # B CELLS
  # =========================================================================
  "CD19",       "B-cell",             "Adaptive-immune",      "B cell marker",
  "MS4A1",      "B-cell",             "Adaptive-immune",      "CD20, B cell marker",
  
  # =========================================================================
  # NEUTROPHILS
  # =========================================================================
  "FCGR3B",     "Neutrophil",         "Context-dependent",    "CD16b, neutrophil marker",
  "CSF3R",      "Neutrophil",         "Context-dependent",    "G-CSF receptor",
  
  # =========================================================================
  # IMMUNE CHECKPOINTS & LIGANDS
  # =========================================================================
  "CD274",      "Checkpoint",         "Immunosuppressive",    "PD-L1, checkpoint ligand",
  "PDCD1LG2",   "Checkpoint",         "Immunosuppressive",    "PD-L2, checkpoint ligand",
  "CD276",      "Checkpoint",         "Immunosuppressive",    "B7-H3, immune checkpoint",
  
  # =========================================================================
  # CYTOKINES & CHEMOKINES
  # =========================================================================
  # Pro-inflammatory (anti-tumor potential)
  "IL6",        "Cytokine",           "Pro-inflammatory",     "IL-6, pleiotropic cytokine",
  "IL12B",      "Cytokine",           "Anti-tumor",           "IL-12p40, Th1 polarization",
  "IL15",       "Cytokine",           "Anti-tumor",           "IL-15, NK/CD8 T cell support",
  "TNF",        "Cytokine",           "Pro-inflammatory",     "TNF-alpha, pro-inflammatory",
  
  # Anti-inflammatory (immunosuppressive)
  "IL4",        "Cytokine",           "Immunosuppressive",    "IL-4, Th2/M2 polarization",
  
  # Chemokines
  "CCL5",       "Chemokine",          "Anti-tumor",           "RANTES, T cell recruitment",
  "CCL2",       "Chemokine",          "Immunosuppressive",    "MCP-1, monocyte/MDSC recruitment",
  
  # =========================================================================
  # LIPID METABOLISM (STREAMLINED - 7 KEY GENES)
  # =========================================================================
  "DGAT1",      "Lipid-LD",           "Metabolic",            "Lipid droplet formation, TAM suppression",
  "SREBF1",     "Lipid-Synthesis",    "Metabolic",            "SREBP-1, master lipogenesis regulator",
  "SCAP",       "Lipid-Regulation",   "Metabolic",            "SREBP escort protein",
  "ACACA",      "Lipid-Synthesis",    "Metabolic",            "ACC, rate-limiting FA synthesis",
  "FASN",       "Lipid-Synthesis",    "Metabolic",            "Fatty acid synthase, de novo FA synthesis",
  "SCD",        "Lipid-Synthesis",    "Metabolic",            "SCD1, MUFA production",
  "HMGCR",      "Lipid-Synthesis",    "Metabolic",            "HMG-CoA reductase, cholesterol synthesis",
  
  # =========================================================================
  # HYPOXIA
  # =========================================================================
  "HIF1A",      "Hypoxia",            "Metabolic",            "HIF-1alpha, metabolic rewiring",
  
  # =========================================================================
  # ANTIGEN PRESENTATION
  # =========================================================================
  "HLA-A",      "MHC-I",              "Anti-tumor",           "MHC class I, antigen presentation",
  "HLA-DRA",    "MHC-II",             "Anti-tumor",           "MHC class II, antigen presentation",
  "B2M",        "MHC-I",              "Anti-tumor",           "Beta-2-microglobulin, MHC-I component",
  "TAP1",       "Antigen-Process",    "Anti-tumor",           "Peptide transporter, antigen processing"
)

# =============================================================================
# COLOR PALETTES FOR VISUALIZATION
# =============================================================================

# Cell Type Colors (biological classification)
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

# Immune Function Colors (functional classification)
FUNCTION_COLORS <- c(
  "Anti-tumor" = "#228B22",            # Dark green - good guys
  "Immunosuppressive" = "#B22222",     # Firebrick red - bad guys
  "Dysfunctional" = "#808080",         # Gray - broken
  "Pro-inflammatory" = "#FF8C00",      # Dark orange - context-dependent
  "Pan-Myeloid" = "#D3D3D3",          # Light gray - neutral markers
  "Adaptive-immune" = "#4169E1",       # Royal blue - B cells
  "Context-dependent" = "#DAA520",     # Goldenrod - depends on context
  "Metabolic" = "#8B4513"              # Saddle brown - metabolism
)

# =============================================================================
# SUMMARY STATISTICS
# =============================================================================

cat("\n=== STREAMLINED MARKER SET SUMMARY ===\n\n")
cat("Total markers:", nrow(CURATED_MARKERS), "\n\n")

cat("Markers by Cell Type:\n")
celltype_summary <- CURATED_MARKERS %>%
  count(CellType, sort = TRUE)
print(celltype_summary)

cat("\n\nMarkers by Immune Function:\n")
function_summary <- CURATED_MARKERS %>%
  count(ImmuneFunction, sort = TRUE)
print(function_summary)

cat("\n=== FUNCTIONAL BREAKDOWN ===\n\n")

# Anti-tumor markers
anti_tumor <- CURATED_MARKERS %>% 
  filter(ImmuneFunction == "Anti-tumor") %>% 
  nrow()
cat("Anti-tumor markers:", anti_tumor, "\n")

# Immunosuppressive markers
immunosupp <- CURATED_MARKERS %>% 
  filter(ImmuneFunction == "Immunosuppressive") %>% 
  nrow()
cat("Immunosuppressive markers:", immunosupp, "\n")

# Dysfunctional markers
dysfunc <- CURATED_MARKERS %>% 
  filter(ImmuneFunction == "Dysfunctional") %>% 
  nrow()
cat("Dysfunctional markers:", dysfunc, "\n")

# Metabolic markers
metabolic <- CURATED_MARKERS %>% 
  filter(ImmuneFunction == "Metabolic") %>% 
  nrow()
cat("Metabolic markers:", metabolic, "\n\n")

cat("=== KEY INSIGHTS ===\n\n")
cat("Ratio of Anti-tumor : Immunosuppressive =", 
    anti_tumor, ":", immunosupp, 
    "(", round(anti_tumor/immunosupp, 2), ":1 )\n\n")

# Save marker set
write_csv(CURATED_MARKERS, "curated_immune_metabolic_markers.csv")

cat("Marker set saved to: curated_immune_metabolic_markers.csv\n")
cat("Color palettes defined for dual plotting strategy\n\n")

cat("==========================================================\n")
cat("DUAL PLOTTING STRATEGY:\n")
cat("  Plot 1: Color by Cell Type (biological context)\n")
cat("  Plot 2: Color by Immune Function (pro/anti-tumor)\n")
cat("==========================================================\n\n")

# =============================================================================
# EXPORT FOR INTEGRATION
# =============================================================================

# Create simplified list for easy integration
MARKER_LIST <- list(
  markers = CURATED_MARKERS,
  celltype_colors = CELLTYPE_COLORS,
  function_colors = FUNCTION_COLORS
)

saveRDS(MARKER_LIST, "curated_markers_with_colors.rds")

cat("Complete marker set with color palettes saved to:\n")
cat("  - curated_immune_metabolic_markers.csv\n")
cat("  - curated_markers_with_colors.rds\n\n")
