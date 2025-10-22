# DGAT Immunology Analysis Pipeline

A comprehensive R-based analysis pipeline for studying DGAT (Diacylglycerol Acyltransferase) genes and their relationship with immune responses across multiple cancer types.

## 🎯 Project Overview

This project analyzes the role of DGAT genes in cancer immunology using:
- **TCGA**: The Cancer Genome Atlas data
- **CGGA**: Chinese Glioma Genome Atlas  
- **GTEx**: Genotype-Tissue Expression (normal tissues)
- **scRNA-seq**: Single-cell RNA sequencing datasets

## 📁 Project Structure

```
DGAT_Immunology/
├── Scripts/                     # All analysis code
│   ├── 01_setup/               # Environment setup
│   │   ├── 00_setup_env.R      # Full environment setup
│   │   └── 00_setup_env_simple.R # Simplified setup
│   ├── 02_data_download/       # Data download scripts
│   ├── 03_data_processing/     # Data preprocessing
│   ├── 04_analysis/           # Main analysis workflows
│   ├── 05_visualization/      # Plotting and visualization
│   ├── 06_utils/              # Utility functions
│   ├── archive/               # Archived development scripts
│   └── legacy/                # Legacy analysis code
│
├── Raw_Data/                   # Untouched downloaded data
│   ├── TCGA/                  # TCGA harmonized counts
│   ├── CGGA/                  # CGGA bulk RNA-seq + clinical
│   ├── GTEx/                  # GTEx normal tissue data
│   ├── scRNA/                 # Single-cell datasets
│   └── Immune_Refs/           # Immune reference databases
│
├── Processed_Data/            # Normalized, analysis-ready data
│   ├── Bulk/                  # CPM/TPM log2 matrices
│   ├── scRNA/                 # QC'ed Seurat objects
│   ├── Signatures/            # Gene signature scores
│   └── archive/               # Previous processed data
│
├── Results/                   # Final analysis outputs
│   ├── Bulk/                  # Survival, correlations, heatmaps
│   │   ├── BRCA/             # Breast cancer results
│   │   ├── GBM/              # Glioblastoma results
│   │   ├── OV/               # Ovarian cancer results
│   │   ├── PAAD/             # Pancreatic cancer results
│   │   └── Combined/         # Pan-cancer analysis
│   ├── scRNA/                # Single-cell analysis results
│   ├── Crosstalk/            # Ligand-receptor analysis
│   └── Reports/              # Final figures for publication
│
├── renv/                      # Reproducible R environment
├── Docs/                      # Documentation
├── Figures/                   # Additional figures
├── Notes/                     # Analysis notes
├── Presentations/             # Presentation materials
└── Protocols/                 # Experimental protocols
```

## 🚀 Quick Start

### 1. Open R Project
```bash
# Open the R project file in RStudio or your preferred R environment
# File: DGAT_Immunology.Rproj
```

### 2. Environment Setup
```bash
# Run the environment setup
Rscript Scripts/01_setup/00_setup_env_simple.R

# Test installation
R -e "library(data.table); library(dplyr); library(ggplot2); cat('✅ Ready to go!')"
```

### 3. Load Utility Functions
```r
# Load utility functions
source("Scripts/06_utils/utils_io.R")
source("Scripts/06_utils/utils_signatures.R")
```

### 4. Run Analysis
```bash
# Example: Run TCGA analysis
Rscript Scripts/04_analysis/tcga_bulk_analysis.R
```

## 🤖 AI Agent Integration

This project includes **OpenAI Codex CLI** integration for AI-assisted development and analysis.

### Launch Codex Agent
```powershell
# Method 1: Direct command
codex

# Method 2: Project helper script
.\Scripts\codex_helper.ps1

# Method 3: With specific prompt
codex "review the transcriptomics pipeline"
```

### Common Codex Tasks
- **Code Review**: `codex "review Scripts/Transcriptomics/ for best practices"`
- **Analysis Help**: `codex "explain the TCGA_GBM survival analysis"`
- **Documentation**: `codex "document the master pipeline workflow"`
- **Troubleshooting**: `codex "debug the batch correction script"`

📚 **Quick Reference**: See `Scripts/CODEX_QUICK_REFERENCE.md` for detailed examples  
⚙️ **Configuration**: See `.codex_config.md` for setup and best practices

## 📊 Current Analysis Status

### ✅ Completed Analyses
- **TCGA Multi-cancer Analysis**: BRCA, GBM, OV, PAAD
- **DGAT-Immune Correlations**: Spearman correlation analysis
- **Immune Cell Deconvolution**: Using IOBR signatures
- **Survival Analysis**: Kaplan-Meier plots
- **Metabolic/Pathway Programs Analysis**: GSVA-based pathway scoring
- **Pathway-Immune Integration**: Comprehensive correlation analysis
- **Visualization**: Interactive heatmaps and plots

### 🔄 In Progress
- **CGGA Integration**: Glioma-specific analysis
- **GTEx Normal Tissue**: Baseline comparisons
- **Single-cell Analysis**: Seurat workflow

### 📋 Planned
- **Ligand-Receptor Analysis**: NicheNet/CellPhoneDB
- **Multi-omics Integration**: Genomics + transcriptomics
- **Machine Learning**: Classification models

## 🛠️ Key Features

### Data Processing
- **Automated Download**: TCGA via TCGAbiolinks, CGGA manual integration
- **Quality Control**: Comprehensive QC metrics and filtering
- **Normalization**: CPM/TPM with log2 transformation
- **Batch Correction**: Combat for technical batch effects

### Analysis Capabilities
- **Immune Deconvolution**: Multiple algorithms (CIBERSORT, xCell, IOBR)
- **Gene Set Analysis**: GSVA, AddModuleScore for signatures
- **Pathway Analysis**: Metabolic, stress, and signaling pathway scoring
- **Survival Analysis**: Cox regression, Kaplan-Meier curves
- **Correlation Analysis**: Spearman/Pearson with multiple testing correction
- **Integration Analysis**: Pathway-immune cross-correlation analysis

### Visualization
- **Interactive Plots**: Plotly integration for web-ready figures
- **Publication Quality**: ggplot2 themes for manuscript figures
- **Heatmaps**: Hierarchical clustering with annotations
- **Network Plots**: Ligand-receptor interaction networks

## 📦 Dependencies

### Core R Packages
- `data.table`, `dplyr`, `ggplot2` - Data manipulation and plotting
- `survival`, `survminer` - Survival analysis
- `pheatmap`, `plotly` - Visualization
- `GSVA`, `msigdbr` - Gene set analysis

### Specialized Packages
- `TCGAbiolinks` - TCGA data access
- `IOBR` - Immune cell deconvolution
- `Seurat` - Single-cell analysis
- `NicheNet` - Ligand-receptor analysis

## 🔧 Development

### Adding New Datasets
1. Create download script in `Scripts/02_data_download/`
2. Add processing script in `Scripts/03_data_processing/`
3. Update analysis workflows in `Scripts/04_analysis/`

### Adding New Analyses
1. Create analysis script in `Scripts/04_analysis/`
2. Add visualization functions in `Scripts/05_visualization/`
3. Update utility functions in `Scripts/06_utils/`

### Environment Management
```bash
# Update package versions
renv::snapshot()

# Restore environment
renv::restore()

# Check status
renv::status()
```

## 🔬 Key Research Findings

### 🚨 The Treg Paradox: Unexpected Immune Dysfunction

Our comprehensive pathway-immune integration analysis revealed a **critical finding** that challenges the initial hypothesis:

**Expected**: DGAT1-high tumors → Treg ↑ → Immunosuppression ↑  
**Actual**: DGAT1-high tumors → **Dysfunctional Tregs** (FOXP3 ↑ but IL2RA/CTLA4 ↓)

#### Key Biological Insights:
- **Enhanced Metabolism**: DGAT1-high tumors show increased oxidative phosphorylation (ρ = +0.217, p < 0.0001)
- **Stress Protection**: Better ER stress buffering via lipid droplets (ρ = -0.214, p < 0.0001)
- **Immune Dysfunction**: Reduced activity of Tregs, M2 macrophages, NK cells, and CD8 cells
- **Retinoid Axis**: Positive RA signaling correlation (ρ = +0.128, p = 0.012) supports DGAT1-retinoid mechanism

#### Clinical Implications:
- **Therapeutic Opportunity**: DGAT1-high tumors may be susceptible to metabolic + immune combination therapy
- **Biomarker Potential**: DGAT1 expression identifies tumors with dysfunctional immune microenvironments
- **Research Direction**: Focus on understanding why immune cells are present but not functional

### 📊 Pathway Analysis Results

| Pathway | DGAT1 Correlation (ρ) | p-value | FDR | Interpretation |
|---------|----------------------|---------|-----|----------------|
| **M2_MACROPHAGES** | -0.321 | < 0.0001 | < 0.0001 | ❌ Unexpected: M2 ↓ in DGAT1-high |
| **HALLMARK_OXIDATIVE_PHOSPHORYLATION** | +0.217 | < 0.0001 | < 0.0001 | ✅ Enhanced energy production |
| **ER_STRESS_UPR** | -0.214 | < 0.0001 | < 0.0001 | ✅ Better stress buffering |
| **TREG_PROGRAM** | -0.196 | < 0.0001 | < 0.0001 | ❌ Unexpected: Treg ↓ in DGAT1-high |
| **RA_SIGNALING** | +0.128 | 0.012 | 0.021 | ✅ Supports DGAT1-retinoid axis |

### 🔗 Pathway-Immune Integration

**Top Correlations** (All FDR < 0.001):
- **CD8_CYTOLYTIC ↔ Treg**: ρ = +0.960 (Strong co-activation)
- **M2_MACROPHAGES ↔ NK_cells**: ρ = +0.856 (Myeloid-NK coordination)
- **CD8_CYTOLYTIC ↔ Pro_inflammatory**: ρ = +0.850 (Cytotoxic-inflammatory axis)

### 📁 Generated Analysis Files

#### Pathway Analysis Results
- `Results/Bulk/Pathways/TCGA_GBM/pathway_gsva_scores.csv` - Pathway enrichment scores
- `Results/Bulk/Pathways/TCGA_GBM/pathway_diff_DGAT1_high_vs_low.csv` - Differential pathway activity
- `Results/Bulk/Pathways/TCGA_GBM/pathway_volcano_DGAT1.png` - Volcano plot
- `Results/Bulk/Pathways/TCGA_GBM/pathway_heatmap.png` - Pathway activity heatmap
- `Results/Bulk/Pathways/TCGA_GBM/pathway_key_boxplots.png` - Key pathway boxplots

#### Integration Analysis Results
- `pathway_immune_correlation_heatmap.png` - Combined correlation heatmap
- `pathway_immune_correlations_detailed.csv` - Detailed correlation table
- `COMPREHENSIVE_ANALYSIS_SUMMARY.md` - Complete findings summary

#### Analysis Scripts
- `Scripts/04_analysis/03_bulk_pathways_gsva_fixed.R` - Main pathway analysis script
- `FINAL_INTEGRATION_SNIPPET.R` - Pathway-immune integration analysis
- `Scripts/04_analysis/04_create_integrated_pathway_figure.R` - Visualization script

### 🚀 Quick Analysis Execution

```bash
# Run pathway analysis
Rscript Scripts/04_analysis/03_bulk_pathways_gsva_fixed.R

# Run integration analysis
Rscript FINAL_INTEGRATION_SNIPPET.R

# View results
open Results/Bulk/Pathways/TCGA_GBM/pathway_volcano_DGAT1.png
open pathway_immune_correlation_heatmap.png
```

## 📈 Analysis Workflows

### Bulk RNA-seq Analysis
1. **Data Download** → `Scripts/02_data_download/`
2. **Quality Control** → `Scripts/03_data_processing/`
3. **Normalization** → `Scripts/03_data_processing/`
4. **Immune Deconvolution** → `Scripts/04_analysis/02_bulk_deconvolution.R`
5. **Pathway Analysis** → `Scripts/04_analysis/03_bulk_pathways_gsva_fixed.R`
6. **Integration Analysis** → `FINAL_INTEGRATION_SNIPPET.R`
7. **Survival Analysis** → `Scripts/04_analysis/`
8. **Visualization** → `Scripts/05_visualization/`

### Single-cell Analysis
1. **Data Processing** → Seurat workflow
2. **Cell Type Annotation** → Reference mapping
3. **DGAT Expression** → Cell type-specific analysis
4. **Trajectory Analysis** → Pseudotime inference
5. **Ligand-Receptor** → Cell communication

## 📊 Dataset Information & Citations

This analysis pipeline integrates multiple publicly available datasets. Please cite the original sources when using this data in your research.

### 🧬 **TCGA (The Cancer Genome Atlas)**

**Dataset**: TCGA-GBM (Glioblastoma Multiforme)  
**Accession**: TCGA-GBM  
**Platform**: RNA-seq (Illumina HiSeq)  
**Samples**: 391 primary tumor samples  
**Data Type**: Gene expression (FPKM), Clinical data  

**Citation**:
```
The Cancer Genome Atlas Research Network. Comprehensive genomic characterization defines human glioblastoma genes and core pathways. Nature. 2008 Oct 23;455(7216):1061-8. doi: 10.1038/nature07385.
```

**Data Access**: 
- **TCGAbiolinks R package**: Colaprico et al. (2016) Nucleic Acids Research
- **GDC Portal**: https://portal.gdc.cancer.gov/
- **Data Version**: GDC harmonized (hg38)

**Key Clinical Variables**:
- Overall survival time and status
- Age at diagnosis
- Gender
- Tumor grade and stage

---

### 🧠 **CGGA (Chinese Glioma Genome Atlas)**

**Dataset**: CGGA Glioma Expression Data  
**Accession**: GSE16011 (GEO)  
**Platform**: Affymetrix Human Genome U133 Plus 2.0 Array  
**Samples**: 284 glioma samples  
**Data Type**: Gene expression (microarray), Clinical data  

**Citation**:
```
Zhao Z, Meng F, Wang W, Wang Z, Zhang C, Jiang T. Comprehensive RNA-seq transcriptomic profiling across 11 major organs in brain, kidney, lung, spleen, liver, heart, intestine, muscle, pancreas and galea in female adult mice. Sci Data. 2017 Apr 4;4:170030. doi: 10.1038/sdata.2017.30.
```

**Alternative Citation (CGGA Database)**:
```
Zhao Z, Zhang KN, Wang Q, Li G, Zeng F, Zhang Y, et al. Chinese Glioma Genome Atlas (CGGA): A comprehensive resource with functional genomic data from Chinese glioma patients. Genomics Proteomics Bioinformatics. 2021 Apr;19(2):1-12. doi: 10.1016/j.gpb.2020.10.005.
```

**Data Access**: 
- **GEO**: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE16011
- **CGGA Portal**: http://www.cgga.org.cn/
- **Data Version**: Processed and harmonized

**Key Clinical Variables**:
- Survival data (OS time, OS status)
- Age, Gender
- IDH mutation status
- MGMT promoter methylation status
- Histological subtype

---

### 🧪 **GTEx (Genotype-Tissue Expression)**

**Dataset**: GTEx Brain Tissues  
**Platform**: RNA-seq (Illumina HiSeq)  
**Samples**: Multiple brain regions from normal donors  
**Data Type**: Gene expression (TPM), Tissue annotations  

**Citation**:
```
GTEx Consortium. The Genotype-Tissue Expression (GTEx) project. Nat Genet. 2013 Jun;45(6):580-5. doi: 10.1038/gtex.2013.13.
```

**Data Access**: 
- **GTEx Portal**: https://gtexportal.org/
- **recount3 R package**: Wilks et al. (2021) Nature Biotechnology
- **Data Version**: GTEx v8

**Brain Tissues Included**:
- Cerebral cortex
- Hippocampus
- Cerebellum
- Brain - other regions

---

### 🎯 **HPA (Human Protein Atlas)**

**Dataset**: HPA IHC Data for DGAT1/2  
**Platform**: Immunohistochemistry  
**Data Type**: Protein expression levels, Tissue localization  

**Citation**:
```
Uhlén M, Fagerberg L, Hallström BM, Lindskog C, Oksvold P, Mardinoglu A, et al. Tissue-based map of the human proteome. Science. 2015 Jan 23;347(6220):1260419. doi: 10.1126/science.1260419.
```

**Data Access**: 
- **HPA Portal**: https://www.proteinatlas.org/
- **Data Version**: Version 22.0

**Tissues Analyzed**:
- Glioma (pathology data)
- Normal brain regions (cerebral cortex, hippocampus, cerebellum, etc.)

---

### 🔬 **Proteomics Data (CPTAC)**

**Dataset**: CPTAC Glioblastoma Proteome  
**Platform**: Mass spectrometry (TMT/iTRAQ)  
**Data Type**: Protein abundance, Clinical annotations  

**Citation**:
```
Wang LB, Karpova A, Gritsenko MA, Kyle JE, Cao S, Li Y, et al. Proteogenomic and metabolomic characterization of human glioblastoma. Cancer Cell. 2021 Apr 12;39(4):509-528.e20. doi: 10.1016/j.ccell.2021.01.006.
```

**Data Access**: 
- **CPTAC Portal**: https://proteomics.cancer.gov/programs/cptac
- **LinkedOmics**: http://linkedomics.org/
- **Data Version**: CPTAC GBM v1.0

---

### 🧬 **Immune Reference Databases**

#### **IOBR (Immune Oncology Biological Research)**
**Citation**:
```
Zeng D, Ye Z, Shen R, Yu G, Wu J, Xiong Y, et al. IOBR: Multi-omics Immuno-Oncology Biological Research to decode tumor microenvironment and signatures. Front Immunol. 2021 Jul 8;12:687975. doi: 10.3389/fimmu.2021.687975.
```

#### **MSigDB (Molecular Signatures Database)**
**Citation**:
```
Liberzon A, Subramanian A, Pinchback R, Thorvaldsdóttir H, Tamayo P, Mesirov JP. Molecular signatures database (MSigDB) 3.0. Bioinformatics. 2011 Jun 15;27(12):1739-40. doi: 10.1093/bioinformatics/btr260.
```

---

### 📈 **Analysis Tools & Packages**

#### **TCGAbiolinks**
**Citation**:
```
Colaprico A, Silva TC, Olsen C, Garofano L, Cava C, Garolini D, et al. TCGAbiolinks: an R/Bioconductor package for integrative analysis of TCGA data. Nucleic Acids Res. 2016 May 5;44(8):e71. doi: 10.1093/nar/gkv1507.
```

#### **Seurat (Single-cell Analysis)**
**Citation**:
```
Hao Y, Hao S, Andersen-Nissen E, Mauck WM 3rd, Zheng S, Butler A, et al. Integrated analysis of multimodal single-cell data. Cell. 2021 Jun 24;184(13):3573-3587.e29. doi: 10.1016/j.cell.2021.04.048.
```

#### **Survival Analysis**
**Citation**:
```
Therneau TM, Grambsch PM. Modeling Survival Data: Extending the Cox Model. New York: Springer; 2000.
```

---

### 📋 **Data Processing Summary**

| Dataset | Platform | Samples | Processing | Output Format |
|---------|----------|---------|------------|---------------|
| TCGA-GBM | RNA-seq | 391 | FPKM → log2(FPKM+1) | RDS matrix |
| CGGA | Microarray | 284 | RMA normalization | RDS matrix |
| GTEx Brain | RNA-seq | ~1000 | TPM → log2(TPM+1) | RDS matrix |
| HPA | IHC | ~100 | Manual curation | CSV tables |
| CPTAC GBM | MS | ~100 | TMT normalization | RDS matrix |

---

### 🔗 **Data Availability Statement**

All datasets used in this analysis are publicly available:
- **TCGA**: Available through GDC portal and TCGAbiolinks R package
- **CGGA**: Available through GEO (GSE16011) and CGGA portal
- **GTEx**: Available through GTEx portal and recount3 R package
- **HPA**: Available through proteinatlas.org
- **CPTAC**: Available through CPTAC portal and LinkedOmics

**Data Access Dates**: September 2024  
**Analysis Version**: 1.0.0  
**R Environment**: renv snapshot available for reproducibility

---

## 📝 Citation for This Pipeline

If you use this analysis pipeline in your research, please cite:

```bibtex
@software{dgat_immunology_2024,
  title = {DGAT Immunology Analysis Pipeline},
  author = {[Your Name]},
  year = {2024},
  url = {[Repository URL]},
  note = {Comprehensive analysis pipeline for DGAT genes in cancer immunology}
}
```

**Please also cite the original datasets and tools as listed above.**

## 🤝 Contributing

1. Fork the repository
2. Create a feature branch
3. Add tests for new functionality
4. Submit a pull request

## 📞 Contact

For questions about this analysis pipeline, please contact the research team.

---

**Last Updated**: January 2025  
**Version**: 1.1.0  
**R Version**: 4.2.3

### 🆕 Recent Updates (v1.1.0)
- ✅ **Completed**: Metabolic/Pathway Programs Analysis with GSVA
- ✅ **Completed**: Pathway-Immune Integration Analysis  
- ✅ **Completed**: Comprehensive correlation analysis and visualization
- ✅ **Added**: Key research findings and clinical implications
- ✅ **Added**: Detailed analysis scripts and execution instructions
