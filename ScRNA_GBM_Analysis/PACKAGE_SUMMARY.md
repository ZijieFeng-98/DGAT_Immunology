# GBM Single-Cell DGAT1 Analysis Pipeline - Complete Package

## 📦 Package Contents

This package contains everything you need to perform comprehensive single-cell RNA-seq analysis of glioblastoma samples with focus on DGAT1 and tumor immunology.

### Core Pipeline
- **sc_glioblastoma_dgat1_pipeline.py** (Main analysis script)
  - Complete workflow from raw counts to annotated cell types
  - DGAT1-focused lipid metabolism analysis
  - CNV-based malignancy inference
  - Best-practice QC, normalization, and batch correction

### Downstream Analysis
- **downstream_analysis.py** (Advanced analyses)
  - Trajectory inference (T cell differentiation/exhaustion)
  - Cell-cell communication networks (LIANA)
  - Advanced metabolic profiling
  - Survival correlation analysis

### Documentation
- **README.md** (Complete documentation)
  - Installation instructions
  - Detailed parameter descriptions
  - Pipeline step explanations
  - Troubleshooting guide
  - Cell type markers reference
  - Output file descriptions

- **QUICKSTART.md** (Quick reference)
  - 5-minute setup guide
  - Copy-paste commands
  - Common workflows
  - Quick fixes for common problems

### Configuration & Setup
- **requirements.txt** (Python dependencies)
  - All required packages with versions
  - Optional packages for advanced features

- **config_example.yaml** (Configuration template)
  - Parameter settings template
  - Easy-to-edit format for batch runs

- **test_installation.py** (Installation verification)
  - Tests all dependencies
  - Verifies scanpy functionality
  - Checks visualization setup

## 🚀 Quick Start (3 Steps)

### 1. Install
```bash
conda create -n gbm_analysis python=3.9 -y
conda activate gbm_analysis
pip install -r requirements.txt
python test_installation.py  # Verify installation
```

### 2. Run Analysis
```bash
python sc_glioblastoma_dgat1_pipeline.py \
    --tumor_paths data/tumor1/ data/tumor2/ \
    --normal_paths data/normal1/ \
    --output_dir results/
```

### 3. Downstream Analyses
```bash
python downstream_analysis.py \
    --adata results/processed_adata.h5ad \
    --output_dir results/downstream/
```

## 📊 What You Get

### Main Pipeline Outputs
```
results/
├── processed_adata.h5ad           # Complete analyzed dataset
├── analysis_summary.csv           # Summary statistics
├── qc_metrics.png                 # Quality control plots
├── umap_overview.png              # Visual overview
├── cell_type_annotation.png       # Cell type assignments
├── cnv_malignancy.png             # Malignancy classification
├── dgat1_expression.png           # DGAT1 patterns
├── lipid_metabolism_umap.png      # Metabolic pathway scores
├── de_tumor_vs_normal.csv         # Differential expression results
└── ...                            # Additional plots and tables
```

### Downstream Analysis Outputs
```
results/downstream/
├── trajectory/
│   ├── paga_trajectory.png        # Cell state transitions
│   ├── tcell_pseudotime.png       # T cell differentiation
│   └── dgat1_along_pseudotime.png # DGAT1 dynamics
├── communication/
│   ├── tumor_immune_interactions.png
│   └── liana_results.csv
├── metabolism/
│   ├── metabolic_profile_heatmap.png
│   └── pathway_correlations.png
└── survival/
    ├── km_dgat1.png               # Kaplan-Meier curves
    └── cox_results.csv            # Survival correlations
```

## 🔬 Scientific Features

### Quality Control
- Multi-metric filtering (genes, counts, mitochondrial %)
- Scrublet doublet detection
- Permissive thresholds following best practices

### Batch Correction
- **Harmony integration** (recommended method based on benchmarking)
- Preserves biological variation
- Handles multiple samples/batches

### Cell Type Annotation
- Marker-based identification
- Comprehensive immune and CNS cell types
- Includes rare populations (OPCs, pericytes, exhausted T cells)

### Malignancy Classification
- CNV-based inference using InferCNVpy
- Uses normal cells as reference
- Distinguishes tumor from normal populations

### DGAT1 Analysis
- Expression across cell types
- Tumor vs. normal comparison
- Correlation with malignancy status

### Lipid Metabolism
- 5 pathway scoring modules:
  - Fatty acid synthesis
  - Fatty acid elongation  
  - Fatty acid oxidation
  - Lipid storage (DGAT1/2, PLINs)
  - Cholesterol metabolism

### Advanced Analyses
- Pseudotime trajectory inference
- Cell-cell communication networks
- Metabolic flux estimation
- Survival correlation (with clinical data)

## 📖 Detailed Documentation

### For First-Time Users
👉 Start with **QUICKSTART.md**
- 5-minute setup
- Basic workflow
- Example commands

### For Detailed Information
👉 Read **README.md**
- Complete parameter guide
- Step-by-step explanations
- Troubleshooting
- Cell type markers
- Literature citations

### For Advanced Users
👉 Explore **downstream_analysis.py**
- Trajectory inference
- Communication analysis
- Custom metabolic scoring
- Integration examples

## 🧪 Key Methods & Citations

### Primary Tools
- **Scanpy**: Wolf et al. (2018) Genome Biology
- **Harmony**: Korsunsky et al. (2019) Nature Methods
- **InferCNVpy**: Virshup et al. (2023) Nature Biotechnology
- **Scrublet**: Wolock et al. (2019) Cell Systems

### Best Practices
- Batch correction benchmarking: Tran et al. (2020) Genome Biology
- Single-cell best practices: Heumos et al. (2023) Nature Reviews Genetics
- QC recommendations: https://www.sc-best-practices.org/

### DGAT1 Biology
- Refer to project files:
  - DGAT1_in_Tumor_Immunity__Experimental_Evidence_and_Research_Gaps.md
  - DGAT1_in_Glioblastoma_Immunology__Comprehensive_Bibliography_and_APA_Citation_Guide.md

## 🎯 Use Cases

### 1. Basic GBM Characterization
Identify cell types, assess malignancy, characterize immune landscape

### 2. DGAT1-Focused Investigation  
Determine DGAT1 expression patterns, correlate with cell states and outcomes

### 3. Immunometabolism Study
Compare metabolic profiles between tumor and immune cells

### 4. Trajectory Analysis
Track T cell exhaustion or myeloid polarization over pseudotime

### 5. Therapeutic Target Discovery
Identify communication networks and metabolic vulnerabilities

### 6. Integration Study
Combine with public GBM datasets or spatial transcriptomics

## ⚙️ Customization

### Adjust Parameters
Edit command-line arguments or modify config_example.yaml

### Add Cell Types
Edit `marker_dict` in sc_glioblastoma_dgat1_pipeline.py

### Custom Pathways
Edit `lipid_genes` dictionary to add metabolic gene sets

### Subset Analysis
Load processed_adata.h5ad and filter to cell types of interest

## 🔧 System Requirements

### Minimum
- 16 GB RAM
- 4 CPU cores
- 10 GB disk space

### Recommended
- 32 GB+ RAM (for large datasets)
- 8+ CPU cores
- 50 GB disk space (for results)

### Runtime
- Small dataset (5K cells): ~30 minutes
- Medium dataset (50K cells): ~2 hours
- Large dataset (200K cells): ~6-8 hours

## 📝 Reproducibility

All analyses are reproducible with:
1. Fixed random seeds (default: 42)
2. Version-controlled code
3. Documented parameters
4. FAIR-compliant data management

## 🆘 Getting Help

1. **Check documentation**
   - QUICKSTART.md for basics
   - README.md for details
   - Test with test_installation.py

2. **Common issues**
   - Memory errors → reduce n_top_genes, n_pcs
   - Poor clustering → adjust resolution
   - Batch effects → verify Harmony installation

3. **Resources**
   - Scanpy docs: https://scanpy.readthedocs.io/
   - Best practices: https://www.sc-best-practices.org/
   - Harmony: https://portals.broadinstitute.org/harmony

## 📄 License & Citation

This pipeline is provided for research use. When publishing results:

1. **Cite the tools used** (see README.md for full citations)
2. **Cite Harmony** for batch correction (Korsunsky et al. 2019)
3. **Follow best practices** referenced in documentation
4. **Share code and data** following FAIR principles

## 🎉 You're Ready!

This complete package provides:
✓ Publication-quality single-cell analysis pipeline
✓ DGAT1-focused immunometabolism investigation  
✓ Best-practice workflows with proper citations
✓ Comprehensive documentation
✓ Downstream analysis tools
✓ Reproducible, FAIR-compliant methods

**Next step:** Run `python test_installation.py` to verify your setup!

---

Package created: 2025
For GBM DGAT1 tumor immunology research
