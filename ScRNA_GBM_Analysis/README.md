# Single-Cell RNA-seq Analysis for GBM (Klemm Dataset)

## Project Overview

This project implements a reproducible single-cell RNA-seq analysis pipeline for glioblastoma (GBM) samples, focusing on the **Klemm et al., 2020** dataset which uniquely includes both tumour and adult normal brain tissue. The pipeline is designed to study immune patterns and lipid metabolism (particularly **DGAT1**) in the GBM tumor microenvironment.

## Key Features

- **Comprehensive QC**: Quality control with MAD-based outlier detection
- **Doublet Detection**: Scrublet-based doublet removal
- **Batch Correction**: Harmony integration (best-performing method)
- **CNV Inference**: InferCNVpy for malignant cell identification
- **Cell-Type Annotation**: Marker-based and reference mapping
- **DGAT1 Analysis**: Lipid metabolism gene expression profiling
- **Reproducible**: FAIR-compliant with documented parameters

## Project Structure

```
ScRNA_GBM_Analysis/
├── scripts/              # Analysis scripts
│   ├── sc_pipeline.py    # Main pipeline
│   ├── step_by_step.py   # Interactive step-by-step runner
│   └── utils.py          # Helper functions
├── data/
│   ├── raw/             # Raw count matrices (10x format)
│   └── processed/       # Processed AnnData objects
├── results/             # Analysis outputs and tables
├── figures/             # Generated plots and visualizations
├── logs/                # Pipeline logs
├── notebooks/           # Jupyter notebooks for exploration
├── README.md            # This file
└── requirements.txt     # Python dependencies

```

## Prerequisites

### Python Environment

Python 3.8+ with the following packages:

- scanpy >= 1.9.0
- scrublet >= 0.2.3
- harmonypy >= 0.0.9
- infercnvpy >= 0.4.0
- numpy >= 1.21.0
- pandas >= 1.3.0
- matplotlib >= 3.5.0
- seaborn >= 0.11.0

Install all dependencies:
```bash
pip install -r requirements.txt
```

## Data Requirements

The pipeline expects **10x Genomics formatted** data:
- `matrix.mtx.gz` or `matrix.mtx`
- `features.tsv.gz` or `genes.tsv`
- `barcodes.tsv.gz` or `barcodes.tsv`

Organize data as:
```
data/raw/
├── tumour/
│   ├── matrix.mtx.gz
│   ├── features.tsv.gz
│   └── barcodes.tsv.gz
└── normal/
    ├── matrix.mtx.gz
    ├── features.tsv.gz
    └── barcodes.tsv.gz
```

## Usage

### Option 1: Run Complete Pipeline

```bash
python scripts/sc_pipeline.py \
    --tumour_path data/raw/tumour/ \
    --normal_path data/raw/normal/ \
    --output_dir results/ \
    --min_genes 200 \
    --max_genes 2500 \
    --max_mito 10 \
    --resolution 0.5
```

### Option 2: Step-by-Step Interactive Mode

```bash
python scripts/step_by_step.py
```

This allows you to:
- Run each analysis step individually
- Inspect intermediate results
- Adjust parameters between steps
- Generate QC plots after each step

## Pipeline Steps

### 1. Data Loading
- Load tumour and normal 10x matrices
- Concatenate samples with metadata
- Harmonize gene names

### 2. Quality Control
- Calculate QC metrics (n_counts, n_genes, pct_mt)
- Filter cells: 200-2500 genes, <10% MT
- Remove low-quality cells

### 3. Doublet Detection
- Run Scrublet per sample
- Remove predicted doublets

### 4. Normalization & Integration
- Normalize to 10,000 counts per cell
- Log-transform
- Identify 3000 highly variable genes
- Scale data
- PCA (50 components)
- **Harmony** batch correction

### 5. Clustering
- KNN graph construction
- UMAP embedding
- Leiden clustering

### 6. CNV Inference
- InferCNVpy analysis
- Identify malignant cells
- Calculate CNV scores

### 7. Cell-Type Annotation
- Marker-based annotation
- Cell types: Myeloid, T cell, Astrocyte, Oligodendrocyte, Endothelial
- Optional: Reference mapping with SingleR/Azimuth

### 8. DGAT1 Expression Analysis
- Compute DGAT1 expression per cell type
- Compare tumour vs normal
- Generate violin plots

## Output Files

### Data Objects
- `processed_adata.h5ad` - Final annotated AnnData object

### Figures
- `umap_sample.png` - UMAP colored by sample
- `umap_clusters.png` - UMAP colored by clusters
- `umap_malignancy.png` - UMAP colored by malignancy status
- `DGAT1_violin.png` - DGAT1 expression by cell type

### Results Tables
- `cell_counts_summary.csv` - Cell counts per sample/cluster
- `dgat1_expression_summary.csv` - DGAT1 stats by cell type
- `qc_metrics.csv` - Quality control metrics

## Quality Control Thresholds

Based on best practices from sc-best-practices.org and recent publications:

| Parameter | Default Value | Rationale |
|-----------|--------------|-----------|
| Min genes | 200 | Remove empty droplets |
| Max genes | 2500 | Remove likely doublets |
| Max MT% | 10% | Remove dying cells |
| Doublet rate | 5% | Typical for 10x |
| HVG count | 3000 | Balance noise vs signal |
| PCA components | 50 | Capture biological variation |
| Leiden resolution | 0.5 | Balanced cluster granularity |

## Customization

### Modify Cell-Type Markers

Edit the marker dictionary in `sc_pipeline.py`:

```python
marker_dict = {
    'Myeloid': ['PTPRC', 'CSF1R', 'ITGAM', 'CD68', 'CD163'],
    'T cell': ['CD3D', 'CD3E', 'TRAC', 'CD8A'],
    'Microglia': ['P2RY12', 'TMEM119', 'CX3CR1'],
    # Add more cell types...
}
```

### Add Lipid Metabolism Genes

Extend the analysis to include additional genes:

```python
lipid_genes = ['DGAT1', 'DGAT2', 'FASN', 'ACLY', 'ELOVL2', 'FABP4', 'FABP5']
```

## Validation

- Orthogonal validation with IHC/flow cytometry
- Integration with TCGA bulk RNA-seq
- Cross-validation with independent GBM datasets

## References

1. **Klemm et al., 2020** - Original dataset (Nature, PMID: 32066951)
2. **Harmony** - Best batch correction method (Genome Biology)
3. **InferCNVpy** - CNV inference documentation
4. **sc-best-practices.org** - Quality control guidelines
5. **DGAT_Immunology Pipeline Guide** - See `Protocols/scRNA_Klemm_Pipeline_Guide.md`

## Troubleshooting

### Common Issues

**Issue**: Out of memory during Harmony
- **Solution**: Reduce `n_pcs` or process samples separately

**Issue**: Too few/many clusters
- **Solution**: Adjust `resolution` parameter (0.1-2.0)

**Issue**: High doublet rate
- **Solution**: Check `expected_doublet_rate`, visualize doublet scores

**Issue**: Poor batch correction
- **Solution**: Try different `use_rep` or increase Harmony iterations

## Contact

For questions about this pipeline, refer to:
- Main project: `DGAT_Immunology`
- Protocol guide: `Protocols/scRNA_Klemm_Pipeline_Guide.md`

## Citation

If you use this pipeline, please cite:
- Klemm et al., 2020 (original dataset)
- Scanpy, Harmony, Scrublet, InferCNVpy (method papers)
- Your analysis code repository

---

**Project Status**: Setup Complete  
**Last Updated**: October 21, 2025  
**Version**: 1.0

