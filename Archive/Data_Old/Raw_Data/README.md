# Raw_Data Directory Structure

This directory contains all raw, unprocessed datasets used in the DGAT immunology analysis.

## Directory Organization

```
Raw_Data/
├── TCGA/                    # The Cancer Genome Atlas
│   ├── GBM/                 # Glioblastoma (391 samples)
│   ├── BRCA/                # Breast cancer
│   ├── OV/                  # Ovarian cancer
│   └── PAAD/                # Pancreatic cancer
├── CGGA/                    # Chinese Glioma Genome Atlas
│   ├── Expression/          # Gene expression data (284 samples)
│   └── Clinical/            # Clinical annotations
├── GTEx/                    # Genotype-Tissue Expression
│   └── Brain/               # Normal brain tissues
├── HPA/                     # Human Protein Atlas
│   ├── Pathology/           # Cancer tissue IHC data
│   └── Normal_Tissue/       # Normal tissue IHC data
├── Proteome/                # Proteomics data
│   └── CPTAC/               # CPTAC GBM proteome
├── Immune_Refs/             # Immune reference databases
│   ├── MSigDB/              # Gene signatures
│   └── IOBR/                # Immune deconvolution references
├── scRNA/                   # Single-cell RNA-seq data
│   ├── Glioma/              # Glioma single-cell datasets
│   └── Normal_Brain/        # Normal brain single-cell data
└── DATASET_INVENTORY.csv    # File inventory and metadata
```

## Data Sources

- **TCGA**: Downloaded via TCGAbiolinks R package
- **CGGA**: Downloaded from GEO (GSE16011) and CGGA portal
- **GTEx**: Downloaded via recount3 R package
- **HPA**: Downloaded from proteinatlas.org
- **CPTAC**: Downloaded from CPTAC portal

## File Formats

- **Expression data**: .tsv (tab-separated values)
- **Clinical data**: .csv (comma-separated values)
- **IHC data**: .csv (comma-separated values)
- **Proteomics**: .tsv/.csv (mass spectrometry data)

## Data Processing

Raw data in this directory is processed by scripts in `Scripts/03_data_processing/`
and stored in `Processed_Data/` for analysis.

## Last Updated

Generated: 2025-09-29 16:20:46
Pipeline Version: 1.0.0

