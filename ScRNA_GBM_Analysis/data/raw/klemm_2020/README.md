# Klemm et al., 2020 GBM Dataset

## Dataset Information

**Study**: Interrogation of the Microenvironmental Landscape in Brain Tumors Reveals Disease-Specific Alterations of Immune Cells

**Authors**: Klemm et al.  
**Journal**: Nature  
**Year**: 2020  
**DOI**: 10.1038/s41586-020-1959-y  
**GEO Accession**: GSE163108  
**Single Cell Portal**: SCP1290  

## Samples

- **40 IDH-wild-type glioblastomas** (tumor samples)
- **6 adult non-tumor brain specimens** (normal controls)
- **~100,000 cells** total (estimated)

## Data Type

- Single-cell RNA-seq (scRNA-seq)
- 10x Genomics Chromium platform
- Flow cytometry validation
- Bulk RNA-seq for comparison

## Download Instructions

### Option 1: GEO Database
1. Visit: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE163108
2. Download supplementary files
3. Extract to this directory

### Option 2: Single Cell Portal
1. Visit: https://singlecell.broadinstitute.org/single_cell/study/SCP1290
2. Download processed data
3. Extract to this directory

### Option 3: Contact Authors
- Email corresponding author for data access
- May require data use agreement

## Expected Files

After download, you should have:
```
data/raw/klemm_2020/
├── matrix.mtx.gz or .h5ad files
├── features.tsv.gz or genes.tsv
├── barcodes.tsv.gz
└── metadata.csv or similar
```

## Next Steps

After downloading:
1. Convert to 10x format if needed (see conversion scripts)
2. Organize into tumor/ and normal/ directories
3. Run the analysis pipeline

## Citation

If you use this data, please cite:

Klemm F, et al. (2020). Interrogation of the microenvironmental landscape 
in brain tumors reveals disease-specific alterations of immune cells. 
Nature, 584(7820), 120-124. doi: 10.1038/s41586-020-1959-y

## Data Location

Downloaded: 2025-10-21
Path: D:\DGAT_Immunology\ScRNA_GBM_Analysis\data\raw\klemm_2020
