# Data Download Guide - Klemm et al. 2020 GBM Dataset

## 🎯 Dataset Information

**Study**: Klemm et al., "Interrogation of the Microenvironmental Landscape in Brain Tumors Reveals Disease-Specific Alterations of Immune Cells"
- **Journal**: Nature (2020)
- **DOI**: 10.1038/s41586-020-1959-y
- **GEO Accession**: GSE163108
- **Cell Count**: ~100,000 cells (tumour + normal)
- **Samples**: 40 GBM patients + 6 normal brains

---

## 📥 Download Options

### Option 1: GEO Database (NCBI) - Recommended

#### Step 1: Visit GEO
https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE163108

#### Step 2: Download Files
1. Scroll to "Supplementary file" section
2. Look for:
   - Count matrices (.mtx, .txt, or .h5 files)
   - Metadata files (sample information)
   - Feature files (gene names)

3. Download all supplementary files

#### Step 3: Extract Data
```powershell
# Extract compressed files
Expand-Archive -Path GSE163108_*.zip -DestinationPath data\raw\geo_data
```

---

### Option 2: Single Cell Portal (Broad Institute)

#### Step 1: Visit Portal
https://singlecell.broadinstitute.org/single_cell/study/SCP1290

#### Step 2: Download Data
1. Click "Download" tab
2. Select:
   - ☑ Expression matrices
   - ☑ Metadata
   - ☑ Cell annotations

3. Choose format: **10x Genomics** (if available) or **Dense matrix**

#### Step 3: Extract to Project
```powershell
# Move files to project
Move-Item downloaded_files\* data\raw\
```

---

### Option 3: Alternative GBM Datasets

If Klemm dataset is not accessible, use these alternatives:

#### **Neftel et al., 2019** - GSE131928
- 28 GBM samples
- Single-cell RNA-seq
- Well-characterized cellular states

#### **Darmanis et al., 2017** - GSE84465
- Already in your project!
- Location: `D:\DGAT_Immunology\Raw_Data\scRNA\Glioma\Darmanis_GSE84465_scRNA.rds`
- Needs conversion from RDS to 10x format

#### **Richards et al., 2021** - GSE182109
- 17 GBM samples
- Includes spatial transcriptomics

---

## 🔄 Data Format Conversion

### If Data is in 10x Format (Ideal)
Already correct! Files should be:
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

### If Data is in H5AD Format
```python
# Run conversion script
python scripts/convert_h5ad_to_10x.py --input downloaded.h5ad --output data/raw/tumour/
```

### If Data is in RDS Format (R object)
```R
# In R console
library(Seurat)
library(Matrix)

# Load RDS
data <- readRDS("path/to/file.rds")

# Extract tumour cells
tumour <- subset(data, subset = sample_type == "tumour")
writeMM(tumour@assays$RNA@counts, "data/raw/tumour/matrix.mtx")
write.table(rownames(tumour), "data/raw/tumour/features.tsv", quote=F, row.names=F, col.names=F)
write.table(colnames(tumour), "data/raw/tumour/barcodes.tsv", quote=F, row.names=F, col.names=F)

# Compress
gzip data/raw/tumour/*.mtx
gzip data/raw/tumour/*.tsv

# Repeat for normal cells
```

### If Data is in Dense Matrix (CSV/TSV)
```python
# Run conversion script
python scripts/convert_dense_to_10x.py --input matrix.csv --genes genes.txt --cells barcodes.txt --output data/raw/tumour/
```

---

## 🧪 Quick Test with Demo Data

Don't have real data yet? Use our generated demo data:

```powershell
python scripts/sc_pipeline.py \
    --tumour_path data/raw/demo_tumour/ \
    --normal_path data/raw/demo_normal/ \
    --output_dir results/demo/ \
    --max_genes 5000
```

This will run the complete pipeline on synthetic data (~5 minutes).

---

## 📋 Data Organization Checklist

After downloading, verify your structure:

```powershell
# Check tumour data
dir data\raw\tumour

# Should show:
# matrix.mtx.gz (or .mtx)
# features.tsv.gz (or genes.tsv.gz)
# barcodes.tsv.gz

# Check normal data
dir data\raw\normal

# Should show same files
```

### Verify File Integrity

```powershell
# Check if files are readable
python -c "import gzip; f=gzip.open('data/raw/tumour/matrix.mtx.gz','rt'); print(f.readline()); f.close()"

# Should print: %%MatrixMarket matrix coordinate integer general
```

---

## 🔍 Alternative Data Sources

### Public scRNA-seq Databases

1. **Human Cell Atlas**
   - https://data.humancellatlas.org/
   - Search: "glioblastoma" or "brain tumor"

2. **CellxGene**
   - https://cellxgene.cziscience.com/
   - Browse: Cancer > Glioblastoma

3. **UCSC Cell Browser**
   - https://cells.ucsc.edu/
   - Search: GBM datasets

4. **ArrayExpress (EMBL-EBI)**
   - https://www.ebi.ac.uk/arrayexpress/
   - Search: glioblastoma single-cell

---

## 💾 Expected Data Sizes

| Dataset | Compressed | Uncompressed | Cells | Genes |
|---------|------------|--------------|-------|-------|
| Klemm et al. | ~500 MB | ~2 GB | ~100K | ~20K |
| Neftel et al. | ~300 MB | ~1.5 GB | ~24K | ~20K |
| Darmanis et al. | ~100 MB | ~500 MB | ~3.5K | ~20K |
| Demo data | ~5 MB | ~20 MB | 500 | 2K |

---

## 🆘 Troubleshooting

### "Access denied" on GEO
- Some datasets require dbGaP approval for patient privacy
- Apply for access: https://www.ncbi.nlm.nih.gov/gap/
- Or use demo data for testing

### "Files are in wrong format"
- Check file extensions (.mtx.gz, not .mtx.txt)
- Verify matrix market format (first line: %%MatrixMarket)
- Use conversion scripts provided

### "Download is very slow"
- GEO mirrors may be slow
- Try direct download from publication supplements
- Use Single Cell Portal (usually faster)
- Download during off-peak hours

### "Can't find the data files"
- Check publication supplementary materials
- Contact authors directly
- Use alternative datasets listed above

---

## 📧 Getting Help

1. **Check GEO accession page**: Full instructions usually provided
2. **Read publication**: Supplementary methods explain data access
3. **Contact authors**: Email corresponding author for data access
4. **Use alternatives**: Many similar GBM scRNA-seq datasets available

---

## ✅ Data Download Checklist

- [ ] Downloaded tumour sample data
- [ ] Downloaded normal sample data
- [ ] Files in correct 10x format
- [ ] Files compressed (.gz extension)
- [ ] Organized in data/raw/ directories
- [ ] Verified files are readable
- [ ] Ready to run pipeline!

---

**Once data is ready, run:**
```powershell
python scripts\step_by_step.py
```

Or test with demo data first:
```powershell
python scripts\create_demo_data.py
python scripts\sc_pipeline.py --tumour_path data/raw/demo_tumour/ --normal_path data/raw/demo_normal/ --output_dir results/demo/
```

