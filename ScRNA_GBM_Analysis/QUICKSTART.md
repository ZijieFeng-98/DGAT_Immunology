# Quick Start Guide

## 🚀 Get Started in 5 Minutes

### 1. Install Dependencies

```bash
# Create a new conda environment (recommended)
conda create -n scRNA python=3.9
conda activate scRNA

# Install required packages
pip install -r requirements.txt
```

### 2. Prepare Your Data

Place your 10x Genomics formatted data in the following structure:

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

### 3. Run the Analysis

#### **Option A: Interactive Step-by-Step** (Recommended for first-time users)

```bash
python scripts/step_by_step.py
```

This will show you an interactive menu where you can:
- Run each step individually
- Inspect results between steps
- Adjust parameters on the fly
- Save/load checkpoints

#### **Option B: Automatic Pipeline**

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

#### **Option C: Jupyter Notebook** (For exploration)

```bash
jupyter notebook notebooks/
```

---

## 📊 What Happens in Each Step?

### Step 1: Load Data (1-2 min)
- Reads 10x matrices
- Concatenates tumour + normal
- Creates initial QC metrics

**Output**: `01_loaded.h5ad`

### Step 2: Quality Control (2-3 min)
- Filters low-quality cells
- Removes high MT% cells
- Generates QC plots

**Output**: `02_qc_filtered.h5ad`, QC plots

### Step 3: Doublet Detection (3-5 min)
- Runs Scrublet per sample
- Removes predicted doublets

**Output**: `03_doublets_removed.h5ad`

### Step 4: Normalization & Integration (5-10 min)
- Normalizes counts
- Finds highly variable genes
- Harmony batch correction
- PCA

**Output**: `04_normalized_integrated.h5ad`, integration plots

### Step 5: Clustering (2-3 min)
- Builds KNN graph
- Computes UMAP
- Leiden clustering

**Output**: `05_clustered.h5ad`, UMAP plots

### Step 6: CNV Inference (10-15 min) ⚠️ SLOW
- Infers copy-number variations
- Identifies malignant cells

**Output**: `06_cnv_inferred.h5ad`, CNV plots

### Step 7: Cell-Type Annotation (1-2 min)
- Assigns cell types using markers
- Validates annotations

**Output**: `07_annotated.h5ad`, annotation plots

### Step 8: DGAT1 Analysis (2-3 min)
- Analyzes lipid metabolism genes
- Generates expression plots
- Exports summary tables

**Output**: `08_final.h5ad`, expression plots, CSV files

---

## 🎯 Example Workflows

### Quick Test Run (Small Dataset)

```bash
# Use lower resolution for faster clustering
python scripts/sc_pipeline.py \
    --tumour_path data/raw/tumour/ \
    --normal_path data/raw/normal/ \
    --output_dir results/test_run/ \
    --max_genes 3000 \
    --resolution 0.3
```

### High-Quality Analysis

```bash
# More stringent QC, higher resolution
python scripts/sc_pipeline.py \
    --tumour_path data/raw/tumour/ \
    --normal_path data/raw/normal/ \
    --output_dir results/high_quality/ \
    --min_genes 500 \
    --max_genes 2000 \
    --max_mito 5 \
    --resolution 0.8
```

### Focus on Specific Cell Types

```python
# In Python/Jupyter
from scripts.step_by_step import StepByStepPipeline

pipeline = StepByStepPipeline(
    tumour_path='data/raw/tumour/',
    normal_path='data/raw/normal/',
    output_dir='results/immune_focused/'
)

# Run steps 1-5
pipeline.step1_load_data()
pipeline.step2_quality_control()
pipeline.step3_doublet_detection()
pipeline.step4_normalize_integrate()
pipeline.step5_clustering()

# Load checkpoint and subset to immune cells
pipeline.load_checkpoint('05_clustered')

# Filter to immune cells only
immune_cells = pipeline.adata[pipeline.adata.obs['cell_type'].isin(['Myeloid', 'T cell', 'NK cell'])]
pipeline.adata = immune_cells

# Continue with CNV and analysis
pipeline.step6_cnv_inference()
pipeline.step8_dgat1_analysis()
```

---

## 🔧 Troubleshooting

### Common Issues

**1. Out of Memory Error**
```bash
# Reduce number of HVGs or PCs
python scripts/sc_pipeline.py ... --n_top_genes 2000
```

**2. Too Many/Few Clusters**
```bash
# Adjust resolution (0.1-2.0)
python scripts/sc_pipeline.py ... --resolution 1.0  # More clusters
python scripts/sc_pipeline.py ... --resolution 0.2  # Fewer clusters
```

**3. Poor Batch Correction**
```bash
# Check sample mixing in UMAP
# Adjust Harmony parameters in script if needed
```

**4. CNV Step is Too Slow**
```bash
# Skip CNV if not needed for your analysis
# Or run on a subset of cells first
```

---

## 📁 Output File Structure

After running the pipeline, you'll have:

```
results/
├── processed_adata.h5ad                 # Final AnnData object
├── checkpoints/                         # Intermediate saves
│   ├── 01_loaded.h5ad
│   ├── 02_qc_filtered.h5ad
│   ├── ...
│   └── 08_final.h5ad
├── qc_plots/                           # Quality control plots
│   ├── _01_highest_genes.png
│   ├── _02_qc_violin.png
│   ├── _04_umap_before_harmony.png
│   ├── _05_umap_clusters.png
│   ├── _06_umap_malignancy.png
│   ├── _07_umap_cell_types.png
│   └── _08_umap_genes.png
├── DGAT1_expression_summary.csv        # Gene expression tables
├── DGAT1_violin.png                    # Expression plots
└── cell_counts_summary.csv             # Cell type counts
```

---

## 🎓 Next Steps

After completing the basic pipeline:

1. **Explore Your Data**
   - Load `processed_adata.h5ad` in Python/R
   - Generate custom plots
   - Examine marker genes

2. **Advanced Analyses**
   - Trajectory inference (Monocle, Slingshot)
   - RNA velocity (scVelo)
   - Cell-cell communication (CellChat, LIANA)
   - Metabolic pathway analysis (scMetabolism)

3. **Validation**
   - Compare with bulk RNA-seq (TCGA)
   - Orthogonal validation (flow cytometry, IHC)
   - Cross-dataset validation

4. **Publication**
   - Refer to pipeline guide: `../Protocols/scRNA_Klemm_Pipeline_Guide.md`
   - Follow FAIR principles
   - Share code and data

---

## 💡 Tips

- **Always save checkpoints** - Steps 4-6 take the longest
- **Check QC plots** - Bad QC = bad results
- **Start with default parameters** - Optimize later
- **Use interactive mode** - Best for learning
- **Document parameter changes** - For reproducibility

---

## 📚 Additional Resources

- Full pipeline documentation: `README.md`
- Detailed protocol: `../Protocols/scRNA_Klemm_Pipeline_Guide.md`
- Scanpy tutorials: https://scanpy.tutorials.readthedocs.io/
- sc-best-practices: https://www.sc-best-practices.org/

---

**Need Help?** Check the main README.md or the protocol guide for detailed explanations.

**Ready to start?** Run `python scripts/step_by_step.py` and follow the prompts!

