# Guide to Completing the Single-Cell GBM Analysis Workflow

## Overview

This guide outlines the steps required to complete a full single-cell RNA-seq (scRNA-seq) analysis for glioblastoma (GBM) samples with matched tumour and normal tissue. It reflects the recommended practices from high-impact publications and scRNA-seq best-practice guidelines and aligns with the pipeline implemented in the `sc_pipeline.py` script provided in this repository. The workflow emphasises reproducibility, careful quality control (QC), appropriate batch correction and integrative analyses.

---

## 1. Data Acquisition and Metadata

### Gather data
Obtain the raw 10x Genomics matrices (or equivalent count matrices) for each tumour and normal sample. Ensure you also have any sample-level metadata (patient ID, region, IDH status, treatment history, sequencing chemistry) to enable later statistical modelling.

### Structure sample metadata
Create a table (e.g. a CSV file) summarising key variables for each sample. Variables that may influence gene expression include patient age, sex, tumour region (core vs peritumour), recurrence status and therapy regimen.

### Ensure FAIR compliance
Follow the FAIR principles by keeping raw FASTQ files (if available), processed count matrices and metadata well-organised, assign unique identifiers and plan to deposit the processed data and code in public repositories upon publication.

---

## 2. Pre-processing and Quality Control (QC)

Pre-processing aims to remove low-quality cells while retaining as much biological variability as possible. Best-practice guidelines recommend assessing three QC metrics—the total number of counts per cell, the number of genes detected and the fraction of mitochondrial reads—because cells with low gene detection and high mitochondrial content often reflect broken membranes or dying cells [[sc-best-practices.org](https://sc-best-practices.org)]. However, these covariates should be considered jointly to avoid discarding biologically relevant cells [[sc-best-practices.org](https://sc-best-practices.org)].

### Load data
Use Scanpy's `read_10x_mtx` or a similar function to read each sample and concatenate the tumour and normal matrices. Align gene names across samples using the intersection of gene lists.

### Compute QC metrics
With Scanpy, call `calculate_qc_metrics` to compute `n_genes_by_counts`, `total_counts` and `pct_counts_mt`. Annotate mitochondrial genes via the prefix "MT-" or "mt-" depending on the species [[sc-best-practices.org](https://sc-best-practices.org)].

### Filter low-quality cells
Apply thresholding on the QC metrics. Typical starting thresholds are:

- **Minimum genes per cell**: ≥200
- **Maximum genes per cell**: ≤2,500–7,000 (to remove potential doublets)
- **Maximum mitochondrial fraction**: ≤10%

If tumour cells are known to have elevated mitochondrial content, consider a higher threshold or use median absolute deviation (MAD)-based outlier detection to be more permissive [[sc-best-practices.org](https://sc-best-practices.org)].

Visualise the distributions of QC metrics (e.g. scatter plots of total counts vs. number of genes, coloured by mitochondrial fraction) to inform threshold choice. The single-cell best-practices guide emphasises being permissive and reassessing thresholds after cell annotation to avoid removing rare populations [[sc-best-practices.org](https://sc-best-practices.org)].

### Remove doublets
Use **Scrublet** (Python) or **DoubletFinder** (R) to identify multiplets. Run the algorithm separately on each sample, as doublet rates depend on cell loading. Filter out cells predicted as doublets.

### Filter genes
Remove genes expressed in fewer than a minimum number of cells (e.g. 20) to reduce noise and runtime.

---

## 3. Normalisation and Batch Integration

### Normalise counts
Normalise the filtered count matrix to a constant library size (e.g. 10,000 counts per cell) using `sc.pp.normalize_total` and log-transform with `sc.pp.log1p`.

### Select highly variable genes (HVGs)
Identify HVGs using the Seurat v3 method (implemented in Scanpy's `highly_variable_genes`) with batch correction based on the sample label. Restrict the analysis to the top ~3,000 HVGs.

### Scale the data
Center and scale each gene to unit variance (e.g. `sc.pp.scale` with `max_value=10`).

### Perform principal component analysis (PCA)
Calculate PCs on the scaled HVG matrix (`sc.tl.pca`).

### Batch correction with Harmony
Batch effects arising from sequencing runs or patient differences can obscure biological signals. A recent benchmarking study found that many batch-correction methods (MNN, SCVI, LIGER, Combat, ComBat-seq, BBKNN and Seurat) introduce artifacts, whereas **Harmony consistently performs well and is the only method recommended** for scRNA-seq batch correction [[PMC](https://pmc.ncbi.nlm.nih.gov)]. Use `harmonypy` to integrate samples by specifying the sample column as the batch key and storing the corrected PCs (e.g. `X_pca_harmony`).

---

## 4. Dimensionality Reduction and Clustering

### Neighbour graph
Build a k-nearest neighbour graph using the Harmony-corrected PCs (`sc.pp.neighbors`). Typical values: `n_neighbors=10–30`, `n_pcs=30–50`.

### UMAP/t-SNE embedding
Compute a two- or three-dimensional representation using UMAP (`sc.tl.umap`) for visualisation.

### Clustering
Apply the **Leiden** or **Louvain** algorithm to identify clusters (`sc.tl.leiden`). The resolution parameter controls the granularity of clusters; start with 0.5 and adjust as needed.

---

## 5. Malignancy Inference via CNV

Large-scale copy-number variations (CNVs) distinguish malignant cells from normal cells. Tools such as **infercnvpy**, **CopyKAT** or **sciCNV** infer CNV profiles by comparing expression across consecutive genomic windows relative to a reference set of normal cells. 

### Procedure:
1. Use the normal samples as the reference category
2. Compute smoothed CNV scores (e.g. 100-gene windows)
3. Summarise CNV profiles per cell (variance or segmental deviation)
4. Designate cells with high CNV scores as malignant
5. Adjust the threshold based on the distribution of CNV scores in the normal population
6. Visualise CNV heatmaps to validate the inferred gains and losses

---

## 6. Cell-Type Annotation

### Marker-based annotation
Assign cell types by scoring clusters for canonical marker genes. For GBM, common markers include:

#### **Glial Cells**
- **Astrocytes**: GFAP, AQP4
- **Oligodendrocytes**: MOG, MBP, OLIG1/2
- **OPCs**: PDGFRA, CSPG4
- **Neurons**: SNAP25, NEUN (RBFOX3)

#### **Immune Cells**
- **Microglia**: P2RY12, CX3CR1, TREM2
- **Monocyte-derived macrophages**: CD14, LILRB4, CCL2
- **T cells**: CD3D, CD3E, TRAC
  - Cytotoxic: GZMB, NKG7
  - Regulatory: FOXP3
  - Exhausted: PDCD1, LAG3

### Reference mapping
Optionally, use reference-based annotation tools such as **SingleR**, **scMap** or **Azimuth** to map your clusters to annotated reference datasets (e.g. the Klemm et al. adult brain dataset). Ensure the reference contains both tumour and normal immune populations.

### Cross-validation
Compare the annotation across modalities (RNA, ATAC) and confirm malignant vs. non-malignant status using CNV profiles and marker genes (e.g. EPCAM, SOX2 for tumour cells).

---

## 7. DGAT1 and Lipid-Metabolism Analysis

### Gene expression
After annotation, examine DGAT1 expression across cell types using violin or box plots. Determine whether DGAT1 is enriched in malignant cells, myeloid cells or other populations.

### Pathway scoring
To understand lipid-metabolism context, score cells for gene sets involved in fatty-acid synthesis, elongation and oxidation (e.g. FASN, ACLY, ELOVL2/4/5, CPT1C). Use single-cell gene set variation analysis (GSVA) or AUCell to compute per-cell pathway scores.

### scMetabolism/METAFlux
For a more comprehensive view, employ metabolic inference tools such as **scMetabolism** or **METAFlux**. METAFlux couples gene expression with flux balance analysis to estimate fluxes and can reveal metabolic competition in immune–tumour interactions. Compare metabolic states between tumour and immune cells and evaluate whether DGAT1-high cells exhibit altered fatty-acid metabolism.

### Clinical correlations
If clinical outcomes are available, correlate DGAT1 expression and lipid-metabolism scores with patient survival or therapy response using Cox regression or Kaplan–Meier analysis. Adjust for confounders (age, treatment) using mixed-effects models.

---

## 8. Additional Analyses

### Trajectory inference
Use pseudotime algorithms (**Monocle 3**, **Slingshot**) or RNA velocity (**scVelo**) to model differentiation trajectories of malignant and immune populations. Investigate whether DGAT1 expression or lipid-metabolism scores change along these trajectories.

### Cell–cell communication
Tools like **CellChat**, **LIANA** or **NicheNet** infer ligand–receptor interactions. Focus on:
- Immunosuppressive signalling (TGF-β, IL-6, CSF1, TREM2)
- Metabolic interactions (e.g. APOE, FABP4/5)

Visualise interactions using chord diagrams or heatmaps.

### Spatial context and multi-omics
If spatial transcriptomics or ATAC data are available (e.g. from the Science Advances 2024 dataset), integrate them to map DGAT1 and lipid-metabolism programmes to specific tissue regions. Use **Signac** or **ArchR** to link ATAC peaks to gene expression and identify regulatory elements controlling DGAT1.

### Batch-aware differential expression
When comparing expression or pathway scores between tumour and normal samples, use linear models or mixed-effect models that include patient ID as a random effect to account for inter-patient variability and avoid false positives.

---

## 9. Running the Provided Pipeline

The `sc_pipeline.py` script automates many of the steps above. To run it:

### Install Dependencies
Install the necessary Python packages in your environment (scanpy, scrublet, harmonypy, infercnvpy, seaborn, matplotlib):

```bash
pip install -r requirements.txt
```

### Prepare Data
Prepare your data directories so that `--tumour_path` and `--normal_path` point to the 10x matrices (containing `barcodes.tsv`, `features.tsv` or `genes.tsv`, and `matrix.mtx`).

### Execute the Script
Execute the script from the command line. For example:

```bash
python sc_pipeline.py \
    --tumour_path path/to/tumour/ \
    --normal_path path/to/normal/ \
    --output_dir results/ \
    --min_genes 200 \
    --max_genes 2500 \
    --max_mito 10 \
    --resolution 0.5 \
    --gene_of_interest DGAT1
```

### Inspect Outputs
Inspect the output directory for:
- UMAP plots coloured by sample, cluster and malignancy status
- Violin plot of DGAT1 expression
- `processed_adata.h5ad` containing all QC metrics and annotations

### Customize
Modify the marker gene dictionary in the script (`marker_dict`) to refine cell-type annotation or add additional cell types.

### Downstream Analyses
For downstream analyses (e.g. trajectory inference, cell–cell communication, metabolic scoring), load `processed_adata.h5ad` into Scanpy or other tools and follow the steps described above.

---

## 10. Reproducibility and FAIR Compliance

### Version control
Keep your analysis scripts under version control (e.g. Git) and document all changes. Use tags or releases to mark stable versions of your pipeline.

### Random seeds
Set random seeds for any stochastic processes (e.g. Leiden clustering, Harmony integration) to ensure reproducibility.

### Documentation
Maintain a README describing data sources, software versions and parameters used. Provide citations to relevant literature.

### Data sharing
Upon publication, deposit raw and processed data (count matrices, metadata, AnnData objects) to public repositories (e.g. GEO, EGA) with appropriate access controls. Provide DOIs and metadata to ensure findability and reuse.

---

## Summary

By following these steps and employing the provided pipeline, you will perform a comprehensive, reproducible single-cell analysis of GBM tumour and normal tissues. The workflow emphasises:

- ✅ Careful QC
- ✅ Robust batch correction using Harmony [[PMC](https://pmc.ncbi.nlm.nih.gov)]
- ✅ Reliable CNV-based malignancy inference
- ✅ Thorough downstream analyses

Adapt the thresholds and parameters to your specific data and research questions, and validate key findings with independent datasets or orthogonal assays.

---

## Quick Reference Commands

### Basic Analysis
```bash
# Activate environment
.\venv\Scripts\Activate.ps1

# Run basic pipeline
py scripts\sc_pipeline.py \
    --tumour_path data\raw\tumour\ \
    --normal_path data\raw\normal\ \
    --output_dir results\basic\
```

### Advanced Analysis
```bash
# Run advanced pipeline
py scripts\sc_pipeline_advanced.py \
    --tumour_path data\raw\tumour\ \
    --normal_path data\raw\normal\ \
    --output_dir results\advanced\

# Run downstream analyses
py scripts\downstream_analysis.py \
    --input results\advanced\processed_adata.h5ad \
    --output results\downstream\
```

### Interactive Mode
```bash
# Step-by-step with checkpoints
py scripts\step_by_step.py
```

---

## References

This protocol is based on:

1. **Batch correction**: Harmony benchmarking study [[PMC](https://pmc.ncbi.nlm.nih.gov)]
2. **Quality control**: sc-best-practices.org guidelines
3. **CNV inference**: InferCNVpy methodology
4. **Cell annotation**: Klemm et al., 2020 GBM dataset
5. **Metabolic analysis**: METAFlux and scMetabolism frameworks

For detailed methodology, also see:
- `../Protocols/scRNA_Klemm_Pipeline_Guide.md` - Scientific background
- `PACKAGE_SUMMARY.md` - Script overview
- `README.md` - Technical documentation

---

**Document Version**: 1.0  
**Last Updated**: October 21, 2025  
**Project**: DGAT_Immunology / ScRNA_GBM_Analysis  
**Status**: Complete Analysis Protocol

