# Single-Cell Analysis Pipeline for the Klemm et al., 2020 GBM Dataset

## Overview

The Klemm et al. (2020) study used flow cytometry, bulk RNA-seq and limited single-cell RNA-seq to profile immune and stromal populations in 40 IDH-wild-type glioblastomas and six adult non-tumour brain specimens [PMC](https://pmc.ncbi.nlm.nih.gov). Because it uniquely includes true adult normal brain samples alongside tumour tissue, the dataset is ideal for studying immune patterns across tumour and normal compartments, even though the single-cell component is modest. The following reproducible pipeline is designed to analyse the tumour-versus-normal single-cell data, integrate immune and stromal cell types, and interrogate lipid-metabolism genes (e.g., DGAT1) within glioblastoma immunology. The workflow is built around modern standards (FAIR data sharing, transparent code and reproducibility) and synthesises best practices from recent high-impact studies (2023–2025). Each step is accompanied by references that document why the step is important and what parameters are commonly used.

---

## 1. Data Acquisition and Metadata

### Download raw and processed data

- Obtain the raw FASTQ files or count matrices (if available) for all tumour and normal brain specimens. Use vendor pipelines (e.g., 10x Genomics CellRanger) to align reads and generate gene–cell count matrices. For the Klemm dataset the single-cell data may be small; still, re-running mapping ensures reproducibility.

- Collect metadata: patient ID, tissue source (tumour vs non-tumour cortex), age, sex, sampling site, sequencing chemistry, run IDs and library preparation. Accurate metadata enable downstream mixed-model analyses and reproducibility.

- Ensure FAIR compliance by assigning unique identifiers (e.g., DOIs), depositing raw and processed data to repositories (GEO/EGA) and publishing analysis scripts with version control.

### Sample size considerations

- Although the Klemm dataset contains fewer cells than large atlases, include all available tumour and normal specimens. Recent GBM studies recommend analysing ≥10 patients per condition to detect effects across clinical subgroups [[nature.com](https://nature.com)]. The Klemm dataset meets this requirement for tumour samples but has only six normal brains, so caution is warranted when interpreting subtle differences.

### Complementary datasets (optional)

- To increase statistical power or validate findings, integrate additional normal-brain single-cell datasets (e.g., from GSE126836 or GSE140231) or tumour datasets that include paired peritumoral tissue (e.g., the 2024 Science Advances GBM dataset). Integration methods described below can correct batch effects and allow cross-study comparisons.

---

## 2. Pre-processing and Quality Control

### Load count matrices

- Use Python (Scanpy/Anndata) or R (Seurat) to import the gene-cell count matrix for each sample.

- Annotate mitochondrial (MT-), ribosomal (RPS, RPL), and hemoglobin (HB) genes to compute QC metrics [[sc-best-practices.org](https://sc-best-practices.org)].

### Calculate QC metrics

For each cell, compute:

1. **Number of UMI counts** (`total_counts`)
2. **Number of detected genes** (`n_genes_by_counts`)
3. **Fraction of mitochondrial RNA** (`pct_counts_mt`)

These metrics are key for identifying dying cells: low UMI counts and gene numbers combined with high mitochondrial content often indicate cells with broken membranes [[sc-best-practices.org](https://sc-best-practices.org)]. Use `scanpy.pp.calculate_qc_metrics` to derive these metrics.

### Filter low-quality cells

- Apply thresholding **jointly** across QC metrics rather than on a single metric to avoid discarding rare but valid populations [[sc-best-practices.org](https://sc-best-practices.org)].

- As a starting guideline, retain cells with:
  - **200–2,500 detected genes**
  - **300–15,000 UMI counts**
  - **<10% mitochondrial reads** [[biostate.ai](https://biostate.ai)]
  
  These thresholds remove damaged or stressed cells while preserving biological variability.

- Use automatic outlier detection via **Median Absolute Deviation (MAD)** to flag cells deviating more than 3–5 MADs from the median for each metric [[sc-best-practices.org](https://sc-best-practices.org)]. Adjust thresholds after visualising QC metric distributions (violin plots and scatter plots).

- Remove cells with extreme gene counts (>7,000 genes) that likely represent doublets [[nature.com](https://nature.com)].

### Doublet detection

- Implement doublet detection using **Scrublet** (Python) or **DoubletFinder** (R). Scrublet was used by a recent tumour–normal single-cell atlas to remove doublets [[nature.com](https://nature.com)].

- Visualise doublet scores and filter cells above an empirically determined threshold.

### Ambient RNA and background correction

- Use tools such as **SoupX** or **CellBender** to remove ambient RNA contamination. This is especially important for droplet-based scRNA-seq where ambient RNA can inflate gene expression signals. Adjust parameters based on the library complexity.

### Normalisation and feature selection

- Log-normalise UMI counts per cell (e.g., `log1p(10^4 * counts / total_counts)` in Scanpy or Seurat).

- Identify **highly variable genes (HVGs)** per sample. Use HVGs for downstream integration and dimensionality reduction.

---

## 3. Data Integration and Batch Correction

### Concatenate samples

- Combine tumour and normal samples into a single Anndata/Seurat object after QC and HVG selection. Tag each cell with sample ID and condition (tumour vs normal).

### Batch effect correction using Harmony

- A recent benchmarking study comparing eight batch-correction methods found that most algorithms (MNN, SCVI, LIGER, ComBat, BBKNN, Seurat) introduce artifacts, whereas **Harmony is the only method that consistently performs well** [[genome.cshlp.org](https://genome.cshlp.org)].

- Apply Harmony on principal components, using sample or sequencing batch as the batch variable. Integrate tumour and normal cells simultaneously to prevent over-correction.

### Dimensionality reduction and clustering

- Perform **PCA** on the Harmony-corrected data and choose an appropriate number of PCs (e.g., 30–50) based on elbow plots or variance explained.

- Use **UMAP** or **t-SNE** for two-dimensional visualisation.

- Cluster cells using the **Louvain** or **Leiden** algorithm at multiple resolution parameters. Evaluate cluster stability and interpretability.

---

## 4. Cell-Type Annotation and Identification of Malignant Cells

### Annotation using canonical markers and reference mapping

- Examine expression of established marker genes:
  - **CD68** for macrophages
  - **GFAP** for astrocytes
  - **OLIG1/2** for oligodendrocytes
  - **CD3D** for T cells

- Use reference mapping tools (**SingleR**, **Azimuth**) to map cells onto curated brain-immune references.

- Where applicable, cross-validate annotations across multiple methods.

### Distinguish malignant vs. non-malignant cells with CNV inference

- Use **inferCNVpy** (Python) or **InferCNV** (R) to infer large-scale copy-number variations. CNV inference requires specifying a set of "normal" reference cells; immune cells often serve as references [[PMC](https://pmc.ncbi.nlm.nih.gov)].

- InferCNV sorts genes by genomic position and uses a moving average to smooth expression; CNV profiles are computed by subtracting normal reference expression from tumour cells [[PMC](https://pmc.ncbi.nlm.nih.gov)].

- In the tumour–normal meta-atlas study, malignant cells were defined by forming separate CNV clusters and having higher CNV scores compared with normal cells [[nature.com](https://nature.com)]. Apply similar criteria to identify cancer cells in the Klemm dataset.

### Refine cell states

- For malignant cells, classify into:
  - **Neural-progenitor-like (NPC-like)**
  - **Oligodendrocyte-progenitor-like (OPC-like)**
  - **Astrocyte-like (AC-like)**
  - **Mesenchymal-like (MES-like)** states
  
  as described by Neftel et al. These states recapitulate GBM heterogeneity and are enriched but not defined by specific genetic aberrations [[PMC](https://pmc.ncbi.nlm.nih.gov)]. Use published gene signatures for each state.

### Quantify immune cell subtypes

- Distinguish **microglia** versus **monocyte-derived macrophages (TAMs)** using ontogenic markers:
  - **P2RY12** and **TMEM119** for microglia
  - **CD163** and **LYZ** for monocyte-derived macrophages

- Identify lymphocyte subsets (CD8 T cells, CD4 T cells, regulatory T cells, NK cells) and activation/exhaustion states based on canonical markers (e.g., **PDCD1/PD-1**, **TOX**, **LAG3**).

---

## 5. Downstream Analyses

### Differential expression and pathway enrichment

- Perform differential expression analyses between tumour and normal immune cells or between malignant and non-malignant cells. Use non-parametric tests (Wilcoxon rank-sum) and adjust P-values via the Benjamini–Hochberg method.

- Employ **Gene Set Variation Analysis (GSVA)** or **AUCell** to compute pathway scores for lipid metabolism, glycolysis, fatty-acid oxidation and immune signalling pathways.

- Highlight expression of **DGAT1** and other lipid metabolism genes (e.g., **FASN**, **ACLY**, **ELOVL2**). Integrate metabolic gene signatures from prior single-cell metabolic analyses (e.g., Saurty-Seerunghen et al.).

### Trajectory and lineage inference

- Use **Monocle 3** or **Slingshot** to construct pseudotime trajectories for tumour cells and immune populations. For RNA velocity, apply **scVelo** on spliced/unspliced counts to infer directional lineage relationships.

- Investigate whether immune cells transition from cytotoxic to exhausted states or whether TAMs differentiate from infiltrating monocytes.

### Cell–cell communication

- Infer ligand–receptor interactions using **CellChat** or **LIANA**. CellChat calculates the probability of communication between cell types based on differentially expressed transmitter and receiver genes and the law of mass action [[PMC](https://pmc.ncbi.nlm.nih.gov)]; significance is assessed via permutation [[PMC](https://pmc.ncbi.nlm.nih.gov)].

- Evaluate immunosuppressive pathways (e.g., **TGF-β**, **IL-6**, **CSF1**) and lipid-related signalling (e.g., **APOE**, **LPL**, **FABP4/5**). Visualise interactions through chord diagrams or heatmaps.

### Metabolic flux analysis

- Apply tools like **scMetabolism**, **scFEA** or **METAFlux** to estimate metabolic pathway activities at the single-cell level. METAFlux uses network-based flux balance analysis to quantify metabolic heterogeneity; high-impact studies have used it to understand metabolic competition in immune cells.

- Compare **DGAT1-high** versus **DGAT1-low** immune cells for differences in fatty-acid synthesis and oxidation pathways.

### Statistical modelling and clinical association

- Use linear mixed models or generalised linear models to assess differences in cell proportions or gene expression between tumour and normal samples while accounting for patient ID as a random effect.

- Although the Klemm dataset lacks survival data, integrate expression signatures with available clinical variables (e.g., IDH status, sex). Validation in external cohorts (e.g., TCGA bulk GBM) can be performed using Cox regression or Kaplan–Meier analysis.

---

## 6. Validation and Reproducibility

### Orthogonal validation

- Confirm key findings (e.g., high DGAT1 expression in TAMs) using immunohistochemistry or flow cytometry on independent GBM samples. Lipid droplet staining and targeted lipidomics can verify metabolic reprogramming.

### Multi-omic integration

- Incorporate **single-cell ATAC-seq** (if available) to link chromatin accessibility to gene expression. Methods such as **Signac** (R) or **ArchR** (R) can identify putative enhancers controlling DGAT1 and lipid metabolism genes.

- Integrate **spatial transcriptomics** (if available) to map metabolic heterogeneity across tumour regions. For example, the tumour–normal meta-atlas integrated spatial data and found that recurrent GBMs shift toward a mesenchymal phenotype [[nature.com](https://nature.com)].

### Code sharing and FAIR compliance

- Publish analysis scripts (Jupyter notebooks or R Markdown) with explicit package versions and random seeds.

- Document all filtering thresholds, parameters (e.g., Harmony settings, differential expression cut-offs) and software versions.

- Deposit processed matrices, cell-type annotations and metadata in public repositories (GEO) and assign DOIs. Include a README that describes the analysis workflow and how to reproduce each figure.

---

## 7. Summary

This pipeline leverages the Klemm et al., 2020 GBM dataset's unique inclusion of adult normal brain tissue to uncover immune patterns across tumour and normal compartments. It integrates rigorous quality control, state-of-the-art batch correction (Harmony) [[genome.cshlp.org](https://genome.cshlp.org)], CNV-based malignant cell identification [[PMC](https://pmc.ncbi.nlm.nih.gov)], and modern downstream analyses such as trajectory inference, cell–cell communication using CellChat [[PMC](https://pmc.ncbi.nlm.nih.gov)] and metabolic flux estimation. By adhering to FAIR principles and sharing reproducible code, this workflow aligns with expectations of high-impact journals and provides a blueprint for future studies that aim to dissect metabolic–immune interactions—such as DGAT1-dependent lipid metabolism—within glioblastoma and other cancers.

---

## Key Quality Control Parameters Summary

| **Parameter** | **Recommended Value** | **Reference** |
|---------------|----------------------|---------------|
| Detected genes per cell | 200–2,500 | sc-best-practices.org |
| UMI counts per cell | 300–15,000 | biostate.ai |
| Mitochondrial RNA % | <10% | biostate.ai |
| Doublet detection | Scrublet or DoubletFinder | nature.com |
| Batch correction | Harmony | genome.cshlp.org |
| PCA components | 30–50 | Standard practice |
| MAD threshold | 3–5 MADs | sc-best-practices.org |

---

## Recommended Tools

### Python Stack
- **Scanpy** - Single-cell analysis
- **Anndata** - Data structures
- **Scrublet** - Doublet detection
- **inferCNVpy** - CNV inference
- **scVelo** - RNA velocity
- **CellBender** - Ambient RNA removal
- **Harmony** - Batch correction

### R Stack
- **Seurat** - Single-cell analysis
- **DoubletFinder** - Doublet detection
- **InferCNV** - CNV inference
- **Monocle 3** - Trajectory inference
- **Slingshot** - Pseudotime analysis
- **CellChat** - Cell-cell communication
- **SoupX** - Ambient RNA removal
- **Signac/ArchR** - ATAC-seq integration

---

## References

This pipeline is based on best practices from recent high-impact publications (2023-2025) in single-cell genomics, glioblastoma immunology, and metabolic reprogramming. Key methodological references include:

1. **Batch correction benchmarking**: genome.cshlp.org
2. **Quality control standards**: sc-best-practices.org, biostate.ai
3. **CNV inference methods**: pmc.ncbi.nlm.nih.gov
4. **Cell-cell communication**: CellChat methodology (pmc.ncbi.nlm.nih.gov)
5. **GBM single-cell studies**: nature.com, Science Advances datasets
6. **Neftel et al. GBM states**: pmc.ncbi.nlm.nih.gov

---

**Document Version:** 1.0  
**Last Updated:** October 21, 2025  
**Project:** DGAT_Immunology  
**Analysis Focus:** Single-cell immune profiling and lipid metabolism in GBM

