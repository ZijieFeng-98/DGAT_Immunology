#!/usr/bin/env python3
"""
Single-Cell Glioblastoma (GBM) Analysis Pipeline with DGAT1 Focus
==================================================================

This script implements a comprehensive scRNA-seq analysis workflow for GBM
samples with matched tumor and normal tissue, focusing on DGAT1 and lipid
metabolism analysis. It follows best practices for batch correction, 
malignancy inference, and cell-type annotation.

Author: Generated for GBM DGAT1 Analysis
Date: 2025
"""

import os
import argparse
import warnings
import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
from typing import Dict, List, Optional, Tuple
import logging

# Suppress warnings
warnings.filterwarnings('ignore')

# Set up logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)


class GBMSingleCellPipeline:
    """
    Complete single-cell RNA-seq analysis pipeline for GBM with DGAT1 focus.
    """
    
    def __init__(
        self,
        tumor_paths: List[str],
        normal_paths: List[str],
        output_dir: str,
        random_seed: int = 42
    ):
        """
        Initialize the pipeline.
        
        Parameters
        ----------
        tumor_paths : list of str
            Paths to tumor sample 10x directories
        normal_paths : list of str
            Paths to normal sample 10x directories
        output_dir : str
            Directory for output files
        random_seed : int
            Random seed for reproducibility
        """
        self.tumor_paths = tumor_paths
        self.normal_paths = normal_paths
        self.output_dir = Path(output_dir)
        self.random_seed = random_seed
        
        # Create output directory
        self.output_dir.mkdir(parents=True, exist_ok=True)
        
        # Set random seeds
        np.random.seed(random_seed)
        
        # Cell type markers
        self.marker_dict = {
            'Astrocytes': ['GFAP', 'AQP4', 'SLC1A2'],
            'Oligodendrocytes': ['MOG', 'MBP', 'OLIG1', 'OLIG2'],
            'OPCs': ['PDGFRA', 'CSPG4', 'VCAN'],
            'Neurons': ['SNAP25', 'RBFOX3', 'SYT1'],
            'Microglia': ['P2RY12', 'CX3CR1', 'TREM2', 'TMEM119'],
            'Macrophages': ['CD14', 'LILRB4', 'CCL2', 'CD68'],
            'T_cells': ['CD3D', 'CD3E', 'TRAC'],
            'Cytotoxic_T': ['GZMB', 'NKG7', 'PRF1'],
            'Regulatory_T': ['FOXP3', 'IL2RA', 'CTLA4'],
            'Exhausted_T': ['PDCD1', 'LAG3', 'HAVCR2'],
            'NK_cells': ['NCAM1', 'KLRB1', 'NKG7'],
            'B_cells': ['CD79A', 'CD79B', 'MS4A1'],
            'Endothelial': ['PECAM1', 'VWF', 'CDH5'],
            'Pericytes': ['RGS5', 'PDGFRB', 'ACTA2'],
            'Malignant': ['SOX2', 'OLIG2', 'EGFR']
        }
        
        # Lipid metabolism gene sets
        self.lipid_genes = {
            'FA_synthesis': ['FASN', 'ACLY', 'ACACA', 'SCD'],
            'FA_elongation': ['ELOVL2', 'ELOVL4', 'ELOVL5', 'ELOVL6'],
            'FA_oxidation': ['CPT1A', 'CPT1C', 'ACOX1', 'HADHA'],
            'Lipid_storage': ['DGAT1', 'DGAT2', 'PLIN2', 'PLIN3'],
            'Cholesterol': ['HMGCR', 'SQLE', 'DHCR7', 'LDLR']
        }
        
        self.adata = None
        
    def load_data(self) -> sc.AnnData:
        """
        Load and concatenate tumor and normal samples.
        
        Returns
        -------
        adata : AnnData
            Concatenated AnnData object
        """
        logger.info("Loading data...")
        
        adatas = []
        
        # Load tumor samples
        for i, path in enumerate(self.tumor_paths):
            logger.info(f"Loading tumor sample {i+1}: {path}")
            adata_temp = sc.read_10x_mtx(path, var_names='gene_symbols', cache=True)
            adata_temp.obs['sample_id'] = f'tumor_{i+1}'
            adata_temp.obs['tissue_type'] = 'tumor'
            adata_temp.obs['patient_id'] = f'patient_{i+1}'
            adatas.append(adata_temp)
        
        # Load normal samples
        for i, path in enumerate(self.normal_paths):
            logger.info(f"Loading normal sample {i+1}: {path}")
            adata_temp = sc.read_10x_mtx(path, var_names='gene_symbols', cache=True)
            adata_temp.obs['sample_id'] = f'normal_{i+1}'
            adata_temp.obs['tissue_type'] = 'normal'
            adata_temp.obs['patient_id'] = f'patient_{len(self.tumor_paths)+i+1}'
            adatas.append(adata_temp)
        
        # Concatenate all samples
        logger.info("Concatenating samples...")
        self.adata = sc.concat(adatas, join='inner', index_unique='_')
        
        # Make variable names unique
        self.adata.var_names_make_unique()
        
        logger.info(f"Loaded {self.adata.n_obs} cells and {self.adata.n_vars} genes")
        
        return self.adata
    
    def quality_control(
        self,
        min_genes: int = 200,
        max_genes: int = 6000,
        max_mito: float = 15.0,
        min_cells: int = 20
    ) -> sc.AnnData:
        """
        Perform quality control filtering.
        
        Parameters
        ----------
        min_genes : int
            Minimum number of genes per cell
        max_genes : int
            Maximum number of genes per cell
        max_mito : float
            Maximum mitochondrial fraction (%)
        min_cells : int
            Minimum cells expressing a gene
            
        Returns
        -------
        adata : AnnData
            Filtered AnnData object
        """
        logger.info("Performing quality control...")
        
        # Identify mitochondrial genes
        self.adata.var['mt'] = self.adata.var_names.str.startswith(('MT-', 'mt-'))
        
        # Calculate QC metrics
        sc.pp.calculate_qc_metrics(
            self.adata,
            qc_vars=['mt'],
            percent_top=None,
            log1p=False,
            inplace=True
        )
        
        # Plot QC metrics before filtering
        fig, axes = plt.subplots(2, 2, figsize=(12, 10))
        
        sc.pl.violin(
            self.adata, ['n_genes_by_counts', 'total_counts', 'pct_counts_mt'],
            jitter=0.4, multi_panel=True, ax=axes[0, 0], show=False
        )
        
        sc.pl.scatter(
            self.adata, x='total_counts', y='n_genes_by_counts',
            color='pct_counts_mt', ax=axes[0, 1], show=False
        )
        
        # Store initial cell count
        n_cells_before = self.adata.n_obs
        
        # Apply filters
        logger.info(f"Filtering cells: min_genes={min_genes}, max_genes={max_genes}, max_mito={max_mito}%")
        
        sc.pp.filter_cells(self.adata, min_genes=min_genes)
        self.adata = self.adata[self.adata.obs['n_genes_by_counts'] < max_genes, :]
        self.adata = self.adata[self.adata.obs['pct_counts_mt'] < max_mito, :]
        
        logger.info(f"Filtering genes: min_cells={min_cells}")
        sc.pp.filter_genes(self.adata, min_cells=min_cells)
        
        # Plot after filtering
        sc.pl.violin(
            self.adata, ['n_genes_by_counts', 'total_counts', 'pct_counts_mt'],
            jitter=0.4, multi_panel=True, ax=axes[1, 0], show=False
        )
        
        sc.pl.scatter(
            self.adata, x='total_counts', y='n_genes_by_counts',
            color='pct_counts_mt', ax=axes[1, 1], show=False
        )
        
        plt.tight_layout()
        plt.savefig(self.output_dir / 'qc_metrics.png', dpi=300, bbox_inches='tight')
        plt.close()
        
        n_cells_after = self.adata.n_obs
        logger.info(f"Removed {n_cells_before - n_cells_after} cells ({100*(n_cells_before - n_cells_after)/n_cells_before:.1f}%)")
        logger.info(f"Retained {n_cells_after} cells and {self.adata.n_vars} genes")
        
        return self.adata
    
    def remove_doublets(self) -> sc.AnnData:
        """
        Remove doublets using Scrublet.
        
        Returns
        -------
        adata : AnnData
            AnnData with doublets removed
        """
        logger.info("Identifying and removing doublets with Scrublet...")
        
        try:
            import scrublet as scr
            
            # Run Scrublet on each sample separately
            doublet_scores = []
            predicted_doublets = []
            
            for sample in self.adata.obs['sample_id'].unique():
                logger.info(f"Processing sample: {sample}")
                
                sample_mask = self.adata.obs['sample_id'] == sample
                adata_sample = self.adata[sample_mask, :]
                
                scrub = scr.Scrublet(adata_sample.X, random_state=self.random_seed)
                scores, preds = scrub.scrub_doublets(
                    min_counts=2,
                    min_cells=3,
                    min_gene_variability_pctl=85,
                    n_prin_comps=30
                )
                
                doublet_scores.extend(scores)
                predicted_doublets.extend(preds)
            
            self.adata.obs['doublet_score'] = doublet_scores
            self.adata.obs['predicted_doublet'] = predicted_doublets
            
            n_doublets = sum(predicted_doublets)
            logger.info(f"Identified {n_doublets} doublets ({100*n_doublets/len(predicted_doublets):.1f}%)")
            
            # Filter doublets
            self.adata = self.adata[~self.adata.obs['predicted_doublet'], :]
            
        except ImportError:
            logger.warning("Scrublet not installed. Skipping doublet removal.")
        
        return self.adata
    
    def normalize_and_process(self, n_top_genes: int = 3000) -> sc.AnnData:
        """
        Normalize, identify HVGs, and scale data.
        
        Parameters
        ----------
        n_top_genes : int
            Number of highly variable genes to select
            
        Returns
        -------
        adata : AnnData
            Processed AnnData object
        """
        logger.info("Normalizing and processing data...")
        
        # Store raw counts
        self.adata.layers['counts'] = self.adata.X.copy()
        
        # Normalize to 10,000 counts per cell
        sc.pp.normalize_total(self.adata, target_sum=1e4)
        
        # Log transform
        sc.pp.log1p(self.adata)
        
        # Store normalized data
        self.adata.layers['log_normalized'] = self.adata.X.copy()
        
        # Identify highly variable genes
        logger.info(f"Selecting {n_top_genes} highly variable genes...")
        sc.pp.highly_variable_genes(
            self.adata,
            n_top_genes=n_top_genes,
            batch_key='sample_id',
            flavor='seurat_v3',
            subset=False
        )
        
        # Plot HVGs
        sc.pl.highly_variable_genes(self.adata, show=False)
        plt.savefig(self.output_dir / 'highly_variable_genes.png', dpi=300, bbox_inches='tight')
        plt.close()
        
        # Scale data
        logger.info("Scaling data...")
        sc.pp.scale(self.adata, max_value=10)
        
        return self.adata
    
    def batch_correction(self, n_pcs: int = 50) -> sc.AnnData:
        """
        Perform PCA and batch correction with Harmony.
        
        Parameters
        ----------
        n_pcs : int
            Number of principal components
            
        Returns
        -------
        adata : AnnData
            AnnData with batch-corrected PCs
        """
        logger.info("Performing PCA...")
        
        # PCA on HVGs
        sc.tl.pca(
            self.adata,
            n_comps=n_pcs,
            use_highly_variable=True,
            svd_solver='arpack',
            random_state=self.random_seed
        )
        
        # Plot PCA variance
        sc.pl.pca_variance_ratio(self.adata, n_pcs=50, show=False)
        plt.savefig(self.output_dir / 'pca_variance.png', dpi=300, bbox_inches='tight')
        plt.close()
        
        # Harmony batch correction
        logger.info("Running Harmony batch correction...")
        
        try:
            import harmonypy as hm
            
            # Run Harmony
            ho = hm.run_harmony(
                self.adata.obsm['X_pca'],
                self.adata.obs,
                'sample_id',
                random_state=self.random_seed
            )
            
            # Store corrected PCs
            self.adata.obsm['X_pca_harmony'] = ho.Z_corr.T
            
            logger.info("Harmony batch correction completed")
            
        except ImportError:
            logger.warning("harmonypy not installed. Using uncorrected PCs.")
            self.adata.obsm['X_pca_harmony'] = self.adata.obsm['X_pca'].copy()
        
        return self.adata
    
    def clustering(
        self,
        n_neighbors: int = 15,
        n_pcs: int = 40,
        resolution: float = 0.5
    ) -> sc.AnnData:
        """
        Build neighbor graph, compute UMAP, and cluster cells.
        
        Parameters
        ----------
        n_neighbors : int
            Number of neighbors for KNN graph
        n_pcs : int
            Number of PCs to use
        resolution : float
            Leiden clustering resolution
            
        Returns
        -------
        adata : AnnData
            AnnData with clustering results
        """
        logger.info("Computing neighbor graph and UMAP...")
        
        # Build neighbor graph on Harmony-corrected PCs
        sc.pp.neighbors(
            self.adata,
            n_neighbors=n_neighbors,
            n_pcs=n_pcs,
            use_rep='X_pca_harmony',
            random_state=self.random_seed
        )
        
        # UMAP
        sc.tl.umap(self.adata, random_state=self.random_seed)
        
        # Leiden clustering
        logger.info(f"Leiden clustering with resolution={resolution}...")
        sc.tl.leiden(
            self.adata,
            resolution=resolution,
            random_state=self.random_seed
        )
        
        # Plot UMAP
        fig, axes = plt.subplots(1, 3, figsize=(18, 5))
        
        sc.pl.umap(self.adata, color='leiden', ax=axes[0], show=False, title='Leiden Clusters')
        sc.pl.umap(self.adata, color='sample_id', ax=axes[1], show=False, title='Sample ID')
        sc.pl.umap(self.adata, color='tissue_type', ax=axes[2], show=False, title='Tissue Type')
        
        plt.tight_layout()
        plt.savefig(self.output_dir / 'umap_overview.png', dpi=300, bbox_inches='tight')
        plt.close()
        
        return self.adata
    
    def infer_malignancy(self, cnv_window: int = 100) -> sc.AnnData:
        """
        Infer malignancy using CNV profiles.
        
        Parameters
        ----------
        cnv_window : int
            Window size for CNV smoothing (number of genes)
            
        Returns
        -------
        adata : AnnData
            AnnData with malignancy annotations
        """
        logger.info("Inferring malignancy using CNV profiles...")
        
        try:
            import infercnvpy as cnv
            
            # Use normal cells as reference
            normal_cells = self.adata.obs['tissue_type'] == 'normal'
            
            if normal_cells.sum() == 0:
                logger.warning("No normal cells found. Skipping CNV inference.")
                self.adata.obs['malignant_cnv'] = 'unknown'
                return self.adata
            
            # Infer CNV
            logger.info("Running infercnvpy...")
            cnv.tl.infercnv(
                self.adata,
                reference_key='tissue_type',
                reference_cat='normal',
                window_size=cnv_window,
                step=10
            )
            
            # Calculate CNV score (variance across genome)
            cnv_scores = np.var(self.adata.obsm['X_cnv'], axis=1)
            self.adata.obs['cnv_score'] = cnv_scores
            
            # Determine threshold (e.g., 95th percentile of normal cells)
            normal_cnv = cnv_scores[normal_cells]
            threshold = np.percentile(normal_cnv, 95)
            
            # Classify cells
            self.adata.obs['malignant_cnv'] = (cnv_scores > threshold).map({
                True: 'malignant',
                False: 'non-malignant'
            })
            
            logger.info(f"CNV threshold: {threshold:.4f}")
            logger.info(f"Malignant cells: {(self.adata.obs['malignant_cnv']=='malignant').sum()}")
            
            # Plot CNV scores
            fig, axes = plt.subplots(1, 2, figsize=(12, 5))
            
            sc.pl.umap(self.adata, color='cnv_score', ax=axes[0], show=False, title='CNV Score')
            sc.pl.umap(self.adata, color='malignant_cnv', ax=axes[1], show=False, title='Malignancy')
            
            plt.tight_layout()
            plt.savefig(self.output_dir / 'cnv_malignancy.png', dpi=300, bbox_inches='tight')
            plt.close()
            
        except ImportError:
            logger.warning("infercnvpy not installed. Skipping CNV inference.")
            self.adata.obs['malignant_cnv'] = 'unknown'
        
        return self.adata
    
    def annotate_cell_types(self) -> sc.AnnData:
        """
        Annotate cell types using marker genes.
        
        Returns
        -------
        adata : AnnData
            AnnData with cell type annotations
        """
        logger.info("Annotating cell types...")
        
        # Calculate marker gene scores for each cell type
        for cell_type, markers in self.marker_dict.items():
            # Filter markers that exist in dataset
            available_markers = [g for g in markers if g in self.adata.var_names]
            
            if len(available_markers) > 0:
                sc.tl.score_genes(
                    self.adata,
                    available_markers,
                    score_name=f'{cell_type}_score',
                    use_raw=False
                )
        
        # Rank clusters by marker scores
        cluster_annotations = {}
        
        for cluster in self.adata.obs['leiden'].unique():
            cluster_cells = self.adata.obs['leiden'] == cluster
            
            # Calculate mean score for each cell type
            scores = {}
            for cell_type in self.marker_dict.keys():
                score_col = f'{cell_type}_score'
                if score_col in self.adata.obs.columns:
                    scores[cell_type] = self.adata.obs.loc[cluster_cells, score_col].mean()
            
            # Assign cell type with highest score
            if scores:
                best_type = max(scores, key=scores.get)
                cluster_annotations[cluster] = best_type
        
        # Map annotations to cells
        self.adata.obs['cell_type'] = self.adata.obs['leiden'].map(cluster_annotations)
        
        # Plot cell types
        fig, axes = plt.subplots(1, 2, figsize=(14, 5))
        
        sc.pl.umap(self.adata, color='cell_type', ax=axes[0], show=False, legend_loc='on data')
        sc.pl.umap(self.adata, color='cell_type', ax=axes[1], show=False, legend_loc='right margin')
        
        plt.tight_layout()
        plt.savefig(self.output_dir / 'cell_type_annotation.png', dpi=300, bbox_inches='tight')
        plt.close()
        
        # Plot heatmap of marker genes
        logger.info("Generating marker gene heatmap...")
        
        # Select representative markers
        all_markers = []
        for markers in self.marker_dict.values():
            all_markers.extend([g for g in markers if g in self.adata.var_names])
        all_markers = list(set(all_markers))[:50]  # Limit to 50 genes
        
        if len(all_markers) > 0:
            sc.pl.matrixplot(
                self.adata,
                all_markers,
                groupby='cell_type',
                dendrogram=True,
                cmap='RdBu_r',
                show=False
            )
            plt.savefig(self.output_dir / 'marker_heatmap.png', dpi=300, bbox_inches='tight')
            plt.close()
        
        return self.adata
    
    def analyze_dgat1(self) -> sc.AnnData:
        """
        Analyze DGAT1 expression across cell types and conditions.
        
        Returns
        -------
        adata : AnnData
            AnnData with DGAT1 analysis results
        """
        logger.info("Analyzing DGAT1 expression...")
        
        if 'DGAT1' not in self.adata.var_names:
            logger.warning("DGAT1 not found in dataset")
            return self.adata
        
        # Plot DGAT1 expression
        fig, axes = plt.subplots(2, 2, figsize=(14, 12))
        
        sc.pl.umap(
            self.adata,
            color='DGAT1',
            use_raw=False,
            ax=axes[0, 0],
            show=False,
            title='DGAT1 Expression'
        )
        
        sc.pl.violin(
            self.adata,
            keys='DGAT1',
            groupby='cell_type',
            rotation=45,
            use_raw=False,
            ax=axes[0, 1],
            show=False
        )
        
        sc.pl.violin(
            self.adata,
            keys='DGAT1',
            groupby='tissue_type',
            use_raw=False,
            ax=axes[1, 0],
            show=False
        )
        
        if 'malignant_cnv' in self.adata.obs.columns:
            sc.pl.violin(
                self.adata,
                keys='DGAT1',
                groupby='malignant_cnv',
                use_raw=False,
                ax=axes[1, 1],
                show=False
            )
        
        plt.tight_layout()
        plt.savefig(self.output_dir / 'dgat1_expression.png', dpi=300, bbox_inches='tight')
        plt.close()
        
        # Statistical comparison
        logger.info("Computing DGAT1 expression statistics...")
        
        stats_df = self.adata.obs.groupby('cell_type').agg({
            'DGAT1': ['mean', 'std', 'count']
        }).round(3)
        
        stats_df.to_csv(self.output_dir / 'dgat1_statistics_by_celltype.csv')
        
        return self.adata
    
    def analyze_lipid_metabolism(self) -> sc.AnnData:
        """
        Analyze lipid metabolism pathways.
        
        Returns
        -------
        adata : AnnData
            AnnData with lipid metabolism scores
        """
        logger.info("Analyzing lipid metabolism pathways...")
        
        # Score each lipid pathway
        for pathway_name, genes in self.lipid_genes.items():
            available_genes = [g for g in genes if g in self.adata.var_names]
            
            if len(available_genes) > 0:
                logger.info(f"Scoring {pathway_name}: {len(available_genes)} genes available")
                sc.tl.score_genes(
                    self.adata,
                    available_genes,
                    score_name=f'{pathway_name}_score',
                    use_raw=False
                )
        
        # Plot pathway scores
        pathway_scores = [f'{p}_score' for p in self.lipid_genes.keys() 
                          if f'{p}_score' in self.adata.obs.columns]
        
        if len(pathway_scores) > 0:
            # UMAP plots
            n_pathways = len(pathway_scores)
            n_cols = 3
            n_rows = (n_pathways + n_cols - 1) // n_cols
            
            fig, axes = plt.subplots(n_rows, n_cols, figsize=(15, 5*n_rows))
            axes = axes.flatten() if n_rows > 1 else [axes]
            
            for i, score in enumerate(pathway_scores):
                sc.pl.umap(
                    self.adata,
                    color=score,
                    ax=axes[i],
                    show=False,
                    title=score.replace('_score', '')
                )
            
            # Hide unused subplots
            for i in range(len(pathway_scores), len(axes)):
                axes[i].axis('off')
            
            plt.tight_layout()
            plt.savefig(self.output_dir / 'lipid_metabolism_umap.png', dpi=300, bbox_inches='tight')
            plt.close()
            
            # Violin plots by cell type
            fig, axes = plt.subplots(n_rows, n_cols, figsize=(15, 5*n_rows))
            axes = axes.flatten() if n_rows > 1 else [axes]
            
            for i, score in enumerate(pathway_scores):
                sc.pl.violin(
                    self.adata,
                    keys=score,
                    groupby='cell_type',
                    rotation=45,
                    ax=axes[i],
                    show=False
                )
            
            # Hide unused subplots
            for i in range(len(pathway_scores), len(axes)):
                axes[i].axis('off')
            
            plt.tight_layout()
            plt.savefig(self.output_dir / 'lipid_metabolism_violin.png', dpi=300, bbox_inches='tight')
            plt.close()
        
        return self.adata
    
    def differential_expression(self) -> pd.DataFrame:
        """
        Perform differential expression analysis.
        
        Returns
        -------
        de_results : DataFrame
            Differential expression results
        """
        logger.info("Performing differential expression analysis...")
        
        # Tumor vs normal
        if 'tissue_type' in self.adata.obs.columns:
            logger.info("DE: Tumor vs Normal")
            sc.tl.rank_genes_groups(
                self.adata,
                groupby='tissue_type',
                method='wilcoxon',
                use_raw=False,
                key_added='de_tissue'
            )
            
            # Plot top genes
            sc.pl.rank_genes_groups(
                self.adata,
                n_genes=20,
                sharey=False,
                key='de_tissue',
                show=False
            )
            plt.savefig(self.output_dir / 'de_tumor_vs_normal.png', dpi=300, bbox_inches='tight')
            plt.close()
            
            # Export results
            de_results = sc.get.rank_genes_groups_df(self.adata, group=None, key='de_tissue')
            de_results.to_csv(self.output_dir / 'de_tumor_vs_normal.csv', index=False)
        
        # Cell type markers
        logger.info("DE: Cell type markers")
        sc.tl.rank_genes_groups(
            self.adata,
            groupby='cell_type',
            method='wilcoxon',
            use_raw=False,
            key_added='de_celltype'
        )
        
        sc.pl.rank_genes_groups(
            self.adata,
            n_genes=20,
            sharey=False,
            key='de_celltype',
            show=False
        )
        plt.savefig(self.output_dir / 'de_celltype_markers.png', dpi=300, bbox_inches='tight')
        plt.close()
        
        # Export cell type markers
        de_celltype = sc.get.rank_genes_groups_df(self.adata, group=None, key='de_celltype')
        de_celltype.to_csv(self.output_dir / 'de_celltype_markers.csv', index=False)
        
        return de_celltype
    
    def save_results(self) -> None:
        """Save the processed AnnData object and summary statistics."""
        logger.info("Saving results...")
        
        # Save AnnData
        output_path = self.output_dir / 'processed_adata.h5ad'
        self.adata.write(output_path)
        logger.info(f"Saved AnnData to {output_path}")
        
        # Save summary statistics
        summary = {
            'n_cells': self.adata.n_obs,
            'n_genes': self.adata.n_vars,
            'n_clusters': len(self.adata.obs['leiden'].unique()),
            'cell_type_counts': self.adata.obs['cell_type'].value_counts().to_dict(),
            'tissue_type_counts': self.adata.obs['tissue_type'].value_counts().to_dict()
        }
        
        if 'malignant_cnv' in self.adata.obs.columns:
            summary['malignancy_counts'] = self.adata.obs['malignant_cnv'].value_counts().to_dict()
        
        summary_df = pd.DataFrame([summary])
        summary_df.to_csv(self.output_dir / 'analysis_summary.csv', index=False)
        
        logger.info("Results saved successfully")
    
    def run_full_pipeline(
        self,
        min_genes: int = 200,
        max_genes: int = 6000,
        max_mito: float = 15.0,
        n_top_genes: int = 3000,
        n_pcs: int = 50,
        n_neighbors: int = 15,
        resolution: float = 0.5
    ) -> sc.AnnData:
        """
        Run the complete analysis pipeline.
        
        Parameters
        ----------
        min_genes : int
            Minimum genes per cell
        max_genes : int
            Maximum genes per cell
        max_mito : float
            Maximum mitochondrial percentage
        n_top_genes : int
            Number of HVGs
        n_pcs : int
            Number of PCs
        n_neighbors : int
            KNN neighbors
        resolution : float
            Leiden resolution
            
        Returns
        -------
        adata : AnnData
            Fully processed AnnData object
        """
        logger.info("=" * 80)
        logger.info("Starting GBM Single-Cell Analysis Pipeline")
        logger.info("=" * 80)
        
        # 1. Load data
        self.load_data()
        
        # 2. Quality control
        self.quality_control(min_genes, max_genes, max_mito)
        
        # 3. Remove doublets
        self.remove_doublets()
        
        # 4. Normalize and process
        self.normalize_and_process(n_top_genes)
        
        # 5. Batch correction
        self.batch_correction(n_pcs)
        
        # 6. Clustering
        self.clustering(n_neighbors, n_pcs, resolution)
        
        # 7. Infer malignancy
        self.infer_malignancy()
        
        # 8. Annotate cell types
        self.annotate_cell_types()
        
        # 9. DGAT1 analysis
        self.analyze_dgat1()
        
        # 10. Lipid metabolism
        self.analyze_lipid_metabolism()
        
        # 11. Differential expression
        self.differential_expression()
        
        # 12. Save results
        self.save_results()
        
        logger.info("=" * 80)
        logger.info("Pipeline completed successfully!")
        logger.info(f"Results saved to: {self.output_dir}")
        logger.info("=" * 80)
        
        return self.adata


def main():
    """Main entry point for the pipeline."""
    parser = argparse.ArgumentParser(
        description='Single-cell GBM analysis pipeline with DGAT1 focus'
    )
    
    parser.add_argument(
        '--tumor_paths',
        nargs='+',
        required=True,
        help='Paths to tumor sample 10x directories'
    )
    
    parser.add_argument(
        '--normal_paths',
        nargs='+',
        required=True,
        help='Paths to normal sample 10x directories'
    )
    
    parser.add_argument(
        '--output_dir',
        type=str,
        required=True,
        help='Output directory'
    )
    
    parser.add_argument(
        '--min_genes',
        type=int,
        default=200,
        help='Minimum genes per cell'
    )
    
    parser.add_argument(
        '--max_genes',
        type=int,
        default=6000,
        help='Maximum genes per cell'
    )
    
    parser.add_argument(
        '--max_mito',
        type=float,
        default=15.0,
        help='Maximum mitochondrial percentage'
    )
    
    parser.add_argument(
        '--n_top_genes',
        type=int,
        default=3000,
        help='Number of highly variable genes'
    )
    
    parser.add_argument(
        '--n_pcs',
        type=int,
        default=50,
        help='Number of principal components'
    )
    
    parser.add_argument(
        '--n_neighbors',
        type=int,
        default=15,
        help='Number of neighbors for KNN'
    )
    
    parser.add_argument(
        '--resolution',
        type=float,
        default=0.5,
        help='Leiden clustering resolution'
    )
    
    parser.add_argument(
        '--random_seed',
        type=int,
        default=42,
        help='Random seed for reproducibility'
    )
    
    args = parser.parse_args()
    
    # Initialize pipeline
    pipeline = GBMSingleCellPipeline(
        tumor_paths=args.tumor_paths,
        normal_paths=args.normal_paths,
        output_dir=args.output_dir,
        random_seed=args.random_seed
    )
    
    # Run pipeline
    pipeline.run_full_pipeline(
        min_genes=args.min_genes,
        max_genes=args.max_genes,
        max_mito=args.max_mito,
        n_top_genes=args.n_top_genes,
        n_pcs=args.n_pcs,
        n_neighbors=args.n_neighbors,
        resolution=args.resolution
    )


if __name__ == '__main__':
    main()
