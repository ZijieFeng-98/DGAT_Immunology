"""
step_by_step.py
===============

Interactive step-by-step runner for the single-cell GBM analysis pipeline.
This script allows you to execute each analysis step individually, inspect
intermediate results, and adjust parameters as needed.

Usage:
    python scripts/step_by_step.py
"""

import os
import sys
import pickle
from typing import Optional, Dict, List

import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt
import seaborn as sns

# Import functions from the main pipeline
from sc_pipeline import (
    load_and_concatenate,
    qc_filtering,
    detect_doublets,
    normalize_and_integrate,
    dimensionality_reduction_and_clustering,
    infer_malignancy,
    annotate_cell_types,
    plot_results
)


class StepByStepPipeline:
    """Interactive step-by-step single-cell analysis pipeline."""
    
    def __init__(self, tumour_path: str, normal_path: str, output_dir: str):
        """Initialize the pipeline with data paths.
        
        Parameters
        ----------
        tumour_path : str
            Path to 10x tumour matrix directory
        normal_path : str
            Path to 10x normal matrix directory
        output_dir : str
            Directory to save results
        """
        self.tumour_path = tumour_path
        self.normal_path = normal_path
        self.output_dir = output_dir
        self.adata = None
        self.current_step = 0
        self.checkpoints = {}
        
        # Create output directories
        os.makedirs(output_dir, exist_ok=True)
        os.makedirs(os.path.join(output_dir, 'checkpoints'), exist_ok=True)
        os.makedirs(os.path.join(output_dir, 'qc_plots'), exist_ok=True)
        
        # Configure scanpy
        sc.settings.verbosity = 3
        sc.settings.figdir = os.path.join(output_dir, 'qc_plots')
        sc.settings.set_figure_params(dpi=80, facecolor='white')
    
    def save_checkpoint(self, step_name: str) -> None:
        """Save a checkpoint of the current AnnData object."""
        checkpoint_path = os.path.join(self.output_dir, 'checkpoints', f'{step_name}.h5ad')
        self.adata.write(checkpoint_path)
        self.checkpoints[step_name] = checkpoint_path
        print(f"✓ Checkpoint saved: {checkpoint_path}")
    
    def load_checkpoint(self, step_name: str) -> None:
        """Load a previously saved checkpoint."""
        if step_name in self.checkpoints:
            self.adata = sc.read_h5ad(self.checkpoints[step_name])
            print(f"✓ Loaded checkpoint: {step_name}")
        else:
            print(f"✗ Checkpoint '{step_name}' not found")
    
    def print_stats(self) -> None:
        """Print current dataset statistics."""
        if self.adata is None:
            print("No data loaded yet")
            return
        
        print("\n" + "="*60)
        print("CURRENT DATASET STATISTICS")
        print("="*60)
        print(f"Number of cells: {self.adata.n_obs:,}")
        print(f"Number of genes: {self.adata.n_vars:,}")
        
        if 'sample' in self.adata.obs.columns:
            print(f"\nCells per sample:")
            print(self.adata.obs['sample'].value_counts())
        
        if 'leiden' in self.adata.obs.columns:
            print(f"\nNumber of clusters: {len(self.adata.obs['leiden'].cat.categories)}")
        
        if 'cell_type' in self.adata.obs.columns:
            print(f"\nCells per type:")
            print(self.adata.obs['cell_type'].value_counts())
        
        print("="*60 + "\n")
    
    # ====== STEP 1: DATA LOADING ======
    
    def step1_load_data(self) -> None:
        """Step 1: Load and concatenate tumour and normal data."""
        print("\n" + "="*60)
        print("STEP 1: LOADING DATA")
        print("="*60)
        print(f"Tumour path: {self.tumour_path}")
        print(f"Normal path: {self.normal_path}")
        
        self.adata = load_and_concatenate(self.tumour_path, self.normal_path)
        self.current_step = 1
        
        print(f"\n✓ Data loaded successfully")
        self.print_stats()
        self.save_checkpoint('01_loaded')
        
        # Basic visualization
        print("Generating initial QC plots...")
        sc.pl.highest_expr_genes(self.adata, n_top=20, save='_01_highest_genes.png')
    
    # ====== STEP 2: QUALITY CONTROL ======
    
    def step2_quality_control(self, min_genes: int = 200, max_genes: int = 2500, 
                              max_mito: float = 10.0, plot: bool = True) -> None:
        """Step 2: Calculate QC metrics and filter cells.
        
        Parameters
        ----------
        min_genes : int
            Minimum genes per cell
        max_genes : int
            Maximum genes per cell
        max_mito : float
            Maximum mitochondrial percentage
        plot : bool
            Whether to generate QC plots
        """
        print("\n" + "="*60)
        print("STEP 2: QUALITY CONTROL")
        print("="*60)
        print(f"Parameters:")
        print(f"  Min genes: {min_genes}")
        print(f"  Max genes: {max_genes}")
        print(f"  Max mito: {max_mito}%")
        
        n_cells_before = self.adata.n_obs
        
        # Apply QC filtering
        self.adata = qc_filtering(self.adata, min_genes, max_genes, max_mito)
        self.current_step = 2
        
        n_cells_after = self.adata.n_obs
        n_removed = n_cells_before - n_cells_after
        pct_removed = (n_removed / n_cells_before) * 100
        
        print(f"\n✓ QC filtering complete")
        print(f"Cells before: {n_cells_before:,}")
        print(f"Cells after: {n_cells_after:,}")
        print(f"Cells removed: {n_removed:,} ({pct_removed:.1f}%)")
        
        self.print_stats()
        self.save_checkpoint('02_qc_filtered')
        
        if plot:
            print("\nGenerating QC plots...")
            # Violin plots of QC metrics
            sc.pl.violin(self.adata, ['n_genes', 'n_counts', 'pct_counts_mt'],
                        jitter=0.4, multi_panel=True, save='_02_qc_violin.png')
            
            # Scatter plots
            sc.pl.scatter(self.adata, x='n_counts', y='pct_counts_mt', save='_02_counts_vs_mito.png')
            sc.pl.scatter(self.adata, x='n_counts', y='n_genes', save='_02_counts_vs_genes.png')
    
    # ====== STEP 3: DOUBLET DETECTION ======
    
    def step3_doublet_detection(self, expected_doublet_rate: float = 0.05, 
                                plot: bool = True) -> None:
        """Step 3: Detect and remove doublets with Scrublet.
        
        Parameters
        ----------
        expected_doublet_rate : float
            Expected fraction of doublets
        plot : bool
            Whether to generate doublet score plots
        """
        print("\n" + "="*60)
        print("STEP 3: DOUBLET DETECTION")
        print("="*60)
        print(f"Expected doublet rate: {expected_doublet_rate:.1%}")
        
        n_cells_before = self.adata.n_obs
        
        # Detect doublets
        self.adata = detect_doublets(self.adata, expected_doublet_rate)
        self.current_step = 3
        
        n_cells_after = self.adata.n_obs
        n_removed = n_cells_before - n_cells_after
        pct_removed = (n_removed / n_cells_before) * 100
        
        print(f"\n✓ Doublet detection complete")
        print(f"Cells before: {n_cells_before:,}")
        print(f"Cells after: {n_cells_after:,}")
        print(f"Doublets removed: {n_removed:,} ({pct_removed:.1f}%)")
        
        self.print_stats()
        self.save_checkpoint('03_doublets_removed')
    
    # ====== STEP 4: NORMALIZATION & INTEGRATION ======
    
    def step4_normalize_integrate(self, batch_key: str = 'sample', 
                                  plot: bool = True) -> None:
        """Step 4: Normalize, find HVGs, and integrate with Harmony.
        
        Parameters
        ----------
        batch_key : str
            Column in .obs for batch correction
        plot : bool
            Whether to generate PCA plots
        """
        print("\n" + "="*60)
        print("STEP 4: NORMALIZATION & INTEGRATION")
        print("="*60)
        print(f"Batch key: {batch_key}")
        
        # Normalize and integrate
        self.adata = normalize_and_integrate(self.adata, batch_key)
        self.current_step = 4
        
        print(f"\n✓ Normalization and Harmony integration complete")
        print(f"Highly variable genes: {self.adata.n_vars:,}")
        print(f"PCs computed: {self.adata.obsm['X_pca'].shape[1]}")
        
        self.save_checkpoint('04_normalized_integrated')
        
        if plot:
            print("\nGenerating integration plots...")
            # PCA before and after Harmony
            sc.pl.pca(self.adata, color='sample', save='_04_pca_before_harmony.png')
            
            # Create temporary UMAP on uncorrected PCs for comparison
            sc.pp.neighbors(self.adata, use_rep='X_pca')
            sc.tl.umap(self.adata)
            sc.pl.umap(self.adata, color='sample', save='_04_umap_before_harmony.png')
    
    # ====== STEP 5: CLUSTERING ======
    
    def step5_clustering(self, resolution: float = 0.5, plot: bool = True) -> None:
        """Step 5: Build neighbor graph, compute UMAP, and cluster.
        
        Parameters
        ----------
        resolution : float
            Leiden clustering resolution
        plot : bool
            Whether to generate UMAP plots
        """
        print("\n" + "="*60)
        print("STEP 5: DIMENSIONALITY REDUCTION & CLUSTERING")
        print("="*60)
        print(f"Resolution: {resolution}")
        
        # Cluster
        self.adata = dimensionality_reduction_and_clustering(
            self.adata, use_rep='X_pca_harmony', resolution=resolution
        )
        self.current_step = 5
        
        print(f"\n✓ Clustering complete")
        print(f"Number of clusters: {len(self.adata.obs['leiden'].cat.categories)}")
        
        self.print_stats()
        self.save_checkpoint('05_clustered')
        
        if plot:
            print("\nGenerating clustering plots...")
            sc.pl.umap(self.adata, color=['sample', 'leiden'], save='_05_umap_clusters.png')
            
            # Cluster composition
            cluster_sample = pd.crosstab(self.adata.obs['leiden'], self.adata.obs['sample'])
            print("\nCluster composition by sample:")
            print(cluster_sample)
    
    # ====== STEP 6: CNV INFERENCE ======
    
    def step6_cnv_inference(self, reference_category: str = 'normal', 
                           plot: bool = True) -> None:
        """Step 6: Infer copy-number variations to identify malignant cells.
        
        Parameters
        ----------
        reference_category : str
            Sample category to use as normal reference
        plot : bool
            Whether to generate CNV plots
        """
        print("\n" + "="*60)
        print("STEP 6: CNV INFERENCE & MALIGNANCY DETECTION")
        print("="*60)
        print(f"Reference: {reference_category}")
        
        # Infer CNV
        self.adata = infer_malignancy(self.adata, reference_category)
        self.current_step = 6
        
        n_malignant = self.adata.obs['malignant'].sum()
        pct_malignant = (n_malignant / self.adata.n_obs) * 100
        
        print(f"\n✓ CNV inference complete")
        print(f"Malignant cells: {n_malignant:,} ({pct_malignant:.1f}%)")
        print(f"Normal cells: {self.adata.n_obs - n_malignant:,} ({100-pct_malignant:.1f}%)")
        
        self.save_checkpoint('06_cnv_inferred')
        
        if plot:
            print("\nGenerating CNV plots...")
            sc.pl.umap(self.adata, color=['malignant', 'cnv_score'], 
                      save='_06_umap_malignancy.png')
            
            # CNV score distribution
            fig, ax = plt.subplots(figsize=(10, 4))
            sns.violinplot(data=self.adata.obs, x='sample', y='cnv_score', ax=ax)
            plt.title('CNV Score Distribution')
            plt.tight_layout()
            plt.savefig(os.path.join(self.output_dir, 'qc_plots', 'cnv_scores.png'))
            plt.close()
    
    # ====== STEP 7: CELL-TYPE ANNOTATION ======
    
    def step7_annotate_cell_types(self, custom_markers: Optional[Dict[str, List[str]]] = None,
                                  plot: bool = True) -> None:
        """Step 7: Annotate cell types based on marker genes.
        
        Parameters
        ----------
        custom_markers : dict, optional
            Custom marker dictionary
        plot : bool
            Whether to generate annotation plots
        """
        print("\n" + "="*60)
        print("STEP 7: CELL-TYPE ANNOTATION")
        print("="*60)
        
        # Default marker dict
        marker_dict = {
            'Myeloid': ['PTPRC', 'CSF1R', 'ITGAM', 'CD68', 'CD163'],
            'Microglia': ['P2RY12', 'TMEM119', 'CX3CR1'],
            'T cell': ['CD3D', 'CD3E', 'TRAC', 'CD8A'],
            'Astrocyte': ['AQP4', 'GFAP', 'SLC1A2'],
            'Oligodendrocyte': ['MOG', 'MBP', 'PLP1'],
            'OPC': ['PDGFRA', 'CSPG4', 'SOX10'],
            'Endothelial': ['PECAM1', 'VWF', 'CDH5'],
            'Neuron': ['RBFOX3', 'SYT1', 'SNAP25'],
        }
        
        if custom_markers:
            marker_dict.update(custom_markers)
        
        print("Using markers:")
        for cell_type, markers in marker_dict.items():
            print(f"  {cell_type}: {', '.join(markers[:3])}...")
        
        # Annotate
        self.adata = annotate_cell_types(self.adata, marker_dict)
        self.current_step = 7
        
        print(f"\n✓ Cell-type annotation complete")
        self.print_stats()
        self.save_checkpoint('07_annotated')
        
        if plot:
            print("\nGenerating annotation plots...")
            sc.pl.umap(self.adata, color=['cell_type', 'leiden'], 
                      save='_07_umap_cell_types.png')
            
            # Marker gene dot plot
            marker_genes = []
            for markers in marker_dict.values():
                marker_genes.extend([m for m in markers if m in self.adata.var_names])
            
            if marker_genes:
                sc.pl.dotplot(self.adata, marker_genes[:30], groupby='cell_type',
                             save='_07_markers_dotplot.png')
    
    # ====== STEP 8: DGAT1 ANALYSIS ======
    
    def step8_dgat1_analysis(self, genes_of_interest: Optional[List[str]] = None,
                            plot: bool = True) -> None:
        """Step 8: Analyze DGAT1 and lipid metabolism gene expression.
        
        Parameters
        ----------
        genes_of_interest : list, optional
            List of genes to analyze (default: DGAT1)
        plot : bool
            Whether to generate expression plots
        """
        print("\n" + "="*60)
        print("STEP 8: DGAT1 & LIPID METABOLISM ANALYSIS")
        print("="*60)
        
        if genes_of_interest is None:
            genes_of_interest = ['DGAT1', 'DGAT2', 'FASN', 'ACLY', 'FABP4', 'FABP5']
        
        available_genes = [g for g in genes_of_interest if g in self.adata.var_names]
        
        print(f"Analyzing genes: {', '.join(available_genes)}")
        
        # Calculate expression summaries
        for gene in available_genes:
            gene_expr = self.adata[:, gene].X.A1 if hasattr(self.adata[:, gene].X, 'A1') else self.adata[:, gene].X.flatten()
            summary = pd.DataFrame({
                'cell_type': self.adata.obs['cell_type'],
                'sample': self.adata.obs['sample'],
                'expression': gene_expr
            })
            
            summary_stats = summary.groupby(['cell_type', 'sample'])['expression'].agg(['mean', 'median', 'std', 'count'])
            print(f"\n{gene} expression by cell type and sample:")
            print(summary_stats)
            
            # Save to file
            summary_stats.to_csv(os.path.join(self.output_dir, f'{gene}_expression_summary.csv'))
        
        self.current_step = 8
        self.save_checkpoint('08_final')
        
        if plot and available_genes:
            print("\nGenerating expression plots...")
            # UMAP with gene expression
            sc.pl.umap(self.adata, color=available_genes, save='_08_umap_genes.png', 
                      cmap='viridis', vmax='p99')
            
            # Violin plots
            for gene in available_genes:
                plot_results(self.adata, self.output_dir, gene)
            
            # Dot plot
            sc.pl.dotplot(self.adata, available_genes, groupby='cell_type',
                         save='_08_lipid_genes_dotplot.png')
    
    # ====== CONVENIENCE METHODS ======
    
    def run_all_steps(self, **kwargs) -> None:
        """Run all steps sequentially with default or provided parameters."""
        print("\n" + "="*60)
        print("RUNNING COMPLETE PIPELINE")
        print("="*60)
        
        self.step1_load_data()
        self.step2_quality_control(**kwargs.get('qc_params', {}))
        self.step3_doublet_detection(**kwargs.get('doublet_params', {}))
        self.step4_normalize_integrate(**kwargs.get('normalize_params', {}))
        self.step5_clustering(**kwargs.get('cluster_params', {}))
        self.step6_cnv_inference(**kwargs.get('cnv_params', {}))
        self.step7_annotate_cell_types(**kwargs.get('annotation_params', {}))
        self.step8_dgat1_analysis(**kwargs.get('dgat1_params', {}))
        
        print("\n" + "="*60)
        print("✓ PIPELINE COMPLETE!")
        print("="*60)
        self.print_stats()
    
    def interactive_menu(self) -> None:
        """Display an interactive menu for step-by-step execution."""
        while True:
            print("\n" + "="*60)
            print("SINGLE-CELL GBM ANALYSIS - INTERACTIVE MENU")
            print("="*60)
            print(f"Current step: {self.current_step}/8")
            print("\nAvailable steps:")
            print("  1. Load data")
            print("  2. Quality control")
            print("  3. Doublet detection")
            print("  4. Normalization & integration")
            print("  5. Clustering")
            print("  6. CNV inference")
            print("  7. Cell-type annotation")
            print("  8. DGAT1 analysis")
            print("\nOther options:")
            print("  9. Run all remaining steps")
            print("  s. Show current statistics")
            print("  l. Load checkpoint")
            print("  q. Quit")
            print("="*60)
            
            choice = input("\nEnter your choice: ").strip().lower()
            
            if choice == '1':
                self.step1_load_data()
            elif choice == '2':
                self.step2_quality_control()
            elif choice == '3':
                self.step3_doublet_detection()
            elif choice == '4':
                self.step4_normalize_integrate()
            elif choice == '5':
                self.step5_clustering()
            elif choice == '6':
                self.step6_cnv_inference()
            elif choice == '7':
                self.step7_annotate_cell_types()
            elif choice == '8':
                self.step8_dgat1_analysis()
            elif choice == '9':
                self.run_all_steps()
            elif choice == 's':
                self.print_stats()
            elif choice == 'l':
                print("\nAvailable checkpoints:")
                for i, name in enumerate(self.checkpoints.keys(), 1):
                    print(f"  {i}. {name}")
                checkpoint_choice = input("Enter checkpoint name: ").strip()
                self.load_checkpoint(checkpoint_choice)
            elif choice == 'q':
                print("Exiting...")
                break
            else:
                print("Invalid choice. Please try again.")


def main():
    """Main entry point for interactive step-by-step analysis."""
    import argparse
    
    parser = argparse.ArgumentParser(description='Interactive step-by-step scRNA-seq analysis')
    parser.add_argument('--tumour_path', type=str, default='data/raw/tumour/',
                       help='Path to 10x tumour matrix')
    parser.add_argument('--normal_path', type=str, default='data/raw/normal/',
                       help='Path to 10x normal matrix')
    parser.add_argument('--output_dir', type=str, default='results/',
                       help='Output directory')
    parser.add_argument('--auto', action='store_true',
                       help='Run all steps automatically without prompts')
    
    args = parser.parse_args()
    
    # Initialize pipeline
    pipeline = StepByStepPipeline(
        tumour_path=args.tumour_path,
        normal_path=args.normal_path,
        output_dir=args.output_dir
    )
    
    if args.auto:
        # Run all steps automatically
        pipeline.run_all_steps()
    else:
        # Interactive menu
        pipeline.interactive_menu()


if __name__ == '__main__':
    main()

