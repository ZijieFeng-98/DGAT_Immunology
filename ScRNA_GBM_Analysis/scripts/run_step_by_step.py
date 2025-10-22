"""
Simple non-interactive step-by-step runner
Runs each step automatically with progress reporting
"""

import os
import sys

# Add parent directory to path
sys.path.insert(0, os.path.dirname(__file__))

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

import scanpy as sc
import numpy as np
import pandas as pd

# Configure scanpy
sc.settings.verbosity = 3
sc.settings.set_figure_params(dpi=100, facecolor='white')

def print_banner(step_num, step_name):
    """Print a nice banner for each step"""
    print("\n" + "="*70)
    print(f"   STEP {step_num}: {step_name}")
    print("="*70 + "\n")

def print_stats(adata, step_name):
    """Print statistics after each step"""
    print(f"\n[{step_name}] Statistics:")
    print(f"  - Cells: {adata.n_obs:,}")
    print(f"  - Genes: {adata.n_vars:,}")
    if 'sample' in adata.obs.columns:
        print(f"  - Samples:")
        for sample, count in adata.obs['sample'].value_counts().items():
            print(f"    * {sample}: {count:,} cells")
    print()

def main():
    """Run the complete pipeline step by step with reporting"""
    
    # Configuration
    tumour_path = "data/raw/demo_tumour/"
    normal_path = "data/raw/demo_normal/"
    output_dir = "results/step_by_step/"
    
    # Create output directory
    os.makedirs(output_dir, exist_ok=True)
    sc.settings.figdir = os.path.join(output_dir, 'figures')
    os.makedirs(sc.settings.figdir, exist_ok=True)
    
    print("\n" + "#"*70)
    print("#" + " "*68 + "#")
    print("#" + "  Single-Cell GBM Analysis - Step-by-Step Execution".center(68) + "#")
    print("#" + " "*68 + "#")
    print("#"*70)
    
    # ========== STEP 1: LOAD DATA ==========
    print_banner(1, "LOADING DATA")
    print(f"Tumour data: {tumour_path}")
    print(f"Normal data: {normal_path}")
    
    adata = load_and_concatenate(tumour_path, normal_path)
    print_stats(adata, "After Loading")
    
    # Save checkpoint
    adata.write(os.path.join(output_dir, '01_loaded.h5ad'))
    print("[OK] Checkpoint saved: 01_loaded.h5ad")
    
    # ========== STEP 2: QUALITY CONTROL ==========
    print_banner(2, "QUALITY CONTROL")
    print("Parameters:")
    print("  - Min genes: 200")
    print("  - Max genes: 5000")
    print("  - Max MT%: 10")
    
    n_before = adata.n_obs
    adata = qc_filtering(adata, min_genes=200, max_genes=5000, max_mito=10.0)
    n_after = adata.n_obs
    
    print(f"\nFiltering results:")
    print(f"  - Cells before: {n_before:,}")
    print(f"  - Cells after: {n_after:,}")
    print(f"  - Removed: {n_before - n_after:,} ({(n_before-n_after)/n_before*100:.1f}%)")
    
    print_stats(adata, "After QC")
    
    # QC plots
    sc.pl.violin(adata, ['n_genes', 'n_counts', 'pct_counts_mt'],
                jitter=0.4, multi_panel=True, save='_02_qc_metrics.png')
    
    adata.write(os.path.join(output_dir, '02_qc_filtered.h5ad'))
    print("[OK] Checkpoint saved: 02_qc_filtered.h5ad")
    
    # ========== STEP 3: DOUBLET DETECTION ==========
    print_banner(3, "DOUBLET DETECTION")
    
    n_before = adata.n_obs
    adata = detect_doublets(adata, expected_doublet_rate=0.05)
    n_after = adata.n_obs
    
    if n_after < n_before:
        print(f"\nDoublet removal:")
        print(f"  - Cells before: {n_before:,}")
        print(f"  - Cells after: {n_after:,}")
        print(f"  - Doublets removed: {n_before - n_after:,} ({(n_before-n_after)/n_before*100:.1f}%)")
    else:
        print("Doublet detection skipped (scrublet not available)")
    
    print_stats(adata, "After Doublet Detection")
    adata.write(os.path.join(output_dir, '03_doublets_removed.h5ad'))
    print("[OK] Checkpoint saved: 03_doublets_removed.h5ad")
    
    # ========== STEP 4: NORMALIZATION & INTEGRATION ==========
    print_banner(4, "NORMALIZATION & BATCH CORRECTION")
    print("Processing:")
    print("  - Normalizing to 10,000 counts per cell")
    print("  - Finding highly variable genes (3000)")
    print("  - Scaling data")
    print("  - Computing PCA")
    print("  - Running Harmony batch correction...")
    
    adata = normalize_and_integrate(adata, batch_key='sample')
    
    print(f"\n[OK] Normalization complete")
    print(f"  - HVGs selected: {adata.n_vars:,}")
    print(f"  - PCs computed: {adata.obsm['X_pca'].shape[1]}")
    print(f"  - Harmony PCs: {adata.obsm['X_pca_harmony'].shape[1]}")
    
    print_stats(adata, "After Normalization")
    adata.write(os.path.join(output_dir, '04_normalized.h5ad'))
    print("[OK] Checkpoint saved: 04_normalized.h5ad")
    
    # ========== STEP 5: CLUSTERING ==========
    print_banner(5, "DIMENSIONALITY REDUCTION & CLUSTERING")
    print("Computing:")
    print("  - Neighbor graph")
    print("  - UMAP embedding")
    print("  - Leiden clustering (resolution=0.5)")
    
    adata = dimensionality_reduction_and_clustering(adata, use_rep='X_pca_harmony', resolution=0.5)
    
    n_clusters = len(adata.obs['leiden'].cat.categories)
    print(f"\n[OK] Clustering complete")
    print(f"  - Number of clusters: {n_clusters}")
    print(f"\nCluster sizes:")
    for cluster, count in adata.obs['leiden'].value_counts().sort_index().items():
        print(f"  - Cluster {cluster}: {count:,} cells ({count/adata.n_obs*100:.1f}%)")
    
    # Generate UMAP plots
    sc.pl.umap(adata, color='sample', save='_05_sample.png')
    sc.pl.umap(adata, color='leiden', save='_05_clusters.png')
    
    print_stats(adata, "After Clustering")
    adata.write(os.path.join(output_dir, '05_clustered.h5ad'))
    print("[OK] Checkpoint saved: 05_clustered.h5ad")
    
    # ========== STEP 6: CNV INFERENCE ==========
    print_banner(6, "MALIGNANCY INFERENCE (CNV)")
    print("This step may take longer...")
    print("Inferring copy-number variations...")
    
    try:
        adata = infer_malignancy(adata, reference_category='normal', cnv_key='cnv')
        
        n_malignant = adata.obs['malignant'].sum()
        print(f"\n[OK] CNV inference complete")
        print(f"  - Malignant cells: {n_malignant:,} ({n_malignant/adata.n_obs*100:.1f}%)")
        print(f"  - Normal cells: {adata.n_obs - n_malignant:,} ({(adata.n_obs-n_malignant)/adata.n_obs*100:.1f}%)")
        
        sc.pl.umap(adata, color='malignant', save='_06_malignancy.png')
    except Exception as e:
        print(f"\n[WARNING] CNV inference failed: {e}")
        print("Continuing without malignancy detection...")
        adata.obs['malignant'] = False
    
    print_stats(adata, "After CNV Inference")
    adata.write(os.path.join(output_dir, '06_cnv_inferred.h5ad'))
    print("[OK] Checkpoint saved: 06_cnv_inferred.h5ad")
    
    # ========== STEP 7: CELL-TYPE ANNOTATION ==========
    print_banner(7, "CELL-TYPE ANNOTATION")
    
    marker_dict = {
        'Myeloid': ['PTPRC', 'CSF1R', 'ITGAM', 'CD68'],
        'T cell': ['CD3D', 'CD3E', 'TRAC'],
        'Astrocyte': ['AQP4', 'GFAP'],
        'Oligodendrocyte': ['MOG', 'MBP'],
        'Endothelial': ['PECAM1', 'VWF'],
    }
    
    print("Using marker genes:")
    for cell_type, markers in marker_dict.items():
        print(f"  - {cell_type}: {', '.join(markers)}")
    
    adata = annotate_cell_types(adata, marker_dict)
    
    print(f"\n[OK] Annotation complete")
    print(f"\nCell type distribution:")
    for cell_type, count in adata.obs['cell_type'].value_counts().items():
        print(f"  - {cell_type}: {count:,} cells ({count/adata.n_obs*100:.1f}%)")
    
    sc.pl.umap(adata, color='cell_type', save='_07_cell_types.png')
    
    print_stats(adata, "After Annotation")
    adata.write(os.path.join(output_dir, '07_annotated.h5ad'))
    print("[OK] Checkpoint saved: 07_annotated.h5ad")
    
    # ========== STEP 8: DGAT1 ANALYSIS ==========
    print_banner(8, "DGAT1 & LIPID METABOLISM ANALYSIS")
    
    genes_of_interest = ['DGAT1', 'DGAT2', 'FASN', 'CD68', 'GFAP']
    available_genes = [g for g in genes_of_interest if g in adata.var_names]
    
    print(f"Analyzing genes: {', '.join(available_genes)}")
    
    # Calculate expression summaries
    for gene in available_genes:
        gene_expr = adata[:, gene].X.A1 if hasattr(adata[:, gene].X, 'A1') else adata[:, gene].X.flatten()
        summary = pd.DataFrame({
            'cell_type': adata.obs['cell_type'],
            'sample': adata.obs['sample'],
            'expression': gene_expr
        })
        
        summary_stats = summary.groupby(['cell_type', 'sample'])['expression'].agg(['mean', 'median', 'std', 'count'])
        
        print(f"\n{gene} expression by cell type:")
        print(summary_stats.to_string())
        
        # Save to file
        summary_stats.to_csv(os.path.join(output_dir, f'{gene}_expression_summary.csv'))
    
    # Generate plots
    if available_genes:
        sc.pl.umap(adata, color=available_genes, save='_08_genes.png', 
                  cmap='viridis', ncols=2)
    
    for gene in available_genes[:2]:  # Plot top 2 genes
        plot_results(adata, output_dir, gene)
    
    print_stats(adata, "Final Dataset")
    adata.write(os.path.join(output_dir, '08_final.h5ad'))
    print("[OK] Final dataset saved: 08_final.h5ad")
    
    # ========== SUMMARY ==========
    print("\n" + "#"*70)
    print("#" + " "*68 + "#")
    print("#" + "  ANALYSIS COMPLETE!".center(68) + "#")
    print("#" + " "*68 + "#")
    print("#"*70)
    
    print(f"\n[SUCCESS] All 8 steps completed!")
    print(f"\nResults saved to: {output_dir}")
    print(f"\nOutput files:")
    print(f"  - Final dataset: 08_final.h5ad")
    print(f"  - Checkpoints: 01_loaded.h5ad through 08_final.h5ad")
    print(f"  - Figures: {sc.settings.figdir}/")
    print(f"  - Expression summaries: *_expression_summary.csv")
    
    print(f"\nFinal dataset contains:")
    print(f"  - {adata.n_obs:,} cells")
    print(f"  - {adata.n_vars:,} genes")
    print(f"  - {len(adata.obs['leiden'].cat.categories)} clusters")
    print(f"  - {len(adata.obs['cell_type'].unique())} cell types")
    
    print("\n" + "="*70)
    print("Next steps:")
    print("  1. View plots in:", sc.settings.figdir)
    print("  2. Load data: import scanpy as sc; adata = sc.read_h5ad('results/step_by_step/08_final.h5ad')")
    print("  3. Run downstream analyses: py scripts/downstream_analysis.py")
    print("="*70 + "\n")

if __name__ == "__main__":
    main()

