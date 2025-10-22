"""
Analyze ALL 21 Samples from GSE222520 (Gupta et al. 2024)

Complete analysis of:
- 3 Normal brain (NGB1-3)
- 4 IDH-mutant primary (IMP1-4)
- 6 IDH-mutant recurrent (IMR1-6)
- 4 IDH-wildtype primary (IWP1-4)
- 4 IDH-wildtype recurrent (IWR1-4)

Total: 21 samples, ~144,000 cells
"""

import os
import numpy as np
import pandas as pd
import scanpy as sc
import harmonypy
from pathlib import Path
import warnings
warnings.filterwarnings('ignore')

def load_all_samples(base_dir='data/raw/gse222520'):
    """
    Load all 21 samples from GSE222520
    """
    print("\n" + "="*70)
    print("LOADING ALL 21 SAMPLES FROM GSE222520")
    print("="*70)
    
    # Define sample groups
    samples = {
        'NGB': ['GSM6925378_NGB1', 'GSM6925379_NGB2', 'GSM6925380_NGB4'],
        'IMP': ['GSM6925381_IMP1', 'GSM6925382_IMP2', 'GSM6925383_IMP3', 'GSM6925384_IMP4'],
        'IMR': ['GSM6925385_IMR1', 'GSM6925386_IMR2', 'GSM6925387_IMR3', 
                'GSM6925388_IMR4', 'GSM6925389_IMR5', 'GSM6925390_IMR6'],
        'IWP': ['GSM6925391_IWP1', 'GSM6925392_IWP2', 'GSM6925393_IWP3', 'GSM6925394_IWP4'],
        'IWR': ['GSM6925395_IWR1', 'GSM6925396_IWR2', 'GSM6925397_IWR3', 'GSM6925398_IWR4']
    }
    
    adatas = []
    sample_info = []
    
    for group, sample_list in samples.items():
        for sample_id in sample_list:
            # Extract sample name (e.g., NGB1, IMP1)
            sample_name = sample_id.split('_')[1]
            
            # Construct path
            if group == 'NGB' and sample_name == 'NGB4':
                # Special case: NGB4 folder is named NGB3
                sample_path = os.path.join(base_dir, sample_id, 'NGB3', 'filtered_feature_bc_matrix')
            else:
                sample_path = os.path.join(base_dir, sample_id, sample_name, 'filtered_feature_bc_matrix')
            
            if not os.path.exists(sample_path):
                print(f"  [WARNING] Skipping {sample_name} - path not found: {sample_path}")
                continue
            
            try:
                # Load 10x data
                adata = sc.read_10x_mtx(sample_path, var_names='gene_symbols', cache=True)
                
                # Add metadata
                adata.obs['sample'] = sample_name
                adata.obs['group'] = group
                adata.obs['sample_id'] = sample_id
                
                # Add detailed group info
                if group == 'NGB':
                    adata.obs['tissue_type'] = 'Normal'
                    adata.obs['idh_status'] = 'N/A'
                    adata.obs['tumor_status'] = 'Normal'
                elif group == 'IMP':
                    adata.obs['tissue_type'] = 'Tumor'
                    adata.obs['idh_status'] = 'IDH-mutant'
                    adata.obs['tumor_status'] = 'Primary'
                elif group == 'IMR':
                    adata.obs['tissue_type'] = 'Tumor'
                    adata.obs['idh_status'] = 'IDH-mutant'
                    adata.obs['tumor_status'] = 'Recurrent'
                elif group == 'IWP':
                    adata.obs['tissue_type'] = 'Tumor'
                    adata.obs['idh_status'] = 'IDH-wildtype'
                    adata.obs['tumor_status'] = 'Primary'
                elif group == 'IWR':
                    adata.obs['tissue_type'] = 'Tumor'
                    adata.obs['idh_status'] = 'IDH-wildtype'
                    adata.obs['tumor_status'] = 'Recurrent'
                
                adatas.append(adata)
                sample_info.append({
                    'sample': sample_name,
                    'group': group,
                    'n_cells': adata.n_obs,
                    'n_genes': adata.n_vars
                })
                
                print(f"  [OK] Loaded {sample_name} ({group}): {adata.n_obs:,} cells, {adata.n_vars:,} genes")
                
            except Exception as e:
                print(f"  [ERROR] Failed to load {sample_name}: {e}")
    
    if len(adatas) == 0:
        raise ValueError("No samples loaded successfully!")
    
    # Combine all samples
    print(f"\n[COMBINING] Merging {len(adatas)} samples...")
    adata_combined = adatas[0].concatenate(adatas[1:], batch_key='sample', 
                                          batch_categories=[a.obs['sample'][0] for a in adatas],
                                          index_unique='-')
    
    # Create summary table
    summary_df = pd.DataFrame(sample_info)
    
    print(f"\n[SUCCESS] Combined dataset:")
    print(f"  - Total cells: {adata_combined.n_obs:,}")
    print(f"  - Total genes: {adata_combined.n_vars:,}")
    print(f"  - Samples: {len(adatas)}")
    
    print(f"\n[SUMMARY] By group:")
    group_summary = adata_combined.obs.groupby('group').size()
    for group, count in group_summary.items():
        pct = count / adata_combined.n_obs * 100
        print(f"  - {group}: {count:,} cells ({pct:.1f}%)")
    
    return adata_combined, summary_df

def quality_control(adata, min_genes=200, max_genes=7000, max_mito=15):
    """
    Quality control filtering
    """
    print("\n" + "="*70)
    print("QUALITY CONTROL")
    print("="*70)
    
    # Calculate QC metrics
    adata.var['mt'] = adata.var_names.str.startswith('MT-')
    sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)
    
    print(f"\n[BEFORE] QC:")
    print(f"  - Cells: {adata.n_obs:,}")
    print(f"  - Genes: {adata.n_vars:,}")
    
    # Filter cells
    adata = adata[
        (adata.obs['n_genes_by_counts'] >= min_genes) &
        (adata.obs['n_genes_by_counts'] <= max_genes) &
        (adata.obs['pct_counts_mt'] < max_mito)
    ].copy()
    
    print(f"\n[AFTER] QC:")
    print(f"  - Cells: {adata.n_obs:,}")
    print(f"  - Genes: {adata.n_vars:,}")
    print(f"  - Filtered: {adata.n_obs} cells passed QC")
    
    return adata

def normalize_and_process(adata, n_top_genes=3000):
    """
    Normalization, HVG selection, scaling, PCA
    """
    print("\n" + "="*70)
    print("NORMALIZATION & PROCESSING")
    print("="*70)
    
    # Normalize
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    
    # Find highly variable genes
    sc.pp.highly_variable_genes(adata, flavor='seurat', batch_key='sample', 
                                n_top_genes=n_top_genes)
    
    print(f"  [OK] Found {adata.var['highly_variable'].sum()} highly variable genes")
    
    # Subset to HVGs and scale
    adata = adata[:, adata.var['highly_variable']].copy()
    sc.pp.scale(adata, max_value=10)
    
    # PCA
    sc.tl.pca(adata, svd_solver='arpack')
    print(f"  [OK] PCA complete ({adata.obsm['X_pca'].shape[1]} components)")
    
    return adata

def batch_correction(adata, batch_key='sample'):
    """
    Harmony batch correction
    """
    print("\n" + "="*70)
    print("BATCH CORRECTION (Harmony)")
    print("="*70)
    
    # Run Harmony
    harmony_out = harmonypy.run_harmony(
        adata.obsm['X_pca'],
        adata.obs,
        batch_key,
        max_iter_harmony=20
    )
    
    adata.obsm['X_pca_harmony'] = harmony_out.Z_corr.T
    
    # Get number of iterations (objective_harmony is a list)
    n_iter = len(harmony_out.objective_harmony) if hasattr(harmony_out, 'objective_harmony') else 'unknown'
    print(f"  [OK] Harmony converged in {n_iter} iterations")
    
    return adata

def clustering_and_umap(adata, resolution=0.8):
    """
    Clustering and UMAP
    """
    print("\n" + "="*70)
    print("CLUSTERING & UMAP")
    print("="*70)
    
    # Neighbors and UMAP
    sc.pp.neighbors(adata, use_rep='X_pca_harmony', n_pcs=50)
    sc.tl.umap(adata)
    
    # Clustering
    sc.tl.leiden(adata, resolution=resolution)
    
    n_clusters = adata.obs['leiden'].nunique()
    print(f"  [OK] Found {n_clusters} clusters (resolution={resolution})")
    
    print(f"\n[CLUSTERS] Size distribution:")
    for cluster, count in adata.obs['leiden'].value_counts().sort_index().items():
        pct = count / adata.n_obs * 100
        print(f"  - Cluster {cluster}: {count:,} cells ({pct:.1f}%)")
    
    return adata

def annotate_cell_types(adata):
    """
    Annotate cell types using marker genes
    """
    print("\n" + "="*70)
    print("CELL TYPE ANNOTATION")
    print("="*70)
    
    # Extended marker dictionary
    markers = {
        # Myeloid
        'Microglia': ['CX3CR1', 'P2RY12', 'TMEM119', 'GPR34'],
        'Macrophage': ['CD68', 'CD14', 'LYZ', 'C1QA'],
        'DC': ['CLEC9A', 'CD1C', 'FCER1A'],
        'MDSC': ['S100A8', 'S100A9', 'LILRB2'],
        # Lymphoid
        'CD8_T': ['CD3D', 'CD8A', 'GZMB', 'PRF1'],
        'CD4_T': ['CD3D', 'CD4', 'IL7R'],
        'NK': ['NKG7', 'KLRD1', 'NCAM1'],
        'B_cell': ['MS4A1', 'CD79A', 'CD79B'],
        'Plasma': ['MZB1', 'IGHG1', 'IGHA1'],
        # Other
        'Astrocyte': ['GFAP', 'AQP4', 'SLC1A3'],
        'Oligodendrocyte': ['MOG', 'MBP', 'PLP1'],
        'Endothelial': ['PECAM1', 'VWF', 'CDH5'],
        'Cycling': ['MKI67', 'TOP2A', 'PCNA']
    }
    
    # Annotate
    for cluster in adata.obs['leiden'].cat.categories:
        idx = adata.obs['leiden'] == cluster
        scores = {}
        
        for cell_type, marker_genes in markers.items():
            valid = [g for g in marker_genes if g in adata.var_names]
            if len(valid) == 0:
                scores[cell_type] = -np.inf
            else:
                scores[cell_type] = np.array(adata[idx, valid].X.mean())
        
        assigned = max(scores, key=scores.get)
        adata.obs.loc[idx, 'cell_type'] = assigned
    
    adata.obs['cell_type'] = adata.obs['cell_type'].astype('category')
    
    print(f"\n[CELL TYPES] Identified:")
    for ct, count in adata.obs['cell_type'].value_counts().items():
        pct = count / adata.n_obs * 100
        print(f"  - {ct}: {count:,} cells ({pct:.1f}%)")
    
    return adata

def main():
    """Main pipeline"""
    import argparse
    
    parser = argparse.ArgumentParser(description='Analyze all 21 samples from GSE222520')
    parser.add_argument('--data_dir', type=str, default='data/raw/gse222520',
                       help='Base directory with sample data')
    parser.add_argument('--output', type=str, default='results/all_samples',
                       help='Output directory')
    parser.add_argument('--min_genes', type=int, default=200,
                       help='Minimum genes per cell')
    parser.add_argument('--max_genes', type=int, default=7000,
                       help='Maximum genes per cell')
    parser.add_argument('--max_mito', type=float, default=15,
                       help='Maximum mitochondrial %')
    parser.add_argument('--resolution', type=float, default=0.8,
                       help='Clustering resolution')
    parser.add_argument('--n_top_genes', type=int, default=3000,
                       help='Number of highly variable genes')
    
    args = parser.parse_args()
    
    print("="*70)
    print("COMPREHENSIVE ANALYSIS: ALL 21 SAMPLES")
    print("Gupta et al. 2024 (GSE222520)")
    print("="*70)
    print(f"\nOutput: {args.output}")
    
    os.makedirs(args.output, exist_ok=True)
    
    # Load all samples
    adata, summary = load_all_samples(args.data_dir)
    
    # Save sample summary
    summary.to_csv(os.path.join(args.output, 'sample_summary.csv'), index=False)
    
    # QC
    adata = quality_control(adata, args.min_genes, args.max_genes, args.max_mito)
    
    # Save raw counts
    adata.raw = adata
    
    # Process
    adata = normalize_and_process(adata, args.n_top_genes)
    
    # Batch correction
    adata = batch_correction(adata, batch_key='sample')
    
    # Clustering
    adata = clustering_and_umap(adata, args.resolution)
    
    # Annotate
    adata = annotate_cell_types(adata)
    
    # Save results
    output_file = os.path.join(args.output, 'all_samples_processed.h5ad')
    adata.write(output_file)
    
    print("\n" + "="*70)
    print("[SUCCESS] ANALYSIS COMPLETE!")
    print("="*70)
    print(f"\nResults saved to: {args.output}/")
    print(f"  - all_samples_processed.h5ad ({adata.n_obs:,} cells)")
    print(f"  - sample_summary.csv")
    print("\nNext steps:")
    print("  1. Run Gupta-style figures on all samples")
    print("  2. Analyze DGAT1 expression across all groups")
    print("  3. Compare primary vs recurrent tumors")
    print("  4. IDH-mutant vs IDH-wildtype analysis")

if __name__ == "__main__":
    main()

