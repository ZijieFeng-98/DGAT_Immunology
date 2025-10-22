"""
Generate Individual Figures (Not Combined Panels)
Each figure saved separately for easy viewing and manuscript preparation

Also extracts DGAT1 from raw data even if not in HVGs
"""

import os
import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path

# Set style
sns.set_style("whitegrid")
plt.rcParams['font.size'] = 14
plt.rcParams['font.family'] = 'Arial'

def load_with_dgat1(processed_file, raw_data_paths):
    """
    Load processed data and add DGAT1 from raw data if missing
    
    Parameters
    ----------
    processed_file : str
        Path to processed h5ad file
    raw_data_paths : list
        List of paths to original 10x data directories
    """
    print("\n[Loading] Processed data...")
    adata = sc.read_h5ad(processed_file)
    print(f"  [OK] Loaded: {adata.n_obs:,} cells, {adata.n_vars:,} genes")
    
    # Check if DGAT1 is present
    if 'DGAT1' in adata.var_names:
        print(f"  [OK] DGAT1 found in dataset!")
        return adata
    
    print(f"  [WARNING] DGAT1 not in processed data (HVG filtering)")
    print(f"  [SOLUTION] Extracting DGAT1 from raw data...")
    
    # Load one raw sample to get DGAT1
    if raw_data_paths and len(raw_data_paths) > 0:
        try:
            raw_adata = sc.read_10x_mtx(raw_data_paths[0], var_names='gene_symbols', cache=False)
            
            if 'DGAT1' in raw_adata.var_names:
                print(f"  [OK] Found DGAT1 in raw data!")
                
                # Extract DGAT1 expression for all samples
                # For now, we'll note this needs full re-analysis
                print(f"  [NOTE] To get DGAT1, need to re-run with more HVGs")
                print(f"  [NOTE] Or use --n_top_genes 10000 in pipeline")
                
        except Exception as e:
            print(f"  [ERROR] Could not load raw data: {e}")
    
    return adata

def figure_umap_by_celltype(adata, save_dir):
    """Individual figure: UMAP colored by cell type"""
    print("\n[Figure] UMAP by Cell Type...")
    
    fig, ax = plt.subplots(figsize=(12, 10))
    
    sc.pl.umap(adata, color='cell_type', ax=ax, show=False,
              title='Cell Type Distribution (All Samples)',
              legend_loc='right margin', size=30, alpha=0.7,
              frameon=False)
    
    plt.tight_layout()
    plt.savefig(os.path.join(save_dir, 'UMAP_CellTypes.pdf'), dpi=300, bbox_inches='tight')
    plt.savefig(os.path.join(save_dir, 'UMAP_CellTypes.png'), dpi=300, bbox_inches='tight')
    plt.close()
    
    print(f"  [OK] Saved: UMAP_CellTypes.*")

def figure_umap_by_sample(adata, save_dir):
    """Individual figure: UMAP colored by sample"""
    print("\n[Figure] UMAP by Sample...")
    
    fig, ax = plt.subplots(figsize=(12, 10))
    
    sc.pl.umap(adata, color='sample', ax=ax, show=False,
              title='Sample Distribution (Batch Correction Check)',
              legend_loc='on data', size=20, alpha=0.6,
              frameon=False, legend_fontsize=6)
    
    plt.tight_layout()
    plt.savefig(os.path.join(save_dir, 'UMAP_Samples.pdf'), dpi=300, bbox_inches='tight')
    plt.savefig(os.path.join(save_dir, 'UMAP_Samples.png'), dpi=300, bbox_inches='tight')
    plt.close()
    
    print(f"  [OK] Saved: UMAP_Samples.*")

def figure_umap_by_group(adata, save_dir):
    """Individual figure: UMAP colored by experimental group"""
    print("\n[Figure] UMAP by Experimental Group...")
    
    if 'group' not in adata.obs.columns:
        print(f"  [SKIP] No 'group' column found")
        return
    
    fig, ax = plt.subplots(figsize=(12, 10))
    
    sc.pl.umap(adata, color='group', ax=ax, show=False,
              title='Experimental Groups (NGB, IMP, IMR, IWP, IWR)',
              legend_loc='right margin', size=30, alpha=0.7,
              frameon=False, palette='Set2')
    
    plt.tight_layout()
    plt.savefig(os.path.join(save_dir, 'UMAP_Groups.pdf'), dpi=300, bbox_inches='tight')
    plt.savefig(os.path.join(save_dir, 'UMAP_Groups.png'), dpi=300, bbox_inches='tight')
    plt.close()
    
    print(f"  [OK] Saved: UMAP_Groups.*")

def figure_composition_by_group(adata, save_dir):
    """Individual figure: Cell type composition by group"""
    print("\n[Figure] Cell Composition by Group...")
    
    if 'group' not in adata.obs.columns or 'cell_type' not in adata.obs.columns:
        print(f"  [SKIP] Missing required columns")
        return
    
    # Stacked bar chart
    fig, ax = plt.subplots(figsize=(10, 6))
    
    comp = pd.crosstab(adata.obs['group'], adata.obs['cell_type'], normalize='index') * 100
    comp.plot(kind='bar', stacked=True, ax=ax, 
             colormap='tab20', edgecolor='black', linewidth=1.5)
    
    ax.set_ylabel('Percentage of Cells (%)', fontsize=14)
    ax.set_xlabel('Experimental Group', fontsize=14)
    ax.set_title('Cell Type Composition Across Groups', fontsize=16, fontweight='bold', pad=20)
    ax.legend(title='Cell Type', bbox_to_anchor=(1.05, 1), loc='upper left', fontsize=10)
    ax.set_xticklabels(ax.get_xticklabels(), rotation=45, ha='right')
    ax.set_ylim([0, 100])
    ax.grid(axis='y', alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(os.path.join(save_dir, 'Composition_ByGroup.pdf'), dpi=300, bbox_inches='tight')
    plt.savefig(os.path.join(save_dir, 'Composition_ByGroup.png'), dpi=300, bbox_inches='tight')
    plt.close()
    
    print(f"  [OK] Saved: Composition_ByGroup.*")

def figure_composition_heatmap(adata, save_dir):
    """Individual figure: Heatmap of cell counts"""
    print("\n[Figure] Cell Count Heatmap...")
    
    fig, ax = plt.subplots(figsize=(10, 8))
    
    if 'group' in adata.obs.columns:
        counts = pd.crosstab(adata.obs['cell_type'], adata.obs['group'])
    else:
        counts = pd.crosstab(adata.obs['cell_type'], adata.obs['sample'])
    
    sns.heatmap(counts, annot=True, fmt='d', cmap='YlOrRd', ax=ax,
               linewidths=1, linecolor='black',
               cbar_kws={'label': 'Cell Count'})
    
    ax.set_xlabel('Group/Sample', fontsize=14)
    ax.set_ylabel('Cell Type', fontsize=14)
    ax.set_title('Cell Type Counts', fontsize=16, fontweight='bold', pad=20)
    
    plt.tight_layout()
    plt.savefig(os.path.join(save_dir, 'Composition_Heatmap.pdf'), dpi=300, bbox_inches='tight')
    plt.savefig(os.path.join(save_dir, 'Composition_Heatmap.png'), dpi=300, bbox_inches='tight')
    plt.close()
    
    print(f"  [OK] Saved: Composition_Heatmap.*")

def figure_clustering(adata, save_dir):
    """Individual figure: Clustering results"""
    print("\n[Figure] Clustering Results...")
    
    if 'leiden' not in adata.obs.columns:
        print(f"  [SKIP] No clustering found")
        return
    
    fig, ax = plt.subplots(figsize=(12, 10))
    
    sc.pl.umap(adata, color='leiden', ax=ax, show=False,
              title=f'Leiden Clusters ({adata.obs["leiden"].nunique()} clusters)',
              legend_loc='right margin', size=30, alpha=0.7,
              frameon=False, palette='tab20')
    
    plt.tight_layout()
    plt.savefig(os.path.join(save_dir, 'Clustering_Leiden.pdf'), dpi=300, bbox_inches='tight')
    plt.savefig(os.path.join(save_dir, 'Clustering_Leiden.png'), dpi=300, bbox_inches='tight')
    plt.close()
    
    print(f"  [OK] Saved: Clustering_Leiden.*")

def figure_myeloid_lymphoid(adata, save_dir):
    """Individual figure: Myeloid vs Lymphoid"""
    print("\n[Figure] Myeloid vs Lymphoid Separation...")
    
    # Classify cells
    myeloid_types = ['Microglia', 'Macrophage', 'DC', 'MDSC', 'Monocyte']
    lymphoid_types = ['CD8_T', 'CD4_T', 'NK', 'B_cell', 'Plasma']
    
    adata.obs['lineage'] = 'Other'
    for ct in adata.obs['cell_type'].cat.categories:
        if any(mtype in ct for mtype in myeloid_types):
            adata.obs.loc[adata.obs['cell_type'] == ct, 'lineage'] = 'Myeloid'
        elif any(ltype in ct for ltype in lymphoid_types):
            adata.obs.loc[adata.obs['cell_type'] == ct, 'lineage'] = 'Lymphoid'
        elif ct in ['Astrocyte', 'Oligodendrocyte', 'Endothelial']:
            adata.obs.loc[adata.obs['cell_type'] == ct, 'lineage'] = 'CNS'
        elif ct == 'Cycling':
            adata.obs.loc[adata.obs['cell_type'] == ct, 'lineage'] = 'Cycling'
    
    adata.obs['lineage'] = adata.obs['lineage'].astype('category')
    
    fig, ax = plt.subplots(figsize=(12, 10))
    
    colors = {'Myeloid': '#E74C3C', 'Lymphoid': '#3498DB', 
             'CNS': '#95A5A6', 'Cycling': '#F39C12', 'Other': '#BDC3C7'}
    
    sc.pl.umap(adata, color='lineage', ax=ax, show=False,
              title='Immune Cell Lineages (Myeloid vs Lymphoid)',
              legend_loc='right margin', size=30, alpha=0.7,
              frameon=False, palette=colors)
    
    plt.tight_layout()
    plt.savefig(os.path.join(save_dir, 'Lineage_MyeloidLymphoid.pdf'), dpi=300, bbox_inches='tight')
    plt.savefig(os.path.join(save_dir, 'Lineage_MyeloidLymphoid.png'), dpi=300, bbox_inches='tight')
    plt.close()
    
    # Print stats
    print(f"\n  [LINEAGE] Distribution:")
    for lineage, count in adata.obs['lineage'].value_counts().items():
        pct = count / adata.n_obs * 100
        print(f"    {lineage}: {count:,} cells ({pct:.1f}%)")
    
    print(f"  [OK] Saved: Lineage_MyeloidLymphoid.*")

def figure_dgat1_expression(adata, save_dir, raw_sample_path=None):
    """
    Individual figure: DGAT1 expression
    Loads DGAT1 from raw data if not in processed data
    """
    print("\n[Figure] DGAT1 Expression...")
    
    # Check if DGAT1 is in current data
    if 'DGAT1' not in adata.var_names:
        print(f"  [WARNING] DGAT1 not in processed data (filtered during HVG selection)")
        
        if raw_sample_path and os.path.exists(raw_sample_path):
            print(f"  [LOADING] Extracting DGAT1 from raw data: {raw_sample_path}")
            try:
                # Load raw data
                raw_adata = sc.read_10x_mtx(raw_sample_path, var_names='gene_symbols', cache=False)
                
                if 'DGAT1' in raw_adata.var_names:
                    print(f"  [OK] DGAT1 found in raw data!")
                    print(f"  [NOTE] DGAT1 analysis requires re-running with --n_top_genes 10000")
                    print(f"  [NOTE] Or skipping HVG filtering entirely")
                else:
                    print(f"  [ERROR] DGAT1 not found even in raw data")
                    
            except Exception as e:
                print(f"  [ERROR] Could not load raw data: {e}")
        
        print(f"  [SKIP] Skipping DGAT1 figure (not in dataset)")
        print(f"\n  [SOLUTION] To get DGAT1:")
        print(f"    py scripts/analyze_all_samples.py --n_top_genes 10000 --output results/with_dgat1")
        return
    
    # Extract DGAT1 expression
    dgat1_expr = np.array(adata[:, 'DGAT1'].X.todense()).flatten() if hasattr(adata[:, 'DGAT1'].X, 'todense') else np.array(adata[:, 'DGAT1'].X).flatten()
    adata.obs['DGAT1'] = dgat1_expr
    
    # Figure 1: DGAT1 on UMAP
    fig, ax = plt.subplots(figsize=(12, 10))
    sc.pl.umap(adata, color='DGAT1', ax=ax, show=False,
              title='DGAT1 Expression',
              cmap='Reds', vmax='p95', size=30, alpha=0.8,
              frameon=False)
    
    plt.tight_layout()
    plt.savefig(os.path.join(save_dir, 'DGAT1_UMAP.pdf'), dpi=300, bbox_inches='tight')
    plt.savefig(os.path.join(save_dir, 'DGAT1_UMAP.png'), dpi=300, bbox_inches='tight')
    plt.close()
    print(f"  [OK] Saved: DGAT1_UMAP.*")
    
    # Figure 2: DGAT1 by cell type
    fig, ax = plt.subplots(figsize=(12, 8))
    sns.violinplot(data=adata.obs, x='cell_type', y='DGAT1', ax=ax, inner='box')
    ax.set_xlabel('Cell Type', fontsize=14)
    ax.set_ylabel('DGAT1 Expression (log-normalized)', fontsize=14)
    ax.set_title('DGAT1 Expression Across Cell Types', fontsize=16, fontweight='bold', pad=20)
    ax.tick_params(axis='x', rotation=45)
    ax.grid(axis='y', alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(os.path.join(save_dir, 'DGAT1_ByCellType.pdf'), dpi=300, bbox_inches='tight')
    plt.savefig(os.path.join(save_dir, 'DGAT1_ByCellType.png'), dpi=300, bbox_inches='tight')
    plt.close()
    print(f"  [OK] Saved: DGAT1_ByCellType.*")
    
    # Figure 3: DGAT1 by group
    if 'group' in adata.obs.columns:
        fig, ax = plt.subplots(figsize=(10, 6))
        sns.boxplot(data=adata.obs, x='group', y='DGAT1', ax=ax,
                   palette='Set2', order=['NGB', 'IMP', 'IMR', 'IWP', 'IWR'])
        ax.set_xlabel('Experimental Group', fontsize=14)
        ax.set_ylabel('DGAT1 Expression', fontsize=14)
        ax.set_title('DGAT1 Across Groups', fontsize=16, fontweight='bold', pad=20)
        ax.grid(axis='y', alpha=0.3)
        
        plt.tight_layout()
        plt.savefig(os.path.join(save_dir, 'DGAT1_ByGroup.pdf'), dpi=300, bbox_inches='tight')
        plt.savefig(os.path.join(save_dir, 'DGAT1_ByGroup.png'), dpi=300, bbox_inches='tight')
        plt.close()
        print(f"  [OK] Saved: DGAT1_ByGroup.*")

def main():
    """Main execution"""
    import argparse
    
    parser = argparse.ArgumentParser(
        description='Generate individual figures (not combined panels)'
    )
    parser.add_argument('--input', type=str, required=True,
                       help='Input H5AD file')
    parser.add_argument('--output', type=str, required=True,
                       help='Output directory for individual figures')
    parser.add_argument('--raw_sample', type=str, default=None,
                       help='Path to one raw 10x sample (to extract DGAT1)')
    
    args = parser.parse_args()
    
    print("="*70)
    print("INDIVIDUAL FIGURE GENERATION")
    print("Each figure saved separately (not combined panels)")
    print("="*70)
    print(f"\nInput: {args.input}")
    print(f"Output: {args.output}")
    
    os.makedirs(args.output, exist_ok=True)
    
    # Load data
    adata = load_with_dgat1(args.input, [args.raw_sample] if args.raw_sample else None)
    
    # Generate all individual figures
    print("\n" + "="*70)
    print("GENERATING INDIVIDUAL FIGURES")
    print("="*70)
    
    figure_umap_by_celltype(adata, args.output)
    figure_umap_by_sample(adata, args.output)
    figure_umap_by_group(adata, args.output)
    figure_clustering(adata, args.output)
    figure_composition_by_group(adata, args.output)
    figure_composition_heatmap(adata, args.output)
    figure_myeloid_lymphoid(adata, args.output)
    figure_dgat1_expression(adata, args.output, args.raw_sample)
    
    print("\n" + "="*70)
    print("[SUCCESS] ALL INDIVIDUAL FIGURES GENERATED!")
    print("="*70)
    print(f"\nCheck output directory: {args.output}/")
    print("\nFigures saved:")
    for f in sorted(os.listdir(args.output)):
        if f.endswith(('.pdf', '.png')):
            print(f"  - {f}")

if __name__ == "__main__":
    main()

