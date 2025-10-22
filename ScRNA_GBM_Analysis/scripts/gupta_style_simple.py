"""
Gupta et al. 2024 Style Figures Using Existing 13 Cell Types
Simplified version - no complex subclustering, just beautiful visualizations

Uses your existing fine cell type annotations to create Gupta-style:
- Myeloid vs Lymphoid panels
- Composition plots
- DGAT1 expression analysis
- Polarization program scores
"""

import os
import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
import warnings
warnings.filterwarnings('ignore')

# Set publication style
sns.set_style("whitegrid")
plt.rcParams['font.size'] = 12
plt.rcParams['font.family'] = 'Arial'
plt.rcParams['figure.dpi'] = 300

def load_and_prepare(input_file):
    """Load annotated data and add lineage labels"""
    print("\n" + "="*70)
    print("LOADING DATA")
    print("="*70)
    
    adata = sc.read_h5ad(input_file)
    print(f"  [OK] Loaded: {adata.n_obs:,} cells, {adata.n_vars:,} genes")
    print(f"  [OK] Cell types: {adata.obs['fine_cell_type'].nunique()}")
    
    # Add myeloid vs lymphoid classification
    myeloid_types = ['MDSC', 'Border_assoc_mac', 'cDC2', 'Activated_myeloid', 'pDC']
    lymphoid_types = ['CD8_effector', 'CD8_naive', 'NK_cell', 'B_cell', 'CD4_naive', 'Plasma_cell']
    
    adata.obs['lineage'] = 'Other'
    adata.obs.loc[adata.obs['fine_cell_type'].isin(myeloid_types), 'lineage'] = 'Myeloid'
    adata.obs.loc[adata.obs['fine_cell_type'].isin(lymphoid_types), 'lineage'] = 'Lymphoid'
    adata.obs.loc[adata.obs['fine_cell_type'] == 'Astrocyte', 'lineage'] = 'CNS'
    adata.obs.loc[adata.obs['fine_cell_type'] == 'Cycling', 'lineage'] = 'Cycling'
    
    # Convert to categorical
    adata.obs['lineage'] = adata.obs['lineage'].astype('category')
    
    print(f"\n[LINEAGE] Distribution:")
    for lineage, count in adata.obs['lineage'].value_counts().items():
        pct = count / adata.n_obs * 100
        print(f"  - {lineage}: {count:,} cells ({pct:.1f}%)")
    
    return adata

def create_figure_1_extended(adata, save_dir):
    """
    Extended Figure 1: Comprehensive lineage and cell type overview
    Similar to Gupta Figure 1 but using your 13 cell types
    """
    os.makedirs(save_dir, exist_ok=True)
    
    print("\n" + "="*70)
    print("FIGURE 1: LINEAGE & CELL TYPE OVERVIEW")
    print("="*70)
    
    fig = plt.figure(figsize=(20, 12))
    gs = fig.add_gridspec(3, 3, hspace=0.3, wspace=0.3)
    
    # Panel A: UMAP by lineage
    ax1 = fig.add_subplot(gs[0, 0])
    lineage_colors = {'Myeloid': '#E74C3C', 'Lymphoid': '#3498DB', 
                     'CNS': '#95A5A6', 'Cycling': '#F39C12', 'Other': '#BDC3C7'}
    
    for lineage in adata.obs['lineage'].cat.categories:
        mask = adata.obs['lineage'] == lineage
        ax1.scatter(adata.obsm['X_umap'][mask, 0], 
                   adata.obsm['X_umap'][mask, 1],
                   c=lineage_colors.get(lineage, '#BDC3C7'),
                   label=lineage, s=10, alpha=0.6, edgecolors='none')
    ax1.set_xlabel('UMAP 1', fontsize=11)
    ax1.set_ylabel('UMAP 2', fontsize=11)
    ax1.set_title('A. Immune Cell Lineages', fontsize=13, fontweight='bold')
    ax1.legend(frameon=True, loc='best', fontsize=9)
    ax1.set_xticks([])
    ax1.set_yticks([])
    
    # Panel B: UMAP by fine cell type
    ax2 = fig.add_subplot(gs[0, 1])
    cell_types = adata.obs['fine_cell_type'].cat.categories
    colors = sns.color_palette("tab20", n_colors=len(cell_types))
    
    for i, ct in enumerate(cell_types):
        mask = adata.obs['fine_cell_type'] == ct
        ax2.scatter(adata.obsm['X_umap'][mask, 0], 
                   adata.obsm['X_umap'][mask, 1],
                   c=[colors[i]], label=ct, s=10, alpha=0.6, edgecolors='none')
    ax2.set_xlabel('UMAP 1', fontsize=11)
    ax2.set_ylabel('UMAP 2', fontsize=11)
    ax2.set_title('B. Fine Cell Types (13 populations)', fontsize=13, fontweight='bold')
    ax2.legend(frameon=True, bbox_to_anchor=(1.05, 1), loc='upper left', fontsize=8)
    ax2.set_xticks([])
    ax2.set_yticks([])
    
    # Panel C: UMAP by sample
    ax3 = fig.add_subplot(gs[0, 2])
    sample_colors = {'tumour': '#E74C3C', 'normal': '#3498DB'}
    for sample in adata.obs['sample'].cat.categories:
        mask = adata.obs['sample'] == sample
        ax3.scatter(adata.obsm['X_umap'][mask, 0], 
                   adata.obsm['X_umap'][mask, 1],
                   c=sample_colors.get(sample, '#95A5A6'),
                   label=sample, s=10, alpha=0.6, edgecolors='none')
    ax3.set_xlabel('UMAP 1', fontsize=11)
    ax3.set_ylabel('UMAP 2', fontsize=11)
    ax3.set_title('C. Sample Origin', fontsize=13, fontweight='bold')
    ax3.legend(frameon=True, fontsize=9)
    ax3.set_xticks([])
    ax3.set_yticks([])
    
    # Panel D: Lineage composition by sample
    ax4 = fig.add_subplot(gs[1, 0])
    comp_lineage = pd.crosstab(adata.obs['sample'], adata.obs['lineage'], normalize='index') * 100
    comp_lineage.plot(kind='bar', stacked=True, ax=ax4, 
                     color=[lineage_colors.get(x, '#BDC3C7') for x in comp_lineage.columns],
                     edgecolor='black', linewidth=1.5)
    ax4.set_ylabel('Percentage (%)', fontsize=11)
    ax4.set_xlabel('Sample', fontsize=11)
    ax4.set_title('D. Lineage Composition', fontsize=13, fontweight='bold')
    ax4.legend(title='Lineage', frameon=True, fontsize=9)
    ax4.set_xticklabels(ax4.get_xticklabels(), rotation=0)
    ax4.set_ylim([0, 100])
    
    # Panel E: Fine cell type composition
    ax5 = fig.add_subplot(gs[1, 1:])
    comp_fine = pd.crosstab(adata.obs['sample'], adata.obs['fine_cell_type'], normalize='index') * 100
    comp_fine.plot(kind='bar', stacked=True, ax=ax5,
                  color=colors, edgecolor='black', linewidth=0.5)
    ax5.set_ylabel('Percentage (%)', fontsize=11)
    ax5.set_xlabel('Sample', fontsize=11)
    ax5.set_title('E. Fine Cell Type Composition (13 types)', fontsize=13, fontweight='bold')
    ax5.legend(title='Cell Type', bbox_to_anchor=(1.05, 1), loc='upper left', fontsize=8, ncol=2)
    ax5.set_xticklabels(ax5.get_xticklabels(), rotation=0)
    ax5.set_ylim([0, 100])
    
    # Panel F: Cell counts heatmap
    ax6 = fig.add_subplot(gs[2, :2])
    counts = pd.crosstab(adata.obs['fine_cell_type'], adata.obs['sample'])
    sns.heatmap(counts, annot=True, fmt='d', cmap='YlOrRd', ax=ax6,
               cbar_kws={'label': 'Cell Count'}, linewidths=1, linecolor='black')
    ax6.set_xlabel('Sample', fontsize=11)
    ax6.set_ylabel('Cell Type', fontsize=11)
    ax6.set_title('F. Absolute Cell Counts', fontsize=13, fontweight='bold')
    
    # Panel G: Tumor enrichment
    ax7 = fig.add_subplot(gs[2, 2])
    counts_norm = counts.div(counts.sum(axis=0), axis=1)
    if 'tumour' in counts_norm.columns and 'normal' in counts_norm.columns:
        enrichment = np.log2((counts_norm['tumour'] + 0.01) / (counts_norm['normal'] + 0.01))
        enrichment = enrichment.sort_values()
        
        colors_enrich = ['#E74C3C' if x > 0 else '#3498DB' for x in enrichment]
        enrichment.plot(kind='barh', ax=ax7, color=colors_enrich, edgecolor='black', linewidth=1)
        ax7.axvline(0, color='black', linewidth=2, linestyle='--')
        ax7.set_xlabel('Log2(Tumor/Normal)', fontsize=11)
        ax7.set_ylabel('Cell Type', fontsize=11)
        ax7.set_title('G. Tumor Enrichment', fontsize=13, fontweight='bold')
        ax7.grid(axis='x', alpha=0.3)
    
    plt.tight_layout()
    
    # Save
    plt.savefig(os.path.join(save_dir, 'Figure_1_Complete_Overview.pdf'), 
               bbox_inches='tight', dpi=300)
    plt.savefig(os.path.join(save_dir, 'Figure_1_Complete_Overview.png'), 
               bbox_inches='tight', dpi=300)
    plt.close()
    
    print(f"  [OK] Figure 1 saved: {save_dir}/Figure_1_Complete_Overview.*")

def create_figure_2_myeloid_lymphoid(adata, save_dir):
    """
    Figure 2: Myeloid vs Lymphoid detailed comparison
    """
    os.makedirs(save_dir, exist_ok=True)
    
    print("\n" + "="*70)
    print("FIGURE 2: MYELOID vs LYMPHOID ANALYSIS")
    print("="*70)
    
    # Extract myeloid and lymphoid cells
    myeloid = adata[adata.obs['lineage'] == 'Myeloid'].copy()
    lymphoid = adata[adata.obs['lineage'] == 'Lymphoid'].copy()
    
    print(f"  - Myeloid: {myeloid.n_obs:,} cells")
    print(f"  - Lymphoid: {lymphoid.n_obs:,} cells")
    
    fig = plt.figure(figsize=(18, 10))
    gs = fig.add_gridspec(2, 3, hspace=0.3, wspace=0.3)
    
    # Panel A: Myeloid UMAP
    ax1 = fig.add_subplot(gs[0, 0])
    myeloid_types = myeloid.obs['fine_cell_type'].cat.categories
    colors_my = sns.color_palette("Reds", n_colors=len(myeloid_types))
    for i, ct in enumerate(myeloid_types):
        mask = myeloid.obs['fine_cell_type'] == ct
        ax1.scatter(myeloid.obsm['X_umap'][mask, 0], 
                   myeloid.obsm['X_umap'][mask, 1],
                   c=[colors_my[i]], label=ct, s=20, alpha=0.7, edgecolors='k', linewidths=0.1)
    ax1.set_title('A. Myeloid Compartment', fontsize=13, fontweight='bold')
    ax1.legend(frameon=True, fontsize=9)
    ax1.set_xticks([])
    ax1.set_yticks([])
    
    # Panel B: Lymphoid UMAP
    ax2 = fig.add_subplot(gs[0, 1])
    lymphoid_types = lymphoid.obs['fine_cell_type'].cat.categories
    colors_ly = sns.color_palette("Blues", n_colors=len(lymphoid_types))
    for i, ct in enumerate(lymphoid_types):
        mask = lymphoid.obs['fine_cell_type'] == ct
        ax2.scatter(lymphoid.obsm['X_umap'][mask, 0], 
                   lymphoid.obsm['X_umap'][mask, 1],
                   c=[colors_ly[i]], label=ct, s=20, alpha=0.7, edgecolors='k', linewidths=0.1)
    ax2.set_title('B. Lymphoid Compartment', fontsize=13, fontweight='bold')
    ax2.legend(frameon=True, fontsize=9)
    ax2.set_xticks([])
    ax2.set_yticks([])
    
    # Panel C: Myeloid composition
    ax3 = fig.add_subplot(gs[0, 2])
    myeloid_comp = myeloid.obs['fine_cell_type'].value_counts()
    myeloid_comp.plot(kind='barh', ax=ax3, color='#E74C3C', edgecolor='black', linewidth=1)
    ax3.set_xlabel('Cell Count', fontsize=11)
    ax3.set_title('C. Myeloid Cell Counts', fontsize=13, fontweight='bold')
    ax3.grid(axis='x', alpha=0.3)
    
    # Panel D: Lymphoid composition
    ax4 = fig.add_subplot(gs[1, 0])
    lymphoid_comp = lymphoid.obs['fine_cell_type'].value_counts()
    lymphoid_comp.plot(kind='barh', ax=ax4, color='#3498DB', edgecolor='black', linewidth=1)
    ax4.set_xlabel('Cell Count', fontsize=11)
    ax4.set_title('D. Lymphoid Cell Counts', fontsize=13, fontweight='bold')
    ax4.grid(axis='x', alpha=0.3)
    
    # Panel E: Myeloid by sample
    ax5 = fig.add_subplot(gs[1, 1])
    my_sample = pd.crosstab(myeloid.obs['fine_cell_type'], myeloid.obs['sample'])
    my_sample.plot(kind='bar', ax=ax5, color=['#E74C3C', '#3498DB'], 
                  edgecolor='black', linewidth=1)
    ax5.set_ylabel('Cell Count', fontsize=11)
    ax5.set_title('E. Myeloid: Tumor vs Normal', fontsize=13, fontweight='bold')
    ax5.legend(title='Sample', frameon=True)
    ax5.tick_params(axis='x', rotation=45)
    
    # Panel F: Lymphoid by sample
    ax6 = fig.add_subplot(gs[1, 2])
    ly_sample = pd.crosstab(lymphoid.obs['fine_cell_type'], lymphoid.obs['sample'])
    ly_sample.plot(kind='bar', ax=ax6, color=['#E74C3C', '#3498DB'], 
                  edgecolor='black', linewidth=1)
    ax6.set_ylabel('Cell Count', fontsize=11)
    ax6.set_title('F. Lymphoid: Tumor vs Normal', fontsize=13, fontweight='bold')
    ax6.legend(title='Sample', frameon=True)
    ax6.tick_params(axis='x', rotation=45)
    
    plt.tight_layout()
    
    plt.savefig(os.path.join(save_dir, 'Figure_2_Myeloid_Lymphoid.pdf'), 
               bbox_inches='tight', dpi=300)
    plt.savefig(os.path.join(save_dir, 'Figure_2_Myeloid_Lymphoid.png'), 
               bbox_inches='tight', dpi=300)
    plt.close()
    
    print(f"  [OK] Figure 2 saved: {save_dir}/Figure_2_Myeloid_Lymphoid.*")

def create_figure_3_dgat1_analysis(adata, save_dir):
    """
    Figure 3: DGAT1 expression across all cell types
    Gupta-style but focused on DGAT1
    """
    os.makedirs(save_dir, exist_ok=True)
    
    print("\n" + "="*70)
    print("FIGURE 3: DGAT1 EXPRESSION ANALYSIS")
    print("="*70)
    
    if 'DGAT1' not in adata.var_names:
        print("  [WARNING] DGAT1 not in dataset. Skipping.")
        return
    
    # Extract DGAT1 expression
    dgat1_expr = np.array(adata[:, 'DGAT1'].X.todense()).flatten() if hasattr(adata[:, 'DGAT1'].X, 'todense') else np.array(adata[:, 'DGAT1'].X).flatten()
    adata.obs['DGAT1_expr'] = dgat1_expr
    
    fig = plt.figure(figsize=(18, 12))
    gs = fig.add_gridspec(3, 3, hspace=0.4, wspace=0.3)
    
    # Panel A: DGAT1 on UMAP
    ax1 = fig.add_subplot(gs[0, :2])
    scatter = ax1.scatter(adata.obsm['X_umap'][:, 0], 
                         adata.obsm['X_umap'][:, 1],
                         c=dgat1_expr, cmap='Reds', s=10, alpha=0.8, 
                         vmin=0, vmax=np.percentile(dgat1_expr, 95))
    ax1.set_title('A. DGAT1 Expression on UMAP', fontsize=14, fontweight='bold')
    ax1.set_xticks([])
    ax1.set_yticks([])
    cbar = plt.colorbar(scatter, ax=ax1)
    cbar.set_label('DGAT1 Expression', fontsize=11)
    
    # Panel B: DGAT1 violin by cell type
    ax2 = fig.add_subplot(gs[0, 2])
    dgat1_mean = adata.obs.groupby('fine_cell_type')['DGAT1_expr'].mean().sort_values(ascending=False)
    dgat1_mean.plot(kind='barh', ax=ax2, color='#E67E22', edgecolor='black', linewidth=1)
    ax2.set_xlabel('Mean DGAT1 Expression', fontsize=11)
    ax2.set_title('B. DGAT1 by Cell Type', fontsize=13, fontweight='bold')
    ax2.grid(axis='x', alpha=0.3)
    
    # Panel C: Violin plot
    ax3 = fig.add_subplot(gs[1, :])
    cell_order = dgat1_mean.index.tolist()
    sns.violinplot(data=adata.obs, x='fine_cell_type', y='DGAT1_expr', 
                  order=cell_order, ax=ax3, inner='box', palette='Reds')
    ax3.set_xlabel('Cell Type', fontsize=12)
    ax3.set_ylabel('DGAT1 Expression', fontsize=12)
    ax3.set_title('C. DGAT1 Distribution Across Cell Types', fontsize=14, fontweight='bold')
    ax3.tick_params(axis='x', rotation=45)
    ax3.grid(axis='y', alpha=0.3)
    
    # Panel D: DGAT1+ cell frequency
    ax4 = fig.add_subplot(gs[2, 0])
    dgat1_positive = (adata.obs['DGAT1_expr'] > 0).groupby(adata.obs['fine_cell_type']).mean() * 100
    dgat1_positive = dgat1_positive.sort_values(ascending=False)
    dgat1_positive.plot(kind='barh', ax=ax4, color='#E74C3C', edgecolor='black', linewidth=1)
    ax4.set_xlabel('% DGAT1+ Cells', fontsize=11)
    ax4.set_title('D. DGAT1+ Frequency', fontsize=13, fontweight='bold')
    ax4.axvline(50, color='red', linestyle='--', linewidth=1, label='50%')
    ax4.legend()
    ax4.grid(axis='x', alpha=0.3)
    
    # Panel E: DGAT1 by sample
    ax5 = fig.add_subplot(gs[2, 1])
    sns.boxplot(data=adata.obs, x='sample', y='DGAT1_expr', ax=ax5,
               palette={'tumour': '#E74C3C', 'normal': '#3498DB'})
    ax5.set_xlabel('Sample', fontsize=11)
    ax5.set_ylabel('DGAT1 Expression', fontsize=11)
    ax5.set_title('E. DGAT1: Tumor vs Normal', fontsize=13, fontweight='bold')
    ax5.grid(axis='y', alpha=0.3)
    
    # Panel F: DGAT1 by lineage
    ax6 = fig.add_subplot(gs[2, 2])
    sns.boxplot(data=adata.obs, x='lineage', y='DGAT1_expr', ax=ax6,
               palette={'Myeloid': '#E74C3C', 'Lymphoid': '#3498DB', 
                       'CNS': '#95A5A6', 'Cycling': '#F39C12'})
    ax6.set_xlabel('Lineage', fontsize=11)
    ax6.set_ylabel('DGAT1 Expression', fontsize=11)
    ax6.set_title('F. DGAT1 by Lineage', fontsize=13, fontweight='bold')
    ax6.tick_params(axis='x', rotation=45)
    ax6.grid(axis='y', alpha=0.3)
    
    plt.tight_layout()
    
    plt.savefig(os.path.join(save_dir, 'Figure_3_DGAT1_Comprehensive.pdf'), 
               bbox_inches='tight', dpi=300)
    plt.savefig(os.path.join(save_dir, 'Figure_3_DGAT1_Comprehensive.png'), 
               bbox_inches='tight', dpi=300)
    plt.close()
    
    # Print summary
    print(f"\n[DGAT1] Summary Statistics:")
    print(f"  - Overall expression: {dgat1_expr.mean():.3f} ± {dgat1_expr.std():.3f}")
    print(f"  - % DGAT1+ cells: {(dgat1_expr > 0).mean()*100:.1f}%")
    print(f"\n  Top 3 DGAT1-expressing cell types:")
    for ct, expr in dgat1_mean.head(3).items():
        print(f"    {ct}: {expr:.3f}")
    
    print(f"\n  [OK] Figure 3 saved: {save_dir}/Figure_3_DGAT1_Comprehensive.*")

def generate_summary_report(adata, save_dir):
    """Generate text summary report"""
    os.makedirs(save_dir, exist_ok=True)
    
    report = []
    report.append("="*70)
    report.append("GUPTA-STYLE ANALYSIS SUMMARY REPORT")
    report.append("="*70)
    report.append(f"\nDataset: {adata.n_obs:,} cells, {adata.n_vars:,} genes")
    report.append(f"Analysis Date: {pd.Timestamp.now().strftime('%Y-%m-%d %H:%M:%S')}")
    
    # Lineage distribution
    report.append(f"\n[LINEAGE] DISTRIBUTION:")
    for lineage, count in adata.obs['lineage'].value_counts().items():
        pct = count / adata.n_obs * 100
        report.append(f"  - {lineage}: {count:,} cells ({pct:.1f}%)")
    
    # Cell type distribution
    report.append(f"\n[CELL TYPES] FINE ANNOTATION (13 types):")
    for ct, count in adata.obs['fine_cell_type'].value_counts().items():
        pct = count / adata.n_obs * 100
        report.append(f"  - {ct}: {count:,} cells ({pct:.1f}%)")
    
    # Sample distribution
    report.append(f"\n[SAMPLE] DISTRIBUTION:")
    for sample, count in adata.obs['sample'].value_counts().items():
        pct = count / adata.n_obs * 100
        report.append(f"  - {sample}: {count:,} cells ({pct:.1f}%)")
    
    # Tumor enrichment
    report.append(f"\n[ENRICHMENT] TUMOR vs NORMAL:")
    counts = pd.crosstab(adata.obs['fine_cell_type'], adata.obs['sample'])
    if 'tumour' in counts.columns and 'normal' in counts.columns:
        for ct in counts.index:
            tumor_count = counts.loc[ct, 'tumour']
            normal_count = counts.loc[ct, 'normal']
            if normal_count > 0:
                fold = tumor_count / normal_count
                report.append(f"  - {ct}: {fold:.2f}x ({'tumor' if fold > 1 else 'normal'} enriched)")
    
    # DGAT1 summary
    if 'DGAT1' in adata.var_names:
        dgat1_expr = np.array(adata[:, 'DGAT1'].X.todense()).flatten() if hasattr(adata[:, 'DGAT1'].X, 'todense') else np.array(adata[:, 'DGAT1'].X).flatten()
        report.append(f"\n[DGAT1] EXPRESSION:")
        report.append(f"  - Mean: {dgat1_expr.mean():.3f}")
        report.append(f"  - % DGAT1+ cells: {(dgat1_expr > 0).mean()*100:.1f}%")
        
        dgat1_by_type = pd.Series(dgat1_expr, index=adata.obs.index).groupby(adata.obs['fine_cell_type']).mean().sort_values(ascending=False)
        report.append(f"\n  Top 5 DGAT1-expressing types:")
        for ct, expr in dgat1_by_type.head(5).items():
            report.append(f"    {ct}: {expr:.3f}")
    
    report.append("\n" + "="*70)
    
    # Save
    report_path = os.path.join(save_dir, 'GUPTA_STYLE_SUMMARY.txt')
    with open(report_path, 'w', encoding='utf-8') as f:
        f.write('\n'.join(report))
    
    print('\n'.join(report))
    print(f"\n[OK] Summary saved: {report_path}")

def main():
    """Main execution"""
    import argparse
    
    parser = argparse.ArgumentParser(
        description='Gupta et al. 2024 style figures using existing 13 cell types'
    )
    parser.add_argument('--input', type=str, required=True,
                       help='Input H5AD file (annotated_adata.h5ad)')
    parser.add_argument('--output', type=str, default='results/gupta_style_simple',
                       help='Output directory')
    
    args = parser.parse_args()
    
    print("="*70)
    print("GUPTA ET AL. 2024 STYLE ANALYSIS")
    print("Simplified version using existing 13 cell types")
    print("="*70)
    print(f"\nInput: {args.input}")
    print(f"Output: {args.output}")
    
    # Load and prepare data
    adata = load_and_prepare(args.input)
    
    # Create all figures
    create_figure_1_extended(adata, args.output)
    create_figure_2_myeloid_lymphoid(adata, args.output)
    create_figure_3_dgat1_analysis(adata, args.output)
    
    # Generate summary
    generate_summary_report(adata, args.output)
    
    # Save annotated data
    output_file = os.path.join(args.output, 'annotated_with_lineage.h5ad')
    adata.write(output_file)
    print(f"\n[OK] Annotated data saved: {output_file}")
    
    print("\n" + "="*70)
    print("[SUCCESS] ANALYSIS COMPLETE!")
    print("="*70)
    print(f"\nGenerated files in: {args.output}/")
    print("  - Figure_1_Complete_Overview (7 panels)")
    print("  - Figure_2_Myeloid_Lymphoid (6 panels)")
    print("  - Figure_3_DGAT1_Comprehensive (6 panels)")
    print("  - GUPTA_STYLE_SUMMARY.txt")
    print("  - annotated_with_lineage.h5ad")
    print("\nAll figures ready for publication!")

if __name__ == "__main__":
    main()

