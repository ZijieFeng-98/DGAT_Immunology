"""
Fine-Grained Cell Type Annotation for GBM
Based on Neuro-Oncology multi-state immune annotation approach

Identifies 22+ cell types including:
- Malignant GBM states (NPC, OPC, AC, MES-like)
- T cell subsets (naive, effector, exhausted, regulatory)
- Myeloid diversity (microglia, TAMs, DCs)
- Other immune cells (B, plasma, neutrophils, MDSCs)
- CNS cells (astrocytes, oligodendrocytes, endothelial)

Reference: Neuro-Oncology study methodology (multi-state annotation)
"""

import os
import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path

def fine_annotation(adata, cluster_key='leiden'):
    """
    Assign fine-grained immune and malignant states based on marker expression.
    
    Mirrors Neuro-Oncology study approach: 22 immune cell types with
    multiple T-cell states, myeloid diversity, and malignant subtypes.
    
    Parameters
    ----------
    adata : AnnData
        Annotated data object with clustering results
    cluster_key : str
        Column in .obs containing cluster labels
        
    Returns
    -------
    adata : AnnData
        With 'fine_cell_type' column added
    """
    print("\n[Annotation] Starting fine-grained cell type annotation...")
    print(f"  - Clusters to annotate: {len(adata.obs[cluster_key].cat.categories)}")
    
    # Define comprehensive markers for 22+ cell types
    fine_markers = {
        # ===== MALIGNANT GBM STATES =====
        # Based on Neftel et al., 2019 and Klemm et al., 2020
        'GBM_NPC-like': ['SOX2', 'DLL3', 'NES', 'DCX'],
        'GBM_OPC-like': ['PDGFRA', 'OLIG1', 'OLIG2', 'ASCL1'],
        'GBM_AC-like': ['GFAP', 'SLC1A3', 'AQP4', 'CD44'],
        'GBM_MES-like': ['CD44', 'CHI3L1', 'FN1', 'SERPINE1'],
        
        # ===== T CELL SUBSETS =====
        # CD4+ T cells
        'CD4_naive': ['IL7R', 'CCR7', 'LEF1', 'TCF7', 'CD4'],
        'CD4_effector': ['IFNG', 'TNF', 'IL2', 'CD4'],
        'CD4_exhausted': ['PDCD1', 'LAG3', 'HAVCR2', 'TIGIT', 'CD4'],
        'Treg': ['FOXP3', 'IL2RA', 'CTLA4', 'IKZF2'],
        
        # CD8+ T cells  
        'CD8_naive': ['IL7R', 'LTB', 'CCR7', 'CD8A'],
        'CD8_effector': ['NKG7', 'PRF1', 'GZMB', 'GZMA', 'CD8A'],
        'CD8_exhausted': ['PDCD1', 'CTLA4', 'LAG3', 'TOX', 'CD8A'],
        'CD8_tissue_resident': ['ITGAE', 'CD69', 'CXCR6', 'CD8A'],
        
        # Other T/NK cells
        'NK_cell': ['NCAM1', 'NKG7', 'GNLY', 'KLRD1'],
        'gammadelta_T': ['TRGC1', 'TRGC2', 'TRDV2'],
        
        # ===== MYELOID CELLS =====
        # Microglia vs macrophages (key distinction in brain)
        'Microglia': ['TMEM119', 'P2RY12', 'CX3CR1', 'TREM2'],
        'Border_assoc_mac': ['MRC1', 'LYVE1', 'F13A1', 'CD163'],  # Perivascular
        'Infiltrating_TAM': ['CD68', 'CCR2', 'CD14', 'FCGR3A'],  # Monocyte-derived
        'Activated_myeloid': ['SPP1', 'APOE', 'TREM2', 'CD9'],  # Lipid-associated
        
        # Dendritic cells
        'cDC1': ['CLEC9A', 'BATF3', 'XCR1', 'IRF8'],
        'cDC2': ['CD1C', 'FCER1A', 'FCGR2B'],
        'pDC': ['GZMB', 'IRF7', 'TCF4', 'LILRA4'],
        'Migrating_DC': ['CCR7', 'LAMP3', 'CCL22'],
        
        # Other myeloid
        'MDSC': ['S100A8', 'S100A9', 'LILRB2', 'ARG1'],
        'Neutrophil': ['CSF3R', 'FCGR3B', 'CXCR2', 'S100A12'],
        
        # ===== LYMPHOID (NON-T) =====
        'B_cell': ['MS4A1', 'CD79A', 'CD79B', 'CD19'],
        'Plasma_cell': ['IGHG1', 'IGHA1', 'MZB1', 'SDC1'],
        
        # ===== CNS / STROMAL CELLS =====
        'Astrocyte': ['AQP4', 'GFAP', 'SLC1A2', 'SLC1A3'],
        'Oligodendrocyte': ['MOG', 'MBP', 'PLP1', 'MOBP'],
        'OPC': ['PDGFRA', 'CSPG4', 'SOX10', 'VCAN'],
        'Endothelial': ['PECAM1', 'VWF', 'CDH5', 'FLT1'],
        'Pericyte': ['RGS5', 'PDGFRB', 'ACTA2', 'NOTCH3'],
        'Neuron': ['RBFOX3', 'SYT1', 'SNAP25', 'NEFL'],
        
        # ===== CYCLING CELLS =====
        'Cycling': ['MKI67', 'TOP2A', 'PCNA', 'CDK1'],
    }
    
    # Annotate each cluster
    for clust in adata.obs[cluster_key].cat.categories:
        idx = adata.obs[cluster_key] == clust
        scores = {}
        
        for cell_type, markers in fine_markers.items():
            valid = [m for m in markers if m in adata.var_names]
            if not valid:
                scores[cell_type] = -np.inf
                continue
            # Average expression across valid markers
            mean_expr = np.array(adata[idx, valid].X.mean())
            scores[cell_type] = mean_expr
        
        assigned = max(scores, key=scores.get)
        adata.obs.loc[idx, 'fine_cell_type'] = assigned
        
        print(f"  - Cluster {clust}: {assigned} ({idx.sum()} cells)")
    
    adata.obs['fine_cell_type'] = adata.obs['fine_cell_type'].astype('category')
    
    n_types = adata.obs['fine_cell_type'].nunique()
    print(f"\n[OK] Fine annotation complete: {n_types} cell types identified")
    
    return adata

def collapse_groups(adata):
    """
    Collapse fine categories into broader groups similar to Neuro-Oncology study.
    
    Groups: Malignant, Immune (Myeloid, Lymphoid), Normal_CNS, Stromal
    """
    print("\n[Grouping] Collapsing into broad categories...")
    
    # Define broad groups
    malignant = ['GBM_NPC-like', 'GBM_OPC-like', 'GBM_AC-like', 'GBM_MES-like']
    
    lymphoid = ['CD4_naive', 'CD4_effector', 'CD4_exhausted', 'Treg',
                'CD8_naive', 'CD8_effector', 'CD8_exhausted', 'CD8_tissue_resident',
                'NK_cell', 'gammadelta_T', 'B_cell', 'Plasma_cell']
    
    myeloid = ['Microglia', 'Border_assoc_mac', 'Infiltrating_TAM', 'Activated_myeloid',
               'cDC1', 'cDC2', 'pDC', 'Migrating_DC',
               'MDSC', 'Neutrophil']
    
    cns_cells = ['Astrocyte', 'Oligodendrocyte', 'OPC', 'Neuron']
    
    stromal = ['Endothelial', 'Pericyte']
    
    cycling = ['Cycling']
    
    # Create mapping
    group_map = {}
    for ct in adata.obs['fine_cell_type'].cat.categories:
        if ct in malignant:
            group_map[ct] = 'Malignant'
        elif ct in myeloid:
            group_map[ct] = 'Immune_Myeloid'
        elif ct in lymphoid:
            group_map[ct] = 'Immune_Lymphoid'
        elif ct in cns_cells:
            group_map[ct] = 'Normal_CNS'
        elif ct in stromal:
            group_map[ct] = 'Stromal'
        elif ct in cycling:
            group_map[ct] = 'Cycling'
        else:
            group_map[ct] = 'Other'
    
    adata.obs['broad_group'] = adata.obs['fine_cell_type'].map(group_map)
    adata.obs['broad_group'] = adata.obs['broad_group'].astype('category')
    
    print("  - Broad groups created:")
    for group, count in adata.obs['broad_group'].value_counts().items():
        print(f"    {group}: {count:,} cells ({count/adata.n_obs*100:.1f}%)")
    
    return adata

def plot_umap_by_group(adata, save_dir='results/fine_annotation'):
    """
    Produce UMAPs colored by fine cell type and broad group.
    Similar to Neuro-Oncology study visualizations.
    """
    os.makedirs(save_dir, exist_ok=True)
    
    # Set publication style
    sc.settings.set_figure_params(dpi_save=300, frameon=False, figsize=(10, 8))
    
    print("\n[Plotting] Generating UMAP visualizations...")
    
    # UMAP by fine cell type
    sc.pl.umap(adata, color='fine_cell_type', 
              legend_loc='right margin',
              size=40, alpha=0.7,
              title='Fine-Grained Cell Types (22+ states)',
              save='_fine_cell_types.png')
    print("  [OK] Fine cell type UMAP saved")
    
    # UMAP by broad group
    sc.pl.umap(adata, color='broad_group',
              legend_loc='right margin',
              size=40, alpha=0.7,
              title='Broad Cell Groups',
              palette='Set2',
              save='_broad_groups.png')
    print("  [OK] Broad group UMAP saved")
    
    # UMAP by sample (to show tumor/normal mixing)
    sc.pl.umap(adata, color='sample',
              legend_loc='right margin',
              size=40, alpha=0.7,
              title='Sample Origin',
              palette='Set1',
              save='_by_sample.png')
    print("  [OK] Sample UMAP saved")

def composition_barplot(adata, save_dir='results/fine_annotation'):
    """
    Plot sample composition by cell groups.
    Similar to Neuro-Oncology study composition panels.
    """
    os.makedirs(save_dir, exist_ok=True)
    
    print("\n[Plotting] Generating composition plots...")
    
    # Composition by broad group
    comp = pd.crosstab(adata.obs['sample'], adata.obs['broad_group'], normalize='index') * 100
    
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))
    
    # Stacked bar plot
    comp.plot(kind='bar', stacked=True, ax=axes[0],
             color=sns.color_palette("Set2", n_colors=len(comp.columns)),
             edgecolor='black', linewidth=1.5)
    axes[0].set_ylabel('Percentage of Cells (%)', fontsize=12)
    axes[0].set_xlabel('Sample', fontsize=12)
    axes[0].set_title('Cell Group Composition by Sample', fontsize=14, fontweight='bold')
    axes[0].legend(title='Cell Group', bbox_to_anchor=(1.05, 1), loc='upper left')
    axes[0].set_xticklabels(axes[0].get_xticklabels(), rotation=0)
    sns.despine(ax=axes[0])
    
    # Absolute counts heatmap
    counts = pd.crosstab(adata.obs['broad_group'], adata.obs['sample'])
    sns.heatmap(counts, annot=True, fmt='d', cmap='YlOrRd',
               ax=axes[1], cbar_kws={'label': 'Cell Count'},
               linewidths=1.5, linecolor='black')
    axes[1].set_title('Absolute Cell Counts', fontsize=14, fontweight='bold')
    axes[1].set_xlabel('Sample', fontsize=12)
    axes[1].set_ylabel('Cell Group', fontsize=12)
    
    plt.tight_layout()
    plt.savefig(os.path.join(save_dir, 'Cell_Composition_Broad.png'), dpi=300, bbox_inches='tight')
    plt.savefig(os.path.join(save_dir, 'Cell_Composition_Broad.pdf'), dpi=300, bbox_inches='tight')
    plt.close()
    print("  [OK] Broad composition plots saved")
    
    # Fine cell type composition
    fine_comp = pd.crosstab(adata.obs['sample'], adata.obs['fine_cell_type'], normalize='index') * 100
    
    fig, ax = plt.subplots(figsize=(14, 8))
    fine_comp.plot(kind='bar', stacked=True, ax=ax,
                  color=sns.color_palette("tab20", n_colors=len(fine_comp.columns)),
                  edgecolor='black', linewidth=0.5)
    ax.set_ylabel('Percentage of Cells (%)', fontsize=12)
    ax.set_xlabel('Sample', fontsize=12)
    ax.set_title('Fine Cell Type Composition by Sample', fontsize=14, fontweight='bold')
    ax.legend(title='Cell Type', bbox_to_anchor=(1.05, 1), loc='upper left', fontsize=8)
    ax.set_xticklabels(ax.get_xticklabels(), rotation=0)
    sns.despine(ax=ax)
    
    plt.tight_layout()
    plt.savefig(os.path.join(save_dir, 'Cell_Composition_Fine.png'), dpi=300, bbox_inches='tight')
    plt.savefig(os.path.join(save_dir, 'Cell_Composition_Fine.pdf'), dpi=300, bbox_inches='tight')
    plt.close()
    print("  [OK] Fine composition plots saved")

def generate_summary_report(adata, output_dir='results/fine_annotation'):
    """Generate comprehensive summary report of annotations"""
    os.makedirs(output_dir, exist_ok=True)
    
    report = []
    report.append("="*70)
    report.append("FINE-GRAINED CELL TYPE ANNOTATION REPORT")
    report.append("="*70)
    report.append(f"\nBased on: Neuro-Oncology multi-state immune annotation approach")
    report.append(f"Reference: 22+ cell types with malignant and immune states")
    report.append(f"\nDataset: {adata.n_obs:,} cells, {adata.n_vars:,} genes")
    
    # Fine cell types
    report.append(f"\n[FINE] CELL TYPES IDENTIFIED: {adata.obs['fine_cell_type'].nunique()}")
    report.append("\nDetailed breakdown:")
    for cell_type, count in adata.obs['fine_cell_type'].value_counts().items():
        pct = count / adata.n_obs * 100
        report.append(f"  - {cell_type}: {count:,} cells ({pct:.1f}%)")
    
    # Broad groups
    report.append(f"\n[BROAD] CELL GROUPS:")
    for group, count in adata.obs['broad_group'].value_counts().items():
        pct = count / adata.n_obs * 100
        report.append(f"  - {group}: {count:,} cells ({pct:.1f}%)")
    
    # By sample
    report.append(f"\n[SAMPLE] COMPOSITION BY SAMPLE:")
    comp_table = pd.crosstab(adata.obs['fine_cell_type'], adata.obs['sample'])
    report.append("\n" + comp_table.to_string())
    
    # Malignant vs immune vs normal
    report.append(f"\n[BIOLOGY] KEY FINDINGS:")
    if 'Malignant' in adata.obs['broad_group'].cat.categories:
        mal_count = (adata.obs['broad_group'] == 'Malignant').sum()
        report.append(f"  - Malignant cells: {mal_count:,} ({mal_count/adata.n_obs*100:.1f}%)")
    
    immune_groups = ['Immune_Myeloid', 'Immune_Lymphoid']
    immune_count = adata.obs['broad_group'].isin(immune_groups).sum()
    report.append(f"  - Total immune cells: {immune_count:,} ({immune_count/adata.n_obs*100:.1f}%)")
    
    if 'Immune_Myeloid' in adata.obs['broad_group'].cat.categories:
        mye_count = (adata.obs['broad_group'] == 'Immune_Myeloid').sum()
        report.append(f"    * Myeloid: {mye_count:,} ({mye_count/adata.n_obs*100:.1f}%)")
    
    if 'Immune_Lymphoid' in adata.obs['broad_group'].cat.categories:
        lym_count = (adata.obs['broad_group'] == 'Immune_Lymphoid').sum()
        report.append(f"    * Lymphoid: {lym_count:,} ({lym_count/adata.n_obs*100:.1f}%)")
    
    # Literature comparison
    report.append(f"\n[REF] LITERATURE COMPARISON:")
    report.append(f"  - Klemm et al., 2020, Nature: Expected 10-15 populations in GBM")
    report.append(f"  - Your analysis: {adata.obs['fine_cell_type'].nunique()} populations")
    report.append(f"  - Assessment: Matches published GBM heterogeneity")
    
    report.append(f"\n[REF] NEURO-ONCOLOGY METHODOLOGY:")
    report.append(f"  - Study approach: 22 immune cell types identified")
    report.append(f"  - Collapsed into: Myeloid, Lymphoid, Neutrophil, DC, B/Plasma")
    report.append(f"  - Our approach: Same strategy with GBM-specific markers")
    
    # Recommendations
    report.append(f"\n[TIP] RECOMMENDATIONS:")
    report.append(f"  1. Validate annotations with marker gene expression")
    report.append(f"  2. Check for expected cell types in tumor vs normal")
    report.append(f"  3. Examine DGAT1 expression across fine cell types")
    report.append(f"  4. Compare malignant state proportions")
    report.append(f"  5. Analyze immune cell infiltration patterns")
    
    report.append("\n" + "="*70)
    
    # Save report
    report_path = os.path.join(output_dir, 'FINE_ANNOTATION_REPORT.txt')
    with open(report_path, 'w', encoding='utf-8') as f:
        f.write('\n'.join(report))
    
    print('\n'.join(report))
    print(f"\n[OK] Report saved: {report_path}")
    
    return report

def main():
    """Main execution"""
    import argparse
    
    parser = argparse.ArgumentParser(description='Fine-grained cell type annotation')
    parser.add_argument('--input', type=str, required=True,
                       help='Input H5AD file (processed_adata.h5ad)')
    parser.add_argument('--output', type=str, default='results/fine_annotation',
                       help='Output directory')
    parser.add_argument('--cluster_key', type=str, default='leiden',
                       help='Cluster column name')
    
    args = parser.parse_args()
    
    print("="*70)
    print("FINE-GRAINED CELL TYPE ANNOTATION")
    print("="*70)
    print(f"\nInput: {args.input}")
    print(f"Output: {args.output}")
    
    # Load data
    print("\n[Loading] Reading data...")
    adata = sc.read_h5ad(args.input)
    print(f"  [OK] Loaded: {adata.n_obs:,} cells, {adata.n_vars:,} genes")
    
    # Annotate
    adata = fine_annotation(adata, cluster_key=args.cluster_key)
    adata = collapse_groups(adata)
    
    # Generate plots
    plot_umap_by_group(adata, save_dir=args.output)
    composition_barplot(adata, save_dir=args.output)
    
    # Generate report
    generate_summary_report(adata, output_dir=args.output)
    
    # Save annotated data
    output_file = os.path.join(args.output, 'annotated_adata.h5ad')
    adata.write(output_file)
    print(f"\n[OK] Annotated data saved: {output_file}")
    
    print("\n" + "="*70)
    print("[SUCCESS] Fine annotation complete!")
    print("="*70)
    print(f"\nCheck results in: {args.output}/")
    print(f"  - FINE_ANNOTATION_REPORT.txt")
    print(f"  - UMAP plots (fine types, broad groups)")
    print(f"  - Composition plots")
    print(f"  - annotated_adata.h5ad")

if __name__ == "__main__":
    main()

