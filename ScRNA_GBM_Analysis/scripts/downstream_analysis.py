#!/usr/bin/env python3
"""
Downstream Analysis Examples for GBM Single-Cell Pipeline
==========================================================

This script demonstrates advanced analyses that can be performed
after running the main pipeline, including:
- Trajectory inference
- Cell-cell communication
- Survival analysis correlation
- Advanced metabolic profiling

Author: Generated for GBM DGAT1 Project
"""

import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
import warnings
warnings.filterwarnings('ignore')

# Set plotting style
sns.set_style('whitegrid')
sc.settings.set_figure_params(dpi=150, facecolor='white')


def load_processed_data(adata_path: str):
    """Load the processed AnnData from main pipeline."""
    print(f"Loading processed data from {adata_path}...")
    adata = sc.read_h5ad(adata_path)
    print(f"Loaded {adata.n_obs} cells and {adata.n_vars} genes")
    return adata


def trajectory_analysis(adata, output_dir='results/trajectory/'):
    """
    Perform trajectory inference on immune cells.
    
    Useful for understanding T cell differentiation/exhaustion
    or myeloid cell polarization trajectories.
    """
    print("\n" + "="*80)
    print("TRAJECTORY ANALYSIS")
    print("="*80)
    
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Example: T cell trajectory
    print("Analyzing T cell differentiation trajectory...")
    
    # Subset to T cells
    t_cell_types = ['T_cells', 'Cytotoxic_T', 'Regulatory_T', 'Exhausted_T']
    t_mask = adata.obs['cell_type'].isin(t_cell_types)
    
    if t_mask.sum() < 50:
        print("Not enough T cells for trajectory analysis")
        return
    
    adata_tcells = adata[t_mask, :].copy()
    
    # Recompute neighbors and UMAP for subset
    sc.pp.neighbors(adata_tcells, n_neighbors=15, n_pcs=30)
    sc.tl.umap(adata_tcells)
    
    # PAGA (Partition-based graph abstraction)
    print("Computing PAGA trajectory...")
    sc.tl.paga(adata_tcells, groups='cell_type')
    
    sc.pl.paga(adata_tcells, show=False)
    plt.savefig(output_dir / 'paga_trajectory.png', dpi=300, bbox_inches='tight')
    plt.close()
    
    # Diffusion pseudotime
    print("Computing diffusion pseudotime...")
    
    # Set root as naive/cytotoxic T cells
    root_cluster = adata_tcells.obs.groupby('leiden')['cell_type'].apply(
        lambda x: (x == 'Cytotoxic_T').sum()
    ).idxmax()
    
    root_cells = adata_tcells.obs['leiden'] == root_cluster
    adata_tcells.uns['iroot'] = np.flatnonzero(root_cells)[0]
    
    sc.tl.dpt(adata_tcells)
    
    # Plot pseudotime
    fig, axes = plt.subplots(1, 3, figsize=(18, 5))
    
    sc.pl.umap(adata_tcells, color='dpt_pseudotime', ax=axes[0], show=False)
    sc.pl.umap(adata_tcells, color='cell_type', ax=axes[1], show=False)
    sc.pl.umap(adata_tcells, color='DGAT1', ax=axes[2], show=False, use_raw=False)
    
    plt.tight_layout()
    plt.savefig(output_dir / 'tcell_pseudotime.png', dpi=300, bbox_inches='tight')
    plt.close()
    
    # DGAT1 expression along pseudotime
    plt.figure(figsize=(10, 6))
    
    # Create bins for pseudotime
    n_bins = 50
    pseudotime_bins = pd.cut(adata_tcells.obs['dpt_pseudotime'], bins=n_bins)
    
    dgat1_expr = adata_tcells[:, 'DGAT1'].X.toarray().flatten() if hasattr(adata_tcells[:, 'DGAT1'].X, 'toarray') else adata_tcells[:, 'DGAT1'].X.flatten()
    
    bin_means = pd.DataFrame({
        'pseudotime': adata_tcells.obs['dpt_pseudotime'],
        'DGAT1': dgat1_expr
    }).groupby(pseudotime_bins)['DGAT1'].mean()
    
    plt.plot(range(len(bin_means)), bin_means.values, linewidth=2)
    plt.xlabel('Pseudotime (binned)')
    plt.ylabel('DGAT1 Expression')
    plt.title('DGAT1 Expression Along T Cell Trajectory')
    plt.tight_layout()
    plt.savefig(output_dir / 'dgat1_along_pseudotime.png', dpi=300, bbox_inches='tight')
    plt.close()
    
    print(f"Trajectory analysis complete. Results saved to {output_dir}")


def cellcell_communication(adata, output_dir='results/communication/'):
    """
    Infer cell-cell communication networks.
    
    Focuses on tumor-immune interactions and metabolic signaling.
    """
    print("\n" + "="*80)
    print("CELL-CELL COMMUNICATION ANALYSIS")
    print("="*80)
    
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    try:
        import liana as li
        
        print("Running LIANA for ligand-receptor analysis...")
        
        # Run LIANA
        li.mt.rank_aggregate(
            adata,
            groupby='cell_type',
            expr_prop=0.1,
            use_raw=False,
            verbose=True,
            resource_name='consensus'  # Use consensus database
        )
        
        # Plot interactions
        print("Generating interaction plots...")
        
        # Tumor to immune interactions
        if 'Malignant' in adata.obs['cell_type'].values:
            immune_cells = ['T_cells', 'Cytotoxic_T', 'Regulatory_T', 
                           'Microglia', 'Macrophages', 'NK_cells']
            
            present_immune = [ct for ct in immune_cells if ct in adata.obs['cell_type'].values]
            
            if present_immune:
                li.pl.dotplot(
                    adata,
                    source_labels=['Malignant'],
                    target_labels=present_immune,
                    top_n=20,
                    orderby='magnitude_rank',
                    orderby_ascending=True,
                    size='magnitude',
                    show=False
                )
                plt.savefig(output_dir / 'tumor_immune_interactions.png', 
                           dpi=300, bbox_inches='tight')
                plt.close()
        
        # Aggregate communication by cell type
        li.pl.aggregate(
            adata,
            source_labels=adata.obs['cell_type'].unique(),
            target_labels=adata.obs['cell_type'].unique(),
            how='weight_norm',
            show=False
        )
        plt.savefig(output_dir / 'communication_heatmap.png', 
                   dpi=300, bbox_inches='tight')
        plt.close()
        
        # Export results
        liana_res = adata.uns['liana_res']
        liana_res.to_csv(output_dir / 'liana_results.csv', index=False)
        
        print(f"Communication analysis complete. Results saved to {output_dir}")
        
    except ImportError:
        print("LIANA not installed. Install with: pip install liana-py")
    except Exception as e:
        print(f"Error in communication analysis: {str(e)}")


def metabolic_analysis(adata, output_dir='results/metabolism/'):
    """
    Advanced metabolic profiling focusing on lipid pathways.
    
    Compares metabolic states between cell types and tumor/normal.
    """
    print("\n" + "="*80)
    print("ADVANCED METABOLIC ANALYSIS")
    print("="*80)
    
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Define additional lipid metabolism gene sets
    metabolism_genes = {
        'Glycolysis': ['HK2', 'PFKP', 'PKM', 'LDHA'],
        'OXPHOS': ['COX5B', 'COX6C', 'NDUFA4', 'ATP5A1'],
        'PPP': ['G6PD', 'PGD', 'TKT', 'TALDO1'],
        'Lipogenesis': ['FASN', 'ACACA', 'ACLY', 'SCD'],
        'FA_Uptake': ['CD36', 'FABP4', 'FABP5', 'LPL'],
        'FA_Oxidation': ['CPT1A', 'CPT1C', 'HADHA', 'ACOX1'],
        'Cholesterol_synth': ['HMGCR', 'SQLE', 'FDFT1', 'DHCR24']
    }
    
    # Score pathways
    print("Scoring metabolic pathways...")
    for pathway, genes in metabolism_genes.items():
        available = [g for g in genes if g in adata.var_names]
        if available:
            sc.tl.score_genes(
                adata,
                available,
                score_name=f'{pathway}_score',
                use_raw=False
            )
    
    # Create metabolic profile heatmap
    print("Creating metabolic profile heatmap...")
    
    pathway_scores = [f'{p}_score' for p in metabolism_genes.keys() 
                     if f'{p}_score' in adata.obs.columns]
    
    if pathway_scores:
        # Average by cell type
        metabolic_profile = adata.obs.groupby('cell_type')[pathway_scores].mean()
        
        plt.figure(figsize=(10, 8))
        sns.heatmap(
            metabolic_profile.T,
            cmap='RdBu_r',
            center=0,
            cbar_kws={'label': 'Average Score'},
            xticklabels=True,
            yticklabels=[s.replace('_score', '') for s in pathway_scores]
        )
        plt.title('Metabolic Profile by Cell Type')
        plt.xlabel('Cell Type')
        plt.ylabel('Metabolic Pathway')
        plt.tight_layout()
        plt.savefig(output_dir / 'metabolic_profile_heatmap.png', 
                   dpi=300, bbox_inches='tight')
        plt.close()
        
        # Save numerical data
        metabolic_profile.to_csv(output_dir / 'metabolic_profiles.csv')
    
    # Correlation analysis
    print("Analyzing pathway correlations...")
    
    if len(pathway_scores) > 1 and 'DGAT1' in adata.var_names:
        # Add DGAT1 expression to correlation
        dgat1_expr = adata[:, 'DGAT1'].X.toarray().flatten() if hasattr(adata[:, 'DGAT1'].X, 'toarray') else adata[:, 'DGAT1'].X.flatten()
        
        corr_data = adata.obs[pathway_scores].copy()
        corr_data['DGAT1_expr'] = dgat1_expr
        
        # Compute correlations
        corr_matrix = corr_data.corr()
        
        plt.figure(figsize=(10, 9))
        sns.heatmap(
            corr_matrix,
            cmap='coolwarm',
            center=0,
            annot=True,
            fmt='.2f',
            square=True,
            cbar_kws={'label': 'Correlation'}
        )
        plt.title('Metabolic Pathway Correlations')
        plt.tight_layout()
        plt.savefig(output_dir / 'pathway_correlations.png', 
                   dpi=300, bbox_inches='tight')
        plt.close()
        
        corr_matrix.to_csv(output_dir / 'pathway_correlations.csv')
    
    print(f"Metabolic analysis complete. Results saved to {output_dir}")


def survival_correlation(adata, clinical_data_path=None, output_dir='results/survival/'):
    """
    Correlate molecular features with patient outcomes.
    
    Requires clinical data with survival information.
    """
    print("\n" + "="*80)
    print("SURVIVAL CORRELATION ANALYSIS")
    print("="*80)
    
    if clinical_data_path is None:
        print("No clinical data provided. Skipping survival analysis.")
        print("To run this analysis, provide a CSV with columns:")
        print("  - patient_id: matching adata.obs['patient_id']")
        print("  - survival_time: time in months")
        print("  - event: 1 for death/progression, 0 for censored")
        return
    
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    try:
        from lifelines import KaplanMeierFitter, CoxPHFitter
        
        # Load clinical data
        clinical = pd.read_csv(clinical_data_path)
        
        # Calculate per-patient DGAT1 expression
        if 'DGAT1' in adata.var_names:
            dgat1_by_patient = adata.obs.groupby('patient_id').apply(
                lambda x: adata[x.index, 'DGAT1'].X.mean()
            )
            
            clinical['dgat1_high'] = clinical['patient_id'].map(
                lambda x: dgat1_by_patient.get(x, np.nan) > dgat1_by_patient.median()
            )
            
            # Kaplan-Meier analysis
            kmf = KaplanMeierFitter()
            
            fig, ax = plt.subplots(figsize=(10, 6))
            
            for group in [True, False]:
                mask = clinical['dgat1_high'] == group
                kmf.fit(
                    clinical.loc[mask, 'survival_time'],
                    clinical.loc[mask, 'event'],
                    label=f'DGAT1 {"High" if group else "Low"}'
                )
                kmf.plot_survival_function(ax=ax)
            
            plt.xlabel('Time (months)')
            plt.ylabel('Survival Probability')
            plt.title('Kaplan-Meier Survival by DGAT1 Expression')
            plt.tight_layout()
            plt.savefig(output_dir / 'km_dgat1.png', dpi=300, bbox_inches='tight')
            plt.close()
            
            # Cox regression
            cph = CoxPHFitter()
            cox_data = clinical[['survival_time', 'event']].copy()
            cox_data['dgat1_expr'] = clinical['patient_id'].map(dgat1_by_patient)
            cox_data = cox_data.dropna()
            
            cph.fit(cox_data, duration_col='survival_time', event_col='event')
            
            print("\nCox Proportional Hazards Results:")
            print(cph.summary)
            
            cph.summary.to_csv(output_dir / 'cox_results.csv')
        
        print(f"Survival analysis complete. Results saved to {output_dir}")
        
    except ImportError:
        print("lifelines not installed. Install with: pip install lifelines")
    except Exception as e:
        print(f"Error in survival analysis: {str(e)}")


def main():
    """Run all downstream analyses."""
    import argparse
    
    parser = argparse.ArgumentParser(
        description='Downstream analyses for GBM single-cell data'
    )
    parser.add_argument(
        '--adata',
        type=str,
        required=True,
        help='Path to processed AnnData h5ad file'
    )
    parser.add_argument(
        '--output_dir',
        type=str,
        default='results/downstream/',
        help='Output directory'
    )
    parser.add_argument(
        '--clinical_data',
        type=str,
        default=None,
        help='Path to clinical data CSV (optional)'
    )
    parser.add_argument(
        '--analyses',
        nargs='+',
        default=['trajectory', 'communication', 'metabolism'],
        choices=['trajectory', 'communication', 'metabolism', 'survival'],
        help='Which analyses to run'
    )
    
    args = parser.parse_args()
    
    # Load data
    adata = load_processed_data(args.adata)
    
    output_base = Path(args.output_dir)
    
    # Run requested analyses
    if 'trajectory' in args.analyses:
        trajectory_analysis(adata, output_base / 'trajectory')
    
    if 'communication' in args.analyses:
        cellcell_communication(adata, output_base / 'communication')
    
    if 'metabolism' in args.analyses:
        metabolic_analysis(adata, output_base / 'metabolism')
    
    if 'survival' in args.analyses:
        survival_correlation(adata, args.clinical_data, output_base / 'survival')
    
    print("\n" + "="*80)
    print("ALL DOWNSTREAM ANALYSES COMPLETE!")
    print(f"Results saved to: {output_base}")
    print("="*80)


if __name__ == '__main__':
    main()
