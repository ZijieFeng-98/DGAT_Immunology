"""Quick script to examine analysis results"""
import scanpy as sc
import sys

if len(sys.argv) > 1:
    input_file = sys.argv[1]
else:
    input_file = 'results/real_pilot/processed_adata.h5ad'

print("="*70)
print("REAL GBM DATA ANALYSIS RESULTS")
print("="*70)

adata = sc.read_h5ad(input_file)

print(f"\n[DATA] Dataset Summary:")
print(f"  - Total cells: {adata.n_obs:,}")
print(f"  - Total genes: {adata.n_vars:,}")

print(f"\n[SAMPLE] Sample Composition:")
for sample, count in adata.obs['sample'].value_counts().items():
    print(f"  - {sample}: {count:,} cells ({count/adata.n_obs*100:.1f}%)")

if 'leiden' in adata.obs.columns:
    n_clusters = len(adata.obs['leiden'].cat.categories)
    print(f"\n[CLUSTER] Clustering Results:")
    print(f"  - Number of clusters: {n_clusters}")
    for cluster, count in adata.obs['leiden'].value_counts().sort_index().items():
        print(f"    Cluster {cluster}: {count:,} cells ({count/adata.n_obs*100:.1f}%)")

if 'cell_type' in adata.obs.columns:
    n_types = adata.obs['cell_type'].nunique()
    print(f"\n[CELL TYPE] Cell Type Annotation:")
    print(f"  - Cell types identified: {n_types}")
    for cell_type, count in adata.obs['cell_type'].value_counts().items():
        print(f"    {cell_type}: {count:,} cells ({count/adata.n_obs*100:.1f}%)")

# Check for DGAT1 and lipid genes
lipid_genes = ['DGAT1', 'DGAT2', 'FASN', 'ACLY', 'FABP4', 'FABP5', 'APOE', 'CD68', 'CD3D']
available = [g for g in lipid_genes if g in adata.var_names]

print(f"\n[GENE] Key Genes Available:")
print(f"  - DGAT1: {'DGAT1' in adata.var_names}")
print(f"  - Available genes: {', '.join(available) if available else 'None from list'}")
print(f"  - Total genes: {adata.n_vars:,}")

print(f"\n[OK] Analysis file loaded successfully!")
print(f"[FILE] Location: {input_file}")

print("\n" + "="*70)

