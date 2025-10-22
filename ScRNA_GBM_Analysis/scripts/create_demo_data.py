"""Generate small demo dataset for testing"""
import numpy as np
import pandas as pd
import scipy.sparse as sp
from pathlib import Path
import gzip

def create_demo_10x_data(output_dir, n_genes=2000, n_cells=500, prefix="demo"):
    """Create a small demo 10x format dataset"""
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)
    
    # Generate random count matrix
    np.random.seed(42)
    density = 0.1  # 10% non-zero
    data = np.random.negative_binomial(5, 0.3, size=int(n_genes * n_cells * density))
    
    # Create sparse matrix
    row_indices = np.random.randint(0, n_genes, size=len(data))
    col_indices = np.random.randint(0, n_cells, size=len(data))
    matrix = sp.coo_matrix((data, (row_indices, col_indices)), shape=(n_genes, n_cells))
    
    # Save matrix.mtx.gz
    with gzip.open(output_path / "matrix.mtx.gz", 'wt') as f:
        f.write(f"%%MatrixMarket matrix coordinate integer general\n")
        f.write(f"%\n")
        f.write(f"{n_genes} {n_cells} {len(matrix.data)}\n")
        for i, j, v in zip(matrix.row + 1, matrix.col + 1, matrix.data):
            f.write(f"{i} {j} {int(v)}\n")
    
    # Generate gene names
    gene_names = [f"GENE{i:04d}" for i in range(n_genes)]
    # Add some realistic gene names
    real_genes = ['DGAT1', 'DGAT2', 'FASN', 'CD68', 'CD3D', 'GFAP', 'MBP', 'PECAM1', 
                  'MT-CO1', 'MT-CO2', 'MT-ND1', 'PTPRC', 'CSF1R', 'ITGAM']
    for idx, gene in enumerate(real_genes[:min(len(real_genes), n_genes)]):
        gene_names[idx] = gene
    
    # Save features.tsv.gz
    with gzip.open(output_path / "features.tsv.gz", 'wt') as f:
        for gene in gene_names:
            f.write(f"{gene}\t{gene}\tGene Expression\n")
    
    # Save barcodes.tsv.gz
    with gzip.open(output_path / "barcodes.tsv.gz", 'wt') as f:
        for i in range(n_cells):
            f.write(f"{prefix}-{i:04d}-1\n")
    
    print(f"[OK] Created demo dataset: {output_dir}")
    print(f"  - {n_genes} genes")
    print(f"  - {n_cells} cells")

if __name__ == "__main__":
    # Create demo tumour data
    create_demo_10x_data("data/raw/demo_tumour", n_genes=2000, n_cells=300, prefix="TUMOUR")
    
    # Create demo normal data
    create_demo_10x_data("data/raw/demo_normal", n_genes=2000, n_cells=200, prefix="NORMAL")
    
    print("\n[OK] Demo datasets created!")
    print("\nTo use demo data, run:")
    print("  py scripts/sc_pipeline.py --tumour_path data/raw/demo_tumour/ --normal_path data/raw/demo_normal/ --output_dir results/demo/")

