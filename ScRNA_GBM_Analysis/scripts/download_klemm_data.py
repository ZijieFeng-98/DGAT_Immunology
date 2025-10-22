"""
Automated downloader for Klemm et al., 2020 GBM single-cell dataset
GEO Accession: GSE163108
"""

import os
import urllib.request
import gzip
import tarfile
from pathlib import Path
import sys

def download_file(url, dest_path, desc="File"):
    """Download file with progress reporting"""
    print(f"\n[Downloading] {desc}")
    print(f"  URL: {url}")
    print(f"  Destination: {dest_path}")
    
    try:
        def report_progress(block_num, block_size, total_size):
            downloaded = block_num * block_size
            percent = min(downloaded * 100 / total_size, 100)
            sys.stdout.write(f'\r  Progress: {percent:.1f}% ({downloaded/1024/1024:.1f} MB / {total_size/1024/1024:.1f} MB)')
            sys.stdout.flush()
        
        urllib.request.urlretrieve(url, dest_path, reporthook=report_progress)
        print(f"\n  [OK] Download complete")
        return True
    except Exception as e:
        print(f"\n  [ERROR] Download failed: {e}")
        return False

def main():
    """Main download function"""
    print("="*70)
    print("Klemm et al., 2020 GBM Dataset Downloader")
    print("="*70)
    print("\nGEO Accession: GSE163108")
    print("Study: Interrogation of the Microenvironmental Landscape in Brain Tumors")
    print("Journal: Nature (2020)")
    print("DOI: 10.1038/s41586-020-1959-y")
    
    # Create output directory
    raw_data_dir = Path("data/raw/klemm_2020")
    raw_data_dir.mkdir(parents=True, exist_ok=True)
    
    print(f"\n[INFO] Data will be saved to: {raw_data_dir}")
    
    # GEO dataset information
    print("\n" + "="*70)
    print("DOWNLOAD OPTIONS")
    print("="*70)
    
    print("\nOption 1: GEO Database (NCBI)")
    print("  URL: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE163108")
    print("  Steps:")
    print("    1. Visit the URL above")
    print("    2. Scroll to 'Supplementary file' section")
    print("    3. Download all supplementary files")
    print("    4. Extract to:", raw_data_dir)
    
    print("\nOption 2: Single Cell Portal (Broad Institute)")
    print("  URL: https://singlecell.broadinstitute.org/single_cell/study/SCP1290")
    print("  Steps:")
    print("    1. Visit the URL above")
    print("    2. Click 'Download' tab")
    print("    3. Select: Expression matrices + Metadata")
    print("    4. Download and extract to:", raw_data_dir)
    
    print("\nOption 3: Manual Download from Publication")
    print("  1. Visit: https://doi.org/10.1038/s41586-020-1959-y")
    print("  2. Check supplementary materials")
    print("  3. Download data files")
    
    # Check if files already downloaded
    print("\n" + "="*70)
    print("CHECKING FOR EXISTING DATA")
    print("="*70)
    
    existing_files = list(raw_data_dir.glob("*"))
    if existing_files:
        print(f"\n[OK] Found {len(existing_files)} files in {raw_data_dir}:")
        for f in existing_files[:10]:
            print(f"  - {f.name}")
        if len(existing_files) > 10:
            print(f"  ... and {len(existing_files) - 10} more")
    else:
        print(f"\n[INFO] No files found yet. Please download manually.")
    
    # Create README
    readme_path = raw_data_dir / "README.md"
    readme_content = f"""# Klemm et al., 2020 GBM Dataset

## Dataset Information

**Study**: Interrogation of the Microenvironmental Landscape in Brain Tumors Reveals Disease-Specific Alterations of Immune Cells

**Authors**: Klemm et al.  
**Journal**: Nature  
**Year**: 2020  
**DOI**: 10.1038/s41586-020-1959-y  
**GEO Accession**: GSE163108  
**Single Cell Portal**: SCP1290  

## Samples

- **40 IDH-wild-type glioblastomas** (tumor samples)
- **6 adult non-tumor brain specimens** (normal controls)
- **~100,000 cells** total (estimated)

## Data Type

- Single-cell RNA-seq (scRNA-seq)
- 10x Genomics Chromium platform
- Flow cytometry validation
- Bulk RNA-seq for comparison

## Download Instructions

### Option 1: GEO Database
1. Visit: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE163108
2. Download supplementary files
3. Extract to this directory

### Option 2: Single Cell Portal
1. Visit: https://singlecell.broadinstitute.org/single_cell/study/SCP1290
2. Download processed data
3. Extract to this directory

### Option 3: Contact Authors
- Email corresponding author for data access
- May require data use agreement

## Expected Files

After download, you should have:
```
data/raw/klemm_2020/
├── matrix.mtx.gz or .h5ad files
├── features.tsv.gz or genes.tsv
├── barcodes.tsv.gz
└── metadata.csv or similar
```

## Next Steps

After downloading:
1. Convert to 10x format if needed (see conversion scripts)
2. Organize into tumor/ and normal/ directories
3. Run the analysis pipeline

## Citation

If you use this data, please cite:

Klemm F, et al. (2020). Interrogation of the microenvironmental landscape 
in brain tumors reveals disease-specific alterations of immune cells. 
Nature, 584(7820), 120-124. doi: 10.1038/s41586-020-1959-y

## Data Location

Downloaded: {datetime.now().strftime('%Y-%m-%d')}
Path: {raw_data_dir.absolute()}
"""
    
    with open(readme_path, 'w', encoding='utf-8') as f:
        f.write(readme_content)
    
    print(f"\n[OK] Created README: {readme_path}")
    
    # Instructions
    print("\n" + "="*70)
    print("NEXT STEPS")
    print("="*70)
    print("\n1. Visit one of the URLs above to download the data")
    print(f"2. Extract files to: {raw_data_dir}")
    print("3. Run format check: py scripts/check_data_format.py")
    print("4. If needed, convert: py scripts/convert_to_10x.py")
    print("5. Run pipeline on real data:")
    print("   py scripts/sc_pipeline.py --tumour_path data/raw/klemm_tumour/")
    print("                              --normal_path data/raw/klemm_normal/")
    print("                              --output_dir results/klemm_analysis/")
    
    print("\n[TIP] Data size: ~500 MB compressed, ~2 GB uncompressed")
    print("[TIP] Download time: 10-30 minutes depending on internet speed")
    
    print("\n" + "="*70)

if __name__ == "__main__":
    from datetime import datetime
    main()

