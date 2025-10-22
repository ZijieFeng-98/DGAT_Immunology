#!/bin/bash
# Automated GBM Single-Cell Analysis Pipeline
# This script runs the complete analysis workflow

set -e  # Exit on error

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

# Function to print colored messages
print_info() {
    echo -e "${GREEN}[INFO]${NC} $1"
}

print_warning() {
    echo -e "${YELLOW}[WARNING]${NC} $1"
}

print_error() {
    echo -e "${RED}[ERROR]${NC} $1"
}

# Check if required arguments are provided
if [ "$#" -lt 3 ]; then
    echo "Usage: $0 <tumor_path1,tumor_path2,...> <normal_path1,normal_path2,...> <output_dir> [options]"
    echo ""
    echo "Example:"
    echo "  $0 data/tumor1,data/tumor2 data/normal1 results/"
    echo ""
    echo "Optional arguments:"
    echo "  --min-genes <int>       Minimum genes per cell (default: 200)"
    echo "  --max-genes <int>       Maximum genes per cell (default: 6000)"
    echo "  --max-mito <float>      Maximum mitochondrial % (default: 15.0)"
    echo "  --resolution <float>    Clustering resolution (default: 0.5)"
    echo "  --downstream            Run downstream analyses"
    echo "  --skip-main             Skip main pipeline (only downstream)"
    echo "  --clinical-data <path>  Path to clinical data CSV for survival analysis"
    exit 1
fi

# Parse arguments
TUMOR_PATHS=$1
NORMAL_PATHS=$2
OUTPUT_DIR=$3
shift 3

# Default parameters
MIN_GENES=200
MAX_GENES=6000
MAX_MITO=15.0
RESOLUTION=0.5
RUN_DOWNSTREAM=false
SKIP_MAIN=false
CLINICAL_DATA=""

# Parse optional arguments
while [[ $# -gt 0 ]]; do
    case $1 in
        --min-genes)
            MIN_GENES="$2"
            shift 2
            ;;
        --max-genes)
            MAX_GENES="$2"
            shift 2
            ;;
        --max-mito)
            MAX_MITO="$2"
            shift 2
            ;;
        --resolution)
            RESOLUTION="$2"
            shift 2
            ;;
        --downstream)
            RUN_DOWNSTREAM=true
            shift
            ;;
        --skip-main)
            SKIP_MAIN=true
            shift
            ;;
        --clinical-data)
            CLINICAL_DATA="$2"
            shift 2
            ;;
        *)
            print_warning "Unknown option: $1"
            shift
            ;;
    esac
done

# Convert comma-separated paths to space-separated for Python
TUMOR_PATHS_ARRAY=(${TUMOR_PATHS//,/ })
NORMAL_PATHS_ARRAY=(${NORMAL_PATHS//,/ })

print_info "Starting GBM Single-Cell Analysis Pipeline"
echo "=========================================="
echo "Tumor samples: ${TUMOR_PATHS_ARRAY[@]}"
echo "Normal samples: ${NORMAL_PATHS_ARRAY[@]}"
echo "Output directory: $OUTPUT_DIR"
echo "Parameters:"
echo "  - Min genes: $MIN_GENES"
echo "  - Max genes: $MAX_GENES"
echo "  - Max mito: $MAX_MITO%"
echo "  - Resolution: $RESOLUTION"
echo "=========================================="
echo ""

# Check if Python script exists
if [ ! -f "sc_glioblastoma_dgat1_pipeline.py" ]; then
    print_error "Main pipeline script not found!"
    print_error "Make sure you're in the correct directory."
    exit 1
fi

# Create output directory
mkdir -p "$OUTPUT_DIR"

# Run main pipeline
if [ "$SKIP_MAIN" = false ]; then
    print_info "Step 1/2: Running main analysis pipeline..."
    echo ""
    
    python sc_glioblastoma_dgat1_pipeline.py \
        --tumor_paths ${TUMOR_PATHS_ARRAY[@]} \
        --normal_paths ${NORMAL_PATHS_ARRAY[@]} \
        --output_dir "$OUTPUT_DIR" \
        --min_genes $MIN_GENES \
        --max_genes $MAX_GENES \
        --max_mito $MAX_MITO \
        --resolution $RESOLUTION
    
    if [ $? -eq 0 ]; then
        print_info "Main pipeline completed successfully!"
        echo ""
    else
        print_error "Main pipeline failed!"
        exit 1
    fi
else
    print_info "Skipping main pipeline (--skip-main specified)"
    echo ""
fi

# Run downstream analyses
if [ "$RUN_DOWNSTREAM" = true ]; then
    print_info "Step 2/2: Running downstream analyses..."
    echo ""
    
    if [ ! -f "downstream_analysis.py" ]; then
        print_error "Downstream analysis script not found!"
        exit 1
    fi
    
    ADATA_PATH="$OUTPUT_DIR/processed_adata.h5ad"
    
    if [ ! -f "$ADATA_PATH" ]; then
        print_error "Processed data file not found: $ADATA_PATH"
        print_error "Make sure the main pipeline completed successfully."
        exit 1
    fi
    
    DOWNSTREAM_CMD="python downstream_analysis.py \
        --adata $ADATA_PATH \
        --output_dir $OUTPUT_DIR/downstream"
    
    # Add clinical data if provided
    if [ -n "$CLINICAL_DATA" ]; then
        if [ -f "$CLINICAL_DATA" ]; then
            DOWNSTREAM_CMD="$DOWNSTREAM_CMD --clinical_data $CLINICAL_DATA"
            print_info "Including survival analysis with clinical data: $CLINICAL_DATA"
        else
            print_warning "Clinical data file not found: $CLINICAL_DATA"
            print_warning "Continuing without survival analysis..."
        fi
    fi
    
    eval $DOWNSTREAM_CMD
    
    if [ $? -eq 0 ]; then
        print_info "Downstream analyses completed successfully!"
        echo ""
    else
        print_error "Downstream analyses failed!"
        exit 1
    fi
else
    print_info "Skipping downstream analyses (use --downstream to enable)"
    echo ""
fi

# Summary
echo ""
echo "=========================================="
print_info "ANALYSIS COMPLETE!"
echo "=========================================="
echo ""
echo "Results saved to: $OUTPUT_DIR"
echo ""
echo "Key output files:"
echo "  - processed_adata.h5ad           (Main analyzed dataset)"
echo "  - umap_overview.png              (Visual overview)"
echo "  - dgat1_expression.png           (DGAT1 patterns)"
echo "  - cell_type_annotation.png       (Cell types)"
echo "  - de_tumor_vs_normal.csv         (Differential expression)"
echo ""

if [ "$RUN_DOWNSTREAM" = true ]; then
    echo "Downstream analysis results:"
    echo "  - downstream/trajectory/         (Pseudotime analysis)"
    echo "  - downstream/communication/      (Cell-cell interactions)"
    echo "  - downstream/metabolism/         (Metabolic profiling)"
    if [ -n "$CLINICAL_DATA" ] && [ -f "$CLINICAL_DATA" ]; then
        echo "  - downstream/survival/           (Survival correlations)"
    fi
    echo ""
fi

echo "Next steps:"
echo "  1. Review the output plots in $OUTPUT_DIR"
echo "  2. Load processed_adata.h5ad in Python for interactive analysis"
echo "  3. Check the CSV files for detailed statistics"
echo ""

if [ "$RUN_DOWNSTREAM" = false ]; then
    echo "To run downstream analyses later:"
    echo "  python downstream_analysis.py \\"
    echo "    --adata $OUTPUT_DIR/processed_adata.h5ad \\"
    echo "    --output_dir $OUTPUT_DIR/downstream"
    echo ""
fi

print_info "All done! Happy analyzing! 🎉"
