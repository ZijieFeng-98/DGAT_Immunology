# Scripts Overview - ScRNA_GBM_Analysis

## 📁 Complete Script Collection

### 🔬 **Main Analysis Pipelines**

#### 1. **sc_pipeline.py** (Basic Pipeline - TESTED ✓)
**Status**: ✅ **Working and Validated**
- Streamlined 8-step workflow
- Successfully tested on demo data
- Best for: Quick analysis, learning, testing
- Runtime: ~5-10 minutes (demo data)
- Features:
  - Data loading & QC
  - Harmony batch correction
  - Clustering & UMAP
  - CNV inference (optional)
  - Cell-type annotation
  - DGAT1 expression

**Usage:**
```bash
py scripts\sc_pipeline.py \
    --tumour_path data\raw\tumour\ \
    --normal_path data\raw\normal\ \
    --output_dir results\
```

---

#### 2. **sc_pipeline_advanced.py** (Comprehensive Pipeline)
**Status**: ⭐ **Advanced Features**
- Extended analysis capabilities
- Additional downstream analyses
- More sophisticated cell typing
- Best for: Publication-quality analysis
- Runtime: ~30-60 minutes
- Features (All of basic +):
  - Enhanced metabolic pathway analysis
  - Trajectory inference support
  - Cell-cell communication preparation
  - Extended marker sets
  - Advanced QC metrics

**Usage:**
```bash
py scripts\sc_pipeline_advanced.py \
    --tumour_path data\raw\tumour\ \
    --normal_path data\raw\normal\ \
    --output_dir results\advanced\
```

---

### 🎯 **Interactive & Downstream**

#### 3. **step_by_step.py** (Interactive Mode)
**Status**: ✅ **Full Control**
- Run each step individually
- Inspect results between steps
- Adjust parameters on-the-fly
- Save/load checkpoints
- Best for: Exploration, parameter tuning
- Features:
  - Interactive menu
  - Progress statistics
  - Checkpoint management
  - Custom parameter input

**Usage:**
```bash
py scripts\step_by_step.py
```

Then choose from menu:
```
1. Load data
2. Quality control
3. Doublet detection
4. Normalization & integration
...
```

---

#### 4. **downstream_analysis.py** (Advanced Analyses)
**Status**: 🧬 **Extended Features**
- Post-clustering analyses
- Requires completed basic pipeline first
- Best for: Deep biological insights
- Features:
  - Trajectory inference (Monocle-style)
  - Cell-cell communication (LIANA)
  - RNA velocity (if available)
  - Metabolic flux analysis
  - Survival association

**Usage:**
```bash
py scripts\downstream_analysis.py \
    --input results\processed_adata.h5ad \
    --output results\downstream\
```

---

### 🛠️ **Utility Scripts**

#### 5. **create_demo_data.py** (Demo Data Generator)
**Status**: ✅ **Tested**
- Creates synthetic test data
- 10x Genomics format
- Best for: Testing installation
- Output: 500 cells, 2000 genes

**Usage:**
```bash
py scripts\create_demo_data.py
```

---

#### 6. **verify_installation.py** (Setup Verification)
**Status**: ✅ **Ready**
- Checks all dependencies
- Validates environment
- Tests critical functions
- Best for: After installation

**Usage:**
```bash
py scripts\verify_installation.py
```

---

## 🎯 **Which Script Should I Use?**

### For First-Time Users:
1. **Start**: `create_demo_data.py` → Generate test data
2. **Test**: `sc_pipeline.py` → Quick validation
3. **Explore**: `step_by_step.py` → Learn each step
4. **Advance**: `sc_pipeline_advanced.py` → Full analysis

### For Quick Analysis:
→ **sc_pipeline.py** (Basic, tested, fast)

### For Publication:
→ **sc_pipeline_advanced.py** + **downstream_analysis.py**

### For Learning:
→ **step_by_step.py** (Interactive mode)

---

## 📊 **Feature Comparison**

| Feature | Basic | Advanced | Downstream |
|---------|-------|----------|------------|
| Data loading | ✅ | ✅ | - |
| Quality control | ✅ | ✅✅ | - |
| Batch correction | ✅ Harmony | ✅ Harmony+ | - |
| Clustering | ✅ | ✅ | - |
| CNV inference | ✅ | ✅✅ | - |
| Cell annotation | ✅ Basic | ✅✅ Extended | - |
| DGAT1 analysis | ✅ | ✅✅ | ✅✅✅ |
| Trajectory | - | ⚠️ Prep | ✅ |
| Communication | - | ⚠️ Prep | ✅ |
| RNA velocity | - | - | ✅ |
| Metabolic flux | - | ⚠️ Prep | ✅ |
| **Runtime** | 5-10 min | 30-60 min | 20-40 min |

✅ = Included, ✅✅ = Enhanced, ⚠️ = Preparation only

---

## 🚀 **Recommended Workflows**

### Workflow 1: Quick Validation
```bash
# Generate demo data
py scripts\create_demo_data.py

# Run basic pipeline
py scripts\sc_pipeline.py \
    --tumour_path data\raw\demo_tumour\ \
    --normal_path data\raw\demo_normal\ \
    --output_dir results\test\
```

### Workflow 2: Complete Analysis
```bash
# 1. Run advanced pipeline
py scripts\sc_pipeline_advanced.py \
    --tumour_path data\raw\tumour\ \
    --normal_path data\raw\normal\ \
    --output_dir results\complete\

# 2. Run downstream analyses
py scripts\downstream_analysis.py \
    --input results\complete\processed_adata.h5ad \
    --output results\complete\downstream\
```

### Workflow 3: Interactive Exploration
```bash
# Interactive step-by-step
py scripts\step_by_step.py

# Then explore in Jupyter
jupyter notebook
```

---

## 📝 **Script Details**

### File Sizes
- `sc_pipeline.py`: 17 KB (443 lines)
- `sc_pipeline_advanced.py`: 33 KB (~900 lines)
- `downstream_analysis.py`: 15 KB (~400 lines)
- `step_by_step.py`: 21 KB (660 lines)
- `create_demo_data.py`: 3 KB (64 lines)
- `verify_installation.py`: 7 KB (~200 lines)

### Dependencies
All scripts require:
- scanpy, anndata, numpy, pandas
- harmonypy (batch correction)
- infercnvpy (CNV analysis)
- matplotlib, seaborn (plotting)

Advanced/downstream additionally need:
- liana (cell-cell communication)
- scvelo (RNA velocity) - optional
- Additional packages listed in requirements.txt

---

## 🔧 **Customization**

All scripts support:
- Command-line arguments
- Configuration files (config.yaml)
- Custom marker dictionaries
- Adjustable QC thresholds
- Flexible output locations

---

## 📚 **Documentation**

Each script includes:
- Comprehensive docstrings
- Parameter descriptions
- Example usage
- References to methodology

For methodology details, see:
- `../Protocols/scRNA_Klemm_Pipeline_Guide.md`

---

## ✅ **Tested Scripts**

- ✅ **sc_pipeline.py** - Fully tested on demo data
- ✅ **create_demo_data.py** - Tested and working
- ⏳ **sc_pipeline_advanced.py** - Ready to test
- ⏳ **downstream_analysis.py** - Ready to test
- ✅ **step_by_step.py** - Structure tested
- ⏳ **verify_installation.py** - Ready to test

---

## 🎓 **Learning Path**

1. **Beginner**: `create_demo_data.py` → `sc_pipeline.py`
2. **Intermediate**: `step_by_step.py` (interactive)
3. **Advanced**: `sc_pipeline_advanced.py`
4. **Expert**: `downstream_analysis.py`

---

**All scripts are now organized and ready to use!**  
**Start with `sc_pipeline.py` for quick results or `step_by_step.py` for learning!**

