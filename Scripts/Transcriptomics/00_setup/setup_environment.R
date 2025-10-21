# =============================================================================
# 00_setup_env_simple.R — Simplified environment setup for DGAT Immunology
# =============================================================================
# 
# This script sets up the R environment with essential packages only,
# avoiding compilation issues on macOS.
#
# Usage: Rscript Scripts/00_setup_env_simple.R
# =============================================================================

cat("🚀 Setting up DGAT Immunology analysis environment (simplified)...\n")

# Initialize renv if not already done
if (!requireNamespace("renv", quietly = TRUE)) {
  cat("📦 Installing renv...\n")
  install.packages("renv", type = "binary")
}

if (!dir.exists("renv")) {
  cat("🔧 Initializing renv environment...\n")
  renv::init(bare = TRUE)
}

# Essential CRAN packages (install as binaries to avoid compilation)
essential_pkgs <- c(
  "data.table",      # Fast data manipulation
  "dplyr",           # Data manipulation  
  "ggplot2",         # Plotting
  "reshape2",        # Data reshaping
  "pheatmap",        # Heatmaps
  "pbapply",         # Progress bars
  "survival",        # Survival analysis
  "survminer",       # Survival plots
  "scales",          # Scale functions
  "jsonlite",        # JSON handling
  "remotes",         # Install from GitHub
  "BiocManager"      # Bioconductor package manager
)

cat("📥 Installing essential CRAN packages (as binaries)...\n")
for (pkg in essential_pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    cat("  Installing", pkg, "...\n")
    tryCatch({
      install.packages(pkg, type = "binary")
      cat("  ✅ Installed", pkg, "\n")
    }, error = function(e) {
      cat("  ⚠️  Failed to install", pkg, "- trying source...\n")
      tryCatch({
        install.packages(pkg)
        cat("  ✅ Installed", pkg, "from source\n")
      }, error = function(e2) {
        cat("  ❌ Failed to install", pkg, "- skipping\n")
      })
    })
  } else {
    cat("  ✅", pkg, "already installed\n")
  }
}

# Try to install Seurat (essential for single-cell)
cat("📥 Installing Seurat...\n")
if (!requireNamespace("Seurat", quietly = TRUE)) {
  tryCatch({
    install.packages("Seurat", type = "binary")
    cat("✅ Installed Seurat\n")
  }, error = function(e) {
    cat("⚠️  Failed to install Seurat from binary - trying source...\n")
    tryCatch({
      install.packages("Seurat")
      cat("✅ Installed Seurat from source\n")
    }, error = function(e2) {
      cat("❌ Failed to install Seurat - will need to install manually\n")
    })
  })
} else {
  cat("✅ Seurat already installed\n")
}

# Bioconductor packages
cat("📥 Installing Bioconductor packages...\n")

# Install BiocManager first
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager", type = "binary")
}

# Essential Bioconductor packages
bioc_pkgs <- c("TCGAbiolinks", "GSVA", "msigdbr")

for (pkg in bioc_pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    cat("  Installing", pkg, "from Bioconductor...\n")
    tryCatch({
      BiocManager::install(pkg, update = FALSE, ask = FALSE)
      cat("  ✅ Installed", pkg, "\n")
    }, error = function(e) {
      cat("  ⚠️  Failed to install", pkg, "- skipping\n")
    })
  } else {
    cat("  ✅", pkg, "already installed\n")
  }
}

# Try to install IOBR from GitHub
cat("📥 Installing IOBR from GitHub...\n")
if (!requireNamespace("IOBR", quietly = TRUE)) {
  tryCatch({
    remotes::install_github("IOBR/IOBR")
    cat("✅ Installed IOBR\n")
  }, error = function(e) {
    cat("⚠️  Failed to install IOBR - will need to install manually\n")
    cat("   Try: remotes::install_github('IOBR/IOBR')\n")
  })
} else {
  cat("✅ IOBR already installed\n")
}

# Take snapshot of current environment
cat("📸 Taking snapshot of package versions...\n")
tryCatch({
  renv::snapshot()
}, error = function(e) {
  cat("⚠️  Could not take snapshot - continuing anyway\n")
})

cat("\n🎉 Environment setup complete!\n")
cat("📁 Project structure:\n")
cat("  Raw_Data/     - Untouched downloads\n")
cat("  Processed_Data/ - Normalized matrices\n")
cat("  Results/      - Final outputs\n")
cat("  Scripts/      - Analysis code\n")
cat("\n💡 Next steps:\n")
cat("  1. Test package loading with: library(data.table); library(dplyr)\n")
cat("  2. Download data using Scripts/10_download_tcga.R\n")
cat("  3. Run analysis using Scripts/20_bulk_analysis.R\n")
cat("\n⚠️  Note: Some packages may need manual installation if compilation failed\n")
cat("    Check with: installed.packages()[c('Seurat', 'IOBR', 'TCGAbiolinks')]\n")
cat("\n")
