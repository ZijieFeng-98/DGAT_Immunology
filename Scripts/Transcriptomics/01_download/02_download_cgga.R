#!/usr/bin/env Rscript
# Scripts/02_data_download/02_get_cgga.R
options(stringsAsFactors = FALSE); options(menu.graphics = FALSE)

req <- c("GEOquery","data.table","dplyr","stringr")
for(p in req) if(!requireNamespace(p, quietly=TRUE)) install.packages(p, repos="https://cloud.r-project.org")
suppressPackageStartupMessages({ library(GEOquery); library(data.table); library(dplyr); library(stringr) })

dir.create("Raw_Data/CGGA", showWarnings = FALSE, recursive = TRUE)
dir.create("Processed_Data/Bulk", showWarnings = FALSE, recursive = TRUE)
dir.create("Results/Reports", showWarnings = FALSE, recursive = TRUE)

# Try common CGGA mirrors on GEO (this varies by cohort); otherwise instruct manual.
# You can add/adjust accession IDs here if you know the exact CGGA dataset you want.
candidate_geo <- c("GSE16011","GSE83300","GSE86946")  # placeholders; at least one should fetch GBM-like matrices

found <- FALSE
for (acc in candidate_geo){
  message("\nTrying GEO accession: ", acc)
  ok <- try({
    gset <- getGEO(acc, GSEMatrix = TRUE, AnnotGPL = TRUE)
    if (length(gset) > 0){
      eset <- gset[[1]]
      expr <- Biobase::exprs(eset)
      fdat <- Biobase::fData(eset)
      pdat <- Biobase::pData(eset)
      # map gene symbols if present
      if (!"Gene.symbol" %in% colnames(fdat)) {
        # fallbacks (platform dependent)
        symcol <- grep("symbol|SYMBOL|Gene", colnames(fdat), ignore.case = TRUE, value = TRUE)[1]
      } else symcol <- "Gene.symbol"
      if (!is.na(symcol)){
        sym <- fdat[[symcol]]
        keep <- !is.na(sym) & sym!=""
        expr <- expr[keep, , drop=FALSE]
        sym  <- sym[keep]
        # collapse duplicates by mean
        DT <- as.data.table(expr); DT[, SYMBOL := sym]
        agg <- DT[, lapply(.SD, mean, na.rm=TRUE), by=SYMBOL]
        mat <- as.matrix(agg[,-1]); rownames(mat) <- agg$SYMBOL
      } else {
        mat <- expr
      }
      # Save
      expr_rds <- file.path("Processed_Data/Bulk", paste0("CGGA_like_", acc, "_Expr_HGNC.rds"))
      clin_csv  <- file.path("Processed_Data/Bulk", paste0("CGGA_like_", acc, "_Clinical.csv"))
      saveRDS(mat, expr_rds)
      fwrite(pdat, clin_csv)
      writeLines(paste0("CGGA/GEO success: ", acc, "  (", nrow(mat), " genes x ", ncol(mat), " samples)"),
                 con = file.path("Results/Reports", "02_get_cgga.log"))
      message("Saved: ", expr_rds)
      message("Saved: ", clin_csv)
      found <- TRUE
      break
    }
  }, silent = TRUE)
}

if (!found){
  msg <- paste(
    "Auto-download via GEO failed to locate a CGGA mirror.",
    "Please download a CGGA GBM expression matrix + clinical file from https://www.cgga.org.cn/ ",
    "and place them here:",
    "  Raw_Data/CGGA/CGGA_GBM_expression.tsv",
    "  Raw_Data/CGGA/CGGA_GBM_clinical.tsv",
    "Then re-run this script to parse and save harmonized RDS/CSV.",
    sep = "\n")
  writeLines(msg, con = file.path("Results/Reports", "02_get_cgga.log"))
  message(msg)
  # Parse if user already dropped files:
  expr_fp <- "Raw_Data/CGGA/CGGA_GBM_expression.tsv"
  clin_fp <- "Raw_Data/CGGA/CGGA_GBM_clinical.tsv"
  if (file.exists(expr_fp) && file.exists(clin_fp)){
    expr <- fread(expr_fp)
    # expect first column gene symbol, rest samples
    gs_col <- 1
    sym <- expr[[gs_col]]
    mat <- as.matrix(expr[,-gs_col]); rownames(mat) <- sym
    clin <- fread(clin_fp)
    saveRDS(mat, file.path("Processed_Data/Bulk","CGGA_GBM_Expr_HGNC.rds"))
    fwrite(clin, file.path("Processed_Data/Bulk","CGGA_GBM_Clinical.csv"))
    message("Parsed dropped CGGA files and saved harmonized objects.")
  }
}
message("\nCGGA step COMPLETE (either via GEO-like set or manual drop-in).")
