#!/usr/bin/env Rscript
# Scripts/02_data_download/01_get_tcga.R
# Non-interactive TCGA download & prep (GBM + LGG)
options(stringsAsFactors = FALSE)
options(bitmapType = "cairo")
options(menu.graphics = FALSE)

req <- c("TCGAbiolinks","SummarizedExperiment","dplyr","data.table","stringr","tools","BiocParallel","AnnotationDbi","org.Hs.eg.db")
for(p in req) if(!requireNamespace(p, quietly=TRUE)) suppressWarnings(install.packages(p, repos="https://cloud.r-project.org"))
suppressPackageStartupMessages({
  library(TCGAbiolinks)
  library(SummarizedExperiment)
  library(dplyr); library(data.table); library(stringr)
  library(BiocParallel); library(AnnotationDbi); library(org.Hs.eg.db)
})

dir.create("Raw_Data/TCGA", showWarnings = FALSE, recursive = TRUE)
dir.create("Processed_Data/Bulk", showWarnings = FALSE, recursive = TRUE)
dir.create("Results/Reports", showWarnings = FALSE, recursive = TRUE)

cohorts <- c("TCGA-GBM","TCGA-LGG")
profile  <- "HTSeq - FPKM"  # human-friendly scale for quick EDA; swap to counts later if needed

fetch_cohort <- function(project){
  message("\n=== TCGA: ", project, " ===")
  q <- GDCquery(project = project,
                data.category = "Transcriptome Profiling",
                data.type = "Gene Expression Quantification",
                workflow.type = profile)
  GDCdownload(q, method = "api", files.per.chunk = 50)
  se <- GDCprepare(q)
  se
}

harmonize_symbols <- function(se){
  # Gene IDs arrive as Ensembl IDs (with/without version). Map to HGNC symbols.
  expr <- as.data.frame(SummarizedExperiment::assay(se))
  genes <- rownames(expr)
  genes <- sub("\\..*$","",genes) # strip version
  # Map Ensembl -> SYMBOL via org.Hs.eg.db
  suppressWarnings({
    map <- AnnotationDbi::select(org.Hs.eg.db, keys = genes, keytype = "ENSEMBL", columns = c("SYMBOL"))
  })
  map <- distinct(map, ENSEMBL, .keep_all = TRUE)
  expr$ENSEMBL <- genes
  M <- merge(expr, map, by.x = "ENSEMBL", by.y = "ENSEMBL", all.x = TRUE)
  # keep rows with a SYMBOL; collapse duplicates by mean
  M <- M[!is.na(M$SYMBOL) & M$SYMBOL != "", ]
  sym <- M$SYMBOL; M$ENSEMBL <- NULL; M$SYMBOL <- NULL
  E <- as.data.frame(M); rownames(E) <- sym
  # collapse duplicate symbols
  E <- as.data.table(E)[, lapply(.SD, mean, na.rm = TRUE), by = rownames(E)]
  mat <- as.matrix(E[,-1]); rownames(mat) <- E[[1]]
  mat
}

save_objects <- function(project, se, mat){
  raw_rds   <- file.path("Raw_Data/TCGA", paste0(gsub("-","_",project), "_SE_", gsub(" ","_",profile), ".rds"))
  expr_rds  <- file.path("Processed_Data/Bulk", paste0(gsub("-","_",project), "_Expr_HGNC_FPKM.rds"))
  clin_csv  <- file.path("Processed_Data/Bulk", paste0(gsub("-","_",project), "_Clinical.csv"))
  saveRDS(se, raw_rds)
  saveRDS(mat, expr_rds)
  clin <- as.data.frame(colData(se))
  data.table::fwrite(clin, clin_csv)
  message("Saved: ", raw_rds)
  message("Saved: ", expr_rds)
  message("Saved: ", clin_csv)
}

log_summary <- function(project, mat){
  keep <- c("DGAT1","DGAT2","CD8A","GZMB","FOXP3","MRC1","CSF1","TGFB1","CCL2")
  hit  <- intersect(keep, rownames(mat))
  msg  <- paste0(project, ": ", nrow(mat), " genes × ", ncol(mat), " samples. Found markers: ", paste(hit, collapse=", "))
  writeLines(msg, con = file.path("Results/Reports", "01_get_tcga.log"), sep = "\n", useBytes = TRUE)
  message(msg)
}

for (proj in cohorts){
  se <- fetch_cohort(proj)
  mat <- harmonize_symbols(se)
  save_objects(proj, se, mat)
  log_summary(proj, mat)
}

message("\nTCGA download COMPLETE.")
