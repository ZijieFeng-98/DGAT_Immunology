#!/usr/bin/env Rscript
# =====================================================================
# CGGA setup & verification (fixed & hardened)
# - Creates folders for mRNAseq_693 and mRNAseq_325
# - Gives clear manual download instructions (no fake automation)
# - When files exist, verifies format, DGAT1/2 presence (critical),
#   reports optional "bonus" panel, and runs QC checks
# - Writes standardized copies for downstream scripts
# =====================================================================

suppressPackageStartupMessages({
  library(data.table)
})

ROOT <- "Data/CGGA"
COHORTS <- list(
  mRNAseq_693 = list(
    dir_raw  = file.path(ROOT, "mRNAseq_693", "Raw"),
    dir_std  = file.path(ROOT, "mRNAseq_693", "Standard"),
    # Flexible, case-insensitive; supports txt/tsv/csv/xls/xlsx
    expr_pat = "CGGA.*693.*RSEM.*\\.(txt|tsv|csv|xls|xlsx)$",
    clin_pat = "CGGA.*693.*clinical.*\\.(txt|tsv|csv|xls|xlsx)$"
  ),
  mRNAseq_325 = list(
    dir_raw  = file.path(ROOT, "mRNAseq_325", "Raw"),
    dir_std  = file.path(ROOT, "mRNAseq_325", "Standard"),
    expr_pat = "CGGA.*325.*RSEM.*symbol.*\\.(txt|tsv|csv|xls|xlsx)$",
    clin_pat = "CGGA.*325.*clinical.*\\.(txt|tsv|csv|xls|xlsx)$"
  )
)

# ---- Panels: Critical (hard fail) vs Bonus (informative) ----
PROTECTED_CRITICAL <- c("DGAT1", "DGAT2")
PROTECTED_BONUS    <- c("SOAT1","SOAT2","CPT1A","CPT1B","PLIN2","PLIN3",
                        "CD8A","GZMB","IFNG","MRC1","CD163","ARG1","S100A8","S100A9")

say <- function(...) cat(sprintf(...), "\n")
mk  <- function(p) { dir.create(p, recursive = TRUE, showWarnings = FALSE); p }

exists_nonempty <- function(p) file.exists(p) && isTRUE(file.info(p)$size > 0)

# --- Safe readers with clear error messages ---
read_any_table <- function(path) {
  ext <- tolower(tools::file_ext(path))
  tryCatch({
    if (ext %in% c("tsv","txt")) {
      # try tab then comma
      dt <- try(fread(path, sep = "\t", data.table=FALSE, na.strings=c("NA",""," ")), silent=TRUE)
      if (inherits(dt,"try-error")) dt <- fread(path, sep = ",", data.table=FALSE, na.strings=c("NA",""," "))
      dt
    } else if (ext %in% c("csv")) {
      fread(path, sep = ",", data.table=FALSE, na.strings=c("NA",""," "))
    } else if (ext %in% c("xls","xlsx")) {
      if (!requireNamespace("readxl", quietly = TRUE))
        stop("readxl not installed. Please install.packages('readxl')")
      as.data.frame(readxl::read_excel(path))
    } else {
      stop("Unsupported file extension: ", ext)
    }
  }, error = function(e) {
    stop(sprintf("Failed to read file: %s\nReason: %s", basename(path), e$message))
  })
}

# --- Find first file by flexible pattern (case-insensitive) ---
find_first <- function(dir, pattern) {
  if (!dir.exists(dir)) return(NA_character_)
  hits <- list.files(dir, pattern = pattern, full.names = TRUE, ignore.case = TRUE, recursive = FALSE)
  if (length(hits)) hits[1] else NA_character_
}

# --- Heuristics to find gene symbol & sample id columns ---
gene_symbol_col <- function(df){
  cands <- c("GeneSymbol","gene_symbol","Gene","SYMBOL","Hugo_Symbol")
  hit <- cands[cands %in% names(df)]
  if (length(hit)) return(hit[1])
  # fallback: first column if looks non-numeric for most rows
  n1 <- names(df)[1]
  v <- as.character(df[[n1]])
  if (mean(!is.na(v) & !grepl("^-?\\d+(\\.\\d+)?$", v)) > 0.7) return(n1)
  # else try to find a mostly non-numeric column
  for (nm in names(df)) {
    v <- as.character(df[[nm]])
    if (mean(!is.na(v) & !grepl("^-?\\d+(\\.\\d+)?$", v)) > 0.7) return(nm)
  }
  NA_character_
}

sample_id_col <- function(df){
  cands <- c("CGGA_ID","Case","Patient","Sample","Sample_ID","SampleID","ID","BarCode","barcode")
  hit <- cands[cands %in% names(df)]
  if (length(hit)) hit[1] else names(df)[1]
}

# --- Standardize expression into numeric matrix with gene symbols as rownames ---
standardize_expr <- function(expr_df){
  gcol <- gene_symbol_col(expr_df)
  if (is.na(gcol)) stop("Could not detect a gene symbol column. Please check the file format.")
  genes <- as.character(expr_df[[gcol]])
  mat <- expr_df[, setdiff(colnames(expr_df), gcol), drop=FALSE]
  # Drop completely empty columns
  if (ncol(mat) == 0) stop("No expression columns found after removing gene symbol column.")
  # Coerce numeric (safely)
  for (j in seq_len(ncol(mat))) {
    mat[[j]] <- suppressWarnings(as.numeric(mat[[j]]))
  }
  rownames(mat) <- genes
  colnames(mat) <- gsub("\\s+","", colnames(mat))
  as.matrix(mat)
}

standardize_clin <- function(clin_df){
  sid <- sample_id_col(clin_df)
  clin_df[[sid]] <- gsub("\\s+","", as.character(clin_df[[sid]]))
  clin_df
}

# --- QC helpers ---
qc_numeric <- function(expr_mat) {
  vals <- as.numeric(expr_mat)
  vals <- vals[is.finite(vals)]
  c(
    n_values = length(vals),
    min = ifelse(length(vals)>0, min(vals), NA_real_),
    max = ifelse(length(vals)>0, max(vals), NA_real_),
    pct_neg = ifelse(length(vals)>0, mean(vals < 0)*100, NA_real_),
    pct_zero = ifelse(length(vals)>0, mean(vals == 0)*100, NA_real_),
    pct_na = 100 - (length(vals) / (nrow(expr_mat)*ncol(expr_mat))) * 100
  )
}

first_nonempty <- function(x) x[which(!is.na(x) & x!="")[1]]

os_time_col <- function(df){
  cands <- c("OS","OS.time","OS_time","survival_time","Overall_Survival","OS_months","OS.days","OS_days")
  hit <- cands[cands %in% names(df)]
  first_nonempty(hit)
}

os_status_col <- function(df){
  cands <- c("Censor (alive=0; dead=1)","Status","OS.event","OS_event","OSstatus","event","vital_status","Censor","censor")
  hit <- cands[cands %in% names(df)]
  first_nonempty(hit)
}

# =====================================================================
# MAIN
# =====================================================================

cat("\n=== CGGA Data Setup (fixed/hardened) ===\n")

# 1) Scaffold + instructions
for (co in names(COHORTS)) {
  meta <- COHORTS[[co]]
  mk(meta$dir_raw); mk(meta$dir_std)
  say("\n[%s] Folders ready:", co)
  say("  Raw:      %s", meta$dir_raw)
  say("  Standard: %s", meta$dir_std)
  say("  Download from https://www.cgga.org.cn/ (login required) and place:")
  say("   • Expression: file matching pattern (case-insensitive): %s", meta$expr_pat)
  say("   • Clinical:   file matching pattern (case-insensitive): %s", meta$clin_pat)
}

# 2) If files present, verify + standardize + QC
for (co in names(COHORTS)) {
  meta <- COHORTS[[co]]
  expr_path <- find_first(meta$dir_raw, meta$expr_pat)
  clin_path <- find_first(meta$dir_raw, meta$clin_pat)

  say("\n[%s] Verification:", co)
  if (is.na(expr_path) || is.na(clin_path)) {
    say("  ✖ Files not found yet. Place them in %s and rerun.", meta$dir_raw)
    next
  }

  say("  ✓ Found expression: %s", basename(expr_path))
  say("  ✓ Found clinical:   %s", basename(clin_path))

  # Load with robust error handling
  expr_raw <- read_any_table(expr_path)
  clin_raw <- read_any_table(clin_path)

  # Standardize structures
  expr_mat <- standardize_expr(expr_raw)
  clin_std <- standardize_clin(clin_raw)

  # Basic counts
  n_genes <- nrow(expr_mat); n_samp <- ncol(expr_mat)
  say("  • Expression: %d genes × %d samples", n_genes, n_samp)

  # QC: duplicates
  dup_genes <- unique(rownames(expr_mat)[duplicated(rownames(expr_mat))])
  if (length(dup_genes)) {
    say("  ! Duplicate gene symbols: %d (first few: %s)", length(dup_genes), paste(head(dup_genes, 6), collapse=", "))
  } else say("  • No duplicate gene symbols detected.")

  dup_cols <- unique(colnames(expr_mat)[duplicated(colnames(expr_mat))])
  if (length(dup_cols)) {
    say("  ! Duplicate sample columns: %d (first few: %s)", length(dup_cols), paste(head(dup_cols, 6), collapse=", "))
  } else say("  • No duplicate sample columns detected.")

  # QC: numeric sanity
  qn <- qc_numeric(expr_mat)
  say("  • Value sanity: n=%s, min=%.4g, max=%.4g, %%neg=%.1f, %%zero=%.1f, %%NA=%.1f",
      format(qn["n_values"], big.mark=","), qn["min"], qn["max"], qn["pct_neg"], qn["pct_zero"], qn["pct_na"])

  # Panels: critical vs bonus
  genes_upper <- toupper(rownames(expr_mat))
  present_crit <- PROTECTED_CRITICAL[toupper(PROTECTED_CRITICAL) %in% genes_upper]
  missing_crit <- setdiff(PROTECTED_CRITICAL, present_crit)
  present_bonus <- PROTECTED_BONUS[toupper(PROTECTED_BONUS) %in% genes_upper]
  missing_bonus <- setdiff(PROTECTED_BONUS, present_bonus)

  if (length(missing_crit)) {
    say("  ✖ CRITICAL: Missing required gene(s): %s", paste(missing_crit, collapse=", "))
  } else {
    say("  ✓ CRITICAL genes present: %s", paste(present_crit, collapse=", "))
  }
  # Bonus are informative only
  if (length(present_bonus)) say("  • Bonus genes present: %s", paste(present_bonus, collapse=", "))
  if (length(missing_bonus)) say("  • Bonus genes missing (OK): %s", paste(missing_bonus, collapse=", "))

  # Clinical survival fields & overlap
  sid <- sample_id_col(clin_std)
  tcol <- os_time_col(clin_std)
  scol <- os_status_col(clin_std)
  say("  • Clinical sample id column: %s", sid)
  say("  • Clinical OS time/status:   %s / %s", ifelse(is.na(tcol),"NOT FOUND",tcol), ifelse(is.na(scol),"NOT FOUND",scol))

  # Overlap between expression columns & clinical IDs
  expr_ids <- colnames(expr_mat)
  clin_ids <- gsub("\\s+","", as.character(clin_std[[sid]]))
  overlap <- intersect(expr_ids, clin_ids)
  pct_overlap <- ifelse(length(expr_ids)>0, 100*length(overlap)/length(expr_ids), 0)
  say("  • Expr–Clin overlap: %d / %d samples (%.1f%%)", length(overlap), length(expr_ids), pct_overlap)

  # Survival completeness, if fields exist
  if (!is.na(tcol) && !is.na(scol)) {
    tvec <- suppressWarnings(as.numeric(clin_std[[tcol]]))
    sraw <- clin_std[[scol]]
    if (is.character(sraw)) {
      u <- toupper(sraw)
      svec <- ifelse(u %in% c("DEAD","DECEASED","1","TRUE","T"), 1L,
                     ifelse(u %in% c("ALIVE","0","FALSE","F"), 0L, NA_integer_))
    } else {
      st <- suppressWarnings(as.integer(sraw))
      uniq <- sort(unique(st[!is.na(st)]))
      svec <- if (all(uniq %in% c(0,1))) st else if (all(uniq %in% c(1,2))) ifelse(st==2L,1L,0L) else st
    }
    valid_surv <- is.finite(tvec) & !is.na(svec)
    say("  • Survival completeness: %d / %d with valid (time,status) [%.1f%%]", sum(valid_surv), length(valid_surv), 100*mean(valid_surv))
  }

  # Write standardized copies (for downstream cleanup script)
  std_expr_tsv <- file.path(meta$dir_std, paste0(co, "_expr.tsv"))
  std_clin_tsv <- file.path(meta$dir_std, paste0(co, "_clinical.tsv"))
  fwrite(data.frame(GeneSymbol=rownames(expr_mat), expr_mat, check.names=FALSE),
         std_expr_tsv, sep="\t", quote=FALSE, na = "")
  fwrite(clin_std, std_clin_tsv, sep="\t", quote=FALSE, na = "")

  # Readiness report
  rpt <- c(
    "CGGA Cohort Readiness Report (fixed)",
    paste0("Cohort: ", co),
    paste0("Date: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
    "",
    sprintf("Expression: %d genes × %d samples", n_genes, n_samp),
    if (length(dup_genes)) sprintf("Duplicate genes: %d", length(dup_genes)) else "Duplicate genes: 0",
    if (length(dup_cols))  sprintf("Duplicate samples: %d", length(dup_cols)) else "Duplicate samples: 0",
    sprintf("Value sanity: min=%.4g max=%.4g %%neg=%.1f %%zero=%.1f %%NA=%.1f",
            qn["min"], qn["max"], qn["pct_neg"], qn["pct_zero"], qn["pct_na"]),
    "",
    sprintf("CRITICAL genes present: %s", ifelse(length(missing_crit), "NO", "YES")),
    sprintf("  present: %s", ifelse(length(present_crit), paste(present_crit, collapse=", "), "None")),
    sprintf("  missing: %s", ifelse(length(missing_crit), paste(missing_crit, collapse=", "), "None")),
    "",
    sprintf("Bonus present: %s", ifelse(length(present_bonus), paste(present_bonus, collapse=", "), "None")),
    sprintf("Bonus missing (OK): %s", ifelse(length(missing_bonus), paste(missing_bonus, collapse=", "), "None")),
    "",
    sprintf("Clinical sample id: %s", sid),
    sprintf("Clinical OS time/status: %s / %s", ifelse(is.na(tcol),"NA",tcol), ifelse(is.na(scol),"NA",scol)),
    sprintf("Expr–Clin overlap: %d / %d (%.1f%%)", length(overlap), length(expr_ids), pct_overlap)
  )
  writeLines(rpt, file.path(meta$dir_std, "READINESS.txt"))

  say("  ✓ Standardized files:")
  say("    - %s", std_expr_tsv)
  say("    - %s", std_clin_tsv)
  say("  ✓ Wrote readiness report: %s", file.path(meta$dir_std, "READINESS.txt"))

  if (length(missing_crit)) {
    say("  >>> ACTION: At least one CRITICAL gene is missing. Verify you downloaded the correct CGGA expression file for %s.", co)
  }
}

cat("\n=== Done ===\n")
cat("Next:\n")
cat("  • Use the standardized files under Data/CGGA/<cohort>/Standard/ with your CGGA_GBM_cleanup_QC.R.\n")
cat("  • Critical gate is DGAT1/DGAT2 presence; bonus markers are informational only.\n")
