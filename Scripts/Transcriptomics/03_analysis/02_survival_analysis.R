#!/usr/bin/env Rscript
# =============================================================================
# 01_bulk_survival_bestcut_auto.R
# - Auto-detect clinical columns (time, status, age, IDH, MGMT)
# - Best cutoff (survminer::surv_cutpoint) for DGAT1/2
# - KM plots + Cox models (drop-in if covariates missing)
# - Saves PNG/PDF plots, Cox tables, and a run log
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(stringr)
  library(survival)
  library(survminer)
})

dir.create("Results/Bulk/Survival/TCGA_GBM", recursive = TRUE, showWarnings = FALSE)
dir.create("Results/Bulk/Survival/CGGA_GBM", recursive = TRUE, showWarnings = FALSE)
log_fp <- "Results/Bulk/Survival/survival_step1.log"

log_msg <- function(...) {
  msg <- paste0(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), " | ", paste0(..., collapse = ""))
  cat(msg, "\n")
  write(msg, file = log_fp, append = TRUE)
}

# ---------- Helpers: flexible clinical mapping ----------
.standardize_status <- function(x) {
  # convert character/factor to numeric 0/1 where 1 = event/death
  if (is.numeric(x)) {
    # Assume already 0/1, but map any values >1 to 1
    return(ifelse(x > 1, 1, x))
  }
  xv <- tolower(as.character(x))
  y <- ifelse(xv %in% c("dead","deceased","1","event","yes","true"), 1,
       ifelse(xv %in% c("alive","0","no","false","censored"), 0, NA_real_))
  return(y)
}

.pick_first_present <- function(df, candidates) {
  hits <- intersect(tolower(candidates), tolower(colnames(df)))
  if (length(hits) == 0) return(NULL)
  # return original-case column name
  colnames(df)[match(hits[1], tolower(colnames(df)))]
}

.guess_time_status <- function(clin) {
  # time (days or months) common patterns
  time_cands <- c(
    "os.time","os_time","overall_survival_time","overall_survival",
    "days_to_death","days_to_last_follow_up",
    "os_months","pfs.time","pfs_time","pfs_months"
  )
  status_cands <- c("os","status","vital_status","event","os_event","pfs","pfs_event")

  time_col <- .pick_first_present(clin, time_cands)
  status_col <- .pick_first_present(clin, status_cands)

  # Fallback: numeric longest-looking column for time, and a binary-looking column for status
  if (is.null(time_col)) {
    num_cols <- names(which(sapply(clin, is.numeric)))
    if (length(num_cols)) {
      # choose one with largest range as time
      rng <- sapply(clin[num_cols], function(v) diff(range(v, na.rm=TRUE)))
      time_col <- num_cols[which.max(rng)]
    }
  }
  if (is.null(status_col)) {
    # try to find categorical with 2 levels
    two_level <- names(which(sapply(clin, function(v) {
      u <- unique(na.omit(as.character(v))); length(u) %in% c(2,3) # Dead/Alive/NA
    })))
    # prefer columns with keywords
    key <- two_level[str_detect(tolower(two_level), "status|event|vital|os")]
    status_col <- if (length(key)) key[1] else if (length(two_level)) two_level[1] else NULL
  }

  list(time_col = time_col, status_col = status_col)
}

.guess_age <- function(clin) {
  .pick_first_present(clin, c("age","age_at_initial_pathologic_diagnosis","age_at_diagnosis","age_at_index","diagnosis_age"))
}

.guess_idh <- function(clin) {
  .pick_first_present(clin, c("idh.status","idh_status","idh","idh1_status","idh1"))
}

.guess_mgmt <- function(clin) {
  .pick_first_present(clin, c("mgmt.status","mgmt_status","mgmt","mgmt.promoter.status","mgmtp_methylation","mgmt_methylation"))
}

.make_numeric <- function(v) {
  if (is.numeric(v)) return(v)
  suppressWarnings(as.numeric(as.character(v)))
}

.normalize_time <- function(time, colname) {
  # Heuristic: if column name contains "month", convert months -> days (approx 30.44)
  if (is.null(colname)) return(time)
  if (str_detect(tolower(colname), "month")) {
    return(time * 30.44)
  }
  time
}

# ---------- Core runner ----------
run_survival <- function(expr_mat, clin_df, gene, cohort, outdir) {
  if (!(gene %in% rownames(expr_mat))) {
    log_msg(cohort, " | ", gene, ": gene not found — skipping.")
    return(invisible(NULL))
  }

  # sample match
  # try: columns of expr_mat are barcodes; clin has rownames or a sample column
  if (is.null(rownames(clin_df)) || length(intersect(colnames(expr_mat), rownames(clin_df))) < 5) {
    # attempt to set rownames from a likely sample column
    cand_ids <- c("sample","sample_id","patient","patient_id","submitter_id","barcode","case_id")
    idcol <- .pick_first_present(clin_df, cand_ids)
    if (!is.null(idcol)) rownames(clin_df) <- make.unique(as.character(clin_df[[idcol]]))
  }
  samples <- intersect(colnames(expr_mat), rownames(clin_df))
  if (length(samples) < 30) {
    log_msg(cohort, " | ", gene, ": too few matched samples (", length(samples), ") — skipping.")
    return(invisible(NULL))
  }
  expr <- expr_mat[gene, samples, drop=FALSE]
  clin <- clin_df[samples, , drop=FALSE]

  # detect key columns
  ts <- .guess_time_status(clin)
  age_col  <- .guess_age(clin)
  idh_col  <- .guess_idh(clin)
  mgmt_col <- .guess_mgmt(clin)

  if (is.null(ts$time_col) || is.null(ts$status_col)) {
    log_msg(cohort, " | ", gene, ": could not detect time/status columns — skipping.")
    return(invisible(NULL))
  }

  time  <- .make_numeric(clin[[ts$time_col]])
  time  <- .normalize_time(time, ts$time_col)
  status <- .standardize_status(clin[[ts$status_col]])

  df <- data.frame(
    expr = as.numeric(expr),
    time = time,
    status = status,
    stringsAsFactors = FALSE
  )
  if (!is.null(age_col))  df$age  <- .make_numeric(clin[[age_col]])
  if (!is.null(idh_col))  df$idh  <- as.character(clin[[idh_col]])
  if (!is.null(mgmt_col)) df$mgmt <- as.character(clin[[mgmt_col]])

  # drop rows with missing essential fields
  df <- df %>% filter(!is.na(time), !is.na(status), !is.na(expr))
  if (nrow(df) < 30) {
    log_msg(cohort, " | ", gene, ": insufficient complete cases after QC — skipping.")
    return(invisible(NULL))
  }

  # best cutoff (ensure reasonable group sizes with minprop)
  cut <- surv_cutpoint(df, time = "time", event = "status",
                       variables = "expr", minprop = 0.1) # at least 10% in each group
  group <- surv_categorize(cut)
  df$group <- group$expr

  # Save cutoff & group sizes
  cut_val <- cut$cutpoint$cutpoint[1]
  sizes <- table(df$group)
  info <- data.frame(
    cohort = cohort, gene = gene,
    cutoff = cut_val,
    n_low = as.integer(sizes["low"]), n_high = as.integer(sizes["high"]),
    time_col = ts$time_col, status_col = ts$status_col,
    age_col = age_col %||% NA_character_,
    idh_col = idh_col %||% NA_character_,
    mgmt_col = mgmt_col %||% NA_character_
  )
  fwrite(info, file.path(outdir, paste0(gene, "_bestcut_info.csv")))

  # KM
  fit <- survfit(Surv(time, status) ~ group, data = df)
  g <- ggsurvplot(
    fit, data = df, pval = TRUE, risk.table = TRUE,
    title = paste0(cohort, " - ", gene, " (best cutoff=", signif(cut_val, 4), ")"),
    legend.title = gene,
    legend.labs = c("Low", "High")
  )
  ggsave(file.path(outdir, paste0(gene, "_KM.png")), g$plot, width = 6, height = 5, dpi = 300)
  ggsave(file.path(outdir, paste0(gene, "_KM.pdf")), g$plot, width = 6, height = 5)

  # Cox models: start with base + add covariates if present
  form <- as.formula("Surv(time, status) ~ expr")
  if ("age" %in% names(df) && sum(!is.na(df$age)) > 0) form <- update(form, . ~ . + age)
  if ("idh" %in% names(df) && length(na.omit(unique(df$idh))) > 1) form <- update(form, . ~ . + idh)
  if ("mgmt" %in% names(df) && length(na.omit(unique(df$mgmt))) > 1) form <- update(form, . ~ . + mgmt)

  cox <- coxph(form, data = df)
  cs  <- summary(cox)

  # Tidy-ish export without extra deps
  coefs <- as.data.frame(cs$coefficients)
  coefs$`HR` <- exp(coefs$coef)
  if (!is.null(cs$conf.int)) {
    ci <- as.data.frame(cs$conf.int)[, c("lower .95","upper .95")]
    names(ci) <- c("HR_lower95","HR_upper95")
    coefs <- cbind(Variable = rownames(coefs), coefs[, c("coef","exp(coef)","se(coef)","z","Pr(>|z|)")], ci)
    names(coefs)[2:5] <- c("coef","HR","se","z")
    names(coefs)[6]   <- c("p")
  } else {
    coefs <- cbind(Variable = rownames(coefs), coefs)
  }
  rownames(coefs) <- NULL
  fwrite(coefs, file.path(outdir, paste0(gene, "_Cox.csv")))

  log_msg(cohort, " | ", gene, ": OK (n=", nrow(df), ", cutoff=", signif(cut_val,4),
          ", groups: low=", sizes["low"], ", high=", sizes["high"], ")")
}

`%||%` <- function(a, b) if (!is.null(a)) a else b

# ---------- Load data ----------
# TCGA
tcga_expr_fp <- "Processed_Data/Bulk/TCGA_GBM_Expr_HGNC_FPKM.rds"
tcga_clin_fp <- "Processed_Data/Bulk/TCGA_GBM_Clinical.csv"
# CGGA (via GEO lookalike you downloaded)
cgga_expr_fp <- "Processed_Data/Bulk/CGGA_like_GSE16011_Expr_HGNC.rds"
cgga_clin_fp <- "Processed_Data/Bulk/CGGA_like_GSE16011_Clinical.csv"

stopifnot(file.exists(tcga_expr_fp), file.exists(tcga_clin_fp))
stopifnot(file.exists(cgga_expr_fp), file.exists(cgga_clin_fp))

tcga_expr <- readRDS(tcga_expr_fp)
tcga_clin <- as.data.frame(fread(tcga_clin_fp))
rownames(tcga_clin) <- tcga_clin$sample_id

cgga_expr <- readRDS(cgga_expr_fp)
cgga_clin <- as.data.frame(fread(cgga_clin_fp))

# If expression is samples x genes, transpose to genes x samples
if (nrow(tcga_expr) < ncol(tcga_expr)) tcga_expr <- t(tcga_expr)
if (nrow(cgga_expr) < ncol(cgga_expr)) cgga_expr <- t(cgga_expr)

# ---------- Run ----------
genes <- c("DGAT1","DGAT2")

for (g in genes) {
  run_survival(tcga_expr, tcga_clin, g, "TCGA_GBM", "Results/Bulk/Survival/TCGA_GBM")
  run_survival(cgga_expr, cgga_clin, g, "CGGA_GBM", "Results/Bulk/Survival/CGGA_GBM")
}

log_msg("DONE: Survival analysis saved under Results/Bulk/Survival/")
cat("\n🎉 Survival analysis complete. See Results/Bulk/Survival/ and ", log_fp, "\n")
