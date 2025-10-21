#!/usr/bin/env Rscript
# =====================================================================
# CGGA mRNAseq_693 — GBM-only cleanup + QC (publication-grade)
# - Filters to GBM (WHO IV / "GBM")
# - Deduplicates to 1 sample per patient (Primary > else DGAT1 > SOAT1)
# - Robust OS status/time detection + flip check
# - Event-rate + median OS validation (GBM expectations)
# - Extracts IDH, MGMT, Age, Sex, Grade
# - Protects DGAT1/DGAT2/SOAT1 from low-expression filtering
# =====================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(survival)
  library(survminer)
  library(ggplot2)
})

# -------------------- USER SETTINGS --------------------
EXPR_FILE <- "../Raw_Data/CGGA/CGGA_693_Expression_Cleaned.tsv"
CLIN_FILE <- "../Raw_Data/CGGA/CGGA_693_Clinical_Cleaned.tsv"
OUT_DIR   <- "../Processed_Data/Clean_CGGA_Cohort/mRNAseq_693_GBM"
LOW_EXPR_MIN_PROP <- 0.10     # keep genes detected (>0) in ≥10% of samples
PROTECT_GENES     <- c("DGAT1","DGAT2","SOAT1")  # never drop these
TIE_GENES         <- c("DGAT1","SOAT1")          # for dedup tie-break
EXPECTED_EVENT_RANGE <- c(0.60, 0.85)            # GBM OS expected event rate
set.seed(1)
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

say <- function(...) cat(sprintf(...), "\n")
`%||%` <- function(a,b) if (is.null(a) || length(a)==0 || all(is.na(a))) b else a

# -------------------- HELPERS --------------------
read_tab <- function(p){
  stopifnot(file.exists(p))
  ext <- tolower(tools::file_ext(p))
  if (ext %in% c("tsv","txt")) fread(p, sep="\t", data.table=FALSE)
  else if (ext=="csv") fread(p, sep=",", data.table=FALSE)
  else stop("Unsupported extension for: ", p)
}

gene_symbol_col <- function(df){
  cands <- c("GeneSymbol","gene_symbol","Gene","SYMBOL","Hugo_Symbol")
  hit <- cands[cands %in% names(df)]
  if (length(hit)) hit[1] else names(df)[1]
}
sample_id_col <- function(df){
  cands <- c("CGGA_ID","Case","Patient","Sample","Sample_ID","SampleID","ID")
  hit <- cands[cands %in% names(df)]
  if (length(hit)) hit[1] else names(df)[1]
}
grade_col <- function(df){
  cands <- c("Grade","WHO.grade","WHO_Grade","WHOgrade","WHO")
  hit <- cands[cands %in% names(df)]
  if (length(hit)) hit[1] else NA_character_
}
sample_type_col <- function(df){
  cands <- c("Sample_Type","SampleType","Type","Tumor.type","Tumor_Type")
  hit <- cands[cands %in% names(df)]
  if (length(hit)) hit[1] else NA_character_
}
idh_col <- function(df){
  cands <- c("IDH.status","IDH_status","IDH","IDH1.status","IDH1_status","IDH1")
  cands[cands %in% names(df)][1] %||% NA_character_
}
mgmt_col <- function(df){
  cands <- c("MGMT.status","MGMT_status","MGMT.methylation.status","MGMT","MGMTp")
  cands[cands %in% names(df)][1] %||% NA_character_
}
age_col <- function(df){
  cands <- c("Age","age","Age.at.diagnosis","Age_at_diagnosis")
  cands[cands %in% names(df)][1] %||% NA_character_
}
sex_col <- function(df){
  cands <- c("Gender","Sex","gender","sex")
  cands[cands %in% names(df)][1] %||% NA_character_
}
os_time_col <- function(df){
  cands <- c("OS","OS.time","OS_time","Overall_Survival","survival_time","OS_days","OS.days","OS_months")
  cands[cands %in% names(df)][1] %||% NA_character_
}
os_status_col <- function(df){
  # Many CGGA tables: "Censor (alive=0; dead=1)" or "OS.event"
  cands <- c("Censor (alive=0; dead=1)","OS.event","OS_event","OSstatus","Status","vital_status","event")
  cands[cands %in% names(df)][1] %||% NA_character_
}

std_expr <- function(df){
  gcol <- gene_symbol_col(df)
  genes <- as.character(df[[gcol]])
  mat <- as.matrix(df[, setdiff(names(df), gcol), drop=FALSE])
  storage.mode(mat) <- "double"
  rownames(mat) <- genes
  colnames(mat) <- gsub("\\s+","", colnames(mat))
  mat
}

# -------------------- LOAD --------------------
say("\nLoading CGGA mRNAseq_693 ...")
expr_raw <- read_tab(EXPR_FILE)
clin_raw <- read_tab(CLIN_FILE)

expr <- std_expr(expr_raw)
clin <- clin_raw
sid <- sample_id_col(clin)
clin[[sid]] <- gsub("\\s+","", as.character(clin[[sid]]))

say("Expression: %d genes × %d samples", nrow(expr), ncol(expr))
say("Clinical rows: %d; ID column: %s", nrow(clin), sid)

# -------------------- ALIGN --------------------
common_ids <- intersect(colnames(expr), clin[[sid]])
if (length(common_ids) < 100) stop("Very low expr/clin overlap: ", length(common_ids))
expr <- expr[, common_ids, drop=FALSE]
clin <- clin[match(common_ids, clin[[sid]]), , drop=FALSE]

# -------------------- GBM FILTER --------------------
gcol <- grade_col(clin)
if (is.na(gcol)) stop("No grade column found in clinical.")
grade_txt <- toupper(gsub("\\s+","", as.character(clin[[gcol]])))
is_gbm <- grepl("GBM|IV|GRADEIV|GLIOBLASTOMA", grade_txt)
if (!any(is_gbm)) stop("No GBM samples detected by grade text.")
expr <- expr[, is_gbm, drop=FALSE]
clin <- clin[is_gbm, , drop=FALSE]
say("GBM subset: %d samples", ncol(expr))

# -------------------- DEDUP (PATIENT-LEVEL) --------------------
# Prefer Sample_Type 'Primary' if available; else tie by DGAT1, then SOAT1
stype_col <- sample_type_col(clin)
patient_key <- sid  # CGGA_ID is typically per-patient; if not, adjust here

# Compute tie vectors
tie_mats <- lapply(TIE_GENES, function(g) if (g %in% rownames(expr)) expr[g,] else rep(-Inf, ncol(expr)))
tie_agg <- Reduce(function(a,b) ifelse(is.finite(a)|is.finite(b), pmax(a,b,na.rm=TRUE), -Inf), tie_mats)

# group by patient_key (CGGA_ID). If CGGA_ID is per sample, this still works 1:1
idx_keep <- unlist(lapply(split(seq_len(ncol(expr)), clin[[patient_key]]), function(ix){
  if (length(ix)==1) return(ix)
  # prefer primary if stype available
  if (!is.na(stype_col)) {
    st <- toupper(as.character(clin[[stype_col]][ix]))
    prim_ix <- ix[grepl("PRIMARY|UNTREATED|INITIAL", st)]
    if (length(prim_ix)==1) return(prim_ix)
    if (length(prim_ix) > 1) ix <- prim_ix
  }
  # tie-break by DGAT1 -> SOAT1
  ix[which.max(tie_agg[ix])]
}))
idx_keep <- sort(unique(idx_keep))

expr <- expr[, idx_keep, drop=FALSE]
clin <- clin[idx_keep, , drop=FALSE]
say("After patient-level dedup: %d patients", ncol(expr))

# -------------------- LOW-EXPR FILTER (protect key genes) --------------------
det_prop <- rowMeans(expr > 0, na.rm=TRUE)
keep <- det_prop >= LOW_EXPR_MIN_PROP
keep[match(PROTECT_GENES, rownames(expr), nomatch=0)] <- TRUE
expr <- expr[keep, , drop=FALSE]
say("Genes after low-expression filter (protected=%s): %d",
    paste(PROTECT_GENES, collapse=","), nrow(expr))

# -------------------- SURVIVAL ENCODING --------------------
tcol <- os_time_col(clin)
scol <- os_status_col(clin)
if (is.na(tcol) || is.na(scol)) stop("Cannot find OS time/status columns.")

time_raw <- suppressWarnings(as.numeric(clin[[tcol]]))
# Detect units: days vs months vs years (heuristic)
med_t <- median(time_raw, na.rm=TRUE)
time_years <- if (is.na(med_t)) rep(NA_real_, length(time_raw)) else {
  if (med_t > 365) time_raw/365.25 else if (med_t > 24) time_raw/12 else time_raw
}

# Map status to 0/1
sraw <- clin[[scol]]
if (is.character(sraw)) {
  u <- toupper(sraw)
  status <- ifelse(u %in% c("DEAD","DECEASED","1","TRUE","T"), 1L,
                   ifelse(u %in% c("ALIVE","0","FALSE","F"), 0L, NA_integer_))
} else {
  st <- suppressWarnings(as.integer(sraw))
  uniq <- sort(unique(st[!is.na(st)]))
  status <- if (all(uniq %in% c(0,1))) st else if (all(uniq %in% c(1,2))) ifelse(st==2L,1L,0L) else st
}

valid <- is.finite(time_years) & !is.na(status) & time_years > 0
expr <- expr[, valid, drop=FALSE]
clin <- clin[valid, , drop=FALSE]
time_years <- time_years[valid]; status <- status[valid]

# Flip detection: dead should have shorter median time
med_alive <- median(time_years[status==0], na.rm=TRUE)
med_dead  <- median(time_years[status==1], na.rm=TRUE)
flipped <- FALSE
if (!is.na(med_alive) && !is.na(med_dead) && med_dead > med_alive) {
  status <- 1L - status
  flipped <- TRUE
  # recompute medians
  med_alive <- median(time_years[status==0], na.rm=TRUE)
  med_dead  <- median(time_years[status==1], na.rm=TRUE)
}

# -------------------- QC STATS --------------------
evt_rate <- mean(status==1)
med_os_events <- median(time_years[status==1], na.rm=TRUE)
expect_ok <- !is.na(evt_rate) && evt_rate >= EXPECTED_EVENT_RANGE[1] && evt_rate <= EXPECTED_EVENT_RANGE[2]

say("\nQC SUMMARY:")
say("  Final patients: %d", ncol(expr))
say("  Events: %d (%.1f%%) | Expected GBM: %d–%d%% | %s",
    sum(status==1), 100*evt_rate, round(100*EXPECTED_EVENT_RANGE[1]), round(100*EXPECTED_EVENT_RANGE[2]),
    ifelse(expect_ok,"PASS","WARN"))
say("  Median OS (events): %.2f years", med_os_events)
say("  Status flipped?: %s", ifelse(flipped,"YES","NO"))

# -------------------- KEY CLINICAL FIELDS --------------------
g_keep <- grade_txt[valid]  # from earlier
clin$Grade_clean <- gsub("\\s+"," ", as.character(clin[[gcol]]))

icol <- idh_col(clin)
mcol <- mgmt_col(clin)
acol <- age_col(clin)
xcol <- sex_col(clin)

clin$IDH_status  <- if (!is.na(icol)) as.character(clin[[icol]]) else NA
clin$MGMT_status <- if (!is.na(mcol)) as.character(clin[[mcol]]) else NA
clin$Age         <- suppressWarnings(as.numeric(if (!is.na(acol)) clin[[acol]] else NA))
clin$Sex         <- if (!is.na(xcol)) as.character(clin[[xcol]]) else NA

# -------------------- SAVE --------------------
OUT_EXPR <- file.path(OUT_DIR, "CGGA_mRNAseq_693_GBM_expr_clean.rds")
OUT_CLIN <- file.path(OUT_DIR, "CGGA_mRNAseq_693_GBM_clin_clean.csv")
OUT_TRACK<- file.path(OUT_DIR, "CGGA_mRNAseq_693_GBM_sample_tracking.csv")
OUT_QC   <- file.path(OUT_DIR, "CGGA_mRNAseq_693_GBM_QC.txt")

saveRDS(expr, OUT_EXPR)

clin_out <- data.frame(
  CGGA_ID = clin[[sid]],
  Grade   = clin$Grade_clean,
  IDH_status = clin$IDH_status,
  MGMT_status= clin$MGMT_status,
  Age = clin$Age,
  Sex = clin$Sex,
  OS_time_years = time_years,
  OS_status = status,
  stringsAsFactors = FALSE
)
fwrite(clin_out, OUT_CLIN)

fwrite(data.frame(
  sample = colnames(expr),
  CGGA_ID = clin[[sid]],
  OS_time_years = time_years,
  OS_status = status
), OUT_TRACK)

qc_lines <- c(
  "CGGA mRNAseq_693 — GBM Cleanup + QC",
  sprintf("Patients: %d", ncol(expr)),
  sprintf("Events: %d (%.1f%%) | Expected: %d–%d%% -> %s",
          sum(status==1), 100*evt_rate, round(100*EXPECTED_EVENT_RANGE[1]),
          round(100*EXPECTED_EVENT_RANGE[2]), ifelse(expect_ok,"PASS","WARN")),
  sprintf("Median OS (events): %.2f years", med_os_events),
  sprintf("Status flipped: %s", ifelse(flipped,"YES","NO")),
  sprintf("Genes kept: %d (low-expr min prop=%.2f; protected=%s)",
          nrow(expr), LOW_EXPR_MIN_PROP, paste(PROTECT_GENES, collapse=",")),
  sprintf("Dedup: Primary preference=%s; tie genes=%s",
          ifelse(!is.na(stype_col),"YES","NO"), paste(TIE_GENES, collapse=" → ")),
  sprintf("Key fields present: IDH=%s, MGMT=%s, Age=%s, Sex=%s",
          !is.na(icol), !is.na(mcol), !is.na(acol), !is.na(xcol))
)
writeLines(qc_lines, OUT_QC)

say("\nSaved:")
say("  %s", OUT_EXPR)
say("  %s", OUT_CLIN)
say("  %s", OUT_TRACK)
say("  %s", OUT_QC)

# -------------------- OPTIONAL QUICK KMs --------------------
quick_km <- function(gene, file_stub){
  if (!(gene %in% rownames(expr))) return(invisible(NULL))
  df <- data.frame(expr=as.numeric(expr[gene,]), time=time_years, status=status)
  df <- df[is.finite(df$expr) & is.finite(df$time), ]
  scp <- try(survminer::surv_cutpoint(df, time="time", event="status", variables="expr", minprop=0.1), silent=TRUE)
  if (inherits(scp,"try-error")) return(invisible(NULL))
  cut <- scp$cutpoint$cutpoint[1]
  grp <- factor(ifelse(df$expr >= cut, "High","Low"), levels=c("Low","High"))
  fit <- survfit(Surv(df$time, df$status) ~ grp)
  sdif <- survdiff(Surv(df$time, df$status) ~ grp)
  p_lr <- pchisq(sdif$chisq, 1, lower.tail=FALSE)
  plt <- ggsurvplot(
    fit, data=df, palette=c("#3B4CC0","#E8601C"), conf.int=FALSE,
    risk.table=TRUE, legend.title=paste0(gene," (RNA)"),
    xlab="Time (years)", ylab="Survival probability",
    title=paste0("CGGA mRNAseq_693 GBM — ", gene, " vs OS"),
    ggtheme=theme_classic()
  )
  plt$plot <- plt$plot + annotate("text", x=0.1, y=0.15,
                                  label=sprintf("cut=%.2f  p=%.3f", cut, p_lr),
                                  hjust=0, vjust=0, size=3.5, fontface="bold", color="gray20")
  out_png <- file.path(OUT_DIR, paste0("KM_", file_stub, ".png"))
  ggsave(out_png, plt$plot, width=8, height=6, dpi=300)
  say("  KM plot saved: %s", out_png)
}
say("\nQuick KMs (optional):")
quick_km("SOAT1","SOAT1")
quick_km("DGAT1","DGAT1")
