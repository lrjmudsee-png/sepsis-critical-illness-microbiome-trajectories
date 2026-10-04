# ============================================================
# Sepsis V2 - Step 84B4B
# PRJNA578267 FULL chimera removal with minParentAbundance = 8
#
# Rationale:
#   Step84B4A showed:
#     207 samples
#     187,235 pre-chimera ASVs
#     17,116,620 merged reads
#   and the batched sequence table was built in ~2 minutes.
#
#   Therefore the overnight bottleneck is consistent with
#   removeBimeraDenovo() on a very large pre-chimera feature set.
#
#   The DADA2 maintainer has specifically recommended
#   minParentAbundance=8 for a nearly identical ~179k-ASV
#   chimera-removal performance problem.
#
# IMPORTANT:
#   - This does NOT abundance-filter the ASV table.
#   - Low-abundance ASVs are NOT deleted before chimera calling.
#   - minParentAbundance only restricts which sequences can serve
#     as candidate chimera PARENTS.
#   - method remains "consensus".
#   - Windows: multithread=FALSE.
#
# R 4.4.0 / Windows
# ============================================================

options(stringsAsFactors = FALSE)

if (!requireNamespace("dada2", quietly=TRUE)) {
  stop("Package 'dada2' is required.")
}
if (!requireNamespace("readr", quietly=TRUE)) {
  stop("Package 'readr' is required.")
}
if (!requireNamespace("dplyr", quietly=TRUE)) {
  stop("Package 'dplyr' is required.")
}
if (!requireNamespace("tibble", quietly=TRUE)) {
  stop("Package 'tibble' is required.")
}

suppressPackageStartupMessages({
  library(dada2)
  library(readr)
  library(dplyr)
  library(tibble)
})

PROJECT <- "PRJNA578267"
ROOT <- "E:/sepsis_project"

PWORK <- file.path(
  ROOT, "data", "_V2_ANALYSIS_READY",
  "01_PROJECTS", PROJECT,
  "05_work", "step84B_full"
)

PRECHIM <- file.path(
  PWORK,
  "PRJNA578267_seqtab_prechim_step84B4A_BATCHED.rds"
)

STEP84B3 <- file.path(
  ROOT,
  "results",
  "V2_24B3_PRJNA578267_SAMPLEWISE"
)

MERGE_TRACK <- file.path(
  STEP84B3,
  "PRJNA578267_samplewise_merge_tracking.csv"
)

OUT <- file.path(
  ROOT,
  "results",
  "V2_24B4B_PRJNA578267_CHIMERA"
)

dir.create(OUT, recursive=TRUE, showWarnings=FALSE)

LOG <- file.path(
  OUT,
  "_STEP84B4B_runtime.log"
)

ERR <- file.path(
  OUT,
  "_STEP84B4B_FATAL_ERROR.txt"
)

if (file.exists(ERR)) unlink(ERR)

logmsg <- function(x) {
  z <- paste0(Sys.time(), " | ", x)
  cat(z, "\n")
  cat(z, "\n", file=LOG, append=TRUE)
}

safe_csv <- function(p) {
  suppressMessages(
    read_csv(
      p,
      show_col_types=FALSE,
      progress=FALSE,
      name_repair="unique"
    )
  )
}

valid_rds <- function(p) {
  if (!file.exists(p)) return(FALSE)
  tryCatch({
    readRDS(p)
    TRUE
  }, error=function(e) FALSE)
}

atomic_saveRDS <- function(obj, p) {
  tmp <- paste0(p, ".tmp")
  if (file.exists(tmp)) unlink(tmp)
  saveRDS(obj, tmp)
  if (file.exists(p)) unlink(p)
  if (!file.rename(tmp, p)) {
    stop("Could not atomically save: ", p)
  }
}

main <- function() {

  logmsg("STEP84B4B STARTED")

  if (!valid_rds(PRECHIM)) {
    stop("Valid Step84B4A pre-chimera sequence table not found.")
  }

  if (!file.exists(MERGE_TRACK)) {
    stop("Step84B3 samplewise merge tracking file not found.")
  }

  seqtab <- readRDS(PRECHIM)

  if (
    !is.matrix(seqtab) ||
    nrow(seqtab) != 207 ||
    ncol(seqtab) != 187235
  ) {
    stop(
      "Pre-chimera hard guard failed. Expected 207 x 187235; observed ",
      nrow(seqtab), " x ", ncol(seqtab)
    )
  }

  if (sum(seqtab) != 17116620) {
    stop(
      "Pre-chimera read-count hard guard failed. Expected 17116620; observed ",
      sum(seqtab)
    )
  }

  pre_abund <- colSums(seqtab)

  pre_diag <- tibble(
    metric=c(
      "samples",
      "pre_chimera_ASVs",
      "pre_chimera_reads",
      "ASVs_total_abundance_1",
      "ASVs_total_abundance_le_2",
      "ASVs_total_abundance_le_5",
      "ASVs_total_abundance_lt_8",
      "ASVs_total_abundance_ge_8"
    ),
    value=c(
      nrow(seqtab),
      ncol(seqtab),
      sum(seqtab),
      sum(pre_abund==1),
      sum(pre_abund<=2),
      sum(pre_abund<=5),
      sum(pre_abund<8),
      sum(pre_abund>=8)
    )
  )

  write_excel_csv(
    pre_diag,
    file.path(
      OUT,
      "PRJNA578267_PRECHIMERA_parent_candidate_diagnostic.csv"
    ),
    na=""
  )

  logmsg(
    paste0(
      "PRECHIM LOADED: 207 samples x 187235 ASVs; reads=",
      sum(seqtab),
      "; globally abundant >=8 ASVs=",
      sum(pre_abund>=8)
    )
  )

  # ----------------------------------------------------------
  # Chimera removal
  # ----------------------------------------------------------
  start_time <- Sys.time()

  logmsg(
    "removeBimeraDenovo START: method=consensus; minParentAbundance=8; multithread=FALSE"
  )

  seqtab_nochim <- removeBimeraDenovo(
    seqtab,
    method="consensus",
    minParentAbundance=8,
    multithread=FALSE,
    verbose=TRUE
  )

  end_time <- Sys.time()

  elapsed_min <- as.numeric(
    difftime(
      end_time,
      start_time,
      units="mins"
    )
  )

  logmsg(
    paste0(
      "removeBimeraDenovo COMPLETE in ",
      round(elapsed_min,2),
      " minutes"
    )
  )

  if (
    !is.matrix(seqtab_nochim) ||
    nrow(seqtab_nochim) == 0 ||
    ncol(seqtab_nochim) == 0
  ) {
    stop("Non-chimeric sequence table is invalid or empty.")
  }

  NOCHIM <- file.path(
    PWORK,
    "PRJNA578267_seqtab_nochim_step84B4B_minParent8.rds"
  )

  atomic_saveRDS(
    seqtab_nochim,
    NOCHIM
  )

  # ----------------------------------------------------------
  # Sample-level tracking
  # ----------------------------------------------------------
  merge_track <- safe_csv(MERGE_TRACK)

  if (
    nrow(merge_track) != 207 ||
    anyDuplicated(merge_track$run_id)
  ) {
    stop("Step84B3 merge tracking hard guard failed.")
  }

  nonchim <- rowSums(seqtab_nochim)

  aligned <- setNames(
    rep(0, nrow(merge_track)),
    merge_track$run_id
  )

  if (!is.null(names(nonchim))) {
    overlap <- intersect(
      names(nonchim),
      names(aligned)
    )
    aligned[overlap] <- nonchim[overlap]
  }

  final_track <- merge_track |>
    mutate(
      nonchim=as.numeric(
        aligned[run_id]
      ),
      nonchim_of_merged_pct=round(
        ifelse(
          merged>0,
          100*nonchim/merged,
          NA_real_
        ),
        2
      )
    )

  write_excel_csv(
    final_track,
    file.path(
      OUT,
      "PRJNA578267_FINAL_tracking_minParent8.csv"
    ),
    na=""
  )

  # ----------------------------------------------------------
  # ASV QC
  # ----------------------------------------------------------
  lens <- nchar(
    colnames(seqtab_nochim)
  )

  asv_qc <- tibble(
    asv_sequence=colnames(seqtab_nochim),
    asv_length=lens,
    total_abundance=colSums(seqtab_nochim),
    within_expected_length=
      lens>=350 &
      lens<=500
  )

  write_excel_csv(
    asv_qc,
    file.path(
      OUT,
      "PRJNA578267_FINAL_ASV_QC_minParent8.csv"
    ),
    na=""
  )

  total_merged <- sum(
    final_track$merged,
    na.rm=TRUE
  )

  total_nonchim <- sum(
    final_track$nonchim,
    na.rm=TRUE
  )

  read_retention <- if (total_merged>0) {
    100*total_nonchim/total_merged
  } else {
    0
  }

  zero_nonchim <- sum(
    final_track$nonchim==0 |
    is.na(final_track$nonchim)
  )

  length_ok <- round(
    100*mean(
      asv_qc$within_expected_length,
      na.rm=TRUE
    ),
    2
  )

  summary <- tibble(
    project=PROJECT,
    method="consensus",
    minParentAbundance=8,
    multithread=FALSE,
    elapsed_minutes=round(elapsed_min,2),

    samples=nrow(seqtab_nochim),

    prechim_ASVs=ncol(seqtab),
    nochim_ASVs=ncol(seqtab_nochim),

    prechim_reads=sum(seqtab),
    nochim_reads=total_nonchim,

    read_retention_after_chimera_pct=round(
      read_retention,
      2
    ),

    zero_nonchim_samples=zero_nonchim,

    median_ASV_length=median(lens),

    ASVs_350_500bp_pct=length_ok,

    chimera_read_retention_pass=
      read_retention>=50,

    zero_sample_pass=
      zero_nonchim==0,

    ASV_length_pass=
      length_ok>=80,

    full_chimera_QC_pass=
      read_retention>=50 &&
      zero_nonchim==0 &&
      length_ok>=80
  )

  write_excel_csv(
    summary,
    file.path(
      OUT,
      "PRJNA578267_FINAL_CHIMERA_SUMMARY_minParent8.csv"
    ),
    na=""
  )

  # ----------------------------------------------------------
  # Promote only if chimera QC passes
  # ----------------------------------------------------------
  if (isTRUE(summary$full_chimera_QC_pass[1])) {

    selected_dir <- file.path(
      ROOT, "data", "_V2_ANALYSIS_READY",
      "01_PROJECTS", PROJECT,
      "02_asv_selected"
    )

    dir.create(
      selected_dir,
      recursive=TRUE,
      showWarnings=FALSE
    )

    dst <- file.path(
      selected_dir,
      "PRJNA578267_seqtab_final_step84B4B_minParent8.rds"
    )

    if (file.exists(dst)) unlink(dst)

    ok <- file.link(
      NOCHIM,
      dst
    )

    if (!isTRUE(ok)) {
      stop("Final chimera-cleaned ASV passed QC but promotion hard-link failed.")
    }

    file.copy(
      file.path(
        OUT,
        "PRJNA578267_FINAL_tracking_minParent8.csv"
      ),
      file.path(
        selected_dir,
        "PRJNA578267_step84B4B_tracking.csv"
      ),
      overwrite=TRUE
    )

    file.copy(
      file.path(
        OUT,
        "PRJNA578267_FINAL_ASV_QC_minParent8.csv"
      ),
      file.path(
        selected_dir,
        "PRJNA578267_step84B4B_ASV_QC.csv"
      ),
      overwrite=TRUE
    )

    status_text <- "PASS_PROMOTED_TO_SELECTED_ASV"

  } else {

    status_text <- "COMPLETE_BUT_QC_REVIEW_REQUIRED"
  }

  writeLines(
    c(
      paste0("Completed: ", Sys.time()),
      paste0("Status: ", status_text),
      paste0("Method: consensus"),
      paste0("minParentAbundance: 8"),
      paste0("Elapsed minutes: ", round(elapsed_min,2)),
      paste0("Prechim ASVs: ", ncol(seqtab)),
      paste0("Nochim ASVs: ", ncol(seqtab_nochim)),
      paste0(
        "Read retention after chimera removal: ",
        round(read_retention,2),
        "%"
      ),
      paste0(
        "Zero nonchim samples: ",
        zero_nonchim
      ),
      paste0(
        "ASVs 350-500 bp: ",
        length_ok,
        "%"
      )
    ),
    file.path(
      OUT,
      "_STEP84B4B_COMPLETE.ok"
    )
  )

  logmsg(
    paste0(
      "STEP84B4B COMPLETE | status=",
      status_text
    )
  )

  cat("\n============================================================\n")
  cat("PRJNA578267 CHIMERA REMOVAL COMPLETE\n")
  cat("============================================================\n")
  print(summary,n=Inf,width=Inf)
  cat("\nOutput:\n",OUT,"\n",sep="")
}

tryCatch(
  main(),
  error=function(e) {

    msg <- c(
      paste0(
        "STEP84B4B FATAL ERROR: ",
        Sys.time()
      ),
      paste0(
        "Message: ",
        conditionMessage(e)
      ),
      paste0(
        "Call: ",
        paste(
          deparse(
            conditionCall(e)
          ),
          collapse=" "
        )
      )
    )

    writeLines(
      msg,
      ERR,
      useBytes=TRUE
    )

    message(
      paste(
        msg,
        collapse="\n"
      )
    )

    quit(
      save="no",
      status=1,
      runLast=FALSE
    )
  }
)
