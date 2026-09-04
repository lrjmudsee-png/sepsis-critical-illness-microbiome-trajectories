# ============================================================
# Sepsis V2 - Step 84B
# FULL-COHORT DADA2 for all pilot-approved Illumina cohorts
#
# STRICT INPUT:
#   E:/sepsis_project/data/_V2_ANALYSIS_READY
#
# Projects:
#   PRJEB82425
#   PRJNA516701
#   PRJNA578267
#   PRJNA851469
#   PRJNA1166732
#
# PRJNA1166732 override from Step84A2:
#   minOverlap = 8
#   maxMismatch = 0
#
# IMPORTANT
#   - Cohorts are processed independently.
#   - No broad filesystem scan.
#   - No HOLD / PacBio / Ion Torrent / FASTA-only cohort is touched.
#   - A project-level COMPLETE marker enables safe resume.
#   - Final ASV is promoted to 02_asv_selected only if full-cohort QC passes.
#
# Usage:
#   Run all projects:
#     Rscript 84B_run_full_cohort_DADA2.R
#
#   Run one project only:
#     Rscript 84B_run_full_cohort_DADA2.R PRJNA578267
#
# R 4.4.0 / Windows
# ============================================================

options(stringsAsFactors = FALSE)

# ------------------------------------------------------------
# Packages
# ------------------------------------------------------------
cran_pkgs <- c("readr","dplyr","stringr","purrr","tibble")
missing_cran <- cran_pkgs[
  !vapply(cran_pkgs, requireNamespace, logical(1), quietly=TRUE)
]

if (length(missing_cran)) {
  install.packages(missing_cran, repos="https://cloud.r-project.org")
}

if (!requireNamespace("dada2", quietly=TRUE)) {
  if (!requireNamespace("BiocManager", quietly=TRUE)) {
    install.packages("BiocManager", repos="https://cloud.r-project.org")
  }
  BiocManager::install("dada2", ask=FALSE, update=FALSE)
}

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(stringr)
  library(purrr)
  library(tibble)
  library(dada2)
})

# ------------------------------------------------------------
# Paths
# ------------------------------------------------------------
ROOT <- "E:/sepsis_project"

CANON <- file.path(
  ROOT,
  "data",
  "_V2_ANALYSIS_READY"
)

PROJECT_ROOT <- file.path(
  CANON,
  "01_PROJECTS"
)

FREEZE_ROOT <- file.path(
  CANON,
  "00_FREEZE"
)

STEP83B_ROOT <- file.path(
  ROOT,
  "results",
  "V2_23B_READ_QUALITY_PRIMER_AUDIT"
)

STEP84A_ROOT <- file.path(
  ROOT,
  "results",
  "V2_24A_DADA2_PILOT"
)

STEP84A2_ROOT <- file.path(
  ROOT,
  "results",
  "V2_24A2_PRJNA1166732_MERGE_DIAGNOSTIC"
)

OUT <- file.path(
  ROOT,
  "results",
  "V2_24B_FULL_DADA2"
)

dir.create(OUT, recursive=TRUE, showWarnings=FALSE)

LOG <- file.path(OUT, "_STEP84B_runtime_checkpoints.txt")
ERR <- file.path(OUT, "_STEP84B_FATAL_ERROR.txt")
if (file.exists(ERR)) unlink(ERR)

ck <- function(x) {
  cat(
    paste0(x, ": ", Sys.time(), "\n"),
    file=LOG,
    append=TRUE
  )
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

extract_run <- function(x) {
  str_extract(
    basename(x),
    regex("(SRR|ERR|DRR|CRR)[0-9]+", ignore_case=TRUE)
  )
}

getN <- function(x) {
  sum(getUniques(x))
}

primer_len <- function(x) {
  lens <- c(
    "341F" = 17L,
    "785R" = 21L,
    "806R" = 20L,
    "515F" = 19L
  )

  if (is.na(x) || !(x %in% names(lens))) return(0L)
  unname(lens[[x]])
}

length_range <- function(project) {

  switch(
    project,

    "PRJEB82425"   = c(350L,500L),
    "PRJNA578267"  = c(350L,500L),
    "PRJNA851469"  = c(200L,330L),
    "PRJNA1166732" = c(200L,330L),
    "PRJNA516701"  = c(150L,500L),

    c(100L,600L)
  )
}

safe_remove <- function(paths) {
  paths <- paths[file.exists(paths)]
  if (length(paths)) unlink(paths, recursive=TRUE, force=TRUE)
}

# ------------------------------------------------------------
# Exact pair table from canonical inventory
# ------------------------------------------------------------
build_pair_table <- function(project, inventory, manifest) {

  z <- inventory |>
    filter(
      project==!!project,
      mate %in% c("R1","R2")
    ) |>
    mutate(
      run_id=extract_run(file_name)
    ) |>
    filter(!is.na(run_id))

  r1 <- z |>
    filter(mate=="R1") |>
    select(run_id, R1=file_path) |>
    distinct(run_id,.keep_all=TRUE)

  r2 <- z |>
    filter(mate=="R2") |>
    select(run_id, R2=file_path) |>
    distinct(run_id,.keep_all=TRUE)

  inner_join(r1,r2,by="run_id") |>
    inner_join(
      manifest |>
        filter(project==!!project) |>
        select(
          project,
          run_id,
          patient_id,
          sample_id,
          time_raw,
          time_day
        ),
      by="run_id"
    ) |>
    arrange(
      patient_id,
      time_day,
      time_raw,
      run_id
    )
}

# ------------------------------------------------------------
# Promote only QC-passed ASV to canonical selected directory
# ------------------------------------------------------------
promote_final_asv <- function(project, source_rds, source_track, source_qc) {

  selected_dir <- file.path(
    PROJECT_ROOT,
    project,
    "02_asv_selected"
  )

  dir.create(
    selected_dir,
    recursive=TRUE,
    showWarnings=FALSE
  )

  dst_rds <- file.path(
    selected_dir,
    paste0(project,"_seqtab_final_step84B.rds")
  )

  dst_track <- file.path(
    selected_dir,
    paste0(project,"_step84B_tracking.csv")
  )

  dst_qc <- file.path(
    selected_dir,
    paste0(project,"_step84B_ASV_length_QC.csv")
  )

  # Do not silently overwrite a different prior object.
  if (file.exists(dst_rds)) {
    unlink(dst_rds)
  }

  ok <- file.link(
    source_rds,
    dst_rds
  )

  if (!isTRUE(ok)) {
    stop(
      paste0(
        project,
        ": failed to hard-link final ASV into 02_asv_selected."
      )
    )
  }

  file.copy(
    source_track,
    dst_track,
    overwrite=TRUE
  )

  file.copy(
    source_qc,
    dst_qc,
    overwrite=TRUE
  )

  dst_rds
}

# ------------------------------------------------------------
# One full cohort
# ------------------------------------------------------------
run_full_project <- function(
  project,
  param,
  inventory,
  manifest
) {

  message("\n============================================================")
  message("STEP84B FULL DADA2: ", project)
  message("============================================================")

  pwork <- file.path(
    PROJECT_ROOT,
    project,
    "05_work",
    "step84B_full"
  )

  filt_dir <- file.path(
    pwork,
    "filtered"
  )

  dir.create(
    filt_dir,
    recursive=TRUE,
    showWarnings=FALSE
  )

  complete_marker <- file.path(
    pwork,
    "_STEP84B_FULL_COMPLETE.ok"
  )

  if (file.exists(complete_marker)) {

    summary_path <- file.path(
      OUT,
      paste0(project,"_full_summary.csv")
    )

    if (file.exists(summary_path)) {
      message(project, ": existing COMPLETE marker found; skipping.")
      return(
        safe_csv(summary_path)
      )
    }
  }

  pairs <- build_pair_table(
    project,
    inventory,
    manifest
  )

  expected_n <- nrow(
    manifest |>
      filter(project==!!project)
  )

  if (nrow(pairs)!=expected_n) {
    stop(
      paste0(
        project,
        ": expected ",
        expected_n,
        " frozen Runs, but found ",
        nrow(pairs),
        " complete canonical R1/R2 pairs."
      )
    )
  }

  if (
    anyDuplicated(pairs$run_id) ||
    anyDuplicated(pairs$R1) ||
    anyDuplicated(pairs$R2)
  ) {
    stop(
      paste0(
        project,
        ": duplicate Run or paired-file mapping detected."
      )
    )
  }

  write_excel_csv(
    pairs,
    file.path(
      OUT,
      paste0(project,"_full_run_manifest.csv")
    ),
    na=""
  )

  fnFs <- pairs$R1
  fnRs <- pairs$R2
  names(fnFs) <- pairs$run_id
  names(fnRs) <- pairs$run_id

  filtFs <- file.path(
    filt_dir,
    paste0(pairs$run_id,"_R1_filt.fastq.gz")
  )

  filtRs <- file.path(
    filt_dir,
    paste0(pairs$run_id,"_R2_filt.fastq.gz")
  )

  names(filtFs) <- pairs$run_id
  names(filtRs) <- pairs$run_id

  # ----------------------------------------------------------
  # Parameters
  # ----------------------------------------------------------
  trimLeftF <- 0L
  trimLeftR <- 0L

  if (
    str_detect(
      param$primer_handling,
      "^CUTADAPT_REMOVE_R1_"
    )
  ) {
    trimLeftF <- primer_len(param$primer_R1)
    trimLeftR <- primer_len(param$primer_R2)
  }

  if (project=="PRJNA1166732") {
    trimLeftF <- primer_len("806R")
    trimLeftR <- primer_len("515F")
  }

  truncF <- as.integer(param$truncLen_R1)
  truncR <- as.integer(param$truncLen_R2)

  if (is.na(truncF)) truncF <- 0L
  if (is.na(truncR)) truncR <- 0L

  maxEEF <- as.numeric(param$maxEE_R1)
  maxEER <- as.numeric(param$maxEE_R2)

  if (is.na(maxEEF)) maxEEF <- 2
  if (is.na(maxEER)) maxEER <- 2

  minLenUse <- as.integer(param$minLen)
  if (is.na(minLenUse)) minLenUse <- 20L

  minOverlapUse <- if (project=="PRJNA1166732") {
    8L
  } else {
    12L
  }

  maxMismatchUse <- 0L

  param_record <- tibble(
    project=project,
    expected_runs=expected_n,
    trimLeft_R1=trimLeftF,
    trimLeft_R2=trimLeftR,
    truncLen_R1=truncF,
    truncLen_R2=truncR,
    maxEE_R1=maxEEF,
    maxEE_R2=maxEER,
    minLen=minLenUse,
    minOverlap=minOverlapUse,
    maxMismatch=maxMismatchUse
  )

  write_excel_csv(
    param_record,
    file.path(
      OUT,
      paste0(project,"_full_parameters.csv")
    ),
    na=""
  )

  # ----------------------------------------------------------
  # Clean only this project's prior incomplete Step84B outputs.
  # ----------------------------------------------------------
  safe_remove(
    c(
      list.files(
        filt_dir,
        full.names=TRUE,
        pattern="\\.fastq\\.gz$"
      ),
      file.path(pwork,paste0(project,"_full_errF.rds")),
      file.path(pwork,paste0(project,"_full_errR.rds")),
      file.path(pwork,paste0(project,"_seqtab_prechim_step84B.rds")),
      file.path(pwork,paste0(project,"_seqtab_nochim_step84B.rds"))
    )
  )

  # ----------------------------------------------------------
  # Filter
  # ----------------------------------------------------------
  ck(paste0(project," FILTER STARTED"))

  out_filter <- filterAndTrim(
    fnFs,
    filtFs,
    fnRs,
    filtRs,
    truncLen=c(truncF,truncR),
    trimLeft=c(trimLeftF,trimLeftR),
    maxN=0,
    maxEE=c(maxEEF,maxEER),
    truncQ=2,
    minLen=minLenUse,
    rm.phix=TRUE,
    compress=TRUE,
    multithread=FALSE,
    verbose=TRUE
  )

  out_df <- as.data.frame(out_filter)

  if (nrow(out_df)!=expected_n) {
    stop(
      paste0(
        project,
        ": filterAndTrim row count mismatch."
      )
    )
  }

  rownames(out_df) <- pairs$run_id

  if (!all(file.exists(filtFs)) || !all(file.exists(filtRs))) {
    stop(
      paste0(
        project,
        ": one or more filtered FASTQ outputs missing."
      )
    )
  }

  ck(paste0(project," FILTER COMPLETE"))

  # ----------------------------------------------------------
  # Learn cohort-specific errors from the FULL cohort
  # ----------------------------------------------------------
  ck(paste0(project," ERROR LEARNING STARTED"))

  errF <- learnErrors(
    filtFs,
    nbases=1e8,
    randomize=TRUE,
    multithread=FALSE,
    verbose=TRUE
  )

  errR <- learnErrors(
    filtRs,
    nbases=1e8,
    randomize=TRUE,
    multithread=FALSE,
    verbose=TRUE
  )

  errF_path <- file.path(
    pwork,
    paste0(project,"_full_errF.rds")
  )

  errR_path <- file.path(
    pwork,
    paste0(project,"_full_errR.rds")
  )

  saveRDS(errF,errF_path)
  saveRDS(errR,errR_path)

  ck(paste0(project," ERROR LEARNING COMPLETE"))

  # ----------------------------------------------------------
  # Derep / DADA
  # ----------------------------------------------------------
  ck(paste0(project," DENOISING STARTED"))

  derepFs <- derepFastq(
    filtFs,
    verbose=FALSE
  )

  derepRs <- derepFastq(
    filtRs,
    verbose=FALSE
  )

  names(derepFs) <- pairs$run_id
  names(derepRs) <- pairs$run_id

  dadaFs <- dada(
    derepFs,
    err=errF,
    multithread=FALSE,
    pool=FALSE
  )

  dadaRs <- dada(
    derepRs,
    err=errR,
    multithread=FALSE,
    pool=FALSE
  )

  names(dadaFs) <- pairs$run_id
  names(dadaRs) <- pairs$run_id

  ck(paste0(project," DENOISING COMPLETE"))

  # ----------------------------------------------------------
  # Merge
  # ----------------------------------------------------------
  ck(paste0(project," MERGING STARTED"))

  mergers <- mergePairs(
    dadaFs,
    derepFs,
    dadaRs,
    derepRs,
    minOverlap=minOverlapUse,
    maxMismatch=maxMismatchUse,
    verbose=FALSE
  )

  names(mergers) <- pairs$run_id

  ck(paste0(project," MERGING COMPLETE"))

  # ----------------------------------------------------------
  # Sequence table / chimera
  # ----------------------------------------------------------
  seqtab <- makeSequenceTable(
    mergers
  )

  if (
    !is.matrix(seqtab) ||
    nrow(seqtab)==0 ||
    ncol(seqtab)==0
  ) {
    stop(
      paste0(
        project,
        ": invalid or empty pre-chimera sequence table."
      )
    )
  }

  seqtab_nochim <- removeBimeraDenovo(
    seqtab,
    method="consensus",
    multithread=FALSE,
    verbose=TRUE
  )

  if (
    !is.matrix(seqtab_nochim) ||
    nrow(seqtab_nochim)==0 ||
    ncol(seqtab_nochim)==0
  ) {
    stop(
      paste0(
        project,
        ": invalid or empty non-chimeric sequence table."
      )
    )
  }

  seqtab_path <- file.path(
    pwork,
    paste0(project,"_seqtab_prechim_step84B.rds")
  )

  nochim_path <- file.path(
    pwork,
    paste0(project,"_seqtab_nochim_step84B.rds")
  )

  saveRDS(seqtab,seqtab_path)
  saveRDS(seqtab_nochim,nochim_path)

  ck(paste0(project," SEQTAB AND CHIMERA COMPLETE"))

  # ----------------------------------------------------------
  # Tracking
  # ----------------------------------------------------------
  denF <- sapply(dadaFs,getN)
  denR <- sapply(dadaRs,getN)
  merged <- sapply(mergers,getN)

  nonchim <- rowSums(seqtab_nochim)

  nonchim_aligned <- setNames(
    rep(0,length(pairs$run_id)),
    pairs$run_id
  )

  if (length(nonchim)) {
    nonchim_aligned[names(nonchim)] <- nonchim
  }

  track <- tibble(
    project=project,
    run_id=pairs$run_id,
    patient_id=pairs$patient_id,
    sample_id=pairs$sample_id,
    time_raw=pairs$time_raw,

    input=as.numeric(
      out_df[,1]
    ),

    filtered=as.numeric(
      out_df[,2]
    ),

    denoised_F=as.numeric(
      denF[pairs$run_id]
    ),

    denoised_R=as.numeric(
      denR[pairs$run_id]
    ),

    merged=as.numeric(
      merged[pairs$run_id]
    ),

    nonchim=as.numeric(
      nonchim_aligned[pairs$run_id]
    )
  ) |>
    mutate(
      filter_retention_pct=round(
        ifelse(input>0,100*filtered/input,NA_real_),
        2
      ),

      merge_retention_pct=round(
        ifelse(filtered>0,100*merged/filtered,NA_real_),
        2
      ),

      nonchim_of_merged_pct=round(
        ifelse(merged>0,100*nonchim/merged,NA_real_),
        2
      ),

      final_retention_pct=round(
        ifelse(input>0,100*nonchim/input,NA_real_),
        2
      )
    )

  track_path <- file.path(
    OUT,
    paste0(project,"_full_tracking.csv")
  )

  write_excel_csv(
    track,
    track_path,
    na=""
  )

  # ----------------------------------------------------------
  # ASV length QC
  # ----------------------------------------------------------
  lr <- length_range(project)

  lens <- nchar(
    colnames(seqtab_nochim)
  )

  abund <- colSums(
    seqtab_nochim
  )

  asv_len <- tibble(
    project=project,
    asv_sequence=colnames(seqtab_nochim),
    asv_length=lens,
    total_abundance=abund,
    within_expected_length=
      lens>=lr[1] &
      lens<=lr[2]
  )

  asv_len_path <- file.path(
    OUT,
    paste0(project,"_full_ASV_length_QC.csv")
  )

  write_excel_csv(
    asv_len,
    asv_len_path,
    na=""
  )

  # ----------------------------------------------------------
  # Full-cohort QC
  # ----------------------------------------------------------
  median_filter <- median(
    track$filter_retention_pct,
    na.rm=TRUE
  )

  median_merge <- median(
    track$merge_retention_pct,
    na.rm=TRUE
  )

  zero_nonchim <- sum(
    track$nonchim==0 |
    is.na(track$nonchim)
  )

  total_merged <- sum(
    track$merged,
    na.rm=TRUE
  )

  total_nonchim <- sum(
    track$nonchim,
    na.rm=TRUE
  )

  chimera_read_retention <- if (total_merged>0) {
    100*total_nonchim/total_merged
  } else {
    0
  }

  length_ok_pct <- round(
    100*mean(
      asv_len$within_expected_length,
      na.rm=TRUE
    ),
    2
  )

  filter_pass <- is.finite(median_filter) &&
    median_filter>=60

  merge_pass <- is.finite(median_merge) &&
    median_merge>=50

  zero_pass <- zero_nonchim==0

  chimera_pass <- is.finite(chimera_read_retention) &&
    chimera_read_retention>=50

  length_pass <- length_ok_pct>=80

  full_qc_pass <-
    filter_pass &&
    merge_pass &&
    zero_pass &&
    chimera_pass &&
    length_pass

  summary <- tibble(
    project=project,
    expected_runs=expected_n,
    processed_runs=nrow(track),

    trimLeft_R1=trimLeftF,
    trimLeft_R2=trimLeftR,
    truncLen_R1=truncF,
    truncLen_R2=truncR,
    maxEE_R1=maxEEF,
    maxEE_R2=maxEER,
    minOverlap=minOverlapUse,
    maxMismatch=maxMismatchUse,

    total_input_reads=sum(track$input,na.rm=TRUE),
    total_filtered_reads=sum(track$filtered,na.rm=TRUE),
    total_merged_reads=total_merged,
    total_nonchim_reads=total_nonchim,

    median_filter_retention_pct=round(
      median_filter,
      2
    ),

    median_merge_retention_pct=round(
      median_merge,
      2
    ),

    total_nonchim_of_merged_pct=round(
      chimera_read_retention,
      2
    ),

    zero_nonchim_samples=zero_nonchim,

    ASVs_prechim=ncol(seqtab),
    ASVs_nochim=ncol(seqtab_nochim),

    median_ASV_length=median(lens),

    expected_ASV_length_min=lr[1],
    expected_ASV_length_max=lr[2],

    ASVs_within_expected_length_pct=length_ok_pct,

    filter_pass=filter_pass,
    merge_pass=merge_pass,
    zero_sample_pass=zero_pass,
    chimera_read_retention_pass=chimera_pass,
    ASV_length_pass=length_pass,

    full_qc_pass=full_qc_pass,

    final_status=ifelse(
      full_qc_pass,
      "FULL_DADA2_PASS_PROMOTED_TO_SELECTED_ASV",
      "FULL_DADA2_COMPLETE_BUT_QC_REVIEW_REQUIRED"
    )
  )

  summary_path <- file.path(
    OUT,
    paste0(project,"_full_summary.csv")
  )

  write_excel_csv(
    summary,
    summary_path,
    na=""
  )

  # ----------------------------------------------------------
  # Promote only if all full-cohort QC passes.
  # ----------------------------------------------------------
  if (full_qc_pass) {

    promoted <- promote_final_asv(
      project,
      nochim_path,
      track_path,
      asv_len_path
    )

    writeLines(
      c(
        paste0("Project: ",project),
        paste0("Final ASV: ",promoted),
        paste0("Created: ",Sys.time()),
        "Status: FULL_DADA2_PASS_PROMOTED_TO_SELECTED_ASV"
      ),
      file.path(
        pwork,
        "_STEP84B_FULL_COMPLETE.ok"
      )
    )

  } else {

    writeLines(
      c(
        paste0("Project: ",project),
        paste0("Created: ",Sys.time()),
        "Status: FULL_DADA2_COMPLETE_BUT_QC_REVIEW_REQUIRED",
        "Final ASV was NOT promoted to 02_asv_selected."
      ),
      file.path(
        pwork,
        "_STEP84B_FULL_QC_REVIEW_REQUIRED.txt"
      )
    )
  }

  summary
}

# ============================================================
# Main
# ============================================================
main <- function() {

  ck("STEP84B STARTED")

  param_path <- file.path(
    FREEZE_ROOT,
    "V2_STEP83C_queue_DADA2_PILOT.csv"
  )

  inventory_path <- file.path(
    STEP83B_ROOT,
    "V2_STEP83B_canonical_raw_inventory.csv"
  )

  manifest_path <- file.path(
    FREEZE_ROOT,
    "V2_FINAL_FROZEN_16S_sequence_manifest.csv"
  )

  retry_decision_path <- file.path(
    STEP84A2_ROOT,
    "PRJNA1166732_FINAL_retry_decision.csv"
  )

  needed <- c(
    param_path,
    inventory_path,
    manifest_path,
    retry_decision_path
  )

  if (!all(file.exists(needed))) {
    stop(
      "Required Step83C/83B/84A2/frozen inputs are missing."
    )
  }

  params <- safe_csv(param_path)
  inventory <- safe_csv(inventory_path)
  manifest <- safe_csv(manifest_path)
  retry <- safe_csv(retry_decision_path)

  expected_projects <- c(
    "PRJEB82425",
    "PRJNA516701",
    "PRJNA578267",
    "PRJNA851469",
    "PRJNA1166732"
  )

  if (
    nrow(params)!=5 ||
    !setequal(params$project,expected_projects)
  ) {
    stop(
      "Step83C DADA2 pilot queue no longer matches the five approved projects."
    )
  }

  if (
    nrow(retry)!=1 ||
    retry$project[1]!="PRJNA1166732" ||
    retry$retry_status[1]!="PASS" ||
    retry$selected_minOverlap[1]!=8
  ) {
    stop(
      "PRJNA1166732 Step84A2 minOverlap=8 PASS guard failed."
    )
  }

  if (nrow(manifest)!=1941) {
    stop(
      "Frozen 16S manifest guard failed."
    )
  }

  # Optional one-project mode.
  args <- commandArgs(trailingOnly=TRUE)

  target_projects <- expected_projects

  if (length(args)>=1 && toupper(args[1])!="ALL") {

    if (!(args[1] %in% expected_projects)) {
      stop(
        paste0(
          "Unknown Step84B project argument: ",
          args[1]
        )
      )
    }

    target_projects <- args[1]
  }

  ck(
    paste0(
      "TARGET PROJECTS: ",
      paste(target_projects,collapse=",")
    )
  )

  # Software versions
  versions <- tibble(
    software=c("R","dada2"),
    version=c(
      R.version.string,
      as.character(packageVersion("dada2"))
    )
  )

  write_excel_csv(
    versions,
    file.path(
      OUT,
      "V2_STEP84B_software_versions.csv"
    ),
    na=""
  )

  # Project-independent processing with error isolation.
  completed <- list()
  errors <- list()

  for (project in target_projects) {

    ck(paste0(project," FULL STARTED"))

    param <- params |>
      filter(project==!!project) |>
      slice(1)

    z <- tryCatch(
      run_full_project(
        project,
        param,
        inventory,
        manifest
      ),
      error=function(e) {

        errors[[
          length(errors)+1
        ]] <<- tibble(
          project=project,
          error_message=conditionMessage(e),
          error_call=paste(
            deparse(conditionCall(e)),
            collapse=" "
          )
        )

        ck(paste0(project," FULL FAILED"))

        NULL
      }
    )

    if (!is.null(z)) {
      completed[[length(completed)+1]] <- z
      ck(paste0(project," FULL COMPLETE"))
    }

    # Try to release memory between cohorts.
    gc(verbose=FALSE)
  }

  combined <- if (length(completed)) {
    bind_rows(completed)
  } else {
    tibble()
  }

  error_df <- if (length(errors)) {
    bind_rows(errors)
  } else {
    tibble(
      project=character(),
      error_message=character(),
      error_call=character()
    )
  }

  write_excel_csv(
    combined,
    file.path(
      OUT,
      "V2_STEP84B_FULL_DADA2_SUMMARY.csv"
    ),
    na=""
  )

  write_excel_csv(
    error_df,
    file.path(
      OUT,
      "V2_STEP84B_project_errors.csv"
    ),
    na=""
  )

  pass_queue <- if (nrow(combined)) {
    combined |>
      filter(full_qc_pass)
  } else {
    tibble()
  }

  review_queue <- if (nrow(combined)) {
    combined |>
      filter(!full_qc_pass)
  } else {
    tibble()
  }

  write_excel_csv(
    pass_queue,
    file.path(
      OUT,
      "V2_STEP84B_queue_TAXONOMY_READY.csv"
    ),
    na=""
  )

  write_excel_csv(
    review_queue,
    file.path(
      OUT,
      "V2_STEP84B_queue_QC_REVIEW.csv"
    ),
    na=""
  )

  global <- tibble(
    metric=c(
      "Target_projects",
      "Completed_projects",
      "Runtime_failed_projects",
      "Full_QC_pass_projects",
      "Full_QC_review_projects"
    ),
    value=c(
      length(target_projects),
      nrow(combined),
      nrow(error_df),
      if (nrow(combined)) sum(combined$full_qc_pass) else 0,
      if (nrow(combined)) sum(!combined$full_qc_pass) else 0
    )
  )

  write_excel_csv(
    global,
    file.path(
      OUT,
      "V2_STEP84B_global_summary.csv"
    ),
    na=""
  )

  readme <- c(
    "SEPSIS V2 - STEP84B FULL-COHORT DADA2",
    paste0("Created: ",Sys.time()),
    "",
    "INPUT ROOT",
    CANON,
    "",
    "PROCESSED COHORTS",
    paste(expected_projects,collapse=", "),
    "",
    "PRJNA1166732",
    "Step84A2 demonstrated that minOverlap 12 and 10 produced zero accepted merged reads, whereas minOverlap 8 and 6 passed.",
    "The most conservative passing threshold, minOverlap=8 with maxMismatch=0, is used here.",
    "",
    "FULL-COHORT QC REQUIRED FOR PROMOTION",
    "- median filter retention >=60%",
    "- median merge retention >=50%",
    "- zero final non-chimeric samples",
    "- >=50% merged reads remain after chimera removal",
    "- >=80% ASVs within broad project-specific expected length range",
    "",
    "DATA SAFETY",
    "Full outputs are first written under each project's 05_work/step84B_full.",
    "Only cohorts passing all full QC are hard-linked into 02_asv_selected.",
    "",
    "RESUME",
    "A completed/passed project gets _STEP84B_FULL_COMPLETE.ok and will be skipped when the script is rerun.",
    "",
    "NEXT",
    "Step85 will assign taxonomy to all selected final ASV tables using one documented reference database/version, while keeping cohorts separate at ASV level."
  )

  writeLines(
    readme,
    file.path(
      OUT,
      "README_STEP84B_FULL_DADA2.txt"
    ),
    useBytes=TRUE
  )

  ck("STEP84B COMPLETE")

  cat("\n============================================================\n")
  cat("SEPSIS V2 - STEP84B FULL DADA2 COMPLETE\n")
  cat("============================================================\n\n")

  print(global,n=Inf,width=Inf)

  if (nrow(combined)) {
    cat("\nFULL PROJECT RESULTS:\n")
    print(combined,n=Inf,width=Inf)
  }

  if (nrow(error_df)) {
    cat("\nRUNTIME ERRORS:\n")
    print(error_df,n=Inf,width=Inf)
  }

  cat("\nOutput:\n")
  cat(OUT,"\n")
  cat("============================================================\n")
}

tryCatch(
  main(),
  error=function(e) {

    msg <- c(
      paste0(
        "STEP84B FATAL ERROR: ",
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

    ck("STEP84B FAILED")

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
