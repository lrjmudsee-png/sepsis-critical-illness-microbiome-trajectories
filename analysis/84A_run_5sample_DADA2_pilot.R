# ============================================================
# Sepsis V2 - Step 84A
# 5-sample-per-cohort DADA2 PILOT
#
# STRICT INPUT:
#   E:/sepsis_project/data/_V2_ANALYSIS_READY
#
# Processes ONLY cohorts marked READY_FOR_PILOT in Step83C.
# DOES NOT process HOLD / ASV-reuse / FASTA-special cohorts.
#
# Per project:
#   1) select 5 representative paired Runs
#   2) trim fixed primer lengths when Step83C requires it
#   3) filterAndTrim
#   4) learnErrors
#   5) dereplicate
#   6) dada
#   7) mergePairs
#   8) makeSequenceTable
#   9) removeBimeraDenovo
#  10) evaluate filtering, merging, chimera retention,
#      zero-output samples and ASV-length plausibility
#
# Pilot outputs go to:
#   E:/sepsis_project/data/_V2_ANALYSIS_READY/
#     01_PROJECTS/<PROJECT>/05_work/step84A_pilot/
#
# Summaries go to:
#   E:/sepsis_project/results/V2_24A_DADA2_PILOT/
#
# R 4.4.0 / Windows
# ============================================================

options(stringsAsFactors = FALSE)

# ------------------------------------------------------------
# Package setup
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

OUT <- file.path(
  ROOT,
  "results",
  "V2_24A_DADA2_PILOT"
)

dir.create(OUT, recursive=TRUE, showWarnings=FALSE)

LOG <- file.path(OUT, "_STEP84A_runtime_checkpoints.txt")
ERR <- file.path(OUT, "_STEP84A_FATAL_ERROR.txt")
if (file.exists(ERR)) unlink(ERR)

ck <- function(x) {
  cat(paste0(x, ": ", Sys.time(), "\n"), file=LOG, append=TRUE)
}

clean_chr <- function(x) {
  x <- trimws(as.character(x))
  x[
    is.na(x) |
    x=="" |
    tolower(x) %in% c("na","nan","n/a","null","none")
  ] <- NA_character_
  x
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

pick_evenly <- function(x, n=5) {
  x <- unique(x)
  if (length(x)<=n) return(x)
  idx <- unique(round(seq(1,length(x),length.out=n)))
  x[idx]
}

getN <- function(x) {
  sum(getUniques(x))
}

primer_len <- function(x) {

  if (is.na(x) || x=="") return(0L)

  # Fixed lengths of primer labels frozen in Step83C.
  lens <- c(
    "341F"=17L,
    "785R"=21L,
    "806R"=20L,
    "515F"=19L
  )

  if (x %in% names(lens)) {
    return(unname(lens[[x]]))
  }

  0L
}

# ------------------------------------------------------------
# Plausible post-primer merged ASV length ranges
# Broad enough for pilot QC, not taxonomy inference.
# ------------------------------------------------------------
length_range <- function(project) {

  switch(
    project,

    "PRJEB82425"  = c(350L,500L),
    "PRJNA578267" = c(350L,500L),
    "PRJNA851469" = c(200L,330L),
    "PRJNA1166732"= c(200L,330L),
    "PRJNA516701" = c(150L,500L),

    c(100L,600L)
  )
}

# ------------------------------------------------------------
# Match canonical R1/R2 files by Run accession.
# Uses ONLY Step83B canonical inventory.
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
    select(run_id, R1=file_path)

  r2 <- z |>
    filter(mate=="R2") |>
    select(run_id, R2=file_path)

  pairs <- inner_join(
    r1,
    r2,
    by="run_id"
  ) |>
    distinct(run_id,.keep_all=TRUE) |>
    inner_join(
      manifest |>
        select(
          project,
          run_id,
          patient_id,
          sample_id,
          time_raw,
          time_day
        ),
      by="run_id"
    )

  pairs
}

# ------------------------------------------------------------
# One-project pilot
# ------------------------------------------------------------
run_project_pilot <- function(project, param, inventory, manifest) {

  message("\n==============================")
  message("STEP84A PILOT: ",project)
  message("==============================")

  pwork <- file.path(
    PROJECT_ROOT,
    project,
    "05_work",
    "step84A_pilot"
  )

  filt_dir <- file.path(pwork,"filtered")
  dir.create(filt_dir, recursive=TRUE, showWarnings=FALSE)

  # Clean only previous Step84A pilot filtered files for reproducibility.
  old <- list.files(
    filt_dir,
    full.names=TRUE,
    pattern="\\.fastq\\.gz$"
  )
  if (length(old)) unlink(old)

  pairs <- build_pair_table(
    project,
    inventory,
    manifest
  )

  if (nrow(pairs)<5) {
    stop(
      paste0(
        project,
        ": fewer than 5 complete canonical R1/R2 pairs (",
        nrow(pairs),
        ")."
      )
    )
  }

  # Sort using time and Run, then take evenly distributed Runs.
  pairs <- pairs |>
    arrange(
      patient_id,
      time_day,
      time_raw,
      run_id
    )

  chosen_ids <- pick_evenly(
    pairs$run_id,
    n=5
  )

  pilot <- pairs |>
    filter(run_id %in% chosen_ids) |>
    arrange(match(run_id,chosen_ids))

  if (nrow(pilot)!=5) {
    stop(
      paste0(
        project,
        ": pilot selection did not yield exactly 5 Runs."
      )
    )
  }

  write_excel_csv(
    pilot,
    file.path(
      OUT,
      paste0(project,"_pilot_selected_runs.csv")
    ),
    na=""
  )

  fnFs <- pilot$R1
  fnRs <- pilot$R2
  names(fnFs) <- pilot$run_id
  names(fnRs) <- pilot$run_id

  filtFs <- file.path(
    filt_dir,
    paste0(pilot$run_id,"_R1_filt.fastq.gz")
  )

  filtRs <- file.path(
    filt_dir,
    paste0(pilot$run_id,"_R2_filt.fastq.gz")
  )

  names(filtFs) <- pilot$run_id
  names(filtRs) <- pilot$run_id

  # ----------------------------------------------------------
  # Primer removal via fixed trimLeft lengths.
  # Step83B detected these primers at the read starts.
  # For cohorts classified as already trimmed, trimLeft=0.
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

  if (
    param$project=="PRJNA1166732"
  ) {
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

  # ----------------------------------------------------------
  # Filter and trim
  # ----------------------------------------------------------
  out <- filterAndTrim(
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

  out_df <- as.data.frame(out)

  # Ensure order is tied to chosen Run IDs.
  if (nrow(out_df)!=5) {
    stop(paste0(project,": filterAndTrim did not return 5 rows."))
  }

  rownames(out_df) <- pilot$run_id

  # Check that filtered outputs actually exist.
  keep <- file.exists(filtFs) & file.exists(filtRs)

  if (!all(keep)) {
    stop(
      paste0(
        project,
        ": one or more filtered pilot files were not created."
      )
    )
  }

  # ----------------------------------------------------------
  # Learn errors from pilot reads
  # ----------------------------------------------------------
  errF <- learnErrors(
    filtFs,
    nbases=5e6,
    randomize=TRUE,
    multithread=FALSE,
    verbose=TRUE
  )

  errR <- learnErrors(
    filtRs,
    nbases=5e6,
    randomize=TRUE,
    multithread=FALSE,
    verbose=TRUE
  )

  saveRDS(
    errF,
    file.path(pwork,paste0(project,"_pilot_errF.rds"))
  )
  saveRDS(
    errR,
    file.path(pwork,paste0(project,"_pilot_errR.rds"))
  )

  # ----------------------------------------------------------
  # Dereplicate and infer ASVs
  # ----------------------------------------------------------
  derepFs <- derepFastq(
    filtFs,
    verbose=FALSE
  )
  derepRs <- derepFastq(
    filtRs,
    verbose=FALSE
  )

  names(derepFs) <- pilot$run_id
  names(derepRs) <- pilot$run_id

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

  # ----------------------------------------------------------
  # Merge pairs
  # ----------------------------------------------------------
  mergers <- mergePairs(
    dadaFs,
    derepFs,
    dadaRs,
    derepRs,
    minOverlap=12,
    maxMismatch=0,
    verbose=TRUE
  )

  seqtab <- makeSequenceTable(
    mergers
  )

  seqtab_nochim <- removeBimeraDenovo(
    seqtab,
    method="consensus",
    multithread=FALSE,
    verbose=TRUE
  )

  saveRDS(
    seqtab,
    file.path(
      pwork,
      paste0(project,"_pilot_seqtab_prechim.rds")
    )
  )

  saveRDS(
    seqtab_nochim,
    file.path(
      pwork,
      paste0(project,"_pilot_seqtab_nochim.rds")
    )
  )

  # ----------------------------------------------------------
  # Tracking
  # ----------------------------------------------------------
  denF <- sapply(dadaFs,getN)
  denR <- sapply(dadaRs,getN)
  merged <- sapply(mergers,getN)

  nochim <- rowSums(seqtab_nochim)

  # makeSequenceTable retains named samples; align explicitly.
  nochim2 <- setNames(
    rep(0,length(pilot$run_id)),
    pilot$run_id
  )

  if (length(nochim)) {
    nochim2[names(nochim)] <- nochim
  }

  track <- tibble(
    project=project,
    run_id=pilot$run_id,
    patient_id=pilot$patient_id,
    sample_id=pilot$sample_id,
    time_raw=pilot$time_raw,

    input=as.numeric(out_df[,1]),
    filtered=as.numeric(out_df[,2]),

    denoised_F=as.numeric(
      denF[pilot$run_id]
    ),

    denoised_R=as.numeric(
      denR[pilot$run_id]
    ),

    merged=as.numeric(
      merged[pilot$run_id]
    ),

    nonchim=as.numeric(
      nochim2[pilot$run_id]
    )
  ) |>
    mutate(
      filter_retention_pct=round(
        100*filtered/input,
        2
      ),

      merge_retention_pct=round(
        ifelse(
          filtered>0,
          100*merged/filtered,
          NA_real_
        ),
        2
      ),

      nonchim_of_merged_pct=round(
        ifelse(
          merged>0,
          100*nonchim/merged,
          NA_real_
        ),
        2
      ),

      final_retention_pct=round(
        ifelse(
          input>0,
          100*nonchim/input,
          NA_real_
        ),
        2
      )
    )

  write_excel_csv(
    track,
    file.path(
      OUT,
      paste0(project,"_pilot_tracking.csv")
    ),
    na=""
  )

  # ----------------------------------------------------------
  # ASV length QC
  # ----------------------------------------------------------
  lens <- nchar(colnames(seqtab_nochim))
  abund <- if (ncol(seqtab_nochim)) colSums(seqtab_nochim) else numeric()

  lr <- length_range(project)

  asv_len <- tibble(
    project=project,
    asv_sequence=colnames(seqtab_nochim),
    asv_length=lens,
    total_abundance=abund,
    within_expected_length=
      lens>=lr[1] &
      lens<=lr[2]
  )

  write_excel_csv(
    asv_len,
    file.path(
      OUT,
      paste0(project,"_pilot_ASV_length_QC.csv")
    ),
    na=""
  )

  # ----------------------------------------------------------
  # Pilot pass/fail rules
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

  n_asv <- ncol(seqtab_nochim)

  length_ok_pct <- if (n_asv>0) {
    round(
      100*mean(
        asv_len$within_expected_length,
        na.rm=TRUE
      ),
      2
    )
  } else {
    0
  }

  median_len <- if (n_asv>0) {
    median(lens)
  } else {
    NA_real_
  }

  filter_pass <- is.finite(median_filter) &&
    median_filter>=60

  merge_pass <- is.finite(median_merge) &&
    median_merge>=50

  zero_pass <- zero_nonchim==0

  length_pass <- n_asv>0 &&
    length_ok_pct>=80

  overall_pass <-
    filter_pass &&
    merge_pass &&
    zero_pass &&
    length_pass

  project_summary <- tibble(
    project=project,

    pilot_runs=5L,

    trimLeft_R1=trimLeftF,
    trimLeft_R2=trimLeftR,

    truncLen_R1=truncF,
    truncLen_R2=truncR,

    maxEE_R1=maxEEF,
    maxEE_R2=maxEER,

    input_reads=sum(
      track$input,
      na.rm=TRUE
    ),

    filtered_reads=sum(
      track$filtered,
      na.rm=TRUE
    ),

    merged_reads=sum(
      track$merged,
      na.rm=TRUE
    ),

    nonchim_reads=sum(
      track$nonchim,
      na.rm=TRUE
    ),

    median_filter_retention_pct=round(
      median_filter,
      2
    ),

    median_merge_retention_pct=round(
      median_merge,
      2
    ),

    zero_nonchim_samples=zero_nonchim,

    seqtab_ASVs_prechim=ncol(seqtab),
    seqtab_ASVs_nochim=n_asv,

    median_ASV_length=median_len,

    expected_ASV_length_min=lr[1],
    expected_ASV_length_max=lr[2],

    ASVs_within_expected_length_pct=
      length_ok_pct,

    filter_pass=filter_pass,
    merge_pass=merge_pass,
    zero_sample_pass=zero_pass,
    length_pass=length_pass,

    pilot_pass=overall_pass,

    recommended_next_action=case_when(
      overall_pass ~
        "ALLOW_FULL_COHORT_DADA2_WITH_THESE_PARAMETERS",

      !filter_pass ~
        "ADJUST_FILTER_OR_TRUNCATION_PARAMETERS_AND_REPEAT_PILOT",

      filter_pass & !merge_pass ~
        "ADJUST_TRUNCLEN_OR_PRIMER_HANDLING_TO_IMPROVE_OVERLAP_AND_REPEAT_PILOT",

      filter_pass & merge_pass & !zero_pass ~
        "INVESTIGATE_ZERO_OUTPUT_SAMPLE_AND_REPEAT_PILOT",

      filter_pass & merge_pass & zero_pass & !length_pass ~
        "CHECK_PRIMER_REGION_OR_MERGE_LENGTH_DISTRIBUTION_BEFORE_FULL_RUN",

      TRUE ~
        "MANUAL_REVIEW"
    )
  )

  write_excel_csv(
    project_summary,
    file.path(
      OUT,
      paste0(project,"_pilot_summary.csv")
    ),
    na=""
  )

  project_summary
}

# ============================================================
# Main
# ============================================================
main <- function() {

  ck("STEP84A STARTED")

  params_path <- file.path(
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

  needed <- c(
    params_path,
    inventory_path,
    manifest_path
  )

  if (!all(file.exists(needed))) {
    stop(
      "Required Step83C/83B/frozen manifest inputs are missing."
    )
  }

  params <- safe_csv(params_path)
  inventory <- safe_csv(inventory_path)
  manifest <- safe_csv(manifest_path)

  if (nrow(params)!=5) {
    stop(
      paste0(
        "Expected 5 READY_FOR_PILOT projects; observed ",
        nrow(params)
      )
    )
  }

  if (nrow(manifest)!=1941) {
    stop("Frozen 16S manifest guard failed.")
  }

  ck("INPUTS GUARDED")

  # ----------------------------------------------------------
  # Version record
  # ----------------------------------------------------------
  versions <- tibble(
    software=c(
      "R",
      "dada2"
    ),
    version=c(
      R.version.string,
      as.character(
        packageVersion("dada2")
      )
    )
  )

  write_excel_csv(
    versions,
    file.path(
      OUT,
      "V2_STEP84A_software_versions.csv"
    ),
    na=""
  )

  # ----------------------------------------------------------
  # Run each project independently so one failed pilot does not
  # destroy results from already completed cohorts.
  # ----------------------------------------------------------
  project_results <- list()
  project_errors <- list()

  for (project in params$project) {

    ck(paste0(project," PILOT STARTED"))

    param <- params |>
      filter(project==!!project) |>
      slice(1)

    result <- tryCatch(
      run_project_pilot(
        project,
        param,
        inventory,
        manifest
      ),
      error=function(e) {

        project_errors[[
          length(project_errors)+1
        ]] <<- tibble(
          project=project,
          error_message=conditionMessage(e),
          error_call=paste(
            deparse(conditionCall(e)),
            collapse=" "
          )
        )

        ck(paste0(project," PILOT FAILED"))

        NULL
      }
    )

    if (!is.null(result)) {
      project_results[[
        length(project_results)+1
      ]] <- result

      ck(paste0(project," PILOT COMPLETE"))
    }
  }

  combined <- if (length(project_results)) {
    bind_rows(project_results)
  } else {
    tibble()
  }

  errors <- if (length(project_errors)) {
    bind_rows(project_errors)
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
      "V2_STEP84A_DADA2_PILOT_SUMMARY.csv"
    ),
    na=""
  )

  write_excel_csv(
    errors,
    file.path(
      OUT,
      "V2_STEP84A_project_errors.csv"
    ),
    na=""
  )

  # ----------------------------------------------------------
  # Final decision queue
  # ----------------------------------------------------------
  if (nrow(combined)) {

    pass_queue <- combined |>
      filter(pilot_pass)

    repeat_queue <- combined |>
      filter(!pilot_pass)

  } else {

    pass_queue <- tibble()
    repeat_queue <- tibble()
  }

  write_excel_csv(
    pass_queue,
    file.path(
      OUT,
      "V2_STEP84A_queue_FULL_DADA2_ALLOWED.csv"
    ),
    na=""
  )

  write_excel_csv(
    repeat_queue,
    file.path(
      OUT,
      "V2_STEP84A_queue_REPEAT_PILOT.csv"
    ),
    na=""
  )

  # ----------------------------------------------------------
  # Summary
  # ----------------------------------------------------------
  summary <- tibble(
    metric=c(
      "Projects_requested",
      "Projects_completed",
      "Projects_failed_runtime",
      "Projects_pilot_pass",
      "Projects_repeat_pilot"
    ),

    value=c(
      nrow(params),
      nrow(combined),
      nrow(errors),
      if (nrow(combined)) sum(combined$pilot_pass) else 0,
      if (nrow(combined)) sum(!combined$pilot_pass) else 0
    )
  )

  write_excel_csv(
    summary,
    file.path(
      OUT,
      "V2_STEP84A_global_summary.csv"
    ),
    na=""
  )

  # ----------------------------------------------------------
  # README
  # ----------------------------------------------------------
  readme <- c(
    "SEPSIS V2 - STEP84A 5-SAMPLE DADA2 PILOT",
    paste0("Created: ",Sys.time()),
    "",
    "STRICT INPUT",
    CANON,
    "",
    "ONLY Step83C READY_FOR_PILOT cohorts are processed.",
    "",
    "PIPELINE",
    "filterAndTrim -> learnErrors -> derepFastq -> dada -> mergePairs -> makeSequenceTable -> removeBimeraDenovo",
    "",
    "WINDOWS",
    "multithread=FALSE is intentionally used.",
    "",
    "PILOT PASS RULES",
    "- median filtering retention >= 60%",
    "- median merge retention >= 50%",
    "- zero samples with no non-chimeric reads",
    "- >=80% of inferred non-chimeric ASVs in broad region-appropriate length range",
    "",
    "IMPORTANT",
    "Primer-positive cohorts use fixed trimLeft equal to the observed primer length during this pilot.",
    "This is appropriate only because Step83B found the frozen primer directly at the start of reads.",
    "No HOLD, PacBio, Ion Torrent anomaly, FASTA-only, static-deferred, or existing-ASV cohort is processed here.",
    "",
    "NEXT",
    "Projects in V2_STEP84A_queue_FULL_DADA2_ALLOWED.csv may proceed to full-cohort DADA2.",
    "Projects in V2_STEP84A_queue_REPEAT_PILOT.csv must have parameters adjusted before full processing."
  )

  writeLines(
    readme,
    file.path(
      OUT,
      "README_STEP84A_DADA2_PILOT.txt"
    ),
    useBytes=TRUE
  )

  ck("STEP84A COMPLETE")

  cat("\n============================================================\n")
  cat("SEPSIS V2 - STEP84A DADA2 PILOT COMPLETE\n")
  cat("============================================================\n\n")

  print(summary,n=Inf,width=Inf)

  if (nrow(combined)) {
    cat("\nPILOT RESULTS:\n")
    print(combined,n=Inf,width=Inf)
  }

  if (nrow(errors)) {
    cat("\nRUNTIME ERRORS:\n")
    print(errors,n=Inf,width=Inf)
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
        "STEP84A FATAL ERROR: ",
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

    ck("STEP84A FAILED")

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
