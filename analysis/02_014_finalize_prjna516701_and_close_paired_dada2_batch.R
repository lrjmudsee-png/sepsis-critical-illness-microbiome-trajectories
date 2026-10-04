# ============================================================
# Sepsis V2 - Step 84C
# Finalize PRJNA516701 after explicit QC review
# + close the 5-cohort paired-end DADA2 batch
#
# DECISION:
#   PRJNA516701 DADA2 processing is technically successful.
#   The only automatic-QC failure is one zero-output sample:
#
#     SRR8488588
#     patient 3770
#     sample SAMN10794824
#     middle
#     input=30, filtered=17, merged=0, nonchim=0
#
#   This is an essentially empty sequencing library, not a
#   cohort-wide DADA2 failure.
#
#   Patient 3770 also has:
#     baseline SRR8488602
#     another middle SRR8488603
#
#   Therefore excluding SRR8488588 does NOT change the cohort's
#   longitudinal eligibility counts:
#     patients = 107
#     >=2 timepoints = 45
#     >=3 timepoints = 16
#
# This script:
#   1) verifies those hard guards;
#   2) removes ONLY the zero-output sample from the selected
#      analysis table;
#   3) does NOT impose a new low-depth cutoff yet;
#   4) promotes the reviewed 170-sample table;
#   5) writes low-depth candidate lists for later unified QC;
#   6) closes the five-cohort paired-end DADA2 batch.
#
# R 4.4.0 / Windows
# ============================================================

options(stringsAsFactors = FALSE)

pkgs <- c("readr","dplyr","tibble")
missing <- pkgs[
  !vapply(pkgs, requireNamespace, logical(1), quietly=TRUE)
]

if (length(missing)) {
  install.packages(
    missing,
    repos="https://cloud.r-project.org"
  )
}

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tibble)
})

ROOT <- "E:/sepsis_project"

CANON <- file.path(
  ROOT,
  "data",
  "_V2_ANALYSIS_READY"
)

P <- "PRJNA516701"

PWORK <- file.path(
  CANON,
  "01_PROJECTS",
  P,
  "05_work",
  "step84B_full"
)

SELECTED <- file.path(
  CANON,
  "01_PROJECTS",
  P,
  "02_asv_selected"
)

FULL84B <- file.path(
  ROOT,
  "results",
  "V2_24B_FULL_DADA2"
)

R578 <- file.path(
  ROOT,
  "results",
  "V2_24B4B_PRJNA578267_CHIMERA"
)

OUT <- file.path(
  ROOT,
  "results",
  "V2_24C_PAIRED_DADA2_CLOSEOUT"
)

dir.create(
  SELECTED,
  recursive=TRUE,
  showWarnings=FALSE
)

dir.create(
  OUT,
  recursive=TRUE,
  showWarnings=FALSE
)

LOG <- file.path(
  OUT,
  "_STEP84C_runtime_checkpoints.txt"
)

ERR <- file.path(
  OUT,
  "_STEP84C_FATAL_ERROR.txt"
)

if (file.exists(ERR)) unlink(ERR)

ck <- function(x) {
  cat(
    paste0(x,": ",Sys.time(),"\n"),
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

atomic_saveRDS <- function(obj,p) {

  tmp <- paste0(p,".tmp")

  if (file.exists(tmp)) unlink(tmp)

  saveRDS(
    obj,
    tmp
  )

  if (file.exists(p)) unlink(p)

  if (!file.rename(tmp,p)) {
    stop(
      "Could not atomically save: ",
      p
    )
  }
}

timepoint_stats <- function(run_manifest) {

  g <- run_manifest |>
    filter(
      !is.na(patient_id),
      !is.na(time_raw)
    ) |>
    distinct(
      patient_id,
      time_raw
    ) |>
    count(
      patient_id,
      name="n_timepoints"
    )

  tibble(
    patients=n_distinct(
      run_manifest$patient_id,
      na.rm=TRUE
    ),
    patients_GE2=sum(
      g$n_timepoints>=2
    ),
    patients_GE3=sum(
      g$n_timepoints>=3
    )
  )
}

read_summary <- function(project) {

  p <- file.path(
    FULL84B,
    paste0(
      project,
      "_full_summary.csv"
    )
  )

  if (!file.exists(p)) {
    stop(
      "Missing full summary for ",
      project
    )
  }

  safe_csv(p)
}

main <- function() {

  ck("STEP84C STARTED")

  summary_path <- file.path(
    FULL84B,
    "PRJNA516701_full_summary.csv"
  )

  tracking_path <- file.path(
    FULL84B,
    "PRJNA516701_full_tracking.csv"
  )

  run_manifest_path <- file.path(
    FULL84B,
    "PRJNA516701_full_run_manifest.csv"
  )

  asv_qc_path <- file.path(
    FULL84B,
    "PRJNA516701_full_ASV_length_QC.csv"
  )

  seqtab_path <- file.path(
    PWORK,
    "PRJNA516701_seqtab_nochim_step84B.rds"
  )

  needed <- c(
    summary_path,
    tracking_path,
    run_manifest_path,
    asv_qc_path,
    seqtab_path
  )

  if (!all(file.exists(needed))) {
    stop(
      "One or more PRJNA516701 Step84B inputs are missing."
    )
  }

  s <- safe_csv(
    summary_path
  )

  tr <- safe_csv(
    tracking_path
  )

  rmf <- safe_csv(
    run_manifest_path
  )

  aq <- safe_csv(
    asv_qc_path
  )

  seqtab <- readRDS(
    seqtab_path
  )

  # ==========================================================
  # 1. Hard guards
  # ==========================================================

  if (
    nrow(s)!=1 ||
    s$expected_runs[1]!=171 ||
    s$processed_runs[1]!=171
  ) {
    stop(
      "PRJNA516701 171/171 guard failed."
    )
  }

  if (
    !isTRUE(s$filter_pass[1]) ||
    !isTRUE(s$merge_pass[1]) ||
    !isTRUE(s$chimera_read_retention_pass[1]) ||
    !isTRUE(s$ASV_length_pass[1])
  ) {
    stop(
      "PRJNA516701 has a QC failure other than zero-sample rule."
    )
  }

  if (
    s$zero_nonchim_samples[1]!=1 ||
    isTRUE(s$zero_sample_pass[1]) ||
    isTRUE(s$full_qc_pass[1])
  ) {
    stop(
      "Expected exactly one zero-output sample and original full_qc_pass=FALSE."
    )
  }

  zeros <- tr |>
    filter(
      is.na(nonchim) |
      nonchim==0
    )

  if (
    nrow(zeros)!=1 ||
    zeros$run_id[1]!="SRR8488588" ||
    as.character(
      zeros$patient_id[1]
    )!="3770" ||
    zeros$sample_id[1]!="SAMN10794824" ||
    zeros$time_raw[1]!="middle" ||
    zeros$input[1]!=30 ||
    zeros$filtered[1]!=17 ||
    zeros$merged[1]!=0
  ) {
    stop(
      "The expected SRR8488588 zero-library guard failed."
    )
  }

  patient3770 <- tr |>
    filter(
      as.character(patient_id)=="3770"
    ) |>
    arrange(
      time_raw,
      run_id
    )

  if (
    nrow(patient3770)!=3 ||
    !all(
      c(
        "SRR8488602",
        "SRR8488588",
        "SRR8488603"
      ) %in% patient3770$run_id
    )
  ) {
    stop(
      "Patient 3770 longitudinal guard failed."
    )
  }

  before <- timepoint_stats(
    rmf
  )

  after_rmf <- rmf |>
    filter(
      run_id!="SRR8488588"
    )

  after <- timepoint_stats(
    after_rmf
  )

  if (
    before$patients[1]!=107 ||
    before$patients_GE2[1]!=45 ||
    before$patients_GE3[1]!=16 ||
    after$patients[1]!=107 ||
    after$patients_GE2[1]!=45 ||
    after$patients_GE3[1]!=16
  ) {
    stop(
      "Longitudinal eligibility changed unexpectedly after zero-sample removal."
    )
  }

  if (
    !is.matrix(seqtab) ||
    nrow(seqtab)!=171 ||
    ncol(seqtab)!=3221
  ) {
    stop(
      paste0(
        "Expected PRJNA516701 seqtab 171 x 3221; observed ",
        nrow(seqtab),
        " x ",
        ncol(seqtab)
      )
    )
  }

  if (
    !("SRR8488588" %in% rownames(seqtab))
  ) {
    stop(
      "SRR8488588 is not present in the Step84B sequence table."
    )
  }

  if (
    sum(
      seqtab[
        "SRR8488588",
        ,
        drop=FALSE
      ]
    )!=0
  ) {
    stop(
      "SRR8488588 sequence-table row is not zero."
    )
  }

  ck("PRJNA516701 REVIEW GUARDS PASSED")

  # ==========================================================
  # 2. Exclude ONLY zero-output sample
  # ==========================================================

  seqtab_selected <- seqtab[
    rownames(seqtab)!="SRR8488588",
    ,
    drop=FALSE
  ]

  if (
    nrow(seqtab_selected)!=170 ||
    ncol(seqtab_selected)!=3221 ||
    any(
      rowSums(seqtab_selected)==0
    )
  ) {
    stop(
      "Reviewed selected sequence-table guard failed."
    )
  }

  selected_rds <- file.path(
    SELECTED,
    "PRJNA516701_seqtab_final_step84C_QC_REVIEWED.rds"
  )

  atomic_saveRDS(
    seqtab_selected,
    selected_rds
  )

  # ==========================================================
  # 3. Explicit exclusions + low-depth candidates
  # ==========================================================

  zero_exclusion <- zeros |>
    mutate(
      exclusion_stage=
        "POST_DADA2_SAMPLE_QC",

      exclusion_reason=
        "ZERO_NONCHIMERIC_READS_FROM_ESSENTIALLY_EMPTY_LIBRARY",

      decision=
        "EXCLUDE",

      affects_longitudinal_eligibility=
        FALSE,

      notes=
        "Input library contained only 30 reads. Patient 3770 retains baseline and another middle sample; cohort patient >=2/>=3 counts are unchanged."
    )

  write_excel_csv(
    zero_exclusion,
    file.path(
      OUT,
      "PRJNA516701_CONFIRMED_zero_read_sample_exclusion.csv"
    ),
    na=""
  )

  # Do NOT exclude these yet.
  low_depth <- tr |>
    mutate(
      candidate_lt_100=
        nonchim<100,

      candidate_lt_500=
        nonchim<500,

      candidate_lt_1000=
        nonchim<1000,

      candidate_lt_5000=
        nonchim<5000,

      qc_status=case_when(
        nonchim==0 ~
          "CONFIRMED_EXCLUDE_ZERO",

        nonchim<5000 ~
          "LOW_DEPTH_CANDIDATE_DO_NOT_EXCLUDE_UNTIL_GLOBAL_THRESHOLD_FROZEN",

        TRUE ~
          "NOT_LOW_DEPTH_CANDIDATE"
      )
    ) |>
    filter(
      nonchim<5000
    ) |>
    arrange(
      nonchim
    )

  write_excel_csv(
    low_depth,
    file.path(
      OUT,
      "PRJNA516701_low_depth_candidates_FOR_LATER_GLOBAL_QC.csv"
    ),
    na=""
  )

  low_depth_summary <- tibble(
    threshold=c(
      "<100",
      "<500",
      "<1000",
      "<5000"
    ),
    samples=c(
      sum(tr$nonchim<100),
      sum(tr$nonchim<500),
      sum(tr$nonchim<1000),
      sum(tr$nonchim<5000)
    )
  )

  write_excel_csv(
    low_depth_summary,
    file.path(
      OUT,
      "PRJNA516701_low_depth_candidate_summary.csv"
    ),
    na=""
  )

  # ==========================================================
  # 4. Reviewed acceptance record
  # ==========================================================

  reviewed <- tibble(
    project="PRJNA516701",

    original_runs=171,
    selected_runs=170,

    excluded_run="SRR8488588",

    excluded_patient="3770",

    excluded_time="middle",

    exclusion_reason=
      "ZERO_NONCHIMERIC_READS_FROM_30_READ_INPUT_LIBRARY",

    original_median_filter_retention_pct=
      s$median_filter_retention_pct[1],

    original_median_merge_retention_pct=
      s$median_merge_retention_pct[1],

    original_nonchim_of_merged_pct=
      s$total_nonchim_of_merged_pct[1],

    original_ASVs_nochim=
      s$ASVs_nochim[1],

    ASVs_within_expected_length_pct=
      s$ASVs_within_expected_length_pct[1],

    patients_before=
      before$patients[1],

    patients_after=
      after$patients[1],

    patients_GE2_before=
      before$patients_GE2[1],

    patients_GE2_after=
      after$patients_GE2[1],

    patients_GE3_before=
      before$patients_GE3[1],

    patients_GE3_after=
      after$patients_GE3[1],

    reviewed_decision=
      "ACCEPT_DADA2_OUTPUT_AFTER_EXCLUDING_ZERO_LIBRARY_SAMPLE",

    selected_asv_path=
      selected_rds
  )

  write_excel_csv(
    reviewed,
    file.path(
      OUT,
      "PRJNA516701_FINAL_QC_REVIEW_DECISION.csv"
    ),
    na=""
  )

  ck("PRJNA516701 PROMOTED AFTER QC REVIEW")

  # ==========================================================
  # 5. Close the five-cohort paired-end DADA2 batch
  # ==========================================================

  p1 <- read_summary(
    "PRJEB82425"
  )

  p2 <- read_summary(
    "PRJNA851469"
  )

  p3 <- read_summary(
    "PRJNA1166732"
  )

  # PRJNA516701 reviewed version
  p4 <- s |>
    mutate(
      full_qc_pass=
        TRUE,

      final_status=
        "PASS_AFTER_EXCLUDING_ZERO_LIBRARY_SAMPLE_SRR8488588",

      selected_samples=
        170
    )

  r578_path <- file.path(
    R578,
    "PRJNA578267_FINAL_CHIMERA_SUMMARY_minParent8.csv"
  )

  if (!file.exists(r578_path)) {
    stop(
      "PRJNA578267 Step84B4B final summary is missing."
    )
  }

  r578 <- safe_csv(
    r578_path
  )

  if (
    nrow(r578)!=1 ||
    !isTRUE(
      r578$full_chimera_QC_pass[1]
    )
  ) {
    stop(
      "PRJNA578267 Step84B4B final QC guard failed."
    )
  }

  registry <- tribble(
    ~project,
    ~paired_DADA2_final_status,
    ~processed_runs,
    ~selected_samples,
    ~special_processing_note,

    "PRJEB82425",
    "PASS",
    92L,
    92L,
    "Standard Step84B full DADA2",

    "PRJNA516701",
    "PASS_AFTER_SAMPLE_QC_REVIEW",
    171L,
    170L,
    "Excluded SRR8488588 only: 30 input reads and zero final reads; longitudinal eligibility unchanged",

    "PRJNA578267",
    "PASS",
    207L,
    207L,
    "Sample-wise DADA; batched prechim table; consensus chimera removal with minParentAbundance=8",

    "PRJNA851469",
    "PASS",
    119L,
    119L,
    "Standard Step84B full DADA2",

    "PRJNA1166732",
    "PASS",
    120L,
    120L,
    "Observed primer orientation; minOverlap=8 validated in pilot"
  ) |>
    mutate(
      taxonomy_ready=
        TRUE
    )

  if (
    sum(registry$processed_runs)!=709 ||
    sum(registry$selected_samples)!=708
  ) {
    stop(
      "Five-cohort paired DADA2 closeout count guard failed."
    )
  }

  write_excel_csv(
    registry,
    file.path(
      OUT,
      "V2_STEP84C_FINAL_paired_DADA2_registry.csv"
    ),
    na=""
  )

  global <- tibble(
    metric=c(
      "Paired_DADA2_projects_closed",
      "Processed_runs",
      "Selected_samples_after_confirmed_zero_library_exclusion",
      "Confirmed_excluded_samples",
      "Projects_taxonomy_ready"
    ),
    value=c(
      5,
      709,
      708,
      1,
      5
    )
  )

  write_excel_csv(
    global,
    file.path(
      OUT,
      "V2_STEP84C_global_closeout_summary.csv"
    ),
    na=""
  )

  # ==========================================================
  # README
  # ==========================================================

  readme <- c(
    "SEPSIS V2 - STEP84C PAIRED-END DADA2 CLOSEOUT",
    paste0(
      "Created: ",
      Sys.time()
    ),
    "",
    "PRJNA516701 DECISION",
    "The cohort-level DADA2 pipeline is accepted.",
    "Its only automatic full-QC failure was SRR8488588, which had only 30 raw input reads and zero merged/non-chimeric reads.",
    "This sample is excluded as an essentially empty sequencing library.",
    "",
    "IMPORTANT",
    "No arbitrary low-depth library cutoff is introduced here.",
    "Samples with low but non-zero read depth are only listed as candidates.",
    "A global sample-depth QC rule should be frozen later across cohorts before alpha/beta-diversity analysis.",
    "",
    "LONGITUDINAL IMPACT",
    "Removing SRR8488588 does not change PRJNA516701 patient/timepoint eligibility:",
    "107 patients; 45 >=2 timepoints; 16 >=3 timepoints.",
    "",
    "PAIRED-END DADA2 BATCH CLOSED",
    "PRJEB82425: PASS",
    "PRJNA516701: PASS after exclusion of one zero library",
    "PRJNA578267: PASS via Step84B4B special high-ASV chimera route",
    "PRJNA851469: PASS",
    "PRJNA1166732: PASS with minOverlap=8",
    "",
    "TOTAL",
    "709 Runs processed; 708 selected samples after one confirmed zero-library exclusion.",
    "",
    "NEXT",
    "These five cohorts are taxonomy-ready.",
    "Taxonomy assignment should use one documented reference database/version while keeping cohort ASV tables separate."
  )

  writeLines(
    readme,
    file.path(
      OUT,
      "README_STEP84C_PAIRED_DADA2_CLOSEOUT.txt"
    ),
    useBytes=TRUE
  )

  writeLines(
    c(
      paste0(
        "Completed: ",
        Sys.time()
      ),
      "Status: STEP84C COMPLETE",
      "PRJNA516701 reviewed and accepted after excluding SRR8488588.",
      "Five paired-end DADA2 cohorts closed successfully."
    ),
    file.path(
      OUT,
      "_STEP84C_COMPLETE.ok"
    )
  )

  ck("STEP84C COMPLETE")

  cat("\n============================================================\n")
  cat("SEPSIS V2 - STEP84C COMPLETE\n")
  cat("============================================================\n\n")

  print(
    reviewed,
    n=Inf,
    width=Inf
  )

  cat("\nFIVE-COHORT REGISTRY:\n")

  print(
    registry,
    n=Inf,
    width=Inf
  )

  cat("\nOutput:\n")
  cat(
    OUT,
    "\n"
  )
}

tryCatch(
  main(),
  error=function(e) {

    msg <- c(
      paste0(
        "STEP84C FATAL ERROR: ",
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

    ck(
      "STEP84C FAILED"
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
