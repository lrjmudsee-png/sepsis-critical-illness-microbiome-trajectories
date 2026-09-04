# ============================================================
# Sepsis V2 - Step 81C
# CORRECTED FINAL cohort freeze
#
# Fixes a semantic bug in Step81B:
# - Legacy Step73 16S cohorts use data_modality = "16"
# - Rescued cohorts use data_modality = "16S"
# - PRJNA884103 had missing modality
#
# Step81B therefore excluded legacy 16S cohorts and accidentally
# included PRJNA884103 in the primary 16S layer.
#
# Step81C:
# 1. Loads the successfully generated Step81B frozen master
# 2. Normalizes modality based on verified project identity
# 3. Recomputes primary 16S inclusion
# 4. Recomputes sequence longitudinal eligibility
# 5. Regenerates all summaries/manifests
# 6. Creates a separate shotgun external manifest
#
# R 4.4.0 / Windows
# ============================================================

options(stringsAsFactors = FALSE)

pkgs <- c("readr","dplyr","tidyr","stringr","purrr","tibble")
missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly=TRUE)]
if (length(missing)) install.packages(missing, repos="https://cloud.r-project.org")

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(purrr)
  library(tibble)
})

ROOT <- "E:/sepsis_project"

IN <- file.path(
  ROOT,
  "results",
  "V2_21B_FINAL_COHORT_FREEZE"
)

OUT <- file.path(
  ROOT,
  "results",
  "V2_21C_FINAL_COHORT_FREEZE_CORRECTED"
)

dir.create(
  OUT,
  recursive=TRUE,
  showWarnings=FALSE
)

LOG <- file.path(
  OUT,
  "_STEP81C_runtime_checkpoints.txt"
)

ERR <- file.path(
  OUT,
  "_STEP81C_FATAL_ERROR.txt"
)

if (file.exists(ERR)) unlink(ERR)

ck <- function(x) {
  cat(
    paste0(x, ": ", Sys.time(), "\n"),
    file=LOG,
    append=TRUE
  )
}

clean_chr <- function(x) {
  x <- trimws(as.character(x))
  x[
    is.na(x) |
    x=="" |
    tolower(x) %in% c(
      "na","nan","n/a","null","none"
    )
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

find_latest <- function(folder, pattern) {
  x <- list.files(
    folder,
    pattern=pattern,
    full.names=TRUE
  )

  if (!length(x)) return(NA_character_)

  x[
    which.max(
      file.info(x)$mtime
    )
  ]
}

count_sequence_files <- function(project) {

  roots <- c(
    file.path(ROOT,"data",project)
  )

  if (project=="PRJNA533528") {
    roots <- c(
      roots,
      "F:/sepsis/PRJNA533528"
    )
  }

  files <- character()

  for (r in roots) {
    if (!dir.exists(r)) next

    files <- c(
      files,
      list.files(
        r,
        recursive=TRUE,
        full.names=TRUE
      )
    )
  }

  files <- unique(files)

  sum(
    str_detect(
      basename(files),
      regex(
        "\\.(fastq|fq|fa|fasta)\\.gz$",
        ignore_case=TRUE
      )
    )
  )
}

main <- function() {

  ck("STEP81C STARTED")

  master_path <- find_latest(
    IN,
    "^V2_FROZEN_master_metadata_[0-9]{8}_[0-9]{6}\\.csv$"
  )

  provenance_path <- file.path(
    IN,
    "V2_STEP81B_input_provenance.csv"
  )

  if (
    is.na(master_path) ||
    !file.exists(master_path)
  ) {
    stop(
      "Step81B frozen master not found."
    )
  }

  full <- safe_csv(master_path)

  ck("STEP81B MASTER READ")

  # ==========================================================
  # Hard guards
  # ==========================================================

  if (nrow(full) != 2693) {
    stop(
      paste0(
        "Unexpected Step81B row count: ",
        nrow(full),
        "; expected 2693."
      )
    )
  }

  if (
    n_distinct(full$project) != 18
  ) {
    stop(
      paste0(
        "Unexpected project count: ",
        n_distinct(full$project),
        "; expected 18."
      )
    )
  }

  if (
    n_distinct(
      full$run_uid[
        !is.na(full$run_uid)
      ]
    ) !=
    sum(!is.na(full$run_uid))
  ) {
    stop(
      "Duplicate run_uid detected before Step81C correction."
    )
  }

  ck("INPUT GUARDS PASSED")

  # ==========================================================
  # Verified modality normalization
  #
  # Current 18-project frozen registry:
  # - PRJNA884103 is shotgun
  # - All other projects are 16S amplicon cohorts
  #
  # PRJNA1125274 is 16S modality but sequence-deferred.
  # ==========================================================

  full <- full |>
    mutate(
      data_modality_original=
        clean_chr(data_modality),

      data_modality=case_when(
        project=="PRJNA884103" ~
          "SHOTGUN",

        TRUE ~
          "16S"
      )
    )

  modality_audit <- full |>
    distinct(
      project,
      data_modality_original,
      data_modality
    ) |>
    arrange(project)

  write_excel_csv(
    modality_audit,
    file.path(
      OUT,
      "V2_STEP81C_modality_correction_audit.csv"
    ),
    na=""
  )

  # ==========================================================
  # Recompute primary sequence inclusion
  # ==========================================================

  # Explicit FALSE from rescued cohorts must be respected.
  seq_flag_false <- rep(FALSE,nrow(full))

  if (
    "sequence_analysis_include" %in%
    names(full)
  ) {

    x <- tolower(
      clean_chr(
        full$sequence_analysis_include
      )
    )

    seq_flag_false <-
      x %in% c(
        "false","f","0","no","n"
      )
  }

  full <- full |>
    mutate(
      primary_16s_sequence_include=
        data_modality=="16S" &
        project!="PRJNA1125274" &
        !is.na(run_id) &
        !seq_flag_false,

      shotgun_external_include=
        data_modality=="SHOTGUN",

      freeze_status=case_when(
        project=="PRJNA1125274" ~
          "FROZEN_METADATA_ONLY_PROVISIONAL_SEQUENCE_DEFERRED",

        data_modality=="SHOTGUN" ~
          "FROZEN_SHOTGUN_EXTERNAL_SEPARATE_LAYER",

        data_modality=="16S" &
        primary_16s_sequence_include ~
          "FROZEN_16S_SEQUENCE_LAYER",

        data_modality=="16S" &
        !primary_16s_sequence_include ~
          "FROZEN_METADATA_SEQUENCE_UNAVAILABLE_OR_UNMAPPED",

        TRUE ~
          "FROZEN_METADATA"
      )
    )

  # ==========================================================
  # Recompute metadata time key
  # ==========================================================

  full <- full |>
    mutate(
      time_key=coalesce(
        clean_chr(time_class),
        clean_chr(time_raw)
      )
    )

  # Metadata eligibility
  meta_pt <- full |>
    filter(
      !is.na(patient_uid),
      !is.na(time_key)
    ) |>
    distinct(
      patient_uid,
      time_key
    ) |>
    count(
      patient_uid,
      name="n_distinct_timepoints_new"
    )

  full <- full |>
    select(
      -any_of(c(
        "n_distinct_timepoints",
        "longitudinal_ge2",
        "longitudinal_ge3"
      ))
    ) |>
    left_join(
      meta_pt,
      by="patient_uid"
    ) |>
    rename(
      n_distinct_timepoints=
        n_distinct_timepoints_new
    ) |>
    mutate(
      n_distinct_timepoints=coalesce(
        n_distinct_timepoints,
        0L
      ),

      longitudinal_ge2=
        n_distinct_timepoints>=2,

      longitudinal_ge3=
        n_distinct_timepoints>=3
    )

  # Sequence-available eligibility
  seq_pt <- full |>
    filter(
      primary_16s_sequence_include,
      !is.na(patient_uid),
      !is.na(time_key)
    ) |>
    distinct(
      patient_uid,
      time_key
    ) |>
    count(
      patient_uid,
      name="sequence_n_distinct_timepoints_new"
    )

  full <- full |>
    select(
      -any_of(c(
        "sequence_n_distinct_timepoints",
        "sequence_longitudinal_ge2",
        "sequence_longitudinal_ge3"
      ))
    ) |>
    left_join(
      seq_pt,
      by="patient_uid"
    ) |>
    rename(
      sequence_n_distinct_timepoints=
        sequence_n_distinct_timepoints_new
    ) |>
    mutate(
      sequence_n_distinct_timepoints=coalesce(
        sequence_n_distinct_timepoints,
        0L
      ),

      sequence_longitudinal_ge2=
        sequence_n_distinct_timepoints>=2,

      sequence_longitudinal_ge3=
        sequence_n_distinct_timepoints>=3
    )

  ck("MODALITY AND SEQUENCE ELIGIBILITY CORRECTED")

  # ==========================================================
  # QC: run duplicates/conflicts remain zero
  # ==========================================================

  run_conflicts <- full |>
    filter(
      !is.na(run_uid)
    ) |>
    group_by(run_uid) |>
    summarise(
      n_rows=n(),

      patients=n_distinct(
        patient_id,
        na.rm=TRUE
      ),

      samples=n_distinct(
        sample_id,
        na.rm=TRUE
      ),

      timepoints=n_distinct(
        time_key,
        na.rm=TRUE
      ),

      .groups="drop"
    ) |>
    filter(
      patients>1 |
      samples>1 |
      timepoints>1
    )

  dup_run <- full |>
    filter(
      !is.na(run_uid)
    ) |>
    count(
      run_uid,
      name="n"
    ) |>
    filter(n>1)

  write_excel_csv(
    run_conflicts,
    file.path(
      OUT,
      "V2_FROZEN_run_conflict_QC.csv"
    ),
    na=""
  )

  write_excel_csv(
    dup_run,
    file.path(
      OUT,
      "V2_FROZEN_duplicate_run_QC.csv"
    ),
    na=""
  )

  if (
    nrow(run_conflicts)>0 ||
    nrow(dup_run)>0
  ) {
    stop(
      "Run conflict/duplicate appeared after correction."
    )
  }

  # ==========================================================
  # Patient summary
  # ==========================================================

  patient_summary <- full |>
    filter(
      !is.na(patient_uid)
    ) |>
    group_by(
      project,
      patient_uid,
      patient_id,
      analysis_module
    ) |>
    summarise(
      n_rows=n(),

      n_samples=n_distinct(
        sample_uid,
        na.rm=TRUE
      ),

      n_runs=n_distinct(
        run_uid,
        na.rm=TRUE
      ),

      n_timepoints=n_distinct(
        time_key[
          !is.na(time_key)
        ]
      ),

      n_sequence_timepoints=n_distinct(
        time_key[
          primary_16s_sequence_include &
          !is.na(time_key)
        ]
      ),

      metadata_longitudinal_ge2=
        n_timepoints>=2,

      metadata_longitudinal_ge3=
        n_timepoints>=3,

      sequence_longitudinal_ge2=
        n_sequence_timepoints>=2,

      sequence_longitudinal_ge3=
        n_sequence_timepoints>=3,

      .groups="drop"
    )

  # ==========================================================
  # Cohort summary
  # ==========================================================

  cohort_summary <- full |>
    group_by(
      project,
      analysis_module,
      data_modality
    ) |>
    summarise(
      rows=n(),

      patients=n_distinct(
        patient_uid,
        na.rm=TRUE
      ),

      samples=n_distinct(
        sample_uid,
        na.rm=TRUE
      ),

      mapped_runs=n_distinct(
        run_uid,
        na.rm=TRUE
      ),

      missing_run_rows=sum(
        is.na(run_id)
      ),

      metadata_patients_ge2=n_distinct(
        patient_uid[
          longitudinal_ge2
        ],
        na.rm=TRUE
      ),

      metadata_patients_ge3=n_distinct(
        patient_uid[
          longitudinal_ge3
        ],
        na.rm=TRUE
      ),

      sequence_patients_ge2=n_distinct(
        patient_uid[
          sequence_longitudinal_ge2
        ],
        na.rm=TRUE
      ),

      sequence_patients_ge3=n_distinct(
        patient_uid[
          sequence_longitudinal_ge3
        ],
        na.rm=TRUE
      ),

      primary_16s_rows=sum(
        primary_16s_sequence_include
      ),

      shotgun_rows=sum(
        shotgun_external_include
      ),

      outcome_nonmissing=sum(
        !is.na(outcome)
      ),

      infection_group_nonmissing=sum(
        !is.na(infection_group)
      ),

      .groups="drop"
    )

  # ==========================================================
  # Module summary
  # ==========================================================

  module_summary <- full |>
    group_by(
      analysis_module
    ) |>
    summarise(
      projects=n_distinct(project),

      patients=n_distinct(
        patient_uid,
        na.rm=TRUE
      ),

      samples=n_distinct(
        sample_uid,
        na.rm=TRUE
      ),

      mapped_runs=n_distinct(
        run_uid,
        na.rm=TRUE
      ),

      metadata_patients_ge2=n_distinct(
        patient_uid[
          longitudinal_ge2
        ],
        na.rm=TRUE
      ),

      sequence_patients_ge2=n_distinct(
        patient_uid[
          sequence_longitudinal_ge2
        ],
        na.rm=TRUE
      ),

      primary_16s_rows=sum(
        primary_16s_sequence_include
      ),

      shotgun_rows=sum(
        shotgun_external_include
      ),

      .groups="drop"
    )

  # ==========================================================
  # Registry
  # ==========================================================

  registry <- cohort_summary |>
    mutate(
      freeze_decision=case_when(
        project=="PRJNA1125274" ~
          "METADATA_ONLY_DEFERRED",

        data_modality=="SHOTGUN" ~
          "SHOTGUN_EXTERNAL_SEPARATE_LAYER",

        analysis_module=="STATIC_SUPPORT" ~
          "STATIC_SUPPORT",

        TRUE ~
          "INCLUDE"
      ),

      next_sequence_action=case_when(
        project=="PRJNA1125274" ~
          "Finish correct 289-Run raw download/mapping before DADA2",

        project=="PRJNA884103" ~
          "Shotgun external layer; do not enter DADA2/16S workflow",

        project=="CRA002354" ~
          "Process cohort-specific 16S reads / ASV table",

        project=="PRJNA851469" ~
          "Use only 119 strict exact-mapped samples in 16S sequence layer",

        TRUE ~
          "Audit existing ASV and raw FASTQ; process cohort separately"
      )
    )

  # ==========================================================
  # Primary 16S manifests
  # ==========================================================

  manifest16s <- full |>
    filter(
      primary_16s_sequence_include
    ) |>
    select(
      project,
      patient_id,
      patient_uid,
      sample_id,
      sample_uid,
      biosample,
      experiment_id,
      run_id,
      run_uid,
      time_raw,
      time_class,
      time_day,
      time_order,
      analysis_module,
      cohort_role,
      sequence_status,
      freeze_status
    ) |>
    arrange(
      project,
      patient_id,
      time_order,
      time_day,
      sample_id
    )

  manifest16s_longitudinal <- full |>
    filter(
      primary_16s_sequence_include,
      sequence_longitudinal_ge2
    ) |>
    select(
      project,
      patient_id,
      patient_uid,
      sample_id,
      sample_uid,
      run_id,
      run_uid,
      time_raw,
      time_class,
      time_day,
      time_order,
      analysis_module,
      cohort_role,
      sequence_n_distinct_timepoints
    ) |>
    arrange(
      project,
      patient_id,
      time_order,
      time_day,
      sample_id
    )

  # ==========================================================
  # Shotgun manifest (separate layer)
  # ==========================================================

  shotgun_manifest <- full |>
    filter(
      shotgun_external_include
    ) |>
    select(
      project,
      patient_id,
      patient_uid,
      sample_id,
      sample_uid,
      biosample,
      run_id,
      run_uid,
      time_raw,
      time_class,
      time_day,
      time_order,
      analysis_module,
      cohort_role,
      sequence_status,
      freeze_status
    ) |>
    arrange(
      project,
      patient_id,
      time_order,
      time_day,
      sample_id
    )

  # ==========================================================
  # Deferred registry
  # ==========================================================

  deferred <- full |>
    filter(
      project=="PRJNA1125274" |
      (
        data_modality=="16S" &
        !primary_16s_sequence_include
      ) |
      data_modality=="SHOTGUN"
    ) |>
    select(
      project,
      patient_id,
      sample_id,
      run_id,
      time_raw,
      analysis_module,
      data_modality,
      sequence_status,
      freeze_status,
      notes
    )

  # ==========================================================
  # Local sequence readiness
  # ==========================================================

  projects <- sort(
    unique(
      full$project[
        !is.na(full$project)
      ]
    )
  )

  seq_readiness <- tibble(
    project=projects,

    local_sequence_files=map_int(
      projects,
      count_sequence_files
    )
  ) |>
    left_join(
      cohort_summary |>
        select(
          project,
          data_modality,
          mapped_runs,
          missing_run_rows,
          primary_16s_rows,
          shotgun_rows
        ),
      by="project"
    ) |>
    mutate(
      local_sequence_file_status=case_when(
        project=="PRJNA1125274" ~
          "LOCAL_FILES_EXIST_BUT_PROJECT_IDENTITY/289_RUN_REPAIR_NOT_FINAL",

        data_modality=="SHOTGUN" ~
          "SHOTGUN_SEPARATE_LAYER_RAW_DEFERRED",

        local_sequence_files>0 ~
          "LOCAL_SEQUENCE_FILES_FOUND",

        primary_16s_rows>0 ~
          "NO_RAW_FILES_DETECTED_CHECK_EXISTING_ASV_OR_DOWNLOAD",

        TRUE ~
          "NO_SEQUENCE_ACTION_IN_PRIMARY_LAYER"
      )
    )

  # ==========================================================
  # Global summary
  # ==========================================================

  global_summary <- tibble(
    metric=c(
      "Rows",
      "Projects",
      "Patients",
      "Samples",
      "Mapped_Runs",
      "Metadata_Patients_GE2",
      "Metadata_Patients_GE3",
      "Primary_16S_Sequence_Rows",
      "Primary_16S_Patients_GE2",
      "Primary_16S_Patients_GE3",
      "Shotgun_Metadata_Rows",
      "Deferred_PRJNA1125274_Rows",
      "Strict_Unmapped_16S_Rows",
      "Missing_Run_Rows",
      "Duplicate_Run_UID",
      "Run_Conflict_Rows"
    ),

    value=c(
      nrow(full),

      n_distinct(
        full$project
      ),

      n_distinct(
        full$patient_uid,
        na.rm=TRUE
      ),

      n_distinct(
        full$sample_uid,
        na.rm=TRUE
      ),

      n_distinct(
        full$run_uid,
        na.rm=TRUE
      ),

      n_distinct(
        full$patient_uid[
          full$longitudinal_ge2
        ],
        na.rm=TRUE
      ),

      n_distinct(
        full$patient_uid[
          full$longitudinal_ge3
        ],
        na.rm=TRUE
      ),

      sum(
        full$primary_16s_sequence_include
      ),

      n_distinct(
        full$patient_uid[
          full$sequence_longitudinal_ge2
        ],
        na.rm=TRUE
      ),

      n_distinct(
        full$patient_uid[
          full$sequence_longitudinal_ge3
        ],
        na.rm=TRUE
      ),

      sum(
        full$shotgun_external_include
      ),

      sum(
        full$project=="PRJNA1125274"
      ),

      sum(
        full$data_modality=="16S" &
        full$project!="PRJNA1125274" &
        !full$primary_16s_sequence_include
      ),

      sum(
        is.na(full$run_id)
      ),

      nrow(dup_run),

      nrow(run_conflicts)
    )
  )

  # Hard expected corrected values
  expected <- c(
    Rows=2693,
    Projects=18,
    Primary_16S_Sequence_Rows=1941,
    Primary_16S_Patients_GE2=548,
    Primary_16S_Patients_GE3=210,
    Shotgun_Metadata_Rows=462,
    Deferred_PRJNA1125274_Rows=286,
    Strict_Unmapped_16S_Rows=4,
    Duplicate_Run_UID=0,
    Run_Conflict_Rows=0
  )

  for (nm in names(expected)) {

    observed <- global_summary$value[
      global_summary$metric==nm
    ]

    if (
      length(observed)!=1 ||
      observed != expected[[nm]]
    ) {
      stop(
        paste0(
          "Corrected freeze guard failed for ",
          nm,
          ": observed=",
          paste(observed,collapse=","),
          " expected=",
          expected[[nm]]
        )
      )
    }
  }

  ck("CORRECTED GLOBAL GUARDS PASSED")

  # ==========================================================
  # Freeze version + write
  # ==========================================================

  stamp <- format(
    Sys.time(),
    "%Y%m%d_%H%M%S"
  )

  freeze_version <- paste0(
    "STEP81C_CORRECTED_",
    stamp
  )

  full$freeze_version <- freeze_version

  full_path <- file.path(
    OUT,
    paste0(
      "V2_FINAL_FROZEN_master_metadata_",
      stamp,
      ".csv"
    )
  )

  write_excel_csv(
    full,
    full_path,
    na=""
  )

  write_excel_csv(
    patient_summary,
    file.path(
      OUT,
      paste0(
        "V2_FINAL_FROZEN_patient_summary_",
        stamp,
        ".csv"
      )
    ),
    na=""
  )

  write_excel_csv(
    cohort_summary,
    file.path(
      OUT,
      paste0(
        "V2_FINAL_FROZEN_cohort_summary_",
        stamp,
        ".csv"
      )
    ),
    na=""
  )

  write_excel_csv(
    global_summary,
    file.path(
      OUT,
      paste0(
        "V2_FINAL_FROZEN_global_summary_",
        stamp,
        ".csv"
      )
    ),
    na=""
  )

  write_excel_csv(
    module_summary,
    file.path(
      OUT,
      "V2_FINAL_FROZEN_analysis_module_summary.csv"
    ),
    na=""
  )

  write_excel_csv(
    registry,
    file.path(
      OUT,
      "V2_FINAL_FROZEN_project_registry.csv"
    ),
    na=""
  )

  write_excel_csv(
    manifest16s,
    file.path(
      OUT,
      "V2_FINAL_FROZEN_16S_sequence_manifest.csv"
    ),
    na=""
  )

  write_excel_csv(
    manifest16s_longitudinal,
    file.path(
      OUT,
      "V2_FINAL_FROZEN_16S_longitudinal_manifest.csv"
    ),
    na=""
  )

  write_excel_csv(
    shotgun_manifest,
    file.path(
      OUT,
      "V2_FINAL_FROZEN_SHOTGUN_external_manifest.csv"
    ),
    na=""
  )

  write_excel_csv(
    deferred,
    file.path(
      OUT,
      "V2_FINAL_FROZEN_deferred_registry.csv"
    ),
    na=""
  )

  write_excel_csv(
    seq_readiness,
    file.path(
      OUT,
      "V2_FINAL_FROZEN_sequence_readiness_by_project.csv"
    ),
    na=""
  )

  # Carry provenance forward
  if (file.exists(provenance_path)) {
    file.copy(
      provenance_path,
      file.path(
        OUT,
        "V2_FINAL_FROZEN_input_provenance_from_STEP81B.csv"
      ),
      overwrite=TRUE
    )
  }

  # ==========================================================
  # README
  # ==========================================================

  readme <- c(
    "SEPSIS V2 - STEP81C CORRECTED FINAL COHORT FREEZE",
    paste0("Freeze version: ",freeze_version),
    paste0("Created: ",Sys.time()),
    "",
    "IMPORTANT",
    "Step81B executed successfully but contained a sequence-layer modality normalization bug.",
    "Legacy Step73 16S cohorts were labeled '16', rescued cohorts were labeled '16S', and PRJNA884103 had missing modality.",
    "Step81C corrected this before final freeze.",
    "",
    "FINAL VERIFIED MODALITY RULE",
    "- PRJNA884103 = SHOTGUN external layer.",
    "- All other currently frozen projects = 16S amplicon.",
    "- PRJNA1125274 is a 16S project but remains sequence-deferred/provisional.",
    "",
    "FINAL 16S SEQUENCE RULE",
    "- Include 16S rows with a mapped Run.",
    "- Exclude PRJNA1125274 until correct 289-Run raw set is finalized.",
    "- PRJNA851469: include only 119 strict exact-mapped samples.",
    "- Shotgun data never enter the DADA2/16S manifest.",
    "",
    "EXPECTED VERIFIED COUNTS",
    "- Rows = 2693",
    "- Projects = 18",
    "- Primary 16S sequence rows = 1941",
    "- Primary 16S patients >=2 timepoints = 548",
    "- Primary 16S patients >=3 timepoints = 210",
    "- Shotgun metadata rows = 462",
    "- PRJNA1125274 deferred rows = 286",
    "- Strict unmapped eligible 16S rows = 4",
    "- Duplicate Run UID = 0",
    "- Run conflict rows = 0",
    "",
    "UNIQUE DOWNSTREAM METADATA ENTRY POINT",
    full_path,
    "",
    "NEXT STEP",
    "Step82 must use V2_FINAL_FROZEN_16S_sequence_manifest.csv, not the Step81B manifest.",
    "Step82 will determine for each 16S cohort whether an existing ASV table is reusable or DADA2/raw-read processing is required."
  )

  writeLines(
    readme,
    file.path(
      OUT,
      "README_STEP81C_CORRECTED_FINAL_FREEZE.txt"
    ),
    useBytes=TRUE
  )

  ck("STEP81C COMPLETE")

  cat("\n============================================================\n")
  cat("SEPSIS V2 - STEP81C CORRECTED FINAL FREEZE COMPLETE\n")
  cat("============================================================\n\n")

  print(
    global_summary,
    n=Inf,
    width=Inf
  )

  cat("\nOutput:\n")
  cat(OUT,"\n")
  cat("============================================================\n")
}

tryCatch(
  main(),
  error=function(e) {

    msg <- c(
      paste0(
        "STEP81C FATAL ERROR: ",
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

    ck("STEP81C FAILED")

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
