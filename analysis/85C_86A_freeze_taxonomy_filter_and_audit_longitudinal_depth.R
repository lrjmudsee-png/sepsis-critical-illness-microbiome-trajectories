# ============================================================
# Sepsis V2 - Step 85C / 86A
# Freeze taxonomy feature rule + audit longitudinal impact
# of candidate sample-depth thresholds
#
# IMPORTANT
#   - NO broad filesystem scan.
#   - Reads only:
#       E:/sepsis_project/data/_V2_ANALYSIS_READY/00_FREEZE
#       E:/sepsis_project/results/V2_25B2_TAXONOMY_QC_AND_DEPTH_AUDIT_FIXED
#       per-project known canonical folders
#   - Does NOT yet delete low-depth samples.
#
# Taxonomy rule frozen here:
#   KEEP:
#     SILVA Kingdom == Bacteria or Archaea
#   REMOVE:
#     chloroplast / mitochondria annotations
#     other kingdoms
#     Kingdom-unassigned features
#
# Why this is acceptable from Step85B2:
#   36,681 ASVs -> 36,004 retained
#   66,428,902 reads -> 66,413,963 retained
#   ~99.98% reads retained globally
#
# Step86A then audits candidate minimum library sizes:
#   100 / 500 / 1000 / 2000 / 5000 / 10000 reads
# and quantifies impact on:
#   - samples
#   - patients
#   - patients with >=2 timepoints
#   - patients with >=3 timepoints
#
# R 4.4.0 / Windows
# ============================================================

options(stringsAsFactors = FALSE)

pkgs <- c("readr","dplyr","stringr","tibble","purrr")
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
  library(stringr)
  library(tibble)
  library(purrr)
})

ROOT <- "E:/sepsis_project"

CANON <- file.path(
  ROOT,
  "data",
  "_V2_ANALYSIS_READY"
)

FREEZE_DIR <- file.path(
  CANON,
  "00_FREEZE"
)

PROJECT_ROOT <- file.path(
  CANON,
  "01_PROJECTS"
)

STEP85B2 <- file.path(
  ROOT,
  "results",
  "V2_25B2_TAXONOMY_QC_AND_DEPTH_AUDIT_FIXED"
)

OUT <- file.path(
  ROOT,
  "results",
  "V2_25C_86A_TAXONOMY_FILTER_AND_LONGITUDINAL_DEPTH_AUDIT"
)

dir.create(
  OUT,
  recursive=TRUE,
  showWarnings=FALSE
)

LOG <- file.path(
  OUT,
  "_STEP85C86A_runtime_checkpoints.txt"
)

ERR <- file.path(
  OUT,
  "_STEP85C86A_FATAL_ERROR.txt"
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
  saveRDS(obj,tmp)
  if (file.exists(p)) unlink(p)
  if (!file.rename(tmp,p)) {
    stop("Could not atomically save: ",p)
  }
}

pick_col <- function(df, candidates, required=TRUE) {

  exact <- candidates[
    candidates %in% names(df)
  ]

  if (length(exact)) {
    return(exact[1])
  }

  lower_names <- tolower(names(df))
  lower_candidates <- tolower(candidates)

  idx <- match(
    lower_candidates,
    lower_names
  )

  idx <- idx[
    !is.na(idx)
  ]

  if (length(idx)) {
    return(names(df)[idx[1]])
  }

  if (required) {
    stop(
      paste0(
        "Could not find required column among: ",
        paste(candidates,collapse=", ")
      )
    )
  }

  NA_character_
}

registry <- tribble(
  ~project, ~expected_samples,

  "PRJNA430161", 26L,
  "PRJNA691455", 49L,
  "PRJNA978257", 13L,
  "PRJEB82425", 92L,
  "PRJNA516701", 170L,
  "PRJNA578267", 207L,
  "PRJNA851469", 119L,
  "PRJNA1166732", 120L
)

thresholds <- c(
  100L,
  500L,
  1000L,
  2000L,
  5000L,
  10000L
)

# ------------------------------------------------------------
# Locate frozen master ONLY inside canonical 00_FREEZE
# ------------------------------------------------------------

master_candidates <- list.files(
  FREEZE_DIR,
  pattern="^V2_FINAL_FROZEN_master_metadata_CLINICAL_REPAIRED.*\\.csv$",
  full.names=TRUE,
  recursive=FALSE
)

if (!length(master_candidates)) {
  stop(
    "No Step81D frozen master found in canonical 00_FREEZE."
  )
}

if (length(master_candidates)>1) {

  info <- file.info(
    master_candidates
  )

  master_path <- rownames(
    info
  )[
    which.max(
      info$mtime
    )
  ]

} else {

  master_path <- master_candidates[1]
}

# ------------------------------------------------------------
# Main
# ------------------------------------------------------------

main <- function() {

  ck("STEP85C86A STARTED")

  # ==========================================================
  # 1. Hard-guard Step85B2 global result
  # ==========================================================

  global85 <- safe_csv(
    file.path(
      STEP85B2,
      "V2_STEP85B_global_summary.csv"
    )
  )

  value_of <- function(metric_name) {

    x <- global85 |>
      filter(
        metric==metric_name
      )

    if (nrow(x)!=1) {
      stop(
        "Missing Step85B2 global metric: ",
        metric_name
      )
    }

    as.numeric(
      x$value[1]
    )
  }

  if (
    value_of("projects")!=8 ||
    value_of("samples")!=796 ||
    value_of("ASVs_before")!=36681 ||
    value_of("reads_before")!=66428902 ||
    value_of("preview_keep_ASVs")!=36004 ||
    value_of("preview_keep_reads")!=66413963 ||
    value_of("samples_zero_after_preview")!=0
  ) {
    stop(
      "Step85B2 global hard guard failed."
    )
  }

  ck("STEP85B2 GLOBAL GUARDS PASSED")

  # ==========================================================
  # 2. Freeze taxonomy feature-filtered sequence tables
  # ==========================================================

  depth_rows <- list()
  freeze_rows <- list()

  for (i in seq_len(nrow(registry))) {

    project <- registry$project[i]

    preview_path <- file.path(
      STEP85B2,
      paste0(
        project,
        "_PREVIEW_prokaryote_nonorganelle_seqtab.rds"
      )
    )

    flags_path <- file.path(
      STEP85B2,
      paste0(
        project,
        "_taxonomy_feature_flags.csv"
      )
    )

    depth_path <- file.path(
      STEP85B2,
      paste0(
        project,
        "_sample_depth_before_after_taxonomy_preview.csv"
      )
    )

    if (
      !file.exists(preview_path) ||
      !file.exists(flags_path) ||
      !file.exists(depth_path)
    ) {
      stop(
        project,
        ": Step85B2 inputs missing."
      )
    }

    seqtab <- readRDS(
      preview_path
    )

    flags <- safe_csv(
      flags_path
    )

    depth <- safe_csv(
      depth_path
    )

    if (
      !is.matrix(seqtab) ||
      nrow(seqtab)!=registry$expected_samples[i] ||
      any(rowSums(seqtab)==0)
    ) {
      stop(
        project,
        ": taxonomy-filtered preview guard failed."
      )
    }

    if (
      nrow(depth)!=registry$expected_samples[i]
    ) {
      stop(
        project,
        ": sample-depth row-count guard failed."
      )
    }

    if (
      !all(
        depth$sample_id %in%
        rownames(seqtab)
      )
    ) {
      stop(
        project,
        ": sample IDs do not match preview sequence table."
      )
    }

    # Freeze into a dedicated work folder, not 02_asv_selected.
    freeze_dir <- file.path(
      PROJECT_ROOT,
      project,
      "05_work",
      "step85C_taxonomy_filtered"
    )

    dir.create(
      freeze_dir,
      recursive=TRUE,
      showWarnings=FALSE
    )

    final_seq_path <- file.path(
      freeze_dir,
      paste0(
        project,
        "_seqtab_taxonomy_filtered_SILVA1382.rds"
      )
    )

    atomic_saveRDS(
      seqtab,
      final_seq_path
    )

    kept_flags <- flags |>
      filter(
        preview_keep
      )

    removed_flags <- flags |>
      filter(
        !preview_keep
      )

    write_excel_csv(
      kept_flags,
      file.path(
        freeze_dir,
        paste0(
          project,
          "_taxonomy_kept_features.csv"
        )
      ),
      na=""
    )

    write_excel_csv(
      removed_flags,
      file.path(
        freeze_dir,
        paste0(
          project,
          "_taxonomy_removed_features.csv"
        )
      ),
      na=""
    )

    freeze_rows[[
      length(freeze_rows)+1
    ]] <- tibble(
      project=project,
      samples=nrow(seqtab),
      retained_ASVs=ncol(seqtab),
      retained_reads=sum(seqtab),
      frozen_seqtab_path=final_seq_path
    )

    depth_rows[[
      length(depth_rows)+1
    ]] <- depth |>
      mutate(
        project=project
      )

    ck(
      paste0(
        project,
        " TAXONOMY FILTER FROZEN"
      )
    )
  }

  freeze_registry <- bind_rows(
    freeze_rows
  )

  depth_all <- bind_rows(
    depth_rows
  )

  if (
    sum(freeze_registry$samples)!=796 ||
    sum(freeze_registry$retained_ASVs)!=36004 ||
    sum(freeze_registry$retained_reads)!=66413963
  ) {
    stop(
      "Frozen taxonomy-filtered global count guard failed."
    )
  }

  write_excel_csv(
    freeze_registry,
    file.path(
      OUT,
      "V2_STEP85C_taxonomy_filter_freeze_registry.csv"
    ),
    na=""
  )

  # ==========================================================
  # 3. Explicitly flag taxonomy-loss outliers
  # ==========================================================

  taxonomy_loss_outliers <- depth_all |>
    filter(
      removed_pct>=5
    ) |>
    arrange(
      desc(
        removed_pct
      )
    )

  write_excel_csv(
    taxonomy_loss_outliers,
    file.path(
      OUT,
      "V2_STEP85C_taxonomy_loss_outliers_GE5pct.csv"
    ),
    na=""
  )

  # Known strong outlier from Step85B2 should be present.
  if (
    !any(
      taxonomy_loss_outliers$project=="PRJNA516701" &
      taxonomy_loss_outliers$sample_id=="SRR8488562" &
      taxonomy_loss_outliers$removed_pct>80
    )
  ) {
    stop(
      "Expected PRJNA516701/SRR8488562 taxonomy-loss outlier not found."
    )
  }

  # ==========================================================
  # 4. Load final frozen master and map Run -> patient/time
  # ==========================================================

  master <- safe_csv(
    master_path
  )

  project_col <- pick_col(
    master,
    c(
      "project",
      "Project",
      "study",
      "Study"
    )
  )

  run_col <- pick_col(
    master,
    c(
      "run_id",
      "Run",
      "run",
      "Run_ID",
      "run_accession"
    )
  )

  patient_col <- pick_col(
    master,
    c(
      "patient_id",
      "Patient_ID",
      "patient",
      "subject_id",
      "Subject_ID"
    )
  )

  sample_col <- pick_col(
    master,
    c(
      "sample_id",
      "Sample_ID",
      "sample",
      "biosample"
    ),
    required=FALSE
  )

  time_raw_col <- pick_col(
    master,
    c(
      "time_raw",
      "Time_raw",
      "timepoint",
      "Timepoint",
      "time_point"
    ),
    required=FALSE
  )

  time_day_col <- pick_col(
    master,
    c(
      "time_day",
      "Time_day",
      "day",
      "Day"
    ),
    required=FALSE
  )

  master_map <- master |>
    transmute(
      project=.data[[project_col]],
      run_id=.data[[run_col]],
      patient_id=.data[[patient_col]],

      source_sample_id=
        if (!is.na(sample_col)) {
          as.character(
            .data[[sample_col]]
          )
        } else {
          NA_character_
        },

      time_raw=
        if (!is.na(time_raw_col)) {
          as.character(
            .data[[time_raw_col]]
          )
        } else {
          NA_character_
        },

      time_day=
        if (!is.na(time_day_col)) {
          suppressWarnings(
            as.numeric(
              .data[[time_day_col]]
            )
          )
        } else {
          NA_real_
        }
    ) |>
    filter(
      project %in%
      registry$project
    ) |>
    distinct(
      project,
      run_id,
      .keep_all=TRUE
    )

  depth_meta <- depth_all |>
    rename(
      run_id=sample_id
    ) |>
    left_join(
      master_map,
      by=c(
        "project",
        "run_id"
      )
    )

  unmapped <- depth_meta |>
    filter(
      is.na(
        patient_id
      ) |
      patient_id==""
    )

  write_excel_csv(
    unmapped,
    file.path(
      OUT,
      "V2_STEP86A_unmapped_selected_runs.csv"
    ),
    na=""
  )

  if (nrow(unmapped)>0) {
    stop(
      paste0(
        "Longitudinal metadata map incomplete: ",
        nrow(unmapped),
        " selected Runs unmapped."
      )
    )
  }

  ck("796 SELECTED RUNS MAPPED TO FROZEN MASTER")

  # ==========================================================
  # 5. Define an analysis timepoint key
  # ==========================================================

  depth_meta <- depth_meta |>
    mutate(
      timepoint_key=case_when(
        !is.na(time_day) ~
          paste0(
            "DAY_",
            format(
              time_day,
              trim=TRUE,
              scientific=FALSE
            )
          ),

        !is.na(time_raw) &
        time_raw!="" ~
          paste0(
            "RAW_",
            time_raw
          ),

        TRUE ~
          paste0(
            "RUN_",
            run_id
          )
      )
    )

  write_excel_csv(
    depth_meta,
    file.path(
      OUT,
      "V2_STEP86A_selected_run_patient_time_map.csv"
    ),
    na=""
  )

  # ==========================================================
  # 6. Candidate depth-threshold longitudinal impact
  # ==========================================================

  impact_rows <- list()
  excluded_rows <- list()

  for (project in registry$project) {

    d <- depth_meta |>
      filter(
        .data$project==project
      )

    for (thr in thresholds) {

      kept <- d |>
        filter(
          reads_after_taxonomy_preview>=thr
        )

      excluded <- d |>
        filter(
          reads_after_taxonomy_preview<thr
        )

      before_patient_time <- d |>
        distinct(
          patient_id,
          timepoint_key
        ) |>
        count(
          patient_id,
          name="n_timepoints"
        )

      after_patient_time <- kept |>
        distinct(
          patient_id,
          timepoint_key
        ) |>
        count(
          patient_id,
          name="n_timepoints"
        )

      impact_rows[[
        length(impact_rows)+1
      ]] <- tibble(
        project=project,
        threshold=thr,

        samples_before=nrow(d),
        samples_after=nrow(kept),
        samples_removed=nrow(excluded),

        patients_before=
          n_distinct(
            d$patient_id
          ),

        patients_after=
          n_distinct(
            kept$patient_id
          ),

        patients_GE2_before=
          sum(
            before_patient_time$n_timepoints>=2
          ),

        patients_GE2_after=
          sum(
            after_patient_time$n_timepoints>=2
          ),

        patients_GE3_before=
          sum(
            before_patient_time$n_timepoints>=3
          ),

        patients_GE3_after=
          sum(
            after_patient_time$n_timepoints>=3
          )
      )

      if (nrow(excluded)) {

        excluded_rows[[
          length(excluded_rows)+1
        ]] <- excluded |>
          mutate(
            threshold=thr
          )
      }
    }
  }

  impact <- bind_rows(
    impact_rows
  ) |>
    mutate(
      patients_lost=
        patients_before-
        patients_after,

      GE2_lost=
        patients_GE2_before-
        patients_GE2_after,

      GE3_lost=
        patients_GE3_before-
        patients_GE3_after
    )

  excluded_candidates <- bind_rows(
    excluded_rows
  )

  write_excel_csv(
    impact,
    file.path(
      OUT,
      "V2_STEP86A_depth_threshold_longitudinal_impact_BY_PROJECT.csv"
    ),
    na=""
  )

  write_excel_csv(
    excluded_candidates,
    file.path(
      OUT,
      "V2_STEP86A_candidate_excluded_runs_ALL_THRESHOLDS.csv"
    ),
    na=""
  )

  global_impact <- impact |>
    group_by(
      threshold
    ) |>
    summarise(
      samples_before=sum(
        samples_before
      ),
      samples_after=sum(
        samples_after
      ),
      samples_removed=sum(
        samples_removed
      ),
      patients_before=sum(
        patients_before
      ),
      patients_after=sum(
        patients_after
      ),
      patients_lost=sum(
        patients_lost
      ),
      GE2_before=sum(
        patients_GE2_before
      ),
      GE2_after=sum(
        patients_GE2_after
      ),
      GE2_lost=sum(
        GE2_lost
      ),
      GE3_before=sum(
        patients_GE3_before
      ),
      GE3_after=sum(
        patients_GE3_after
      ),
      GE3_lost=sum(
        GE3_lost
      ),
      .groups="drop"
    )

  write_excel_csv(
    global_impact,
    file.path(
      OUT,
      "V2_STEP86A_depth_threshold_longitudinal_impact_GLOBAL.csv"
    ),
    na=""
  )

  # ==========================================================
  # 7. Freeze rule documentation, but NOT a depth cutoff
  # ==========================================================

  rule <- tibble(
    rule_component=c(
      "Taxonomy_database",
      "Taxonomy_keep",
      "Taxonomy_remove",
      "Kingdom_unassigned",
      "Sample_depth_cutoff"
    ),

    frozen_value=c(
      "SILVA 138.2 SSU NR99",
      "Bacteria + Archaea",
      "Chloroplast/Mitochondria + other Kingdoms",
      "REMOVE",
      "NOT_YET_FROZEN_PENDING_STEP86A_REVIEW"
    )
  )

  write_excel_csv(
    rule,
    file.path(
      OUT,
      "V2_STEP85C86A_rule_freeze_status.csv"
    ),
    na=""
  )

  readme <- c(
    "SEPSIS V2 - STEP85C / 86A",
    paste0(
      "Created: ",
      Sys.time()
    ),
    "",
    "STEP85C TAXONOMY FEATURE RULE FROZEN",
    "Keep: SILVA Kingdom Bacteria or Archaea.",
    "Remove: chloroplast/mitochondria annotations, other kingdoms, and Kingdom-unassigned features.",
    "",
    "GLOBAL STEP85B2 EVIDENCE",
    "36,681 ASVs before -> 36,004 retained.",
    "66,428,902 reads before -> 66,413,963 retained.",
    "No sample became zero after taxonomy filtering.",
    "",
    "IMPORTANT OUTLIER",
    "PRJNA516701 SRR8488562 loses >80% of reads after target-taxonomy filtering and is explicitly flagged.",
    "",
    "STEP86A",
    "Candidate minimum library-size thresholds are audited only.",
    "No depth threshold is applied by this script.",
    "",
    "NEXT",
    "Review GLOBAL and BY_PROJECT longitudinal impact tables, then freeze one depth threshold for all eight cohorts."
  )

  writeLines(
    readme,
    file.path(
      OUT,
      "README_STEP85C86A.txt"
    ),
    useBytes=TRUE
  )

  writeLines(
    c(
      paste0(
        "Completed: ",
        Sys.time()
      ),
      "Status: STEP85C86A COMPLETE",
      "Taxonomy feature rule frozen.",
      "No sample-depth threshold applied."
    ),
    file.path(
      OUT,
      "_STEP85C86A_COMPLETE.ok"
    )
  )

  ck("STEP85C86A COMPLETE")

  cat("\n============================================================\n")
  cat("SEPSIS V2 - STEP85C / 86A COMPLETE\n")
  cat("============================================================\n\n")

  print(
    freeze_registry,
    n=Inf,
    width=Inf
  )

  cat("\nGLOBAL DEPTH THRESHOLD IMPACT:\n")

  print(
    global_impact,
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
        "STEP85C86A FATAL ERROR: ",
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
      "STEP85C86A FAILED"
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
