# ============================================================
# Sepsis V2 - Step 87A
# Rebuild FULL analysis metadata + time-axis/replicate audit
#
# PURPOSE
#   Step86B produced valid final sequence tables (786 samples),
#   but its compact metadata contains mainly Run/patient/time/read-QC
#   fields. Before formal longitudinal statistics, restore the full
#   Step81D frozen clinical metadata and audit duplicate patient-time
#   observations.
#
# IMPORTANT
#   - NO broad filesystem scan.
#   - NO sample/ASV deletion.
#   - NO alpha/beta statistics yet.
#   - Reads only canonical/frozen known paths.
#   - Preserves both members of any same-patient/same-time replicate.
#
# Expected current state:
#   8 projects
#   786 samples
#   390 project-specific patients
#   258 patients >=2 timepoints
#   132 patients >=3 timepoints
#
# Known item to audit:
#   PRJNA516701 patient 4200 has two retained Runs at DAY_5:
#     SRR8488665
#     SRR8488668
#   This script flags the pair but DOES NOT resolve/drop it.
#
# R 4.4.0 / Windows
# ============================================================

options(stringsAsFactors = FALSE)

pkgs <- c("readr","dplyr","stringr","tibble","purrr","tidyr")
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
  library(tidyr)
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

STEP86B <- file.path(
  ROOT,
  "results",
  "V2_26B_DEPTH2000_FINAL_ANALYSIS_READY"
)

OUT <- file.path(
  ROOT,
  "results",
  "V2_27A_ANALYSIS_METADATA_REBUILD_AND_REPLICATE_AUDIT"
)

dir.create(
  OUT,
  recursive=TRUE,
  showWarnings=FALSE
)

LOG <- file.path(
  OUT,
  "_STEP87A_runtime_checkpoints.txt"
)

ERR <- file.path(
  OUT,
  "_STEP87A_FATAL_ERROR.txt"
)

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

pick_col <- function(df, candidates, required=TRUE) {

  hit <- candidates[
    candidates %in% names(df)
  ]

  if (length(hit)) {
    return(hit[1])
  }

  lower_names <- tolower(names(df))
  idx <- match(
    tolower(candidates),
    lower_names
  )

  idx <- idx[!is.na(idx)]

  if (length(idx)) {
    return(names(df)[idx[1]])
  }

  if (required) {
    stop(
      "Required column not found. Candidates: ",
      paste(candidates, collapse=", ")
    )
  }

  NA_character_
}

nonmissing <- function(x) {
  !(is.na(x) | trimws(as.character(x))=="")
}

# ------------------------------------------------------------
# Frozen project roles
# ------------------------------------------------------------

role_map <- tribble(
  ~project, ~analysis_role,

  "PRJNA691455",
  "CORE_SEPSIS_LONGITUDINAL",

  "PRJEB82425",
  "ICU_INFECTION_LONGITUDINAL_EXTERNAL",

  "PRJNA516701",
  "ICU_BACKGROUND_LONGITUDINAL",

  "PRJNA851469",
  "ICU_BACKGROUND_LONGITUDINAL",

  "PRJNA578267",
  "NONSEPSIS_LONGITUDINAL_CONTROL",

  "PRJNA430161",
  "INTERVENTION_LONGITUDINAL_SUPPORT",

  "PRJNA1166732",
  "INTERVENTION_LONGITUDINAL_SUPPORT",

  "PRJNA978257",
  "STATIC_SUPPORT"
)

# ------------------------------------------------------------
# Locate Step81D master in canonical 00_FREEZE only
# ------------------------------------------------------------

master_candidates <- list.files(
  FREEZE_DIR,
  pattern="^V2_FINAL_FROZEN_master_metadata_CLINICAL_REPAIRED.*\\.csv$",
  full.names=TRUE,
  recursive=FALSE
)

if (!length(master_candidates)) {
  stop(
    "No Step81D final frozen master found in canonical 00_FREEZE."
  )
}

if (length(master_candidates)>1) {

  info <- file.info(
    master_candidates
  )

  master_path <- rownames(info)[
    which.max(info$mtime)
  ]

} else {

  master_path <- master_candidates[1]
}

main <- function() {

  ck("STEP87A STARTED")

  # ==========================================================
  # 1. Load Step86B final 786-sample metadata
  # ==========================================================

  final_meta_path <- file.path(
    STEP86B,
    "V2_STEP86B_FINAL_analysis_ready_metadata_786samples.csv"
  )

  final_reg_path <- file.path(
    STEP86B,
    "V2_STEP86B_FINAL_analysis_ready_registry.csv"
  )

  if (
    !file.exists(final_meta_path) ||
    !file.exists(final_reg_path)
  ) {
    stop(
      "Step86B final metadata/registry missing."
    )
  }

  final_meta <- safe_csv(
    final_meta_path
  )

  final_reg <- safe_csv(
    final_reg_path
  )

  if (
    nrow(final_meta)!=786 ||
    n_distinct(final_meta$project)!=8 ||
    anyDuplicated(
      final_meta[c("project","run_id")]
    )
  ) {
    stop(
      "Step86B 786-sample metadata hard guard failed."
    )
  }

  if (
    nrow(final_reg)!=8 ||
    sum(final_reg$samples_final)!=786 ||
    sum(final_reg$patients_final)!=390 ||
    sum(final_reg$patients_GE2_final)!=258 ||
    sum(final_reg$patients_GE3_final)!=132
  ) {
    stop(
      "Step86B final registry hard guard failed."
    )
  }

  ck("STEP86B FINAL STATE GUARDED")

  # ==========================================================
  # 2. Load full Step81D frozen master
  # ==========================================================

  master <- safe_csv(
    master_path
  )

  project_col <- pick_col(
    master,
    c("project","Project","study","Study")
  )

  run_col <- pick_col(
    master,
    c("run_id","Run","run","Run_ID","run_accession")
  )

  # Standardize the join keys only; retain every original field.
  master_std <- master |>
    mutate(
      .join_project=
        as.character(.data[[project_col]]),

      .join_run=
        as.character(.data[[run_col]])
    )

  # Conflicting duplicate Run mappings are not allowed.
  dup_master <- master_std |>
    filter(
      nonmissing(.join_run)
    ) |>
    count(
      .join_project,
      .join_run,
      name="n"
    ) |>
    filter(
      n>1
    )

  if (nrow(dup_master)>0) {
    write_excel_csv(
      dup_master,
      file.path(
        OUT,
        "V2_STEP87A_FATAL_duplicate_runs_in_frozen_master.csv"
      ),
      na=""
    )

    stop(
      "Frozen master contains duplicated project+Run keys; refusing to guess."
    )
  }

  # ==========================================================
  # 3. Join full clinical master onto final 786 samples
  # ==========================================================

  # Prefix master fields that would collide with Step86B fields.
  key_cols <- c(
    project_col,
    run_col,
    ".join_project",
    ".join_run"
  )

  collision <- intersect(
    names(master_std),
    names(final_meta)
  )

  collision <- setdiff(
    collision,
    key_cols
  )

  master_join <- master_std

  if (length(collision)) {

    rename_map <- setNames(
      paste0(
        "master_",
        collision
      ),
      collision
    )

    master_join <- master_join |>
      rename(
        !!!rename_map
      )
  }

  full <- final_meta |>
    mutate(
      .join_project=as.character(project),
      .join_run=as.character(run_id)
    ) |>
    left_join(
      master_join,
      by=c(
        ".join_project",
        ".join_run"
      )
    )

  if (nrow(full)!=786) {
    stop(
      "Full metadata join changed row count; possible row multiplication."
    )
  }

  # A matched master row must carry at least its original project/run.
  master_project_joined <- if (
    paste0("master_",project_col) %in% names(full)
  ) {
    paste0("master_",project_col)
  } else {
    project_col
  }

  master_run_joined <- if (
    paste0("master_",run_col) %in% names(full)
  ) {
    paste0("master_",run_col)
  } else {
    run_col
  }

  unmatched <- full |>
    filter(
      is.na(
        .data[[master_run_joined]]
      ) |
      as.character(
        .data[[master_run_joined]]
      )==""
    )

  write_excel_csv(
    unmatched,
    file.path(
      OUT,
      "V2_STEP87A_unmatched_final_runs_to_frozen_master.csv"
    ),
    na=""
  )

  if (nrow(unmatched)>0) {
    stop(
      "Not all 786 final Runs map to Step81D frozen master."
    )
  }

  ck("FULL STEP81D CLINICAL METADATA RESTORED")

  # ==========================================================
  # 4. Add project role and explicit time-axis semantics
  # ==========================================================

  full <- full |>
    left_join(
      role_map,
      by="project"
    )

  if (any(is.na(full$analysis_role))) {
    stop(
      "Missing analysis_role for one or more projects."
    )
  }

  full <- full |>
    mutate(
      time_axis_type=case_when(

        project %in%
        c(
          "PRJNA430161",
          "PRJNA516701",
          "PRJNA1166732"
        ) ~
          "NUMERIC_DAY",

        project=="PRJNA691455" &
        !is.na(time_day) ~
          "NUMERIC_DAY",

        project=="PRJNA691455" &
        is.na(time_day) ~
          "STATIC_OR_ONE_TIME",

        project=="PRJNA851469" &
        !is.na(time_day) ~
          "NUMERIC_DAY",

        project=="PRJNA851469" &
        is.na(time_day) ~
          "STATIC_CONTROL",

        project=="PRJNA578267" ~
          "ORDINAL_VISIT",

        project=="PRJEB82425" ~
          "EVENT_STAGE",

        project=="PRJNA978257" ~
          "STATIC",

        TRUE ~
          "UNRESOLVED"
      ),

      analysis_time_order=case_when(

        !is.na(time_day) ~
          as.numeric(time_day),

        project=="PRJNA578267" &
        time_raw=="Time first" ~
          1,

        project=="PRJNA578267" &
        time_raw=="Time second" ~
          2,

        project=="PRJNA578267" &
        time_raw=="Time third" ~
          3,

        # Event order for plotting/model factor ordering only;
        # NOT treated as continuous equal-interval time.
        project=="PRJEB82425" &
        time_raw=="Inclusion" ~
          1,

        project=="PRJEB82425" &
        time_raw=="Infection_D1" ~
          2,

        project=="PRJEB82425" &
        time_raw=="Infection_D5" ~
          3,

        project=="PRJEB82425" &
        time_raw=="Extubation" ~
          4,

        project=="PRJEB82425" &
        time_raw=="Discharge" ~
          5,

        TRUE ~
          NA_real_
      ),

      analysis_time_is_continuous=
        time_axis_type=="NUMERIC_DAY"
    )

  # ==========================================================
  # 5. Duplicate patient-timepoint audit
  # ==========================================================

  dup_groups <- full |>
    count(
      project,
      patient_id,
      timepoint_key,
      name="n_runs_same_patient_time"
    ) |>
    filter(
      n_runs_same_patient_time>1
    )

  dup_detail <- full |>
    inner_join(
      dup_groups,
      by=c(
        "project",
        "patient_id",
        "timepoint_key"
      )
    ) |>
    arrange(
      project,
      patient_id,
      timepoint_key,
      desc(
        reads_after_taxonomy_preview
      )
    )

  write_excel_csv(
    dup_groups,
    file.path(
      OUT,
      "V2_STEP87A_duplicate_patient_time_groups.csv"
    ),
    na=""
  )

  write_excel_csv(
    dup_detail,
    file.path(
      OUT,
      "V2_STEP87A_duplicate_patient_time_detail.csv"
    ),
    na=""
  )

  # Current expected duplicate pair.
  if (
    nrow(dup_groups)!=1 ||
    dup_groups$project[1]!="PRJNA516701" ||
    as.character(
      dup_groups$patient_id[1]
    )!="4200" ||
    dup_groups$timepoint_key[1]!="DAY_5" ||
    dup_groups$n_runs_same_patient_time[1]!=2
  ) {
    stop(
      "Duplicate patient-time audit differs from expected current state."
    )
  }

  expected_dup_runs <- c(
    "SRR8488665",
    "SRR8488668"
  )

  if (
    !setequal(
      dup_detail$run_id,
      expected_dup_runs
    )
  ) {
    stop(
      "PRJNA516701 patient 4200 duplicate Run IDs differ from expected."
    )
  }

  full <- full |>
    left_join(
      dup_groups,
      by=c(
        "project",
        "patient_id",
        "timepoint_key"
      )
    ) |>
    mutate(
      n_runs_same_patient_time=
        coalesce(
          n_runs_same_patient_time,
          1L
        ),

      duplicate_patient_time=
        n_runs_same_patient_time>1,

      replicate_resolution_status=case_when(
        duplicate_patient_time ~
          "UNRESOLVED_RETAIN_BOTH_PENDING_SOURCE_REVIEW",

        TRUE ~
          "NOT_DUPLICATE"
      )
    )

  ck("PATIENT-TIME REPLICATE AUDIT COMPLETE")

  # ==========================================================
  # 6. Timepoint / longitudinal eligibility summaries
  # ==========================================================

  project_time_summary <- full |>
    count(
      project,
      analysis_role,
      time_axis_type,
      time_raw,
      time_day,
      analysis_time_order,
      name="samples"
    ) |>
    arrange(
      project,
      analysis_time_order,
      time_raw
    )

  write_excel_csv(
    project_time_summary,
    file.path(
      OUT,
      "V2_STEP87A_project_timepoint_distribution.csv"
    ),
    na=""
  )

  patient_time <- full |>
    distinct(
      project,
      patient_id,
      timepoint_key
    ) |>
    count(
      project,
      patient_id,
      name="n_distinct_timepoints"
    )

  eligibility <- patient_time |>
    group_by(
      project
    ) |>
    summarise(
      patients=n(),
      patients_GE2=sum(
        n_distinct_timepoints>=2
      ),
      patients_GE3=sum(
        n_distinct_timepoints>=3
      ),
      .groups="drop"
    ) |>
    left_join(
      role_map,
      by="project"
    )

  if (
    sum(eligibility$patients)!=390 ||
    sum(eligibility$patients_GE2)!=258 ||
    sum(eligibility$patients_GE3)!=132
  ) {
    stop(
      "Rebuilt longitudinal eligibility differs from Step86B freeze."
    )
  }

  write_excel_csv(
    eligibility,
    file.path(
      OUT,
      "V2_STEP87A_longitudinal_eligibility_by_project.csv"
    ),
    na=""
  )

  # ==========================================================
  # 7. Full-field completeness audit by project
  # ==========================================================

  # Exclude helper/join keys from field-completeness output.
  audit_cols <- setdiff(
    names(full),
    c(
      ".join_project",
      ".join_run"
    )
  )

  completeness <- map_dfr(
    unique(full$project),
    function(proj) {

      d <- full |>
        filter(
          .data$project==.env$proj
        )

      map_dfr(
        audit_cols,
        function(col) {

          x <- d[[col]]

          tibble(
            project=proj,
            field=col,
            nonmissing_n=sum(
              nonmissing(x)
            ),
            total_n=length(x),
            completeness_pct=round(
              100*mean(
                nonmissing(x)
              ),
              2
            )
          )
        }
      )
    }
  )

  write_excel_csv(
    completeness,
    file.path(
      OUT,
      "V2_STEP87A_full_metadata_field_completeness_by_project.csv"
    ),
    na=""
  )

  # ==========================================================
  # 8. Save rebuilt full analysis metadata
  # ==========================================================

  full_out <- full |>
    select(
      -any_of(
        c(
          ".join_project",
          ".join_run"
        )
      )
    )

  write_excel_csv(
    full_out,
    file.path(
      OUT,
      "V2_STEP87A_FULL_analysis_metadata_786samples.csv"
    ),
    na=""
  )

  # ==========================================================
  # 9. Registry / README
  # ==========================================================

  registry <- eligibility |>
    left_join(
      final_reg |>
        select(
          project,
          samples_final,
          ASVs_final,
          reads_final,
          final_seqtab_path
        ),
      by="project"
    ) |>
    mutate(
      duplicate_patient_time_groups=
        map_int(
          project,
          ~sum(
            dup_groups$project==.x
          )
        ),

      statistical_metadata_ready=
        duplicate_patient_time_groups==0
    )

  # We deliberately mark PRJNA516701 pending until its duplicate
  # same-time samples are source-reviewed.
  write_excel_csv(
    registry,
    file.path(
      OUT,
      "V2_STEP87A_analysis_object_readiness_registry.csv"
    ),
    na=""
  )

  readme <- c(
    "SEPSIS V2 - STEP87A ANALYSIS METADATA REBUILD",
    paste0(
      "Created: ",
      Sys.time()
    ),
    "",
    "SUCCESS",
    "All 786 Step86B final samples were rejoined to the Step81D frozen master metadata.",
    "",
    "CURRENT LONGITUDINAL STATE",
    "390 project-specific patients.",
    "258 patients with >=2 distinct timepoints.",
    "132 patients with >=3 distinct timepoints.",
    "",
    "TIME AXES",
    "Numeric-day cohorts retain numeric time_day.",
    "PRJNA578267 uses ordered visits (first/second/third), not continuous time.",
    "PRJEB82425 uses event-stage order for display/factor ordering; event-stage numbers are NOT equal-interval continuous days.",
    "Static/one-time samples are explicitly labeled.",
    "",
    "REPLICATE AUDIT",
    "One same-patient/same-time group remains:",
    "PRJNA516701 patient 4200, DAY_5, Runs SRR8488665 and SRR8488668.",
    "Both are retained in Step87A. No automatic merge/drop is performed.",
    "",
    "NEXT",
    "Resolve the PRJNA516701 patient 4200 same-time pair from source metadata.",
    "After that, construct cohort-specific analysis objects and begin primary longitudinal alpha/beta diversity and within-patient change analyses.",
    "",
    "IMPORTANT",
    "Do not merge ASV identities across different 16S variable regions.",
    "Cross-cohort synthesis should use harmonized taxonomic levels and/or cohort-level effect sizes."
  )

  writeLines(
    readme,
    file.path(
      OUT,
      "README_STEP87A.txt"
    ),
    useBytes=TRUE
  )

  writeLines(
    c(
      paste0(
        "Completed: ",
        Sys.time()
      ),
      "Status: STEP87A COMPLETE",
      "Full Step81D metadata restored for 786 final samples.",
      "One duplicate patient-time group remains unresolved: PRJNA516701 patient 4200 DAY_5."
    ),
    file.path(
      OUT,
      "_STEP87A_COMPLETE.ok"
    )
  )

  ck("STEP87A COMPLETE")

  cat("\n============================================================\n")
  cat("SEPSIS V2 - STEP87A COMPLETE\n")
  cat("============================================================\n\n")

  cat("LONGITUDINAL ELIGIBILITY:\n")
  print(
    eligibility,
    n=Inf,
    width=Inf
  )

  cat("\nDUPLICATE PATIENT-TIME GROUPS:\n")
  print(
    dup_groups,
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
        "STEP87A FATAL ERROR: ",
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
      "STEP87A FAILED"
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
