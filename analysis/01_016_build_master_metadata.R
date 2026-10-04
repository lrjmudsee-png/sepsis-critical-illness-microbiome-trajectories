# ============================================================
# Sepsis V2 - Step 72
# Build validated public-cohort master metadata
# R 4.4.0 / Windows
#
# This step consolidates the cohorts whose sample-patient-time
# mapping has now been explicitly validated:
#
#   PRJNA516701
#   PRJNA578267
#   PRJNA595346
#   PRJNA884103
#   PRJEB67798
#   PRJEB82425 (gut/rectal only)
#
# PRJEB37289 is retained as STATIC_BACKGROUND and is NOT mixed
# into the longitudinal master table.
#
# This script does not modify source files.
# ============================================================

options(stringsAsFactors = FALSE)

pkgs <- c(
  "readr",
  "dplyr",
  "tidyr",
  "stringr",
  "tibble",
  "purrr"
)

missing <- pkgs[
  !vapply(
    pkgs,
    requireNamespace,
    logical(1),
    quietly = TRUE
  )
]

if (length(missing) > 0) {
  install.packages(
    missing,
    repos = "https://cloud.r-project.org"
  )
}

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(tibble)
  library(purrr)
})

# ------------------------------------------------------------
# Paths
# ------------------------------------------------------------
ROOT68 <- "E:/sepsis_project/results/V2_08_remaining_metadata_resolution"
ROOT70 <- "E:/sepsis_project/results/V2_10_exact_metadata_mapping"
ROOT71B <- "E:/sepsis_project/results/V2_11_PRJEB82425_body_site_final"

OUT_ROOT <- "E:/sepsis_project/results/V2_12_master_metadata"
dir.create(
  OUT_ROOT,
  recursive = TRUE,
  showWarnings = FALSE
)

# ------------------------------------------------------------
# Helpers
# ------------------------------------------------------------
clean_chr <- function(x) {
  x <- trimws(as.character(x))
  x[
    is.na(x) |
    x == "" |
    tolower(x) %in% c(
      "na",
      "nan",
      "n/a",
      "null",
      "none"
    )
  ] <- NA_character_
  x
}

find_latest <- function(folder, pattern) {
  if (!dir.exists(folder)) return(NA_character_)

  x <- list.files(
    folder,
    pattern = pattern,
    full.names = TRUE
  )

  if (length(x) == 0) return(NA_character_)

  x[
    which.max(
      file.info(x)$mtime
    )
  ]
}

safe_csv <- function(path) {
  if (
    length(path) == 0 ||
    is.na(path) ||
    !file.exists(path)
  ) return(NULL)

  suppressMessages(
    read_csv(
      path,
      show_col_types = FALSE,
      progress = FALSE,
      name_repair = "unique"
    )
  )
}

get1 <- function(df, candidates) {
  if (is.null(df)) return(character(0))

  hit <- candidates[
    candidates %in% names(df)
  ]

  if (length(hit) == 0) {
    return(
      rep(
        NA_character_,
        nrow(df)
      )
    )
  }

  clean_chr(
    df[[hit[1]]]
  )
}

get_num <- function(df, candidates) {
  suppressWarnings(
    as.numeric(
      get1(
        df,
        candidates
      )
    )
  )
}

standardize <- function(
  df,
  project,
  cohort_role,
  analysis_tier
) {

  if (is.null(df) || nrow(df) == 0) {
    return(tibble())
  }

  tibble(
    project = project,

    patient_id = get1(
      df,
      c(
        "Patient_ID",
        "patient_id"
      )
    ),

    sample_id = get1(
      df,
      c(
        "Sample_ID",
        "sample_id"
      )
    ),

    biosample = get1(
      df,
      c(
        "BioSample",
        "biosample"
      )
    ),

    run_id = get1(
      df,
      c(
        "Run_ID",
        "run_id"
      )
    ),

    experiment_id = get1(
      df,
      c(
        "Experiment_ID",
        "experiment_id"
      )
    ),

    time_raw = get1(
      df,
      c(
        "Time_Raw",
        "time_raw"
      )
    ),

    time_day = get_num(
      df,
      c(
        "Time_Day",
        "time_day"
      )
    ),

    time_class = get1(
      df,
      c(
        "Time_Class",
        "Time_Label",
        "Study_Period",
        "time_class"
      )
    ),

    collection_date = get1(
      df,
      c(
        "Collection_Date",
        "Collection_Date_Raw",
        "Collection_Date_Shifted",
        "collection_date"
      )
    ),

    body_site = get1(
      df,
      c(
        "Body_Site",
        "body_site"
      )
    ),

    infection_group = get1(
      df,
      c(
        "Infection_Group",
        "infection_group"
      )
    ),

    outcome = get1(
      df,
      c(
        "Outcome",
        "outcome"
      )
    ),

    antibiotics = get1(
      df,
      c(
        "ABx_score_sample",
        "Abx_days",
        "antibiotics"
      )
    ),

    severity = get1(
      df,
      c(
        "SOFA",
        "APACHE",
        "severity"
      )
    ),

    sex = get1(
      df,
      c(
        "Sex",
        "sex"
      )
    ),

    age = get1(
      df,
      c(
        "Age",
        "age"
      )
    ),

    gestational_age = get1(
      df,
      c(
        "Gestational_age",
        "gestational_age"
      )
    ),

    birthweight = get1(
      df,
      c(
        "Birthweight",
        "birthweight"
      )
    ),

    cohort_role = cohort_role,

    analysis_tier = analysis_tier
  )
}

# ------------------------------------------------------------
# 1. Step68 READY map
# ------------------------------------------------------------
ready68_path <- find_latest(
  ROOT68,
  "^V2_READY_map_after68_.*\\.csv$"
)

ready68 <- safe_csv(
  ready68_path
)

if (is.null(ready68)) {
  stop(
    "Cannot find V2_READY_map_after68_*.csv"
  )
}

# ------------------------------------------------------------
# 2. PRJEB67798 exact map
# ------------------------------------------------------------
map677_path <- file.path(
  ROOT70,
  "PRJEB67798_exact_sample_patient_run_map.csv"
)

map677 <- safe_csv(
  map677_path
)

if (is.null(map677)) {
  stop(
    "Cannot find PRJEB67798_exact_sample_patient_run_map.csv"
  )
}

# ------------------------------------------------------------
# 3. PRJEB82425 gut-only final map
# ------------------------------------------------------------
map824_path <- file.path(
  ROOT71B,
  "PRJEB82425_GUT_ONLY_FINAL.csv"
)

map824 <- safe_csv(
  map824_path
)

if (is.null(map824)) {
  stop(
    "Cannot find PRJEB82425_GUT_ONLY_FINAL.csv"
  )
}

# ------------------------------------------------------------
# 4. Standardize cohort by cohort
# ------------------------------------------------------------
parts <- list()

# PRJNA516701
x <- ready68 |>
  filter(Project == "PRJNA516701")

parts[["PRJNA516701"]] <- standardize(
  x,
  project = "PRJNA516701",
  cohort_role = "CRITICAL_ILLNESS_LONGITUDINAL",
  analysis_tier = "SUPPORTING_LONGITUDINAL"
)

# PRJNA578267
x <- ready68 |>
  filter(Project == "PRJNA578267")

parts[["PRJNA578267"]] <- standardize(
  x,
  project = "PRJNA578267",
  cohort_role = "NONSEPSIS_SURGICAL_LONGITUDINAL_CONTROL",
  analysis_tier = "CONTROL_LONGITUDINAL"
)

# PRJNA595346
x <- ready68 |>
  filter(
    Project == "PRJNA595346"
  ) |>
  filter(
    is.na(Body_Site) |
    str_detect(
      Body_Site,
      regex(
        "rectal|stool|fecal|faecal|gut",
        ignore_case = TRUE
      )
    )
  )

parts[["PRJNA595346"]] <- standardize(
  x,
  project = "PRJNA595346",
  cohort_role = "ICU_GUT_LONGITUDINAL_EXTERNAL_VALIDATION",
  analysis_tier = "EXTERNAL_VALIDATION"
)

# PRJNA884103
x <- ready68 |>
  filter(Project == "PRJNA884103")

parts[["PRJNA884103"]] <- standardize(
  x,
  project = "PRJNA884103",
  cohort_role = "NEONATAL_BSI_LONGITUDINAL_EXTERNAL_VALIDATION",
  analysis_tier = "EXTERNAL_VALIDATION_SHOTGUN"
)

# PRJEB67798
parts[["PRJEB67798"]] <- standardize(
  map677,
  project = "PRJEB67798",
  cohort_role = "NONSEPSIS_SURGICAL_LONGITUDINAL_CONTROL",
  analysis_tier = "CONTROL_LONGITUDINAL"
)

# PRJEB82425 rectal only
parts[["PRJEB82425"]] <- standardize(
  map824,
  project = "PRJEB82425",
  cohort_role = "ICU_INFECTION_LONGITUDINAL_EXTERNAL_VALIDATION",
  analysis_tier = "EXTERNAL_VALIDATION"
)

master <- bind_rows(
  parts
) |>
  mutate(
    patient_id = clean_chr(patient_id),
    sample_id = clean_chr(sample_id),
    run_id = clean_chr(run_id),
    time_raw = clean_chr(time_raw),
    body_site = clean_chr(body_site),

    patient_uid = ifelse(
      is.na(patient_id),
      NA_character_,
      paste(
        project,
        patient_id,
        sep = "::"
      )
    ),

    sample_uid = ifelse(
      is.na(sample_id),
      NA_character_,
      paste(
        project,
        sample_id,
        sep = "::"
      )
    ),

    run_uid = ifelse(
      is.na(run_id),
      NA_character_,
      paste(
        project,
        run_id,
        sep = "::"
      )
    )
  )

# ------------------------------------------------------------
# 5. Longitudinal eligibility
# ------------------------------------------------------------
patient_time <- master |>
  filter(
    !is.na(patient_uid),
    !is.na(time_raw)
  ) |>
  distinct(
    patient_uid,
    time_raw
  ) |>
  count(
    patient_uid,
    name = "n_distinct_timepoints"
  )

master <- master |>
  left_join(
    patient_time,
    by = "patient_uid"
  ) |>
  mutate(
    n_distinct_timepoints = ifelse(
      is.na(n_distinct_timepoints),
      0L,
      n_distinct_timepoints
    ),

    longitudinal_ge2 =
      n_distinct_timepoints >= 2,

    longitudinal_ge3 =
      n_distinct_timepoints >= 3
  )

# ------------------------------------------------------------
# 6. Duplicate QC
# ------------------------------------------------------------
duplicate_run <- master |>
  filter(!is.na(run_uid)) |>
  count(
    run_uid,
    name = "n"
  ) |>
  filter(n > 1)

duplicate_sample <- master |>
  filter(!is.na(sample_uid)) |>
  count(
    sample_uid,
    name = "n"
  ) |>
  filter(n > 1)

duplicate_patient_sample <- master |>
  filter(
    !is.na(patient_uid),
    !is.na(sample_uid)
  ) |>
  count(
    patient_uid,
    sample_uid,
    name = "n"
  ) |>
  filter(n > 1)

# ------------------------------------------------------------
# 7. Cohort-level QC
# ------------------------------------------------------------
cohort_qc <- master |>
  group_by(
    project,
    cohort_role,
    analysis_tier
  ) |>
  summarise(
    rows = n(),

    patients = n_distinct(
      patient_uid,
      na.rm = TRUE
    ),

    samples = n_distinct(
      sample_uid,
      na.rm = TRUE
    ),

    runs = n_distinct(
      run_uid,
      na.rm = TRUE
    ),

    missing_patient = sum(
      is.na(patient_id)
    ),

    missing_sample = sum(
      is.na(sample_id)
    ),

    missing_run = sum(
      is.na(run_id)
    ),

    missing_time = sum(
      is.na(time_raw)
    ),

    patients_ge2_timepoints =
      n_distinct(
        patient_uid[
          longitudinal_ge2
        ],
        na.rm = TRUE
      ),

    patients_ge3_timepoints =
      n_distinct(
        patient_uid[
          longitudinal_ge3
        ],
        na.rm = TRUE
      ),

    .groups = "drop"
  )

# ------------------------------------------------------------
# 8. Analysis-ready patient table
# ------------------------------------------------------------
patient_summary <- master |>
  filter(
    !is.na(patient_uid)
  ) |>
  group_by(
    project,
    patient_uid,
    patient_id,
    cohort_role,
    analysis_tier
  ) |>
  summarise(
    n_rows = n(),
    n_samples = n_distinct(
      sample_uid,
      na.rm = TRUE
    ),
    n_runs = n_distinct(
      run_uid,
      na.rm = TRUE
    ),
    n_timepoints = n_distinct(
      time_raw[
        !is.na(time_raw)
      ]
    ),
    longitudinal_ge2 =
      n_timepoints >= 2,
    longitudinal_ge3 =
      n_timepoints >= 3,
    .groups = "drop"
  )

# ------------------------------------------------------------
# 9. Save
# ------------------------------------------------------------
stamp <- format(
  Sys.time(),
  "%Y%m%d_%H%M%S"
)

master_path <- file.path(
  OUT_ROOT,
  paste0(
    "V2_master_metadata_",
    stamp,
    ".csv"
  )
)

qc_path <- file.path(
  OUT_ROOT,
  paste0(
    "V2_master_metadata_QC_",
    stamp,
    ".csv"
  )
)

patient_path <- file.path(
  OUT_ROOT,
  paste0(
    "V2_master_patient_summary_",
    stamp,
    ".csv"
  )
)

dup_run_path <- file.path(
  OUT_ROOT,
  paste0(
    "V2_master_duplicate_run_QC_",
    stamp,
    ".csv"
  )
)

dup_sample_path <- file.path(
  OUT_ROOT,
  paste0(
    "V2_master_duplicate_sample_QC_",
    stamp,
    ".csv"
  )
)

dup_ps_path <- file.path(
  OUT_ROOT,
  paste0(
    "V2_master_duplicate_patient_sample_QC_",
    stamp,
    ".csv"
  )
)

write_excel_csv(
  master,
  master_path,
  na = ""
)

write_excel_csv(
  cohort_qc,
  qc_path,
  na = ""
)

write_excel_csv(
  patient_summary,
  patient_path,
  na = ""
)

write_excel_csv(
  duplicate_run,
  dup_run_path,
  na = ""
)

write_excel_csv(
  duplicate_sample,
  dup_sample_path,
  na = ""
)

write_excel_csv(
  duplicate_patient_sample,
  dup_ps_path,
  na = ""
)

# ------------------------------------------------------------
# 10. Console
# ------------------------------------------------------------
cat(
  "\n============================================================\n"
)

cat(
  "SEPSIS V2 - STEP 72 COMPLETE\n"
)

cat(
  "R version: ",
  R.version.string,
  "\n",
  sep = ""
)

cat(
  "============================================================\n\n"
)

cat(
  "MASTER COHORT QC:\n"
)

print(
  cohort_qc,
  n = Inf,
  width = Inf
)

cat(
  "\nDuplicate QC:\n"
)

cat(
  "duplicate run IDs: ",
  nrow(duplicate_run),
  "\n",
  sep = ""
)

cat(
  "duplicate sample IDs: ",
  nrow(duplicate_sample),
  "\n",
  sep = ""
)

cat(
  "duplicate patient-sample pairs: ",
  nrow(duplicate_patient_sample),
  "\n",
  sep = ""
)

cat(
  "\nOutputs:\n",
  OUT_ROOT,
  "\n",
  sep = ""
)

cat(
  "============================================================\n"
)
