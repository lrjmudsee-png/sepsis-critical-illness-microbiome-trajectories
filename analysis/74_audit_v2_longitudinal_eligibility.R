# ============================================================
# Sepsis V2 - Step 74
# Longitudinal eligibility and analysis-module audit
# R 4.4.0 / Windows
#
# Input:
#   latest V2_FULL_master_metadata_*.csv from Step73
#
# Outputs:
#   - project-level longitudinal eligibility
#   - patient-level longitudinal eligibility
#   - analysis-module assignment
#   - clinical covariate coverage
#   - timepoint dictionary candidates
#   - project directories present in data/ but absent from master
#
# This script does NOT modify metadata or sequencing files.
# ============================================================

options(stringsAsFactors = FALSE)

pkgs <- c(
  "readr",
  "dplyr",
  "tidyr",
  "stringr",
  "purrr",
  "tibble"
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
  library(purrr)
  library(tibble)
})

# ------------------------------------------------------------
# Paths
# ------------------------------------------------------------
PROJECT_ROOT <- "E:/sepsis_project"
DATA_ROOT <- file.path(PROJECT_ROOT, "data")
STEP73_ROOT <- file.path(
  PROJECT_ROOT,
  "results",
  "V2_13_full_master_metadata"
)

OUT_ROOT <- file.path(
  PROJECT_ROOT,
  "results",
  "V2_14_longitudinal_eligibility"
)

dir.create(
  OUT_ROOT,
  recursive = TRUE,
  showWarnings = FALSE
)

# ------------------------------------------------------------
# Helpers
# ------------------------------------------------------------
find_latest <- function(folder, pattern) {
  x <- list.files(
    folder,
    pattern = pattern,
    full.names = TRUE
  )

  if (length(x) == 0) {
    return(NA_character_)
  }

  x[
    which.max(
      file.info(x)$mtime
    )
  ]
}

clean_chr <- function(x) {
  x <- trimws(as.character(x))

  x[
    is.na(x) |
    x == "" |
    tolower(x) %in% c(
      "na","nan","n/a","null","none"
    )
  ] <- NA_character_

  x
}

collapse_values <- function(x, max_n = 12) {
  x <- sort(unique(clean_chr(x)))
  x <- x[!is.na(x)]

  if (length(x) == 0) return("")

  paste(
    head(x, max_n),
    collapse = ";"
  )
}

# ------------------------------------------------------------
# Load full master
# ------------------------------------------------------------
master_path <- find_latest(
  STEP73_ROOT,
  "^V2_FULL_master_metadata_[0-9]{8}_[0-9]{6}\\.csv$"
)

if (is.na(master_path)) {
  stop("Cannot find Step73 V2_FULL_master_metadata file.")
}

master <- suppressMessages(
  read_csv(
    master_path,
    show_col_types = FALSE,
    progress = FALSE
  )
)

message("Loaded: ", master_path)
message("Rows: ", nrow(master))

# ------------------------------------------------------------
# Analysis module assignment
# ------------------------------------------------------------
project_plan <- tibble(
  project = c(
    "PRJEB33360",
    "PRJNA691455",
    "PRJEB82425",
    "PRJNA516701",
    "PRJNA595346",
    "PRJNA578267",
    "PRJEB67798",
    "PRJNA430161",
    "PRJNA1166732",
    "PRJNA912621",
    "PRJNA797231",
    "PRJNA978257",
    "PRJNA1010969",
    "PRJNA884103"
  ),

  analysis_module = c(
    "CORE_SEPSIS_LONGITUDINAL",
    "CORE_SEPSIS_LONGITUDINAL",
    "ICU_INFECTION_LONGITUDINAL_EXTERNAL",
    "ICU_BACKGROUND_LONGITUDINAL",
    "ICU_BACKGROUND_LONGITUDINAL",
    "NONSEPSIS_LONGITUDINAL_CONTROL",
    "NONSEPSIS_LONGITUDINAL_CONTROL",
    "INTERVENTION_LONGITUDINAL_SUPPORT",
    "INTERVENTION_LONGITUDINAL_SUPPORT",
    "ORGAN_DYSFUNCTION_LONGITUDINAL_SUPPORT",
    "ORGAN_DYSFUNCTION_STATIC_SUPPORT",
    "STATIC_SEPSIS_CONTROL_SUPPORT",
    "STATIC_SEPSIS_CONTROL_SUPPORT",
    "SHOTGUN_BSI_LONGITUDINAL_EXTERNAL"
  ),

  primary_use = c(
    "Primary within-patient sepsis trajectory",
    "Primary within-patient sepsis trajectory",
    "External infection-event trajectory validation",
    "General critical-illness longitudinal background",
    "General ICU longitudinal external validation",
    "Non-sepsis surgery temporal control",
    "Non-sepsis surgery temporal control",
    "Treatment-response supportive analysis",
    "Treatment-response supportive analysis after DADA2",
    "Organ-dysfunction trajectory support",
    "Cross-sectional organ-dysfunction support",
    "Cross-sectional sepsis/control support",
    "Cross-sectional external support",
    "Shotgun neonatal BSI external validation"
  ),

  use_in_primary_16S_longitudinal = c(
    TRUE,
    TRUE,
    FALSE,
    FALSE,
    FALSE,
    FALSE,
    FALSE,
    FALSE,
    FALSE,
    FALSE,
    FALSE,
    FALSE,
    FALSE,
    FALSE
  ),

  use_as_longitudinal_validation = c(
    FALSE,
    FALSE,
    TRUE,
    TRUE,
    TRUE,
    TRUE,
    TRUE,
    TRUE,
    TRUE,
    TRUE,
    FALSE,
    FALSE,
    FALSE,
    TRUE
  )
)

master74 <- master |>
  left_join(
    project_plan,
    by = "project"
  )

if (any(is.na(master74$analysis_module))) {
  warning(
    "Some projects in master are missing analysis-module assignment."
  )
}

# ------------------------------------------------------------
# Patient-level eligibility
# ------------------------------------------------------------
patient_elig <- master74 |>
  filter(!is.na(patient_uid)) |>
  group_by(
    project,
    patient_uid,
    patient_id,
    analysis_module,
    primary_use,
    data_modality,
    sequence_status
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
      time_raw[!is.na(time_raw)]
    ),

    phenotype_values =
      collapse_values(phenotype),

    infection_group_values =
      collapse_values(infection_group),

    outcome_values =
      collapse_values(outcome),

    organ_dysfunction_values =
      collapse_values(
        organ_dysfunction_status
      ),

    .groups = "drop"
  ) |>
  mutate(
    longitudinal_ge2 =
      n_timepoints >= 2,

    longitudinal_ge3 =
      n_timepoints >= 3,

    sequence_ready_now = case_when(
      sequence_status ==
        "FASTQ_DOWNLOADED_DADA2_PENDING" ~ FALSE,

      sequence_status ==
        "SHOTGUN_METADATA_READY_RAW_DEFERRED" ~ FALSE,

      TRUE ~ TRUE
    ),

    primary_longitudinal_eligible =
      analysis_module ==
        "CORE_SEPSIS_LONGITUDINAL" &
      longitudinal_ge2 &
      sequence_ready_now,

    validation_longitudinal_eligible =
      analysis_module %in% c(
        "ICU_INFECTION_LONGITUDINAL_EXTERNAL",
        "ICU_BACKGROUND_LONGITUDINAL",
        "NONSEPSIS_LONGITUDINAL_CONTROL",
        "INTERVENTION_LONGITUDINAL_SUPPORT",
        "ORGAN_DYSFUNCTION_LONGITUDINAL_SUPPORT"
      ) &
      longitudinal_ge2 &
      sequence_ready_now,

    deferred_longitudinal_eligible =
      longitudinal_ge2 &
      !sequence_ready_now
  )

write_excel_csv(
  patient_elig,
  file.path(
    OUT_ROOT,
    "V2_longitudinal_eligible_patients.csv"
  ),
  na = ""
)

# ------------------------------------------------------------
# Project-level eligibility
# ------------------------------------------------------------
project_elig <- master74 |>
  group_by(
    project,
    analysis_module,
    primary_use,
    cohort_role,
    analysis_tier,
    data_modality,
    sequence_status
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

    missing_patient =
      sum(is.na(patient_id)),

    missing_time =
      sum(is.na(time_raw)),

    n_unique_time_labels =
      n_distinct(
        time_raw[!is.na(time_raw)]
      ),

    time_labels =
      collapse_values(time_raw),

    phenotype_coverage =
      sum(!is.na(phenotype)),

    infection_group_coverage =
      sum(!is.na(infection_group)),

    outcome_coverage =
      sum(!is.na(outcome)),

    antibiotics_coverage =
      sum(!is.na(antibiotics)),

    severity_coverage =
      sum(!is.na(severity)),

    organ_dysfunction_coverage =
      sum(!is.na(organ_dysfunction_status)),

    .groups = "drop"
  ) |>
  left_join(
    patient_elig |>
      group_by(project) |>
      summarise(
        patients_ge2 =
          sum(longitudinal_ge2),

        patients_ge3 =
          sum(longitudinal_ge3),

        primary_eligible_patients =
          sum(primary_longitudinal_eligible),

        validation_eligible_patients =
          sum(validation_longitudinal_eligible),

        deferred_eligible_patients =
          sum(deferred_longitudinal_eligible),

        .groups = "drop"
      ),
    by = "project"
  ) |>
  mutate(
    project_longitudinal_status =
      case_when(
        primary_eligible_patients > 0 ~
          "PRIMARY_READY",

        validation_eligible_patients > 0 ~
          "VALIDATION_READY",

        deferred_eligible_patients > 0 ~
          "LONGITUDINAL_METADATA_READY_SEQUENCE_DEFERRED",

        patients_ge2 > 0 ~
          "LONGITUDINAL_METADATA_PRESENT_NOT_CURRENTLY_USED",

        TRUE ~
          "STATIC_ONLY"
      )
  )

write_excel_csv(
  project_elig,
  file.path(
    OUT_ROOT,
    "V2_longitudinal_eligibility_by_project.csv"
  ),
  na = ""
)

# ------------------------------------------------------------
# Clinical covariate audit
# ------------------------------------------------------------
clinical_fields <- c(
  "phenotype",
  "infection_group",
  "outcome",
  "antibiotics",
  "severity",
  "organ_dysfunction_status",
  "sex",
  "age"
)

clinical_audit <- map_dfr(
  clinical_fields,
  function(field) {

    master74 |>
      group_by(project) |>
      summarise(
        field = field,

        rows =
          n(),

        nonmissing =
          sum(
            !is.na(
              .data[[field]]
            )
          ),

        coverage_fraction =
          mean(
            !is.na(
              .data[[field]]
            )
          ),

        unique_values =
          collapse_values(
            .data[[field]]
          ),

        .groups = "drop"
      )
  }
)

write_excel_csv(
  clinical_audit,
  file.path(
    OUT_ROOT,
    "V2_clinical_covariate_coverage.csv"
  ),
  na = ""
)

# ------------------------------------------------------------
# Timepoint dictionary candidate table
# ------------------------------------------------------------
time_dictionary <- master74 |>
  group_by(
    project,
    analysis_module,
    time_raw,
    time_day,
    time_class,
    time_anchor
  ) |>
  summarise(
    rows = n(),

    patients =
      n_distinct(
        patient_uid,
        na.rm = TRUE
      ),

    .groups = "drop"
  ) |>
  arrange(
    project,
    time_day,
    time_raw
  )

write_excel_csv(
  time_dictionary,
  file.path(
    OUT_ROOT,
    "V2_timepoint_dictionary_candidates.csv"
  ),
  na = ""
)

# ------------------------------------------------------------
# Analysis architecture summary
# ------------------------------------------------------------
module_summary <- patient_elig |>
  group_by(
    analysis_module
  ) |>
  summarise(
    projects =
      n_distinct(project),

    patients =
      n(),

    patients_ge2 =
      sum(longitudinal_ge2),

    patients_ge3 =
      sum(longitudinal_ge3),

    sequence_ready_ge2 =
      sum(
        longitudinal_ge2 &
        sequence_ready_now
      ),

    .groups = "drop"
  )

write_excel_csv(
  module_summary,
  file.path(
    OUT_ROOT,
    "V2_analysis_module_summary.csv"
  ),
  na = ""
)

# ------------------------------------------------------------
# Scan data directories not represented in full master
# ------------------------------------------------------------
data_projects <- character()

if (dir.exists(DATA_ROOT)) {

  dirs <- list.dirs(
    DATA_ROOT,
    recursive = FALSE,
    full.names = FALSE
  )

  data_projects <- dirs[
    str_detect(
      dirs,
      regex(
        "^(PRJ|ERP|CRA)",
        ignore_case = TRUE
      )
    )
  ]
}

master_projects <- sort(
  unique(
    master74$project
  )
)

not_in_master <- tibble(
  project = setdiff(
    sort(data_projects),
    master_projects
  ),
  status = "PRESENT_IN_DATA_ROOT_BUT_NOT_IN_FULL_MASTER",
  required_action =
    "Review before declaring final cohort set; do not automatically include."
)

write_excel_csv(
  not_in_master,
  file.path(
    OUT_ROOT,
    "V2_projects_present_but_not_in_full_master.csv"
  ),
  na = ""
)

# ------------------------------------------------------------
# Overall summary
# ------------------------------------------------------------
overall <- tibble(
  Metric = c(
    "Full_master_projects",
    "Full_master_runs",
    "Full_master_patients",
    "Patients_GE2",
    "Patients_GE3",
    "Primary_16S_longitudinal_patients_ready",
    "Validation_16S_longitudinal_patients_ready",
    "Longitudinal_patients_sequence_deferred",
    "Projects_in_data_root_not_in_master"
  ),

  Value = c(
    n_distinct(master74$project),
    n_distinct(master74$run_uid, na.rm=TRUE),
    n_distinct(master74$patient_uid, na.rm=TRUE),
    sum(patient_elig$longitudinal_ge2),
    sum(patient_elig$longitudinal_ge3),
    sum(patient_elig$primary_longitudinal_eligible),
    sum(patient_elig$validation_longitudinal_eligible),
    sum(patient_elig$deferred_longitudinal_eligible),
    nrow(not_in_master)
  )
)

write_excel_csv(
  overall,
  file.path(
    OUT_ROOT,
    "V2_step74_overall_summary.csv"
  ),
  na = ""
)

# ------------------------------------------------------------
# Console
# ------------------------------------------------------------
cat(
  "\n============================================================\n"
)

cat(
  "SEPSIS V2 - STEP 74 COMPLETE\n"
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

cat("OVERALL:\n")
print(
  overall,
  n = Inf,
  width = Inf
)

cat("\nANALYSIS MODULE SUMMARY:\n")
print(
  module_summary,
  n = Inf,
  width = Inf
)

cat("\nPROJECTS PRESENT IN DATA ROOT BUT NOT MASTER:\n")
print(
  not_in_master,
  n = Inf,
  width = Inf
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
