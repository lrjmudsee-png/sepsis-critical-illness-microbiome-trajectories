# ============================================================
# Sepsis V2 - Step 73
# Merge validated legacy/V1 cohorts into the V2 master metadata
# R 4.4.0 / Windows
#
# Inputs
# 1) Latest V2 master from Step72
# 2) Frozen reviewed workbook:
#       03_Metadata_Analysis_Ready_v1.xlsx
#
# Legacy workbook sheets used:
#   Core_Case_Control
#   Longitudinal_Support
#   Organ_Dysfunction
#
# Only rows with:
#   Include_Flag == "Include"
#   Review_Status == "Confirmed"
# are imported.
#
# No patient IDs or timepoints are inferred from sample names here.
# ============================================================

options(stringsAsFactors = FALSE)

pkgs <- c(
  "readr",
  "readxl",
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
  message(
    "Installing missing packages: ",
    paste(missing, collapse = ", ")
  )

  install.packages(
    missing,
    repos = "https://cloud.r-project.org"
  )
}

suppressPackageStartupMessages({
  library(readr)
  library(readxl)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(purrr)
  library(tibble)
})


# ------------------------------------------------------------
# 1. Paths
# ------------------------------------------------------------
PROJECT_ROOT <- "E:/sepsis_project"

STEP72_ROOT <- paste0(
  PROJECT_ROOT,
  "/results/V2_12_master_metadata"
)

OUT_ROOT <- paste0(
  PROJECT_ROOT,
  "/results/V2_13_full_master_metadata"
)

dir.create(
  OUT_ROOT,
  recursive = TRUE,
  showWarnings = FALSE
)


# ------------------------------------------------------------
# 2. Helpers
# ------------------------------------------------------------
clean_chr <- function(x) {

  x <- trimws(
    as.character(x)
  )

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


clean_num <- function(x) {

  suppressWarnings(
    as.numeric(
      clean_chr(x)
    )
  )
}


find_latest <- function(
  folder,
  pattern
) {

  if (!dir.exists(folder)) {
    return(NA_character_)
  }

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


find_recursive_exact <- function(
  root,
  filename
) {

  x <- list.files(
    root,
    recursive = TRUE,
    full.names = TRUE
  )

  hit <- x[
    basename(x) == filename
  ]

  if (length(hit) == 0) {
    return(NA_character_)
  }

  # Prefer the newest copy if multiple copies exist.
  hit[
    which.max(
      file.info(hit)$mtime
    )
  ]
}


safe_csv <- function(path) {

  if (
    length(path) == 0 ||
      is.na(path) ||
      !file.exists(path)
  ) {
    return(NULL)
  }

  suppressMessages(
    read_csv(
      path,
      show_col_types = FALSE,
      progress = FALSE,
      name_repair = "unique"
    )
  )
}


col_chr <- function(
  df,
  name
) {

  if (!name %in% names(df)) {
    return(
      rep(
        NA_character_,
        nrow(df)
      )
    )
  }

  clean_chr(
    df[[name]]
  )
}


col_num <- function(
  df,
  name
) {

  if (!name %in% names(df)) {
    return(
      rep(
        NA_real_,
        nrow(df)
      )
    )
  }

  clean_num(
    df[[name]]
  )
}


# ------------------------------------------------------------
# 3. Load Step72 master
# ------------------------------------------------------------
step72_path <- find_latest(
  STEP72_ROOT,
  "^V2_master_metadata_[0-9]{8}_[0-9]{6}\\.csv$"
)

if (is.na(step72_path)) {
  stop(
    paste0(
      "Cannot find Step72 master metadata in: ",
      STEP72_ROOT
    )
  )
}


v2_new <- safe_csv(
  step72_path
)

if (is.null(v2_new)) {
  stop(
    "Failed to read Step72 master metadata."
  )
}


message(
  "Step72 master: ",
  step72_path
)

message(
  "Step72 rows: ",
  nrow(v2_new)
)


# ------------------------------------------------------------
# 4. Locate frozen V1 metadata workbook
# ------------------------------------------------------------
legacy_file <- find_recursive_exact(
  PROJECT_ROOT,
  "03_Metadata_Analysis_Ready_v1.xlsx"
)

if (is.na(legacy_file)) {

  stop(
    paste(
      "Cannot find 03_Metadata_Analysis_Ready_v1.xlsx",
      "anywhere under E:/sepsis_project."
    )
  )
}


message(
  "Legacy workbook: ",
  legacy_file
)


available_sheets <- excel_sheets(
  legacy_file
)


target_sheets <- c(
  "Core_Case_Control",
  "Longitudinal_Support",
  "Organ_Dysfunction"
)


missing_sheets <- setdiff(
  target_sheets,
  available_sheets
)


if (length(missing_sheets) > 0) {

  stop(
    paste0(
      "Legacy workbook missing sheet(s): ",
      paste(
        missing_sheets,
        collapse = ", "
      )
    )
  )
}


# ------------------------------------------------------------
# 5. Read reviewed legacy sheets
# ------------------------------------------------------------
legacy_parts <- list()


for (sheet in target_sheets) {

  message(
    "Reading legacy sheet: ",
    sheet
  )


  x <- suppressMessages(
    read_excel(
      legacy_file,
      sheet = sheet,
      .name_repair = "unique"
    )
  ) |>
    mutate(
      Legacy_Module = sheet
    )


  # Only frozen, confirmed analysis-ready rows.
  if (
    all(
      c(
        "Include_Flag",
        "Review_Status"
      ) %in% names(x)
    )
  ) {

    x <- x |>
      filter(
        Include_Flag == "Include",
        Review_Status == "Confirmed"
      )
  } else {

    stop(
      paste0(
        "Sheet ",
        sheet,
        " lacks Include_Flag or Review_Status."
      )
    )
  }


  legacy_parts[[sheet]] <- x
}


legacy_raw <- bind_rows(
  legacy_parts
)


# ------------------------------------------------------------
# 6. Legacy duplicate-run audit BEFORE deduplication
# ------------------------------------------------------------
legacy_run_modules <- legacy_raw |>
  filter(
    !is.na(Run_ID)
  ) |>
  group_by(
    Project_ID,
    Run_ID
  ) |>
  summarise(
    Legacy_Modules = paste(
      sort(
        unique(
          Legacy_Module
        )
      ),
      collapse = ";"
    ),

    N_Module_Rows = n(),

    Patient_ID_N = n_distinct(
      Patient_ID[
        !is.na(
          Patient_ID
        )
      ]
    ),

    Timepoint_N = n_distinct(
      Timepoint_Raw[
        !is.na(
          Timepoint_Raw
        )
      ]
    ),

    .groups = "drop"
  )


legacy_conflicts <- legacy_run_modules |>
  filter(
    Patient_ID_N > 1 |
      Timepoint_N > 1
  )


write_excel_csv(
  legacy_conflicts,
  file.path(
    OUT_ROOT,
    "V2_legacy_run_conflict_QC.csv"
  ),
  na = ""
)


if (nrow(legacy_conflicts) > 0) {

  stop(
    paste0(
      "Legacy workbook contains ",
      nrow(legacy_conflicts),
      " Run-level patient/time conflicts. ",
      "See V2_legacy_run_conflict_QC.csv."
    )
  )
}


# ------------------------------------------------------------
# 7. Collapse legacy duplicate module rows by Run
# ------------------------------------------------------------
# A baseline Run may legitimately appear in both Core_Case_Control
# and Longitudinal_Support. We keep it once and preserve all module
# memberships in legacy_module.
#
# We select the first non-missing value for each metadata field.
first_nonmissing <- function(x) {

  x <- clean_chr(x)

  hit <- x[
    !is.na(x)
  ]

  if (length(hit) == 0) {
    return(NA_character_)
  }

  hit[1]
}


first_nonmissing_num <- function(x) {

  x <- suppressWarnings(
    as.numeric(x)
  )

  hit <- x[
    !is.na(x)
  ]

  if (length(hit) == 0) {
    return(NA_real_)
  }

  hit[1]
}


legacy <- legacy_raw |>
  group_by(
    Project_ID,
    Run_ID
  ) |>
  summarise(
    BioSample_ID =
      first_nonmissing(
        BioSample_ID
      ),

    Patient_ID =
      first_nonmissing(
        Patient_ID
      ),

    Specimen_ID =
      first_nonmissing(
        Specimen_ID
      ),

    Sample_Alias =
      first_nonmissing(
        Sample_Alias
      ),

    Sample_Type =
      first_nonmissing(
        Sample_Type
      ),

    Phenotype =
      first_nonmissing(
        Phenotype
      ),

    Control_Type =
      first_nonmissing(
        Control_Type
      ),

    Primary_Sample_Flag =
      first_nonmissing(
        Primary_Sample_Flag
      ),

    Timepoint_Raw =
      first_nonmissing(
        Timepoint_Raw
      ),

    Timepoint_Day =
      first_nonmissing_num(
        Timepoint_Day
      ),

    Timepoint_Order =
      first_nonmissing_num(
        Timepoint_Order
      ),

    Time_Anchor =
      first_nonmissing(
        Time_Anchor
      ),

    Intervention =
      first_nonmissing(
        Intervention
      ),

    Intervention_Arm =
      first_nonmissing(
        Intervention_Arm
      ),

    Organ_Dysfunction_Type =
      first_nonmissing(
        Organ_Dysfunction_Type
      ),

    Organ_Dysfunction_Status =
      first_nonmissing(
        Organ_Dysfunction_Status
      ),

    Analysis_Role =
      first_nonmissing(
        Analysis_Role
      ),

    Dataset_Split =
      first_nonmissing(
        Dataset_Split
      ),

    Label_Source =
      first_nonmissing(
        Label_Source
      ),

    Label_Confidence =
      first_nonmissing(
        Label_Confidence
      ),

    Review_Status =
      first_nonmissing(
        Review_Status
      ),

    QC_Status =
      first_nonmissing(
        QC_Status
      ),

    Notes =
      first_nonmissing(
        Notes
      ),

    Legacy_Module = paste(
      sort(
        unique(
          Legacy_Module
        )
      ),
      collapse = ";"
    ),

    .groups = "drop"
  )


# ------------------------------------------------------------
# 8. Assign legacy project role/tier
# ------------------------------------------------------------
legacy <- legacy |>
  mutate(

    cohort_role = case_when(

      Project_ID == "PRJEB33360" ~
        "SEPSIS_ICU_CORE_AND_LONGITUDINAL",

      Project_ID == "PRJNA691455" ~
        "SEPSIS_CORE_AND_LONGITUDINAL",

      Project_ID == "PRJNA978257" ~
        "SEPSIS_HEALTHY_STATIC_DISCOVERY",

      Project_ID == "PRJNA1010969" ~
        "SEPSIS_TRAUMA_HEALTHY_STATIC_EXTERNAL",

      Project_ID == "PRJNA430161" ~
        "INTERVENTION_LONGITUDINAL_SUPPORT",

      Project_ID == "PRJNA1166732" ~
        "INTERVENTION_LONGITUDINAL_SUPPORT",

      Project_ID == "PRJNA797231" ~
        "ORGAN_DYSFUNCTION_STATIC_SUPPORT",

      Project_ID == "PRJNA912621" ~
        "ORGAN_DYSFUNCTION_LONGITUDINAL_SUPPORT",

      TRUE ~
        "LEGACY_REVIEWED_SUPPORT"
    ),


    analysis_tier = case_when(

      Project_ID %in% c(
        "PRJEB33360",
        "PRJNA691455",
        "PRJNA978257"
      ) ~
        "LEGACY_DEVELOPMENT",

      Project_ID == "PRJNA1010969" ~
        "LEGACY_LOCKED_EXTERNAL",

      Project_ID %in% c(
        "PRJNA430161",
        "PRJNA1166732"
      ) ~
        "LEGACY_INTERVENTION_SUPPORT",

      Project_ID %in% c(
        "PRJNA797231",
        "PRJNA912621"
      ) ~
        "LEGACY_ORGAN_SUPPORT",

      TRUE ~
        "LEGACY_SUPPORT"
    ),


    sequence_status = case_when(

      Project_ID == "PRJNA1166732" ~
        "FASTQ_DOWNLOADED_DADA2_PENDING",

      TRUE ~
        "EXISTING_ASV_OR_PROCESSED"
    ),


    data_modality =
      "16S"
  )


# ------------------------------------------------------------
# 9. Convert legacy metadata to expanded V2 schema
# ------------------------------------------------------------
legacy_std <- legacy |>
  transmute(

    project = clean_chr(
      Project_ID
    ),

    patient_id = clean_chr(
      Patient_ID
    ),

    sample_id = coalesce(
      clean_chr(
        Specimen_ID
      ),
      clean_chr(
        Sample_Alias
      ),
      clean_chr(
        BioSample_ID
      )
    ),

    biosample = clean_chr(
      BioSample_ID
    ),

    run_id = clean_chr(
      Run_ID
    ),

    experiment_id =
      NA_character_,

    time_raw = clean_chr(
      Timepoint_Raw
    ),

    time_day = as.numeric(
      Timepoint_Day
    ),

    time_class = clean_chr(
      Timepoint_Raw
    ),

    time_order = as.numeric(
      Timepoint_Order
    ),

    time_anchor = clean_chr(
      Time_Anchor
    ),

    collection_date =
      NA_character_,

    body_site = clean_chr(
      Sample_Type
    ),

    phenotype = clean_chr(
      Phenotype
    ),

    control_type = clean_chr(
      Control_Type
    ),

    infection_group = clean_chr(
      Phenotype
    ),

    outcome =
      NA_character_,

    antibiotics =
      NA_character_,

    severity =
      NA_character_,

    intervention = clean_chr(
      Intervention
    ),

    intervention_arm = clean_chr(
      Intervention_Arm
    ),

    organ_dysfunction_type = clean_chr(
      Organ_Dysfunction_Type
    ),

    organ_dysfunction_status = clean_chr(
      Organ_Dysfunction_Status
    ),

    sex =
      NA_character_,

    age =
      NA_character_,

    gestational_age =
      NA_character_,

    birthweight =
      NA_character_,

    primary_sample_flag = clean_chr(
      Primary_Sample_Flag
    ),

    analysis_role = clean_chr(
      Analysis_Role
    ),

    dataset_split = clean_chr(
      Dataset_Split
    ),

    cohort_role =
      cohort_role,

    analysis_tier =
      analysis_tier,

    sequence_status =
      sequence_status,

    data_modality =
      data_modality,

    metadata_source =
      basename(
        legacy_file
      ),

    label_source = clean_chr(
      Label_Source
    ),

    label_confidence = clean_chr(
      Label_Confidence
    ),

    review_status = clean_chr(
      Review_Status
    ),

    qc_status = clean_chr(
      QC_Status
    ),

    notes = clean_chr(
      Notes
    ),

    legacy_module = clean_chr(
      Legacy_Module
    )
  )


# ------------------------------------------------------------
# 10. Expand Step72 rows to same schema
# ------------------------------------------------------------
v2_new_std <- v2_new |>
  transmute(

    project =
      clean_chr(project),

    patient_id =
      clean_chr(patient_id),

    sample_id =
      clean_chr(sample_id),

    biosample =
      clean_chr(biosample),

    run_id =
      clean_chr(run_id),

    experiment_id =
      clean_chr(experiment_id),

    time_raw =
      clean_chr(time_raw),

    time_day =
      suppressWarnings(
        as.numeric(
          time_day
        )
      ),

    time_class =
      clean_chr(time_class),

    time_order =
      NA_real_,

    time_anchor =
      NA_character_,

    collection_date =
      clean_chr(collection_date),

    body_site =
      clean_chr(body_site),

    phenotype =
      clean_chr(infection_group),

    control_type =
      NA_character_,

    infection_group =
      clean_chr(infection_group),

    outcome =
      clean_chr(outcome),

    antibiotics =
      clean_chr(antibiotics),

    severity =
      clean_chr(severity),

    intervention =
      NA_character_,

    intervention_arm =
      NA_character_,

    organ_dysfunction_type =
      NA_character_,

    organ_dysfunction_status =
      NA_character_,

    sex =
      clean_chr(sex),

    age =
      clean_chr(age),

    gestational_age =
      clean_chr(gestational_age),

    birthweight =
      clean_chr(birthweight),

    primary_sample_flag =
      NA_character_,

    analysis_role =
      NA_character_,

    dataset_split =
      NA_character_,

    cohort_role =
      clean_chr(cohort_role),

    analysis_tier =
      clean_chr(analysis_tier),

    sequence_status = case_when(

      project == "PRJNA884103" ~
        "SHOTGUN_METADATA_READY_RAW_DEFERRED",

      TRUE ~
        "PUBLIC_METADATA_AND_SEQUENCE_AVAILABLE"
    ),

    data_modality = case_when(

      project == "PRJNA884103" ~
        "SHOTGUN",

      TRUE ~
        "16S"
    ),

    metadata_source =
      "V2_public_metadata_rescue_steps_63_to_72",

    label_source =
      NA_character_,

    label_confidence =
      NA_character_,

    review_status =
      "AUTO_VALIDATED",

    qc_status =
      "PASS_CORE_MAPPING",

    notes =
      NA_character_,

    legacy_module =
      NA_character_
  )


# ------------------------------------------------------------
# 11. Merge old + new
# ------------------------------------------------------------
full <- bind_rows(
  legacy_std,
  v2_new_std
)


# ------------------------------------------------------------
# 12. Global UID fields
# ------------------------------------------------------------
full <- full |>
  mutate(

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
# 13. Cross-source duplicate Run check
# ------------------------------------------------------------
dup_run <- full |>
  filter(
    !is.na(run_uid)
  ) |>
  count(
    run_uid,
    name = "n"
  ) |>
  filter(
    n > 1
  )


write_excel_csv(
  dup_run,
  file.path(
    OUT_ROOT,
    "V2_full_master_duplicate_run_QC.csv"
  ),
  na = ""
)


if (nrow(dup_run) > 0) {

  # Do not silently drop conflicting duplicated Runs.
  duplicate_rows <- full |>
    semi_join(
      dup_run,
      by = "run_uid"
    ) |>
    arrange(
      run_uid
    )

  write_excel_csv(
    duplicate_rows,
    file.path(
      OUT_ROOT,
      "V2_full_master_duplicate_run_rows_QC.csv"
    ),
    na = ""
  )

  stop(
    paste0(
      "Duplicate run_uid detected after legacy + V2 merge: ",
      nrow(dup_run),
      ". Inspect duplicate-run QC outputs."
    )
  )
}


# ------------------------------------------------------------
# 14. Longitudinal eligibility
# ------------------------------------------------------------
patient_time <- full |>
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


full <- full |>
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
# 15. Cohort QC
# ------------------------------------------------------------
cohort_qc <- full |>
  group_by(
    project,
    cohort_role,
    analysis_tier,
    data_modality,
    sequence_status
  ) |>
  summarise(

    rows =
      n(),

    patients =
      n_distinct(
        patient_uid,
        na.rm = TRUE
      ),

    samples =
      n_distinct(
        sample_uid,
        na.rm = TRUE
      ),

    runs =
      n_distinct(
        run_uid,
        na.rm = TRUE
      ),

    missing_patient =
      sum(
        is.na(
          patient_id
        )
      ),

    missing_sample =
      sum(
        is.na(
          sample_id
        )
      ),

    missing_run =
      sum(
        is.na(
          run_id
        )
      ),

    missing_time =
      sum(
        is.na(
          time_raw
        )
      ),

    patients_ge2 =
      n_distinct(
        patient_uid[
          longitudinal_ge2
        ],
        na.rm = TRUE
      ),

    patients_ge3 =
      n_distinct(
        patient_uid[
          longitudinal_ge3
        ],
        na.rm = TRUE
      ),

    phenotype_nonmissing =
      sum(
        !is.na(
          phenotype
        )
      ),

    infection_group_nonmissing =
      sum(
        !is.na(
          infection_group
        )
      ),

    outcome_nonmissing =
      sum(
        !is.na(
          outcome
        )
      ),

    antibiotics_nonmissing =
      sum(
        !is.na(
          antibiotics
        )
      ),

    organ_dysfunction_nonmissing =
      sum(
        !is.na(
          organ_dysfunction_status
        )
      ),

    .groups = "drop"
  )


# ------------------------------------------------------------
# 16. Patient summary
# ------------------------------------------------------------
patient_summary <- full |>
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

    n_rows =
      n(),

    n_samples =
      n_distinct(
        sample_uid,
        na.rm = TRUE
      ),

    n_runs =
      n_distinct(
        run_uid,
        na.rm = TRUE
      ),

    n_timepoints =
      n_distinct(
        time_raw[
          !is.na(
            time_raw
          )
        ]
      ),

    phenotype = paste(
      sort(
        unique(
          phenotype[
            !is.na(
              phenotype
            )
          ]
        )
      ),
      collapse = ";"
    ),

    organ_dysfunction_status = paste(
      sort(
        unique(
          organ_dysfunction_status[
            !is.na(
              organ_dysfunction_status
            )
          ]
        )
      ),
      collapse = ";"
    ),

    longitudinal_ge2 =
      n_timepoints >= 2,

    longitudinal_ge3 =
      n_timepoints >= 3,

    .groups = "drop"
  )


# ------------------------------------------------------------
# 17. Global summary
# ------------------------------------------------------------
global_summary <- tibble(

  Metric = c(
    "Rows",
    "Projects",
    "Patients",
    "Samples",
    "Runs",
    "Patients_GE2_timepoints",
    "Patients_GE3_timepoints",
    "16S_runs",
    "Shotgun_runs",
    "Missing_patient_rows",
    "Missing_time_rows",
    "Duplicate_run_uid"
  ),

  Value = c(

    nrow(full),

    n_distinct(
      full$project
    ),

    n_distinct(
      full$patient_uid,
      na.rm = TRUE
    ),

    n_distinct(
      full$sample_uid,
      na.rm = TRUE
    ),

    n_distinct(
      full$run_uid,
      na.rm = TRUE
    ),

    n_distinct(
      full$patient_uid[
        full$longitudinal_ge2
      ],
      na.rm = TRUE
    ),

    n_distinct(
      full$patient_uid[
        full$longitudinal_ge3
      ],
      na.rm = TRUE
    ),

    n_distinct(
      full$run_uid[
        full$data_modality == "16S"
      ],
      na.rm = TRUE
    ),

    n_distinct(
      full$run_uid[
        full$data_modality == "SHOTGUN"
      ],
      na.rm = TRUE
    ),

    sum(
      is.na(
        full$patient_id
      )
    ),

    sum(
      is.na(
        full$time_raw
      )
    ),

    nrow(
      dup_run
    )
  )
)


# ------------------------------------------------------------
# 18. Save
# ------------------------------------------------------------
stamp <- format(
  Sys.time(),
  "%Y%m%d_%H%M%S"
)


full_path <- file.path(
  OUT_ROOT,
  paste0(
    "V2_FULL_master_metadata_",
    stamp,
    ".csv"
  )
)


qc_path <- file.path(
  OUT_ROOT,
  paste0(
    "V2_FULL_master_QC_",
    stamp,
    ".csv"
  )
)


patient_path <- file.path(
  OUT_ROOT,
  paste0(
    "V2_FULL_master_patient_summary_",
    stamp,
    ".csv"
  )
)


summary_path <- file.path(
  OUT_ROOT,
  paste0(
    "V2_FULL_master_global_summary_",
    stamp,
    ".csv"
  )
)


legacy_path <- file.path(
  OUT_ROOT,
  paste0(
    "V2_legacy_metadata_standardized_",
    stamp,
    ".csv"
  )
)


write_excel_csv(
  full,
  full_path,
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
  global_summary,
  summary_path,
  na = ""
)

write_excel_csv(
  legacy_std,
  legacy_path,
  na = ""
)


# ------------------------------------------------------------
# 19. Console
# ------------------------------------------------------------
cat(
  "\n============================================================\n"
)

cat(
  "SEPSIS V2 - STEP 73 COMPLETE\n"
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
  "GLOBAL SUMMARY:\n"
)

print(
  global_summary,
  n = Inf,
  width = Inf
)


cat(
  "\nCOHORT QC:\n"
)

print(
  cohort_qc,
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
