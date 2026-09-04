# ============================================================
# Sepsis V2 - Step 66
# Build Sample -> Patient -> Time candidate mapping
# Designed for R 4.4.0 on Windows
#
# Ready projects:
#   PRJNA516701
#   PRJNA578267
#   PRJNA595346
#   PRJNA884103
#
# Unresolved projects retained without guessing patient IDs:
#   PRJEB37289
#   PRJEB67798
#   PRJEB82425
#
# Input:
#   E:/sepsis_project/data/<PROJECT>/00_metadata/
#
# Output:
#   E:/sepsis_project/results/V2_06_sample_patient_time_map/
#
# This script does not modify source metadata.
# ============================================================

options(stringsAsFactors = FALSE)

# ---------------------------
# 1. Packages
# ---------------------------
pkgs <- c(
  "readr",
  "dplyr",
  "stringr",
  "tibble",
  "purrr",
  "tidyr"
)

missing_pkgs <- pkgs[
  !vapply(
    pkgs,
    requireNamespace,
    logical(1),
    quietly = TRUE
  )
]

if (length(missing_pkgs) > 0) {
  message(
    "Installing missing packages: ",
    paste(missing_pkgs, collapse = ", ")
  )

  install.packages(
    missing_pkgs,
    repos = "https://cloud.r-project.org"
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


# ---------------------------
# 2. Paths
# ---------------------------
DATA_ROOT <- "E:/sepsis_project/data"

OUT_ROOT <- paste0(
  "E:/sepsis_project/results/",
  "V2_06_sample_patient_time_map"
)

dir.create(
  OUT_ROOT,
  recursive = TRUE,
  showWarnings = FALSE
)


READY_PROJECTS <- c(
  "PRJNA516701",
  "PRJNA578267",
  "PRJNA595346",
  "PRJNA884103"
)

UNRESOLVED_PROJECTS <- c(
  "PRJEB37289",
  "PRJEB67798",
  "PRJEB82425"
)


# ---------------------------
# 3. Helper functions
# ---------------------------
clean_chr <- function(x) {

  x <- as.character(x)

  x <- trimws(x)

  x[
    is.na(x) |
    x == "" |
    tolower(x) %in% c(
      "nan",
      "na",
      "n/a",
      "null",
      "none"
    )
  ] <- NA_character_

  x
}


safe_read_csv <- function(path) {

  if (!file.exists(path)) {
    return(NULL)
  }

  tryCatch(
    suppressMessages(
      read_csv(
        path,
        show_col_types = FALSE,
        progress = FALSE,
        name_repair = "unique"
      )
    ),
    error = function(e) {
      message(
        "FAILED CSV: ",
        path,
        " | ",
        conditionMessage(e)
      )

      NULL
    }
  )
}


safe_read_tsv <- function(path) {

  if (!file.exists(path)) {
    return(NULL)
  }

  tryCatch(
    suppressMessages(
      read_tsv(
        path,
        show_col_types = FALSE,
        progress = FALSE,
        name_repair = "unique"
      )
    ),
    error = function(e) {
      message(
        "FAILED TSV: ",
        path,
        " | ",
        conditionMessage(e)
      )

      NULL
    }
  )
}


first_existing <- function(
  df,
  candidates
) {

  if (is.null(df)) {
    return(NA_character_)
  }

  hit <- candidates[
    candidates %in% names(df)
  ]

  if (length(hit) == 0) {
    return(NA_character_)
  }

  hit[1]
}


get_col <- function(
  df,
  candidates,
  default = NA_character_
) {

  if (is.null(df)) {
    return(
      rep(
        default,
        0
      )
    )
  }

  col <- first_existing(
    df,
    candidates
  )

  if (is.na(col)) {
    return(
      rep(
        default,
        nrow(df)
      )
    )
  }

  clean_chr(
    df[[col]]
  )
}


find_file <- function(
  project,
  pattern
) {

  base <- file.path(
    DATA_ROOT,
    project,
    "00_metadata"
  )

  if (!dir.exists(base)) {
    return(NA_character_)
  }

  x <- list.files(
    base,
    recursive = TRUE,
    full.names = TRUE
  )

  hit <- x[
    str_detect(
      basename(x),
      regex(
        pattern,
        ignore_case = TRUE
      )
    )
  ]

  if (length(hit) == 0) {
    return(NA_character_)
  }

  hit[1]
}


read_ena_run <- function(project) {

  path <- find_file(
    project,
    "^01_ENA_read_run_metadata\\.tsv$"
  )

  if (is.na(path)) {

    path <- find_file(
      project,
      paste0(
        "^",
        project,
        "_ENA_run_manifest\\.tsv$"
      )
    )
  }

  if (is.na(path)) {
    return(NULL)
  }

  df <- safe_read_tsv(path)

  if (is.null(df)) {
    return(NULL)
  }

  df |>
    mutate(
      run_accession = get_col(
        df,
        c(
          "run_accession",
          "Run"
        )
      ),

      ena_sample_accession = get_col(
        df,
        c(
          "sample_accession"
        )
      ),

      ena_secondary_sample_accession = get_col(
        df,
        c(
          "secondary_sample_accession"
        )
      ),

      ena_sample_alias = get_col(
        df,
        c(
          "sample_alias"
        )
      ),

      ena_sample_title = get_col(
        df,
        c(
          "sample_title"
        )
      ),

      ena_collection_date = get_col(
        df,
        c(
          "collection_date"
        )
      ),

      ena_body_site = get_col(
        df,
        c(
          "host_body_site",
          "isolation_source"
        )
      )
    ) |>
    select(
      run_accession,
      ena_sample_accession,
      ena_secondary_sample_accession,
      ena_sample_alias,
      ena_sample_title,
      ena_collection_date,
      ena_body_site,
      everything()
    )
}


read_biosample <- function(project) {

  path <- find_file(
    project,
    "^06_NCBI_BioSample_attributes\\.csv$"
  )

  if (is.na(path)) {
    return(NULL)
  }

  df <- safe_read_csv(path)

  if (is.null(df)) {
    return(NULL)
  }

  df
}


# Join BioSample metadata to ENA run metadata.
# ENA may store the NCBI BioSample accession in either
# sample_accession or secondary_sample_accession.
join_biosample_to_ena <- function(
  biosample,
  ena
) {

  if (
    is.null(biosample) ||
    is.null(ena)
  ) {
    return(NULL)
  }

  biosample <- biosample |>
    mutate(
      BioSample_accession = clean_chr(
        BioSample_accession
      )
    )

  ena1 <- ena |>
    select(
      run_accession,
      ena_sample_accession,
      ena_secondary_sample_accession,
      ena_sample_alias,
      ena_sample_title,
      ena_collection_date,
      ena_body_site
    )


  by_primary <- biosample |>
    left_join(
      ena1,
      by = c(
        "BioSample_accession" =
          "ena_sample_accession"
      )
    )


  unmatched_primary <- is.na(
    by_primary$run_accession
  )


  if (any(unmatched_primary)) {

    second <- biosample[
      unmatched_primary,
      ,
      drop = FALSE
    ] |>
      left_join(
        ena1,
        by = c(
          "BioSample_accession" =
            "ena_secondary_sample_accession"
        )
      )


    by_primary[
      unmatched_primary,
      names(second)
    ] <- second
  }


  by_primary
}


standard_map_columns <- function(
  project,
  df,
  patient,
  sample,
  run,
  time_raw,
  time_day = NA_character_,
  collection_date = NA_character_,
  body_site = NA_character_,
  source_file = NA_character_,
  mapping_method = NA_character_,
  extra = list()
) {

  n <- nrow(df)

  out <- tibble(
    Project = rep(
      project,
      n
    ),

    Patient_ID = clean_chr(
      patient
    ),

    Sample_ID = clean_chr(
      sample
    ),

    Run_ID = clean_chr(
      run
    ),

    Time_Raw = clean_chr(
      time_raw
    ),

    Time_Day = suppressWarnings(
      as.numeric(
        clean_chr(
          time_day
        )
      )
    ),

    Collection_Date = clean_chr(
      collection_date
    ),

    Body_Site = clean_chr(
      body_site
    ),

    Source_File = rep(
      source_file,
      n
    ),

    Mapping_Method = rep(
      mapping_method,
      n
    )
  )


  if (length(extra) > 0) {

    for (nm in names(extra)) {

      value <- extra[[nm]]

      if (length(value) == 1) {
        value <- rep(
          value,
          n
        )
      }

      out[[nm]] <- value
    }
  }

  out
}


# ---------------------------
# 4. PRJNA516701
# ---------------------------
build_PRJNA516701 <- function() {

  project <- "PRJNA516701"

  bio <- read_biosample(
    project
  )

  ena <- read_ena_run(
    project
  )

  if (
    is.null(bio) ||
    is.null(ena)
  ) {
    return(NULL)
  }


  joined <- join_biosample_to_ena(
    bio,
    ena
  )


  patient <- get_col(
    joined,
    c(
      "SubjectID",
      "Patient",
      "Subject"
    )
  )

  sample <- get_col(
    joined,
    c(
      "BioSample_accession"
    )
  )

  run <- get_col(
    joined,
    c(
      "run_accession"
    )
  )

  study_day <- get_col(
    joined,
    c(
      "Study_Day",
      "Day"
    )
  )

  study_period <- get_col(
    joined,
    c(
      "Study_Period",
      "timepoint"
    )
  )

  time_raw <- ifelse(
    !is.na(study_period),
    study_period,
    study_day
  )

  collection <- get_col(
    joined,
    c(
      "collection_date",
      "ena_collection_date"
    )
  )


  standard_map_columns(
    project = project,
    df = joined,
    patient = patient,
    sample = sample,
    run = run,
    time_raw = time_raw,
    time_day = study_day,
    collection_date = collection,
    body_site = get_col(
      joined,
      c(
        "ena_body_site"
      )
    ),
    source_file =
      "06_NCBI_BioSample_attributes.csv + 01_ENA_read_run_metadata.tsv",
    mapping_method =
      "BioSample accession joined to ENA sample accession",
    extra = list(
      Study_Period = study_period
    )
  )
}


# ---------------------------
# 5. PRJNA578267
# ---------------------------
build_PRJNA578267 <- function() {

  project <- "PRJNA578267"

  bio <- read_biosample(
    project
  )

  ena <- read_ena_run(
    project
  )

  if (
    is.null(bio) ||
    is.null(ena)
  ) {
    return(NULL)
  }


  joined <- join_biosample_to_ena(
    bio,
    ena
  )


  patient <- get_col(
    joined,
    c(
      "Patient",
      "SubjectID",
      "Subject"
    )
  )

  sample <- get_col(
    joined,
    c(
      "BioSample_accession"
    )
  )

  timepoint <- get_col(
    joined,
    c(
      "timepoint",
      "Timepoint"
    )
  )


  standard_map_columns(
    project = project,
    df = joined,
    patient = patient,
    sample = sample,
    run = get_col(
      joined,
      c(
        "run_accession"
      )
    ),
    time_raw = timepoint,
    collection_date = get_col(
      joined,
      c(
        "collection_date",
        "ena_collection_date"
      )
    ),
    body_site = get_col(
      joined,
      c(
        "isolation_source",
        "ena_body_site"
      )
    ),
    source_file =
      "06_NCBI_BioSample_attributes.csv + 01_ENA_read_run_metadata.tsv",
    mapping_method =
      "BioSample accession joined to ENA sample accession"
  )
}


# ---------------------------
# 6. PRJNA595346
# ---------------------------
build_PRJNA595346 <- function() {

  project <- "PRJNA595346"

  bio <- read_biosample(
    project
  )

  ena <- read_ena_run(
    project
  )

  if (
    is.null(bio) ||
    is.null(ena)
  ) {
    return(NULL)
  }


  joined <- join_biosample_to_ena(
    bio,
    ena
  )


  patient <- get_col(
    joined,
    c(
      "SubjectID",
      "Patient",
      "Subject"
    )
  )

  day <- get_col(
    joined,
    c(
      "Day",
      "Study_Day"
    )
  )

  study_date <- get_col(
    joined,
    c(
      "StudyDate",
      "collection_date"
    )
  )

  sample_type <- get_col(
    joined,
    c(
      "sample_type",
      "tissue",
      "ena_body_site"
    )
  )


  standard_map_columns(
    project = project,
    df = joined,
    patient = patient,
    sample = get_col(
      joined,
      c(
        "BioSample_accession",
        "SampleID_Merge"
      )
    ),
    run = get_col(
      joined,
      c(
        "run_accession"
      )
    ),
    time_raw = day,
    time_day = day,
    collection_date = study_date,
    body_site = sample_type,
    source_file =
      "06_NCBI_BioSample_attributes.csv + 01_ENA_read_run_metadata.tsv",
    mapping_method =
      "BioSample accession joined to ENA sample accession",
    extra = list(
      SampleID_Merge = get_col(
        joined,
        c(
          "SampleID_Merge"
        )
      ),
      Age = get_col(
        joined,
        c(
          "age"
        )
      ),
      Sex = get_col(
        joined,
        c(
          "sex"
        )
      )
    )
  )
}


# ---------------------------
# 7. PRJNA884103
# ---------------------------
build_PRJNA884103 <- function() {

  project <- "PRJNA884103"

  master_path <- find_file(
    project,
    "^MasterMetadata_FINAL_GitHub\\.csv$"
  )

  ena <- read_ena_run(
    project
  )

  if (
    is.na(master_path) ||
    is.null(ena)
  ) {
    return(NULL)
  }


  master <- safe_read_csv(
    master_path
  )

  if (is.null(master)) {
    return(NULL)
  }


  master <- master |>
    mutate(
      Sequence_name_join = get_col(
        master,
        c(
          "Sequence.name",
          "Sequence_name",
          "Sample"
        )
      )
    )


  ena_small <- ena |>
    transmute(
      run_accession,
      ena_sample_accession,
      ena_secondary_sample_accession,
      ena_sample_alias,
      ena_sample_title,
      ena_collection_date,
      ena_body_site
    )


  joined <- master |>
    left_join(
      ena_small,
      by = c(
        "Sequence_name_join" =
          "ena_sample_alias"
      )
    )


  # Secondary attempt using Sample if Sequence.name did not join.
  if ("Sample" %in% names(master)) {

    miss <- is.na(
      joined$run_accession
    )

    if (any(miss)) {

      second <- master[
        miss,
        ,
        drop = FALSE
      ] |>
        mutate(
          Sample_join = clean_chr(
            Sample
          )
        ) |>
        left_join(
          ena_small,
          by = c(
            "Sample_join" =
              "ena_sample_alias"
          )
        )


      common <- intersect(
        names(joined),
        names(second)
      )

      joined[
        miss,
        common
      ] <- second[
        ,
        common,
        drop = FALSE
      ]
    }
  }


  patient <- get_col(
    joined,
    c(
      "Subject",
      "Subject_id"
    )
  )

  dol <- get_col(
    joined,
    c(
      "DOL"
    )
  )

  sample <- get_col(
    joined,
    c(
      "Sequence.name",
      "Sample",
      "Sample_id"
    )
  )


  standard_map_columns(
    project = project,
    df = joined,
    patient = patient,
    sample = sample,
    run = get_col(
      joined,
      c(
        "run_accession"
      )
    ),
    time_raw = dol,
    time_day = dol,
    collection_date = get_col(
      joined,
      c(
        "ena_collection_date"
      )
    ),
    body_site = get_col(
      joined,
      c(
        "sample_type",
        "Specimen..",
        "ena_body_site"
      )
    ),
    source_file =
      "MasterMetadata_FINAL_GitHub.csv + 01_ENA_read_run_metadata.tsv",
    mapping_method =
      "GitHub Sequence.name joined to ENA sample_alias",
    extra = list(

      Infection_Group = get_col(
        joined,
        c(
          "BSI",
          "bacteremia",
          "Bacteremia"
        )
      ),

      Bacteremia_DOL = get_col(
        joined,
        c(
          "Bacteremia_DOL",
          "BSI_DOL"
        )
      ),

      Outcome = get_col(
        joined,
        c(
          "outcome",
          "Outcome"
        )
      ),

      Abx_days = get_col(
        joined,
        c(
          "Abx_days"
        )
      ),

      ABx_days_beforebacteremia = get_col(
        joined,
        c(
          "ABx_days_beforebacteremia"
        )
      ),

      ABx_score_sample = get_col(
        joined,
        c(
          "ABx_score_sample"
        )
      ),

      Sex = get_col(
        joined,
        c(
          "sex"
        )
      ),

      Gestational_age = get_col(
        joined,
        c(
          "gest_age"
        )
      ),

      Birthweight = get_col(
        joined,
        c(
          "birthweight"
        )
      )
    )
  )
}


# ---------------------------
# 8. Build maps
# ---------------------------
builders <- list(
  PRJNA516701 = build_PRJNA516701,
  PRJNA578267 = build_PRJNA578267,
  PRJNA595346 = build_PRJNA595346,
  PRJNA884103 = build_PRJNA884103
)


ready_maps <- list()


for (project in names(builders)) {

  message(
    "Building map: ",
    project
  )

  x <- tryCatch(
    builders[[project]](),
    error = function(e) {

      message(
        "FAILED ",
        project,
        ": ",
        conditionMessage(e)
      )

      NULL
    }
  )


  if (!is.null(x)) {

    ready_maps[[project]] <- x

    write_excel_csv(
      x,
      file.path(
        OUT_ROOT,
        paste0(
          project,
          "_sample_patient_time_map.csv"
        )
      ),
      na = ""
    )
  }
}


combined_ready <- if (
  length(ready_maps) > 0
) {

  bind_rows(
    ready_maps
  )

} else {

  tibble()
}


# ---------------------------
# 9. Unresolved projects
# Keep sample/run/time information,
# but NEVER infer patient IDs from sample names.
# ---------------------------
unresolved_rows <- list()


for (project in UNRESOLVED_PROJECTS) {

  ena <- read_ena_run(
    project
  )

  if (is.null(ena)) {
    next
  }


  unresolved_rows[[project]] <- ena |>
    transmute(
      Project = project,

      Patient_ID = NA_character_,

      Sample_ID = dplyr::coalesce(
        clean_chr(
          ena_sample_alias
        ),
        clean_chr(
          ena_secondary_sample_accession
        ),
        clean_chr(
          ena_sample_accession
        )
      ),

      Run_ID = clean_chr(
        run_accession
      ),

      Time_Raw = clean_chr(
        ena_collection_date
      ),

      Time_Day = NA_real_,

      Collection_Date = clean_chr(
        ena_collection_date
      ),

      Body_Site = clean_chr(
        ena_body_site
      ),

      Patient_Mapping_Status =
        "UNRESOLVED_DO_NOT_INFER",

      Needed_Field =
        "patient_id",

      Recommended_Action =
        "Parse existing supplement/repository before any manual collection"
    )
}


unresolved <- if (
  length(unresolved_rows) > 0
) {

  bind_rows(
    unresolved_rows
  )

} else {

  tibble()
}


# ---------------------------
# 10. QC
# ---------------------------
qc_ready <- combined_ready |>
  group_by(
    Project
  ) |>
  summarise(

    Rows = n(),

    Unique_Patients = n_distinct(
      Patient_ID[
        !is.na(
          Patient_ID
        )
      ]
    ),

    Unique_Samples = n_distinct(
      Sample_ID[
        !is.na(
          Sample_ID
        )
      ]
    ),

    Unique_Runs = n_distinct(
      Run_ID[
        !is.na(
          Run_ID
        )
      ]
    ),

    Missing_Patient = sum(
      is.na(
        Patient_ID
      )
    ),

    Missing_Sample = sum(
      is.na(
        Sample_ID
      )
    ),

    Missing_Run = sum(
      is.na(
        Run_ID
      )
    ),

    Missing_Time = sum(
      is.na(
        Time_Raw
      )
    ),

    Patients_GE2_Samples = {

      tmp <- tibble(
        Patient_ID = Patient_ID,
        Sample_ID = Sample_ID
      ) |>
        filter(
          !is.na(
            Patient_ID
          )
        ) |>
        distinct(
          Patient_ID,
          Sample_ID
        ) |>
        count(
          Patient_ID,
          name = "n_samples"
        )

      sum(
        tmp$n_samples >= 2
      )
    },

    .groups = "drop"
  )


# Add duplicate diagnostics.
dup_patient_sample <- combined_ready |>
  filter(
    !is.na(
      Patient_ID
    ),
    !is.na(
      Sample_ID
    )
  ) |>
  count(
    Project,
    Patient_ID,
    Sample_ID,
    name = "Rows"
  ) |>
  filter(
    Rows > 1
  )


run_unmatched <- combined_ready |>
  filter(
    is.na(
      Run_ID
    )
  )


# ---------------------------
# 11. Save combined outputs
# ---------------------------
stamp <- format(
  Sys.time(),
  "%Y%m%d_%H%M%S"
)


p_ready <- file.path(
  OUT_ROOT,
  paste0(
    "V2_sample_patient_time_map_READY_",
    stamp,
    ".csv"
  )
)


p_unresolved <- file.path(
  OUT_ROOT,
  paste0(
    "V2_patient_mapping_UNRESOLVED_",
    stamp,
    ".csv"
  )
)


p_qc <- file.path(
  OUT_ROOT,
  paste0(
    "V2_sample_patient_time_QC_",
    stamp,
    ".csv"
  )
)


p_dup <- file.path(
  OUT_ROOT,
  paste0(
    "V2_duplicate_patient_sample_QC_",
    stamp,
    ".csv"
  )
)


p_unmatched <- file.path(
  OUT_ROOT,
  paste0(
    "V2_unmatched_run_QC_",
    stamp,
    ".csv"
  )
)


write_excel_csv(
  combined_ready,
  p_ready,
  na = ""
)


write_excel_csv(
  unresolved,
  p_unresolved,
  na = ""
)


write_excel_csv(
  qc_ready,
  p_qc,
  na = ""
)


write_excel_csv(
  dup_patient_sample,
  p_dup,
  na = ""
)


write_excel_csv(
  run_unmatched,
  p_unmatched,
  na = ""
)


# ---------------------------
# 12. Console summary
# ---------------------------
cat(
  "\n============================================================\n"
)

cat(
  "SEPSIS V2 - STEP 66 COMPLETE\n"
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
  "READY PROJECT QC:\n"
)

print(
  qc_ready,
  n = Inf,
  width = Inf
)


cat(
  "\nUNRESOLVED PATIENT-ID PROJECTS:\n"
)

if (
  nrow(
    unresolved
  ) > 0
) {

  print(
    unresolved |>
      count(
        Project,
        name = "Rows"
      ),
    n = Inf
  )
}


cat(
  "\nOutput files:\n"
)

cat(
  p_ready,
  "\n"
)

cat(
  p_unresolved,
  "\n"
)

cat(
  p_qc,
  "\n"
)

cat(
  p_dup,
  "\n"
)

cat(
  p_unmatched,
  "\n"
)

cat(
  "============================================================\n"
)
