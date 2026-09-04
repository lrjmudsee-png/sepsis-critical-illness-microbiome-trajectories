# ============================================================
# Sepsis V2 - Step 80
# Finalize CRA002354 + PRJNA851469 master-ready metadata
# OFFLINE ONLY
#
# What this script does:
# 1. CRA002354:
#    - correctly promotes the real header row in S1/S2
#    - exact Patient ID join
#    - creates 131-sample longitudinal clinical map
#    - maps Sample ID -> GSA Run automatically from local CRA002354.xlsx
#
# 2. PRJNA851469:
#    - reads Supplementary Table 2 sample IDs
#    - parses EXPLICIT PatientXX-Day-Y metadata
#    - separates healthy controls from ICU longitudinal samples
#    - automatically searches existing local ENA metadata for Run mapping
#
# 3. Produces final QC and cohort decision files.
#
# No internet access.
# R 4.4.0 / Windows
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
  !vapply(pkgs, requireNamespace, logical(1), quietly=TRUE)
]

if (length(missing) > 0) {
  install.packages(
    missing,
    repos="https://cloud.r-project.org"
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

# ============================================================
# PATHS
# ============================================================

ROOT <- "E:/sepsis_project"

STEP79 <- file.path(
  ROOT,
  "results",
  "V2_19_offline_metadata_parse"
)

OUT <- file.path(
  ROOT,
  "results",
  "V2_20_finalize_CRA002354_PRJNA851469"
)

dir.create(
  OUT,
  recursive=TRUE,
  showWarnings=FALSE
)

LOG <- file.path(
  OUT,
  "_STEP80_runtime_checkpoints.txt"
)

cat(
  paste0(
    "STEP80 STARTED: ",
    Sys.time(),
    "\n"
  ),
  file=LOG,
  append=TRUE
)

# ============================================================
# HELPERS
# ============================================================

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

norm_id <- function(x) {

  x <- clean_chr(x)

  x <- toupper(x)

  x <- str_replace_all(
    x,
    "[^A-Z0-9]",
    ""
  )

  x
}

safe_read_csv <- function(path) {

  tryCatch(
    suppressMessages(
      read_csv(
        path,
        show_col_types=FALSE,
        progress=FALSE
      )
    ),
    error=function(e) NULL
  )
}

safe_read_tsv <- function(path) {

  tryCatch(
    suppressMessages(
      read_tsv(
        path,
        show_col_types=FALSE,
        progress=FALSE
      )
    ),
    error=function(e) NULL
  )
}

first_nonmissing <- function(x) {

  x <- clean_chr(x)

  x <- x[
    !is.na(x)
  ]

  if (length(x)==0) {
    return(NA_character_)
  }

  x[1]
}

parse_day <- function(x) {

  as.integer(
    str_extract(
      clean_chr(x),
      "(?i)(?<=D)\\d+|(?<=DAY-)\\d+"
    )
  )
}

find_col_name <- function(
  nms,
  patterns
) {

  hit <- which(
    str_detect(
      nms,
      regex(
        paste(
          patterns,
          collapse="|"
        ),
        ignore_case=TRUE
      )
    )
  )

  if (length(hit)==0) {
    return(NA_character_)
  }

  nms[hit[1]]
}

count_regex_values <- function(
  x,
  pattern
) {

  x <- clean_chr(x)

  sum(
    str_detect(
      x,
      regex(
        pattern,
        ignore_case=TRUE
      )
    ),
    na.rm=TRUE
  )
}

# ============================================================
# PART A - CRA002354
# ============================================================

CRA_S1_PATH <- file.path(
  STEP79,
  "CRA002354_S1_sheet_01_sheet1.csv"
)

CRA_S2_PATH <- file.path(
  STEP79,
  "CRA002354_S2_sheet_01_sheet1.csv"
)

if (
  !file.exists(CRA_S1_PATH) ||
  !file.exists(CRA_S2_PATH)
) {
  stop(
    "Step79 CRA S1/S2 exported CSV files not found."
  )
}

cra_s1_raw <- safe_read_csv(
  CRA_S1_PATH
)

cra_s2_raw <- safe_read_csv(
  CRA_S2_PATH
)

if (
  is.null(cra_s1_raw) ||
  is.null(cra_s2_raw)
) {
  stop(
    "Failed to read Step79 CRA S1/S2 CSV files."
  )
}

# ------------------------------------------------------------
# CRA workbook structure:
#
# Original Excel row 1 = table title
# Original Excel row 2 = note
# Original Excel row 3 = TRUE COLUMN HEADERS
#
# Step79 read row 1 as the CSV column header.
# Therefore in the Step79 CSV:
#   data row 1 = note
#   data row 2 = true headers
#   data rows 3+ = actual data
# ------------------------------------------------------------

promote_real_header <- function(
  x,
  header_row=2
) {

  new_names <- clean_chr(
    unlist(
      x[
        header_row,
        ,
        drop=TRUE
      ],
      use.names=FALSE
    )
  )

  bad <- is.na(new_names)

  if (any(bad)) {
    new_names[bad] <- paste0(
      "unnamed_",
      which(bad)
    )
  }

  new_names <- make.unique(
    new_names,
    sep="_"
  )

  dat <- x[
    seq.int(
      header_row + 1,
      nrow(x)
    ),
    ,
    drop=FALSE
  ]

  names(dat) <- new_names

  dat |>
    mutate(
      across(
        everything(),
        clean_chr
      )
    )
}

cra_s1 <- promote_real_header(
  cra_s1_raw,
  header_row=2
)

cra_s2 <- promote_real_header(
  cra_s2_raw,
  header_row=2
)

# Sanity checks from the published supplementary tables.
if (nrow(cra_s1) != 64) {

  stop(
    paste0(
      "CRA S1 expected 64 patient rows, found ",
      nrow(cra_s1)
    )
  )
}

if (nrow(cra_s2) != 131) {

  stop(
    paste0(
      "CRA S2 expected 131 sample rows, found ",
      nrow(cra_s2)
    )
  )
}

required_s1 <- c(
  "Patient ID",
  "Collected day",
  "Sample ID"
)

required_s2 <- c(
  "Patient ID",
  "Collected day",
  "Sample ID"
)

if (
  !all(
    required_s1 %in% names(cra_s1)
  )
) {
  stop(
    "CRA S1 true header parsing failed."
  )
}

if (
  !all(
    required_s2 %in% names(cra_s2)
  )
) {
  stop(
    "CRA S2 true header parsing failed."
  )
}

# ------------------------------------------------------------
# Patient-level clinical table
# ------------------------------------------------------------

cra_patient <- cra_s1 |>
  rename(
    patient_id=`Patient ID`,
    first_collected_day=`Collected day`,
    first_sample_id=`Sample ID`
  ) |>
  mutate(
    project="CRA002354",
    patient_id=clean_chr(patient_id),
    first_sample_id=clean_chr(first_sample_id),
    .before=1
  )

# ------------------------------------------------------------
# Sample-level table
# ------------------------------------------------------------

cra_sample <- cra_s2 |>
  rename(
    sample_id=`Sample ID`,
    time_raw=`Collected day`,
    patient_id=`Patient ID`
  ) |>
  mutate(
    project="CRA002354",
    patient_id=clean_chr(patient_id),
    sample_id=clean_chr(sample_id),
    time_raw=clean_chr(time_raw),
    time_day=parse_day(time_raw),
    time_order=time_day,
    .before=1
  )

# ------------------------------------------------------------
# Exact patient join
# ------------------------------------------------------------

cra_master <- cra_sample |>
  left_join(
    cra_patient,
    by=c(
      "project",
      "patient_id"
    ),
    suffix=c(
      "_sample",
      "_patient"
    )
  )

# ------------------------------------------------------------
# Verify S1 first samples exist EXACTLY in S2.
# ------------------------------------------------------------

cra_first_sample_check <- cra_patient |>
  select(
    patient_id,
    first_sample_id,
    first_collected_day
  ) |>
  left_join(
    cra_sample |>
      select(
        patient_id,
        sample_id,
        time_raw
      ),
    by=c(
      "patient_id",
      "first_sample_id"="sample_id"
    )
  ) |>
  mutate(
    sample_exact_match=
      !is.na(time_raw),

    day_exact_match=
      clean_chr(
        first_collected_day
      ) ==
      clean_chr(
        time_raw
      )
  )

# ============================================================
# CRA GSA RUN MAPPING
# ============================================================

# Search local project tree for non-supplement workbooks/tables.
cra_search_roots <- c(
  file.path(
    ROOT,
    "data",
    "CRA002354"
  ),
  file.path(
    ROOT,
    "results"
  )
)

cra_files <- character()

for (r in cra_search_roots) {

  if (!dir.exists(r)) next

  f <- list.files(
    r,
    recursive=TRUE,
    full.names=TRUE,
    pattern="\\.(csv|tsv|xlsx)$",
    ignore.case=TRUE
  )

  cra_files <- c(
    cra_files,
    f
  )
}

cra_files <- unique(
  cra_files
)

cra_files <- cra_files[
  !str_detect(
    cra_files,
    regex(
      "V2_19_offline_metadata_parse|V2_20_finalize|step78A",
      ignore_case=TRUE
    )
  )
]

cra_ids_norm <- unique(
  norm_id(
    cra_sample$sample_id
  )
)

cra_candidate_tables <- list()
cra_candidate_scores <- list()

inspect_table_for_cra <- function(
  dat,
  source_file,
  source_sheet=NA_character_
) {

  if (
    is.null(dat) ||
    nrow(dat)==0 ||
    ncol(dat)==0
  ) {
    return(NULL)
  }

  dat <- as_tibble(dat)

  overlaps <- map_int(
    dat,
    function(v) {

      vv <- unique(
        norm_id(v)
      )

      length(
        intersect(
          vv[
            !is.na(vv)
          ],
          cra_ids_norm
        )
      )
    }
  )

  best_overlap <- max(
    overlaps,
    na.rm=TRUE
  )

  if (
    !is.finite(best_overlap) ||
    best_overlap == 0
  ) {
    return(NULL)
  }

  best_col <- names(dat)[
    which.max(overlaps)
  ]

  run_scores <- map_int(
    dat,
    count_regex_values,
    pattern="^(CRR|SRR|ERR|DRR)[0-9]+$"
  )

  run_col <- if (
    max(
      run_scores,
      na.rm=TRUE
    ) > 0
  ) {
    names(dat)[
      which.max(
        run_scores
      )
    ]
  } else {
    NA_character_
  }

  score <- tibble(
    source_file=source_file,
    source_sheet=source_sheet,
    rows=nrow(dat),
    columns=ncol(dat),
    sample_match_column=best_col,
    exact_sample_overlap=best_overlap,
    run_column=run_col,
    run_accession_count=
      ifelse(
        is.na(run_col),
        0L,
        max(
          run_scores,
          na.rm=TRUE
        )
      )
  )

  list(
    data=dat,
    score=score
  )
}

for (p in cra_files) {

  ext <- tolower(
    tools::file_ext(p)
  )

  # Avoid giant unrelated tables.
  if (
    file.exists(p) &&
    file.info(p)$size >
      50 * 1024^2
  ) next

  if (ext=="csv") {

    dat <- safe_read_csv(p)

    z <- inspect_table_for_cra(
      dat,
      p
    )

    if (!is.null(z)) {

      cra_candidate_tables[[
        length(
          cra_candidate_tables
        ) + 1
      ]] <- z

      cra_candidate_scores[[
        length(
          cra_candidate_scores
        ) + 1
      ]] <- z$score
    }
  }

  if (ext=="tsv") {

    dat <- safe_read_tsv(p)

    z <- inspect_table_for_cra(
      dat,
      p
    )

    if (!is.null(z)) {

      cra_candidate_tables[[
        length(
          cra_candidate_tables
        ) + 1
      ]] <- z

      cra_candidate_scores[[
        length(
          cra_candidate_scores
        ) + 1
      ]] <- z$score
    }
  }

  if (ext=="xlsx") {

    sheets <- tryCatch(
      excel_sheets(p),
      error=function(e) character()
    )

    for (sh in sheets) {

      dat <- tryCatch(
        read_excel(
          p,
          sheet=sh,
          .name_repair="unique"
        ),
        error=function(e) NULL
      )

      z <- inspect_table_for_cra(
        dat,
        p,
        sh
      )

      if (!is.null(z)) {

        cra_candidate_tables[[
          length(
            cra_candidate_tables
          ) + 1
        ]] <- z

        cra_candidate_scores[[
          length(
            cra_candidate_scores
          ) + 1
        ]] <- z$score
      }
    }
  }
}

cra_run_search <- if (
  length(
    cra_candidate_scores
  ) == 0
) {
  tibble()
} else {
  bind_rows(
    cra_candidate_scores
  ) |>
    arrange(
      desc(
        exact_sample_overlap
      ),
      desc(
        run_accession_count
      )
    )
}

write_excel_csv(
  cra_run_search,
  file.path(
    OUT,
    "CRA002354_local_run_mapping_search.csv"
  ),
  na=""
)

# ------------------------------------------------------------
# Use best candidate table when it has sample overlap AND Runs.
# ------------------------------------------------------------

cra_technical_map <- tibble()

if (
  nrow(
    cra_run_search
  ) > 0
) {

  usable <- cra_run_search |>
    filter(
      exact_sample_overlap > 0,
      run_accession_count > 0
    )

  if (
    nrow(
      usable
    ) > 0
  ) {

    best <- usable[1, ]

    source_index <- which(
      map_lgl(
        cra_candidate_tables,
        function(z) {

          identical(
            z$score$source_file[1],
            best$source_file[1]
          ) &&
          identical(
            z$score$source_sheet[1],
            best$source_sheet[1]
          )
        }
      )
    )[1]

    z <- cra_candidate_tables[[
      source_index
    ]]

    dat <- z$data

    sample_col <-
      z$score$sample_match_column[1]

    run_col <-
      z$score$run_column[1]

    experiment_scores <- map_int(
      dat,
      count_regex_values,
      pattern="^(CRX|SRX|ERX|DRX)[0-9]+$"
    )

    experiment_col <- if (
      max(
        experiment_scores,
        na.rm=TRUE
      ) > 0
    ) {
      names(dat)[
        which.max(
          experiment_scores
        )
      ]
    } else {
      NA_character_
    }

    cra_technical_map <- tibble(
      sample_id=clean_chr(
        dat[[sample_col]]
      ),

      sample_norm=norm_id(
        dat[[sample_col]]
      ),

      run_id=clean_chr(
        dat[[run_col]]
      ),

      experiment_id=
        if (
          !is.na(
            experiment_col
          )
        ) {
          clean_chr(dat[[experiment_col]])
        } else {
          NA_character_
        },

      metadata_source_file=
        z$score$source_file[1],

      metadata_source_sheet=
        z$score$source_sheet[1]
    ) |>
      filter(
        !is.na(sample_id),
        !is.na(run_id)
      ) |>
      distinct()

    # Join by normalized public sample ID.
    cra_master <- cra_master |>
      mutate(
        sample_norm=
          norm_id(
            sample_id
          )
      ) |>
      left_join(
        cra_technical_map |>
          select(
            sample_norm,
            run_id,
            experiment_id,
            metadata_source_file,
            metadata_source_sheet
          ),
        by="sample_norm"
      )
  }
}

# ------------------------------------------------------------
# Standardized CRA clinical variables
# ------------------------------------------------------------

sepsis_col <- find_col_name(
  names(cra_master),
  c(
    "Sepsis/Septic shock"
  )
)

survival_col <- find_col_name(
  names(cra_master),
  c(
    "28 days survival"
  )
)

infection_col <- find_col_name(
  names(cra_master),
  c(
    "^Site of infection"
  )
)

pathogen_col <- find_col_name(
  names(cra_master),
  c(
    "Cultured pathogenic microorganism"
  )
)

gender_col <- find_col_name(
  names(cra_master),
  c(
    "^Gender$"
  )
)

age_col <- find_col_name(
  names(cra_master),
  c(
    "^Age"
  )
)

sofa_patient_col <- find_col_name(
  names(cra_master),
  c(
    "SOFA \\("
  )
)

apache_patient_col <- find_col_name(
  names(cra_master),
  c(
    "APACHE II \\("
  )
)

abx_icu_col <- find_col_name(
  names(cra_master),
  c(
    "Antibiotics use \\(during ICU"
  )
)

infection_label <- c(
  `0`="other",
  `1`="lung",
  `2`="intestinal",
  `3`="abdominal",
  `4`="blood",
  `5`="urinary",
  `6`="brain",
  `7`="surgical_site"
)

cra_master <- cra_master |>
  mutate(
    cohort_role=
      "CORE_SEPSIS_LONGITUDINAL",

    sepsis_status=
      if (
        !is.na(
          sepsis_col
        )
      ) {
        case_when(
          .data[[sepsis_col]]=="1" ~
            "sepsis",

          .data[[sepsis_col]]=="2" ~
            "septic_shock",

          TRUE ~
            NA_character_
        )
      } else {
        NA_character_
      },

    outcome_28d=
      if (
        !is.na(
          survival_col
        )
      ) {
        case_when(
          .data[[survival_col]]=="1" ~
            "survived",

          .data[[survival_col]]=="2" ~
            "dead",

          TRUE ~
            NA_character_
        )
      } else {
        NA_character_
      },

    infection_source=
      if (
        !is.na(
          infection_col
        )
      ) {
        unname(
          infection_label[
            clean_chr(.data[[infection_col]])
          ]
        )
      } else {
        NA_character_
      },

    cultured_pathogen=
      if (
        !is.na(
          pathogen_col
        )
      ) {
        clean_chr(.data[[pathogen_col]])
      } else {
        NA_character_
      },

    sex=
      if (
        !is.na(
          gender_col
        )
      ) {
        clean_chr(.data[[gender_col]])
      } else {
        NA_character_
      },

    age=
      if (
        !is.na(
          age_col
        )
      ) {
        suppressWarnings(
          as.numeric(.data[[age_col]])
        )
      } else {
        NA_real_
      },

    baseline_sofa=
      if (
        !is.na(
          sofa_patient_col
        )
      ) {
        suppressWarnings(
          as.numeric(.data[[sofa_patient_col]])
        )
      } else {
        NA_real_
      },

    baseline_apache_ii=
      if (
        !is.na(
          apache_patient_col
        )
      ) {
        suppressWarnings(
          as.numeric(.data[[apache_patient_col]])
        )
      } else {
        NA_real_
      },

    antibiotics_during_icu=
      if (
        !is.na(
          abx_icu_col
        )
      ) {
        clean_chr(.data[[abx_icu_col]])
      } else {
        NA_character_
      },

    patient_mapping_status=
      "EXPLICIT_PUBLIC_PATIENT_ID",

    time_mapping_status=
      "EXPLICIT_PUBLIC_COLLECTED_DAY",

    clinical_mapping_status=
      "EXPLICIT_PUBLIC_SUPPLEMENT_S1_S2_JOIN",

    metadata_source=
      "CRA002354 Supplementary Tables S1 and S2"
  )

# ============================================================
# CRA QC
# ============================================================

cra_pt <- cra_master |>
  distinct(
    patient_id,
    sample_id
  ) |>
  count(
    patient_id,
    name="n_samples"
  )

cra_qc <- tibble(
  metric=c(
    "samples",
    "patients",
    "patients_GE2",
    "patients_GE3",
    "max_timepoints",
    "S1_first_samples_exactly_found_in_S2",
    "S1_first_sample_day_matches",
    "sepsis_patients",
    "septic_shock_patients",
    "28d_survived_patients",
    "28d_dead_patients",
    "run_id_nonmissing",
    "run_id_missing"
  ),

  value=c(
    nrow(cra_master),

    n_distinct(
      cra_master$patient_id
    ),

    sum(
      cra_pt$n_samples >= 2
    ),

    sum(
      cra_pt$n_samples >= 3
    ),

    max(
      cra_pt$n_samples
    ),

    sum(
      cra_first_sample_check$
        sample_exact_match,
      na.rm=TRUE
    ),

    sum(
      cra_first_sample_check$
        day_exact_match,
      na.rm=TRUE
    ),

    n_distinct(
      cra_master$patient_id[
        cra_master$sepsis_status ==
          "sepsis"
      ]
    ),

    n_distinct(
      cra_master$patient_id[
        cra_master$sepsis_status ==
          "septic_shock"
      ]
    ),

    n_distinct(
      cra_master$patient_id[
        cra_master$outcome_28d ==
          "survived"
      ]
    ),

    n_distinct(
      cra_master$patient_id[
        cra_master$outcome_28d ==
          "dead"
      ]
    ),

    sum(
      !is.na(
        cra_master$run_id
      )
    ),

    sum(
      is.na(
        cra_master$run_id
      )
    )
  )
)

write_excel_csv(
  cra_master,
  file.path(
    OUT,
    "CRA002354_MASTER_READY_metadata.csv"
  ),
  na=""
)

write_excel_csv(
  cra_patient,
  file.path(
    OUT,
    "CRA002354_patient_level_clinical.csv"
  ),
  na=""
)

write_excel_csv(
  cra_first_sample_check,
  file.path(
    OUT,
    "CRA002354_S1_S2_linkage_QC.csv"
  ),
  na=""
)

write_excel_csv(
  cra_qc,
  file.path(
    OUT,
    "CRA002354_FINAL_QC.csv"
  ),
  na=""
)

cat(
  paste0(
    "CRA002354 COMPLETE: ",
    Sys.time(),
    "\n"
  ),
  file=LOG,
  append=TRUE
)

# ============================================================
# PART B - PRJNA851469
# ============================================================

P851_T2 <- file.path(
  STEP79,
  "PRJNA851469_sheet_02_supplementary_table_2.csv"
)

if (!file.exists(P851_T2)) {
  stop(
    "PRJNA851469 Supplementary Table 2 export not found."
  )
}

t2 <- safe_read_csv(
  P851_T2
)

if (
  is.null(t2) ||
  !("Sample ID" %in% names(t2))
) {
  stop(
    "Could not read PRJNA851469 Sample ID column."
  )
}

# ------------------------------------------------------------
# Unique microbiome samples.
#
# IMPORTANT:
# PatientXX-Day-Y is already explicit in the public supplement.
# This is NOT inference from an undocumented naming convention:
# the paper itself labels each microbiome sample this way.
# ------------------------------------------------------------

p851_samples <- t2 |>
  distinct(
    sample_id=`Sample ID`
  ) |>
  mutate(
    sample_id=
      clean_chr(
        sample_id
      ),

    sample_class=
      case_when(
        str_detect(
          sample_id,
          regex(
            "^HealthyVolunteer-",
            ignore_case=TRUE
          )
        ) ~
          "HEALTHY_CONTROL",

        str_detect(
          sample_id,
          regex(
            "^Patient[0-9]+-Day-[0-9]+$",
            ignore_case=TRUE
          )
        ) ~
          "ICU_PATIENT",

        TRUE ~
          "OTHER"
      ),

    patient_id=
      case_when(
        sample_class ==
          "ICU_PATIENT" ~
          str_extract(
            sample_id,
            regex(
              "^Patient[0-9]+",
              ignore_case=TRUE
            )
          ),

        sample_class ==
          "HEALTHY_CONTROL" ~
          sample_id,

        TRUE ~
          NA_character_
      ),

    time_raw=
      case_when(
        sample_class ==
          "ICU_PATIENT" ~
          str_extract(
            sample_id,
            regex(
              "Day-[0-9]+",
              ignore_case=TRUE
            )
          ),

        TRUE ~
          NA_character_
      ),

    time_day=
      case_when(
        sample_class ==
          "ICU_PATIENT" ~
          suppressWarnings(
            as.integer(
              str_extract(
                time_raw,
                "[0-9]+"
              )
            )
          ),

        TRUE ~
          NA_integer_
      ),

    time_order=time_day,

    cohort_role=
      case_when(
        sample_class ==
          "ICU_PATIENT" ~
          "ICU_BACKGROUND_LONGITUDINAL",

        sample_class ==
          "HEALTHY_CONTROL" ~
          "HEALTHY_STATIC_CONTROL",

        TRUE ~
          "UNCLASSIFIED"
      ),

    patient_mapping_status=
      case_when(
        sample_class ==
          "ICU_PATIENT" ~
          "EXPLICIT_PUBLIC_SAMPLE_ID",

        sample_class ==
          "HEALTHY_CONTROL" ~
          "EXPLICIT_PUBLIC_SAMPLE_ID",

        TRUE ~
          NA_character_
      ),

    time_mapping_status=
      case_when(
        sample_class ==
          "ICU_PATIENT" ~
          "EXPLICIT_PUBLIC_SAMPLE_ID_DAY",

        TRUE ~
          "NOT_APPLICABLE"
      ),

    metadata_source=
      "Nature Medicine Supplementary Table 2",

    sample_norm=
      norm_id(
        sample_id
      )
  )

# ============================================================
# Search existing LOCAL metadata for Run IDs
# ============================================================

p851_search_roots <- c(
  file.path(
    ROOT,
    "data",
    "PRJNA851469"
  ),
  file.path(
    ROOT,
    "results"
  )
)

p851_files <- character()

for (r in p851_search_roots) {

  if (!dir.exists(r)) next

  f <- list.files(
    r,
    recursive=TRUE,
    full.names=TRUE,
    pattern="\\.(csv|tsv|xlsx)$",
    ignore.case=TRUE
  )

  p851_files <- c(
    p851_files,
    f
  )
}

p851_files <- unique(
  p851_files
)

# Keep project-related result files preferentially.
p851_files <- p851_files[
  str_detect(
    p851_files,
    regex(
      "PRJNA851469",
      ignore_case=TRUE
    )
  )
]

p851_files <- p851_files[
  !str_detect(
    p851_files,
    regex(
      "V2_19_offline_metadata_parse|V2_20_finalize|Supplementary_Tables_2_17",
      ignore_case=TRUE
    )
  )
]

p851_ids <- unique(
  p851_samples$sample_norm
)

p851_candidate_tables <- list()
p851_candidate_scores <- list()

inspect_table_851 <- function(
  dat,
  source_file,
  source_sheet=NA_character_
) {

  if (
    is.null(dat) ||
    nrow(dat)==0 ||
    ncol(dat)==0
  ) {
    return(NULL)
  }

  dat <- as_tibble(dat)

  # Search exact normalized sample ID overlap.
  overlaps <- map_int(
    dat,
    function(v) {

      vv <- unique(
        norm_id(v)
      )

      length(
        intersect(
          vv[
            !is.na(vv)
          ],
          p851_ids
        )
      )
    }
  )

  # Also identify Run accession columns.
  run_scores <- map_int(
    dat,
    count_regex_values,
    pattern="^(SRR|ERR|DRR|CRR)[0-9]+$"
  )

  best_overlap <- max(
    overlaps,
    na.rm=TRUE
  )

  best_run <- max(
    run_scores,
    na.rm=TRUE
  )

  if (
    (
      !is.finite(
        best_overlap
      ) ||
      best_overlap==0
    ) &&
    (
      !is.finite(
        best_run
      ) ||
      best_run==0
    )
  ) {
    return(NULL)
  }

  best_col <- if (
    best_overlap > 0
  ) {
    names(dat)[
      which.max(
        overlaps
      )
    ]
  } else {
    NA_character_
  }

  run_col <- if (
    best_run > 0
  ) {
    names(dat)[
      which.max(
        run_scores
      )
    ]
  } else {
    NA_character_
  }

  list(
    data=dat,

    score=tibble(
      source_file=source_file,
      source_sheet=source_sheet,
      rows=nrow(dat),
      columns=ncol(dat),
      sample_match_column=best_col,
      exact_sample_overlap=best_overlap,
      run_column=run_col,
      run_accession_count=best_run
    )
  )
}

for (p in p851_files) {

  if (
    file.exists(p) &&
    file.info(p)$size >
      50 * 1024^2
  ) next

  ext <- tolower(
    tools::file_ext(p)
  )

  if (ext=="csv") {

    dat <- safe_read_csv(p)

    z <- inspect_table_851(
      dat,
      p
    )

    if (!is.null(z)) {

      p851_candidate_tables[[
        length(
          p851_candidate_tables
        ) + 1
      ]] <- z

      p851_candidate_scores[[
        length(
          p851_candidate_scores
        ) + 1
      ]] <- z$score
    }
  }

  if (ext=="tsv") {

    dat <- safe_read_tsv(p)

    z <- inspect_table_851(
      dat,
      p
    )

    if (!is.null(z)) {

      p851_candidate_tables[[
        length(
          p851_candidate_tables
        ) + 1
      ]] <- z

      p851_candidate_scores[[
        length(
          p851_candidate_scores
        ) + 1
      ]] <- z$score
    }
  }

  if (ext=="xlsx") {

    sheets <- tryCatch(
      excel_sheets(p),
      error=function(e) character()
    )

    for (sh in sheets) {

      dat <- tryCatch(
        read_excel(
          p,
          sheet=sh,
          .name_repair="unique"
        ),
        error=function(e) NULL
      )

      z <- inspect_table_851(
        dat,
        p,
        sh
      )

      if (!is.null(z)) {

        p851_candidate_tables[[
          length(
            p851_candidate_tables
          ) + 1
        ]] <- z

        p851_candidate_scores[[
          length(
            p851_candidate_scores
          ) + 1
        ]] <- z$score
      }
    }
  }
}

p851_run_search <- if (
  length(
    p851_candidate_scores
  ) == 0
) {
  tibble()
} else {
  bind_rows(
    p851_candidate_scores
  ) |>
    arrange(
      desc(
        exact_sample_overlap
      ),
      desc(
        run_accession_count
      )
    )
}

write_excel_csv(
  p851_run_search,
  file.path(
    OUT,
    "PRJNA851469_local_run_mapping_search.csv"
  ),
  na=""
)

# ------------------------------------------------------------
# Attempt exact run mapping.
# ------------------------------------------------------------

p851_map <- p851_samples

if (
  nrow(
    p851_run_search
  ) > 0
) {

  usable <- p851_run_search |>
    filter(
      exact_sample_overlap > 0,
      run_accession_count > 0
    )

  if (
    nrow(
      usable
    ) > 0
  ) {

    best <- usable[1, ]

    source_index <- which(
      map_lgl(
        p851_candidate_tables,
        function(z) {

          identical(
            z$score$source_file[1],
            best$source_file[1]
          ) &&
          identical(
            z$score$source_sheet[1],
            best$source_sheet[1]
          )
        }
      )
    )[1]

    z <- p851_candidate_tables[[
      source_index
    ]]

    dat <- z$data

    sample_col <-
      z$score$sample_match_column[1]

    run_col <-
      z$score$run_column[1]

    technical <- tibble(
      sample_norm=
        norm_id(dat[[sample_col]]),

      run_id=
        clean_chr(dat[[run_col]]),

      run_metadata_source=
        z$score$source_file[1],

      run_metadata_sheet=
        z$score$source_sheet[1]
    ) |>
      filter(
        !is.na(sample_norm),
        !is.na(run_id)
      ) |>
      distinct()

    p851_map <- p851_map |>
      left_join(
        technical,
        by="sample_norm"
      )
  }
}

# ============================================================
# PRJNA851469 QC
# ============================================================

icu851 <- p851_map |>
  filter(
    sample_class ==
      "ICU_PATIENT"
  )

pt851 <- icu851 |>
  distinct(
    patient_id,
    time_day
  ) |>
  count(
    patient_id,
    name="n_timepoints"
  )

p851_qc <- tibble(
  metric=c(
    "all_microbiome_samples",
    "healthy_control_samples",
    "ICU_microbiome_samples",
    "ICU_patients_with_public_16S_samples",
    "ICU_Day1_samples",
    "ICU_Day3_samples",
    "ICU_Day7_samples",
    "ICU_patients_GE2_timepoints",
    "ICU_patients_GE3_timepoints",
    "study_guide_reported_ICU_cohort_N",
    "difference_study_cohort_vs_public_16S_patient_count",
    "run_id_nonmissing_all_samples",
    "run_id_missing_all_samples"
  ),

  value=c(
    nrow(
      p851_map
    ),

    sum(
      p851_map$sample_class ==
        "HEALTHY_CONTROL"
    ),

    nrow(
      icu851
    ),

    n_distinct(
      icu851$patient_id
    ),

    sum(
      icu851$time_day == 1,
      na.rm=TRUE
    ),

    sum(
      icu851$time_day == 3,
      na.rm=TRUE
    ),

    sum(
      icu851$time_day == 7,
      na.rm=TRUE
    ),

    sum(
      pt851$n_timepoints >= 2
    ),

    sum(
      pt851$n_timepoints >= 3
    ),

    51,

    51 -
      n_distinct(
        icu851$patient_id
      ),

    if (
      "run_id" %in%
        names(
          p851_map
        )
    ) {
      sum(
        !is.na(
          p851_map$run_id
        )
      )
    } else {
      0
    },

    if (
      "run_id" %in%
        names(
          p851_map
        )
    ) {
      sum(
        is.na(
          p851_map$run_id
        )
      )
    } else {
      nrow(
        p851_map
      )
    }
  )
)

write_excel_csv(
  p851_map,
  file.path(
    OUT,
    "PRJNA851469_MASTER_READY_metadata.csv"
  ),
  na=""
)

write_excel_csv(
  icu851,
  file.path(
    OUT,
    "PRJNA851469_ICU_LONGITUDINAL_ONLY.csv"
  ),
  na=""
)

write_excel_csv(
  p851_map |>
    filter(
      sample_class ==
        "HEALTHY_CONTROL"
    ),
  file.path(
    OUT,
    "PRJNA851469_HEALTHY_CONTROLS.csv"
  ),
  na=""
)

write_excel_csv(
  p851_qc,
  file.path(
    OUT,
    "PRJNA851469_FINAL_QC.csv"
  ),
  na=""
)

cat(
  paste0(
    "PRJNA851469 COMPLETE: ",
    Sys.time(),
    "\n"
  ),
  file=LOG,
  append=TRUE
)

# ============================================================
# FINAL DECISION TABLE
# ============================================================

cra_run_n <- sum(
  !is.na(
    cra_master$run_id
  )
)

p851_run_n <- if (
  "run_id" %in%
    names(
      p851_map
    )
) {
  sum(
    !is.na(
      p851_map$run_id
    )
  )
} else {
  0
}

decision <- tibble(
  project=c(
    "CRA002354",
    "PRJNA851469"
  ),

  patient_time_mapping=c(
    "EXPLICIT_CONFIRMED",
    "EXPLICIT_CONFIRMED"
  ),

  longitudinal_patient_n=c(
    sum(
      cra_pt$n_samples >= 2
    ),

    sum(
      pt851$n_timepoints >= 2
    )
  ),

  clinical_depth=c(
    paste(
      "Excellent: sepsis/shock, 28d survival, infection site,",
      "pathogen, SOFA, APACHE II, lactate, antibiotics"
    ),

    paste(
      "Limited in public supplement used here;",
      "patient/time are explicit but individual infection/outcome",
      "labels are not recovered from Supplementary Table 2"
    )
  ),

  run_mapping_status=c(
    ifelse(
      cra_run_n ==
        nrow(
          cra_master
        ),
      "COMPLETE",
      paste0(
        "PARTIAL_",
        cra_run_n,
        "_OF_",
        nrow(
          cra_master
        )
      )
    ),

    ifelse(
      p851_run_n ==
        nrow(
          p851_map
        ),
      "COMPLETE",
      paste0(
        "PARTIAL_",
        p851_run_n,
        "_OF_",
        nrow(
          p851_map
        )
      )
    )
  ),

  recommended_role=c(
    "CORE_SEPSIS_LONGITUDINAL_WITH_OUTCOME_AND_INFECTION_SOURCE",
    "ICU_BACKGROUND_LONGITUDINAL"
  ),

  metadata_decision=c(
    "INCLUDE",
    "INCLUDE"
  )
)

write_excel_csv(
  decision,
  file.path(
    OUT,
    "V2_step80_cohort_decision.csv"
  ),
  na=""
)

cat(
  paste0(
    "STEP80 COMPLETE: ",
    Sys.time(),
    "\n"
  ),
  file=LOG,
  append=TRUE
)

cat("\n")
cat("============================================================\n")
cat("SEPSIS V2 - STEP 80 COMPLETE\n")
cat("============================================================\n")
cat("\nCRA002354 QC:\n")
print(
  cra_qc,
  n=Inf,
  width=Inf
)

cat("\nPRJNA851469 QC:\n")
print(
  p851_qc,
  n=Inf,
  width=Inf
)

cat("\nCOHORT DECISION:\n")
print(
  decision,
  n=Inf,
  width=Inf
)

cat("\nOutput folder:\n")
cat(
  OUT,
  "\n"
)
cat("============================================================\n")
