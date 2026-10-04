# ============================================================
# Sepsis V2 - Step 79
# OFFLINE metadata parsing for:
#   1) PRJNA851469
#   2) CRA002354
#
# Purpose:
# - No internet access
# - Read the validated supplementary XLSX files downloaded in 78A3/78A4B
# - Profile every worksheet
# - Detect patient/sample/time/infection/outcome/antibiotic fields
# - Export candidate metadata tables for exact adjudication
# - Do NOT invent patient/time mappings
#
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
# Paths
# ============================================================

PROJECT_ROOT <- "E:/sepsis_project"

OUT_ROOT <- file.path(
  PROJECT_ROOT,
  "results",
  "V2_19_offline_metadata_parse"
)

dir.create(
  OUT_ROOT,
  recursive=TRUE,
  showWarnings=FALSE
)

LOG <- file.path(
  OUT_ROOT,
  "_STEP79_runtime_checkpoints.txt"
)

cat(
  paste0("STEP79 STARTED: ", Sys.time(), "\n"),
  file=LOG,
  append=TRUE
)

PRJNA851469_XLSX <- file.path(
  PROJECT_ROOT,
  "data",
  "_step78A3_US_fast_download",
  "PRJNA851469",
  "PRJNA851469_Supplementary_Tables_2_17.xlsx"
)

CRA_DIR <- file.path(
  PROJECT_ROOT,
  "data",
  "_step78A3_US_fast_download",
  "CRA002354"
)

CRA_S1_CANDIDATES <- c(
  file.path(CRA_DIR, "CRA002354_Table_S1_64_patients.xlsx"),
  file.path(CRA_DIR, "mmc5.xlsx")
)

CRA_S2_CANDIDATES <- c(
  file.path(CRA_DIR, "CRA002354_Table_S2_131_samples.xlsx"),
  file.path(CRA_DIR, "mmc6.xlsx")
)

first_existing <- function(paths) {
  hit <- paths[file.exists(paths)]
  if (length(hit) == 0) {
    return(NA_character_)
  }
  hit[1]
}

CRA_S1_XLSX <- first_existing(CRA_S1_CANDIDATES)
CRA_S2_XLSX <- first_existing(CRA_S2_CANDIDATES)

# ============================================================
# Helpers
# ============================================================

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

safe_name <- function(x) {
  x <- iconv(
    x,
    from="",
    to="ASCII//TRANSLIT"
  )

  x <- tolower(x)

  x <- str_replace_all(
    x,
    "[^a-z0-9]+",
    "_"
  )

  x <- str_replace_all(
    x,
    "^_+|_+$",
    ""
  )

  ifelse(
    is.na(x) | x == "",
    "unnamed",
    x
  )
}

validate_xlsx <- function(path) {

  if (
    is.na(path) ||
    !file.exists(path)
  ) {
    return(FALSE)
  }

  if (
    file.info(path)$size < 1000
  ) {
    return(FALSE)
  }

  tryCatch(
    {
      sh <- excel_sheets(path)
      length(sh) > 0
    },
    error=function(e) FALSE
  )
}

read_sheet_safe <- function(
  path,
  sheet,
  skip=0
) {

  tryCatch(
    {
      x <- read_excel(
        path,
        sheet=sheet,
        skip=skip,
        .name_repair="unique"
      )

      as_tibble(x)
    },
    error=function(e) {
      NULL
    }
  )
}

nonempty_fraction <- function(x) {

  x <- clean_chr(x)

  if (length(x) == 0) {
    return(0)
  }

  mean(!is.na(x))
}

n_unique_nonmissing <- function(x) {

  x <- clean_chr(x)

  length(
    unique(
      x[
        !is.na(x)
      ]
    )
  )
}

keyword_groups <- list(

  patient = c(
    "patient",
    "subject",
    "participant",
    "individual",
    "case",
    "id"
  ),

  sample = c(
    "sample",
    "biosample",
    "specimen",
    "swab",
    "stool",
    "fecal",
    "faecal",
    "run",
    "accession"
  ),

  time = c(
    "time",
    "day",
    "visit",
    "admission",
    "discharge",
    "baseline",
    "infection",
    "onset",
    "date"
  ),

  infection = c(
    "infection",
    "sepsis",
    "septic",
    "nosocomial",
    "culture",
    "pathogen",
    "bacteremia",
    "bacteraemia",
    "vap",
    "pneumonia"
  ),

  outcome = c(
    "outcome",
    "mortality",
    "death",
    "dead",
    "survival",
    "survivor",
    "discharge",
    "icu_free",
    "hospital_free"
  ),

  antibiotics = c(
    "antibiotic",
    "antimicrobial",
    "abx",
    "carbapenem",
    "vancomycin",
    "meropenem",
    "piperacillin"
  ),

  severity = c(
    "sofa",
    "apache",
    "shock",
    "vasopressor",
    "organ",
    "severity"
  )
)

group_regex <- lapply(
  keyword_groups,
  function(x) {
    regex(
      paste(
        x,
        collapse="|"
      ),
      ignore_case=TRUE
    )
  }
)

column_group_flags <- function(nms) {

  tibble(
    column_name=nms,
    column_name_clean=safe_name(nms),

    patient_field=
      str_detect(
        nms,
        group_regex$patient
      ),

    sample_field=
      str_detect(
        nms,
        group_regex$sample
      ),

    time_field=
      str_detect(
        nms,
        group_regex$time
      ),

    infection_field=
      str_detect(
        nms,
        group_regex$infection
      ),

    outcome_field=
      str_detect(
        nms,
        group_regex$outcome
      ),

    antibiotics_field=
      str_detect(
        nms,
        group_regex$antibiotics
      ),

    severity_field=
      str_detect(
        nms,
        group_regex$severity
      )
  )
}

sheet_profile <- function(
  path,
  project,
  source_label
) {

  sheets <- excel_sheets(path)

  out <- vector(
    "list",
    length(sheets)
  )

  for (i in seq_along(sheets)) {

    sh <- sheets[i]

    x <- read_sheet_safe(
      path,
      sh
    )

    if (is.null(x)) {

      out[[i]] <- tibble(
        project=project,
        source=source_label,
        sheet_index=i,
        sheet=sh,
        read_ok=FALSE,
        rows=NA_integer_,
        columns=NA_integer_,
        patient_fields=NA_integer_,
        sample_fields=NA_integer_,
        time_fields=NA_integer_,
        infection_fields=NA_integer_,
        outcome_fields=NA_integer_,
        antibiotics_fields=NA_integer_,
        severity_fields=NA_integer_,
        keyword_score=NA_integer_
      )

      next
    }

    flags <- column_group_flags(
      names(x)
    )

    score <-
      sum(flags$patient_field) +
      sum(flags$sample_field) +
      sum(flags$time_field) +
      sum(flags$infection_field) +
      sum(flags$outcome_field) +
      sum(flags$antibiotics_field) +
      sum(flags$severity_field)

    out[[i]] <- tibble(
      project=project,
      source=source_label,
      sheet_index=i,
      sheet=sh,
      read_ok=TRUE,
      rows=nrow(x),
      columns=ncol(x),

      patient_fields=
        sum(flags$patient_field),

      sample_fields=
        sum(flags$sample_field),

      time_fields=
        sum(flags$time_field),

      infection_fields=
        sum(flags$infection_field),

      outcome_fields=
        sum(flags$outcome_field),

      antibiotics_fields=
        sum(flags$antibiotics_field),

      severity_fields=
        sum(flags$severity_field),

      keyword_score=score,

      column_names=
        paste(
          names(x),
          collapse=" | "
        )
    )
  }

  bind_rows(out)
}

export_workbook <- function(
  path,
  project,
  prefix
) {

  sheets <- excel_sheets(path)

  sheet_inventory <- list()
  column_inventory <- list()

  for (i in seq_along(sheets)) {

    sh <- sheets[i]

    x <- read_sheet_safe(
      path,
      sh
    )

    if (is.null(x)) next

    # Export every worksheet in full.
    safe_sh <- safe_name(sh)

    out_csv <- file.path(
      OUT_ROOT,
      paste0(
        prefix,
        "_sheet_",
        sprintf("%02d", i),
        "_",
        safe_sh,
        ".csv"
      )
    )

    write_excel_csv(
      x,
      out_csv,
      na=""
    )

    flags <- column_group_flags(
      names(x)
    ) |>
      mutate(
        project=project,
        source_file=basename(path),
        sheet=sh,
        sheet_index=i,
        .before=1
      )

    stats <- tibble(
      column_name=names(x),

      nonempty_fraction=
        map_dbl(
          x,
          nonempty_fraction
        ),

      unique_nonmissing=
        map_int(
          x,
          n_unique_nonmissing
        )
    )

    flags <- flags |>
      left_join(
        stats,
        by="column_name"
      )

    column_inventory[[
      length(column_inventory)+1
    ]] <- flags

    sheet_inventory[[
      length(sheet_inventory)+1
    ]] <- tibble(
      project=project,
      source_file=basename(path),
      sheet=sh,
      sheet_index=i,
      rows=nrow(x),
      columns=ncol(x),
      exported_csv=out_csv
    )
  }

  list(
    sheets=
      bind_rows(
        sheet_inventory
      ),

    columns=
      bind_rows(
        column_inventory
      )
  )
}

extract_candidate_tables <- function(
  path,
  project,
  prefix
) {

  sheets <- excel_sheets(path)

  outputs <- list()

  for (i in seq_along(sheets)) {

    sh <- sheets[i]

    x <- read_sheet_safe(
      path,
      sh
    )

    if (
      is.null(x) ||
      nrow(x) == 0 ||
      ncol(x) == 0
    ) next

    flags <- column_group_flags(
      names(x)
    )

    keep <- flags |>
      filter(
        patient_field |
        sample_field |
        time_field |
        infection_field |
        outcome_field |
        antibiotics_field |
        severity_field
      ) |>
      pull(column_name)

    if (length(keep) == 0) next

    candidate <- x |>
      select(
        all_of(keep)
      ) |>
      mutate(
        Source_Project=project,
        Source_Workbook=basename(path),
        Source_Sheet=sh,
        Source_Row=row_number(),
        .before=1
      )

    safe_sh <- safe_name(sh)

    file_out <- file.path(
      OUT_ROOT,
      paste0(
        prefix,
        "_CANDIDATE_METADATA_",
        sprintf("%02d", i),
        "_",
        safe_sh,
        ".csv"
      )
    )

    write_excel_csv(
      candidate,
      file_out,
      na=""
    )

    outputs[[
      length(outputs)+1
    ]] <- tibble(
      project=project,
      workbook=basename(path),
      sheet_index=i,
      sheet=sh,
      candidate_columns=length(keep),
      candidate_rows=nrow(candidate),
      output_csv=file_out
    )
  }

  if (length(outputs)==0) {
    tibble()
  } else {
    bind_rows(outputs)
  }
}

# ============================================================
# Validate inputs
# ============================================================

input_status <- tibble(
  project=c(
    "PRJNA851469",
    "CRA002354",
    "CRA002354"
  ),

  source=c(
    "Nature_Supplement_Tables_2_17",
    "Supplement_Table_S1_patient_level",
    "Supplement_Table_S2_sample_level"
  ),

  path=c(
    PRJNA851469_XLSX,
    CRA_S1_XLSX,
    CRA_S2_XLSX
  )
) |>
  mutate(
    exists=
      !is.na(path) &
      file.exists(path),

    bytes=
      ifelse(
        exists,
        file.info(path)$size,
        NA_real_
      ),

    valid_xlsx=
      map_lgl(
        path,
        validate_xlsx
      )
  )

write_excel_csv(
  input_status,
  file.path(
    OUT_ROOT,
    "STEP79_input_validation.csv"
  ),
  na=""
)

if (!all(input_status$valid_xlsx)) {

  print(input_status)

  stop(
    paste(
      "At least one required XLSX is missing or invalid.",
      "See STEP79_input_validation.csv."
    )
  )
}

cat(
  paste0(
    "INPUT VALIDATION COMPLETE: ",
    Sys.time(),
    "\n"
  ),
  file=LOG,
  append=TRUE
)

# ============================================================
# PRJNA851469
# ============================================================

profile851 <- sheet_profile(
  PRJNA851469_XLSX,
  "PRJNA851469",
  "Nature_Supplement_Tables_2_17"
)

write_excel_csv(
  profile851,
  file.path(
    OUT_ROOT,
    "PRJNA851469_sheet_profile.csv"
  ),
  na=""
)

export851 <- export_workbook(
  PRJNA851469_XLSX,
  "PRJNA851469",
  "PRJNA851469"
)

write_excel_csv(
  export851$sheets,
  file.path(
    OUT_ROOT,
    "PRJNA851469_sheet_inventory.csv"
  ),
  na=""
)

write_excel_csv(
  export851$columns,
  file.path(
    OUT_ROOT,
    "PRJNA851469_column_inventory.csv"
  ),
  na=""
)

candidate851 <- extract_candidate_tables(
  PRJNA851469_XLSX,
  "PRJNA851469",
  "PRJNA851469"
)

write_excel_csv(
  candidate851,
  file.path(
    OUT_ROOT,
    "PRJNA851469_candidate_table_index.csv"
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
# CRA002354 - Supplement Table S1
# ============================================================

profileCRA_S1 <- sheet_profile(
  CRA_S1_XLSX,
  "CRA002354",
  "Supplement_Table_S1_patient_level"
)

write_excel_csv(
  profileCRA_S1,
  file.path(
    OUT_ROOT,
    "CRA002354_S1_sheet_profile.csv"
  ),
  na=""
)

exportCRA_S1 <- export_workbook(
  CRA_S1_XLSX,
  "CRA002354",
  "CRA002354_S1"
)

write_excel_csv(
  exportCRA_S1$sheets,
  file.path(
    OUT_ROOT,
    "CRA002354_S1_sheet_inventory.csv"
  ),
  na=""
)

write_excel_csv(
  exportCRA_S1$columns,
  file.path(
    OUT_ROOT,
    "CRA002354_S1_column_inventory.csv"
  ),
  na=""
)

candidateCRA_S1 <- extract_candidate_tables(
  CRA_S1_XLSX,
  "CRA002354",
  "CRA002354_S1"
)

write_excel_csv(
  candidateCRA_S1,
  file.path(
    OUT_ROOT,
    "CRA002354_S1_candidate_table_index.csv"
  ),
  na=""
)

cat(
  paste0(
    "CRA002354 S1 COMPLETE: ",
    Sys.time(),
    "\n"
  ),
  file=LOG,
  append=TRUE
)

# ============================================================
# CRA002354 - Supplement Table S2
# ============================================================

profileCRA_S2 <- sheet_profile(
  CRA_S2_XLSX,
  "CRA002354",
  "Supplement_Table_S2_sample_level"
)

write_excel_csv(
  profileCRA_S2,
  file.path(
    OUT_ROOT,
    "CRA002354_S2_sheet_profile.csv"
  ),
  na=""
)

exportCRA_S2 <- export_workbook(
  CRA_S2_XLSX,
  "CRA002354",
  "CRA002354_S2"
)

write_excel_csv(
  exportCRA_S2$sheets,
  file.path(
    OUT_ROOT,
    "CRA002354_S2_sheet_inventory.csv"
  ),
  na=""
)

write_excel_csv(
  exportCRA_S2$columns,
  file.path(
    OUT_ROOT,
    "CRA002354_S2_column_inventory.csv"
  ),
  na=""
)

candidateCRA_S2 <- extract_candidate_tables(
  CRA_S2_XLSX,
  "CRA002354",
  "CRA002354_S2"
)

write_excel_csv(
  candidateCRA_S2,
  file.path(
    OUT_ROOT,
    "CRA002354_S2_candidate_table_index.csv"
  ),
  na=""
)

cat(
  paste0(
    "CRA002354 S2 COMPLETE: ",
    Sys.time(),
    "\n"
  ),
  file=LOG,
  append=TRUE
)

# ============================================================
# CRA S1 vs S2 cross-table linkage diagnostics
#
# This DOES NOT assume which column is the patient ID.
# It only computes exact overlap between plausible ID-like columns
# in S1 and S2 so we can see whether a direct public linkage exists.
# ============================================================

read_all_sheets <- function(path) {

  sheets <- excel_sheets(path)

  out <- list()

  for (i in seq_along(sheets)) {

    x <- read_sheet_safe(
      path,
      sheets[i]
    )

    if (is.null(x)) next

    out[[length(out)+1]] <- list(
      sheet=sheets[i],
      data=x
    )
  }

  out
}

cra_s1_all <- read_all_sheets(
  CRA_S1_XLSX
)

cra_s2_all <- read_all_sheets(
  CRA_S2_XLSX
)

linkage_results <- list()

for (a in cra_s1_all) {

  x1 <- a$data

  flags1 <- column_group_flags(
    names(x1)
  )

  id1 <- flags1 |>
    filter(
      patient_field |
      sample_field
    ) |>
    pull(column_name)

  if (length(id1)==0) next

  for (b in cra_s2_all) {

    x2 <- b$data

    flags2 <- column_group_flags(
      names(x2)
    )

    id2 <- flags2 |>
      filter(
        patient_field |
        sample_field
      ) |>
      pull(column_name)

    if (length(id2)==0) next

    for (c1 in id1) {

      v1 <- unique(
        clean_chr(
          x1[[c1]]
        )
      )

      v1 <- v1[
        !is.na(v1)
      ]

      if (length(v1)==0) next

      for (c2 in id2) {

        v2 <- unique(
          clean_chr(
            x2[[c2]]
          )
        )

        v2 <- v2[
          !is.na(v2)
        ]

        if (length(v2)==0) next

        common <- intersect(
          v1,
          v2
        )

        linkage_results[[
          length(linkage_results)+1
        ]] <- tibble(
          S1_sheet=a$sheet,
          S1_column=c1,
          S1_unique=length(v1),

          S2_sheet=b$sheet,
          S2_column=c2,
          S2_unique=length(v2),

          exact_overlap_n=length(common),

          overlap_fraction_S1=
            length(common) /
            length(v1),

          overlap_fraction_S2=
            length(common) /
            length(v2),

          example_overlap=
            paste(
              head(
                common,
                10
              ),
              collapse=" | "
            )
        )
      }
    }
  }
}

cra_linkage <- if (
  length(linkage_results)==0
) {
  tibble()
} else {
  bind_rows(
    linkage_results
  ) |>
    arrange(
      desc(exact_overlap_n),
      desc(overlap_fraction_S1),
      desc(overlap_fraction_S2)
    )
}

write_excel_csv(
  cra_linkage,
  file.path(
    OUT_ROOT,
    "CRA002354_S1_S2_exact_ID_overlap_diagnostic.csv"
  ),
  na=""
)

# ============================================================
# Summary
# ============================================================

summary <- tibble(
  Project=c(
    "PRJNA851469",
    "CRA002354_S1",
    "CRA002354_S2"
  ),

  Workbook=c(
    basename(PRJNA851469_XLSX),
    basename(CRA_S1_XLSX),
    basename(CRA_S2_XLSX)
  ),

  Sheets=c(
    nrow(profile851),
    nrow(profileCRA_S1),
    nrow(profileCRA_S2)
  ),

  Total_Rows_Across_Sheets=c(
    sum(profile851$rows, na.rm=TRUE),
    sum(profileCRA_S1$rows, na.rm=TRUE),
    sum(profileCRA_S2$rows, na.rm=TRUE)
  ),

  Candidate_Metadata_Tables=c(
    nrow(candidate851),
    nrow(candidateCRA_S1),
    nrow(candidateCRA_S2)
  ),

  Status=c(
    "OFFLINE_PARSE_COMPLETE",
    "OFFLINE_PARSE_COMPLETE",
    "OFFLINE_PARSE_COMPLETE"
  )
)

write_excel_csv(
  summary,
  file.path(
    OUT_ROOT,
    "V2_step79_status.csv"
  ),
  na=""
)

cat(
  paste0(
    "STEP79 COMPLETE: ",
    Sys.time(),
    "\n"
  ),
  file=LOG,
  append=TRUE
)

cat("\n")
cat("============================================================\n")
cat("SEPSIS V2 - STEP 79 COMPLETE\n")
cat("============================================================\n")
print(summary, n=Inf, width=Inf)
cat("\nOutput folder:\n")
cat(OUT_ROOT, "\n")
cat("============================================================\n")
