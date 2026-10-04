# ============================================================
# Sepsis V2 - Step 91A
# Infection-source / infection-site feasibility audit
#
# PURPOSE
# Before any source-specific microbiome analysis, audit the
# frozen Step87B metadata to determine which cohorts contain
# defensible infection-source variables.
#
# This step DOES NOT create a source classification automatically.
# It identifies candidate metadata columns, distributions,
# missingness, patient-level consistency, and feasibility tiers.
#
# Examples of source/site concepts searched:
# pneumonia / respiratory / lung
# urinary / UTI
# abdominal / intra-abdominal
# bloodstream / bacteremia
# catheter
# skin / soft tissue
# CNS / meningitis
# other infection
#
# IMPORTANT:
# - "phenotype", "diagnosis", "disease" etc. are only candidate
#   fields and are NOT assumed to be infection source.
# - Manual/source-paper validation remains mandatory before
#   freezing any source variable.
# ============================================================

options(stringsAsFactors = FALSE)

pkgs <- c(
  "readr",
  "dplyr",
  "tidyr",
  "tibble",
  "purrr",
  "stringr"
)

missing_pkgs <- pkgs[
  !vapply(
    pkgs,
    requireNamespace,
    logical(1),
    quietly = TRUE
  )
]

if (length(missing_pkgs)) {
  install.packages(
    missing_pkgs,
    repos = "https://cloud.r-project.org"
  )
}

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(purrr)
  library(stringr)
})

ROOT <- "E:/sepsis_project"

STEP87B <- file.path(
  ROOT,
  "results",
  "V2_27B_ANALYSIS_OBJECTS_AND_REPLICATE_FREEZE"
)

STEP90B <- file.path(
  ROOT,
  "results",
  "V2_30B_STEP90B_CORRECTED_COHERENCE_INFERENCE"
)

OUT <- file.path(
  ROOT,
  "results",
  "V2_31A_STEP91A_INFECTION_SOURCE_FEASIBILITY_AUDIT"
)

DIR_COL <- file.path(
  OUT,
  "01_COLUMN_INVENTORY"
)

DIR_VAL <- file.path(
  OUT,
  "02_CANDIDATE_VALUES"
)

DIR_PAT <- file.path(
  OUT,
  "03_PATIENT_CONSISTENCY"
)

DIR_REG <- file.path(
  OUT,
  "04_FEASIBILITY_REGISTRY"
)

for (d in c(
  OUT,
  DIR_COL,
  DIR_VAL,
  DIR_PAT,
  DIR_REG
)) {
  dir.create(
    d,
    recursive = TRUE,
    showWarnings = FALSE
  )
}

LOG <- file.path(
  OUT,
  "_STEP91A_runtime_checkpoints.txt"
)

ERR <- file.path(
  OUT,
  "_STEP91A_FATAL_ERROR.txt"
)

if (file.exists(ERR)) {
  unlink(ERR)
}

ck <- function(x) {
  cat(
    paste0(
      x,
      ": ",
      Sys.time(),
      "\n"
    ),
    file = LOG,
    append = TRUE
  )
}

safe_csv <- function(p) {
  suppressMessages(
    read_csv(
      p,
      show_col_types = FALSE,
      progress = FALSE,
      name_repair = "minimal"
    )
  )
}

write_csv_safe <- function(x, p) {
  write_excel_csv(
    x,
    p,
    na = ""
  )
}

clean_chr <- function(x) {

  y <- as.character(x)
  y <- str_trim(y)

  y[
    is.na(y) |
    y == "" |
    tolower(y) %in%
      c(
        "na",
        "nan",
        "null",
        "unknown",
        "not available",
        "not_applicable",
        "n/a"
      )
  ] <- NA_character_

  y
}

is_identifier_like <- function(name) {

  str_detect(
    tolower(name),
    paste0(
      "(^|_)(run|sample|biosample|patient|subject|study|project|",
      "accession|id|file|path|read|library|batch|plate|",
      "time|day|date|visit|order|axis|replicate)($|_)"
    )
  )
}

candidate_name_score <- function(name) {

  x <- tolower(name)

  score <- 0L

  if (
    str_detect(
      x,
      "infection.?source|source.?infection|infection.?site|site.?infection|sepsis.?source|focus"
    )
  ) {
    score <- score + 6L
  }

  if (
    str_detect(
      x,
      "source|site|origin"
    )
  ) {
    score <- score + 3L
  }

  if (
    str_detect(
      x,
      "phenotype|diagnos|disease|condition|clinical|infection|sepsis"
    )
  ) {
    score <- score + 2L
  }

  if (
    str_detect(
      x,
      "pneum|respir|pulmon|lung|urinar|uti|uro|abdom|blood|bacter|catheter|skin|soft|mening|cns"
    )
  ) {
    score <- score + 4L
  }

  if (
    is_identifier_like(name)
  ) {
    score <- score - 5L
  }

  score
}

source_value_regex <- paste0(
  "(",
  "pneumonia|pneumon|respiratory|pulmonary|lung|",
  "urinary|urine|\\buti\\b|urosepsis|urogenital|",
  "abdominal|intra[-_ ]?abdominal|gastrointestinal|periton|biliary|cholang|",
  "bloodstream|blood stream|bacteremia|bacteraemia|septicemia|septicaemia|",
  "catheter|central line|",
  "skin|soft tissue|wound|cellulitis|",
  "mening|central nervous|\\bcns\\b|",
  "other infection|infection source|infection site",
  ")"
)

value_source_fraction <- function(x) {

  y <- clean_chr(x)

  y <- y[
    !is.na(y)
  ]

  if (!length(y)) {
    return(NA_real_)
  }

  mean(
    str_detect(
      tolower(y),
      source_value_regex
    )
  )
}

summarize_candidate_column <- function(
  md,
  project,
  col
) {

  x <- clean_chr(
    md[[col]]
  )

  nonmissing <- !is.na(x)

  n_nonmissing <- sum(nonmissing)

  n_unique <- n_distinct(
    x[nonmissing]
  )

  patient_nonmissing <- md |>
    mutate(
      .value =
        x
    ) |>
    filter(
      !is.na(.value)
    ) |>
    distinct(
      patient_id
    ) |>
    nrow()

  patient_total <- n_distinct(
    md$patient_id
  )

  source_frac <- value_source_fraction(
    x
  )

  tibble(
    project =
      project,
    column =
      col,
    name_score =
      candidate_name_score(
        col
      ),
    n_samples =
      nrow(md),
    n_nonmissing_samples =
      n_nonmissing,
    sample_nonmissing_fraction =
      n_nonmissing /
        nrow(md),
    n_patients_total =
      patient_total,
    n_patients_nonmissing =
      patient_nonmissing,
    patient_nonmissing_fraction =
      patient_nonmissing /
        patient_total,
    n_unique_nonmissing_values =
      n_unique,
    source_like_value_fraction =
      source_frac
  )
}

patient_consistency <- function(
  md,
  project,
  col
) {

  d <- md |>
    mutate(
      .value =
        clean_chr(
          .data[[col]]
        )
    ) |>
    filter(
      !is.na(
        .value
      )
    )

  if (!nrow(d)) {
    return(
      tibble(
        project =
          project,
        column =
          col,
        n_patients_with_value =
          0L,
        n_patients_one_unique_value =
          0L,
        n_patients_multiple_values =
          0L,
        patient_consistency_fraction =
          NA_real_
      )
    )
  }

  p <- d |>
    distinct(
      patient_id,
      .value
    ) |>
    count(
      patient_id,
      name =
        "n_unique_values"
    )

  tibble(
    project =
      project,
    column =
      col,
    n_patients_with_value =
      nrow(p),
    n_patients_one_unique_value =
      sum(
        p$n_unique_values ==
          1
      ),
    n_patients_multiple_values =
      sum(
        p$n_unique_values >
          1
      ),
    patient_consistency_fraction =
      mean(
        p$n_unique_values ==
          1
      )
  )
}

assign_feasibility <- function(d) {

  d |>
    mutate(
      automated_feasibility_tier =
        case_when(
          name_score >=
            5 &
            n_patients_nonmissing >=
              20 &
            patient_nonmissing_fraction >=
              0.70 &
            n_unique_nonmissing_values >=
              2 &
            n_unique_nonmissing_values <=
              12 &
            source_like_value_fraction >=
              0.50 ~
            "HIGH_PRIORITY_MANUAL_VALIDATION",

          (
            name_score >=
              3 |
            source_like_value_fraction >=
              0.30
          ) &
            n_patients_nonmissing >=
              10 &
            n_unique_nonmissing_values >=
              2 &
            n_unique_nonmissing_values <=
              20 ~
            "POTENTIAL_MANUAL_REVIEW",

          TRUE ~
            "LOW_OR_UNCERTAIN"
        ),
      manual_validation_required =
        TRUE
    )
}

main <- function() {

  ck(
    "STEP91A STARTED"
  )

  registry_path <- file.path(
    STEP87B,
    "V2_STEP87B_analysis_object_registry.csv"
  )

  step87_complete <- file.path(
    STEP87B,
    "_STEP87B_COMPLETE.ok"
  )

  step90_complete <- file.path(
    STEP90B,
    "_STEP90B_COMPLETE.ok"
  )

  required <- c(
    registry_path,
    step87_complete,
    step90_complete
  )

  if (
    !all(
      file.exists(
        required
      )
    )
  ) {
    stop(
      paste0(
        "Required upstream files missing:\n",
        paste(
          required[
            !file.exists(
              required
            )
          ],
          collapse = "\n"
        )
      )
    )
  }

  registry <- safe_csv(
    registry_path
  )

  if (
    nrow(registry) !=
      8 ||
    sum(
      registry$primary_samples
    ) !=
      785
  ) {
    stop(
      "Step87B registry guard failed."
    )
  }

  ck(
    "UPSTREAM STATE GUARDED"
  )

  inventory_rows <- list()
  candidate_rows <- list()
  value_rows <- list()
  consistency_rows <- list()

  for (
    proj in registry$project
  ) {

    ck(
      paste0(
        proj,
        " START"
      )
    )

    rr <- registry |>
      filter(
        project ==
          proj
      )

    object_path <- as.character(
      rr$analysis_object_path[1]
    )

    if (
      !file.exists(
        object_path
      )
    ) {
      stop(
        proj,
        ": analysis object missing: ",
        object_path
      )
    }

    obj <- readRDS(
      object_path
    )

    md <- as_tibble(
      obj$metadata
    )

    if (
      !("patient_id" %in%
          names(md))
    ) {
      stop(
        proj,
        ": metadata has no patient_id."
      )
    }

    inv <- tibble(
      project =
        proj,
      column =
        names(md),
      class =
        vapply(
          md,
          function(x)
            paste(
              class(x),
              collapse = ";"
            ),
          character(1)
        ),
      name_score =
        vapply(
          names(md),
          candidate_name_score,
          integer(1)
        ),
      identifier_like =
        vapply(
          names(md),
          is_identifier_like,
          logical(1)
        )
    )

    inventory_rows[[
      length(
        inventory_rows
      ) + 1
    ]] <- inv

    # Candidate columns are retained if either the field name
    # is suggestive OR the values themselves contain source-like
    # concepts. This broad screen is for audit only.
    candidate_cols <- character(0)

    for (
      col in names(md)
    ) {

      if (
        col ==
          "patient_id"
      ) {
        next
      }

      x <- md[[col]]

      # Skip list/matrix columns if ever present.
      if (
        is.list(x) ||
        is.matrix(x)
      ) {
        next
      }

      nscore <- candidate_name_score(
        col
      )

      vfrac <- value_source_fraction(
        x
      )

      if (
        nscore >= 2 ||
        (
          !is.na(vfrac) &&
          vfrac >= 0.15
        )
      ) {
        candidate_cols <- c(
          candidate_cols,
          col
        )
      }
    }

    candidate_cols <- unique(
      candidate_cols
    )

    if (
      length(
        candidate_cols
      )
    ) {

      for (
        col in candidate_cols
      ) {

        sm <- summarize_candidate_column(
          md,
          proj,
          col
        )

        candidate_rows[[
          length(
            candidate_rows
          ) + 1
        ]] <- sm

        consistency_rows[[
          length(
            consistency_rows
          ) + 1
        ]] <- patient_consistency(
          md,
          proj,
          col
        )

        vals <- md |>
          transmute(
            project =
              proj,
            column =
              col,
            patient_id =
              patient_id,
            value =
              clean_chr(
                .data[[col]]
              )
          ) |>
          filter(
            !is.na(value)
          )

        counts <- vals |>
          count(
            value,
            sort = TRUE,
            name =
              "n_samples"
          ) |>
          mutate(
            n_patients =
              vapply(
                value,
                function(v)
                  n_distinct(
                    vals$patient_id[
                      vals$value ==
                        v
                    ]
                  ),
                integer(1)
              ),
            source_like =
              str_detect(
                tolower(value),
                source_value_regex
              ),
            .before = 1
          ) |>
          mutate(
            project =
              proj,
            column =
              col,
            .before = 1
          )

        value_rows[[
          length(
            value_rows
          ) + 1
        ]] <- counts
      }
    }

    ck(
      paste0(
        proj,
        " COMPLETE"
      )
    )
  }

  inventory <- bind_rows(
    inventory_rows
  )

  candidates <- bind_rows(
    candidate_rows
  )

  values <- bind_rows(
    value_rows
  )

  consistency <- bind_rows(
    consistency_rows
  )

  candidates2 <- candidates |>
    left_join(
      consistency,
      by = c(
        "project",
        "column"
      )
    ) |>
    assign_feasibility() |>
    arrange(
      factor(
        automated_feasibility_tier,
        levels = c(
          "HIGH_PRIORITY_MANUAL_VALIDATION",
          "POTENTIAL_MANUAL_REVIEW",
          "LOW_OR_UNCERTAIN"
        )
      ),
      desc(
        name_score
      ),
      desc(
        source_like_value_fraction
      ),
      desc(
        n_patients_nonmissing
      )
    )

  write_csv_safe(
    inventory,
    file.path(
      DIR_COL,
      "V2_STEP91A_all_metadata_column_inventory.csv"
    )
  )

  write_csv_safe(
    candidates2,
    file.path(
      DIR_REG,
      "V2_STEP91A_infection_source_candidate_registry.csv"
    )
  )

  write_csv_safe(
    values,
    file.path(
      DIR_VAL,
      "V2_STEP91A_candidate_column_value_distributions.csv"
    )
  )

  write_csv_safe(
    consistency,
    file.path(
      DIR_PAT,
      "V2_STEP91A_candidate_patient_level_consistency.csv"
    )
  )

  high <- candidates2 |>
    filter(
      automated_feasibility_tier ==
        "HIGH_PRIORITY_MANUAL_VALIDATION"
    )

  potential <- candidates2 |>
    filter(
      automated_feasibility_tier ==
        "POTENTIAL_MANUAL_REVIEW"
    )

  write_csv_safe(
    high,
    file.path(
      DIR_REG,
      "V2_STEP91A_HIGH_PRIORITY_columns_for_manual_validation.csv"
    )
  )

  write_csv_safe(
    potential,
    file.path(
      DIR_REG,
      "V2_STEP91A_POTENTIAL_columns_for_manual_review.csv"
    )
  )

  project_summary <- candidates2 |>
    group_by(
      project
    ) |>
    summarise(
      n_candidate_columns =
        n(),
      n_high_priority =
        sum(
          automated_feasibility_tier ==
            "HIGH_PRIORITY_MANUAL_VALIDATION"
        ),
      n_potential =
        sum(
          automated_feasibility_tier ==
            "POTENTIAL_MANUAL_REVIEW"
        ),
      best_candidate_column =
        ifelse(
          n() > 0,
          first(column),
          ""
        ),
      best_candidate_tier =
        ifelse(
          n() > 0,
          first(
            automated_feasibility_tier
          ),
          ""
        ),
      .groups = "drop"
    )

  all_projects <- registry |>
    select(
      project,
      analysis_role
    )

  project_summary <- all_projects |>
    left_join(
      project_summary,
      by = "project"
    ) |>
    mutate(
      across(
        c(
          n_candidate_columns,
          n_high_priority,
          n_potential
        ),
        ~replace_na(
          .x,
          0L
        )
      ),
      best_candidate_column =
        replace_na(
          best_candidate_column,
          ""
        ),
      best_candidate_tier =
        replace_na(
          best_candidate_tier,
          "NO_CANDIDATE_COLUMN"
        )
    )

  write_csv_safe(
    project_summary,
    file.path(
      DIR_REG,
      "V2_STEP91A_project_level_source_feasibility_summary.csv"
    )
  )

  readme <- c(
    "SEPSIS V2 - STEP91A INFECTION-SOURCE FEASIBILITY AUDIT",
    paste0(
      "Created: ",
      Sys.time()
    ),
    "",
    "PURPOSE",
    "Determine which frozen cohorts contain defensible metadata for infection-source/site analysis before creating any source categories.",
    "",
    "IMPORTANT",
    "Automated feasibility tiers are screening aids only.",
    "No field is accepted as infection source until its values are manually checked against the source paper/data dictionary.",
    "Phenotype, diagnosis, disease, and clinical labels must not be automatically reinterpreted as infection source.",
    "",
    "OUTPUTS",
    "All metadata columns, candidate column distributions, patient-level consistency, and a project-level feasibility summary are exported.",
    "",
    paste0(
      "High-priority candidate columns: ",
      nrow(high)
    ),
    paste0(
      "Potential manual-review columns: ",
      nrow(potential)
    ),
    "",
    "NEXT",
    "Manually validate the high-priority/potential fields against original study definitions.",
    "Only after that should Step91B freeze a harmonized infection-source variable and define eligible cohorts."
  )

  writeLines(
    readme,
    file.path(
      OUT,
      "README_STEP91A.txt"
    ),
    useBytes = TRUE
  )

  writeLines(
    c(
      paste0(
        "Completed: ",
        Sys.time()
      ),
      "Status: STEP91A COMPLETE",
      "Eight Step87B cohort metadata objects audited.",
      paste0(
        "High-priority source/site candidate columns: ",
        nrow(high)
      ),
      paste0(
        "Potential candidate columns: ",
        nrow(potential)
      ),
      "No infection-source variable was automatically frozen.",
      "Manual source-paper/data-dictionary validation remains mandatory."
    ),
    file.path(
      OUT,
      "_STEP91A_COMPLETE.ok"
    ),
    useBytes = TRUE
  )

  ck(
    "STEP91A COMPLETE"
  )

  cat(
    "\n============================================================\n"
  )
  cat(
    "SEPSIS V2 - STEP91A COMPLETE\n"
  )
  cat(
    "============================================================\n\n"
  )

  cat(
    "Project-level feasibility summary:\n"
  )

  print(
    project_summary,
    n = Inf,
    width = Inf
  )

  cat(
    "\nHigh-priority candidate columns:\n"
  )

  print(
    high,
    n = Inf,
    width = Inf
  )

  cat(
    "\nOutput directory:\n"
  )

  cat(
    OUT,
    "\n"
  )
}

tryCatch(
  main(),
  error = function(e) {

    msg <- c(
      paste0(
        "STEP91A FATAL ERROR: ",
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
          collapse = " "
        )
      )
    )

    writeLines(
      msg,
      ERR,
      useBytes = TRUE
    )

    ck(
      "STEP91A FAILED"
    )

    message(
      paste(
        msg,
        collapse = "\n"
      )
    )

    quit(
      save = "no",
      status = 1,
      runLast = FALSE
    )
  }
)
