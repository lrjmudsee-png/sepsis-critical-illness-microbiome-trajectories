# ============================================================
# Sepsis V2 - Step 75
# Audit and auto-fetch metadata for projects present in data/
# but absent from the Step73/74 full master
# R 4.4.0 / Windows
#
# This is a TRIAGE step:
# - fetch public ENA run metadata when possible
# - inspect local metadata / sequence availability
# - detect patient/time/outcome/infection fields
# - assign scientific priority
#
# It does NOT add a project into the master automatically.
# ============================================================

options(stringsAsFactors = FALSE)

pkgs <- c(
  "readr","readxl","dplyr","tidyr","stringr",
  "purrr","tibble","httr2","xml2"
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
  library(readxl)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(purrr)
  library(tibble)
  library(httr2)
  library(xml2)
})

# ------------------------------------------------------------
# Paths
# ------------------------------------------------------------
PROJECT_ROOT <- "E:/sepsis_project"
DATA_ROOT <- file.path(PROJECT_ROOT, "data")

OUT_ROOT <- file.path(
  PROJECT_ROOT,
  "results",
  "V2_15_missing_project_triage"
)

dir.create(
  OUT_ROOT,
  recursive = TRUE,
  showWarnings = FALSE
)

# ------------------------------------------------------------
# Projects found by Step74 but absent from full master
# ------------------------------------------------------------
TARGETS <- c(
  "CRA002354",
  "PRJEB37289",
  "PRJEB68229",
  "PRJNA1125274",
  "PRJNA1143057",
  "PRJNA480143",
  "PRJNA533528",
  "PRJNA641414",
  "PRJNA851469",
  "PRJNA951907"
)

# ------------------------------------------------------------
# Scientific triage based on project design already reviewed
# ------------------------------------------------------------
scientific_plan <- tibble(
  project = TARGETS,

  proposed_role = c(
    "CORE_SEPSIS_LONGITUDINAL_CANDIDATE",
    "STATIC_ICU_SEPSIS_BACKGROUND",
    "ICU_OUTCOME_LONGITUDINAL_CANDIDATE",
    "CORE_SEPSIS_LONGITUDINAL_CANDIDATE",
    "SEPSIS_AKI_OUTCOME_SUPPORT",
    "UNRESOLVED_REVIEW_REQUIRED",
    "SHOTGUN_ICU_BACKGROUND_LONGITUDINAL",
    "STATIC_INFECTION_RISK_SUPPORT",
    "ICU_NOSOCOMIAL_INFECTION_LONGITUDINAL_CANDIDATE",
    "UNRESOLVED_REVIEW_REQUIRED"
  ),

  priority = c(
    "HIGH",
    "LOW",
    "HIGH",
    "VERY_HIGH",
    "MEDIUM",
    "REVIEW",
    "MEDIUM",
    "LOW",
    "HIGH",
    "REVIEW"
  ),

  auto_include = FALSE,

  rationale = c(
    "ICU sepsis/septic shock cohort; recover explicit patient/time mapping before inclusion.",
    "Previously adjudicated as static ICU/sepsis background rather than core longitudinal.",
    "Serial ICU gut microbiome with mortality/outcome relevance; high-value longitudinal validation candidate.",
    "SURVEIL longitudinal rectal swabs at baseline, sepsis onset and discharge; highly aligned with V2 primary question.",
    "Sepsis-associated AKI microbiome study; useful organ-dysfunction/outcome support if sample timing can be resolved.",
    "Do not include until local/public metadata establishes study design and patient/time structure.",
    "Serial ICU faecal shotgun cohort; keep as later shotgun/background validation.",
    "Rectal Klebsiella colonization/infection-risk study; not a core sepsis longitudinal cohort.",
    "Prospective critical-illness rectal-swab longitudinal cohort with nosocomial infection outcomes; public metadata may be partly restricted.",
    "Do not include until local/public metadata establishes study design and patient/time structure."
  )
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
      "na","nan","n/a","null","none"
    )
  ] <- NA_character_
  x
}

safe_download_text <- function(url) {
  tryCatch(
    {
      resp <- request(url) |>
        req_user_agent("SepsisV2-Step75-R44/1.0") |>
        req_timeout(120) |>
        req_perform()

      if (resp_status(resp) >= 300) return(NULL)

      resp_body_string(resp)
    },
    error = function(e) NULL
  )
}

safe_csv <- function(path) {
  tryCatch(
    suppressMessages(
      read_csv(
        path,
        show_col_types = FALSE,
        progress = FALSE,
        n_max = 5000,
        name_repair = "unique"
      )
    ),
    error = function(e) NULL
  )
}

safe_tsv <- function(path) {
  tryCatch(
    suppressMessages(
      read_tsv(
        path,
        show_col_types = FALSE,
        progress = FALSE,
        n_max = 5000,
        name_repair = "unique"
      )
    ),
    error = function(e) NULL
  )
}

safe_excel <- function(path) {
  sheets <- tryCatch(
    excel_sheets(path),
    error = function(e) character()
  )

  out <- list()

  for (sh in sheets) {
    x <- tryCatch(
      suppressMessages(
        read_excel(
          path,
          sheet = sh,
          n_max = 5000,
          .name_repair = "unique"
        )
      ),
      error = function(e) NULL
    )

    if (!is.null(x)) {
      out[[length(out) + 1]] <- list(
        sheet = sh,
        data = x
      )
    }
  }

  out
}

field_patterns <- list(
  patient_id = c(
    "patient","subject","participant",
    "host.?subject","individual.?id"
  ),

  sample_id = c(
    "sample.?id","biosample","sample.?alias",
    "specimen","sequence.?name"
  ),

  run_id = c(
    "run_accession","\\bsrr","\\berr","\\bcrr"
  ),

  time = c(
    "timepoint","time.?point","\\bday\\b",
    "collection.?date","visit","baseline",
    "discharge","admission","onset"
  ),

  sepsis = c(
    "sepsis","septic","shock"
  ),

  infection = c(
    "infection","bacteremia","bacteraemia",
    "bloodstream","nosocomial","pneumonia"
  ),

  outcome = c(
    "mortality","death","survival","outcome",
    "icu.?stay","hospital.?stay"
  ),

  antibiotics = c(
    "antibiotic","antimicrobial","abx",
    "meropenem","carbapenem"
  ),

  severity = c(
    "sofa","apache","severity","vasopressor",
    "organ.?failure"
  ),

  body_site = c(
    "stool","fecal","faecal","rectal",
    "gut","tracheal","body.?site","sample.?type"
  )
)

detect_fields <- function(text) {
  text <- tolower(
    paste(
      text,
      collapse = " "
    )
  )

  vapply(
    field_patterns,
    function(pats) {
      any(
        vapply(
          pats,
          function(p) {
            str_detect(
              text,
              regex(
                p,
                ignore_case = TRUE
              )
            )
          },
          logical(1)
        )
      )
    },
    logical(1)
  )
}

profile_df <- function(df, project, source, sheet = "") {
  if (is.null(df) || ncol(df) == 0) return(tibble())

  det <- detect_fields(
    c(
      names(df),
      unlist(
        lapply(
          df[
            seq_len(
              min(
                ncol(df),
                80
              )
            )
          ],
          function(x) {
            head(
              unique(
                clean_chr(x)
              ),
              15
            )
          }
        )
      )
    )
  )

  tibble(
    project = project,
    source = source,
    sheet = sheet,
    rows = nrow(df),
    columns = ncol(df),

    patient_id = det["patient_id"],
    sample_id = det["sample_id"],
    run_id = det["run_id"],
    time = det["time"],
    sepsis = det["sepsis"],
    infection = det["infection"],
    outcome = det["outcome"],
    antibiotics = det["antibiotics"],
    severity = det["severity"],
    body_site = det["body_site"],

    column_names = paste(
      names(df),
      collapse = ";"
    )
  )
}

# ------------------------------------------------------------
# ENA fetch
# ------------------------------------------------------------
ena_fields <- paste(
  c(
    "run_accession",
    "sample_accession",
    "secondary_sample_accession",
    "sample_alias",
    "sample_title",
    "sample_description",
    "study_accession",
    "secondary_study_accession",
    "experiment_accession",
    "experiment_alias",
    "experiment_title",
    "scientific_name",
    "collection_date",
    "host_sex",
    "host_body_site",
    "library_layout",
    "library_strategy",
    "fastq_ftp",
    "fastq_md5",
    "fastq_bytes"
  ),
  collapse = ","
)

fetch_ena_run <- function(project) {
  url <- paste0(
    "https://www.ebi.ac.uk/ena/portal/api/filereport?",
    "accession=",
    URLencode(project, reserved = TRUE),
    "&result=read_run",
    "&fields=",
    URLencode(ena_fields, reserved = TRUE),
    "&format=tsv",
    "&download=true"
  )

  txt <- safe_download_text(url)

  if (
    is.null(txt) ||
    nchar(txt) < 20 ||
    !str_detect(txt, "run_accession")
  ) {
    return(
      list(
        data = NULL,
        url = url,
        status = "ENA_NOT_AVAILABLE_OR_QUERY_FAILED"
      )
    )
  }

  tf <- tempfile(
    fileext = ".tsv"
  )

  writeLines(
    txt,
    tf,
    useBytes = TRUE
  )

  df <- safe_tsv(tf)

  unlink(tf)

  list(
    data = df,
    url = url,
    status = if (
      is.null(df)
    ) {
      "ENA_PARSE_FAILED"
    } else {
      "ENA_OK"
    }
  )
}

# ------------------------------------------------------------
# Local inventory
# ------------------------------------------------------------
inventory_rows <- list()
profile_rows <- list()
ena_rows <- list()

for (project in TARGETS) {

  message(
    "Auditing ",
    project
  )

  pdir <- file.path(
    DATA_ROOT,
    project
  )

  if (!dir.exists(pdir)) {
    next
  }

  files <- list.files(
    pdir,
    recursive = TRUE,
    full.names = TRUE
  )

  files <- files[
    file.exists(files) &
    !dir.exists(files)
  ]

  ext <- tolower(
    tools::file_ext(files)
  )

  inventory_rows[[project]] <- tibble(
    project = project,

    local_files = length(files),

    fastq_gz = sum(
      str_detect(
        basename(files),
        regex(
          "\\.fastq\\.gz$",
          ignore_case = TRUE
        )
      )
    ),

    fa_gz = sum(
      str_detect(
        basename(files),
        regex(
          "\\.fa\\.gz$",
          ignore_case = TRUE
        )
      )
    ),

    aria2 = sum(
      str_detect(
        basename(files),
        regex(
          "\\.aria2$",
          ignore_case = TRUE
        )
      )
    ),

    csv = sum(
      ext == "csv"
    ),

    tsv = sum(
      ext %in% c(
        "tsv","txt"
      )
    ),

    xlsx = sum(
      ext %in% c(
        "xlsx","xls"
      )
    ),

    rds_rdata = sum(
      ext %in% c(
        "rds","rdata","rda"
      )
    ),

    biom = sum(
      ext == "biom"
    ),

    total_size_GB = round(
      sum(
        file.info(files)$size,
        na.rm = TRUE
      ) / 1024^3,
      3
    )
  )

  # ----------------------------------------------------------
  # Local structured metadata profile
  # ----------------------------------------------------------
  meta_files <- files[
    ext %in% c(
      "csv","tsv","txt","xlsx","xls"
    )
  ]

  # Skip huge taxonomic/feature tables when possible.
  meta_files <- meta_files[
    !str_detect(
      basename(meta_files),
      regex(
        "asv|otu|taxonomy|taxa|feature.?table|abundance",
        ignore_case = TRUE
      )
    )
  ]

  for (path in head(meta_files, 30)) {

    e <- tolower(
      tools::file_ext(path)
    )

    if (e == "csv") {
      df <- safe_csv(path)

      if (!is.null(df)) {
        profile_rows[[length(profile_rows) + 1]] <-
          profile_df(
            df,
            project,
            basename(path)
          )
      }
    }

    if (e %in% c("tsv","txt")) {
      df <- safe_tsv(path)

      if (!is.null(df)) {
        profile_rows[[length(profile_rows) + 1]] <-
          profile_df(
            df,
            project,
            basename(path)
          )
      }
    }

    if (e %in% c("xlsx","xls")) {

      sheets <- safe_excel(
        path
      )

      for (item in sheets) {
        profile_rows[[length(profile_rows) + 1]] <-
          profile_df(
            item$data,
            project,
            basename(path),
            item$sheet
          )
      }
    }
  }

  # ----------------------------------------------------------
  # Public ENA metadata
  # ----------------------------------------------------------
  if (
    str_detect(
      project,
      "^PRJ"
    )
  ) {

    ena <- fetch_ena_run(
      project
    )

    save_dir <- file.path(
      pdir,
      "00_metadata",
      "step75_auto_review"
    )

    dir.create(
      save_dir,
      recursive = TRUE,
      showWarnings = FALSE
    )

    if (!is.null(ena$data)) {

      ena_path <- file.path(
        save_dir,
        paste0(
          project,
          "_ENA_read_run_step75.tsv"
        )
      )

      write_tsv(
        ena$data,
        ena_path,
        na = ""
      )

      profile_rows[[length(profile_rows) + 1]] <-
        profile_df(
          ena$data,
          project,
          basename(ena_path)
        )
    }

    ena_rows[[project]] <- tibble(
      project = project,
      ena_status = ena$status,
      ena_rows = if (
        is.null(
          ena$data
        )
      ) {
        0
      } else {
        nrow(
          ena$data
        )
      },
      ena_url = ena$url
    )
  }
}

local_inventory <- bind_rows(
  inventory_rows
)

field_profile <- if (
  length(profile_rows) > 0
) {
  bind_rows(
    profile_rows
  )
} else {
  tibble()
}

ena_status <- if (
  length(ena_rows) > 0
) {
  bind_rows(
    ena_rows
  )
} else {
  tibble(
    project = character(),
    ena_status = character(),
    ena_rows = integer(),
    ena_url = character()
  )
}

# ------------------------------------------------------------
# Aggregate metadata evidence
# ------------------------------------------------------------
field_summary <- if (
  nrow(field_profile) > 0
) {

  field_profile |>
    group_by(project) |>
    summarise(
      patient_id_detected =
        any(patient_id),

      sample_id_detected =
        any(sample_id),

      run_id_detected =
        any(run_id),

      time_detected =
        any(time),

      sepsis_detected =
        any(sepsis),

      infection_detected =
        any(infection),

      outcome_detected =
        any(outcome),

      antibiotics_detected =
        any(antibiotics),

      severity_detected =
        any(severity),

      body_site_detected =
        any(body_site),

      structured_sources =
        n(),

      .groups = "drop"
    )

} else {

  tibble(
    project = TARGETS
  )
}

# ------------------------------------------------------------
# Final triage
# ------------------------------------------------------------
triage <- scientific_plan |>
  left_join(
    local_inventory,
    by = "project"
  ) |>
  left_join(
    ena_status |>
      select(
        project,
        ena_status,
        ena_rows
      ),
    by = "project"
  ) |>
  left_join(
    field_summary,
    by = "project"
  ) |>
  mutate(

    local_sequence_present =
      coalesce(
        fastq_gz,
        0L
      ) +
      coalesce(
        fa_gz,
        0L
      ) > 0,

    download_incomplete =
      coalesce(
        aria2,
        0L
      ) > 0,

    core_mapping_fields_detected =
      coalesce(
        patient_id_detected,
        FALSE
      ) &
      coalesce(
        sample_id_detected,
        FALSE
      ) &
      coalesce(
        time_detected,
        FALSE
      ),

    recommended_next_action = case_when(

      project == "PRJNA1125274" &
        core_mapping_fields_detected ~
        "BUILD_PATIENT_TIME_MAP_AND_ADD_TO_CORE",

      project == "PRJNA1125274" ~
        "PARSE_PUBLICATION_SUPPLEMENT_FOR_PATIENT_TIME",

      project == "CRA002354" ~
        "PARSE_GSA_SAMPLE_NAMES_WITH_PUBLICATION_MAPPING",

      project == "PRJEB68229" &
        core_mapping_fields_detected ~
        "BUILD_PATIENT_TIME_OUTCOME_MAP",

      project == "PRJEB68229" ~
        "RECOVER_PATIENT_TIME_OUTCOME_FROM_PUBLICATION_METADATA",

      project == "PRJNA851469" &
        core_mapping_fields_detected ~
        "BUILD_PUBLIC_LONGITUDINAL_MAP",

      project == "PRJNA851469" ~
        "PARSE_SUPPLEMENTARY_TABLES; DO_NOT WAIT FOR RESTRICTED_METADATA",

      project == "PRJNA1143057" ~
        "ASSESS_AKI_OUTCOME_MAPPING_FOR_SUPPORT_MODULE",

      project == "PRJNA533528" ~
        "KEEP_SHOTGUN_BACKGROUND; DOWNLOAD_CAN_CONTINUE_IN_PARALLEL",

      project == "PRJEB37289" ~
        "KEEP_STATIC_BACKGROUND_ONLY",

      project == "PRJNA641414" ~
        "LOW_PRIORITY_STATIC_INFECTION_RISK_SUPPORT",

      priority == "REVIEW" ~
        "REVIEW_LOCAL_AND_ENA_METADATA_BEFORE_ANY_INCLUSION",

      TRUE ~
        "NO_IMMEDIATE_ACTION"
    )
  )

# ------------------------------------------------------------
# Save
# ------------------------------------------------------------
write_excel_csv(
  scientific_plan,
  file.path(
    OUT_ROOT,
    "V2_step75_scientific_priority_plan.csv"
  ),
  na = ""
)

write_excel_csv(
  local_inventory,
  file.path(
    OUT_ROOT,
    "V2_step75_local_project_inventory.csv"
  ),
  na = ""
)

write_excel_csv(
  field_profile,
  file.path(
    OUT_ROOT,
    "V2_step75_metadata_file_profiles.csv"
  ),
  na = ""
)

write_excel_csv(
  ena_status,
  file.path(
    OUT_ROOT,
    "V2_step75_ENA_fetch_status.csv"
  ),
  na = ""
)

write_excel_csv(
  triage,
  file.path(
    OUT_ROOT,
    "V2_step75_missing_project_triage.csv"
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
  "SEPSIS V2 - STEP 75 COMPLETE\n"
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

print(
  triage |>
    select(
      project,
      proposed_role,
      priority,
      local_sequence_present,
      download_incomplete,
      ena_rows,
      patient_id_detected,
      time_detected,
      outcome_detected,
      recommended_next_action
    ),
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
