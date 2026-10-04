
# ============================================================
# V2_36E RAW / SUPPLEMENTARY CLINICAL METADATA RECOVERY AUDIT
#
# Goal:
# Search raw project folders and supplementary tabular files for
# clinically useful variables that were NOT mapped into the frozen master.
#
# This is a recovery audit only.
# NO microbiome inferential model is run.
# ============================================================

suppressPackageStartupMessages({
  library(readr)
  library(readxl)
  library(dplyr)
  library(stringr)
  library(purrr)
  library(tidyr)
})

ROOT <- "E:/sepsis_project"
DATA <- file.path(ROOT, "data")
OUT <- file.path(ROOT, "results", "V2_36E_RAW_CLINICAL_METADATA_RECOVERY_AUDIT")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

target_projects <- c(
  "PRJNA691455",
  "PRJEB82425",
  "PRJNA516701",
  "PRJNA851469",
  "PRJNA578267",
  "PRJNA430161",
  "PRJNA1166732",
  "PRJNA978257"
)

# More precise clinical regex than previous keyword audits.
clinical_patterns <- list(
  antibiotic = paste0(
    "(^|[^a-z])(",
    "antibiotic|antimicrobial|antibacterial|anti[-_ ]?infective|",
    "ampicillin|amoxicillin|piperacillin|tazobactam|",
    "vancomycin|meropenem|imipenem|ertapenem|",
    "cefazolin|cefotaxime|ceftriaxone|cefepime|ceftazidime|",
    "gentamicin|amikacin|tobramycin|clindamycin|",
    "metronidazole|ciprofloxacin|levofloxacin|linezolid",
    ")($|[^a-z])"
  ),
  sofa = "(^|[^a-z])sofa($|[^a-z])|sequential.?organ.?failure",
  apache = "(^|[^a-z])apache.?ii($|[^a-z])|apache.?2",
  lactate = "(^|[^a-z])lactate($|[^a-z])",
  mortality = "mortality|28.?day|90.?day|death|dead|survival|survivor",
  icu_los = "icu.?los|length.?of.?stay|hospital.?los|icu.?days",
  vasopressor = "vasopressor|norepinephrine|noradrenaline|epinephrine|vasopressin",
  ventilation = "mechanical.?vent|ventilat",
  infection_source = "infection.?source|source.?of.?infection|infection.?site|pneumonia|abdominal|urinary|bloodstream"
)

id_pattern <- "patient|subject|participant|sample|run|accession|biosample|specimen|individual|study.?id|record.?id"

is_informative <- function(x) {
  y <- str_trim(as.character(x))
  !is.na(x) &
    y != "" &
    !tolower(y) %in% c(
      "na","n/a","nan","null","none","unknown",
      "not available","not_available",
      "not reported","not_reported","missing"
    )
}

safe_examples <- function(x, n=5) {
  y <- unique(str_trim(as.character(x[is_informative(x)])))
  paste(head(y, n), collapse = " || ")
}

detect_project <- function(path) {
  hits <- target_projects[str_detect(path, fixed(target_projects))]
  if (length(hits) == 0) return(NA_character_)
  hits[1]
}

read_table_safe <- function(f) {
  ext <- tolower(tools::file_ext(f))
  size_mb <- file.info(f)$size / 1024^2

  # Avoid giant feature/abundance matrices.
  if (is.na(size_mb) || size_mb > 100) return(NULL)

  tryCatch({
    if (ext == "csv") {
      read_csv(f, n_max = 5000, show_col_types = FALSE, progress = FALSE)
    } else if (ext %in% c("tsv","txt")) {
      # Try TSV first; read_delim handles many metadata text files.
      read_delim(f, delim = "\t", n_max = 5000,
                 show_col_types = FALSE, progress = FALSE)
    } else if (ext %in% c("xlsx","xls")) {
      sh <- excel_sheets(f)
      if (length(sh) == 0) return(NULL)
      read_excel(f, sheet = sh[1], n_max = 5000)
    } else {
      NULL
    }
  }, error = function(e) NULL)
}

# ------------------------------------------------------------
# 1. Candidate file discovery
# ------------------------------------------------------------
all_files <- list.files(
  DATA,
  pattern = "\\.(csv|tsv|txt|xlsx|xls)$",
  recursive = TRUE,
  full.names = TRUE,
  ignore.case = TRUE
)

# Prefer files in exact target-project paths.
candidate_files <- all_files[
  map_lgl(all_files, ~ any(str_detect(.x, fixed(target_projects))))
]

# Exclude known downstream feature/abundance-heavy products when obvious.
candidate_files <- candidate_files[
  !str_detect(
    tolower(candidate_files),
    "feature.?table|asv.?table|otu.?table|taxonomy|bray|distance|ordination|rarefied|genus.?matrix"
  )
]

file_inventory <- tibble(
  project = map_chr(candidate_files, detect_project),
  file = candidate_files,
  size_mb = round(file.info(candidate_files)$size / 1024^2, 3)
)

write_csv(file_inventory, file.path(OUT, "01_candidate_tabular_files.csv"))

# ------------------------------------------------------------
# 2. Value-level clinical field recovery
# ------------------------------------------------------------
scan_one <- function(f) {
  d <- read_table_safe(f)
  if (is.null(d) || ncol(d) == 0 || nrow(d) == 0) return(tibble())

  cn <- names(d)
  low <- tolower(cn)

  id_cols <- cn[str_detect(low, id_pattern)]

  rows <- list()

  for (category in names(clinical_patterns)) {
    hit_cols <- cn[str_detect(low, clinical_patterns[[category]])]
    if (length(hit_cols) == 0) next

    for (v in hit_cols) {
      x <- d[[v]]
      inf <- is_informative(x)
      if (sum(inf) == 0) next

      rows[[length(rows)+1]] <- tibble(
        project = detect_project(f),
        file = f,
        category = category,
        variable = v,
        n_rows_read = nrow(d),
        n_informative = sum(inf),
        coverage_pct = round(100 * mean(inf), 2),
        n_unique_nonmissing = n_distinct(as.character(x[inf])),
        has_variation = n_distinct(as.character(x[inf])) >= 2,
        id_columns_present = paste(id_cols, collapse = " | "),
        has_candidate_linkage_id = length(id_cols) > 0,
        example_values = safe_examples(x)
      )
    }
  }

  bind_rows(rows)
}

recovered <- map_dfr(candidate_files, scan_one)

write_csv(
  recovered,
  file.path(OUT, "02_recovered_clinical_variables_value_level.csv")
)

# ------------------------------------------------------------
# 3. Cohort-level recoverability decision
# ------------------------------------------------------------
if (nrow(recovered) > 0) {
  cohort_summary <- recovered %>%
    group_by(project, category) %>%
    summarise(
      candidate_files = n_distinct(file),
      candidate_variables = n_distinct(variable),
      max_informative_rows = max(n_informative, na.rm = TRUE),
      any_variation = any(has_variation),
      any_linkage_id = any(has_candidate_linkage_id),
      .groups = "drop"
    ) %>%
    mutate(
      recovery_priority = case_when(
        any_linkage_id & any_variation & max_informative_rows >= 10 ~ "HIGH",
        any_linkage_id & any_variation & max_informative_rows >= 5 ~ "MODERATE",
        any_linkage_id ~ "LOW",
        TRUE ~ "UNLINKED"
      )
    )
} else {
  cohort_summary <- tibble(
    project = character(),
    category = character(),
    candidate_files = integer(),
    candidate_variables = integer(),
    max_informative_rows = integer(),
    any_variation = logical(),
    any_linkage_id = logical(),
    recovery_priority = character()
  )
}

write_csv(
  cohort_summary,
  file.path(OUT, "03_cohort_level_recovery_summary.csv")
)

# ------------------------------------------------------------
# 4. Shortlist files worth manual mapping
# ------------------------------------------------------------
shortlist <- recovered %>%
  filter(
    has_candidate_linkage_id,
    has_variation,
    n_informative >= 5
  ) %>%
  arrange(
    project,
    category,
    desc(n_informative)
  )

write_csv(
  shortlist,
  file.path(OUT, "04_manual_mapping_shortlist.csv")
)

# ------------------------------------------------------------
# 5. Decision summary
# ------------------------------------------------------------
high <- cohort_summary %>% filter(recovery_priority == "HIGH")
moderate <- cohort_summary %>% filter(recovery_priority == "MODERATE")

lines <- c(
  "V2 STEP36E RAW CLINICAL METADATA RECOVERY AUDIT",
  "",
  paste0("Target frozen cohorts searched: ", length(target_projects)),
  paste0("Candidate tabular files screened: ", length(candidate_files)),
  paste0("Populated clinical variables recovered: ", nrow(recovered)),
  paste0("HIGH-priority recoverable cohort/category pairs: ", nrow(high)),
  paste0("MODERATE-priority recoverable cohort/category pairs: ", nrow(moderate)),
  ""
)

if (nrow(high) + nrow(moderate) == 0) {
  lines <- c(
    lines,
    "DECISION:",
    "No convincing recoverable clinical metadata were identified locally for the frozen main cohorts.",
    "Do NOT launch antibiotic-, SOFA-, APACHE-, lactate-, or outcome-adjusted microbiome models from current local data.",
    "",
    "Recommended pivot:",
    "Preserve the frozen manuscript and consider non-clinical enhancement only (e.g., trajectory/resilience or carefully framed functional inference),",
    "or perform a targeted publication/supplement web retrieval for specific cohorts before any new clinical analysis."
  )
} else {
  lines <- c(
    lines,
    "RECOVERY CANDIDATES EXIST.",
    "Inspect 04_manual_mapping_shortlist.csv.",
    "Only after patient/sample linkage is manually verified should these variables be merged into the frozen analysis layer."
  )
}

writeLines(lines, file.path(OUT, "05_DECISION_SUMMARY.txt"))
writeLines(
  c(paste0("Completed: ", Sys.time()), "STEP36E COMPLETE"),
  file.path(OUT, "_STEP36E_COMPLETE.txt")
)

cat("STEP36E COMPLETE\n")
cat("Output:", OUT, "\n")
