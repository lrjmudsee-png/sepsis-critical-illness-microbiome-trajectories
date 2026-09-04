
# ============================================================
# V2 Metadata Enhancement Audit
# Purpose:
# Identify whether current public cohorts contain variables for:
# 1. Antibiotic exposure
# 2. Clinical severity (SOFA/APACHE)
# 3. Outcomes (mortality/ICU stay)
# 4. Treatment and infection metadata
#
# NO new statistical analysis.
# Only metadata discovery.
# ============================================================

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(stringr)
  library(purrr)
})

ROOT <- "E:/sepsis_project"
DATA <- file.path(ROOT, "data")
OUT <- file.path(ROOT, "results",
                 "V2_36_METADATA_ENHANCEMENT_AUDIT")

dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

# Search all csv/tsv/xlsx exported metadata files
files <- list.files(
  DATA,
  pattern = "\\.(csv|tsv|txt)$",
  recursive = TRUE,
  full.names = TRUE,
  ignore.case = TRUE
)

keywords <- list(
  antibiotic = c(
    "antibiotic","antimicrobial","drug",
    "meropenem","vancomycin","cef",
    "ampicillin","gentamicin"
  ),
  severity = c(
    "sofa","apache","qsofa",
    "severity","score"
  ),
  outcome = c(
    "mortality","death","survival",
    "outcome","icu_days","length",
    "los"
  ),
  treatment = c(
    "vasopressor","ventilation",
    "mechanical","treatment"
  ),
  infection = c(
    "source","origin","site",
    "pneumonia","lung","abdominal"
  )
)

scan_file <- function(f){

  x <- tryCatch(
    readLines(f, warn = FALSE),
    error=function(e) character()
  )

  txt <- tolower(
    paste(x, collapse=" ")
  )

  tibble(
    file=f,
    n_lines=length(x),
    antibiotic=
      any(str_detect(txt, keywords$antibiotic)),
    severity=
      any(str_detect(txt, keywords$severity)),
    outcome=
      any(str_detect(txt, keywords$outcome)),
    treatment=
      any(str_detect(txt, keywords$treatment)),
    infection=
      any(str_detect(txt, keywords$infection))
  )
}

result <- map_dfr(files, scan_file)

write_csv(
  result,
  file.path(
    OUT,
    "metadata_keyword_audit.csv"
  )
)

summary <- result %>%
  summarise(
    files_scanned=n(),
    antibiotic_files=sum(antibiotic),
    severity_files=sum(severity),
    outcome_files=sum(outcome),
    treatment_files=sum(treatment),
    infection_files=sum(infection)
  )

write_csv(
  summary,
  file.path(
    OUT,
    "metadata_audit_summary.csv"
  )
)

# cohort-level keyword mapping
cohort_pattern <- c(
  "PRJNA691455",
  "PRJNA516701",
  "PRJNA851469",
  "PRJEB82425",
  "PRJNA578267",
  "PRJNA430161",
  "PRJNA1166732",
  "PRJNA978257",
  "PRJNA1010969",
  "CRA002354"
)

cohort_map <- map_dfr(
  cohort_pattern,
  function(cn){
    result %>%
      filter(str_detect(file, fixed(cn))) %>%
      mutate(cohort=cn)
  }
)

write_csv(
  cohort_map,
  file.path(
    OUT,
    "cohort_level_metadata_enhancement_map.csv"
  )
)

writeLines(
  c(
    paste0("Completed: ", Sys.time()),
    "Metadata enhancement audit complete.",
    "Review cohort_level_metadata_enhancement_map.csv before adding analyses."
  ),
  file.path(
    OUT,
    "_AUDIT_COMPLETE.txt"
  )
)

cat("V2_36 METADATA AUDIT COMPLETE\n")
