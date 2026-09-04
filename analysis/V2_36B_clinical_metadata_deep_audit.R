
# ============================================================
# V2_36B Clinical Metadata Deep Audit
# Purpose:
# Confirm whether real clinical metadata support:
# 1. antibiotic exposure analysis
# 2. severity score analysis
# 3. outcome analysis
#
# This script does NOT run new statistics.
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
                 "V2_36B_CLINICAL_METADATA_DEEP_AUDIT")

dir.create(OUT, recursive=TRUE, showWarnings=FALSE)

files <- list.files(
  DATA,
  pattern="\\.(csv|tsv|txt)$",
  recursive=TRUE,
  full.names=TRUE,
  ignore.case=TRUE
)

patterns <- list(
  antibiotic=c(
    "antibiotic","antimicrobial","ampicillin",
    "vancomycin","meropenem","cef",
    "gentamicin","drug","medication"
  ),
  severity=c(
    "sofa","apache","qsofa",
    "severity","score"
  ),
  outcome=c(
    "mortality","death","survival",
    "icu","length.of.stay","los"
  ),
  infection_source=c(
    "source","origin","site",
    "pneumonia","lung","abdominal"
  )
)

scan_columns <- function(f){

  dat <- tryCatch(
    read_csv(f, n_max=5, show_col_types=FALSE),
    error=function(e) NULL
  )

  if(is.null(dat)){
    return(tibble())
  }

  cols <- names(dat)

  tibble(
    file=f,
    n_columns=length(cols),
    columns=paste(cols, collapse=" | "),
    antibiotic=any(str_detect(
      tolower(cols),
      paste(patterns$antibiotic,collapse="|")
    )),
    severity=any(str_detect(
      tolower(cols),
      paste(patterns$severity,collapse="|")
    )),
    outcome=any(str_detect(
      tolower(cols),
      paste(patterns$outcome,collapse="|")
    )),
    infection_source=any(str_detect(
      tolower(cols),
      paste(patterns$infection_source,collapse="|")
    ))
  )
}

inventory <- map_dfr(files, scan_columns)

write_csv(
  inventory,
  file.path(OUT,"clinical_variable_inventory.csv")
)

feasibility <- tibble(
  analysis=c(
    "Antibiotic exposure adjustment",
    "Severity association (SOFA/APACHE)",
    "Outcome association (mortality/LOS)",
    "Infection source stratification"
  ),
  required_variables=c(
    "antibiotic / antimicrobial / medication variables",
    "SOFA/APACHE/qSOFA/severity variables",
    "mortality/death/LOS variables",
    "infection source variables"
  ),
  files_found=c(
    sum(inventory$antibiotic,na.rm=TRUE),
    sum(inventory$severity,na.rm=TRUE),
    sum(inventory$outcome,na.rm=TRUE),
    sum(inventory$infection_source,na.rm=TRUE)
  )
)

write_csv(
  feasibility,
  file.path(OUT,"enhancement_analysis_feasibility.csv")
)

writeLines(
  c(
    paste0("Completed: ",Sys.time()),
    "",
    "Review:",
    "clinical_variable_inventory.csv",
    "enhancement_analysis_feasibility.csv",
    "",
    "No statistical analysis performed."
  ),
  file.path(OUT,"_AUDIT_COMPLETE.txt")
)

cat("V2_36B COMPLETE\n")
