
# ============================================================
# Step93T4
# Search EII from Excel/workbook sources
#
# Previous:
# CSV search found only audit tables.
#
# This step searches xlsx files and sheets.
# Criteria:
# Patient identifier + EII-like numeric column
# ============================================================

options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(readxl)
  library(readr)
  library(dplyr)
  library(stringr)
  library(purrr)
})

ROOT <- "E:/sepsis_project"

OUT <- file.path(
  ROOT,
  "results",
  "V2_33T4_EII_XLSX_PATIENT_TABLE_SEARCH"
)

dir.create(
  OUT,
  recursive=TRUE,
  showWarnings=FALSE
)

files <- list.files(
  file.path(ROOT,"results"),
  recursive=TRUE,
  full.names=TRUE,
  pattern="\\.xlsx?$",
  ignore.case=TRUE
)


scan_sheet <- function(f,s){

  x <- tryCatch(
    read_excel(f, sheet=s),
    error=function(e) NULL
  )

  if(is.null(x))
    return(NULL)

  cols <- names(x)

  patient <- any(
    str_detect(
      cols,
      regex("patient|sample|id",
            ignore_case=TRUE)
    )
  )

  eii <- cols[
    str_detect(
      cols,
      regex("EII|instability|index|ecological",
            ignore_case=TRUE)
    )
  ]

  if(patient && length(eii)>0){

    tibble(
      file=f,
      sheet=s,
      nrow=nrow(x),
      columns=paste(cols,collapse=";"),
      eii_columns=paste(eii,collapse=";")
    )

  } else {
    NULL
  }
}


audit <- map_dfr(
  files,
  function(f){

    sheets <- excel_sheets(f)

    map_dfr(
      sheets,
      ~scan_sheet(f,.x)
    )

  }
)


write_csv(
  audit,
  file.path(
    OUT,
    "EII_xlsx_candidate_audit.csv"
  )
)


if(nrow(audit)==0)
  stop("No Excel EII table found")


# choose candidate with smallest patient-level table
candidate <- audit %>%
  arrange(nrow)


writeLines(
  c(
    paste0("Candidate: ",candidate$file[1]),
    paste0("Sheet: ",candidate$sheet[1]),
    "STEP93T4 COMPLETE"
  ),
  file.path(
    OUT,
    "_STEP93T4_COMPLETE.ok"
  )
)

cat("STEP93T4 COMPLETE\n")
