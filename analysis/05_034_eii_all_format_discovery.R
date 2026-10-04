
# ============================================================
# Step93T5
# Search EII patient table in ALL common formats
#
# Searches:
# csv txt tsv xlsx rds rdata
#
# Goal:
# Find patient-level EII table
# ============================================================

options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(stringr)
  library(purrr)
})

ROOT <- "E:/sepsis_project"

OUT <- file.path(
  ROOT,
  "results",
  "V2_33T5_EII_ALL_FORMAT_DISCOVERY"
)

dir.create(
  OUT,
  recursive=TRUE,
  showWarnings=FALSE
)


files <- list.files(
  file.path(ROOT,"results"),
  recursive=TRUE,
  full.names=TRUE
)


# -----------------------------
# CSV/TSV scan
# -----------------------------

scan_table <- function(x, file){

  if(is.null(x))
    return(NULL)

  cols <- names(x)

  patient <- any(
    str_detect(
      cols,
      regex("patient|sample|subject|id",
            ignore_case=TRUE)
    )
  )

  eii <- cols[
    str_detect(
      cols,
      regex("eii|instability|ecological|index",
            ignore_case=TRUE)
    )
  ]

  if(patient && length(eii)>0){

    tibble(
      file=file,
      rows=nrow(x),
      columns=paste(cols,collapse=";"),
      eii_columns=paste(eii,collapse=";")
    )

  } else NULL

}


csv_candidates <- files[
  str_detect(
    files,
    regex("\\.(csv|txt|tsv)$",
          ignore_case=TRUE)
  )
]


csv_audit <- map_dfr(
  csv_candidates,
  function(f){

    x <- tryCatch(
      read_delim(
        f,
        delim=",",
        show_col_types=FALSE,
        guess_max=5000
      ),
      error=function(e) NULL
    )

    scan_table(x,f)

  }
)


# -----------------------------
# RDS scan
# -----------------------------

rds_candidates <- files[
  str_detect(
    files,
    regex("\\.rds$",
          ignore_case=TRUE)
  )
]


rds_audit <- map_dfr(
  rds_candidates,
  function(f){

    x <- tryCatch(
      readRDS(f),
      error=function(e) NULL
    )

    if(is.data.frame(x))
      scan_table(x,f)
    else NULL

  }
)


audit <- bind_rows(
  csv_audit,
  rds_audit
)


write_csv(
  audit,
  file.path(
    OUT,
    "EII_all_format_candidate_audit.csv"
  )
)


if(nrow(audit)==0){

  writeLines(
    "No patient-level EII table detected.",
    file.path(
      OUT,
      "NO_EII_TABLE_FOUND.txt"
    )
  )

  stop(
    "No patient-level EII table detected."
  )

}


writeLines(
  c(
    paste0("Completed: ",Sys.time()),
    "STEP93T5 COMPLETE"
  ),
  file.path(
    OUT,
    "_STEP93T5_COMPLETE.ok"
  )
)

cat("STEP93T5 COMPLETE\n")
