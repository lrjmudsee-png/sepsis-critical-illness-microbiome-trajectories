
# ============================================================
# Step93T3
# Identify real patient-level EII table
#
# Previous issue:
# metadata readiness/audit tables were selected as EII.
#
# Criteria:
# 1. Must contain patient/sample identifier
# 2. Must contain numeric EII-like column
# 3. Prefer small patient-level tables
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
  "V2_33T3_EII_REAL_TABLE_IDENTIFICATION"
)

dir.create(
  OUT,
  recursive=TRUE,
  showWarnings=FALSE
)


# trajectory cluster
cluster_file <- list.files(
  file.path(ROOT,"results"),
  recursive=TRUE,
  full.names=TRUE,
  pattern="trajectory_clusters\\.csv$"
)[1]

clusters <- read_csv(
  cluster_file,
  show_col_types=FALSE
)

patient_col <- names(clusters)[
  str_detect(names(clusters), regex("patient",ignore_case=TRUE))
][1]

names(clusters)[names(clusters)==patient_col] <- "Patient_true"

clusters <- clusters %>%
  mutate(Patient_true=as.character(Patient_true))


# all csv candidates
files <- list.files(
  file.path(ROOT,"results"),
  recursive=TRUE,
  full.names=TRUE,
  pattern="\\.csv$"
)


scan_one <- function(f){

  x <- tryCatch(
    read_csv(f, show_col_types=FALSE),
    error=function(e) NULL
  )

  if(is.null(x))
    return(NULL)

  cols <- names(x)

  has_patient <- any(
    str_detect(
      cols,
      regex("patient|sample|id",
            ignore_case=TRUE)
    )
  )

  eii_cols <- cols[
    str_detect(
      cols,
      regex("EII|instability|index",
            ignore_case=TRUE)
    )
  ]

  if(!has_patient || length(eii_cols)==0)
    return(NULL)


  numeric_eii <- eii_cols[
    sapply(
      x[eii_cols],
      is.numeric
    )
  ]

  tibble(
    file=f,
    nrow=nrow(x),
    columns=paste(cols,collapse=";"),
    eii_numeric_columns=paste(
      numeric_eii,
      collapse=";"
    )
  )
}


audit <- map_dfr(
  files,
  scan_one
)


write_csv(
  audit,
  file.path(
    OUT,
    "real_EII_candidate_audit.csv"
  )
)


if(nrow(audit)==0)
  stop("No real EII table found")


# select smallest patient-level candidate
candidate <- audit %>%
  filter(
    eii_numeric_columns!=""
  ) %>%
  arrange(nrow)


eii <- read_csv(
  candidate$file[1],
  show_col_types=FALSE
)

patient_col <- names(eii)[
  str_detect(
    names(eii),
    regex("patient|sample|id",
          ignore_case=TRUE)
  )
][1]

names(eii)[names(eii)==patient_col] <- "Patient_true"

eii <- eii %>%
  mutate(
    Patient_true=as.character(Patient_true)
  )


merged <- inner_join(
  clusters,
  eii,
  by="Patient_true"
)


write_csv(
  merged,
  file.path(
    OUT,
    "V2_STEP93T3_EII_TRAJECTORY_INTEGRATED.csv"
  )
)

writeLines(
  c(
    paste0("Completed: ",Sys.time()),
    "STEP93T3 COMPLETE"
  ),
  file.path(
    OUT,
    "_STEP93T3_COMPLETE.ok"
  )
)

cat("STEP93T3 COMPLETE\n")
