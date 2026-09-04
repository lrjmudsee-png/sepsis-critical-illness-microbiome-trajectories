
# ============================================================
# Step93T9
# Run ID -> Patient ID mapping reconstruction
#
# Goal:
# ERR/SRR Run
#       ->
# Sample/BioSample
#       ->
# Patient ID
#
# Used for:
# EII + trajectory integration
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
  "V2_33T9_RUN_TO_PATIENT_MAPPING"
)

dir.create(
  OUT,
  recursive = TRUE,
  showWarnings = FALSE
)


# EII run list
eii_file <- file.path(
  ROOT,
  "results",
  "V2_32B_STEP92B_EII_CONSTRUCTION",
  "V2_STEP92B_patient_timepoint_EII.csv"
)

eii <- read_csv(
  eii_file,
  show_col_types = FALSE
)

run_col <- names(eii)[
  str_detect(
    names(eii),
    regex("ERR|SRR|run",
          ignore_case = TRUE)
  )
][1]

if(is.na(run_col)){
  stop("Cannot detect run column in EII table")
}

eii_runs <- unique(
  as.character(eii[[run_col]])
)


# search metadata
files <- list.files(
  file.path(ROOT,"data"),
  recursive = TRUE,
  full.names = TRUE,
  pattern="\\.(csv|tsv)$",
  ignore.case=TRUE
)


results <- list()

for(f in files){

  x <- tryCatch(
    read_delim(
      f,
      delim=",",
      show_col_types=FALSE,
      guess_max=5000
    ),
    error=function(e) NULL
  )

  if(is.null(x)) next

  cols <- names(x)

  run_exists <- any(
    str_detect(
      cols,
      regex("ERR|SRR|run",
            ignore_case=TRUE)
    )
  )

  patient_exists <- any(
    str_detect(
      cols,
      regex("patient|subject|sample|biosample",
            ignore_case=TRUE)
    )
  )

  if(run_exists && patient_exists){

    results[[length(results)+1]] <- tibble(
      file=f,
      columns=paste(cols,collapse=";"),
      rows=nrow(x)
    )
  }
}


if(length(results)>0){

  bind_rows(results) %>%
    write_csv(
      file.path(
        OUT,
        "RUN_TO_PATIENT_MAPPING_CANDIDATES.csv"
      )
    )

}


writeLines(
  c(
    paste0("Completed: ",Sys.time()),
    "STEP93T9 COMPLETE"
  ),
  file.path(
    OUT,
    "_STEP93T9_COMPLETE.ok"
  )
)

cat("STEP93T9 COMPLETE\n")
