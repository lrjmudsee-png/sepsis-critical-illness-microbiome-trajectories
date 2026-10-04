
# ============================================================
# Sepsis V2 - Step93N2
# Metadata autodiscovery and sample mapping
#
# Automatically searches:
# E:/sepsis_project/data
#
# Supports:
# - SraRunTable.csv
# - runinfo.csv
# - txt
# - tsv
#
# Goal:
# Run -> BioSample -> Sample -> Group -> Patient_ID
# ============================================================

options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(stringr)
  library(purrr)
})

ROOT <- "E:/sepsis_project"

DATA_ROOT <- file.path(
  ROOT,
  "data"
)

OUT <- file.path(
  ROOT,
  "results",
  "V2_33N2_METADATA_AUTODISCOVERY"
)

dir.create(
  OUT,
  recursive = TRUE,
  showWarnings = FALSE
)


# ------------------------------------------------------------
# Find metadata candidates
# ------------------------------------------------------------

all_files <- list.files(
  DATA_ROOT,
  recursive = TRUE,
  full.names = TRUE
)

metadata_files <- all_files[
  str_detect(
    all_files,
    regex(
      "\\.(csv|txt|tsv)$",
      ignore_case = TRUE
    )
  )
]


# ------------------------------------------------------------
# Validate SRA table
# ------------------------------------------------------------

check_run_table <- function(f){

  x <- tryCatch(
    read_csv(
      f,
      show_col_types = FALSE,
      guess_max = 10000
    ),
    error=function(e) NULL
  )

  if(is.null(x))
    return(NULL)

  if("Run" %in% names(x)){
    return(
      tibble(
        file=f,
        rows=nrow(x),
        cols=ncol(x),
        columns=paste(names(x),collapse=";")
      )
    )
  }

  NULL
}


audit <- map_dfr(
  metadata_files,
  check_run_table
)


write_csv(
  audit,
  file.path(
    OUT,
    "V2_STEP93N2_metadata_discovery_audit.csv"
  )
)


if(nrow(audit)==0){

  stop(
    "No valid SRA Run table found."
  )

}


# ------------------------------------------------------------
# Build mapping
# ------------------------------------------------------------

build_mapping <- function(f){

  dat <- read_csv(
    f,
    show_col_types = FALSE,
    guess_max = 10000
  )

  cols <- names(dat)

  get_col <- function(pattern){

    idx <- which(
      str_detect(
        tolower(cols),
        pattern
      )
    )

    if(length(idx)>0)
      cols[idx[1]]
    else
      NA
  }


  sample_col <- get_col(
    "sample_name|sample name|biosample"
  )

  source_col <- get_col(
    "source|title|material|condition|disease"
  )


  result <- tibble(

    Run =
      dat$Run,

    BioSample =
      if("BioSample" %in% cols)
        dat$BioSample
      else
        NA,

    Sample_Name =
      if(!is.na(sample_col))
        dat[[sample_col]]
      else
        NA,

    Source =
      if(!is.na(source_col))
        dat[[source_col]]
      else
        NA,

    metadata_file=f

  )


  result %>%
    mutate(

      Project =
        str_extract(
          f,
          "PRJNA[0-9]+"
        ),

      Group =
        case_when(

          str_detect(
            tolower(
              paste(
                Sample_Name,
                Source
              )
            ),
            "sepsis"
          )
          ~ "Sepsis",

          str_detect(
            tolower(
              paste(
                Sample_Name,
                Source
              )
            ),
            "trauma"
          )
          ~ "Trauma",

          str_detect(
            tolower(
              paste(
                Sample_Name,
                Source
              )
            ),
            "control|healthy|normal"
          )
          ~ "Control",

          TRUE
          ~ "Unknown"
        ),

      Patient_ID =
        str_extract(
          paste(
            Sample_Name,
            Source
          ),
          regex(
            "(sepsis|trauma|control)[-_]?[0-9]+",
            ignore_case = TRUE
          )
        )

    )

}


mapping <- map_dfr(
  audit$file,
  build_mapping
)


write_csv(
  mapping,
  file.path(
    OUT,
    "V2_ALL_PROJECT_SAMPLE_MAPPING.csv"
  )
)


writeLines(
  c(
    paste0(
      "Completed: ",
      Sys.time()
    ),
    "STEP93N2 COMPLETE"
  ),
  file.path(
    OUT,
    "_STEP93N2_COMPLETE.ok"
  )
)


cat(
  "STEP93N2 COMPLETE\n"
)
