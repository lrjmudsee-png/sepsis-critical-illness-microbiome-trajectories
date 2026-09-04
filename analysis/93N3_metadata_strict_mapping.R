
# ============================================================
# Sepsis V2 - Step93N3
# Strict metadata mapping for selected projects only
#
# Projects:
# PRJNA1010969
# PRJNA691455
#
# Avoid contamination from audit/QC files.
#
# Required metadata:
# Run + BioSample + BioProject
# OR Run + BioSample + Sample_Name
# ============================================================

options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(stringr)
  library(purrr)
})

ROOT <- "E:/sepsis_project"

DATA_ROOT <- file.path(ROOT,"data")

OUT <- file.path(
  ROOT,
  "results",
  "V2_33N3_METADATA_STRICT_MAPPING"
)

dir.create(
  OUT,
  recursive=TRUE,
  showWarnings=FALSE
)

projects <- c(
  "PRJNA1010969",
  "PRJNA691455"
)


# ------------------------------------------------------------
# Find only project-specific metadata
# ------------------------------------------------------------

find_files <- function(project){

  project_dir <- file.path(
    DATA_ROOT,
    project
  )

  files <- list.files(
    project_dir,
    recursive=TRUE,
    full.names=TRUE
  )

  files <- files[
    str_detect(
      files,
      regex(
        "\\.(csv|txt|tsv)$",
        ignore_case=TRUE
      )
    )
  ]

  # exclude analysis/QC files
  files <- files[
    !str_detect(
      basename(files),
      regex(
        "audit|repair|before|after|summary|qc|depth|problem|locked",
        ignore_case=TRUE
      )
    )
  ]

  files
}


# ------------------------------------------------------------
# Identify true SRA tables
# ------------------------------------------------------------

validate_sra <- function(f){

  x <- tryCatch(
    read_csv(
      f,
      show_col_types=FALSE,
      guess_max=10000
    ),
    error=function(e) NULL
  )

  if(is.null(x))
    return(NULL)

  cols <- names(x)

  valid <-
    "Run" %in% cols &&
    (
      "BioSample" %in% cols ||
      any(
        str_detect(
          cols,
          regex(
            "sample",
            ignore_case=TRUE
          )
        )
      )
    )

  if(valid){

    tibble(
      file=f,
      rows=nrow(x),
      cols=ncol(x),
      columns=paste(cols,collapse=";")
    )

  } else {
    NULL
  }

}


audit <- map_dfr(
  projects,
  function(project){

    files <- find_files(project)

    map_dfr(
      files,
      validate_sra
    ) %>%
      mutate(
        Project=project
      )
  }
)


write_csv(
  audit,
  file.path(
    OUT,
    "V2_STEP93N3_valid_metadata_candidates.csv"
  )
)


if(nrow(audit)==0){

  stop(
    "No valid SRA metadata found. Check file placement."
  )

}


# ------------------------------------------------------------
# Build final mapping
# ------------------------------------------------------------

build_mapping <- function(f,project){

  dat <- read_csv(
    f,
    show_col_types=FALSE,
    guess_max=10000
  )

  cols <- names(dat)


  sample_col <- cols[
    str_detect(
      cols,
      regex(
        "Sample_Name|SampleName|sample name|sample_name",
        ignore_case=TRUE
      )
    )
  ][1]


  source_col <- cols[
    str_detect(
      cols,
      regex(
        "source|title|material|condition",
        ignore_case=TRUE
      )
    )
  ][1]


  result <- tibble(
    Project=project,
    Run=dat$Run,

    BioSample=
      if("BioSample"%in%cols)
        dat$BioSample
      else NA,

    Sample_Name=
      if(!is.na(sample_col))
        dat[[sample_col]]
      else NA,

    Source=
      if(!is.na(source_col))
        dat[[source_col]]
      else NA
  )


  result %>%
    mutate(

      Group=
        case_when(

          str_detect(
            tolower(
              paste(Sample_Name,Source)
            ),
            "sepsis"
          )
          ~ "Sepsis",

          str_detect(
            tolower(
              paste(Sample_Name,Source)
            ),
            "trauma"
          )
          ~ "Trauma",

          str_detect(
            tolower(
              paste(Sample_Name,Source)
            ),
            "control|healthy|normal"
          )
          ~ "Control",

          TRUE
          ~ "Unknown"
        ),

      Patient_ID=
        str_extract(
          paste(Sample_Name,Source),
          regex(
            "(sepsis|trauma|control)[-_]?[0-9]+",
            ignore_case=TRUE
          )
        )

    )

}


mapping <- map_dfr(
  1:nrow(audit),
  function(i){

    build_mapping(
      audit$file[i],
      audit$Project[i]
    )
  }
)


write_csv(
  mapping,
  file.path(
    OUT,
    "V2_ALL_PROJECT_SAMPLE_MAPPING.csv"
  )
)


mapping %>%
  filter(
    Project=="PRJNA1010969"
  ) %>%
  write_csv(
    file.path(
      OUT,
      "PRJNA1010969_sample_mapping.csv"
    )
  )


mapping %>%
  filter(
    Project=="PRJNA691455"
  ) %>%
  write_csv(
    file.path(
      OUT,
      "PRJNA691455_sample_mapping.csv"
    )
  )


writeLines(
  c(
    paste0("Completed: ",Sys.time()),
    "STEP93N3 COMPLETE"
  ),
  file.path(
    OUT,
    "_STEP93N3_COMPLETE.ok"
  )
)

cat("STEP93N3 COMPLETE\n")
