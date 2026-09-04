
# ============================================================
# Sepsis V2 - Step93N
# Build Run_ID -> sample -> group metadata mapping
#
# Supports:
# PRJNA1010969
# PRJNA691455
#
# Input:
# SraRunTable.csv / *_runinfo.csv
#
# Output:
# Unified sample mapping tables
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
  "V2_33N_STEP93N_SAMPLE_MAPPING"
)

dir.create(
  OUT,
  recursive = TRUE,
  showWarnings = FALSE
)

projects <- c(
  "PRJNA1010969",
  "PRJNA691455"
)

find_run_table <- function(project){

  candidates <- list.files(
    file.path(ROOT,"data",project),
    pattern = "SraRunTable|runinfo|\\.csv$",
    recursive = TRUE,
    full.names = TRUE,
    ignore.case = TRUE
  )

  # exclude html files accidentally downloaded
  candidates <- candidates[
    !str_detect(
      candidates,
      "html"
    )
  ]

  candidates
}


build_mapping <- function(project){

  files <- find_run_table(project)

  if(length(files)==0){
    warning(
      paste("No run table found:",project)
    )
    return(NULL)
  }

  # choose first valid csv with Run column
  dat <- NULL

  for(f in files){

    tmp <- tryCatch(
      read_csv(
        f,
        show_col_types = FALSE
      ),
      error=function(e) NULL
    )

    if(!is.null(tmp) && "Run" %in% names(tmp)){
      dat <- tmp
      break
    }
  }

  if(is.null(dat)){
    stop(
      paste(
        "Cannot identify valid SRA table:",
        project
      )
    )
  }


  names_lower <- tolower(names(dat))


  pick_col <- function(patterns){

    idx <- which(
      Reduce(
        "|",
        lapply(
          patterns,
          function(x)
            str_detect(
              names_lower,
              x
            )
        )
      )
    )

    if(length(idx)>0)
      names(dat)[idx[1]]
    else
      NA
  }


  sample_col <- pick_col(
    c(
      "sample_name",
      "sample name",
      "biosample"
    )
  )

  source_col <- pick_col(
    c(
      "source_material",
      "source",
      "title",
      "sample_title"
    )
  )


  result <- tibble(
    project = project,
    Run = dat$Run,

    BioSample =
      if("BioSample" %in% names(dat))
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
        NA
  )


  result <- result %>%
    mutate(

      Group = case_when(

        str_detect(
          tolower(
            paste(
              Sample_Name,
              Source
            )
          ),
          "sepsis"
        ) ~ "Sepsis",

        str_detect(
          tolower(
            paste(
              Sample_Name,
              Source
            )
          ),
          "trauma"
        ) ~ "Trauma",

        str_detect(
          tolower(
            paste(
              Sample_Name,
              Source
            ),
          ),
          "control|healthy|hc"
        ) ~ "Control",

        TRUE ~ "Unknown"
      ),

      Patient_ID =
        case_when(

          str_detect(
            paste(
              Sample_Name,
              Source
            ),
            regex(
              "sepsis\\d+|trauma\\d+|control\\d+",
              ignore_case=TRUE
            )
          ) ~ str_extract(
            paste(
              Sample_Name,
              Source
            ),
            regex(
              "(sepsis|trauma|control)\\d+",
              ignore_case=TRUE
            )
          ),

          TRUE ~ Sample_Name
        )
    )


  write_csv(
    result,
    file.path(
      OUT,
      paste0(
        project,
        "_sample_mapping.csv"
      )
    )
  )

  result
}


all_map <- map_dfr(
  projects,
  build_mapping
)


write_csv(
  all_map,
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
    "STEP93N COMPLETE"
  ),
  file.path(
    OUT,
    "_STEP93N_COMPLETE.ok"
  )
)

cat("STEP93N COMPLETE\n")
