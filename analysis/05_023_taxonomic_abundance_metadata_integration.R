
# ============================================================
# Sepsis V2 - Step93O
# Taxonomic abundance + metadata integration
#
# Input:
# genus_relative_abundance.csv
# patient/timepoint mapping
#
# Projects:
# PRJNA691455
# PRJNA1010969
#
# Output:
# integrated genus profiles
# matching audit
# ============================================================

options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(stringr)
})

ROOT <- "E:/sepsis_project"

OUT <- file.path(
  ROOT,
  "results",
  "V2_33O_TAXONOMIC_ABUNDANCE_METADATA_INTEGRATION"
)

dir.create(
  OUT,
  recursive=TRUE,
  showWarnings=FALSE
)


projects <- c(
  "PRJNA691455",
  "PRJNA1010969"
)


integrate_project <- function(project){

  abundance_candidates <- list.files(
    file.path(ROOT,"data",project),
    pattern="genus_relative_abundance\\.csv$",
    recursive=TRUE,
    full.names=TRUE,
    ignore.case=TRUE
  )

  if(length(abundance_candidates)==0){
    stop(
      paste(
        "No genus abundance found:",
        project
      )
    )
  }

  abundance_file <- abundance_candidates[1]


  abundance <- read_csv(
    abundance_file,
    show_col_types=FALSE
  )


  mapping_file <- file.path(
    ROOT,
    "results",
    "V2_33N4_PATIENT_TIMEPOINT_RECONSTRUCTION",
    paste0(
      project,
      "_patient_timepoint_mapping.csv"
    )
  )


  if(!file.exists(mapping_file)){
    stop(
      paste(
        "Missing mapping:",
        project
      )
    )
  }


  mapping <- read_csv(
    mapping_file,
    show_col_types=FALSE
  )


  # standardize Run column

  if(!"Run" %in% names(abundance)){

    run_col <- names(abundance)[
      str_detect(
        names(abundance),
        regex(
          "run|sample|id",
          ignore_case=TRUE
        )
      )
    ][1]

    names(abundance)[
      names(abundance)==run_col
    ] <- "Run"
  }


  merged <- mapping %>%
    inner_join(
      abundance,
      by="Run"
    )


  audit <- tibble(
    project=project,
    abundance_rows=nrow(abundance),
    mapping_rows=nrow(mapping),
    matched_rows=nrow(merged),
    match_rate=
      nrow(merged)/nrow(mapping)
  )


  write_csv(
    audit,
    file.path(
      OUT,
      paste0(
        project,
        "_integration_audit.csv"
      )
    )
  )


  write_csv(
    merged,
    file.path(
      OUT,
      paste0(
        project,
        "_integrated_genus_profile.csv"
      )
    )
  )


  merged
}


all <- lapply(
  projects,
  integrate_project
)


write_csv(
  bind_rows(all),
  file.path(
    OUT,
    "V2_ALL_INTEGRATED_GENUS_PROFILE.csv"
  )
)


writeLines(
  c(
    paste0(
      "Completed: ",
      Sys.time()
    ),
    "STEP93O COMPLETE"
  ),
  file.path(
    OUT,
    "_STEP93O_COMPLETE.ok"
  )
)

cat("STEP93O COMPLETE\n")
