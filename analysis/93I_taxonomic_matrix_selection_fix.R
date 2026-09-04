# ============================================================
# Sepsis V2 - Step93I
# Taxonomic matrix selection fix
#
# Problem:
# Step93H selected a paired differential genus table, not an
# abundance trajectory matrix.
#
# This step searches previous outputs for true abundance matrices:
# sample x genus / patient x genus.
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
  "V2_33I_STEP93I_TAXONOMIC_MATRIX_SELECTION_FIX"
)

dir.create(
  OUT,
  recursive=TRUE,
  showWarnings=FALSE
)

files <- list.files(
  file.path(ROOT,"results"),
  pattern="genus|relative|abundance|profile|clr",
  recursive=TRUE,
  full.names=TRUE,
  ignore.case=TRUE
)

audit <- lapply(
  files,
  function(f){

    x <- tryCatch(
      read_csv(
        f,
        show_col_types=FALSE,
        n_max=5
      ),
      error=function(e) NULL
    )

    if(is.null(x)){
      return(NULL)
    }

    cols <- names(x)

    tibble(
      file=f,
      n_preview_cols=length(cols),
      has_patient=
        any(str_detect(tolower(cols),"patient|subject|sample|id")),
      has_time=
        any(str_detect(tolower(cols),"day|time|visit|early|late")),
      has_taxon=
        any(str_detect(tolower(cols),"taxon|genus|species")),
      has_abundance=
        any(str_detect(tolower(cols),"abundance|relative|ra|clr"))
    )
  }
)

audit <- bind_rows(audit)

write_csv(
  audit,
  file.path(
    OUT,
    "V2_STEP93I_taxonomy_matrix_search_audit.csv"
  )
)

candidates <- audit %>%
  filter(
    has_patient,
    has_taxon,
    has_abundance
  )

write_csv(
  candidates,
  file.path(
    OUT,
    "V2_STEP93I_true_trajectory_matrix_candidates.csv"
  )
)

writeLines(
  c(
    paste0("Completed: ",Sys.time()),
    "STEP93I COMPLETE",
    "True abundance matrix selection completed."
  ),
  file.path(
    OUT,
    "_STEP93I_COMPLETE.ok"
  )
)

cat("STEP93I COMPLETE\n")
