# ============================================================
# Sepsis V2 - Step93K
# Genus abundance matrix validation
#
# Purpose:
# Validate true longitudinal genus relative abundance matrices.
#
# Step93J found:
# PRJNA691455_genus_relative_abundance.csv
# PRJNA1010969_genus_relative_abundance.csv
#
# This step inspects structure:
# - metadata columns
# - genus feature columns
# - patient/time availability
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
  "V2_33K_STEP93K_GENUS_MATRIX_VALIDATION"
)

dir.create(
  OUT,
  recursive=TRUE,
  showWarnings=FALSE
)

candidate <- c(
  "E:/sepsis_project/data/PRJNA691455/03_dada2_rerun_v2/PRJNA691455_genus_relative_abundance.csv",
  "E:/sepsis_project/data/PRJNA1010969/03_dada2_rerun_v2/PRJNA1010969_genus_relative_abundance.csv"
)

audit <- lapply(
  candidate,
  function(f){

    if(!file.exists(f)){
      return(
        tibble(
          file=f,
          exists=FALSE
        )
      )
    }

    x <- read_csv(
      f,
      show_col_types=FALSE,
      n_max=5
    )

    tibble(
      file=f,
      exists=TRUE,
      rows_preview=nrow(x),
      cols=ncol(x),
      columns=paste(names(x),collapse=";"),
      has_patient=
        any(str_detect(
          tolower(names(x)),
          "patient|subject|sample|id"
        )),
      has_time=
        any(str_detect(
          tolower(names(x)),
          "day|time|visit"
        ))
    )
  }
) %>%
  bind_rows()

write_csv(
  audit,
  file.path(
    OUT,
    "V2_STEP93K_genus_matrix_validation.csv"
  )
)

writeLines(
  c(
    paste0("Completed: ",Sys.time()),
    "STEP93K COMPLETE"
  ),
  file.path(
    OUT,
    "_STEP93K_COMPLETE.ok"
  )
)

cat("STEP93K COMPLETE\n")
