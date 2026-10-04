# ============================================================
# Sepsis V2 - Step93C
# Taxonomic divergence versus EII conservation
#
# Purpose:
# Quantify whether patients show heterogeneous taxonomic routes
# despite conserved ecological instability.
#
# This step searches available abundance/trajectory objects and
# prepares diversity-of-route analysis.
# ============================================================

options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(stringr)
  library(tidyr)
})

ROOT <- "E:/sepsis_project"

EII_FILE <- file.path(
  ROOT,
  "results",
  "V2_33B_STEP93B_TAXONOMIC_TRAJECTORY_HETEROGENEITY",
  "V2_STEP93B_EII_patient_reference.csv"
)

OUT <- file.path(
  ROOT,
  "results",
  "V2_33C_STEP93C_TAXONOMIC_DIVERGENCE_VS_EII"
)

dir.create(
  OUT,
  recursive=TRUE,
  showWarnings=FALSE
)

if(!file.exists(EII_FILE)){
  stop("Missing EII reference.")
}

eii <- read_csv(
  EII_FILE,
  show_col_types=FALSE
)

write_csv(
  eii,
  file.path(
    OUT,
    "V2_STEP93C_EII_reference_for_taxonomy_comparison.csv"
  )
)

# Search candidate abundance matrices generated in previous steps

candidate_files <- list.files(
  file.path(ROOT,"results"),
  pattern="abundance|relative|genus|ASV|taxa|taxonomy",
  recursive=TRUE,
  full.names=TRUE,
  ignore.case=TRUE
)

candidate_table <- tibble(
  taxonomy_candidate_file=candidate_files
)

write_csv(
  candidate_table,
  file.path(
    OUT,
    "V2_STEP93C_taxonomy_candidate_files.csv"
  )
)

writeLines(
  c(
    paste0(
      "Completed: ",
      Sys.time()
    ),
    "STEP93C COMPLETE",
    "Taxonomic divergence framework prepared."
  ),
  file.path(
    OUT,
    "_STEP93C_COMPLETE.ok"
  )
)

cat("STEP93C COMPLETE\n")
