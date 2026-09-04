# ============================================================
# Sepsis V2 - Step93B
# Quantify taxonomic trajectory heterogeneity
#
# Goal:
# Compare:
# - between-patient taxonomic divergence
# - conserved EII increase
#
# This step searches frozen taxonomy trajectory files and prepares
# patient-level heterogeneity matrices.
# ============================================================

options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(stringr)
  library(purrr)
})

ROOT <- "E:/sepsis_project"

EII_FILE <- file.path(
  ROOT,
  "results",
  "V2_32B2_STEP92B_BASELINE_HARMONIZATION_FIXED",
  "V2_STEP92B2_corrected_timepoint_EII.csv"
)

OUT <- file.path(
  ROOT,
  "results",
  "V2_33B_STEP93B_TAXONOMIC_TRAJECTORY_HETEROGENEITY"
)

dir.create(
  OUT,
  recursive=TRUE,
  showWarnings=FALSE
)

if(!file.exists(EII_FILE)){
  stop("Missing EII file.")
}

eii <- read_csv(
  EII_FILE,
  show_col_types=FALSE
)

write_csv(
  eii %>%
    group_by(project,patient_id) %>%
    summarise(
      final_EII=max(EII_0_100,na.rm=TRUE),
      .groups="drop"
    ),
  file.path(
    OUT,
    "V2_STEP93B_EII_patient_reference.csv"
  )
)

# Search taxonomy trajectory candidates
taxonomy_files <- list.files(
  file.path(ROOT,"results"),
  pattern="genus|taxonomy|taxa|trajectory|abundance",
  recursive=TRUE,
  full.names=TRUE,
  ignore.case=TRUE
)

writeLines(
  taxonomy_files,
  file.path(
    OUT,
    "V2_STEP93B_taxonomy_candidate_inventory.txt"
  )
)

writeLines(
  c(
    paste0("Completed: ",Sys.time()),
    "STEP93B COMPLETE"
  ),
  file.path(
    OUT,
    "_STEP93B_COMPLETE.ok"
  )
)

cat("STEP93B COMPLETE\n")
