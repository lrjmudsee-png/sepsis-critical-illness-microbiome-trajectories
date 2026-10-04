# ============================================================
# Sepsis V2 - Step93A
# Taxonomic heterogeneity versus ecological instability
#
# Goal:
# Demonstrate:
# 1. patients have heterogeneous taxonomic trajectories
# 2. ecological instability is conserved
#
# Inputs:
# Step92B2 EII
# Step90B taxonomy trajectory outputs (if available)
# ============================================================

options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(stringr)
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
  "V2_33A_STEP93_TAXONOMIC_HETEROGENEITY_VS_EII"
)

dir.create(
  OUT,
  recursive=TRUE,
  showWarnings=FALSE
)

if(!file.exists(EII_FILE)){
  stop("Missing EII input.")
}

eii <- read_csv(
  EII_FILE,
  show_col_types=FALSE
)

# ------------------------------------------------------------
# EII conservation summary
# ------------------------------------------------------------

eii_summary <- eii %>%
  group_by(
    project,
    patient_id
  ) %>%
  summarise(
    max_EII=max(
      EII_0_100,
      na.rm=TRUE
    ),
    EII_increase=
      max(EII_0_100,na.rm=TRUE) -
      min(EII_0_100,na.rm=TRUE),
    .groups="drop"
  )

write_csv(
  eii_summary,
  file.path(
    OUT,
    "V2_STEP93A_EII_patient_conservation_summary.csv"
  )
)

# ------------------------------------------------------------
# Optional taxonomy inventory
# ------------------------------------------------------------

taxonomy_candidates <- list.files(
  file.path(
    ROOT,
    "results"
  ),
  pattern="taxonomy|genus|family|feature",
  recursive=TRUE,
  full.names=TRUE,
  ignore.case=TRUE
)

writeLines(
  taxonomy_candidates,
  file.path(
    OUT,
    "V2_STEP93A_taxonomy_file_inventory.txt"
  )
)

writeLines(
  c(
    paste0("Completed: ",Sys.time()),
    "STEP93A COMPLETE",
    "Taxonomic heterogeneity framework prepared."
  ),
  file.path(
    OUT,
    "_STEP93A_COMPLETE.ok"
  )
)

cat("STEP93A COMPLETE\n")
