
# ============================================================
# Sepsis V2 - Step93P
# Longitudinal taxonomic trajectory reconstruction
#
# Input:
# PRJNA691455_integrated_genus_profile.csv
#
# Goal:
# Build patient-level genus trajectories.
#
# Outputs:
# - patient trajectory matrix
# - baseline vs follow-up changes
# ============================================================

options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(stringr)
})

ROOT <- "E:/sepsis_project"

IN <- file.path(
  ROOT,
  "results",
  "V2_33O_TAXONOMIC_ABUNDANCE_METADATA_INTEGRATION",
  "PRJNA691455_integrated_genus_profile.csv"
)

OUT <- file.path(
  ROOT,
  "results",
  "V2_33P_LONGITUDINAL_TAXONOMIC_TRAJECTORY_ANALYSIS"
)

dir.create(
  OUT,
  recursive=TRUE,
  showWarnings=FALSE
)

if(!file.exists(IN)){
  stop("Missing integrated genus profile.")
}

dat <- read_csv(
  IN,
  show_col_types=FALSE
)


# identify genus columns

meta_cols <- c(
  "Project",
  "Run",
  "BioSample",
  "Sample_Name",
  "Source",
  "Group",
  "Patient_ID",
  "Patient_ID_reconstructed",
  "Timepoint_reconstructed"
)

genus_cols <- setdiff(
  names(dat),
  meta_cols
)


# choose reconstructed patient ID

dat <- dat %>%
  mutate(
    Patient_final =
      coalesce(
        Patient_ID_reconstructed,
        Patient_ID,
        Sample_Name
      )
  )


# clean timepoint

dat <- dat %>%
  mutate(
    Time_numeric =
      as.numeric(
        str_extract(
          Timepoint_reconstructed,
          "[0-9]+"
        )
      )
  )


write_csv(
  dat,
  file.path(
    OUT,
    "V2_STEP93P_clean_longitudinal_genus_profile.csv"
  )
)


# only patients with repeated measurements

trajectory_patients <- dat %>%
  group_by(Patient_final) %>%
  summarise(
    n_timepoints =
      n_distinct(Time_numeric),
    .groups="drop"
  ) %>%
  filter(
    n_timepoints >= 2
  )


write_csv(
  trajectory_patients,
  file.path(
    OUT,
    "V2_STEP93P_repeated_measurement_patients.csv"
  )
)


writeLines(
  c(
    paste0("Completed: ",Sys.time()),
    "STEP93P COMPLETE"
  ),
  file.path(
    OUT,
    "_STEP93P_COMPLETE.ok"
  )
)

cat("STEP93P COMPLETE\n")
