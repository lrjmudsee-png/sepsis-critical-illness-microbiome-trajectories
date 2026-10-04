
# ============================================================
# Sepsis V2 - Step93Q
# Fix patient trajectory reconstruction
#
# Problem from Step93P:
# Patient_final was sample-level because PRJNA691455
# sample names contain patient and day information.
#
# Goal:
# Extract true patient ID:
# NY35H13006FFO021 -> NY35H13006
# NY35H13006FFO011 -> NY35H13006
#
# Then calculate repeated longitudinal patients.
# ============================================================

options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(stringr)
})

ROOT <- "E:/sepsis_project"

IN <- file.path(
  ROOT,
  "results",
  "V2_33P_LONGITUDINAL_TAXONOMIC_TRAJECTORY_ANALYSIS",
  "V2_STEP93P_clean_longitudinal_genus_profile.csv"
)

OUT <- file.path(
  ROOT,
  "results",
  "V2_33Q_TRAJECTORY_PATIENT_ID_RECONSTRUCTION_FIX"
)

dir.create(
  OUT,
  recursive=TRUE,
  showWarnings=FALSE
)

dat <- read_csv(
  IN,
  show_col_types=FALSE
)

dat <- dat %>%
  mutate(

    Patient_true =
      case_when(

        Project=="PRJNA691455" ~

          str_extract(
            Sample_Name,
            "^[A-Za-z0-9]+(?=FFO)"
          ),

        TRUE ~ Patient_final
      )

  )


patient_summary <- dat %>%
  group_by(Patient_true) %>%
  summarise(
    n_samples=n(),
    n_timepoints=n_distinct(Time_numeric),
    timepoints=paste(
      sort(unique(Time_numeric)),
      collapse=","
    ),
    .groups="drop"
  )


write_csv(
  dat,
  file.path(
    OUT,
    "V2_STEP93Q_fixed_longitudinal_genus_profile.csv"
  )
)

write_csv(
  patient_summary,
  file.path(
    OUT,
    "V2_STEP93Q_patient_timepoint_summary.csv"
  )
)

write_csv(
  patient_summary %>%
    filter(n_timepoints>=2),
  file.path(
    OUT,
    "V2_STEP93Q_repeated_measurement_patients.csv"
  )
)

writeLines(
  c(
    paste0("Completed: ",Sys.time()),
    "STEP93Q COMPLETE"
  ),
  file.path(
    OUT,
    "_STEP93Q_COMPLETE.ok"
  )
)

cat("STEP93Q COMPLETE\n")
