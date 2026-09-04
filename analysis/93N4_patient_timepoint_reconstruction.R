
# ============================================================
# Sepsis V2 - Step93N4
# Patient and timepoint reconstruction
#
# Input:
# V2_33N3_METADATA_STRICT_MAPPING outputs
#
# Goal:
# Run -> Patient_ID -> Timepoint -> Group
#
# Handles:
# PRJNA691455
# PRJNA1010969
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
  "V2_33N3_METADATA_STRICT_MAPPING",
  "V2_ALL_PROJECT_SAMPLE_MAPPING.csv"
)

OUT <- file.path(
  ROOT,
  "results",
  "V2_33N4_PATIENT_TIMEPOINT_RECONSTRUCTION"
)

dir.create(
  OUT,
  recursive=TRUE,
  showWarnings=FALSE
)

if(!file.exists(IN)){
  stop("Missing Step93N3 mapping file.")
}

dat <- read_csv(
  IN,
  show_col_types=FALSE
)


# ------------------------------------------------------------
# PRJNA691455 reconstruction
# ------------------------------------------------------------

dat <- dat %>%
  mutate(

    Patient_ID_reconstructed =
      case_when(

        Project=="PRJNA691455" ~

          str_extract(
            Sample_Name,
            "^[A-Za-z0-9]+"
          ),

        TRUE ~ NA_character_
      ),


    Timepoint_reconstructed =
      case_when(

        Project=="PRJNA691455" &
          str_detect(
            tolower(Source),
            "day"
          ) ~

          str_extract(
            tolower(Source),
            "day[ ]?[0-9]+"
          ),


        TRUE ~ NA_character_
      )
  )


# ------------------------------------------------------------
# PRJNA1010969 reconstruction
#
# Expected coding:
# HMB = healthy/control
# TMB = trauma
# SMB = sepsis
# ------------------------------------------------------------

dat <- dat %>%
  mutate(

    Patient_ID_reconstructed =
      case_when(

        Project=="PRJNA1010969" &
          str_detect(
            Sample_Name,
            "^[A-Za-z]+[0-9]+"
          ) ~

          str_extract(
            Sample_Name,
            "^[A-Za-z]+[0-9]+"
          ),

        TRUE ~ Patient_ID_reconstructed
      ),


    Group_reconstructed =
      case_when(

        Project=="PRJNA1010969" &
          str_detect(
            Sample_Name,
            "^SMB"
          ) ~ "Sepsis",

        Project=="PRJNA1010969" &
          str_detect(
            Sample_Name,
            "^TMB"
          ) ~ "Trauma",

        Project=="PRJNA1010969" &
          str_detect(
            Sample_Name,
            "^HMB"
          ) ~ "Control",

        TRUE ~ Group
      )

  )


write_csv(
  dat,
  file.path(
    OUT,
    "V2_ALL_PROJECT_PATIENT_TIMEPOINT_MAPPING.csv"
  )
)


write_csv(
  dat %>% filter(Project=="PRJNA691455"),
  file.path(
    OUT,
    "PRJNA691455_patient_timepoint_mapping.csv"
  )
)


write_csv(
  dat %>% filter(Project=="PRJNA1010969"),
  file.path(
    OUT,
    "PRJNA1010969_patient_timepoint_mapping.csv"
  )
)


writeLines(
  c(
    paste0("Completed: ",Sys.time()),
    "STEP93N4 COMPLETE"
  ),
  file.path(
    OUT,
    "_STEP93N4_COMPLETE.ok"
  )
)

cat("STEP93N4 COMPLETE\n")
