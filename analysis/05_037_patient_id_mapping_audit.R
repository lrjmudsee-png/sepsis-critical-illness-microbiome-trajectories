
# ============================================================
# Step93T7
# Patient ID mapping audit
#
# Goal:
# Compare trajectory patient IDs and EII patient IDs.
# Generate candidate mappings.
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
  "V2_33T7_PATIENT_ID_MAPPING_AUDIT"
)

dir.create(
  OUT,
  recursive=TRUE,
  showWarnings=FALSE
)


# trajectory cluster
cluster_file <- list.files(
  file.path(ROOT,"results"),
  recursive=TRUE,
  full.names=TRUE,
  pattern="trajectory_clusters\\.csv$"
)[1]

clusters <- read_csv(
  cluster_file,
  show_col_types=FALSE
)

cluster_id <- names(clusters)[
  str_detect(
    names(clusters),
    regex("patient", ignore_case=TRUE)
  )
][1]

names(clusters)[names(clusters)==cluster_id] <- "trajectory_patient"

clusters <- clusters %>%
  mutate(
    trajectory_patient = as.character(trajectory_patient)
  ) %>%
  distinct(trajectory_patient)


# EII table
eii_file <- file.path(
  ROOT,
  "results",
  "V2_32B_STEP92B_EII_CONSTRUCTION",
  "V2_STEP92B_patient_timepoint_EII.csv"
)

eii <- read_csv(
  eii_file,
  show_col_types=FALSE
)

eii_id <- names(eii)[
  str_detect(
    names(eii),
    regex("patient|sample|subject|id",
          ignore_case=TRUE)
  )
][1]

names(eii)[names(eii)==eii_id] <- "eii_patient"

eii <- eii %>%
  mutate(
    eii_patient=as.character(eii_patient)
  ) %>%
  distinct(eii_patient)


# exact match
exact <- clusters %>%
  inner_join(
    eii,
    by=c(
      "trajectory_patient"="eii_patient"
    )
  ) %>%
  mutate(
    match_type="exact"
  )


# substring candidates
candidate <- expand.grid(
  trajectory_patient=clusters$trajectory_patient,
  eii_patient=eii$eii_patient,
  stringsAsFactors=FALSE
) %>%
  mutate(
    match_type=case_when(
      str_detect(
        eii_patient,
        fixed(trajectory_patient)
      ) ~ "trajectory_in_eii",

      str_detect(
        trajectory_patient,
        fixed(eii_patient)
      ) ~ "eii_in_trajectory",

      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(match_type))


write_csv(
  exact,
  file.path(
    OUT,
    "exact_patient_ID_match.csv"
  )
)

write_csv(
  candidate,
  file.path(
    OUT,
    "patient_ID_mapping_candidates.csv"
  )
)

write_csv(
  clusters,
  file.path(
    OUT,
    "trajectory_patient_ID_list.csv"
  )
)

write_csv(
  eii,
  file.path(
    OUT,
    "EII_patient_ID_list.csv"
  )
)


writeLines(
  c(
    paste0("Completed: ",Sys.time()),
    "STEP93T7 COMPLETE"
  ),
  file.path(
    OUT,
    "_STEP93T7_COMPLETE.ok"
  )
)

cat("STEP93T7 COMPLETE\n")
