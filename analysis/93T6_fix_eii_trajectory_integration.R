
# ============================================================
# Step93T6 FIX
# Integrate patient trajectory cluster with real EII table
#
# Inputs:
# trajectory_clusters.csv
# V2_STEP92B_patient_timepoint_EII.csv
#
# Goal:
# patient-level trajectory subtype + ecological instability
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
  "V2_33T6_FIX_EII_TRAJECTORY_INTEGRATION"
)

dir.create(
  OUT,
  recursive = TRUE,
  showWarnings = FALSE
)


# -----------------------------
# trajectory cluster
# -----------------------------

cluster_file <- list.files(
  file.path(ROOT,"results"),
  recursive = TRUE,
  full.names = TRUE,
  pattern = "trajectory_clusters\\.csv$"
)[1]

clusters <- read_csv(
  cluster_file,
  show_col_types = FALSE
)

cluster_id_col <- names(clusters)[
  str_detect(
    names(clusters),
    regex("patient", ignore_case=TRUE)
  )
][1]

names(clusters)[names(clusters)==cluster_id_col] <- "Patient_true"

clusters <- clusters %>%
  mutate(
    Patient_true = as.character(Patient_true)
  )


# -----------------------------
# real EII table
# -----------------------------

eii_file <- file.path(
  ROOT,
  "results",
  "V2_32B_STEP92B_EII_CONSTRUCTION",
  "V2_STEP92B_patient_timepoint_EII.csv"
)

if(!file.exists(eii_file)){
  stop("Cannot find Step92B patient timepoint EII file")
}

eii <- read_csv(
  eii_file,
  show_col_types = FALSE
)


# detect patient column
eii_patient_col <- names(eii)[
  str_detect(
    names(eii),
    regex("patient|subject|id",
          ignore_case=TRUE)
  )
][1]

names(eii)[names(eii)==eii_patient_col] <- "Patient_true"


eii <- eii %>%
  mutate(
    Patient_true = as.character(Patient_true)
  )


# -----------------------------
# patient-level aggregation
# -----------------------------

eii_numeric <- names(eii)[
  sapply(
    eii,
    is.numeric
  )
]


if(length(eii_numeric)==0){
  stop("No numeric EII column found")
}

eii_col <- eii_numeric[1]


patient_eii <- eii %>%
  group_by(Patient_true) %>%
  summarise(
    EII_mean = mean(
      .data[[eii_col]],
      na.rm=TRUE
    ),
    EII_max = max(
      .data[[eii_col]],
      na.rm=TRUE
    ),
    n_timepoints=n(),
    .groups="drop"
  )


# -----------------------------
# merge
# -----------------------------

merged <- inner_join(
  clusters,
  patient_eii,
  by="Patient_true"
)


write_csv(
  merged,
  file.path(
    OUT,
    "V2_STEP93T6_patient_trajectory_EII_integrated.csv"
  )
)


cluster_summary <- merged %>%
  group_by(cluster) %>%
  summarise(
    n=n(),
    mean_EII=mean(EII_mean,na.rm=TRUE),
    sd_EII=sd(EII_mean,na.rm=TRUE),
    .groups="drop"
  )

write_csv(
  cluster_summary,
  file.path(
    OUT,
    "V2_STEP93T6_cluster_EII_summary.csv"
  )
)


writeLines(
  c(
    paste0("Completed: ",Sys.time()),
    "STEP93T6 FIX COMPLETE"
  ),
  file.path(
    OUT,
    "_STEP93T6_FIX_COMPLETE.ok"
  )
)

cat("STEP93T6 FIX COMPLETE\n")
