
# ============================================================
# Step93T FIX
# Fix Patient_true type mismatch in EII integration
#
# Error:
# clusters$Patient_true character
# eii$Patient_true double
#
# Solution:
# force both identifiers to character before join
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
  "V2_33T_EII_TRAJECTORY_INTEGRATION"
)

dir.create(
  OUT,
  recursive=TRUE,
  showWarnings=FALSE
)


cluster_file <- list.files(
  file.path(ROOT,"results"),
  recursive=TRUE,
  full.names=TRUE,
  pattern="trajectory_clusters\\.csv$"
)[1]


eii_candidates <- list.files(
  file.path(ROOT,"results"),
  recursive=TRUE,
  full.names=TRUE,
  pattern="(EII|eii|instability|instability_index).*\\.csv$"
)

if(length(eii_candidates)==0){
  stop("No EII file found")
}

eii_file <- eii_candidates[1]


clusters <- read_csv(
  cluster_file,
  show_col_types=FALSE
)

eii <- read_csv(
  eii_file,
  show_col_types=FALSE
)


# identify patient columns

cluster_patient <- names(clusters)[
  str_detect(
    names(clusters),
    regex("patient",ignore_case=TRUE)
  )
][1]

names(clusters)[
  names(clusters)==cluster_patient
] <- "Patient_true"


eii_patient <- names(eii)[
  str_detect(
    names(eii),
    regex("patient|sample|id",ignore_case=TRUE)
  )
][1]

names(eii)[
  names(eii)==eii_patient
] <- "Patient_true"


# critical fix
clusters <- clusters %>%
  mutate(
    Patient_true = as.character(Patient_true)
  )

eii <- eii %>%
  mutate(
    Patient_true = as.character(Patient_true)
  )


merged <- inner_join(
  clusters,
  eii,
  by="Patient_true"
)


write_csv(
  merged,
  file.path(
    OUT,
    "V2_STEP93T_EII_TRAJECTORY_CLUSTER_INTEGRATED_FIXED.csv"
  )
)


writeLines(
  c(
    paste0("Completed: ",Sys.time()),
    "STEP93T FIX COMPLETE"
  ),
  file.path(
    OUT,
    "_STEP93T_FIX_COMPLETE.ok"
  )
)

cat("STEP93T FIX COMPLETE\n")
