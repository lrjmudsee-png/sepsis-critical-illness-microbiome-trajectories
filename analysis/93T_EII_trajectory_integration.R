
# ============================================================
# Sepsis V2 - Step93T
# EII and taxonomic trajectory integration
#
# Goal:
# Combine:
# 1. trajectory clusters
# 2. signed trajectory distance
# 3. EII patient-level index
#
# Since EII filename may differ between versions,
# script searches V2 results recursively.
# ============================================================

options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(stringr)
  library(ggplot2)
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


# ----------------------------
# Locate trajectory cluster
# ----------------------------

cluster_file <- list.files(
  file.path(ROOT,"results"),
  recursive=TRUE,
  full.names=TRUE,
  pattern="trajectory_clusters\\.csv$"
)[1]


if(is.na(cluster_file)){
  stop("trajectory_clusters.csv not found")
}


clusters <- read_csv(
  cluster_file,
  show_col_types=FALSE
)


# ----------------------------
# Locate EII file
# ----------------------------

eii_candidates <- list.files(
  file.path(ROOT,"results"),
  recursive=TRUE,
  full.names=TRUE,
  pattern="(EII|eii|instability|instability_index).*\\.csv$"
)


if(length(eii_candidates)==0){

  stop(
    "No EII csv found. Please check previous EII output location."
  )

}


eii_file <- eii_candidates[1]


eii <- read_csv(
  eii_file,
  show_col_types=FALSE
)


# ----------------------------
# Detect patient column
# ----------------------------

patient_col <- names(clusters)[
  str_detect(
    names(clusters),
    regex("patient",ignore_case=TRUE)
  )
][1]


names(clusters)[
  names(clusters)==patient_col
] <- "Patient_true"


eii_patient_col <- names(eii)[
  str_detect(
    names(eii),
    regex("patient|sample|id",ignore_case=TRUE)
  )
][1]


names(eii)[
  names(eii)==eii_patient_col
] <- "Patient_true"


merged <- inner_join(
  clusters,
  eii,
  by="Patient_true"
)


write_csv(
  merged,
  file.path(
    OUT,
    "V2_STEP93T_EII_TRAJECTORY_CLUSTER_INTEGRATED.csv"
  )
)


# cluster EII summary

numeric_cols <- names(merged)[
  sapply(
    merged,
    is.numeric
  )
]


if(length(numeric_cols)>0){

  eii_col <- numeric_cols[1]

  summary <- merged %>%
    group_by(cluster) %>%
    summarise(
      n=n(),
      mean_EII=mean(
        .data[[eii_col]],
        na.rm=TRUE
      ),
      sd_EII=sd(
        .data[[eii_col]],
        na.rm=TRUE
      ),
      .groups="drop"
    )

  write_csv(
    summary,
    file.path(
      OUT,
      "V2_STEP93T_cluster_EII_summary.csv"
    )
  )

}


writeLines(
  c(
    paste0("Completed: ",Sys.time()),
    "STEP93T COMPLETE"
  ),
  file.path(
    OUT,
    "_STEP93T_COMPLETE.ok"
  )
)

cat("STEP93T COMPLETE\n")
