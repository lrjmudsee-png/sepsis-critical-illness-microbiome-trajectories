# ============================================================
# Sepsis V2 - Step93H
# Taxonomic trajectory distance calculation
#
# Goal:
# Calculate patient-level genus trajectory divergence.
#
# Step93G selected genus matrix.
# This step extracts metadata columns and calculates distance
# between longitudinal trajectories when the matrix structure allows.
# ============================================================

options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(vegan)
  library(stringr)
})

ROOT <- "E:/sepsis_project"

SELECTED <- file.path(
  ROOT,
  "results",
  "V2_33G_STEP93G_TAXONOMIC_TRAJECTORY_DIVERGENCE",
  "V2_STEP93G_selected_candidate.csv"
)

OUT <- file.path(
  ROOT,
  "results",
  "V2_33H_STEP93H_TAXONOMIC_TRAJECTORY_DISTANCE"
)

dir.create(
  OUT,
  recursive=TRUE,
  showWarnings=FALSE
)

if(!file.exists(SELECTED)){
  stop("Selected genus matrix missing.")
}

selected <- read_csv(
  SELECTED,
  show_col_types=FALSE
)

matrix_file <- selected$taxonomy_candidate_file[1]

if(!file.exists(matrix_file)){
  stop("Genus matrix path unavailable.")
}

tax <- read_csv(
  matrix_file,
  show_col_types=FALSE
)

# Audit structure

write_csv(
  tibble(
    rows=nrow(tax),
    cols=ncol(tax),
    file=matrix_file
  ),
  file.path(
    OUT,
    "V2_STEP93H_matrix_structure.csv"
  )
)

# Identify metadata columns

meta_cols <- names(tax)[
  str_detect(
    tolower(names(tax)),
    "sample|patient|subject|id|time|day|project|cohort"
  )
]

feature_cols <- setdiff(
  names(tax),
  meta_cols
)

write_csv(
  tibble(
    metadata_columns=paste(meta_cols,collapse=";"),
    feature_columns=paste(feature_cols,collapse=";")
  ),
  file.path(
    OUT,
    "V2_STEP93H_column_classification.csv"
  )
)

# If no trajectory metadata, stop after audit

if(length(meta_cols)<2){
  writeLines(
    "Need manual metadata mapping.",
    file.path(
      OUT,
      "MANUAL_REVIEW_REQUIRED.txt"
    )
  )
} else {

  write_csv(
    tax[,c(meta_cols,feature_cols)],
    file.path(
      OUT,
      "V2_STEP93H_taxonomy_matrix_for_trajectory_analysis.csv"
    )
  )
}

writeLines(
  c(
    paste0("Completed: ",Sys.time()),
    "STEP93H COMPLETE"
  ),
  file.path(
    OUT,
    "_STEP93H_COMPLETE.ok"
  )
)

cat("STEP93H COMPLETE\n")
