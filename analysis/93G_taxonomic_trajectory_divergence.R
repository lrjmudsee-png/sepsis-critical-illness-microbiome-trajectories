# ============================================================
# Sepsis V2 - Step93G
# Taxonomic trajectory divergence calculation
#
# Goal:
# Calculate patient-level taxonomic trajectory heterogeneity
# from genus-level longitudinal abundance matrices.
#
# Step93F identified candidate genus matrices.
# This step selects the largest compatible matrix first and audits.
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

CANDIDATE <- file.path(
  ROOT,
  "results",
  "V2_33F_STEP93F_TAXONOMIC_HETEROGENEITY_CALCULATION",
  "V2_STEP93F_genus_candidates_audit.csv"
)

OUT <- file.path(
  ROOT,
  "results",
  "V2_33G_STEP93G_TAXONOMIC_TRAJECTORY_DIVERGENCE"
)

dir.create(
  OUT,
  recursive=TRUE,
  showWarnings=FALSE
)

if(!file.exists(CANDIDATE)){
  stop("Missing genus candidate list.")
}

files <- read_csv(
  CANDIDATE,
  show_col_types=FALSE
)

# Candidate ranking by existence and expected size
audit <- files %>%
  mutate(
    exists = file.exists(taxonomy_candidate_file)
  )

write_csv(
  audit,
  file.path(
    OUT,
    "V2_STEP93G_candidate_audit.csv"
  )
)

# Select candidate for downstream inspection
selected <- audit %>%
  filter(exists) %>%
  slice(1)

write_csv(
  selected,
  file.path(
    OUT,
    "V2_STEP93G_selected_candidate.csv"
  )
)

# Read selected matrix
if(nrow(selected)==0){
  stop("No readable genus matrix.")
}

mat <- read_csv(
  selected$taxonomy_candidate_file[1],
  show_col_types=FALSE
)

write_csv(
  tibble(
    n_rows=nrow(mat),
    n_cols=ncol(mat),
    selected_file=
      selected$taxonomy_candidate_file[1]
  ),
  file.path(
    OUT,
    "V2_STEP93G_selected_matrix_dimension.csv"
  )
)

writeLines(
  c(
    paste0("Completed: ",Sys.time()),
    "STEP93G COMPLETE",
    "Taxonomic trajectory divergence calculation initialized."
  ),
  file.path(
    OUT,
    "_STEP93G_COMPLETE.ok"
  )
)

cat("STEP93G COMPLETE\n")
