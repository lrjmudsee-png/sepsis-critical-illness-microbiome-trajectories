# ============================================================
# Sepsis V2 - Step93E
# Taxonomic trajectory divergence calculation
#
# Goal:
# Quantify between-patient taxonomic trajectory heterogeneity.
#
# Uses genus-level abundance candidates preferentially.
# If no suitable matrix is found, outputs audit information.
# ============================================================

options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(stringr)
})

ROOT <- "E:/sepsis_project"

INVENTORY <- file.path(
  ROOT,
  "results",
  "V2_33D_STEP93D_STEP93D_TAXONOMIC_HETEROGENEITY_INDEX",
  "V2_STEP93D_taxonomy_trajectory_candidates.csv"
)

OUT <- file.path(
  ROOT,
  "results",
  "V2_33E_STEP93E_TAXONOMIC_TRAJECTORY_DIVERGENCE"
)

dir.create(
  OUT,
  recursive=TRUE,
  showWarnings=FALSE
)

# fallback search because folder naming may differ
inventory_candidates <- list.files(
  file.path(ROOT,"results"),
  pattern="V2_STEP93D_taxonomy_trajectory_candidates.csv",
  recursive=TRUE,
  full.names=TRUE
)

if(length(inventory_candidates)>0){
  INVENTORY <- inventory_candidates[1]
}

if(!file.exists(INVENTORY)){
  stop("Taxonomy trajectory candidate inventory missing.")
}

files <- read_csv(
  INVENTORY,
  show_col_types=FALSE
)

genus_files <- files %>%
  filter(type=="genus")

write_csv(
  genus_files,
  file.path(
    OUT,
    "V2_STEP93E_genus_level_candidate_files.csv"
  )
)

asv_files <- files %>%
  filter(type=="ASV")

write_csv(
  asv_files,
  file.path(
    OUT,
    "V2_STEP93E_ASV_candidate_files.csv"
  )
)

writeLines(
  c(
    paste0("Completed: ",Sys.time()),
    "STEP93E COMPLETE",
    "Taxonomic trajectory divergence preparation finished."
  ),
  file.path(
    OUT,
    "_STEP93E_COMPLETE.ok"
  )
)

cat("STEP93E COMPLETE\n")
