# ============================================================
# Sepsis V2 - Step93D
# Taxonomic heterogeneity index preparation
#
# Goal:
# Build patient-level taxonomic trajectory similarity analysis.
#
# This step discovers available taxonomy matrices and prepares
# input tables for heterogeneity calculation.
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
  "V2_33C_STEP93C_TAXONOMIC_DIVERGENCE_VS_EII",
  "V2_STEP93C_taxonomy_candidate_files.csv"
)

OUT <- file.path(
  ROOT,
  "results",
  "V2_33D_STEP93D_TAXONOMIC_HETEROGENEITY_INDEX"
)

dir.create(
  OUT,
  recursive=TRUE,
  showWarnings=FALSE
)

if(!file.exists(INVENTORY)){
  stop("Taxonomy inventory missing.")
}

files <- read_csv(
  INVENTORY,
  show_col_types=FALSE
)

# Classify possible taxonomy objects

taxonomy_inventory <- files %>%
  mutate(
    type =
      case_when(
        str_detect(
          taxonomy_candidate_file,
          "genus"
        ) ~ "genus",

        str_detect(
          taxonomy_candidate_file,
          "family"
        ) ~ "family",

        str_detect(
          taxonomy_candidate_file,
          "ASV|asv"
        ) ~ "ASV",

        TRUE ~ "unknown"
      )
  )

write_csv(
  taxonomy_inventory,
  file.path(
    OUT,
    "V2_STEP93D_taxonomy_inventory_classified.csv"
  )
)

# Keep candidates for downstream trajectory matrix construction

trajectory_candidates <- taxonomy_inventory %>%
  filter(
    type!="unknown"
  )

write_csv(
  trajectory_candidates,
  file.path(
    OUT,
    "V2_STEP93D_taxonomy_trajectory_candidates.csv"
  )
)

writeLines(
  c(
    paste0("Completed: ",Sys.time()),
    "STEP93D COMPLETE",
    "Taxonomic heterogeneity index preparation finished."
  ),
  file.path(
    OUT,
    "_STEP93D_COMPLETE.ok"
  )
)

cat("STEP93D COMPLETE\n")
