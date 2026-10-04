# ============================================================
# Sepsis V2 - Step92D
# EII stability comparison
#
# Questions:
# 1. Does EII reproduce Bray displacement?
# 2. Is EII trajectory more conserved than taxonomy signatures?
# 3. Leave-one-cohort-out robustness preparation
# ============================================================

options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(stats)
})

ROOT <- "E:/sepsis_project"

EII_FILE <- file.path(
  ROOT,
  "results",
  "V2_32B2_STEP92B_BASELINE_HARMONIZATION_FIXED",
  "V2_STEP92B2_corrected_timepoint_EII.csv"
)

BRAY_FILE <- file.path(
  ROOT,
  "results",
  "V2_28A2_STEP88A_LONGITUDINAL_DIVERSITY_AND_DISPLACEMENT_FIXED",
  "03_BETA_DISPLACEMENT",
  "V2_STEP88A2_ALL_within_patient_bray_displacement.csv"
)

OUT <- file.path(
  ROOT,
  "results",
  "V2_32D_STEP92D_EII_STABILITY_COMPARISON"
)

dir.create(
  OUT,
  recursive=TRUE,
  showWarnings=FALSE
)

if(!file.exists(EII_FILE)){
  stop("Missing EII file.")
}

if(!file.exists(BRAY_FILE)){
  stop("Missing Bray file.")
}

eii <- read_csv(
  EII_FILE,
  show_col_types=FALSE
)

bray <- read_csv(
  BRAY_FILE,
  show_col_types=FALSE
)

# ------------------------------------------------------------
# Patient-level comparison
# ------------------------------------------------------------

eii_patient <- eii %>%
  group_by(
    project,
    patient_id
  ) %>%
  summarise(
    final_EII =
      EII_0_100[
        which.max(
          row_number()
        )
      ][1],
    max_EII =
      max(
        EII_0_100,
        na.rm=TRUE
      ),
    .groups="drop"
  )

bray_patient <- bray %>%
  group_by(
    project,
    patient_id
  ) %>%
  summarise(
    max_bray =
      max(
        bray_from_patient_baseline,
        na.rm=TRUE
      ),
    .groups="drop"
  )

comparison <- eii_patient %>%
  inner_join(
    bray_patient,
    by=c(
      "project",
      "patient_id"
    )
  )

write_csv(
  comparison,
  file.path(
    OUT,
    "V2_STEP92D_EII_vs_Bray_patient_comparison.csv"
  )
)

cor_test <- cor.test(
  comparison$max_EII,
  comparison$max_bray,
  method="spearman",
  exact=FALSE
)

cor_result <- tibble(
  method="Spearman",
  rho=
    unname(cor_test$estimate),
  p_value=
    cor_test$p.value,
  n=
    nrow(comparison)
)

write_csv(
  cor_result,
  file.path(
    OUT,
    "V2_STEP92D_EII_Bray_correlation.csv"
  )
)

# ------------------------------------------------------------
# Leave-one-cohort-out preparation
# ------------------------------------------------------------

projects <- unique(
  comparison$project
)

loo <- lapply(
  projects,
  function(drop){

    remain <- comparison %>%
      filter(
        project != drop
      )

    tibble(
      excluded_cohort=drop,
      n=nrow(remain),
      rho=
        cor(
          remain$max_EII,
          remain$max_bray,
          method="spearman"
        )
    )
  }
) %>%
  bind_rows()

write_csv(
  loo,
  file.path(
    OUT,
    "V2_STEP92D_leave_one_cohort_out_EII_Bray.csv"
  )
)

writeLines(
  c(
    paste0(
      "Completed: ",
      Sys.time()
    ),
    "STEP92D COMPLETE"
  ),
  file.path(
    OUT,
    "_STEP92D_COMPLETE.ok"
  )
)

cat("STEP92D COMPLETE\n")
