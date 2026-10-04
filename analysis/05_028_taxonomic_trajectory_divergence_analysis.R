
# ============================================================
# Sepsis V2 - Step93R
# Taxonomic trajectory divergence analysis
#
# Input:
# V2_STEP93Q_fixed_longitudinal_genus_profile.csv
#
# Goal:
# 1. Build patient-level longitudinal genus trajectories
# 2. Calculate baseline-followup genus change vectors
# 3. Calculate Bray-Curtis distance between patients
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

IN <- file.path(
  ROOT,
  "results",
  "V2_33Q_TRAJECTORY_PATIENT_ID_RECONSTRUCTION_FIX",
  "V2_STEP93Q_fixed_longitudinal_genus_profile.csv"
)

OUT <- file.path(
  ROOT,
  "results",
  "V2_33R_TAXONOMIC_TRAJECTORY_DIVERGENCE_ANALYSIS"
)

dir.create(
  OUT,
  recursive=TRUE,
  showWarnings=FALSE
)

if(!file.exists(IN)){
  stop("Missing Step93Q input.")
}

dat <- read_csv(
  IN,
  show_col_types=FALSE
)

meta_cols <- c(
  "Project",
  "Run",
  "BioSample",
  "Sample_Name",
  "Source",
  "Group",
  "Patient_ID",
  "Patient_ID_reconstructed",
  "Timepoint_reconstructed",
  "Patient_true"
)

genus_cols <- setdiff(
  names(dat),
  meta_cols
)


# ------------------------------------------------------------
# Keep only repeated patients
# ------------------------------------------------------------

patients <- dat %>%
  group_by(Patient_true) %>%
  summarise(
    n_timepoints=n_distinct(Time_numeric),
    .groups="drop"
  ) %>%
  filter(
    n_timepoints >= 2
  )


dat2 <- dat %>%
  filter(
    Patient_true %in% patients$Patient_true
  )


# ------------------------------------------------------------
# Baseline and follow-up selection
# ------------------------------------------------------------

baseline <- dat2 %>%
  group_by(Patient_true) %>%
  filter(
    Time_numeric == min(Time_numeric, na.rm=TRUE)
  ) %>%
  slice(1) %>%
  ungroup()


followup <- dat2 %>%
  group_by(Patient_true) %>%
  filter(
    Time_numeric == max(Time_numeric, na.rm=TRUE)
  ) %>%
  slice(1) %>%
  ungroup()


# ------------------------------------------------------------
# Delta genus trajectory matrix
# ------------------------------------------------------------

delta_list <- lapply(
  seq_len(nrow(baseline)),
  function(i){

    b <- baseline[i, genus_cols]
    f <- followup[i, genus_cols]

    as.numeric(f) - as.numeric(b)
  }
)

delta <- as.data.frame(
  do.call(rbind, delta_list)
)

names(delta) <- genus_cols

delta$Patient_true <- baseline$Patient_true


write_csv(
  delta,
  file.path(
    OUT,
    "V2_STEP93R_patient_genus_delta_matrix.csv"
  )
)


# ------------------------------------------------------------
# Bray-Curtis distance
# ------------------------------------------------------------

delta_matrix <- delta %>%
  select(
    all_of(genus_cols)
  )

rownames(delta_matrix) <- delta$Patient_true


# Convert negative changes to comparable abundance vectors
delta_distance_matrix <- vegdist(
  abs(as.matrix(delta_matrix)),
  method="bray"
)

write.csv(
  as.matrix(delta_distance_matrix),
  file.path(
    OUT,
    "V2_STEP93R_taxonomic_trajectory_bray_distance_matrix.csv"
  )
)


heterogeneity <- tibble(
  mean_pairwise_distance =
    mean(
      as.vector(delta_distance_matrix)
    ),
  n_patients=nrow(delta)
)

write_csv(
  heterogeneity,
  file.path(
    OUT,
    "V2_STEP93R_taxonomic_trajectory_heterogeneity_index.csv"
  )
)


writeLines(
  c(
    paste0("Completed: ",Sys.time()),
    "STEP93R COMPLETE"
  ),
  file.path(
    OUT,
    "_STEP93R_COMPLETE.ok"
  )
)

cat("STEP93R COMPLETE\n")
