# ============================================================
# Sepsis V2 - Step 71B
# Finalize PRJEB82425 body-site classification
# R 4.4.0 / Windows
#
# Rationale:
# - Study contains only two specimen types:
#     rectal swab (V3-V4)
#     endotracheal aspirate / ETA (V1-V2)
# - Step71 FASTQ scan shows two completely separated V3-V4
#   signal clusters.
# - High V3-V4 cluster = rectal.
# - Low V3-V4 cluster = ETA, provided:
#     1) separation is large,
#     2) low-cluster size matches the published ETA count (73),
#     3) total ENA runs are 165, one fewer than published total 166.
#
# This script DOES NOT inspect/download raw FASTQs again.
# ============================================================

options(stringsAsFactors = FALSE)

pkgs <- c(
  "readr",
  "dplyr",
  "stringr",
  "tibble"
)

missing <- pkgs[
  !vapply(
    pkgs,
    requireNamespace,
    logical(1),
    quietly = TRUE
  )
]

if (length(missing) > 0) {
  install.packages(
    missing,
    repos = "https://cloud.r-project.org"
  )
}

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(stringr)
  library(tibble)
})


# ------------------------------------------------------------
# 1. Paths
# ------------------------------------------------------------
STEP71_ROOT <- paste0(
  "E:/sepsis_project/results/",
  "V2_11_PRJEB82425_body_site"
)

STEP70_ROOT <- paste0(
  "E:/sepsis_project/results/",
  "V2_10_exact_metadata_mapping"
)

CLASS_FILE <- file.path(
  STEP71_ROOT,
  "PRJEB82425_FASTQ_body_site_classification.csv"
)

MAP_FILE <- file.path(
  STEP70_ROOT,
  "PRJEB82425_explicit_sample_patient_time_map.csv"
)

OUT_ROOT <- paste0(
  "E:/sepsis_project/results/",
  "V2_11_PRJEB82425_body_site_final"
)

dir.create(
  OUT_ROOT,
  recursive = TRUE,
  showWarnings = FALSE
)


# ------------------------------------------------------------
# 2. Published study counts
# ------------------------------------------------------------
PUBLISHED_RECTAL <- 93L
PUBLISHED_ETA <- 73L
PUBLISHED_TOTAL <- PUBLISHED_RECTAL + PUBLISHED_ETA


# ------------------------------------------------------------
# 3. Load
# ------------------------------------------------------------
if (!file.exists(CLASS_FILE)) {
  stop(
    paste0(
      "Missing Step71 classification: ",
      CLASS_FILE
    )
  )
}

if (!file.exists(MAP_FILE)) {
  stop(
    paste0(
      "Missing Step70 mapping: ",
      MAP_FILE
    )
  )
}


cls <- suppressMessages(
  read_csv(
    CLASS_FILE,
    show_col_types = FALSE
  )
)

meta <- suppressMessages(
  read_csv(
    MAP_FILE,
    show_col_types = FALSE
  )
)


if (!all(
  c(
    "Run_ID",
    "V34_Prop"
  ) %in% names(cls)
)) {
  stop(
    "Step71 classification lacks Run_ID or V34_Prop."
  )
}


# ------------------------------------------------------------
# 4. Detect the largest gap in V3-V4 signal
# ------------------------------------------------------------
signal <- cls |>
  filter(
    !is.na(V34_Prop)
  ) |>
  arrange(
    V34_Prop
  )


if (nrow(signal) < 2) {
  stop(
    "Too few V34_Prop values to identify clusters."
  )
}


gap_tbl <- signal |>
  mutate(
    Next_V34 = lead(V34_Prop),
    Gap = Next_V34 - V34_Prop
  )


gap_row <- gap_tbl |>
  filter(
    !is.na(Gap)
  ) |>
  slice_max(
    Gap,
    n = 1,
    with_ties = FALSE
  )


max_low <- gap_row$V34_Prop[1]
min_high <- gap_row$Next_V34[1]
largest_gap <- gap_row$Gap[1]

threshold <- (
  max_low + min_high
) / 2


# ------------------------------------------------------------
# 5. Preliminary two-cluster assignment
# ------------------------------------------------------------
cls2 <- cls |>
  mutate(
    V34_Cluster = case_when(
      is.na(V34_Prop) ~ "missing",
      V34_Prop > threshold ~ "high_V34",
      TRUE ~ "low_V34"
    )
  )


n_high <- sum(
  cls2$V34_Cluster == "high_V34"
)

n_low <- sum(
  cls2$V34_Cluster == "low_V34"
)


# ------------------------------------------------------------
# 6. Validate before assigning ETA
# ------------------------------------------------------------
# We require:
# - Very large separation between clusters.
# - Low cluster exactly equals published ETA count.
# - High cluster is one fewer than published rectal count,
#   consistent with the 165 deposited runs vs 166 study samples.
SEPARATION_OK <- (
  max_low <= 0.05 &&
  min_high >= 0.50 &&
  largest_gap >= 0.40
)

ETA_COUNT_OK <- (
  n_low == PUBLISHED_ETA
)

RECTAL_COUNT_COMPATIBLE <- (
  n_high == PUBLISHED_RECTAL ||
  n_high == PUBLISHED_RECTAL - 1L
)

TOTAL_COMPATIBLE <- (
  nrow(cls2) == PUBLISHED_TOTAL ||
  nrow(cls2) == PUBLISHED_TOTAL - 1L
)


VALIDATED <- (
  SEPARATION_OK &&
  ETA_COUNT_OK &&
  RECTAL_COUNT_COMPATIBLE &&
  TOTAL_COMPATIBLE
)


if (!VALIDATED) {

  stop(
    paste(
      "Automatic body-site validation FAILED.",
      "Review cluster counts/separation before assigning ETA."
    )
  )
}


# ------------------------------------------------------------
# 7. Final classification
# ------------------------------------------------------------
cls_final <- cls2 |>
  mutate(
    Body_Site_Final = case_when(
      V34_Cluster == "high_V34" ~ "rectal",
      V34_Cluster == "low_V34" ~ "tracheal_ETA",
      TRUE ~ NA_character_
    ),

    Body_Site_Evidence = case_when(

      Body_Site_Final == "rectal" ~
        paste(
          "Strong V3-V4 FASTQ signal;",
          "rectal swabs used V3-V4 in study protocol"
        ),

      Body_Site_Final == "tracheal_ETA" ~
        paste(
          "Low V3-V4 FASTQ cluster;",
          "study contains only rectal(V3-V4) and ETA(V1-V2);",
          "low-cluster count equals published 73 ETA samples"
        ),

      TRUE ~
        NA_character_
    ),

    Body_Site_Confidence = case_when(
      Body_Site_Final == "rectal" ~ "HIGH",
      Body_Site_Final == "tracheal_ETA" ~ "HIGH_BY_DESIGN_AND_CLUSTER",
      TRUE ~ "UNRESOLVED"
    )
  )


write_excel_csv(
  cls_final,
  file.path(
    OUT_ROOT,
    "PRJEB82425_body_site_classification_FINAL.csv"
  ),
  na = ""
)


# ------------------------------------------------------------
# 8. Merge with patient/time metadata
# ------------------------------------------------------------
final_map <- meta |>
  select(
    -any_of(
      c(
        "Body_Site",
        "Body_Site_Mapping_Method"
      )
    )
  ) |>
  left_join(
    cls_final |>
      select(
        Run_ID,
        Body_Site_Final,
        Body_Site_Evidence,
        Body_Site_Confidence,
        V12_Prop,
        V34_Prop,
        Classification_Ratio
      ),
    by = "Run_ID"
  ) |>
  rename(
    Body_Site = Body_Site_Final
  ) |>
  mutate(
    Body_Site_Mapping_Method =
      "FASTQ V3-V4 bimodal cluster + published two-specimen study design"
  )


write_excel_csv(
  final_map,
  file.path(
    OUT_ROOT,
    "PRJEB82425_sample_patient_time_FINAL.csv"
  ),
  na = ""
)


# ------------------------------------------------------------
# 9. Gut / ETA subsets
# ------------------------------------------------------------
gut <- final_map |>
  filter(
    Body_Site == "rectal"
  )


eta <- final_map |>
  filter(
    Body_Site == "tracheal_ETA"
  )


write_excel_csv(
  gut,
  file.path(
    OUT_ROOT,
    "PRJEB82425_GUT_ONLY_FINAL.csv"
  ),
  na = ""
)


write_excel_csv(
  eta,
  file.path(
    OUT_ROOT,
    "PRJEB82425_ETA_ONLY_FINAL.csv"
  ),
  na = ""
)


# ------------------------------------------------------------
# 10. Longitudinal QC helper
# ------------------------------------------------------------
patient_time_qc <- function(df) {

  pt <- df |>
    filter(
      !is.na(Patient_ID),
      !is.na(Time_Raw)
    ) |>
    distinct(
      Patient_ID,
      Time_Raw
    ) |>
    count(
      Patient_ID,
      name = "N_Timepoints"
    )

  tibble(
    Runs = n_distinct(
      df$Run_ID,
      na.rm = TRUE
    ),

    Samples = n_distinct(
      df$Sample_ID,
      na.rm = TRUE
    ),

    Patients = n_distinct(
      df$Patient_ID,
      na.rm = TRUE
    ),

    Patients_GE2_Timepoints = sum(
      pt$N_Timepoints >= 2
    ),

    Patients_GE3_Timepoints = sum(
      pt$N_Timepoints >= 3
    )
  )
}


gut_qc <- patient_time_qc(
  gut
) |>
  mutate(
    Body_Site = "rectal"
  )


eta_qc <- patient_time_qc(
  eta
) |>
  mutate(
    Body_Site = "tracheal_ETA"
  )


long_qc <- bind_rows(
  gut_qc,
  eta_qc
) |>
  select(
    Body_Site,
    everything()
  )


write_excel_csv(
  long_qc,
  file.path(
    OUT_ROOT,
    "PRJEB82425_FINAL_longitudinal_QC.csv"
  ),
  na = ""
)


# ------------------------------------------------------------
# 11. Validation report
# ------------------------------------------------------------
validation <- tibble(
  Metric = c(
    "ENA_runs",
    "Published_total_samples",
    "High_V34_runs",
    "Published_rectal_samples",
    "Low_V34_runs",
    "Published_ETA_samples",
    "Max_low_V34_prop",
    "Min_high_V34_prop",
    "Largest_signal_gap",
    "Cluster_threshold",
    "Validation_passed"
  ),

  Value = c(
    nrow(cls_final),
    PUBLISHED_TOTAL,
    n_high,
    PUBLISHED_RECTAL,
    n_low,
    PUBLISHED_ETA,
    max_low,
    min_high,
    largest_gap,
    threshold,
    VALIDATED
  )
)


write_excel_csv(
  validation,
  file.path(
    OUT_ROOT,
    "PRJEB82425_body_site_validation_report.csv"
  ),
  na = ""
)


# ------------------------------------------------------------
# 12. Missing published sample note
# ------------------------------------------------------------
missing_note <- tibble(
  Study_Published_Total = PUBLISHED_TOTAL,
  ENA_Mapped_Runs = nrow(cls_final),
  Difference = PUBLISHED_TOTAL - nrow(cls_final),
  Published_Rectal = PUBLISHED_RECTAL,
  ENA_Rectal_Classified = n_high,
  Published_ETA = PUBLISHED_ETA,
  ENA_ETA_Classified = n_low,
  Interpretation = paste(
    "ENA contains one fewer run than the 166 samples reported",
    "in the publication. The deficit is in the rectal group:",
    "92 deposited/classified vs 93 reported. Do not invent or",
    "impute the missing sample."
  )
)


write_excel_csv(
  missing_note,
  file.path(
    OUT_ROOT,
    "PRJEB82425_published_vs_ENA_sample_count_note.csv"
  ),
  na = ""
)


# ------------------------------------------------------------
# 13. Console
# ------------------------------------------------------------
cat(
  "\n============================================================\n"
)

cat(
  "SEPSIS V2 - STEP 71B COMPLETE\n"
)

cat(
  "R version: ",
  R.version.string,
  "\n",
  sep = ""
)

cat(
  "============================================================\n\n"
)


cat(
  "FASTQ V3-V4 signal separation:\n"
)

cat(
  "  max low cluster  = ",
  max_low,
  "\n",
  sep = ""
)

cat(
  "  min high cluster = ",
  min_high,
  "\n",
  sep = ""
)

cat(
  "  largest gap      = ",
  largest_gap,
  "\n",
  sep = ""
)

cat(
  "  threshold        = ",
  threshold,
  "\n\n",
  sep = ""
)


cat(
  "Final body-site counts:\n"
)

cat(
  "  rectal       = ",
  nrow(gut),
  "\n",
  sep = ""
)

cat(
  "  tracheal_ETA = ",
  nrow(eta),
  "\n",
  sep = ""
)

cat(
  "  unresolved   = ",
  sum(
    is.na(
      final_map$Body_Site
    )
  ),
  "\n\n",
  sep = ""
)


cat(
  "Longitudinal QC:\n"
)

print(
  long_qc,
  n = Inf,
  width = Inf
)


cat(
  "\nOutput folder:\n",
  OUT_ROOT,
  "\n",
  sep = ""
)

cat(
  "============================================================\n"
)
