# ============================================================
# Sepsis V2 - Step 70B
# Finish Step70 after exact public-ID matching
# R 4.4.0 / Windows
#
# This script does NOT redownload anything.
# It uses Step70 files already produced.
# ============================================================

options(stringsAsFactors = FALSE)

pkgs <- c("readr","dplyr","stringr","tibble")
missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly=TRUE)]

if (length(missing) > 0) {
  install.packages(missing, repos="https://cloud.r-project.org")
}

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(stringr)
  library(tibble)
})

ROOT <- "E:/sepsis_project/results/V2_10_exact_metadata_mapping"

candidate_path <- file.path(
  ROOT,
  "PRJEB67798_exact_ID_match_candidates.csv"
)

ambiguity_path <- file.path(
  ROOT,
  "PRJEB67798_exact_match_ambiguity_QC.csv"
)

map824_path <- file.path(
  ROOT,
  "PRJEB82425_explicit_sample_patient_time_map.csv"
)

if (!file.exists(candidate_path)) {
  stop("Missing: PRJEB67798_exact_ID_match_candidates.csv")
}

if (!file.exists(ambiguity_path)) {
  stop("Missing: PRJEB67798_exact_match_ambiguity_QC.csv")
}

if (!file.exists(map824_path)) {
  stop("Missing: PRJEB82425_explicit_sample_patient_time_map.csv")
}

cand <- suppressMessages(
  read_csv(candidate_path, show_col_types=FALSE)
)

amb <- suppressMessages(
  read_csv(ambiguity_path, show_col_types=FALSE)
)

map824 <- suppressMessages(
  read_csv(map824_path, show_col_types=FALSE)
)

# ------------------------------------------------------------
# PRJEB67798 exact map
# ------------------------------------------------------------
good <- amb |>
  filter(
    N_ENA_Rows == 1,
    N_Runs == 1
  ) |>
  select(
    Public_Sample_ID,
    Public_Patient_ID
  )

map677 <- cand |>
  inner_join(
    good,
    by=c(
      "Public_Sample_ID",
      "Public_Patient_ID"
    )
  ) |>
  filter(
    !is.na(run_accession),
    !is.na(Public_Patient_ID),
    !is.na(Public_Time_Condition),
    Public_Patient_ID != "categorical"
  ) |>
  mutate(
    Time_Label = case_when(
      str_detect(
        sample_description,
        regex("before", ignore_case=TRUE)
      ) ~ "Before_surgery",

      str_detect(
        sample_description,
        regex("after", ignore_case=TRUE)
      ) ~ "After_surgery",

      TRUE ~ Public_Time_Condition
    )
  ) |>
  transmute(
    Project = "PRJEB67798",

    Patient_ID = Public_Patient_ID,

    Sample_ID = Public_Sample_ID,

    Run_ID = run_accession,

    Time_Raw = Public_Time_Condition,

    Time_Label = Time_Label,

    Body_Site = "rectal",

    Sample_Description = sample_description,

    Collection_Date_Raw = as.character(collection_date),

    Crosswalk_Source = Crosswalk_Source,

    ENA_ID_Field = ENA_ID_Field,

    Mapping_Method =
      "Exact public BORIS SampleID -> ENA identifier"
  ) |>
  distinct()

write_excel_csv(
  map677,
  file.path(
    ROOT,
    "PRJEB67798_exact_sample_patient_run_map.csv"
  ),
  na=""
)

# ------------------------------------------------------------
# Longitudinal QC helpers
# ------------------------------------------------------------
long_qc <- function(df, project) {

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
      name="N_Timepoints"
    )

  tibble(
    Project = project,

    Rows = nrow(df),

    Unique_Patients = n_distinct(
      df$Patient_ID,
      na.rm=TRUE
    ),

    Unique_Samples = n_distinct(
      df$Sample_ID,
      na.rm=TRUE
    ),

    Unique_Runs = n_distinct(
      df$Run_ID,
      na.rm=TRUE
    ),

    Missing_Patient = sum(
      is.na(df$Patient_ID)
    ),

    Missing_Time = sum(
      is.na(df$Time_Raw)
    ),

    Missing_Run = sum(
      is.na(df$Run_ID)
    ),

    Patients_GE2_Timepoints = sum(
      pt$N_Timepoints >= 2
    ),

    Patients_GE3_Timepoints = sum(
      pt$N_Timepoints >= 3
    )
  )
}

qc677 <- long_qc(
  map677,
  "PRJEB67798"
)

qc824 <- long_qc(
  map824,
  "PRJEB82425"
) |>
  mutate(
    Body_Site_Resolved = sum(
      !is.na(map824$Body_Site)
    ),

    Body_Site_Unresolved = sum(
      is.na(map824$Body_Site)
    )
  )

qc <- bind_rows(
  qc677,
  qc824
)

write_excel_csv(
  qc,
  file.path(
    ROOT,
    "V2_step70B_mapping_QC.csv"
  ),
  na=""
)

# ------------------------------------------------------------
# Project decision
# ------------------------------------------------------------
decision <- tibble(
  Project = c(
    "PRJEB67798",
    "PRJEB82425"
  ),

  Patient_Mapping = c(
    "COMPLETE_FOR_40_EXACT_PUBLIC_SAMPLES",
    "COMPLETE_FOR_ALL_165_RUNS"
  ),

  Time_Mapping = c(
    "COMPLETE_T1_T2_FOR_20_PATIENTS",
    "COMPLETE_EVENT_TIME_FOR_ALL_165_RUNS"
  ),

  Body_Site = c(
    "RECTAL_EXPLICIT_FROM_SAMPLE_DESCRIPTION",
    "UNRESOLVED_FOR_165_RUNS"
  ),

  Longitudinal_Use = c(
    "YES_SURGICAL_NONSEPSIS_CONTROL",
    "YES_AFTER_GUT_VS_TRACHEAL_CLASSIFICATION"
  ),

  Manual_Metadata_Collection_Now = c(
    "NO",
    "NO"
  ),

  Next_Action = c(
    "Include the 40 exact rectal samples in master metadata.",
    paste(
      "Automatically classify rectal vs tracheal using raw-read",
      "amplicon/primer evidence; do not infer from numeric sample_alias."
    )
  )
)

write_excel_csv(
  decision,
  file.path(
    ROOT,
    "V2_metadata_decision_after70B.csv"
  ),
  na=""
)

cat("\n=============================================\n")
cat("SEPSIS V2 STEP 70B COMPLETE\n")
cat("=============================================\n\n")

cat("PRJEB67798 exact mapping:\n")
print(qc677, n=Inf, width=Inf)

cat("\nPRJEB82425 mapping:\n")
print(qc824, n=Inf, width=Inf)

cat("\nDecision:\n")
print(decision, n=Inf, width=Inf)

cat("\nOutputs:\n")
cat(file.path(ROOT,"PRJEB67798_exact_sample_patient_run_map.csv"),"\n")
cat(file.path(ROOT,"V2_step70B_mapping_QC.csv"),"\n")
cat(file.path(ROOT,"V2_metadata_decision_after70B.csv"),"\n")
