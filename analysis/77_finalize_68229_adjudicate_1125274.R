# ============================================================
# Sepsis V2 - Step 77
# Finalize PRJEB68229 and adjudicate PRJNA1125274
# OFFLINE ONLY - no network access
# R 4.4.0 / Windows
# ============================================================

options(stringsAsFactors = FALSE)

pkgs <- c("readr","dplyr","tidyr","stringr","tibble")

missing <- pkgs[
  !vapply(pkgs, requireNamespace, logical(1), quietly=TRUE)
]

if (length(missing) > 0) {
  install.packages(
    missing,
    repos="https://cloud.r-project.org"
  )
}

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(tibble)
})

PROJECT_ROOT <- "E:/sepsis_project"

IN_ROOT <- file.path(
  PROJECT_ROOT,
  "results",
  "V2_16D_high_value_metadata_rescue_fixed"
)

OUT_ROOT <- file.path(
  PROJECT_ROOT,
  "results",
  "V2_17_finalize_68229_adjudicate_1125274"
)

dir.create(
  OUT_ROOT,
  recursive=TRUE,
  showWarnings=FALSE
)

# ------------------------------------------------------------
# Helpers
# ------------------------------------------------------------
clean_chr <- function(x) {
  x <- trimws(as.character(x))
  x[
    is.na(x) |
    x == "" |
    tolower(x) %in% c("na","nan","n/a","null","none")
  ] <- NA_character_
  x
}

read_required <- function(name) {

  path <- file.path(
    IN_ROOT,
    name
  )

  if (!file.exists(path)) {
    stop(
      paste0(
        "Required file not found: ",
        path
      )
    )
  }

  suppressMessages(
    read_csv(
      path,
      show_col_types=FALSE,
      progress=FALSE
    )
  )
}

# ============================================================
# PART A. PRJEB68229
# ============================================================
map682 <- read_required(
  "PRJEB68229_patient_time_map_CANDIDATE_HIGH_CONFIDENCE.csv"
)

logical682 <- read_required(
  "PRJEB68229_logical_samples_deduplicated.csv"
)

# ------------------------------------------------------------
# 1. Final patient/time mapping
#
# ENA raw labels are D1/D3.
# Public supplementary material defines the study observations
# as:
#   S1 = first stool sample
#   S2 = second stool sample
#
# The dataset contains exactly:
#   96 first observations + 94 second observations = 190
#
# Therefore:
#   D1 -> S1
#   D3 -> S2
#
# IMPORTANT:
# D1/D3 are preserved as raw submission labels.
# They are NOT interpreted as literal ICU day 1/day 3.
# ------------------------------------------------------------
final682 <- map682 |>
  mutate(
    Patient_ID = clean_chr(Patient_ID),

    Time_Raw = clean_chr(Time_Raw),

    Time_Standard = case_when(
      Time_Raw == "D1" ~ "S1_first_stool",
      Time_Raw == "D3" ~ "S2_second_stool",
      TRUE ~ NA_character_
    ),

    Time_Order = case_when(
      Time_Raw == "D1" ~ 1L,
      Time_Raw == "D3" ~ 2L,
      TRUE ~ NA_integer_
    ),

    Time_Anchor = case_when(
      Time_Raw == "D1" ~
        "First_stool_after_ICU_admission",

      Time_Raw == "D3" ~
        "Second_stool_at_least_24h_after_S1",

      TRUE ~
        NA_character_
    ),

    Time_Day = NA_real_,

    Project_Role =
      "ICU_BACKGROUND_LONGITUDINAL",

    Analysis_Tier =
      "LONGITUDINAL_VALIDATION",

    Patient_Level_Mortality_Available =
      FALSE,

    Outcome_Note =
      paste(
        "Study reports aggregate 60-day mortality,",
        "but current public Run-patient map does not identify",
        "which individual patients were survivors/non-survivors."
      ),

    Mapping_Status =
      "CONFIRMED_FOR_PATIENT_AND_LONGITUDINAL_ORDER",

    Mapping_Evidence =
      paste(
        "ENA aliases provide repeated patient structure;",
        "supplement defines S1/S2;",
        "96 patients and 190 longitudinal observations reconcile exactly."
      )
  )

# ------------------------------------------------------------
# 2. Collection-date interval QC
# ------------------------------------------------------------
dates682 <- final682 |>
  mutate(
    Collection_Date_Parsed =
      as.Date(Collection_Date)
  ) |>
  select(
    Patient_ID,
    Time_Raw,
    Collection_Date_Parsed
  ) |>
  distinct() |>
  pivot_wider(
    names_from=Time_Raw,
    values_from=Collection_Date_Parsed
  ) |>
  mutate(
    Interval_Days =
      as.numeric(D3 - D1),

    Collection_Date_QC = case_when(
      is.na(D1) | is.na(D3) ~
        "ONE_TIMEPOINT_ONLY",

      Interval_Days < 0 ~
        "SOURCE_DATE_ORDER_ANOMALY",

      TRUE ~
        "PASS"
    )
  )

final682 <- final682 |>
  left_join(
    dates682 |>
      select(
        Patient_ID,
        Interval_Days,
        Collection_Date_QC
      ),
    by="Patient_ID"
  )

write_excel_csv(
  final682,
  file.path(
    OUT_ROOT,
    "PRJEB68229_FINAL_patient_time_map.csv"
  ),
  na=""
)

write_excel_csv(
  dates682,
  file.path(
    OUT_ROOT,
    "PRJEB68229_collection_date_QC.csv"
  ),
  na=""
)

# ------------------------------------------------------------
# 3. NEC control QC
# ------------------------------------------------------------
nec682 <- logical682 |>
  filter(
    str_detect(
      sample_alias,
      regex("NEC", ignore_case=TRUE)
    )
  )

write_excel_csv(
  nec682,
  file.path(
    OUT_ROOT,
    "PRJEB68229_NEC_controls.csv"
  ),
  na=""
)

# ------------------------------------------------------------
# 4. Longitudinal QC
# ------------------------------------------------------------
pt682 <- final682 |>
  distinct(
    Patient_ID,
    Time_Standard
  ) |>
  count(
    Patient_ID,
    name="N_Timepoints"
  )

qc682 <- tibble(
  Metric=c(
    "Patient_runs",
    "Patients",
    "Patients_S1",
    "Patients_S2",
    "Patients_GE2_timepoints",
    "Patients_GE3_timepoints",
    "NEC_controls",
    "Negative_collection_date_intervals",
    "Median_S1_to_S2_days_valid_pairs",
    "Patient_level_mortality_linkage_available"
  ),

  Value=c(
    nrow(final682),

    n_distinct(
      final682$Patient_ID
    ),

    n_distinct(
      final682$Patient_ID[
        final682$Time_Standard ==
          "S1_first_stool"
      ]
    ),

    n_distinct(
      final682$Patient_ID[
        final682$Time_Standard ==
          "S2_second_stool"
      ]
    ),

    sum(
      pt682$N_Timepoints >= 2
    ),

    sum(
      pt682$N_Timepoints >= 3
    ),

    nrow(nec682),

    sum(
      dates682$Collection_Date_QC ==
        "SOURCE_DATE_ORDER_ANOMALY",
      na.rm=TRUE
    ),

    median(
      dates682$Interval_Days[
        dates682$Interval_Days >= 0
      ],
      na.rm=TRUE
    ),

    FALSE
  )
)

write_excel_csv(
  qc682,
  file.path(
    OUT_ROOT,
    "PRJEB68229_FINAL_QC.csv"
  ),
  na=""
)

# ============================================================
# PART B. PRJNA1125274
# ============================================================
audit112 <- read_required(
  "PRJNA1125274_alias_pattern_audit_NOT_FINAL.csv"
)

recon112 <- read_required(
  "PRJNA1125274_run_reconciliation.csv"
)

# ------------------------------------------------------------
# 5. Raw-data adjudication
# ------------------------------------------------------------
current_runs112 <- recon112 |>
  filter(
    In_Current_ENA %in% TRUE
  ) |>
  select(
    Run_ID,
    In_Current_ENA,
    In_Local_FASTQ,
    In_SRR_Acc_List,
    Status
  ) |>
  arrange(Run_ID)

wrong_local112 <- recon112 |>
  filter(
    In_Local_FASTQ %in% TRUE,
    !(In_Current_ENA %in% TRUE)
  ) |>
  select(
    Run_ID,
    everything()
  ) |>
  arrange(Run_ID)

write_excel_csv(
  current_runs112,
  file.path(
    OUT_ROOT,
    "PRJNA1125274_CORRECT_current_289_Run_list.csv"
  ),
  na=""
)

write_excel_csv(
  wrong_local112,
  file.path(
    OUT_ROOT,
    "PRJNA1125274_LOCAL_308_runs_WRONG_FOR_THIS_PROJECT.csv"
  ),
  na=""
)

# ------------------------------------------------------------
# 6. Baseline-defined patient candidates
#
# This remains a PROVISIONAL planning map:
# patient candidate = site + numeric prefix
#
# Only candidates having an A sample are retained.
# This avoids inventing patient identity for later-only samples.
#
# A/B/C interpretation is strong study-design inference, not an
# explicit ENA patient/time field. It MUST remain labelled as such.
# ------------------------------------------------------------
baseline_patients112 <- audit112 |>
  filter(
    alias_letter == "A",
    !is.na(patient_candidate)
  ) |>
  distinct(
    patient_candidate
  )

strict112 <- audit112 |>
  semi_join(
    baseline_patients112,
    by="patient_candidate"
  ) |>
  mutate(
    Patient_ID_Provisional =
      patient_candidate,

    Time_Raw =
      alias_letter,

    Time_Standard_Candidate = case_when(
      alias_letter == "A" ~
        "T0_baseline_admission_candidate",

      alias_letter == "B" ~
        "T1_sepsis_diagnosis_candidate",

      alias_letter == "C" ~
        "T2_discharge_or_death_candidate",

      TRUE ~
        NA_character_
    ),

    Time_Order_Candidate = case_when(
      alias_letter == "A" ~ 1L,
      alias_letter == "B" ~ 2L,
      alias_letter == "C" ~ 3L,
      TRUE ~ NA_integer_
    ),

    Patient_Mapping_Status =
      "HIGH_CONFIDENCE_INFERENCE_NOT_EXPLICIT",

    Time_Mapping_Status =
      "HIGH_CONFIDENCE_STUDY_DESIGN_INFERENCE_NOT_EXPLICIT",

    Raw_Data_Status =
      "CORRECT_CURRENT_289_RUNS_NOT_PRESENT_LOCALLY",

    Analysis_Status =
      "DO_NOT_RUN_DADA2_UNTIL_CORRECT_289_RUNS_DOWNLOADED"
  )

later_only112 <- audit112 |>
  anti_join(
    baseline_patients112,
    by="patient_candidate"
  )

write_excel_csv(
  strict112,
  file.path(
    OUT_ROOT,
    "PRJNA1125274_STRICT_PROVISIONAL_patient_time_map.csv"
  ),
  na=""
)

write_excel_csv(
  later_only112,
  file.path(
    OUT_ROOT,
    "PRJNA1125274_later_only_samples_EXCLUDED_FROM_STRICT_MAP.csv"
  ),
  na=""
)

# ------------------------------------------------------------
# 7. Collection-date ordering diagnostic
# This is diagnostic only: dates are not edited.
# ------------------------------------------------------------
date112 <- strict112 |>
  mutate(
    Collection_Date_Parsed =
      as.Date(collection_date)
  ) |>
  select(
    Patient_ID_Provisional,
    Time_Raw,
    Collection_Date_Parsed
  ) |>
  distinct() |>
  pivot_wider(
    names_from=Time_Raw,
    values_from=Collection_Date_Parsed
  ) |>
  mutate(
    B_before_A =
      !is.na(A) & !is.na(B) & B < A,

    C_before_A =
      !is.na(A) & !is.na(C) & C < A,

    C_before_B =
      !is.na(B) & !is.na(C) & C < B,

    Any_Date_Order_Flag =
      B_before_A |
      C_before_A |
      C_before_B
  )

write_excel_csv(
  date112 |>
    filter(
      Any_Date_Order_Flag
    ),
  file.path(
    OUT_ROOT,
    "PRJNA1125274_collection_date_order_flags.csv"
  ),
  na=""
)

# ------------------------------------------------------------
# 8. Provisional longitudinal QC
# ------------------------------------------------------------
pt112 <- strict112 |>
  distinct(
    Patient_ID_Provisional,
    Time_Raw
  ) |>
  count(
    Patient_ID_Provisional,
    name="N_Timepoints"
  )

qc112 <- tibble(
  Metric=c(
    "Current_ENA_runs",
    "Local_runs_in_current_ENA",
    "Local_runs_wrong_for_current_project",
    "Strict_provisional_rows",
    "Baseline_defined_patient_candidates",
    "Patients_GE2_timepoints_candidate",
    "Patients_GE3_timepoints_candidate",
    "Later_only_excluded_rows",
    "Collection_date_order_flag_patients",
    "Explicit_ENA_patient_fields",
    "Explicit_ENA_time_fields",
    "Ready_for_DADA2_now"
  ),

  Value=c(
    nrow(current_runs112),

    sum(
      current_runs112$In_Local_FASTQ %in% TRUE
    ),

    nrow(wrong_local112),

    nrow(strict112),

    n_distinct(
      strict112$Patient_ID_Provisional
    ),

    sum(
      pt112$N_Timepoints >= 2
    ),

    sum(
      pt112$N_Timepoints >= 3
    ),

    nrow(later_only112),

    sum(
      date112$Any_Date_Order_Flag,
      na.rm=TRUE
    ),

    sum(
      !is.na(
        audit112$patient_id
      )
    ),

    sum(
      !is.na(
        audit112$timepoint
      )
    ),

    FALSE
  )
)

write_excel_csv(
  qc112,
  file.path(
    OUT_ROOT,
    "PRJNA1125274_FINAL_ADJUDICATION_QC.csv"
  ),
  na=""
)

# ============================================================
# PART C. Decision table
# ============================================================
decision <- tibble(
  Project=c(
    "PRJEB68229",
    "PRJNA1125274"
  ),

  Metadata_Status=c(
    "PATIENT_AND_LONGITUDINAL_ORDER_CONFIRMED",
    "STRICT_PROVISIONAL_MAPPING_ONLY"
  ),

  Sequence_Status=c(
    "EXISTING_SEQUENCE_DATA_CAN_BE_USED_IF_PROJECT_FILES_PRESENT",
    "CORRECT_CURRENT_289_RUNS_NOT_LOCAL"
  ),

  Recommended_Role=c(
    "ICU_BACKGROUND_LONGITUDINAL",
    "CORE_SEPSIS_LONGITUDINAL_CANDIDATE_AFTER_RAW_REPAIR"
  ),

  Outcome_Analysis=c(
    "NO_PATIENT_LEVEL_MORTALITY_LINKAGE_IN_CURRENT_PUBLIC_MAP",
    "NOT_YET_ADJUDICATED"
  ),

  Immediate_Action=c(
    "ADD_TO_METADATA_MASTER_AS_LONGITUDINAL_BACKGROUND",
    "DOWNLOAD_CORRECT_289_RUNS_BEFORE_ANY_DADA2"
  )
)

write_excel_csv(
  decision,
  file.path(
    OUT_ROOT,
    "V2_step77_decision.csv"
  ),
  na=""
)

cat("\n============================================================\n")
cat("SEPSIS V2 - STEP 77 COMPLETE\n")
cat("============================================================\n\n")

cat("PRJEB68229 FINAL QC:\n")
print(qc682, n=Inf, width=Inf)

cat("\nPRJNA1125274 ADJUDICATION QC:\n")
print(qc112, n=Inf, width=Inf)

cat("\nDECISION:\n")
print(decision, n=Inf, width=Inf)

cat(
  "\nOutput folder:\n",
  OUT_ROOT,
  "\n",
  sep=""
)

cat("============================================================\n")
