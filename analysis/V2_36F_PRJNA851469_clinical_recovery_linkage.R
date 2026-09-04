
# ============================================================
# V2_36F PRJNA851469 CLINICAL RECOVERY + LINKAGE QC
#
# Purpose:
#   Recover publicly available patient-level clinical variables
#   from Schlechte et al. Nature Medicine 2023 Supplementary Table 1
#   and link them to the frozen PRJNA851469 longitudinal Bray table.
#
# This step DOES NOT run final inferential outcome models.
# It determines whether the recovered data are sufficiently linked
# for a temporally valid Day-3 landmark clinical analysis and an
# antibiotic-at-admission sensitivity analysis.
# ============================================================

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(stringr)
  library(tidyr)
})

ROOT <- "E:/sepsis_project"
CODE <- file.path(ROOT, "code", "03_data_processing")
RESULTS <- file.path(ROOT, "results")
OUT <- file.path(RESULTS, "V2_36F_PRJNA851469_CLINICAL_RECOVERY_LINKAGE")
dir.create(OUT, recursive=TRUE, showWarnings=FALSE)

MAP_FILE <- file.path(
  CODE,
  "PRJNA851469_SuppTable1_patient_clinical_mapping_CURATED.csv"
)

if (!file.exists(MAP_FILE)) {
  stop("Curated Supplementary Table 1 mapping not found: ", MAP_FILE)
}

clinical <- read_csv(MAP_FILE, show_col_types=FALSE)

# -------------------------------
# Hard publication cross-checks
# -------------------------------
if (nrow(clinical) != 51) stop("Expected 51 published ICU patients.")
if (sum(clinical$nosocomial_infection_30d, na.rm=TRUE) != 28)
  stop("Nosocomial infection count must equal published 28/51.")
if (sum(clinical$mortality_30d, na.rm=TRUE) != 17)
  stop("Mortality count must equal published 17/51.")
if (sum(clinical$antibiotics_at_icu_admission, na.rm=TRUE) != 28)
  stop("Admission antibiotic count must equal published 28/51.")
if (sum(clinical$progressive_enterobacteriaceae_enrichment == "Yes", na.rm=TRUE) != 18)
  stop("Progressive enrichment YES count must equal 18.")
if (sum(clinical$progressive_enterobacteriaceae_enrichment == "No", na.rm=TRUE) != 26)
  stop("Progressive enrichment NO count must equal 26.")
if (sum(is.na(clinical$progressive_enterobacteriaceae_enrichment)) != 7)
  stop("Progressive enrichment NA count must equal 7.")

write_csv(clinical, file.path(OUT, "01_patient_clinical_mapping_curated.csv"))

# ------------------------------------------------------------
# Locate authoritative Step88A2 displacement file
# ------------------------------------------------------------
disp_candidates <- list.files(
  RESULTS,
  pattern="^V2_STEP88A2_ALL_within_patient_bray_displacement\\.csv$",
  recursive=TRUE,
  full.names=TRUE
)

if (length(disp_candidates) == 0) {
  stop("Could not locate V2_STEP88A2_ALL_within_patient_bray_displacement.csv")
}

# Prefer the most recently modified exact file if duplicates exist.
disp_path <- disp_candidates[which.max(file.info(disp_candidates)$mtime)]
disp <- read_csv(disp_path, show_col_types=FALSE)

req <- c("project","patient_id","time_factor","bray_from_patient_baseline")
miss <- setdiff(req, names(disp))
if (length(miss)) stop("Step88A2 displacement missing columns: ", paste(miss, collapse=", "))

p851 <- disp %>%
  filter(project == "PRJNA851469")

if (nrow(p851) == 0) stop("No PRJNA851469 rows in Step88A2 displacement.")

# Only ICU patient identifiers should be mapped. HealthyVolunteer IDs remain unmapped.
extract_patient_no <- function(x) {
  x <- as.character(x)
  is_patient <- str_detect(tolower(x), "patient")
  out <- rep(NA_integer_, length(x))
  out[is_patient] <- suppressWarnings(
    as.integer(str_extract(x[is_patient], "[0-9]+"))
  )
  out
}

p851_linked <- p851 %>%
  mutate(
    published_patient_no = extract_patient_no(patient_id)
  ) %>%
  left_join(
    clinical,
    by="published_patient_no",
    suffix=c("", "_clinical")
  )

write_csv(
  p851_linked,
  file.path(OUT, "02_PRJNA851469_displacement_with_recovered_clinical.csv")
)

unmatched <- p851_linked %>%
  filter(
    str_detect(tolower(as.character(patient_id)), "patient"),
    is.na(admission_group)
  ) %>%
  distinct(patient_id, published_patient_no, .keep_all=TRUE)

write_csv(
  unmatched,
  file.path(OUT, "02B_unmatched_ICU_patient_ids.csv")
)

# ------------------------------------------------------------
# Day-3 landmark dataset
#
# Important:
# Day-3 ecological displacement is NOT allowed to "predict" events
# that occurred on/before day 3. Those patients are excluded from
# the landmark analysis to avoid temporal leakage.
# ------------------------------------------------------------
day3 <- p851_linked %>%
  filter(time_factor == "Day-3") %>%
  filter(!is.na(admission_group)) %>%
  mutate(
    event_before_or_on_day3 =
      !is.na(composite_event_day) & composite_event_day <= 3,

    landmark_eligible =
      !event_before_or_on_day3,

    subsequent_composite_event =
      landmark_eligible &
      !is.na(composite_event_day) &
      composite_event_day > 3 &
      composite_event_day <= 30,

    followup_days_from_day3 =
      case_when(
        !landmark_eligible ~ NA_real_,
        subsequent_composite_event ~ composite_event_day - 3,
        TRUE ~ 27
      ),

    z_day3_bray =
      as.numeric(scale(bray_from_patient_baseline))
  )

write_csv(
  day3,
  file.path(OUT, "03_day3_landmark_candidate.csv")
)

land <- day3 %>% filter(landmark_eligible)

# ------------------------------------------------------------
# Readiness rules
# ------------------------------------------------------------
n_day3 <- nrow(day3)
n_land <- nrow(land)
n_future_events <- sum(land$subsequent_composite_event, na.rm=TRUE)

abx_counts <- land %>%
  count(antibiotics_at_icu_admission, name="n")

n_abx_yes <- sum(land$antibiotics_at_icu_admission %in% TRUE, na.rm=TRUE)
n_abx_no  <- sum(land$antibiotics_at_icu_admission %in% FALSE, na.rm=TRUE)

outcome_ready <- n_land >= 25 && n_future_events >= 8
abx_ready <- n_land >= 20 && n_abx_yes >= 8 && n_abx_no >= 8

readiness <- tibble(
  analysis=c(
    "DAY3_DISPLACEMENT_VS_SUBSEQUENT_INFECTION_OR_DEATH_LANDMARK",
    "DAY3_DISPLACEMENT_ADJUSTED_FOR_ADMISSION_ANTIBIOTICS",
    "PATIENT_LEVEL_SOFA_ASSOCIATION"
  ),
  ready=c(
    outcome_ready,
    outcome_ready && abx_ready,
    FALSE
  ),
  n_day3=n_day3,
  n_landmark_eligible=n_land,
  n_subsequent_events=n_future_events,
  n_antibiotic_yes=n_abx_yes,
  n_antibiotic_no=n_abx_no,
  reason=c(
    ifelse(outcome_ready,
           "Enough temporally eligible Day-3 patients/events for an exploratory landmark model.",
           "Insufficient linked Day-3 patients and/or subsequent events."),
    ifelse(outcome_ready && abx_ready,
           "Outcome landmark analysis has both antibiotic-at-admission groups represented.",
           "Outcome and/or antibiotic-group counts do not meet conservative readiness threshold."),
    "Individual patient SOFA values are not publicly present in Supplementary Table 1; only aggregate/model summaries are public."
  )
)

write_csv(
  readiness,
  file.path(OUT, "04_analysis_readiness_summary.csv")
)

# ------------------------------------------------------------
# QC summary
# ------------------------------------------------------------
qc <- tibble(
  metric=c(
    "step88a2_source_file",
    "PRJNA851469_displacement_rows",
    "linked_ICU_displacement_rows",
    "linked_ICU_patients",
    "unmatched_ICU_patient_ids",
    "Day3_linked_patients",
    "Day3_events_on_or_before_landmark",
    "Day3_landmark_eligible_patients",
    "Day3_subsequent_composite_events",
    "Day3_admission_antibiotic_yes",
    "Day3_admission_antibiotic_no"
  ),
  value=c(
    disp_path,
    as.character(nrow(p851)),
    as.character(sum(!is.na(p851_linked$admission_group))),
    as.character(n_distinct(p851_linked$published_patient_no[!is.na(p851_linked$admission_group)])),
    as.character(nrow(unmatched)),
    as.character(n_day3),
    as.character(sum(day3$event_before_or_on_day3, na.rm=TRUE)),
    as.character(n_land),
    as.character(n_future_events),
    as.character(n_abx_yes),
    as.character(n_abx_no)
  )
)

write_csv(qc, file.path(OUT, "05_linkage_QC_summary.csv"))

decision_lines <- c(
  "V2 STEP36F PRJNA851469 CLINICAL RECOVERY / LINKAGE",
  "",
  paste0("Published patient mapping rows: ", nrow(clinical)),
  "Cross-checks: infection 28/51; mortality 17/51; antibiotics at ICU admission 28/51; progressive enrichment 18 Yes / 26 No / 7 NA.",
  paste0("Linked ICU patients in frozen displacement table: ",
         n_distinct(p851_linked$published_patient_no[!is.na(p851_linked$admission_group)])),
  paste0("Day-3 linked patients: ", n_day3),
  paste0("Events on/before Day 3 excluded from landmark: ",
         sum(day3$event_before_or_on_day3, na.rm=TRUE)),
  paste0("Day-3 landmark eligible: ", n_land),
  paste0("Subsequent infection/death events after Day 3: ", n_future_events),
  "",
  paste0("Outcome landmark ready: ", outcome_ready),
  paste0("Antibiotic-adjusted landmark ready: ", outcome_ready && abx_ready),
  "",
  "SOFA: NOT RECOVERED at individual-patient level from the public Supplementary Table 1.",
  "",
  "IMPORTANT:",
  "The antibiotic variable used here is only 'antibiotics at ICU admission' (binary), manually curated from the published graphical timeline and cross-checked against Table 1 total 28/51.",
  "Exact antibiotic spectrum/duration is not reconstructed in this step.",
  "No final inferential clinical model was run in Step36F."
)

writeLines(decision_lines, file.path(OUT, "06_DECISION_SUMMARY.txt"))
writeLines(
  c(paste0("Completed: ", Sys.time()), "STEP36F COMPLETE"),
  file.path(OUT, "_STEP36F_COMPLETE.txt")
)

cat("STEP36F COMPLETE\n")
cat("Outcome landmark ready:", outcome_ready, "\n")
cat("Antibiotic-adjusted landmark ready:", outcome_ready && abx_ready, "\n")
