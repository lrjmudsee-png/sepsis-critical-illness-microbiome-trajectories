# ============================================================
# Sepsis V2 - Step94B2A
# CRA002354 SOURCE-STRATIFIED MODEL DESIGN PREPARATION
#
# SAFE TO RUN WHILE VSEARCH STEP94B1B IS STILL RUNNING.
# This script reads ONLY frozen metadata.
# ============================================================

options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(stringr)
  library(tibble)
})

ROOT <- "E:/sepsis_project"

MASTER <- file.path(
  ROOT,
  "data",
  "_V2_ANALYSIS_READY",
  "00_FREEZE",
  "V2_FINAL_FROZEN_master_metadata_CLINICAL_REPAIRED_20260814_131620.csv"
)

OUT <- file.path(
  ROOT,
  "results",
  "V2_34B2A_CRA002354_SOURCE_MODEL_DESIGN_PREP"
)

dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

norm_chr <- function(x) {
  z <- str_squish(as.character(x))
  z[z %in% c("", "NA", "N/A", "NULL", "None", "none", ".", "-")] <- NA_character_
  z
}

num <- function(x) suppressWarnings(as.numeric(as.character(x)))

m <- read_csv(
  MASTER,
  show_col_types = FALSE,
  guess_max = 50000,
  name_repair = "unique"
)

required <- c("project","run_id","patient_id","time_day","infection_source_standard")
miss <- setdiff(required, names(m))
if (length(miss) > 0) {
  stop(paste0("Missing required columns: ", paste(miss, collapse = "; ")))
}

cra <- m %>%
  filter(toupper(project) == "CRA002354") %>%
  transmute(
    Run_ID = norm_chr(run_id),
    patient_id = norm_chr(patient_id),
    sample_id = if ("sample_id" %in% names(m)) norm_chr(sample_id) else NA_character_,
    time_day = num(time_day),
    time_raw = if ("time_raw" %in% names(m)) norm_chr(time_raw) else NA_character_,
    infection_source_raw = norm_chr(infection_source_standard),
    baseline_sofa = if ("baseline_sofa" %in% names(m)) num(baseline_sofa) else NA_real_,
    baseline_apache_ii = if ("baseline_apache_ii" %in% names(m)) num(baseline_apache_ii) else NA_real_,
    baseline_lactate = if ("baseline_lactate" %in% names(m)) num(baseline_lactate) else NA_real_,
    outcome_28d = if ("outcome_28d" %in% names(m)) norm_chr(outcome_28d) else NA_character_,
    age = if ("age" %in% names(m)) num(age) else NA_real_,
    sex = if ("sex" %in% names(m)) norm_chr(sex) else NA_character_
  ) %>%
  mutate(
    infection_source_group = case_when(
      tolower(infection_source_raw) == "lung" ~ "RESPIRATORY",
      tolower(infection_source_raw) %in% c("abdominal", "intestinal") ~ "ABDOMINAL_GI",
      tolower(infection_source_raw) == "blood" ~ "BLOODSTREAM",
      tolower(infection_source_raw) == "urinary" ~ "URINARY",
      tolower(infection_source_raw) == "surgical_site" ~ "SKIN_SOFT_TISSUE",
      tolower(infection_source_raw) == "other" ~ "OTHER_UNKNOWN",
      TRUE ~ "UNMAPPED"
    ),
    pulmonary_binary = ifelse(
      infection_source_group == "RESPIRATORY",
      "PULMONARY",
      "NONPULMONARY_RECORDED"
    )
  ) %>%
  filter(!is.na(Run_ID), !is.na(patient_id), !is.na(time_day), !is.na(infection_source_raw)) %>%
  distinct(Run_ID, .keep_all = TRUE)

write_csv(cra, file.path(OUT, "01_CRA002354_SOURCE_ANALYSIS_METADATA.csv"))

patient_design <- cra %>%
  arrange(patient_id, time_day) %>%
  group_by(patient_id) %>%
  summarise(
    pulmonary_binary = first(pulmonary_binary),
    infection_source_group = first(infection_source_group),
    infection_source_raw = first(infection_source_raw),
    n_samples = n(),
    first_day = min(time_day),
    last_day = max(time_day),
    followup_span_days = last_day - first_day,
    has_day1 = any(time_day == 1),
    has_day3 = any(time_day == 3),
    has_day5plus = any(time_day >= 5),
    has_day7plus = any(time_day >= 7),
    baseline_sofa = {z <- na.omit(baseline_sofa); if(length(z)) z[1] else NA_real_},
    baseline_apache_ii = {z <- na.omit(baseline_apache_ii); if(length(z)) z[1] else NA_real_},
    baseline_lactate = {z <- na.omit(baseline_lactate); if(length(z)) z[1] else NA_real_},
    outcome_28d = {z <- na.omit(outcome_28d); if(length(z)) z[1] else NA_character_},
    age = {z <- na.omit(age); if(length(z)) z[1] else NA_real_},
    sex = {z <- na.omit(sex); if(length(z)) z[1] else NA_character_},
    .groups = "drop"
  )

write_csv(patient_design, file.path(OUT, "02_PATIENT_LEVEL_SOURCE_DESIGN.csv"))

group_balance <- patient_design %>%
  group_by(pulmonary_binary) %>%
  summarise(
    patients = n(),
    GE2_samples = sum(n_samples >= 2),
    GE3_samples = sum(n_samples >= 3),
    GE4_samples = sum(n_samples >= 4),
    median_samples = median(n_samples),
    median_followup_span = median(followup_span_days),
    day1_patients = sum(has_day1),
    day3_patients = sum(has_day3),
    day5plus_patients = sum(has_day5plus),
    day7plus_patients = sum(has_day7plus),
    .groups = "drop"
  )

write_csv(group_balance, file.path(OUT, "03_GROUP_BALANCE_AND_LONGITUDINAL_DEPTH.csv"))

time_distribution <- cra %>%
  count(pulmonary_binary, time_day, name = "samples") %>%
  arrange(pulmonary_binary, time_day)

write_csv(time_distribution, file.path(OUT, "04_TIME_DISTRIBUTION_BY_SOURCE.csv"))

patient_time_presence <- cra %>%
  distinct(patient_id, pulmonary_binary, time_day) %>%
  count(pulmonary_binary, time_day, name = "patients")

write_csv(patient_time_presence, file.path(OUT, "05_PATIENT_TIMEPOINT_COUNTS_BY_SOURCE.csv"))

covariates <- c("baseline_sofa","baseline_apache_ii","baseline_lactate","outcome_28d","age","sex")

cov_audit <- bind_rows(lapply(covariates, function(v) {
  x <- patient_design[[v]]
  tibble(
    variable = v,
    patients_total = nrow(patient_design),
    nonmissing = sum(!is.na(x)),
    coverage = mean(!is.na(x)),
    unique_nonmissing = n_distinct(x[!is.na(x)])
  )
}))

write_csv(cov_audit, file.path(OUT, "06_BASELINE_COVARIATE_AVAILABILITY.csv"))

n_pulm <- sum(patient_design$pulmonary_binary == "PULMONARY")
n_non <- sum(patient_design$pulmonary_binary == "NONPULMONARY_RECORDED")
n_pulm_ge2 <- sum(patient_design$pulmonary_binary == "PULMONARY" & patient_design$n_samples >= 2)
n_non_ge2 <- sum(patient_design$pulmonary_binary == "NONPULMONARY_RECORDED" & patient_design$n_samples >= 2)

estimands <- tibble(
  estimand = c(
    "Baseline ecological-state difference",
    "Continuous-time longitudinal interaction",
    "Within-patient displacement slope",
    "Early-to-late paired displacement",
    "Alpha-diversity time interaction",
    "Detailed infection-source multicategory"
  ),
  proposed_model = c(
    "Baseline-only Wilcoxon / PERMANOVA with dispersion check",
    "metric ~ time_day * pulmonary_binary + (1|patient_id)",
    "Bray-from-personal-baseline ~ time_day * pulmonary_binary + (1|patient_id)",
    "Within-patient paired early-vs-late effect, compared by source group",
    "alpha ~ time_day * pulmonary_binary + (1|patient_id)",
    "Exploratory/descriptive only"
  ),
  status = c(
    ifelse(n_pulm >= 30 & n_non >= 15, "PRIMARY", "SUPPORTIVE"),
    ifelse(n_pulm_ge2 >= 20 & n_non_ge2 >= 10, "PRIMARY", "SUPPORTIVE"),
    ifelse(n_pulm_ge2 >= 20 & n_non_ge2 >= 10, "PRIMARY", "SUPPORTIVE"),
    "PRIMARY_OR_SENSITIVITY",
    ifelse(n_pulm_ge2 >= 20 & n_non_ge2 >= 10, "PRIMARY", "SUPPORTIVE"),
    "DO_NOT_USE_FOR_FORMAL_INFERENCE"
  )
)

write_csv(estimands, file.path(OUT, "07_PREDEFINED_SOURCE_ANALYSIS_ESTIMANDS.csv"))

usable_cov <- cov_audit %>%
  filter(coverage >= 0.70, unique_nonmissing >= 2) %>%
  pull(variable)

model_plan <- tibble(
  analysis_tier = c("PRIMARY","ADJUSTED_SENSITIVITY","DEPTH_SENSITIVITY","SOURCE_DETAIL_SENSITIVITY"),
  specification = c(
    "time_day * pulmonary_binary + patient random intercept",
    paste0(
      "Primary model + ",
      ifelse(length(usable_cov) == 0,
             "no baseline covariates pass >=70% coverage rule",
             paste(usable_cov, collapse = " + "))
    ),
    "Repeat key results at >=2000 reads if primary uses >=4000 reads",
    "Exclude OTHER_UNKNOWN from non-pulmonary group and repeat key contrasts"
  )
)

write_csv(model_plan, file.path(OUT, "08_PREDEFINED_MODEL_SENSITIVITY_PLAN.csv"))

guardrails <- tibble(
  claim = c(
    "Pulmonary infection causes a different microbiome trajectory",
    "Pulmonary and recorded non-pulmonary sepsis show different longitudinal ecological trajectories",
    "All non-pulmonary infection sources are biologically equivalent",
    "Detailed anatomical source categories can be formally compared",
    "OTHER_UNKNOWN can be treated as a known anatomical source",
    "Within-cohort source-by-time association can be tested"
  ),
  status = c(
    "NO_CAUSAL_CLAIM",
    "TESTABLE_ASSOCIATION",
    "NO",
    "NO",
    "NO",
    "YES"
  )
)

write_csv(guardrails, file.path(OUT, "09_SOURCE_ANALYSIS_CLAIM_GUARDRAILS.csv"))

summary <- tibble(
  patients = nrow(patient_design),
  pulmonary_patients = n_pulm,
  nonpulmonary_recorded_patients = n_non,
  pulmonary_GE2 = n_pulm_ge2,
  nonpulmonary_GE2 = n_non_ge2,
  usable_adjustment_covariates = paste(usable_cov, collapse = ";"),
  metadata_ready_for_formal_binary_source_model =
    (n_pulm >= 30 && n_non >= 15 && n_pulm_ge2 >= 20 && n_non_ge2 >= 10),
  microbiome_abundance_ready = FALSE,
  safe_to_run_during_VSEARCH = TRUE
)

write_csv(summary, file.path(OUT, "10_STEP94B2A_METADATA_DESIGN_READINESS.csv"))

writeLines(
  c(
    "STEP94B2A SOURCE MODEL DESIGN PREP COMPLETE",
    "",
    paste0("Pulmonary patients: ", n_pulm),
    paste0("Non-pulmonary recorded patients: ", n_non),
    paste0("Pulmonary >=2 samples: ", n_pulm_ge2),
    paste0("Non-pulmonary >=2 samples: ", n_non_ge2),
    "",
    "This step is metadata-only and safe to run while VSEARCH is processing."
  ),
  file.path(OUT, "11_STEP94B2A_INTERPRETATION.txt")
)

writeLines(
  c(
    paste0("Completed: ", Sys.time()),
    "STEP94B2A COMPLETE"
  ),
  file.path(OUT, "_STEP94B2A_COMPLETE.ok")
)

cat("STEP94B2A COMPLETE\n")
