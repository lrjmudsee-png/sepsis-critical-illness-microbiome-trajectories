# ============================================================
# PRJNA1125274 final metadata and statistical-governance freeze
# Step 98H. Downstream statistics only: no FASTQ, DADA2, taxonomy,
# chimera removal, 98D processing, or modification of patient identity.
# ============================================================

options(stringsAsFactors = FALSE)
suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(digest)
})

set.seed(20260907)

ROOT <- "E:/sepsis_project"
UPGRADE <- file.path(ROOT, "results", "V2_UPGRADE_20260907")
OUT <- file.path(ROOT, "results", "V2_UPGRADE_20260910_FINAL_FREEZE")
RECON <- file.path(UPGRADE, "98G_PRJNA1125274_PATIENT_RECONCILIATION_20260909")
BASE98E <- file.path(UPGRADE, "98E_PRJNA1125274_EXTERNAL_VALIDATION")
FREEZE98C <- file.path(UPGRADE, "98C_PRJNA1125274_METADATA_FREEZE")

if (dir.exists(OUT) && length(list.files(OUT, all.files = FALSE, no.. = TRUE)) > 0L) {
  stop("Final-freeze output directory already contains files. Refusing to overwrite: ", OUT)
}
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

write_csv_safe <- function(x, path) readr::write_excel_csv(x, path, na = "")
sha256 <- function(path) digest::digest(file = path, algo = "sha256")
J_paired <- function(n) 1 - 3 / (4 * (n - 1) - 1)
DATE_ANOMALY <- "LETTER_ORDER_TIME_CODE_DATE_ANOMALY"
ORPHAN_IDS <- c("SNO_4", "SNO_34")
ORPHAN_RUNS <- c("SRR29453250", "SRR29453247", "SRR29453297")
SAL8_ID <- "SAL_8"
METRICS <- c("bray_from_T0", "aitchison_from_T0_CZM")
METRIC_LABELS <- c(
  bray_from_T0 = "Bray-Curtis distance from patient-specific T0",
  aitchison_from_T0_CZM = "Aitchison distance from patient-specific T0 (CZM)"
)
ZERO_HANDLING <- c(
  bray_from_T0 = "not applicable",
  aitchison_from_T0_CZM = "CZM multiplicative replacement"
)

path_sap <- file.path(FREEZE98C, "PRJNA1125274_VALIDATION_SAP.md")
path_manifest <- file.path(FREEZE98C, "PRJNA1125274_FINAL_289_run_manifest.csv")
path_map <- file.path(FREEZE98C, "PRJNA1125274_FINAL_patient_time_map.csv")
path_qc <- file.path(FREEZE98C, "PRJNA1125274_metadata_and_run_QC_summary.csv")
path_object <- file.path(ROOT, "data", "PRJNA1125274", "03_dada2", "PRJNA1125274_analysis_object_external_validation.rds")
path_original_disp <- file.path(BASE98E, "PRJNA1125274_external_validation_sample_displacement.csv")
path_original_primary <- file.path(BASE98E, "PRJNA1125274_external_validation_primary_contrast.csv")
path_original_hospital <- file.path(BASE98E, "PRJNA1125274_external_validation_primary_by_hospital.csv")
path_orphan_disp <- file.path(RECON, "98E_ORPHAN_EXCLUSION_SENSITIVITY", "PRJNA1125274_external_validation_sample_displacement.csv")
path_orphan_comparison <- file.path(RECON, "PRJNA1125274_orphan_exclusion_sensitivity_comparison.csv")
path_resolution <- file.path(RECON, "PRJNA1125274_132_vs_134_resolution.md")
path_patient_audit <- file.path(RECON, "PRJNA1125274_patient_reconciliation_table.csv")
path_accounting_audit <- file.path(RECON, "PRJNA1125274_reconciled_population_accounting.csv")
path_object_accounting <- file.path(RECON, "analysis_object_run_accounting.csv")
path_prior_hashes <- file.path(RECON, "input_sha256_after.csv")
path_script <- file.path(ROOT, "code", "03_data_processing", "98H_PRJNA1125274_final_governance_freeze.R")

required <- c(path_sap, path_manifest, path_map, path_qc, path_object,
              path_original_disp, path_original_primary, path_original_hospital,
              path_orphan_disp, path_orphan_comparison, path_resolution,
              path_patient_audit, path_accounting_audit, path_object_accounting,
              path_prior_hashes, path_script)
if (!all(file.exists(required))) stop("Missing required frozen input: ", paste(required[!file.exists(required)], collapse = "; "))

prior_hashes <- read_csv(path_prior_hashes, show_col_types = FALSE)
normalise_path_key <- function(path) tolower(gsub("\\\\", "/", path))
verify_prior_source <- function(source_path) {
  expected <- prior_hashes |>
    filter(normalise_path_key(.data$path) == normalise_path_key(source_path)) |>
    pull(sha256)
  if (length(expected) != 1L) stop("No unique prior SHA256 entry for ", basename(source_path))
  if (!identical(tolower(sha256(source_path)), tolower(expected))) stop("Frozen-source SHA256 mismatch for ", basename(source_path))
}
invisible(lapply(c(path_manifest, path_map, path_qc), verify_prior_source))

manifest <- read_csv(path_manifest, show_col_types = FALSE)
patient_map <- read_csv(path_map, show_col_types = FALSE)
original_disp <- read_csv(path_original_disp, show_col_types = FALSE)
orphan_disp <- read_csv(path_orphan_disp, show_col_types = FALSE)
saved_primary <- read_csv(path_original_primary, show_col_types = FALSE)
saved_orphan <- read_csv(path_orphan_comparison, show_col_types = FALSE)
recon_accounting <- read_csv(path_accounting_audit, show_col_types = FALSE)
object_accounting <- read_csv(path_object_accounting, show_col_types = FALSE)

stopifnot(nrow(manifest) == 289L,
          n_distinct(manifest$patient_id) == 134L,
          n_distinct(manifest$patient_id[manifest$timepoint == "T0"]) == 132L,
          all(ORPHAN_IDS %in% manifest$patient_id),
          !any(orphan_disp$patient_id %in% ORPHAN_IDS))

build_pairs <- function(disp, metric, exclude_date_anomaly = FALSE,
                        exclude_patient_ids = character()) {
  long <- disp |>
    filter(timepoint %in% c("T1", "T2")) |>
    select(patient_id, hospital, timepoint, mapping_status, value = all_of(metric))
  dup <- long |> count(patient_id, timepoint, name = "n") |> filter(n > 1L)
  if (nrow(dup) > 0L) stop("Duplicate patient-timepoint found in displacement table.")
  all_pairs <- long |>
    select(patient_id, hospital, timepoint, value) |>
    pivot_wider(names_from = timepoint, values_from = value) |>
    filter(!is.na(T1), !is.na(T2)) |>
    arrange(patient_id)
  anomaly_ids <- unique(long$patient_id[long$mapping_status == DATE_ANOMALY])
  ids <- all_pairs$patient_id
  if (exclude_date_anomaly) ids <- setdiff(ids, anomaly_ids)
  ids <- setdiff(ids, exclude_patient_ids)
  list(
    pairs = all_pairs |>
      filter(patient_id %in% ids) |>
      mutate(paired_difference = T2 - T1) |>
      arrange(patient_id),
    all_pair_ids = all_pairs$patient_id,
    anomaly_ids = anomaly_ids,
    n_date_anomaly_complete_patients = length(intersect(all_pairs$patient_id, anomaly_ids)),
    n_explicit_exclusions = length(intersect(all_pairs$patient_id, exclude_patient_ids))
  )
}

summarize_contrast <- function(disp, metric, analysis_id, analysis_role,
                               population_label, exclude_date_anomaly = FALSE,
                               exclude_patient_ids = character(),
                               orphan_alias_exclusion = FALSE,
                               feature_count = NA_integer_, input_distance_table) {
  p <- build_pairs(disp, metric, exclude_date_anomaly, exclude_patient_ids)
  delta <- p$pairs$paired_difference
  n <- length(delta)
  if (n < 2L) stop("Too few paired patients for ", analysis_id)
  m <- mean(delta)
  s <- sd(delta)
  se <- s / sqrt(n)
  tcrit <- qt(0.975, n - 1L)
  tstat <- m / se
  dz <- m / s
  gz <- J_paired(n) * dz
  gz_se <- J_paired(n) * sqrt(1 / n + dz^2 / (2 * n))
  tibble(
    analysis_id = analysis_id,
    analysis_role = analysis_role,
    population_label = population_label,
    metric = metric,
    metric_label = unname(METRIC_LABELS[[metric]]),
    zero_handling = unname(ZERO_HANDLING[[metric]]),
    n_patients = n,
    n_complete_before_filter = length(p$all_pair_ids),
    n_date_anomaly_complete_patients = p$n_date_anomaly_complete_patients,
    n_date_anomaly_patients_excluded = ifelse(exclude_date_anomaly, p$n_date_anomaly_complete_patients, 0L),
    explicit_patient_exclusions = paste(exclude_patient_ids, collapse = ";"),
    n_explicit_patient_exclusions = p$n_explicit_exclusions,
    orphan_alias_exclusion = orphan_alias_exclusion,
    orphan_aliases_removed = ifelse(orphan_alias_exclusion, paste(ORPHAN_IDS, collapse = ";"), ""),
    orphan_runs_removed = ifelse(orphan_alias_exclusion, paste(ORPHAN_RUNS, collapse = ";"), ""),
    n_features_after_prevalence_filter = feature_count,
    mean_paired_difference = m,
    sd_paired_difference = s,
    mean_paired_difference_se = se,
    mean_paired_difference_ci_low = m - tcrit * se,
    mean_paired_difference_ci_high = m + tcrit * se,
    median_paired_difference = median(delta),
    cohen_dz = dz,
    hedges_gz = gz,
    hedges_gz_se = gz_se,
    hedges_gz_ci_low = gz - tcrit * gz_se,
    hedges_gz_ci_high = gz + tcrit * gz_se,
    paired_t_statistic = tstat,
    paired_t_df = n - 1L,
    paired_t_p_value = 2 * pt(-abs(tstat), n - 1L),
    wilcoxon_p_value = suppressWarnings(wilcox.test(delta, mu = 0, exact = FALSE, correct = FALSE)$p.value),
    direction = ifelse(m > 0, "INCREASE", ifelse(m < 0, "DECREASE", "NO_CHANGE")),
    discovery_concordant_direction = "INCREASE",
    input_distance_table = basename(input_distance_table),
    sequence_processing_rerun = FALSE,
    patient_mapping_changed = FALSE
  )
}

assert_close <- function(actual, expected, label, tol = 1e-10) {
  if (length(actual) != 1L || length(expected) != 1L || is.na(actual) || is.na(expected) || abs(actual - expected) > tol) {
    stop("Frozen-result verification failed: ", label)
  }
}

primary_rows <- bind_rows(lapply(METRICS, function(metric) {
  summarize_contrast(original_disp, metric,
                     "PRESPECIFIED_PRIMARY_DATE_CLEAN", "PRESPECIFIED_PRIMARY",
                     "24 depth-qualified complete cases without date-order anomaly",
                     exclude_date_anomaly = TRUE, feature_count = 9373L,
                     input_distance_table = path_original_disp)
}))
complete_rows <- bind_rows(lapply(METRICS, function(metric) {
  summarize_contrast(original_disp, metric,
                     "COMPLETE_CASE_SENSITIVITY", "COMPLETE_CASE_SENSITIVITY",
                     "30 depth-qualified complete T0/T1/T2 cases",
                     feature_count = 9373L, input_distance_table = path_original_disp)
}))
orphan_primary_rows <- bind_rows(lapply(METRICS, function(metric) {
  summarize_contrast(orphan_disp, metric,
                     "ORPHAN_ALIAS_EXCLUSION_PRESPECIFIED_PRIMARY", "METADATA_ROBUSTNESS_SENSITIVITY",
                     "24 primary cases after removing all 3 orphan-alias runs before feature filtering",
                     exclude_date_anomaly = TRUE, orphan_alias_exclusion = TRUE,
                     feature_count = 9262L, input_distance_table = path_orphan_disp)
}))
orphan_complete_rows <- bind_rows(lapply(METRICS, function(metric) {
  summarize_contrast(orphan_disp, metric,
                     "ORPHAN_ALIAS_EXCLUSION_COMPLETE_CASE", "METADATA_ROBUSTNESS_SENSITIVITY",
                     "30 complete cases after removing all 3 orphan-alias runs before feature filtering",
                     orphan_alias_exclusion = TRUE, feature_count = 9262L,
                     input_distance_table = path_orphan_disp)
}))
sal8_rows <- bind_rows(lapply(METRICS, function(metric) {
  summarize_contrast(original_disp, metric,
                     "SAL8_EXCLUSION_PRESPECIFIED_PRIMARY", "SOURCE_METADATA_CONFLICT_SENSITIVITY",
                     "23 primary cases after exclusion of SAL_8",
                     exclude_date_anomaly = TRUE, exclude_patient_ids = SAL8_ID,
                     feature_count = 9373L, input_distance_table = path_original_disp)
}))

final_results <- bind_rows(primary_rows, complete_rows, orphan_primary_rows,
                           orphan_complete_rows, sal8_rows)

# Reproduce the frozen 98E 30-person and date-clean rows exactly. This checks
# that only governance labels, not numerical endpoints, have changed.
for (metric in METRICS) {
  saved_30 <- saved_primary |> filter(metric == !!metric, sensitivity_set == "ALL_COMPLETE_T0T1T2")
  saved_24 <- saved_primary |> filter(metric == !!metric, sensitivity_set == "EXCLUDING_DATE_ANOMALY")
  new_30 <- complete_rows |> filter(metric == !!metric)
  new_24 <- primary_rows |> filter(metric == !!metric)
  if (nrow(saved_30) != 1L || nrow(saved_24) != 1L) stop("Missing frozen 98E verification row for ", metric)
  assert_close(new_30$mean_paired_difference, saved_30$mean_paired_diff, paste(metric, "complete mean"))
  assert_close(new_30$hedges_gz, saved_30$hedges_gz, paste(metric, "complete Hedges gz"))
  assert_close(new_30$paired_t_p_value, saved_30$paired_t_p, paste(metric, "complete p"))
  assert_close(new_24$mean_paired_difference, saved_24$mean_paired_diff, paste(metric, "primary mean"))
  assert_close(new_24$hedges_gz, saved_24$hedges_gz, paste(metric, "primary Hedges gz"))
  assert_close(new_24$paired_t_p_value, saved_24$paired_t_p, paste(metric, "primary p"))
  saved_orphan_row <- saved_orphan |> filter(metric == !!metric, sensitivity_set == "EXCLUDING_DATE_ANOMALY")
  new_orphan <- orphan_primary_rows |> filter(metric == !!metric)
  if (nrow(saved_orphan_row) != 1L) stop("Missing archived orphan sensitivity for ", metric)
  assert_close(new_orphan$mean_paired_difference, saved_orphan_row$mean_paired_diff_orphan_exclusion, paste(metric, "orphan mean"))
  assert_close(new_orphan$hedges_gz, saved_orphan_row$hedges_gz_orphan_exclusion, paste(metric, "orphan Hedges gz"))
  assert_close(new_orphan$paired_t_p_value, saved_orphan_row$paired_t_p_orphan_exclusion, paste(metric, "orphan p"))
}

write_csv_safe(final_results, file.path(OUT, "PRJNA1125274_FINAL_primary_and_sensitivity_results.csv"))

sal8_comparison <- sal8_rows |>
  select(metric, metric_label, n_patients, mean_paired_difference, sd_paired_difference,
         mean_paired_difference_se, mean_paired_difference_ci_low,
         mean_paired_difference_ci_high, cohen_dz, hedges_gz, hedges_gz_se,
         hedges_gz_ci_low, hedges_gz_ci_high, paired_t_p_value,
         wilcoxon_p_value, direction) |>
  rename_with(~ paste0(.x, "_SAL8_excluded"), -c(metric, metric_label)) |>
  left_join(
    primary_rows |>
      select(metric, n_patients, mean_paired_difference, sd_paired_difference,
             mean_paired_difference_se, mean_paired_difference_ci_low,
             mean_paired_difference_ci_high, cohen_dz, hedges_gz, hedges_gz_se,
             hedges_gz_ci_low, hedges_gz_ci_high, paired_t_p_value,
             wilcoxon_p_value, direction) |>
      rename_with(~ paste0(.x, "_prespecified_primary"), -metric),
    by = "metric"
  ) |>
  mutate(
    excluded_patient_id = SAL8_ID,
    source_conflict = "SAL_8 formal alias/SUBMITTER_ID/geo support SAL; TITLE is 8C-SNO_S17",
    delta_mean_paired_difference = mean_paired_difference_SAL8_excluded - mean_paired_difference_prespecified_primary,
    delta_hedges_gz = hedges_gz_SAL8_excluded - hedges_gz_prespecified_primary,
    delta_paired_t_p_value = paired_t_p_value_SAL8_excluded - paired_t_p_value_prespecified_primary,
    sequence_processing_rerun = FALSE,
    patient_mapping_changed = FALSE
  )
write_csv_safe(sal8_comparison, file.path(OUT, "PRJNA1125274_FINAL_SAL8_exclusion_sensitivity.csv"))

# Centre estimates use the same 24-person prespecified population. The Welch
# comparison is an exploratory difference of patient-level Delta trajectories,
# not a replacement for the overall confirmatory endpoint.
hospital_effects <- bind_rows(lapply(METRICS, function(metric) {
  pairs <- build_pairs(original_disp, metric, exclude_date_anomaly = TRUE)$pairs
  bind_rows(lapply(sort(unique(pairs$hospital)), function(h) {
    delta <- pairs |> filter(hospital == h) |> pull(paired_difference)
    n <- length(delta); m <- mean(delta); s <- sd(delta); se <- s / sqrt(n)
    tcrit <- qt(0.975, n - 1L); tstat <- m / se; dz <- m / s; gz <- J_paired(n) * dz
    gz_se <- J_paired(n) * sqrt(1 / n + dz^2 / (2 * n))
    tibble(
      record_type = "CENTRE_STRATIFIED_EFFECT",
      metric = metric,
      metric_label = unname(METRIC_LABELS[[metric]]),
      hospital = h,
      comparison_hospital = "",
      n_patients = n,
      n_comparison_patients = NA_integer_,
      mean_paired_difference = m,
      sd_paired_difference = s,
      mean_paired_difference_se = se,
      mean_paired_difference_ci_low = m - tcrit * se,
      mean_paired_difference_ci_high = m + tcrit * se,
      cohen_dz = dz,
      hedges_gz = gz,
      hedges_gz_se = gz_se,
      hedges_gz_ci_low = gz - tcrit * gz_se,
      hedges_gz_ci_high = gz + tcrit * gz_se,
      paired_t_p_value = 2 * pt(-abs(tstat), n - 1L),
      wilcoxon_p_value = suppressWarnings(wilcox.test(delta, mu = 0, exact = FALSE, correct = FALSE)$p.value),
      interaction_difference_in_delta = NA_real_,
      interaction_ci_low = NA_real_,
      interaction_ci_high = NA_real_,
      interaction_welch_df = NA_real_,
      interaction_p_value = NA_real_,
      direction = ifelse(m > 0, "INCREASE", ifelse(m < 0, "DECREASE", "NO_CHANGE")),
      analysis_role = "HETEROGENEITY_SENSITIVITY_NOT_CONFIRMATORY",
      interpretation = "Centre-specific estimate; small Novara n makes this estimate imprecise.",
      sequence_processing_rerun = FALSE,
      patient_mapping_changed = FALSE
    )
  }))
}))

hospital_interactions <- bind_rows(lapply(METRICS, function(metric) {
  pairs <- build_pairs(original_disp, metric, exclude_date_anomaly = TRUE)$pairs
  sal <- pairs |> filter(hospital == "ALESSANDRIA") |> pull(paired_difference)
  sno <- pairs |> filter(hospital == "NOVARA") |> pull(paired_difference)
  if (length(sal) < 2L || length(sno) < 2L) stop("Insufficient centre n for interaction proxy.")
  test <- t.test(sal, sno, var.equal = FALSE)
  tibble(
    record_type = "EXPLORATORY_BETWEEN_CENTRE_DIFFERENCE_OF_DELTAS",
    metric = metric,
    metric_label = unname(METRIC_LABELS[[metric]]),
    hospital = "ALESSANDRIA",
    comparison_hospital = "NOVARA",
    n_patients = length(sal),
    n_comparison_patients = length(sno),
    mean_paired_difference = NA_real_, sd_paired_difference = NA_real_,
    mean_paired_difference_se = NA_real_, mean_paired_difference_ci_low = NA_real_,
    mean_paired_difference_ci_high = NA_real_, cohen_dz = NA_real_,
    hedges_gz = NA_real_, hedges_gz_se = NA_real_, hedges_gz_ci_low = NA_real_,
    hedges_gz_ci_high = NA_real_, paired_t_p_value = NA_real_, wilcoxon_p_value = NA_real_,
    interaction_difference_in_delta = mean(sal) - mean(sno),
    interaction_ci_low = unname(test$conf.int[1]),
    interaction_ci_high = unname(test$conf.int[2]),
    interaction_welch_df = unname(test$parameter),
    interaction_p_value = test$p.value,
    direction = ifelse(mean(sal) > mean(sno), "ALESSANDRIA_GREATER_DELTA", "NOVARA_GREATER_DELTA"),
    analysis_role = "EXPLORATORY_HOSPITAL_X_TRAJECTORY_INTERACTION_PROXY",
    interpretation = "Welch comparison of patient-level Delta trajectories; exploratory and not a confirmatory endpoint.",
    sequence_processing_rerun = FALSE,
    patient_mapping_changed = FALSE
  )
}))
hospital_results <- bind_rows(hospital_effects, hospital_interactions)
write_csv_safe(hospital_results, file.path(OUT, "PRJNA1125274_FINAL_hospital_sensitivity.csv"))

# Population flow separates public run/alias counts from the actual paired
# analysis denominators. It replaces potentially ambiguous legacy labels.
public_complete_ids <- patient_map |>
  filter(timepoint %in% c("T0", "T1", "T2")) |>
  distinct(patient_id, timepoint, .keep_all = TRUE) |>
  count(patient_id, name = "n_timepoints") |>
  filter(n_timepoints == 3L) |>
  pull(patient_id)
complete_ids <- build_pairs(original_disp, "bray_from_T0")$all_pair_ids
primary_ids <- build_pairs(original_disp, "bray_from_T0", exclude_date_anomaly = TRUE)$pairs$patient_id
date_anomaly_complete_ids <- setdiff(complete_ids, primary_ids)
depth_record <- recon_accounting |> filter(stage == "DEPTH_GE_2000")
object_runs <- unique(object_accounting$run_id)
sal78 <- manifest |> filter(run_accession == "SRR29453361")
if (nrow(depth_record) != 1L || nrow(sal78) != 1L || sal78$patient_id != "SAL_78") {
  stop("Population-flow controls could not verify depth stage or documented SAL_78 integrity exclusion.")
}

population_accounting <- tribble(
  ~flow_order, ~stage_id, ~stage_label, ~count_unit, ~n_at_stage, ~n_runs_at_stage, ~n_alias_groups_at_stage, ~n_alias_groups_with_T0, ~n_complete_T0_T1_T2_alias_groups, ~definition, ~governance_interpretation,
  1L, "PUBLIC_RUNS", "Publicly released Run records", "runs", 289L, 289L, 134L, 132L, 37L,
  "Current public PRJNA1125274 manifest", "Public-data denominator; it is not a clinical enrollment roster.",
  2L, "HOSPITAL_RESTRICTED_ALIAS_GROUPS", "Hospital-restricted patient alias groups", "alias groups", 134L, 289L, 134L, 132L, 37L,
  "Unique hospital-scoped IDs parsed from formal aliases", "Includes SNO_4 and SNO_34, which lack T0.",
  3L, "ALIAS_GROUPS_WITH_T0", "Alias groups with a T0 baseline", "alias groups", 132L, 286L, 132L, 132L, 37L,
  "At least one public A/T0 record", "Baseline-observed denominator; individual clinical identities remain unresolved.",
  4L, "METADATA_COMPLETE_T0_T1_T2", "Metadata-complete T0/T1/T2 groups before sequencing eligibility", "alias groups", length(public_complete_ids), 3L * length(public_complete_ids), NA_integer_, NA_integer_, length(public_complete_ids),
  "Public formal alias map contains A, B and C", "Pre-depth availability only; it is not an analysed complete-case denominator.",
  5L, "ANALYSIS_OBJECT_RUNS", "Runs in existing 98D analysis object", "runs", length(object_runs), length(object_runs), 134L, 131L, 36L,
  "SRR29453361 (SAL_78 T0) is absent because its R1 FASTQ was not byte-verifiable", "Data-integrity exclusion, not a patient-identity or metadata correction.",
  6L, "DEPTH_QUALIFIED_COMPLETE_CASES", "Sequencing-depth-qualified complete T0/T1/T2 cases", "patients", depth_record$complete_T0_T1_T2_recomputed, depth_record$samples, depth_record$patient_alias_groups, depth_record$patients_with_T0, depth_record$complete_T0_T1_T2_recomputed,
  "Existing 98E framework, depth >= 2,000 target reads", "Complete-case sensitivity population: n = 30.",
  7L, "PRESPECIFIED_PRIMARY_DATE_CLEAN", "Date-order-anomaly-excluded complete cases", "patients", length(primary_ids), 3L * length(primary_ids), length(primary_ids), length(primary_ids), length(primary_ids),
  "Frozen SAP excludes chronology violations from primary paired analyses", "Prespecified primary external-validation population: n = 24."
)
write_csv_safe(population_accounting, file.path(OUT, "PRJNA1125274_FINAL_population_accounting.csv"))

membership <- original_disp |>
  filter(patient_id %in% complete_ids, timepoint %in% c("T0", "T1", "T2")) |>
  distinct(patient_id, hospital, timepoint, mapping_status, .keep_all = TRUE) |>
  select(patient_id, hospital, timepoint, run_id, collection_date, mapping_status, depth) |>
  mutate(
    in_complete_case_sensitivity = TRUE,
    in_prespecified_primary = patient_id %in% primary_ids,
    excluded_from_primary_reason = ifelse(patient_id %in% date_anomaly_complete_ids, "DATE_ORDER_ANOMALY_PER_FROZEN_SAP", ""),
    SAL8_source_conflict_flag = patient_id == SAL8_ID,
    in_SAL8_exclusion_sensitivity = patient_id %in% setdiff(primary_ids, SAL8_ID),
    patient_mapping_changed = FALSE
  ) |>
  arrange(hospital, patient_id, timepoint)
write_csv_safe(membership, file.path(OUT, "PRJNA1125274_FINAL_population_membership.csv"))

issue_registry <- tribble(
  ~issue_id, ~issue, ~affected_records_or_population, ~affects_patient_time_mapping, ~affects_prespecified_primary, ~sensitivity_analysis_performed, ~final_processing_decision, ~final_status, ~evidence_basis,
  "ISSUE_01", "132 versus 134 public metadata discrepancy", "134 hospital-restricted alias groups; 132 groups with T0", "No", "No direct membership change", "Yes: orphan-alias removal before feature filtering", "Keep original map; do not merge identities without an author-level crosswalk.", "UNRESOLVED_LIKELY_METADATA_BOOKKEEPING", "All aliases parse consistently; SNO_4/SNO_34 lack T0; public clinical roster unavailable; published counts are internally inconsistent.",
  "ISSUE_02", "SNO_4 lacks T0", "4B-SNO / 4C-SNO; SRR29453250 / SRR29453247", "No", "No", "Yes: included in 3-run orphan-alias sensitivity", "Retain source alias and exclude from all T0-anchored analyses.", "UNRESOLVED_ALIAS_IDENTITY", "No public clinical crosswalk or unique supported merge candidate.",
  "ISSUE_03", "SNO_34 lacks T0", "34C-SNO; SRR29453297", "No", "No", "Yes: included in 3-run orphan-alias sensitivity", "Retain source alias and exclude from all T0-anchored analyses.", "UNRESOLVED_ALIAS_IDENTITY", "No public clinical crosswalk or unique supported merge candidate.",
  "ISSUE_04", "Date-order anomalies", "12 public alias groups; 6 depth-qualified complete cases", "No", "Yes", "Yes: 24-person primary versus 30-person complete-case sensitivity", "Keep A/B/C mapping; classify 24 date-clean cases as primary and 30 cases as sensitivity.", "GOVERNANCE_RESOLVED", "The frozen pre-results SAP requires chronology-violation exclusion from primary paired trajectories.",
  "ISSUE_05", "SAL_8 / SNO_8 TITLE hospital conflict", "SAL_8 T2 and SNO_8 T2 source TITLE fields", "No", "SAL_8: Yes; SNO_8: No", "Yes: SAL_8 exclusion from prespecified primary", "Retain formal alias/SUBMITTER_ID/geo mapping; do not exchange samples.", "UNRESOLVED_SOURCE_METADATA_CONFLICT", "Formal alias, SUBMITTER_ID and geo agree, but TITLE suffix conflicts; author adjudication is unavailable.",
  "ISSUE_06", "One public Run absent from existing analysis object", "SRR29453361 / SAL_78 T0", "No", "No", "No new sequence processing is permitted or needed", "Retain prior data-integrity exclusion and disclose 289 public versus 288 object runs.", "DATA_INTEGRITY_EXCLUSION_DOCUMENTED", "Manifest documents repeatedly MD5-mismatched R1 FASTQ after re-download attempts.",
  "ISSUE_07", "Original 98E naming inconsistent with frozen SAP", "24 date-clean versus 30 complete cases", "No", "Yes", "Yes: numerical results reverified with corrected governance labels", "Relabel 24 as prespecified primary and 30 as complete-case sensitivity; do not alter results post hoc.", "GOVERNANCE_RESOLVED", "SAP was frozen before 98D/98E and explicitly required date-anomaly exclusion from primary paired analyses.",
  "ISSUE_08", "Legacy 37-complete label can be misleading", "37 pre-depth public metadata-complete groups", "No", "No", "No additional analysis", "Use the explicit population-flow table for all reporting.", "REPORTING_RESOLVED", "37 is metadata completeness before depth/data-integrity eligibility, not the final analysed n."
)
write_csv_safe(issue_registry, file.path(OUT, "PRJNA1125274_FINAL_metadata_issue_registry.csv"))

fmt_num <- function(x, digits = 3) formatC(x, digits = digits, format = "f")
fmt_p <- function(x) ifelse(x < 0.001, "<0.001", formatC(x, digits = 3, format = "f"))
get_result <- function(analysis_id, metric) {
  result <- final_results |> filter(analysis_id == !!analysis_id, metric == !!metric)
  if (nrow(result) != 1L) stop("Expected one result: ", analysis_id, " / ", metric)
  result
}
primary_bray <- get_result("PRESPECIFIED_PRIMARY_DATE_CLEAN", "bray_from_T0")
primary_ait <- get_result("PRESPECIFIED_PRIMARY_DATE_CLEAN", "aitchison_from_T0_CZM")
complete_bray <- get_result("COMPLETE_CASE_SENSITIVITY", "bray_from_T0")
complete_ait <- get_result("COMPLETE_CASE_SENSITIVITY", "aitchison_from_T0_CZM")
orphan_bray <- get_result("ORPHAN_ALIAS_EXCLUSION_PRESPECIFIED_PRIMARY", "bray_from_T0")
orphan_ait <- get_result("ORPHAN_ALIAS_EXCLUSION_PRESPECIFIED_PRIMARY", "aitchison_from_T0_CZM")
sal8_bray <- get_result("SAL8_EXCLUSION_PRESPECIFIED_PRIMARY", "bray_from_T0")
sal8_ait <- get_result("SAL8_EXCLUSION_PRESPECIFIED_PRIMARY", "aitchison_from_T0_CZM")
sal_n <- hospital_effects |> filter(metric == "bray_from_T0", hospital == "ALESSANDRIA") |> pull(n_patients)
sno_n <- hospital_effects |> filter(metric == "bray_from_T0", hospital == "NOVARA") |> pull(n_patients)

resolution_text <- c(
  "# PRJNA1125274 final resolution of the 132 vs 134 discrepancy",
  "",
  "Freeze date: 2026-09-10. This governance freeze does not rerun or alter Step98D sequence processing, the frozen patient-time map, or original Step98E files.",
  "",
  "## Final classification",
  "",
  "**unresolved metadata discrepancy / likely metadata bookkeeping issue.** Public metadata support 134 hospital-restricted alias groups but exactly 132 groups with a T0 baseline. No publicly available clinical enrollment crosswalk can establish whether the two no-T0 groups are additional patients or administratively distinct sample groups. No patient identity correction is justified.",
  "",
  "## Localized discrepancy",
  "",
  "- `SNO_4`: `4B-SNO` (T1, SRR29453250) and `4C-SNO` (T2, SRR29453247); no T0.",
  "- `SNO_34`: `34C-SNO` (T2, SRR29453297); no T0.",
  "- These three records do not enter the 24-person prespecified primary population or the 30-person complete-case sensitivity population because neither alias group has a personal baseline.",
  "- Removing all three records before feature filtering reduced retained ASVs from 9,373 to 9,262 but did not materially change the date-clean Bray or Aitchison paired results.",
  "",
  "## Evidence chain",
  "",
  "1. All 289 formal aliases parse consistently to hospital-scoped SAL/SNO aliases and A/B/C time letters; no Run, BioSample or patient-timepoint duplication was detected.",
  "2. The public manifest has 134 alias groups, 132 with T0 and 37 with A/B/C before sequencing eligibility. The two groups without T0 are exactly SNO_4 and SNO_34.",
  "3. The formal report describes 132 participants, but its reported age groups (85 plus 49) total 134. No public centre-specific roster or supplementary patient list resolves this inconsistency.",
  "4. Candidate merges for SNO_4/SNO_34 are non-unique and depend only on date/number proximity. No merge was applied.",
  "",
  "## Final decision",
  "",
  "The original patient-time map is retained unchanged (`patient_mapping_changed = FALSE`). SNO_4 and SNO_34 remain in public-metadata accounting and remain ineligible for T0-anchored paired analyses. A corrected map must not be generated unless source authors provide an auditable de-identified clinical crosswalk.",
  "",
  "## Separate run-count issue",
  "",
  "The 289-public-run denominator and the existing 288-run analysis object differ because SRR29453361 (SAL_78 T0) had a repeatedly MD5-mismatched R1 FASTQ and was excluded as not byte-verifiable. This is a data-integrity exclusion, not a metadata identity decision."
)
writeLines(resolution_text, file.path(OUT, "PRJNA1125274_FINAL_132_vs_134_resolution.md"), useBytes = TRUE)

conclusion_text <- c(
  "# PRJNA1125274 final external-validation conclusion",
  "",
  "## Statistical governance",
  "",
  "The SAP was frozen before Step98D/98E and states that chronology violations are excluded from the primary paired trajectory analysis. The final prespecified primary population is therefore 24 depth-qualified complete cases without a date-order anomaly. The 30 depth-qualified complete cases are retained as a complete-case sensitivity population. This corrects reporting governance only; it does not change the endpoint, the patient map, sequence data, feature construction, or saved results.",
  "",
  "## Prespecified primary result",
  "",
  paste0("The contrast was Δ = D(T0,T2) − D(T0,T1). In 24 date-clean complete cases, Bray-Curtis distance increased: mean Δ=", fmt_num(primary_bray$mean_paired_difference, 6),
         ", 95% CI ", fmt_num(primary_bray$mean_paired_difference_ci_low, 6), " to ", fmt_num(primary_bray$mean_paired_difference_ci_high, 6),
         "; Hedges g_z=", fmt_num(primary_bray$hedges_gz, 3), " (95% CI ", fmt_num(primary_bray$hedges_gz_ci_low, 3), " to ", fmt_num(primary_bray$hedges_gz_ci_high, 3),
         "); paired t p=", fmt_p(primary_bray$paired_t_p_value), "."),
  paste0("Aitchison distance with CZM zero replacement also increased: mean Δ=", fmt_num(primary_ait$mean_paired_difference, 6),
         ", 95% CI ", fmt_num(primary_ait$mean_paired_difference_ci_low, 6), " to ", fmt_num(primary_ait$mean_paired_difference_ci_high, 6),
         "; Hedges g_z=", fmt_num(primary_ait$hedges_gz, 3), " (95% CI ", fmt_num(primary_ait$hedges_gz_ci_low, 3), " to ", fmt_num(primary_ait$hedges_gz_ci_high, 3),
         "); paired t p=", fmt_p(primary_ait$paired_t_p_value), "."),
  "",
  "## Manuscript-ready Methods text",
  "",
  "For PRJNA1125274, hospital-restricted public aliases were mapped to T0 (A), T1 (B), and T2 (C) according to the prespecified framework. The primary external-validation population comprised depth-qualified complete trajectories with no date-order anomaly, as specified before sequence processing. The patient-level contrast was D(T0,T2) minus D(T0,T1), calculated separately for Bray-Curtis relative-abundance distance and Aitchison distance after CZM zero replacement. All depth-qualified complete trajectories were analysed as a complete-case sensitivity population. Public aliases lacking a personal T0 baseline were not eligible for T0-anchored paired analyses.",
  "",
  "## Manuscript-ready Results text",
  "",
  paste0("In the prespecified date-clean external-validation population (n=24), both Bray-Curtis and Aitchison distances from the patient-specific baseline increased from T1 to T2 (Bray-Curtis mean Δ=", fmt_num(primary_bray$mean_paired_difference, 3),
         ", 95% CI ", fmt_num(primary_bray$mean_paired_difference_ci_low, 3), " to ", fmt_num(primary_bray$mean_paired_difference_ci_high, 3),
         ", paired t p=", fmt_p(primary_bray$paired_t_p_value), "; Aitchison mean Δ=", fmt_num(primary_ait$mean_paired_difference, 3),
         ", 95% CI ", fmt_num(primary_ait$mean_paired_difference_ci_low, 3), " to ", fmt_num(primary_ait$mean_paired_difference_ci_high, 3),
         ", paired t p=", fmt_p(primary_ait$paired_t_p_value), "). Prespecified external replication supported continued ecological displacement from the patient-specific baseline across the clinical trajectory."),
  "",
  "## Sensitivity and heterogeneity interpretation",
  "",
  paste0("The 30-person complete-case sensitivity was directionally concordant (Bray Hedges g_z=", fmt_num(complete_bray$hedges_gz, 3),
         ", p=", fmt_p(complete_bray$paired_t_p_value), "; Aitchison Hedges g_z=", fmt_num(complete_ait$hedges_gz, 3),
         ", p=", fmt_p(complete_ait$paired_t_p_value), ")."),
  paste0("Removal of all three no-T0 orphan-alias runs before feature filtering retained the primary direction (Bray Hedges g_z=", fmt_num(orphan_bray$hedges_gz, 3),
         ", p=", fmt_p(orphan_bray$paired_t_p_value), "; Aitchison Hedges g_z=", fmt_num(orphan_ait$hedges_gz, 3),
         ", p=", fmt_p(orphan_ait$paired_t_p_value), ")."),
  paste0("Excluding SAL_8, the one primary-case alias with an unresolved source TITLE conflict, retained the positive direction (Bray Hedges g_z=", fmt_num(sal8_bray$hedges_gz, 3),
         ", p=", fmt_p(sal8_bray$paired_t_p_value), "; Aitchison Hedges g_z=", fmt_num(sal8_ait$hedges_gz, 3),
         ", p=", fmt_p(sal8_ait$paired_t_p_value), ")."),
  paste0("Centre-stratified estimates are heterogeneity evidence only (Alessandria n=", sal_n, "; Novara n=", sno_n,
         "). They should not be interpreted as confirmation by each hospital because the smaller Novara estimate is imprecise. Recommended wording: ‘External replication was supported at the overall SURVEIL-cohort level, while centre-stratified estimates suggested potential heterogeneity and were imprecise in the smaller Novara subset.’"),
  "",
  "This analysis supports an ecological-distance trajectory in the analysed cohort. It does not establish universal validation, validation at each hospital separately, causality, clinical deterioration, or an author-adjudicated 132-person clinical roster."
)
writeLines(conclusion_text, file.path(OUT, "PRJNA1125274_FINAL_external_validation_conclusion.md"), useBytes = TRUE)

summary_text <- c(
  "# PRJNA1125274 final governance freeze: concise summary",
  "",
  "1. **132 versus 134:** unresolved metadata discrepancy / likely metadata bookkeeping issue. Public aliases yield 134 groups, but 132 have T0; no clinical crosswalk supports a merge or proves two extra patients.",
  "2. **Patient-time mapping:** unchanged. No identity correction, patient merge or sample exchange was made.",
  "3. **SNO_4 and SNO_34:** neither enters a T0-anchored confirmatory analysis. Their three runs remain covered by the full orphan-alias exclusion sensitivity.",
  "4. **Prespecified primary n:** 24 depth-qualified complete cases without date-order anomaly.",
  "5. **Why 24 rather than 30:** the SAP frozen before results requires date-order anomalies be excluded from primary paired trajectories. The 30 complete cases are sensitivity only.",
  paste0("6. **Primary Bray:** mean Δ=", fmt_num(primary_bray$mean_paired_difference, 6), ", Hedges g_z=", fmt_num(primary_bray$hedges_gz, 3),
         " (95% CI ", fmt_num(primary_bray$hedges_gz_ci_low, 3), " to ", fmt_num(primary_bray$hedges_gz_ci_high, 3), "), p=", fmt_p(primary_bray$paired_t_p_value), "."),
  paste0("7. **Primary Aitchison:** mean Δ=", fmt_num(primary_ait$mean_paired_difference, 6), ", Hedges g_z=", fmt_num(primary_ait$hedges_gz, 3),
         " (95% CI ", fmt_num(primary_ait$hedges_gz_ci_low, 3), " to ", fmt_num(primary_ait$hedges_gz_ci_high, 3), "), p=", fmt_p(primary_ait$paired_t_p_value), "."),
  paste0("8. **30-person complete-case sensitivity:** same positive direction (Bray g_z=", fmt_num(complete_bray$hedges_gz, 3),
         "; Aitchison g_z=", fmt_num(complete_ait$hedges_gz, 3), ")."),
  paste0("9. **SAL_8 exclusion:** stable positive direction after excluding SAL_8 (Bray g_z=", fmt_num(sal8_bray$hedges_gz, 3),
         "; Aitchison g_z=", fmt_num(sal8_ait$hedges_gz, 3), ")."),
  paste0("10. **Centre heterogeneity:** Alessandria and Novara are sensitivity-only estimates (n=", sal_n, " and ", sno_n,
         "). Small Novara n makes its estimate imprecise; a non-significant centre estimate cannot establish validation failure."),
  "11. **Could metadata or statistics overturn replication?** No. Date-clean primary, complete-case, full orphan-alias exclusion and SAL_8 exclusion all retain positive Bray and Aitchison trajectories. Unresolved metadata requires disclosure, not speculative remapping.",
  "12. **Freeze decision:** PRJNA1125274 can be formally frozen for the present V2 manuscript. Reopen it only if authors supply a de-identified clinical crosswalk or corrected source metadata, using a new versioned adjudication."
)
writeLines(summary_text, file.path(OUT, "PRJNA1125274_FINAL_two_page_summary.md"), useBytes = TRUE)

# Preserve this exact final-freeze code alongside outputs, then hash all source
# and output evidence. The source folders remain untouched.
dir.create(file.path(OUT, "scripts"), recursive = TRUE, showWarnings = FALSE)
copied_script <- file.path(OUT, "scripts", basename(path_script))
if (!file.copy(path_script, copied_script, overwrite = FALSE)) stop("Could not copy final-governance script into output package.")
if (!identical(sha256(path_script), sha256(copied_script))) stop("Copied final-governance script hash mismatch.")

writeLines(capture.output(sessionInfo()), file.path(OUT, "98H_R_sessionInfo.txt"))
writeLines(c(
  paste0("started and completed: ", Sys.time()),
  "scope: existing metadata/displacement tables and archived orphan sensitivity only",
  "sequence_processing_rerun: FALSE",
  "patient_mapping_changed: FALSE"
), file.path(OUT, "_STEP98H_runtime.txt"))

source_items <- c(
  "Frozen SAP" = path_sap,
  "Frozen 289-run manifest" = path_manifest,
  "Frozen patient-time map" = path_map,
  "Frozen metadata/run QC summary" = path_qc,
  "Existing 98D analysis object" = path_object,
  "Existing 98E displacement" = path_original_disp,
  "Existing 98E primary table" = path_original_primary,
  "Existing 98E hospital table" = path_original_hospital,
  "Archived orphan-exclusion displacement" = path_orphan_disp,
  "Archived orphan-exclusion comparison" = path_orphan_comparison,
  "Prior 132-vs-134 resolution" = path_resolution,
  "Prior patient reconciliation table" = path_patient_audit,
  "Prior population accounting" = path_accounting_audit,
  "Prior source-hash inventory" = path_prior_hashes,
  "Final governance script" = path_script
)
source_register <- tibble(
  record_type = "SOURCE_INPUT",
  item = names(source_items),
  path = unname(source_items),
  role = c(
    "Prespecified governance rule", "Public denominator source", "Source patient-time map",
    "Source metadata QC", "Read-only existing downstream object", "Frozen distance input",
    "Frozen original-result verification", "Frozen original-centre-result verification",
    "Completed 3-run full-exclusion sensitivity input", "Completed 3-run sensitivity verification",
    "Prior reconciliation evidence", "Prior patient-level evidence", "Prior accounting evidence",
    "Prior immutability evidence", "Reproducible final-freeze code"
  ),
  sequence_processing_rerun = FALSE,
  patient_mapping_changed = FALSE
) |>
  rowwise() |>
  mutate(bytes = as.numeric(file.info(path)$size), sha256 = sha256(path),
         verification_status = "READ_ONLY_SOURCE_VERIFIED") |>
  ungroup()

final_paths <- list.files(OUT, recursive = TRUE, full.names = TRUE)
final_paths <- final_paths[!file.info(final_paths)$isdir]
final_paths <- final_paths[basename(final_paths) != "V2_UPGRADE_FINAL_evidence_freeze.csv"]
final_paths <- final_paths[basename(final_paths) != "OUTPUT_SHA256_MANIFEST.csv"]
output_register <- tibble(
  record_type = "FINAL_OUTPUT",
  item = basename(final_paths),
  path = final_paths,
  role = "Final governance-freeze deliverable",
  sequence_processing_rerun = FALSE,
  patient_mapping_changed = FALSE
) |>
  rowwise() |>
  mutate(bytes = as.numeric(file.info(path)$size), sha256 = sha256(path),
         verification_status = "CREATED_AND_HASHED") |>
  ungroup()
write_csv_safe(bind_rows(source_register, output_register),
               file.path(OUT, "V2_UPGRADE_FINAL_evidence_freeze.csv"))

# The package manifest excludes itself by design, avoiding a self-referential
# hash while covering all other final-freeze files (including the evidence CSV).
manifest_paths <- list.files(OUT, recursive = TRUE, full.names = TRUE)
manifest_paths <- manifest_paths[!file.info(manifest_paths)$isdir]
manifest_paths <- manifest_paths[basename(manifest_paths) != "OUTPUT_SHA256_MANIFEST.csv"]
output_manifest <- tibble(
  relative_path = substring(manifest_paths, nchar(OUT) + 2L),
  bytes = as.numeric(file.info(manifest_paths)$size),
  sha256 = vapply(manifest_paths, sha256, character(1))
) |>
  arrange(relative_path)
write_csv_safe(output_manifest, file.path(OUT, "OUTPUT_SHA256_MANIFEST.csv"))

cat("DONE 98H final governance freeze\n")
print(final_results |>
  select(analysis_id, metric, n_patients, mean_paired_difference, hedges_gz,
         paired_t_p_value, direction))
