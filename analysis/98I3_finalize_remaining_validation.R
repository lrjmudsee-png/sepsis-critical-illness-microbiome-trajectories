# Finalize Step98I validation reports from completed downstream results.
# No high-dimensional distances are recomputed here.

options(stringsAsFactors = FALSE)
suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(digest)
})

ROOT <- "E:/sepsis_project"
OUT <- file.path(ROOT, "results", "V2_UPGRADE_20260910_REMAINING_VALIDATION")
SAVED <- file.path(ROOT, "results", "V2_UPGRADE_20260910_FINAL_FREEZE",
                   "PRJNA1125274_FINAL_primary_and_sensitivity_results.csv")
DISP <- file.path(ROOT, "results", "V2_UPGRADE_20260907",
                  "98E_PRJNA1125274_EXTERNAL_VALIDATION",
                  "PRJNA1125274_external_validation_sample_displacement.csv")
SCRIPT98I <- file.path(ROOT, "code", "03_data_processing", "98I_remaining_validation_no_fastq.R")
SCRIPT98I2 <- file.path(ROOT, "code", "03_data_processing", "98I2_CZM_5_2_diagnostics.R")
SCRIPT98I3 <- file.path(ROOT, "code", "03_data_processing", "98I3_finalize_remaining_validation.R")

paths <- c(
  comparison = file.path(OUT, "PRJNA1125274_sensitivity_comparison_to_executed_primary.csv"),
  sensitivity = file.path(OUT, "PRJNA1125274_non_target_and_threshold_sensitivity.csv"),
  feature_qc = file.path(OUT, "PRJNA1125274_feature_taxonomy_QC.csv"),
  selection = file.path(OUT, "PRJNA1125274_complete_case_selection_audit.csv"),
  selection_summary = file.path(OUT, "PRJNA1125274_complete_case_selection_summary_and_tests.csv"),
  s1 = file.path(OUT, "TAXONOMIC_CONCORDANCE_S1_PANEL_AUDIT.csv"),
  patient_diag = file.path(OUT, "PRJNA1125274_CZM_5_2_patient_diagnostics.csv"),
  run_diag = file.path(OUT, "PRJNA1125274_CZM_5_2_nonfinite_run_diagnostics.csv"),
  checklist = file.path(OUT, "REMAINING_VALIDATION_CHECKLIST.csv")
)
required <- c(paths, SAVED, DISP, SCRIPT98I, SCRIPT98I2, SCRIPT98I3)
if (!all(file.exists(required))) stop("Missing finalization input")

write_csv_safe <- function(x, p) write_excel_csv(x, p, na = "")
sha256 <- function(path) digest(file = path, algo = "sha256", serialize = FALSE)
DATE_ANOMALY <- "LETTER_ORDER_TIME_CODE_DATE_ANOMALY"

# Correct stability labels: preserve statistical direction/inference separately
# from whether CZM retained the exact same patient set.
comparison <- read_csv(paths[["comparison"]], show_col_types = FALSE)
reference_n <- comparison |>
  filter(scenario_id == "EXECUTED_ALL_FEATURES_T20_P3") |>
  select(population_id, metric, reference_n_patients = n_patients)
comparison <- comparison |>
  select(-any_of(c("same_n_as_reference", "directional_inference_stable",
                   "strict_same_population_stability", "stability_interpretation",
                   "reference_n_patients"))) |>
  left_join(reference_n, by = c("population_id", "metric")) |>
  mutate(
    same_n_as_reference = n_patients == reference_n_patients,
    directional_inference_stable = direction == "INCREASE" &
      paired_t_p_value < 0.05 & wilcoxon_p_value < 0.05 &
      mean_ci_low > 0,
    strict_same_population_stability = directional_inference_stable & same_n_as_reference,
    stability_interpretation = case_when(
      scenario_id == "EXECUTED_ALL_FEATURES_T20_P3" ~ "REFERENCE_NUMERICALLY_REPRODUCED",
      directional_inference_stable & same_n_as_reference ~ "DIRECTION_EFFECT_AND_INFERENCE_STABLE",
      directional_inference_stable & !same_n_as_reference ~ "DIRECTION_AND_INFERENCE_STABLE_WITH_CZM_ROW_LOSS",
      TRUE ~ "RESULT_NOT_STABLE"
    )
  )
write_csv_safe(comparison, paths[["comparison"]])

# Correct the S1 finding: the displayed pairwise summaries are source-identical,
# although PRJNA516701 has small differences in eligible within-cohort taxa.
s1 <- read_csv(paths[["s1"]], show_col_types = FALSE) |>
  select(-any_of(c("display_relevant", "row_interpretation",
                   "overall_figure_integrity_status", "informativeness"))) |>
  mutate(
    display_relevant = check_type == "PAIRWISE_SUMMARY",
    row_interpretation = case_when(
      check_type == "PAIRWISE_SUMMARY" & exact_match ~ "DISPLAYED_PAIRWISE_SOURCE_METRICS_IDENTICAL",
      check_type == "WITHIN_COHORT_EFFECT_VECTOR" & exact_match ~ "WITHIN_COHORT_EFFECT_VECTOR_IDENTICAL",
      check_type == "WITHIN_COHORT_EFFECT_VECTOR" & !exact_match ~ "WITHIN_COHORT_ELIGIBLE_TAXA_OR_EFFECTS_DIFFER_BUT_DISPLAYED_PAIRWISE_SUMMARY_UNCHANGED",
      TRUE ~ "REVIEW_REQUIRED"
    ),
    overall_figure_integrity_status = ifelse(
      all(exact_match[check_type == "PAIRWISE_SUMMARY"]),
      "PASS_SOURCE_VALUES_IDENTICAL_NOT_PANEL_REUSE", "FAIL_SOURCE_MISMATCH"),
    informativeness = "LIMITED_BECAUSE_DISPLAYED_PRIMARY_AND_SENSITIVITY_SUMMARIES_ARE_IDENTICAL"
  )
write_csv_safe(s1, paths[["s1"]])

# Deterministic numeric reproduction comparison to frozen 98H rows.
sens <- read_csv(paths[["sensitivity"]], show_col_types = FALSE)
saved <- read_csv(SAVED, show_col_types = FALSE)
rerun <- sens |>
  filter(scenario_id == "EXECUTED_ALL_FEATURES_T20_P3") |>
  select(population_id, metric, rerun_n_patients = n_patients,
         rerun_mean = mean_paired_difference, rerun_gz = hedges_gz,
         rerun_p = paired_t_p_value)
original <- saved |>
  filter(analysis_id %in% c("PRESPECIFIED_PRIMARY_DATE_CLEAN", "COMPLETE_CASE_SENSITIVITY")) |>
  select(population_id = analysis_id, metric,
         original_n_patients = n_patients,
         original_mean = mean_paired_difference,
         original_gz = hedges_gz,
         original_p = paired_t_p_value)
repro <- original |>
  inner_join(rerun, by = c("population_id", "metric")) |>
  mutate(
    n_match = original_n_patients == rerun_n_patients,
    abs_diff_mean = abs(original_mean - rerun_mean),
    abs_diff_gz = abs(original_gz - rerun_gz),
    abs_diff_p = abs(original_p - rerun_p),
    numeric_tolerance = 1e-10,
    status = ifelse(n_match & pmax(abs_diff_mean, abs_diff_gz, abs_diff_p) <= numeric_tolerance,
                    "MATCH_WITHIN_1E-10", "MISMATCH")
  )
write_csv_safe(repro, file.path(OUT, "PRJNA1125274_frozen_numeric_reproducibility.csv"))
if (!all(repro$status == "MATCH_WITHIN_1E-10")) stop("Frozen numeric reproduction failed")

# Paired-t distribution check for the archived primary 24-person contrasts.
disp <- read_csv(DISP, show_col_types = FALSE)
anomaly_ids <- unique(disp$patient_id[disp$mapping_status == DATE_ANOMALY])
assumptions <- bind_rows(lapply(c("bray_from_T0", "aitchison_from_T0_CZM"), function(metric) {
  pair <- disp |>
    filter(!patient_id %in% anomaly_ids, timepoint %in% c("T1", "T2")) |>
    select(patient_id, timepoint, value = all_of(metric)) |>
    distinct(patient_id, timepoint, .keep_all = TRUE) |>
    pivot_wider(names_from = timepoint, values_from = value) |>
    filter(!is.na(T1), !is.na(T2)) |>
    mutate(delta = T2 - T1)
  sw <- shapiro.test(pair$delta)
  tibble(metric = metric, n_patients = nrow(pair),
         check = "Shapiro-Wilk on paired differences",
         statistic = unname(sw$statistic), p_value = sw$p.value,
         interpretation = ifelse(sw$p.value >= 0.05,
                                 "NO_STRONG_NORMALITY_DEPARTURE_DETECTED",
                                 "NORMALITY_DEPARTURE_DETECTED_USE_WILCOXON_AS_ROBUSTNESS"))
}))
write_csv_safe(assumptions, file.path(OUT, "PRJNA1125274_primary_assumption_checks.csv"))

# Explicit 11/11 statistical fallacy scan required by validation protocol.
fallacy <- tribble(
  ~fallacy_id, ~fallacy, ~severity, ~finding, ~affected_scope,
  1L, "Simpson's paradox", "CAUTION", "Overall effects are positive, but Novara centre estimates are near zero or negative; this is heterogeneity, not a demonstrated full reversal across every stratum.", "External-validation generalisation",
  2L, "Ecological fallacy", "NOTE", "Primary contrasts are patient-level. Cohort-level results must not be used to infer that every patient deteriorates.", "Individual-level wording",
  3L, "Berkson's paradox", "CAUTION", "Public ICU/sepsis cohorts are selected clinical samples; selection mechanisms cannot be reconstructed from available metadata.", "Transportability",
  4L, "Collider bias", "CAUTION", "Antibiotics, nutrition, severity and outcome are unavailable, so conditioning/collider structures cannot be assessed.", "Clinical covariate interpretation",
  5L, "Base-rate neglect", "NOTE_NOT_APPLICABLE", "No diagnostic sensitivity, specificity, PPV or NPV claim is made.", "Not applicable",
  6L, "Regression to the mean", "NOTE", "Participants were not selected by extreme microbiome displacement, but baseline-anchored change without an equivalent clinical control cannot establish disease-specific progression.", "Trajectory interpretation",
  7L, "Survivorship bias", "CAUTION", "Only 30 of 116 patients with a depth-qualified T0 had complete T0/T1/T2, and 24 entered the date-clean primary analysis.", "Complete-case external validation",
  8L, "Look-elsewhere effect", "CAUTION", "Multiple geometries and sensitivities were examined. These are robustness checks, not independent discoveries; no multiplicity-adjusted discovery claim is warranted.", "Sensitivity analyses",
  9L, "Garden of forking paths", "CAUTION", "The historical addendum says 5/2 while executed code used 20/3. Directional results are robust, but local files are not immutable public preregistration.", "Feature-filter governance",
  10L, "Correlation is not causation", "CAUTION", "All cohorts are observational secondary analyses; ecological displacement cannot be described as caused by sepsis or treatment.", "Causal language",
  11L, "Reverse causality", "CAUTION", "Temporal sampling establishes order but does not separate disease evolution from treatment, nutrition or other time-varying exposures.", "Mechanistic interpretation"
)
write_csv_safe(fallacy, file.path(OUT, "STATISTICAL_VALIDATION_FALLACY_SCAN.csv"))

# Update final checklist with qualified rather than binary conclusions.
checklist <- read_csv(paths[["checklist"]], show_col_types = FALSE) |>
  mutate(
    status = case_when(
      check_id == "V05" ~ "PASS_DIRECTIONAL_INFERENCE_WITH_CZM_ROW_LOSS",
      check_id == "V07" ~ "PASS_SOURCE_VALUES_NOT_PANEL_REUSE",
      TRUE ~ status
    ),
    remaining_action = case_when(
      check_id == "V05" ~ "Report executed 20/3 primary; report 5/2 Aitchison row loss (SAL_4 and SAL_51) and prokaryote-only stability",
      check_id == "V07" ~ "Explain identical displayed summaries and limited incremental visual information",
      TRUE ~ remaining_action
    )
  )
write_csv_safe(checklist, paths[["checklist"]])

feature <- read_csv(paths[["feature_qc"]], show_col_types = FALSE)
selection <- read_csv(paths[["selection"]], show_col_types = FALSE)
selection_tests <- read_csv(paths[["selection_summary"]], show_col_types = FALSE) |>
  filter(record_type == "TEST")
patient_diag <- read_csv(paths[["patient_diag"]], show_col_types = FALSE)
lost <- patient_diag |> filter(!complete_aitchison_pair)

get_result <- function(scenario, metric, pop = "PRESPECIFIED_PRIMARY_DATE_CLEAN") {
  sens |> filter(scenario_id == .env$scenario, .data$metric == .env$metric,
                 population_id == .env$pop)
}
format_result <- function(scenario, metric) {
  z <- get_result(scenario, metric)
  paste0("n=", z$n_patients, ", mean delta=", sprintf("%.3f", z$mean_paired_difference),
         ", 95% CI ", sprintf("%.3f", z$mean_ci_low), " to ", sprintf("%.3f", z$mean_ci_high),
         ", Hedges g_z=", sprintf("%.3f", z$hedges_gz),
         " (95% CI ", sprintf("%.3f", z$hedges_gz_ci_low), " to ", sprintf("%.3f", z$hedges_gz_ci_high), ")",
         ", paired-t p=", format.pval(z$paired_t_p_value, digits = 3, eps = 0.001),
         ", Wilcoxon p=", format.pval(z$wilcoxon_p_value, digits = 3, eps = 0.001))
}

fisher_p <- selection_tests$p_value[selection_tests$variable == "hospital"]
depth_p <- selection_tests$p_value[selection_tests$variable == "t0_depth"]
complete_n <- sum(selection$complete_depth_T0_T1_T2)
baseline_n <- sum(selection$has_depth_T0)
primary_n <- sum(selection$in_prespecified_primary)
euk <- feature |> filter(scenario_id == "EXECUTED_ALL_FEATURES_T20_P3")
max_repro_diff <- max(repro$abs_diff_mean, repro$abs_diff_gz, repro$abs_diff_p)
s1_pass <- all(s1$exact_match[s1$check_type == "PAIRWISE_SUMMARY"])

report <- c(
  "## Material Passport",
  "",
  "- Origin Skill: academic-research-suite / experiment-agent",
  "- Origin Mode: validate",
  paste0("- Origin Date: ", format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z")),
  "- Verification Status: VERIFIED",
  "- Version Label: validation_v1",
  "",
  "# Remaining validation report",
  "",
  "- Source: frozen PRJNA1125274 Step98D object, Step98E displacement output, Step98H results, and frozen cross-cohort taxonomic concordance tables",
  "- Overall Confidence: CAUTION",
  "- Scope: downstream validation only; no FASTQ preprocessing, DADA2, patient mapping, or frozen 98H file was changed",
  "",
  "## Statistical findings",
  "",
  "| Analysis | Result | Confidence |",
  "|---|---|---|",
  paste0("| Executed 20/3 Bray primary | ", format_result("EXECUTED_ALL_FEATURES_T20_P3", "bray_from_T0"), " | CAUTION: medium effect; paired differences are non-normal, but Wilcoxon corroborates |"),
  paste0("| Executed 20/3 Aitchison primary | ", format_result("EXECUTED_ALL_FEATURES_T20_P3", "aitchison_from_T0_CZM"), " | SOLID within analysed 24-person population; medium effect |"),
  paste0("| Exclude Eukaryota Bray | ", format_result("EXCLUDE_EUKARYOTA_T20_P3", "bray_from_T0"), " | CAUTION sensitivity with paired-t/Wilcoxon concordance |"),
  paste0("| Exclude Eukaryota Aitchison | ", format_result("EXCLUDE_EUKARYOTA_T20_P3", "aitchison_from_T0_CZM"), " | SOLID sensitivity |"),
  paste0("| Bacteria/Archaea-only Bray | ", format_result("PROKARYOTE_ONLY_T20_P3", "bray_from_T0"), " | CAUTION sensitivity with paired-t/Wilcoxon concordance |"),
  paste0("| Bacteria/Archaea-only Aitchison | ", format_result("PROKARYOTE_ONLY_T20_P3", "aitchison_from_T0_CZM"), " | SOLID sensitivity |"),
  paste0("| Historical-document 5/2 Aitchison | ", format_result("DOCUMENTED_ALL_FEATURES_T5_P2", "aitchison_from_T0_CZM"), " | CAUTION: one date-clean patient lost numerically |"),
  "",
  "## Technical findings",
  "",
  paste0("- The executed 20/3 table contained ", euk$n_features, " ASVs; ", euk$eukaryota_features,
         " were labelled Eukaryota and accounted for ", sprintf("%.2f%%", 100 * euk$eukaryota_read_fraction_within_retained), " of retained reads."),
  "- ASV IDs, count columns, full-seqtab sequences, and taxonomy rows aligned exactly.",
  "- Removing Eukaryota or retaining only Bacteria/Archaea preserved the positive direction, 95% CI exclusion of zero, paired-t inference, and Wilcoxon inference for both primary metrics.",
  "- Bray paired differences departed from normality (Shapiro-Wilk W=0.8884, p=0.0123); Aitchison paired differences did not show evidence of departure (W=0.9594, p=0.4262). The Bray conclusion is retained because the paired Wilcoxon result is concordant.",
  paste0("- The old 5/2 rule retained 22,031 ASVs. CZM-Aitchison retained 23/24 date-clean primary patients and 28/30 complete cases. The affected complete patients were ", paste(lost$patient_id, collapse = " and "), "."),
  "- SAL_51 T1 and SAL_4 T0 produced non-finite CLR rows after CZM in the very sparse 5/2 table. The Bacteria/Archaea-only 5/2 sensitivity restored 24/24 and 30/30 while retaining positive inference.",
  "- Therefore, 20/3 remains the executed primary rule; 5/2 is a documentation-concordant sensitivity with a disclosed numerical limitation, not a replacement primary analysis.",
  "",
  "## Complete-case selection",
  "",
  paste0("- Depth-qualified T0 patients: ", baseline_n, "; complete T0/T1/T2: ", complete_n,
         "; date-clean primary: ", primary_n, "."),
  paste0("- Complete versus baseline-incomplete hospital distribution: Fisher p=", format.pval(fisher_p, digits = 3),
         "; T0 sequencing depth: Wilcoxon p=", format.pval(depth_p, digits = 3), "."),
  "- Complete-case patients necessarily had more qualifying longitudinal runs; that comparison is structural rather than evidence against selection bias.",
  "- Age, sex, severity, antibiotics, nutrition, comorbidity and outcome are absent from the analysis object, so clinical attrition bias cannot be tested.",
  "",
  "## Supplementary Figure S1",
  "",
  paste0("- Integrity verdict: ", ifelse(s1_pass, "PASS", "FAIL"), ". The six displayed pairwise primary values and six displayed common-anchor values are source-identical; the figure code reads distinct rows and does not reuse a panel."),
  "- PRJNA516701 has small differences in eligible within-cohort taxa/effects, but these do not change the pairwise summary values used by the figure. The identical panels are therefore real but add limited visual information.",
  "",
  "## Reproducibility",
  "",
  "- Method: deterministic numeric rerun from the frozen count/taxonomy/metadata object using the same seed and functions.",
  paste0("- Verdict: REPRODUCIBLE. Four frozen primary/complete metric rows matched in patient count and all compared numeric fields within 1e-10; maximum absolute difference=", format(max_repro_diff, scientific = TRUE), "."),
  "- File-level byte identity is not expected because Step98I uses a new output schema; the frozen Step98H files were not rewritten.",
  "",
  "## Warnings",
  "",
  "| Type | Detail |",
  "|---|---|",
  "| Survivorship/attrition | Only 30/116 depth-qualified-baseline patients had complete trajectories; 24/116 entered the primary analysis. |",
  "| Centre heterogeneity | Overall validation must not be restated as successful replication in each hospital. |",
  "| Metadata reconciliation | The public 134-group versus paper 132-person discrepancy still requires an author-adjudicated roster. |",
  "| Governance | Local timestamps identify 20/3 as the executed rule, but do not constitute immutable public preregistration. |",
  "| Multiplicity | Sensitivity p-values are robustness diagnostics, not independent discovery tests. |",
  "| Distributional assumption | Bray paired differences were non-normal; interpret the paired t-test together with the concordant Wilcoxon sensitivity. |",
  "",
  "## Fallacy scan",
  "",
  "- Coverage: 11/11 statistical fallacy types checked.",
  "- Detailed findings: STATISTICAL_VALIDATION_FALLACY_SCAN.csv.",
  "",
  "## Final verdict",
  "",
  "The main PRJNA1125274 ecological-trajectory conclusion survives every same-population non-target-feature sensitivity and the prokaryote-only 5/2 sensitivity. Confidence remains CAUTION rather than unrestricted SOLID because complete-case attrition, centre heterogeneity, unavailable clinical covariates, non-public threshold governance, and the unresolved 132-versus-134 roster issue limit generalisation.",
  "",
  "The remaining unresolved items cannot be answered by further computation on the current files: an author-adjudicated 132-person clinical roster and clinical covariates for attrition/confounding assessment."
)
writeLines(report, file.path(OUT, "REMAINING_VALIDATION_FINAL_REPORT.md"), useBytes = TRUE)

code_paths <- c(SCRIPT98I, SCRIPT98I2, SCRIPT98I3)
code_manifest <- tibble(
  file = basename(code_paths),
  bytes = file.info(code_paths)$size,
  sha256 = vapply(code_paths, sha256, character(1))
)
write_csv_safe(code_manifest, file.path(OUT, "VALIDATION_CODE_SHA256.csv"))

# Rebuild the output manifest after every final report/table update.
manifest_path <- file.path(OUT, "OUTPUT_SHA256_MANIFEST.csv")
targets <- list.files(OUT, full.names = TRUE, recursive = TRUE)
targets <- targets[normalizePath(targets, winslash = "/", mustWork = TRUE) !=
                     normalizePath(manifest_path, winslash = "/", mustWork = FALSE)]
out_norm <- normalizePath(OUT, winslash = "/", mustWork = TRUE)
manifest <- tibble(
  relative_path = substring(normalizePath(targets, winslash = "/", mustWork = TRUE),
                            nchar(out_norm) + 2L),
  bytes = file.info(targets)$size,
  sha256 = vapply(targets, sha256, character(1))
) |>
  arrange(relative_path)
write_csv_safe(manifest, manifest_path)

cat("STEP98I FINALIZATION COMPLETE\n")
