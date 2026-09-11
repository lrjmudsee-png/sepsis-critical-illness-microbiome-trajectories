# ============================================================
# Step 98I: remaining technical validation
# - no FASTQ preprocessing
# - no DADA2 rerun
# - no patient-map modification
# - reads the frozen Step98D object and prior frozen results only
# ============================================================

options(stringsAsFactors = FALSE)
suppressPackageStartupMessages({
  library(tidyr)
  library(tibble)
  library(readr)
  library(purrr)
  library(vegan)
  library(zCompositions)
  library(digest)
  # Load dplyr last so dplyr::select masks MASS::select attached by vegan.
  library(dplyr)
})

set.seed(20260907)

ROOT <- "E:/sepsis_project"
OUT <- file.path(ROOT, "results", "V2_UPGRADE_20260910_REMAINING_VALIDATION")
OBJ <- file.path(ROOT, "data", "PRJNA1125274", "03_dada2",
                 "PRJNA1125274_analysis_object_external_validation.rds")
SEQTAB <- file.path(ROOT, "data", "PRJNA1125274", "03_dada2",
                    "PRJNA1125274_full_seqtab_se.rds")
SAVED <- file.path(ROOT, "results", "V2_UPGRADE_20260910_FINAL_FREEZE",
                   "PRJNA1125274_FINAL_primary_and_sensitivity_results.csv")
ADDENDUM <- file.path(ROOT, "results", "V2_UPGRADE_20260907",
                     "98C_PRJNA1125274_METADATA_FREEZE",
                     "PRJNA1125274_SE_R1_strategy_addendum.md")
SCRIPT98E <- file.path(ROOT, "code", "03_data_processing",
                      "98E_PRJNA1125274_external_validation.R")
RUNTIME98E <- file.path(ROOT, "results", "V2_UPGRADE_20260907",
                       "98E_PRJNA1125274_EXTERNAL_VALIDATION",
                       "_STEP98E_runtime.txt")
EFFECTS <- file.path(ROOT, "results",
                     "V2_36I3_CROSS_COHORT_TAXONOMIC_REPRODUCIBILITY_FINAL",
                     "02_HARMONIZED_WITHIN_COHORT_EFFECTS.csv")
PAIRWISE <- file.path(ROOT, "results",
                      "V2_36I3_CROSS_COHORT_TAXONOMIC_REPRODUCIBILITY_FINAL",
                      "03_PRIMARY3_PAIRWISE_CONCORDANCE.csv")
FIGDATA <- file.path(ROOT, "results",
                     "V2_96E_FIGURE3D_CROSS_COHORT_REPRODUCIBILITY_FINAL",
                     "02_SUPPLEMENTARY_SENSITIVITY_DATA.csv")
THIS_SCRIPT <- file.path(ROOT, "code", "03_data_processing",
                         "98I_remaining_validation_no_fastq.R")
CACHE <- file.path(ROOT, "results", ".98I_remaining_validation_cache_full268")
OBSOLETE_CACHE <- file.path(ROOT, "results", ".98I_remaining_validation_cache")

required <- c(OBJ, SEQTAB, SAVED, ADDENDUM, SCRIPT98E, RUNTIME98E,
              EFFECTS, PAIRWISE, FIGDATA, THIS_SCRIPT)
if (!all(file.exists(required))) {
  stop("Missing required input: ", paste(required[!file.exists(required)], collapse = "; "))
}
if (dir.exists(OUT) && length(list.files(OUT, all.files = FALSE, no.. = TRUE)) > 0L) {
  stop("Output directory already contains files; refusing to overwrite: ", OUT)
}
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
dir.create(CACHE, recursive = TRUE, showWarnings = FALSE)

write_csv_safe <- function(x, path) write_excel_csv(x, path, na = "")
sha256 <- function(path) digest(file = path, algo = "sha256", serialize = FALSE)
J_paired <- function(n) 1 - 3 / (4 * (n - 1) - 1)
DATE_ANOMALY <- "LETTER_ORDER_TIME_CODE_DATE_ANOMALY"

runtime_log <- file.path(OUT, "_STEP98I_runtime.txt")
stamp <- function(x) cat(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), x, "\n",
                         file = runtime_log, append = TRUE)
stamp("START")

# -------------------------------------------------------------------------
# 1. Load and verify frozen analysis object and taxonomy/sequence alignment
# -------------------------------------------------------------------------
obj <- readRDS(OBJ)
seqtab <- readRDS(SEQTAB)
counts <- obj$counts
meta <- as.data.frame(obj$metadata)
tax <- as.data.frame(obj$taxonomy)

stopifnot(nrow(counts) == nrow(meta),
          ncol(counts) == nrow(tax),
          ncol(seqtab) == nrow(tax),
          all(rownames(counts) == as.character(meta$run_accession)),
          all(colnames(counts) == tax$ASV_ID),
          all(colnames(seqtab) == tax$ASV_sequence),
          !anyDuplicated(tax$ASV_ID),
          !anyDuplicated(tax$ASV_sequence))

meta$run_id <- as.character(meta$run_accession)
meta$patient_id <- as.character(meta$patient_id)
meta$depth <- rowSums(counts[meta$run_id, , drop = FALSE])

framework <- meta |>
  filter(timepoint %in% c("T0", "T1", "T2"), !is.na(patient_id), depth >= 2000) |>
  select(run_id, patient_id, hospital, timepoint, alias_letter,
         collection_date, mapping_status, complete_T0_T1_T2, depth)

stopifnot(nrow(framework) == 268L,
          n_distinct(framework$patient_id) == 131L)

scenario_registry <- tribble(
  ~scenario_id, ~taxonomy_policy, ~min_total_count, ~min_prevalence,
  "EXECUTED_ALL_FEATURES_T20_P3", "ALL_FEATURES_AS_EXECUTED", 20L, 3L,
  "DOCUMENTED_ALL_FEATURES_T5_P2", "ALL_FEATURES_DOCUMENTED_ADDENDUM", 5L, 2L,
  "EXCLUDE_EUKARYOTA_T20_P3", "EXCLUDE_EUKARYOTA_RETAIN_UNCLASSIFIED", 20L, 3L,
  "PROKARYOTE_ONLY_T20_P3", "BACTERIA_OR_ARCHAEA_ONLY", 20L, 3L,
  "EXCLUDE_EUKARYOTA_T5_P2", "EXCLUDE_EUKARYOTA_RETAIN_UNCLASSIFIED", 5L, 2L,
  "PROKARYOTE_ONLY_T5_P2", "BACTERIA_OR_ARCHAEA_ONLY", 5L, 2L
)

taxonomy_mask <- function(policy) {
  kingdom <- as.character(tax$Kingdom)
  if (policy %in% c("ALL_FEATURES_AS_EXECUTED", "ALL_FEATURES_DOCUMENTED_ADDENDUM")) {
    rep(TRUE, nrow(tax))
  } else if (policy == "EXCLUDE_EUKARYOTA_RETAIN_UNCLASSIFIED") {
    is.na(kingdom) | kingdom == "" | kingdom != "Eukaryota"
  } else if (policy == "BACTERIA_OR_ARCHAEA_ONLY") {
    !is.na(kingdom) & kingdom %in% c("Bacteria", "Archaea")
  } else {
    stop("Unknown taxonomy policy: ", policy)
  }
}

clr_matrix <- function(x) log(x) - rowMeans(log(x))

summarize_delta <- function(disp, metric, population_id, scenario) {
  d <- disp |> filter(timepoint %in% c("T1", "T2"))
  if (population_id == "PRESPECIFIED_PRIMARY_DATE_CLEAN") {
    anomaly_ids <- unique(d$patient_id[d$mapping_status == DATE_ANOMALY])
    d <- d |> filter(!patient_id %in% anomaly_ids)
  }
  pair <- d |>
    select(patient_id, timepoint, value = all_of(metric)) |>
    distinct(patient_id, timepoint, .keep_all = TRUE) |>
    pivot_wider(names_from = timepoint, values_from = value) |>
    filter(!is.na(T1), !is.na(T2)) |>
    mutate(delta = T2 - T1)
  x <- pair$delta
  n <- length(x)
  m <- mean(x)
  s <- sd(x)
  se <- s / sqrt(n)
  tcrit <- qt(0.975, n - 1L)
  dz <- m / s
  gz <- J_paired(n) * dz
  gz_se <- J_paired(n) * sqrt(1 / n + dz^2 / (2 * n))
  tibble(
    scenario_id = scenario$scenario_id,
    taxonomy_policy = scenario$taxonomy_policy,
    min_total_count = scenario$min_total_count,
    min_prevalence = scenario$min_prevalence,
    population_id = population_id,
    metric = metric,
    n_patients = n,
    n_features = scenario$n_features,
    retained_read_fraction = scenario$retained_read_fraction,
    mean_paired_difference = m,
    sd_paired_difference = s,
    mean_ci_low = m - tcrit * se,
    mean_ci_high = m + tcrit * se,
    cohen_dz = dz,
    hedges_gz = gz,
    hedges_gz_ci_low = gz - tcrit * gz_se,
    hedges_gz_ci_high = gz + tcrit * gz_se,
    paired_t_p_value = 2 * pt(-abs(m / se), n - 1L),
    wilcoxon_p_value = suppressWarnings(wilcox.test(x, mu = 0, exact = FALSE,
                                                    correct = FALSE)$p.value),
    direction = ifelse(m > 0, "INCREASE", ifelse(m < 0, "DECREASE", "NO_CHANGE")),
    sequence_processing_rerun = FALSE,
    patient_mapping_changed = FALSE
  )
}

run_scenario <- function(i) {
  sc <- scenario_registry[i, ]
  cache_file <- file.path(CACHE, paste0(sc$scenario_id, ".rds"))
  if (file.exists(cache_file)) {
    stamp(paste("SCENARIO_CACHE_HIT", sc$scenario_id))
    return(readRDS(cache_file))
  }
  stamp(paste("SCENARIO_START", sc$scenario_id))
  raw <- counts[framework$run_id, , drop = FALSE]
  policy_keep <- taxonomy_mask(sc$taxonomy_policy)
  candidate_all <- raw[, policy_keep, drop = FALSE]
  keep <- colSums(candidate_all) >= sc$min_total_count &
    colSums(candidate_all > 0) >= sc$min_prevalence
  cnt_all <- candidate_all[, keep, drop = FALSE]
  # CZM multiplicative replacement depends on the full composition matrix.
  # Retain all 268 depth-qualified runs exactly as in the archived Step98E.
  cnt <- cnt_all
  if (ncol(cnt) < 2L || any(rowSums(cnt) == 0)) stop("Invalid feature table for ", sc$scenario_id)

  kept_ids <- colnames(cnt_all)
  tx <- tax[match(kept_ids, tax$ASV_ID), , drop = FALSE]
  weights <- colSums(cnt_all)
  kingdom <- as.character(tx$Kingdom)
  euk <- !is.na(kingdom) & kingdom == "Eukaryota"
  unclassified <- is.na(kingdom) | kingdom == ""
  feature_summary <- tibble(
    scenario_id = sc$scenario_id,
    taxonomy_policy = sc$taxonomy_policy,
    min_total_count = sc$min_total_count,
    min_prevalence = sc$min_prevalence,
    n_depth_qualified_samples = nrow(cnt_all),
    n_distance_samples = nrow(cnt),
    n_features = ncol(cnt_all),
    retained_reads = sum(cnt_all),
    retained_read_fraction = sum(cnt_all) / sum(raw),
    min_retained_reads_per_sample = min(rowSums(cnt_all)),
    median_retained_reads_per_sample = median(rowSums(cnt_all)),
    eukaryota_features = sum(euk),
    eukaryota_read_fraction_within_retained = sum(weights[euk]) / sum(weights),
    unclassified_kingdom_features = sum(unclassified),
    unclassified_kingdom_read_fraction_within_retained = sum(weights[unclassified]) / sum(weights),
    min_sequence_length = min(nchar(tx$ASV_sequence)),
    max_sequence_length = max(nchar(tx$ASV_sequence))
  )

  rel <- cnt / rowSums(cnt)
  bray <- as.matrix(vegdist(rel, method = "bray"))
  imp <- zCompositions::cmultRepl(cnt, label = 0, method = "CZM",
                                 output = "p-counts", z.warning = 1,
                                 z.delete = FALSE, suppress.print = TRUE)
  aitch <- as.matrix(dist(clr_matrix(imp), method = "euclidean"))

  base <- framework |>
    filter(timepoint == "T0") |>
    distinct(patient_id, .keep_all = TRUE) |>
    select(patient_id, base_run = run_id)
  disp <- framework |>
    left_join(base, by = "patient_id") |>
    filter(!is.na(base_run)) |>
    rowwise() |>
    mutate(
      bray_from_T0 = bray[run_id, base_run],
      aitchison_from_T0_CZM = aitch[run_id, base_run]
    ) |>
    ungroup()

  sc$n_features <- ncol(cnt_all)
  sc$retained_read_fraction <- sum(cnt_all) / sum(raw)
  results <- bind_rows(
    map_dfr(c("bray_from_T0", "aitchison_from_T0_CZM"),
            ~ summarize_delta(disp, .x, "PRESPECIFIED_PRIMARY_DATE_CLEAN", sc)),
    map_dfr(c("bray_from_T0", "aitchison_from_T0_CZM"),
            ~ summarize_delta(disp, .x, "COMPLETE_CASE_SENSITIVITY", sc))
  )
  stamp(paste("SCENARIO_COMPLETE", sc$scenario_id, "FEATURES", ncol(cnt_all)))
  gc()
  ans <- list(summary = feature_summary, results = results)
  saveRDS(ans, cache_file)
  ans
}

scenario_runs <- lapply(seq_len(nrow(scenario_registry)), run_scenario)
feature_qc <- bind_rows(lapply(scenario_runs, `[[`, "summary"))
sensitivity <- bind_rows(lapply(scenario_runs, `[[`, "results"))

# Exact numerical verification of the executed scenario against frozen 98H.
saved <- read_csv(SAVED, show_col_types = FALSE)
for (pop in c("PRESPECIFIED_PRIMARY_DATE_CLEAN", "COMPLETE_CASE_SENSITIVITY")) {
  for (metric in c("bray_from_T0", "aitchison_from_T0_CZM")) {
    old <- saved |> filter(analysis_id == .env$pop, .data$metric == .env$metric)
    new <- sensitivity |>
      filter(scenario_id == "EXECUTED_ALL_FEATURES_T20_P3",
             population_id == .env$pop, .data$metric == .env$metric)
    if (nrow(old) != 1L || nrow(new) != 1L) stop("Frozen comparison row not unique")
    for (nm in c("mean_paired_difference", "hedges_gz")) {
      if (abs(old[[nm]] - new[[nm]]) > 1e-10) stop("Frozen mismatch: ", pop, " ", metric, " ", nm)
    }
    if (abs(old$paired_t_p_value - new$paired_t_p_value) > 1e-10) {
      stop("Frozen p-value mismatch: ", pop, " ", metric)
    }
  }
}

reference <- sensitivity |>
  filter(scenario_id == "EXECUTED_ALL_FEATURES_T20_P3") |>
  select(population_id, metric,
         reference_mean = mean_paired_difference,
         reference_gz = hedges_gz,
         reference_p = paired_t_p_value)

comparison <- sensitivity |>
  left_join(reference, by = c("population_id", "metric")) |>
  mutate(
    delta_mean_vs_executed = mean_paired_difference - reference_mean,
    delta_gz_vs_executed = hedges_gz - reference_gz,
    direction_preserved = direction == "INCREASE",
    nominal_significance_preserved = paired_t_p_value < 0.05,
    mean_ci_excludes_zero = mean_ci_low > 0 | mean_ci_high < 0,
    result_stable = direction_preserved & nominal_significance_preserved & mean_ci_excludes_zero
  )

write_csv_safe(feature_qc, file.path(OUT, "PRJNA1125274_feature_taxonomy_QC.csv"))
write_csv_safe(sensitivity, file.path(OUT, "PRJNA1125274_non_target_and_threshold_sensitivity.csv"))
write_csv_safe(comparison, file.path(OUT, "PRJNA1125274_sensitivity_comparison_to_executed_primary.csv"))

# Per-sample non-target read fractions under the executed 20/3 filter.
raw_fw <- counts[framework$run_id, , drop = FALSE]
k_exec <- colSums(raw_fw) >= 20 & colSums(raw_fw > 0) >= 3
exec_ids <- colnames(raw_fw)[k_exec]
exec_tax <- tax[match(exec_ids, tax$ASV_ID), , drop = FALSE]
exec_euk <- !is.na(exec_tax$Kingdom) & exec_tax$Kingdom == "Eukaryota"
exec_uncl <- is.na(exec_tax$Kingdom) | exec_tax$Kingdom == ""
sample_non_target <- framework |>
  mutate(
    retained_reads = rowSums(raw_fw[, k_exec, drop = FALSE]),
    eukaryota_reads = rowSums(raw_fw[, k_exec, drop = FALSE][, exec_euk, drop = FALSE]),
    unclassified_kingdom_reads = rowSums(raw_fw[, k_exec, drop = FALSE][, exec_uncl, drop = FALSE]),
    eukaryota_read_fraction = eukaryota_reads / retained_reads,
    unclassified_kingdom_read_fraction = unclassified_kingdom_reads / retained_reads
  )
write_csv_safe(sample_non_target,
               file.path(OUT, "PRJNA1125274_sample_non_target_read_fraction.csv"))

# -------------------------------------------------------------------------
# 2. Complete-case selection audit using all public fields available in obj
# -------------------------------------------------------------------------
mapped <- meta |>
  filter(timepoint %in% c("T0", "T1", "T2"), !is.na(patient_id)) |>
  mutate(depth_pass = depth >= 2000)

patient_selection <- mapped |>
  group_by(patient_id) |>
  summarise(
    hospital = first(hospital),
    n_mapped_runs = n(),
    n_depth_qualified_runs = sum(depth_pass),
    mapped_timepoints = paste(sort(unique(timepoint)), collapse = ";"),
    depth_qualified_timepoints = paste(sort(unique(timepoint[depth_pass])), collapse = ";"),
    has_depth_T0 = any(timepoint == "T0" & depth_pass),
    has_depth_T1 = any(timepoint == "T1" & depth_pass),
    has_depth_T2 = any(timepoint == "T2" & depth_pass),
    t0_depth = ifelse(any(timepoint == "T0" & depth_pass),
                      depth[which(timepoint == "T0" & depth_pass)[1]], NA_real_),
    date_anomaly = any(mapping_status == DATE_ANOMALY),
    .groups = "drop"
  ) |>
  mutate(
    complete_depth_T0_T1_T2 = has_depth_T0 & has_depth_T1 & has_depth_T2,
    in_prespecified_primary = complete_depth_T0_T1_T2 & !date_anomaly,
    selection_group = case_when(
      in_prespecified_primary ~ "PRESPECIFIED_PRIMARY_DATE_CLEAN",
      complete_depth_T0_T1_T2 & date_anomaly ~ "COMPLETE_CASE_DATE_ANOMALY",
      has_depth_T0 ~ "BASELINE_ELIGIBLE_INCOMPLETE",
      TRUE ~ "NO_DEPTH_QUALIFIED_T0"
    )
  ) |>
  arrange(hospital, patient_id)

stopifnot(nrow(patient_selection) == 134L,
          sum(patient_selection$complete_depth_T0_T1_T2) == 30L,
          sum(patient_selection$in_prespecified_primary) == 24L,
          sum(patient_selection$has_depth_T0) == 116L)

write_csv_safe(patient_selection,
               file.path(OUT, "PRJNA1125274_complete_case_selection_audit.csv"))

baseline <- patient_selection |> filter(has_depth_T0) |>
  mutate(complete_group = ifelse(complete_depth_T0_T1_T2,
                                 "COMPLETE_T0_T1_T2", "INCOMPLETE_AFTER_BASELINE"))

desc_numeric <- function(variable) {
  bind_rows(lapply(sort(unique(baseline$complete_group)), function(g) {
    x <- baseline[baseline$complete_group == g, variable, drop = TRUE]
    tibble(record_type = "DESCRIPTIVE", comparison = "COMPLETE_VS_BASELINE_INCOMPLETE",
           variable = variable, level = g, n = sum(!is.na(x)),
           value_1 = median(x, na.rm = TRUE), value_1_label = "median",
           value_2 = quantile(x, 0.25, na.rm = TRUE), value_2_label = "q1",
           value_3 = quantile(x, 0.75, na.rm = TRUE), value_3_label = "q3",
           test = "", statistic = NA_real_, p_value = NA_real_)
  }))
}

hospital_desc <- baseline |>
  count(complete_group, hospital, name = "n") |>
  group_by(complete_group) |>
  mutate(proportion = n / sum(n)) |>
  ungroup() |>
  transmute(record_type = "DESCRIPTIVE", comparison = "COMPLETE_VS_BASELINE_INCOMPLETE",
            variable = "hospital", level = paste(complete_group, hospital, sep = ":"),
            n = n, value_1 = proportion, value_1_label = "proportion",
            value_2 = NA_real_, value_2_label = "", value_3 = NA_real_, value_3_label = "",
            test = "", statistic = NA_real_, p_value = NA_real_)

hospital_tab <- table(baseline$complete_group, baseline$hospital)
hospital_test <- fisher.test(hospital_tab)
depth_test <- wilcox.test(t0_depth ~ complete_group, data = baseline, exact = FALSE)
run_test <- wilcox.test(n_depth_qualified_runs ~ complete_group, data = baseline, exact = FALSE)

tests <- bind_rows(
  tibble(record_type = "TEST", comparison = "COMPLETE_VS_BASELINE_INCOMPLETE",
         variable = "hospital", level = "", n = nrow(baseline), value_1 = NA_real_,
         value_1_label = "", value_2 = NA_real_, value_2_label = "", value_3 = NA_real_,
         value_3_label = "", test = "Fisher exact", statistic = NA_real_,
         p_value = hospital_test$p.value),
  tibble(record_type = "TEST", comparison = "COMPLETE_VS_BASELINE_INCOMPLETE",
         variable = "t0_depth", level = "", n = nrow(baseline), value_1 = NA_real_,
         value_1_label = "", value_2 = NA_real_, value_2_label = "", value_3 = NA_real_,
         value_3_label = "", test = "Wilcoxon rank sum", statistic = unname(depth_test$statistic),
         p_value = depth_test$p.value),
  tibble(record_type = "TEST", comparison = "COMPLETE_VS_BASELINE_INCOMPLETE",
         variable = "n_depth_qualified_runs", level = "", n = nrow(baseline), value_1 = NA_real_,
         value_1_label = "", value_2 = NA_real_, value_2_label = "", value_3 = NA_real_,
         value_3_label = "", test = "Wilcoxon rank sum", statistic = unname(run_test$statistic),
         p_value = run_test$p.value)
)

selection_summary <- bind_rows(desc_numeric("t0_depth"),
                               desc_numeric("n_depth_qualified_runs"),
                               hospital_desc, tests)
write_csv_safe(selection_summary,
               file.path(OUT, "PRJNA1125274_complete_case_selection_summary_and_tests.csv"))

# -------------------------------------------------------------------------
# 3. Supplementary Figure S1: determine whether identical panels are real
# -------------------------------------------------------------------------
effects <- read_csv(EFFECTS, show_col_types = FALSE)
pairwise <- read_csv(PAIRWISE, show_col_types = FALSE)
figdata <- read_csv(FIGDATA, show_col_types = FALSE)
primary_projects <- c("PRJNA691455", "PRJNA851469", "PRJNA516701")
ranks <- c("GENUS", "FAMILY")

effect_vector_checks <- bind_rows(lapply(primary_projects, function(prj) {
  bind_rows(lapply(ranks, function(rk) {
    z <- effects |>
      filter(project == prj, tax_rank == rk, eligible, named_taxon,
             analysis_role %in% c("PRIMARY", "SENSITIVITY")) |>
      group_by(analysis_role, canonical_taxon) |>
      summarise(effect = mean(effect, na.rm = TRUE), .groups = "drop") |>
      pivot_wider(names_from = analysis_role, values_from = effect)
    tibble(
      check_type = "WITHIN_COHORT_EFFECT_VECTOR",
      project_a = prj, project_b = "", rank = rk,
      n_primary = sum(!is.na(z$PRIMARY)), n_sensitivity = sum(!is.na(z$SENSITIVITY)),
      n_shared = sum(!is.na(z$PRIMARY) & !is.na(z$SENSITIVITY)),
      n_primary_only = sum(!is.na(z$PRIMARY) & is.na(z$SENSITIVITY)),
      n_sensitivity_only = sum(is.na(z$PRIMARY) & !is.na(z$SENSITIVITY)),
      max_abs_difference = max(abs(z$PRIMARY - z$SENSITIVITY), na.rm = TRUE),
      exact_match = all(!is.na(z$PRIMARY) & !is.na(z$SENSITIVITY)) &&
        all(z$PRIMARY == z$SENSITIVITY)
    )
  }))
}))

pair_checks <- pairwise |>
  filter(named_only, project_a %in% primary_projects, project_b %in% primary_projects,
         analysis_role %in% c("PRIMARY", "SENSITIVITY")) |>
  select(project_a, project_b, rank, analysis_role, n_shared, spearman_rho,
         sign_agreement, top_jaccard) |>
  pivot_wider(names_from = analysis_role,
              values_from = c(n_shared, spearman_rho, sign_agreement, top_jaccard)) |>
  transmute(
    check_type = "PAIRWISE_SUMMARY",
    project_a, project_b, rank,
    n_primary = n_shared_PRIMARY, n_sensitivity = n_shared_SENSITIVITY,
    n_shared = pmin(n_shared_PRIMARY, n_shared_SENSITIVITY),
    n_primary_only = NA_integer_, n_sensitivity_only = NA_integer_,
    max_abs_difference = pmax(abs(spearman_rho_PRIMARY - spearman_rho_SENSITIVITY),
                              abs(sign_agreement_PRIMARY - sign_agreement_SENSITIVITY),
                              abs(top_jaccard_PRIMARY - top_jaccard_SENSITIVITY)),
    exact_match = n_shared_PRIMARY == n_shared_SENSITIVITY &
      spearman_rho_PRIMARY == spearman_rho_SENSITIVITY &
      sign_agreement_PRIMARY == sign_agreement_SENSITIVITY &
      top_jaccard_PRIMARY == top_jaccard_SENSITIVITY
  )

fig_rows_ok <- nrow(figdata) == 12L &&
  sum(figdata$analysis_role == "Primary") == 6L &&
  sum(figdata$analysis_role == "Common-anchor sensitivity") == 6L

s1_audit <- bind_rows(effect_vector_checks, pair_checks) |>
  mutate(
    figure_data_has_expected_12_rows = fig_rows_ok,
    conclusion = ifelse(exact_match & figure_data_has_expected_12_rows,
                        "IDENTICAL_BY_SOURCE_DATA_NOT_PANEL_REUSE",
                        "DIFFERENCE_OR_SOURCE_PROBLEM_REQUIRES_REVIEW")
  )
write_csv_safe(s1_audit, file.path(OUT, "TAXONOMIC_CONCORDANCE_S1_PANEL_AUDIT.csv"))

# -------------------------------------------------------------------------
# 4. Governance resolution and final validation report
# -------------------------------------------------------------------------
mtime <- function(p) format(file.info(p)$mtime, "%Y-%m-%d %H:%M:%S %z")
runtime_lines <- readLines(RUNTIME98E, warn = FALSE)

gov <- c(
  "# PRJNA1125274 ASV filter governance resolution",
  "",
  "## Finding",
  "",
  "The Step98C strategy addendum states total count >=5 and prevalence >=2 samples, whereas the archived Step98E code actually executed total count >=20 and prevalence >=3 samples.",
  "",
  "## Local chronology evidence",
  "",
  paste0("- Strategy addendum last modified: ", mtime(ADDENDUM)),
  paste0("- Executed Step98E script last modified: ", mtime(SCRIPT98E)),
  paste0("- Archived Step98E runtime log last modified: ", mtime(RUNTIME98E)),
  paste0("- Runtime entries: ", paste(runtime_lines, collapse = " | ")),
  "",
  "The archived code predates the completed archived results, so 20/3 is the identifiable executed rule. The local files do not constitute immutable public preregistration and cannot prove when the analytical preference for 20/3 was formed.",
  "",
  "## Resolution",
  "",
  "- Keep the original addendum unchanged as historical evidence.",
  "- Report 20/3 as the executed primary feature filter.",
  "- Report 5/2 as a documentation-concordant sensitivity analysis.",
  "- Describe this as a documentation/version-control discrepancy, not as a patient-mapping or sequence-processing error.",
  "- Do not state that the 20/3 threshold was publicly preregistered."
)
writeLines(gov, file.path(OUT, "PRJNA1125274_ASV_filter_governance_resolution.md"), useBytes = TRUE)

primary_cmp <- comparison |>
  filter(population_id == "PRESPECIFIED_PRIMARY_DATE_CLEAN") |>
  arrange(metric, scenario_id)
all_stable <- all(primary_cmp$result_stable)
taxonomy_stable <- all(primary_cmp$result_stable[
  primary_cmp$scenario_id %in% c("EXCLUDE_EUKARYOTA_T20_P3", "PROKARYOTE_ONLY_T20_P3")])
threshold_stable <- all(primary_cmp$result_stable[
  primary_cmp$scenario_id == "DOCUMENTED_ALL_FEATURES_T5_P2"])
s1_ok <- all(s1_audit$exact_match) && fig_rows_ok

sel_counts <- patient_selection |> count(selection_group, name = "n")
get_n <- function(g) sel_counts$n[match(g, sel_counts$selection_group)]
euk_exec <- feature_qc |>
  filter(scenario_id == "EXECUTED_ALL_FEATURES_T20_P3")

fmt_result <- function(metric, scenario) {
  z <- primary_cmp |> filter(.data$metric == .env$metric,
                             scenario_id == .env$scenario)
  paste0("mean delta=", sprintf("%.3f", z$mean_paired_difference),
         ", 95% CI ", sprintf("%.3f", z$mean_ci_low), " to ", sprintf("%.3f", z$mean_ci_high),
         ", Hedges g_z=", sprintf("%.3f", z$hedges_gz),
         ", p=", format.pval(z$paired_t_p_value, digits = 3, eps = 0.001))
}

report <- c(
  "# Remaining validation final report",
  "",
  "Scope: downstream validation from the frozen count/taxonomy/metadata object. No FASTQ preprocessing, DADA2, patient mapping, or frozen 98H files were changed.",
  "",
  "## 1. Feature and taxonomy validation",
  "",
  paste0("The executed 20/3 table retained ", euk_exec$n_features, " ASVs. Eukaryota-labelled ASVs comprised ",
         euk_exec$eukaryota_features, " features and ", sprintf("%.2f%%", 100 * euk_exec$eukaryota_read_fraction_within_retained),
         " of reads within the retained table. Sequence-to-ASV-to-taxonomy alignment passed exactly."),
  paste0("Taxonomy exclusion stability: ", ifelse(taxonomy_stable, "PASSED", "FAILED"), "."),
  paste0("- Exclude Eukaryota, Bray: ", fmt_result("bray_from_T0", "EXCLUDE_EUKARYOTA_T20_P3")),
  paste0("- Exclude Eukaryota, Aitchison: ", fmt_result("aitchison_from_T0_CZM", "EXCLUDE_EUKARYOTA_T20_P3")),
  paste0("- Bacteria/Archaea only, Bray: ", fmt_result("bray_from_T0", "PROKARYOTE_ONLY_T20_P3")),
  paste0("- Bacteria/Archaea only, Aitchison: ", fmt_result("aitchison_from_T0_CZM", "PROKARYOTE_ONLY_T20_P3")),
  "",
  "## 2. ASV threshold validation",
  "",
  paste0("The executed rule is 20/3; the historical addendum says 5/2. Documentation-concordant 5/2 sensitivity: ",
         ifelse(threshold_stable, "PASSED", "FAILED"), "."),
  paste0("- 5/2 Bray: ", fmt_result("bray_from_T0", "DOCUMENTED_ALL_FEATURES_T5_P2")),
  paste0("- 5/2 Aitchison: ", fmt_result("aitchison_from_T0_CZM", "DOCUMENTED_ALL_FEATURES_T5_P2")),
  "",
  "## 3. Complete-case selection",
  "",
  paste0("Mapped patient groups=134; depth-qualified T0=116; complete T0/T1/T2=30; prespecified date-clean primary=24; baseline-eligible but incomplete=",
         get_n("BASELINE_ELIGIBLE_INCOMPLETE"), "."),
  paste0("Hospital distribution Fisher p=", format.pval(hospital_test$p.value, digits = 3),
         "; T0 depth Wilcoxon p=", format.pval(depth_test$p.value, digits = 3), "."),
  "Only hospital, run availability, collection date, mapping status, and sequencing depth are available in the analysis object. Clinical selection bias by age, sex, severity, antibiotics, nutrition, comorbidity, and outcome remains untestable from these files.",
  "",
  "## 4. Supplementary Figure S1",
  "",
  paste0("S1 source audit: ", ifelse(s1_ok, "PASSED", "FAILED"), ". Primary and common-anchor panels are numerically identical because the eligible within-cohort effect vectors themselves are identical for all three selected cohorts at genus and family levels. This is not a plotting-panel reuse error."),
  "",
  "## Final decision",
  "",
  paste0("All newly testable numerical sensitivities passed: ", ifelse(all_stable, "YES", "NO"), "."),
  "The external trajectory conclusion may be retained if the manuscript explicitly reports the executed 20/3 rule, the 5/2 sensitivity, the non-target feature sensitivity, centre heterogeneity, and the unresolved 132-versus-134 bookkeeping issue.",
  "The remaining limitations requiring external information are the author-adjudicated 132-person roster and unavailable clinical covariates; neither can be resolved by additional computation on the current files."
)
writeLines(report, file.path(OUT, "REMAINING_VALIDATION_FINAL_REPORT.md"), useBytes = TRUE)

checklist <- tribble(
  ~check_id, ~check, ~status, ~evidence_file, ~remaining_action,
  "V01", "Frozen 98H numerical reproduction under executed 20/3 rule", "PASS", "PRJNA1125274_non_target_and_threshold_sensitivity.csv", "None",
  "V02", "ASV ID, sequence, taxonomy alignment", "PASS", "PRJNA1125274_feature_taxonomy_QC.csv", "None",
  "V03", "Exclude Eukaryota sensitivity", ifelse(taxonomy_stable, "PASS", "FAIL"), "PRJNA1125274_sensitivity_comparison_to_executed_primary.csv", "Report in supplement",
  "V04", "Bacteria/Archaea-only sensitivity", ifelse(taxonomy_stable, "PASS", "FAIL"), "PRJNA1125274_sensitivity_comparison_to_executed_primary.csv", "Report in supplement",
  "V05", "Documented 5/2 versus executed 20/3 threshold", ifelse(threshold_stable, "PASS_WITH_DOCUMENTATION_DISCREPANCY", "FAIL"), "PRJNA1125274_ASV_filter_governance_resolution.md", "Correct manuscript methods; retain historical addendum",
  "V06", "Complete-case selection using available metadata", "PASS_WITH_LIMITATION", "PRJNA1125274_complete_case_selection_summary_and_tests.csv", "State unavailable clinical-covariate limitation",
  "V07", "Supplementary Figure S1 identical panels", ifelse(s1_ok, "PASS_EXPECTED_IDENTICAL_INPUT", "FAIL"), "TAXONOMIC_CONCORDANCE_S1_PANEL_AUDIT.csv", "Explain why panels are identical",
  "V08", "132 versus 134 author-adjudicated clinical roster", "UNRESOLVED_EXTERNAL_INFORMATION", "Prior Step98G reconciliation", "Request original clinical roster or author confirmation",
  "V09", "Cross-centre generalisability", "LIMITATION_NOT_FAILURE", "Prior Step98H hospital sensitivity", "Retain heterogeneity wording; do not claim hospital-level replication",
  "V10", "No raw sequence preprocessing rerun", "CONFIRMED", "_STEP98I_runtime.txt", "Not required for these downstream checks"
)
write_csv_safe(checklist, file.path(OUT, "REMAINING_VALIDATION_CHECKLIST.csv"))

inputs <- tibble(
  path = normalizePath(required, winslash = "/", mustWork = TRUE),
  bytes = file.info(required)$size,
  modified = format(file.info(required)$mtime, "%Y-%m-%d %H:%M:%S %z"),
  sha256 = vapply(required, sha256, character(1))
)
write_csv_safe(inputs, file.path(OUT, "INPUT_SHA256.csv"))

writeLines(capture.output(sessionInfo()), file.path(OUT, "98I_R_sessionInfo.txt"), useBytes = TRUE)
stamp("ANALYSIS_COMPLETE")
stamp("DONE")

manifest_targets <- list.files(OUT, full.names = TRUE, recursive = TRUE)
manifest_targets <- manifest_targets[basename(manifest_targets) != "OUTPUT_SHA256_MANIFEST.csv"]
manifest <- tibble(
  relative_path = substring(normalizePath(manifest_targets, winslash = "/", mustWork = TRUE),
                            nchar(normalizePath(OUT, winslash = "/", mustWork = TRUE)) + 2L),
  bytes = file.info(manifest_targets)$size,
  sha256 = vapply(manifest_targets, sha256, character(1))
) |> arrange(relative_path)
write_csv_safe(manifest, file.path(OUT, "OUTPUT_SHA256_MANIFEST.csv"))
unlink(CACHE, recursive = TRUE, force = TRUE)
if (dir.exists(OBSOLETE_CACHE)) unlink(OBSOLETE_CACHE, recursive = TRUE, force = TRUE)

cat("STEP98I COMPLETE\n")
cat("OUTPUT:", OUT, "\n")
