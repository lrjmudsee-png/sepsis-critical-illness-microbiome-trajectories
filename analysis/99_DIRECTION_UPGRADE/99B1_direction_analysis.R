#!/usr/bin/env Rscript

# Stage-B direction analysis. Reads frozen stage-A tables and existing analysis
# objects/results. Does not rerun sequence processing or modify any prior file.

options(stringsAsFactors = FALSE, warn = 1)

args <- commandArgs(trailingOnly = TRUE)
arg_value <- function(flag) {
  idx <- match(flag, args)
  if (is.na(idx) || idx == length(args)) stop("Missing argument: ", flag)
  args[[idx + 1L]]
}
out_dir <- normalizePath(arg_value("--output"), winslash = "/", mustWork = TRUE)
dir.create(file.path(out_dir, "figures"), showWarnings = FALSE, recursive = TRUE)

read_csv <- function(path) {
  read.csv(path, check.names = FALSE, stringsAsFactors = FALSE, fileEncoding = "UTF-8-BOM")
}
write_csv_once <- function(x, name) {
  path <- file.path(out_dir, name)
  if (file.exists(path)) stop("Refusing to overwrite: ", path)
  if (!is.data.frame(x) || nrow(x) == 0L) stop("No rows for required output: ", name)
  write.csv(x, path, row.names = FALSE, na = "", fileEncoding = "UTF-8")
}
write_lines_once <- function(x, name) {
  path <- file.path(out_dir, name)
  if (file.exists(path)) stop("Refusing to overwrite: ", path)
  writeLines(x, path, useBytes = TRUE)
}
as_flag <- function(x) {
  if (is.logical(x)) return(!is.na(x) & x)
  toupper(trimws(as.character(x))) %in% c("TRUE", "T", "1", "YES", "Y")
}
num <- function(x) suppressWarnings(as.numeric(x))
safe_value <- function(x, digits = 4L) {
  ifelse(is.finite(x), formatC(x, digits = digits, format = "f"), "NA")
}
seed_from_label <- function(label) {
  20260913L + (sum(utf8ToInt(label)) %% 100000L)
}

A_DIR <- "E:/sepsis_project/results/V2_UPGRADE_20260913_DIRECTION/A_FEASIBILITY"
a04_path <- file.path(A_DIR, "A04_patient_time_join_audit.csv")
a05_path <- file.path(A_DIR, "A05_pair_eligibility.csv")
a07_path <- file.path(A_DIR, "A07_taxonomy_mapping_proposal.csv")
a08_path <- file.path(A_DIR, "A08_feature_and_zero_coverage.csv")
sap_path <- "E:/sepsis_project/code/03_data_processing/99_DIRECTION_UPGRADE/B00_FIXED_SECONDARY_ANALYSIS_SPECIFICATION.md"
external_population_path <- "E:/sepsis_project/results/V2_UPGRADE_20260910_FINAL_FREEZE/PRJNA1125274_FINAL_population_membership.csv"
prjna851_object_path <- "E:/sepsis_project/data/_V2_ANALYSIS_READY/01_PROJECTS/PRJNA851469/05_work/step87B_analysis_object/PRJNA851469_analysis_object_step87B.rds"

required <- c(a04_path, a05_path, a07_path, a08_path, sap_path, external_population_path, prjna851_object_path)
if (any(!file.exists(required))) stop("Missing required stage-B input: ", paste(required[!file.exists(required)], collapse = " | "))

samples <- read_csv(a04_path)
pairs <- read_csv(a05_path)
taxonomy_map <- read_csv(a07_path)
coverage <- read_csv(a08_path)
external_population <- read_csv(external_population_path)

logical_columns <- c(
  "primary_patient_time_include", "run_present_exactly_in_counts",
  "metadata_run_id_duplicated", "counts_run_id_duplicated",
  "patient_in_frozen_pair_population", "in_external_complete_case_sensitivity",
  "in_external_prespecified_primary", "P_zero", "C_zero", "P_plus_C_double_zero"
)
for (column in intersect(logical_columns, names(samples))) samples[[column]] <- as_flag(samples[[column]])
pair_logical_columns <- c(
  "anchor_input_available", "early_late_inputs_available", "object_inputs_available",
  "in_frozen_primary_pair_population", "in_frozen_sensitivity_pair_population",
  "A_stage_eligible", "any_P_or_C_zero_at_early_or_late",
  "any_P_plus_C_double_zero_at_early_or_late"
)
for (column in intersect(pair_logical_columns, names(pairs))) pairs[[column]] <- as_flag(pairs[[column]])
for (column in c("total_reads", "P_reads", "C_reads")) samples[[column]] <- num(samples[[column]])

balance_value <- function(p, c, pseudocount) {
  result <- rep(NA_real_, length(p))
  valid <- is.finite(p) & is.finite(c) & (p + c > 0)
  result[valid] <- log((p[valid] + pseudocount) / (c[valid] + pseudocount))
  result
}
samples$B_pc0_5 <- balance_value(samples$P_reads, samples$C_reads, 0.5)
samples$B_pc0_1 <- balance_value(samples$P_reads, samples$C_reads, 0.1)
samples$B_pc1 <- balance_value(samples$P_reads, samples$C_reads, 1)
samples$B_no_pseudocount <- ifelse(
  samples$P_reads > 0 & samples$C_reads > 0,
  log(samples$P_reads / samples$C_reads),
  NA_real_
)
samples$P_relative_abundance <- samples$P_reads / samples$total_reads
samples$C_relative_abundance <- samples$C_reads / samples$total_reads
samples$OTHER_reads <- samples$total_reads - samples$P_reads - samples$C_reads
samples$OTHER_relative_abundance <- samples$OTHER_reads / samples$total_reads

if (any(duplicated(paste(samples$project, samples$run_id, sep = "::")))) {
  stop("Sample balance input contains duplicated project+run_id keys")
}

lookup_sample_value <- function(project, run_id, column) {
  idx <- match(paste(project, run_id, sep = "::"), paste(samples$project, samples$run_id, sep = "::"))
  result <- rep(NA_real_, length(idx))
  valid <- !is.na(idx)
  result[valid] <- samples[[column]][idx[valid]]
  result
}
for (column in c("B_pc0_5", "B_pc0_1", "B_pc1", "B_no_pseudocount")) {
  early_name <- paste0("early_", column)
  late_name <- paste0("late_", column)
  delta_name <- paste0("delta_", column)
  pairs[[early_name]] <- lookup_sample_value(pairs$project, pairs$early_run_id, column)
  pairs[[late_name]] <- lookup_sample_value(pairs$project, pairs$late_run_id, column)
  pairs[[delta_name]] <- pairs[[late_name]] - pairs[[early_name]]
}
pairs$main_balance_exclusion_reason <- ifelse(
  !pairs$early_late_inputs_available,
  "missing_or_ambiguous_early_late_input",
  ifelse(
    !is.finite(pairs$delta_B_pc0_5),
    "P_plus_C_double_zero_at_early_or_late",
    ""
  )
)

# Exact join to frozen Bray displacement results.
displacement_paths <- c(
  PRJNA691455 = "E:/sepsis_project/results/V2_UPGRADE_20260907/98A_AITCHISON_LONGITUDINAL_ROBUSTNESS/01_SAMPLE_DISPLACEMENT/PRJNA691455_within_patient_displacement_bray_aitchison.csv",
  PRJNA851469 = "E:/sepsis_project/results/V2_UPGRADE_20260907/98A_AITCHISON_LONGITUDINAL_ROBUSTNESS/01_SAMPLE_DISPLACEMENT/PRJNA851469_within_patient_displacement_bray_aitchison.csv",
  PRJNA516701 = "E:/sepsis_project/results/V2_UPGRADE_20260907/98A_AITCHISON_LONGITUDINAL_ROBUSTNESS/01_SAMPLE_DISPLACEMENT/PRJNA516701_within_patient_displacement_bray_aitchison.csv",
  PRJNA1125274 = "E:/sepsis_project/results/V2_UPGRADE_20260907/98E_PRJNA1125274_EXTERNAL_VALIDATION/PRJNA1125274_external_validation_sample_displacement.csv"
)
displacement_rows <- list()
for (project in names(displacement_paths)) {
  table <- read_csv(displacement_paths[[project]])
  if (project == "PRJNA1125274") {
    table <- data.frame(
      project = project,
      run_id = table$run_id,
      patient_id = table$patient_id,
      displacement = num(table$bray_from_T0),
      stringsAsFactors = FALSE
    )
  } else {
    table <- data.frame(
      project = project,
      run_id = table$run_id,
      patient_id = table$patient_id,
      displacement = num(table$bray_from_patient_baseline),
      stringsAsFactors = FALSE
    )
  }
  displacement_rows[[project]] <- table
}
displacement <- do.call(rbind, displacement_rows)
if (any(duplicated(paste(displacement$project, displacement$run_id, sep = "::")))) {
  stop("Displacement source contains duplicated project+run_id keys")
}
displacement_lookup <- function(project, run_id) {
  idx <- match(paste(project, run_id, sep = "::"), paste(displacement$project, displacement$run_id, sep = "::"))
  result <- rep(NA_real_, length(idx))
  result[!is.na(idx)] <- displacement$displacement[idx[!is.na(idx)]]
  result
}
pairs$early_D <- displacement_lookup(pairs$project, pairs$early_run_id)
pairs$late_D <- displacement_lookup(pairs$project, pairs$late_run_id)
pairs$delta_D <- pairs$late_D - pairs$early_D

# Healthy-reference distance within the frozen PRJNA851469 feature space.
obj851 <- readRDS(prjna851_object_path)
counts851 <- obj851$counts
meta851 <- obj851$metadata
run851 <- as.character(meta851$run_id)
if (!identical(run851, rownames(counts851))) {
  idx <- match(run851, rownames(counts851))
  if (anyNA(idx)) stop("PRJNA851469 metadata cannot be strictly ordered to counts")
  counts851 <- counts851[idx, , drop = FALSE]
}
candidate_columns <- intersect(
  c("source_sample_id", "sample_id", "patient_id", "phenotype", "control_type", "sample_class", "analysis_role"),
  names(meta851)
)
metadata_text <- apply(meta851[, candidate_columns, drop = FALSE], 1L, function(row) paste(row, collapse = "|"))
healthy_flag <- grepl("healthy|volunteer|non[-_ ]?icu|control", metadata_text, ignore.case = TRUE)
if (sum(healthy_flag) != 13L) stop("Expected 13 frozen healthy PRJNA851469 samples; observed ", sum(healthy_flag))
relative851 <- counts851 / rowSums(counts851)
healthy_rel <- relative851[healthy_flag, , drop = FALSE]
patient_rel <- relative851[!healthy_flag, , drop = FALSE]
bray_to_health <- matrix(
  NA_real_, nrow = nrow(patient_rel), ncol = nrow(healthy_rel),
  dimnames = list(rownames(patient_rel), rownames(healthy_rel))
)
for (j in seq_len(nrow(healthy_rel))) {
  reference <- healthy_rel[j, ]
  numerator <- rowSums(abs(sweep(patient_rel, 2L, reference, "-")))
  denominator <- rowSums(sweep(patient_rel, 2L, reference, "+"))
  bray_to_health[, j] <- numerator / denominator
}
healthy_sample_metrics <- data.frame(
  project = "PRJNA851469",
  run_id = rownames(patient_rel),
  patient_id = as.character(meta851$patient_id[match(rownames(patient_rel), run851)]),
  source_sample_id = as.character(meta851$source_sample_id[match(rownames(patient_rel), run851)]),
  time_raw = as.character(meta851$time_raw[match(rownames(patient_rel), run851)]),
  n_healthy_references = nrow(healthy_rel),
  mean_bray_to_healthy = rowMeans(bray_to_health),
  median_bray_to_healthy = apply(bray_to_health, 1L, median),
  min_bray_to_healthy = apply(bray_to_health, 1L, min),
  max_bray_to_healthy = apply(bray_to_health, 1L, max),
  stringsAsFactors = FALSE
)
h_lookup <- setNames(healthy_sample_metrics$mean_bray_to_healthy, healthy_sample_metrics$run_id)
pairs$early_H <- ifelse(
  pairs$project == "PRJNA851469" & pairs$early_run_id %in% names(h_lookup),
  h_lookup[pairs$early_run_id], NA_real_
)
pairs$late_H <- ifelse(
  pairs$project == "PRJNA851469" & pairs$late_run_id %in% names(h_lookup),
  h_lookup[pairs$late_run_id], NA_real_
)
pairs$delta_H <- pairs$late_H - pairs$early_H

effect_summary <- function(values, project, analysis, population, metric, pseudocount = NA_real_, n_expected = NA_integer_) {
  values <- as.numeric(values)
  values <- values[is.finite(values)]
  n <- length(values)
  result <- data.frame(
    project = project, analysis = analysis, population = population, metric = metric,
    pseudocount = pseudocount, n_expected = n_expected, n_analyzed = n,
    n_missing = if (is.finite(n_expected)) n_expected - n else NA_integer_,
    mean = NA_real_, sd = NA_real_, se = NA_real_, t_df = NA_real_,
    t_statistic = NA_real_, t_p = NA_real_, ci_low = NA_real_, ci_high = NA_real_,
    median = NA_real_, wilcoxon_p = NA_real_, cohen_dz = NA_real_,
    bootstrap_ci_low = NA_real_, bootstrap_ci_high = NA_real_, shapiro_p = NA_real_,
    positive_fraction = NA_real_, direction = "NOT_ESTIMABLE", status = "NOT_ESTIMABLE",
    stringsAsFactors = FALSE
  )
  if (n == 0L) return(result)
  result$mean <- mean(values)
  result$median <- median(values)
  result$positive_fraction <- mean(values > 0)
  result$direction <- ifelse(result$mean > 0, "POSITIVE", ifelse(result$mean < 0, "NEGATIVE", "ZERO"))
  if (n >= 2L) {
    result$sd <- sd(values)
    result$se <- result$sd / sqrt(n)
    result$t_df <- n - 1L
    if (is.finite(result$sd) && result$sd > 0) {
      result$t_statistic <- result$mean / result$se
      result$t_p <- 2 * pt(-abs(result$t_statistic), df = result$t_df)
      critical <- qt(0.975, df = result$t_df)
      result$ci_low <- result$mean - critical * result$se
      result$ci_high <- result$mean + critical * result$se
      result$cohen_dz <- result$mean / result$sd
    }
    result$wilcoxon_p <- suppressWarnings(wilcox.test(values, mu = 0, paired = FALSE, exact = FALSE)$p.value)
    set.seed(seed_from_label(paste(project, analysis, population, metric, pseudocount, sep = "::")))
    boot <- replicate(5000L, mean(sample(values, n, replace = TRUE)))
    result$bootstrap_ci_low <- unname(quantile(boot, 0.025, type = 7L))
    result$bootstrap_ci_high <- unname(quantile(boot, 0.975, type = 7L))
    if (n >= 3L && n <= 5000L) result$shapiro_p <- shapiro.test(values)$p.value
    result$status <- ifelse(is.finite(result$t_p), "ESTIMATED", "ZERO_VARIANCE_OR_NONFINITE")
  }
  result
}

natural_projects <- c("PRJNA691455", "PRJNA851469", "PRJNA516701")
primary_rows <- list()
for (project in natural_projects) {
  selected <- pairs$project == project & pairs$in_frozen_primary_pair_population
  primary_rows[[project]] <- effect_summary(
    pairs$delta_B_pc0_5[selected], project, "PRIMARY_COHORT", "FROZEN_ALL_PAIRED",
    "delta_B", 0.5, sum(selected)
  )
}
external_primary_selected <- pairs$project == "PRJNA1125274" & pairs$A_stage_eligible
primary_rows[["PRJNA1125274"]] <- effect_summary(
  pairs$delta_B_pc0_5[external_primary_selected], "PRJNA1125274", "EXTERNAL_PRIMARY",
  "PRESPECIFIED_DATE_CLEAN_24", "delta_B", 0.5, sum(external_primary_selected)
)
cohort_effects <- do.call(rbind, primary_rows)

# Sensitivity estimates.
sensitivity_rows <- list()
add_sensitivity <- function(project, analysis, population, metric_column, pc, selected, expected = sum(selected)) {
  sensitivity_rows[[length(sensitivity_rows) + 1L]] <<- effect_summary(
    pairs[[metric_column]][selected], project, analysis, population, metric_column, pc, expected
  )
}
for (project in natural_projects) {
  base <- pairs$project == project & pairs$in_frozen_primary_pair_population
  add_sensitivity(project, "PSEUDOCOUNT_SENSITIVITY", "FROZEN_ALL_PAIRED", "delta_B_pc0_1", 0.1, base)
  add_sensitivity(project, "PSEUDOCOUNT_SENSITIVITY", "FROZEN_ALL_PAIRED", "delta_B_pc1", 1, base)
  add_sensitivity(project, "NONZERO_ONLY_SENSITIVITY", "FROZEN_ALL_PAIRED_NONZERO", "delta_B_no_pseudocount", 0, base)
  anchor <- base & pairs$anchor_input_available
  add_sensitivity(project, "COMMON_ANCHOR_SENSITIVITY", "STRICT_ANCHOR_AND_ALL_PAIRED", "delta_B_pc0_5", 0.5, anchor)
}
external_all24 <- external_primary_selected
add_sensitivity("PRJNA1125274", "PSEUDOCOUNT_SENSITIVITY", "PRESPECIFIED_DATE_CLEAN_24", "delta_B_pc0_1", 0.1, external_all24)
add_sensitivity("PRJNA1125274", "PSEUDOCOUNT_SENSITIVITY", "PRESPECIFIED_DATE_CLEAN_24", "delta_B_pc1", 1, external_all24)
add_sensitivity("PRJNA1125274", "NONZERO_ONLY_SENSITIVITY", "PRESPECIFIED_DATE_CLEAN_24_NONZERO", "delta_B_no_pseudocount", 0, external_all24)
external_complete <- pairs$project == "PRJNA1125274" & pairs$in_frozen_sensitivity_pair_population & pairs$early_late_inputs_available
add_sensitivity("PRJNA1125274", "COMPLETE_CASE_SENSITIVITY", "COMPLETE_T0_T1_T2_30", "delta_B_pc0_5", 0.5, external_complete)

external_population$in_SAL8_exclusion_sensitivity <- as_flag(external_population$in_SAL8_exclusion_sensitivity)
sal8_patients <- unique(external_population$patient_id[external_population$in_SAL8_exclusion_sensitivity])
external_sal8 <- external_all24 & pairs$patient_id %in% sal8_patients
add_sensitivity("PRJNA1125274", "SAL8_EXCLUSION_SENSITIVITY", "PRESPECIFIED_PRIMARY_WITH_SAL8_RULE", "delta_B_pc0_5", 0.5, external_sal8)
for (center in sort(unique(pairs$hospital[external_all24]))) {
  center_selected <- external_all24 & pairs$hospital == center
  add_sensitivity("PRJNA1125274", "HOSPITAL_STRATIFIED_PRIMARY", center, "delta_B_pc0_5", 0.5, center_selected)
}
sensitivity_effects <- do.call(rbind, sensitivity_rows)

# REML random-effects meta-analysis with conservative Hartung-Knapp SE.
reml_meta <- function(effect_table, analysis_label, note = "") {
  usable <- is.finite(effect_table$mean) & is.finite(effect_table$sd) & effect_table$sd > 0 & effect_table$n_analyzed >= 2
  table <- effect_table[usable, , drop = FALSE]
  k <- nrow(table)
  result <- data.frame(
    analysis = analysis_label, k = k, pooled_mean = NA_real_, tau2_REML = NA_real_,
    conventional_se = NA_real_, hk_scale_q = NA_real_, hk_raw_se = NA_real_,
    used_se = NA_real_, df = ifelse(k > 0, k - 1L, NA_integer_),
    ci_low = NA_real_, ci_high = NA_real_, t_statistic = NA_real_, p_value = NA_real_,
    Q = NA_real_, Q_df = ifelse(k > 0, k - 1L, NA_integer_), Q_p = NA_real_, I2_percent = NA_real_,
    direction = "NOT_ESTIMABLE", note = note, stringsAsFactors = FALSE
  )
  weights <- data.frame()
  if (k < 2L) return(list(summary = result, weights = weights))
  yi <- table$mean
  vi <- (table$sd ^ 2) / table$n_analyzed
  if (any(!is.finite(vi) | vi <= 0)) return(list(summary = result, weights = weights))
  objective <- function(tau2) {
    w <- 1 / (vi + tau2)
    mu <- sum(w * yi) / sum(w)
    0.5 * (sum(log(vi + tau2)) + log(sum(w)) + sum(w * (yi - mu) ^ 2))
  }
  upper <- max(1, var(yi) * 100, max(vi) * 100)
  optimized <- optimize(objective, interval = c(0, upper), tol = 1e-12)
  tau2 <- if (objective(0) <= optimized$objective + 1e-10) 0 else optimized$minimum
  w <- 1 / (vi + tau2)
  mu <- sum(w * yi) / sum(w)
  se_conventional <- sqrt(1 / sum(w))
  q_hk <- sum(w * (yi - mu) ^ 2) / (k - 1L)
  se_hk <- sqrt(q_hk / sum(w))
  se_used <- max(se_conventional, se_hk)
  critical <- qt(0.975, df = k - 1L)
  t_value <- mu / se_used
  p_value <- 2 * pt(-abs(t_value), df = k - 1L)
  w_fixed <- 1 / vi
  mu_fixed <- sum(w_fixed * yi) / sum(w_fixed)
  Q <- sum(w_fixed * (yi - mu_fixed) ^ 2)
  I2 <- ifelse(Q > 0, max(0, (Q - (k - 1L)) / Q) * 100, 0)
  result$pooled_mean <- mu
  result$tau2_REML <- tau2
  result$conventional_se <- se_conventional
  result$hk_scale_q <- q_hk
  result$hk_raw_se <- se_hk
  result$used_se <- se_used
  result$ci_low <- mu - critical * se_used
  result$ci_high <- mu + critical * se_used
  result$t_statistic <- t_value
  result$p_value <- p_value
  result$Q <- Q
  result$Q_p <- pchisq(Q, df = k - 1L, lower.tail = FALSE)
  result$I2_percent <- I2
  result$direction <- ifelse(mu > 0, "POSITIVE", ifelse(mu < 0, "NEGATIVE", "ZERO"))
  weights <- data.frame(
    analysis = analysis_label,
    project = table$project,
    yi = yi,
    vi = vi,
    random_weight = w,
    random_weight_percent = 100 * w / sum(w),
    stringsAsFactors = FALSE
  )
  list(summary = result, weights = weights)
}

meta_results <- list()
meta_weights <- list()
primary_natural <- cohort_effects[cohort_effects$project %in% natural_projects, , drop = FALSE]
meta_results[["PRIMARY_PC0_5"]] <- reml_meta(
  primary_natural, "PRIMARY_PC0_5",
  "Three natural-history cohorts; frozen ALL_PAIRED; double-zero pairs are missing"
)

meta_from_sensitivity <- function(metric_column, analysis_label, filter_analysis = NULL) {
  rows <- list()
  for (project in natural_projects) {
    base <- pairs$project == project & pairs$in_frozen_primary_pair_population
    if (identical(filter_analysis, "COMMON_ANCHOR")) base <- base & pairs$anchor_input_available
    rows[[project]] <- effect_summary(
      pairs[[metric_column]][base], project, analysis_label, "NATURAL_HISTORY", metric_column,
      ifelse(metric_column == "delta_B_pc0_1", 0.1, ifelse(metric_column == "delta_B_pc1", 1, ifelse(metric_column == "delta_B_no_pseudocount", 0, 0.5))),
      sum(base)
    )
  }
  reml_meta(do.call(rbind, rows), analysis_label, "Prespecified meta-analysis sensitivity")
}
meta_results[["SENSITIVITY_PC0_1"]] <- meta_from_sensitivity("delta_B_pc0_1", "SENSITIVITY_PC0_1")
meta_results[["SENSITIVITY_PC1"]] <- meta_from_sensitivity("delta_B_pc1", "SENSITIVITY_PC1")
meta_results[["SENSITIVITY_NONZERO_ONLY"]] <- meta_from_sensitivity("delta_B_no_pseudocount", "SENSITIVITY_NONZERO_ONLY")
meta_results[["SENSITIVITY_COMMON_ANCHOR"]] <- meta_from_sensitivity("delta_B_pc0_5", "SENSITIVITY_COMMON_ANCHOR", "COMMON_ANCHOR")
for (excluded in natural_projects) {
  table <- primary_natural[primary_natural$project != excluded, , drop = FALSE]
  label <- paste0("LEAVE_ONE_OUT_EXCLUDE_", excluded)
  meta_results[[label]] <- reml_meta(table, label, "k=2 instability diagnostic; not a primary significance result")
}
for (name in names(meta_results)) {
  meta_weights[[name]] <- meta_results[[name]]$weights
}
meta_summary <- do.call(rbind, lapply(meta_results, `[[`, "summary"))
meta_weight_table <- do.call(rbind, meta_weights)

# Independent REML implementation cross-check. Prefer metafor when installed;
# otherwise use a second base-R optimizer with the same restricted-likelihood
# objective. This checks the optimizer boundary/solution without adding a
# runtime package dependency.
metafor_check <- data.frame(
  package_available = requireNamespace("metafor", quietly = TRUE),
  comparison = "PRIMARY_PC0_5_REML_TAU2",
  internal_tau2 = meta_summary$tau2_REML[meta_summary$analysis == "PRIMARY_PC0_5"],
  crosscheck_method = NA_character_, crosscheck_tau2 = NA_real_,
  absolute_difference = NA_real_, status = "NOT_AVAILABLE",
  stringsAsFactors = FALSE
)
yi <- primary_natural$mean
vi <- (primary_natural$sd ^ 2) / primary_natural$n_analyzed
if (metafor_check$package_available) {
  fit <- metafor::rma.uni(yi = yi, vi = vi, method = "REML")
  metafor_check$crosscheck_method <- "metafor_rma_uni_REML"
  metafor_check$crosscheck_tau2 <- as.numeric(fit$tau2)
} else {
  objective_crosscheck <- function(tau2) {
    w <- 1 / (vi + tau2)
    mu <- sum(w * yi) / sum(w)
    0.5 * (sum(log(vi + tau2)) + log(sum(w)) + sum(w * (yi - mu) ^ 2))
  }
  cross_upper <- max(1, var(yi) * 100, max(vi) * 100)
  cross_fit <- nlminb(
    start = max(1e-8, metafor_check$internal_tau2),
    objective = objective_crosscheck,
    lower = 0,
    upper = cross_upper,
    control = list(eval.max = 1000L, iter.max = 1000L, rel.tol = 1e-12)
  )
  boundary_value <- objective_crosscheck(0)
  optimized_value <- objective_crosscheck(cross_fit$par)
  metafor_check$crosscheck_method <- "base_R_nlminb_REML_objective"
  metafor_check$crosscheck_tau2 <- ifelse(boundary_value <= optimized_value + 1e-10, 0, cross_fit$par)
}
metafor_check$absolute_difference <- abs(metafor_check$internal_tau2 - metafor_check$crosscheck_tau2)
metafor_check$status <- ifelse(metafor_check$absolute_difference <= 1e-7, "PASS", "REVIEW")

# Primary Holm family.
primary_meta <- meta_summary[meta_summary$analysis == "PRIMARY_PC0_5", , drop = FALSE]
external_effect <- cohort_effects[cohort_effects$project == "PRJNA1125274", , drop = FALSE]
primary_family <- data.frame(
  test_id = c("NATURAL_HISTORY_REML_HK", "PRJNA1125274_EXTERNAL_24"),
  estimate = c(primary_meta$pooled_mean, external_effect$mean),
  ci_low = c(primary_meta$ci_low, external_effect$ci_low),
  ci_high = c(primary_meta$ci_high, external_effect$ci_high),
  raw_p = c(primary_meta$p_value, external_effect$t_p),
  stringsAsFactors = FALSE
)
primary_family$holm_p <- p.adjust(primary_family$raw_p, method = "holm")
primary_family$pointwise_ci_not_multiplicity_adjusted <- TRUE
both_positive <- all(primary_family$estimate > 0)
both_holm <- all(primary_family$holm_p < 0.05)
primary_conclusion <- if (both_positive && both_holm) {
  "STRONG_COMMON_DIRECTION_SUPPORT"
} else if (both_positive) {
  "DIRECTIONALLY_CONCORDANT_BUT_NOT_BOTH_HOLM_SIGNIFICANT"
} else if (all(primary_family$estimate < 0)) {
  "CONCORDANT_OPPOSITE_TO_HYPOTHESIS"
} else {
  "DIRECTIONALLY_HETEROGENEOUS"
}
primary_family$conclusion_rule_result <- primary_conclusion

# Healthy-reference paired result with two-level bootstrap.
h_selected <- pairs$project == "PRJNA851469" & pairs$in_frozen_primary_pair_population & is.finite(pairs$delta_H)
healthy_effect <- effect_summary(
  pairs$delta_H[h_selected], "PRJNA851469", "HEALTHY_REFERENCE_SECONDARY",
  "FROZEN_ALL_PAIRED", "delta_H", NA_real_, sum(pairs$project == "PRJNA851469" & pairs$in_frozen_primary_pair_population)
)
h_pair_rows <- pairs[h_selected, , drop = FALSE]
early_matrix <- bray_to_health[h_pair_rows$early_run_id, , drop = FALSE]
late_matrix <- bray_to_health[h_pair_rows$late_run_id, , drop = FALSE]
set.seed(seed_from_label("PRJNA851469_TWO_LEVEL_HEALTHY_BOOTSTRAP"))
two_level_boot <- replicate(5000L, {
  h_idx <- sample(seq_len(ncol(early_matrix)), ncol(early_matrix), replace = TRUE)
  p_idx <- sample(seq_len(nrow(early_matrix)), nrow(early_matrix), replace = TRUE)
  patient_delta <- rowMeans(late_matrix[, h_idx, drop = FALSE]) - rowMeans(early_matrix[, h_idx, drop = FALSE])
  mean(patient_delta[p_idx])
})
healthy_effect$two_level_bootstrap_ci_low <- unname(quantile(two_level_boot, 0.025, type = 7L))
healthy_effect$two_level_bootstrap_ci_high <- unname(quantile(two_level_boot, 0.975, type = 7L))
healthy_effect$n_healthy_references <- nrow(healthy_rel)
healthy_effect$independence_status <- "SAME_COHORT_SECONDARY_REFERENCE"

# Same-source association summaries.
spearman_summary <- function(df, x, y, project, analysis) {
  keep <- is.finite(df[[x]]) & is.finite(df[[y]])
  n <- sum(keep)
  rho <- p <- NA_real_
  if (n >= 3L && length(unique(df[[x]][keep])) > 1L && length(unique(df[[y]][keep])) > 1L) {
    test <- suppressWarnings(cor.test(df[[x]][keep], df[[y]][keep], method = "spearman", exact = FALSE))
    rho <- unname(test$estimate)
    p <- test$p.value
  }
  data.frame(
    project = project, analysis = analysis, x = x, y = y, n = n,
    spearman_rho = rho, p_value = p,
    multiplicity_status = "EXPLORATORY_NOT_IN_PRIMARY_HOLM_FAMILY",
    independence_status = "SAME_16S_SOURCE_CONSISTENCY_ANALYSIS",
    stringsAsFactors = FALSE
  )
}
association_rows <- list()
for (project in c(natural_projects, "PRJNA1125274")) {
  selected <- pairs$project == project & if (project == "PRJNA1125274") pairs$A_stage_eligible else pairs$in_frozen_primary_pair_population
  association_rows[[length(association_rows) + 1L]] <- spearman_summary(
    pairs[selected, , drop = FALSE], "delta_D", "delta_B_pc0_5", project, "DELTA_D_VS_DELTA_B"
  )
}
association_rows[[length(association_rows) + 1L]] <- spearman_summary(
  pairs[h_selected, , drop = FALSE], "delta_D", "delta_H", "PRJNA851469", "DELTA_D_VS_DELTA_H"
)
association_rows[[length(association_rows) + 1L]] <- spearman_summary(
  pairs[h_selected, , drop = FALSE], "delta_B_pc0_5", "delta_H", "PRJNA851469", "DELTA_B_VS_DELTA_H"
)
associations <- do.call(rbind, association_rows)
associations$fdr_within_exploratory_family <- p.adjust(associations$p_value, method = "BH")

# Descriptive transition quadrants.
main_pair <- (pairs$project %in% natural_projects & pairs$in_frozen_primary_pair_population) |
  (pairs$project == "PRJNA1125274" & pairs$A_stage_eligible)
transition_patient <- pairs[main_pair & is.finite(pairs$delta_D) & is.finite(pairs$delta_B_pc0_5),
                            c("project", "patient_id", "hospital", "delta_D", "delta_B_pc0_5"), drop = FALSE]
transition_patient$transition_quadrant <- ifelse(
  transition_patient$delta_D > 0 & transition_patient$delta_B_pc0_5 > 0,
  "GREATER_DISPLACEMENT__P_RELATIVE_INCREASE",
  ifelse(
    transition_patient$delta_D > 0 & transition_patient$delta_B_pc0_5 < 0,
    "GREATER_DISPLACEMENT__C_RELATIVE_INCREASE",
    ifelse(
      transition_patient$delta_D < 0 & transition_patient$delta_B_pc0_5 > 0,
      "LESS_DISPLACEMENT__P_RELATIVE_INCREASE",
      ifelse(
        transition_patient$delta_D < 0 & transition_patient$delta_B_pc0_5 < 0,
        "LESS_DISPLACEMENT__C_RELATIVE_INCREASE", "AXIS_TIE"
      )
    )
  )
)
transition_summary <- aggregate(
  list(n_patients = transition_patient$patient_id),
  by = list(project = transition_patient$project, transition_quadrant = transition_patient$transition_quadrant),
  FUN = length
)
project_totals <- aggregate(n_patients ~ project, data = transition_summary, FUN = sum)
transition_summary$fraction <- transition_summary$n_patients / project_totals$n_patients[match(transition_summary$project, project_totals$project)]

# Exclusion ledger for all frozen relevant populations.
relevant_frozen <- (pairs$project %in% natural_projects & pairs$in_frozen_primary_pair_population) |
  (pairs$project == "PRJNA1125274" & pairs$in_frozen_sensitivity_pair_population)
exclusion_ledger <- pairs[relevant_frozen, c(
  "project", "patient_id", "hospital", "in_frozen_primary_pair_population",
  "in_frozen_sensitivity_pair_population", "A_stage_eligible", "anchor_input_available",
  "early_run_id", "late_run_id", "early_P_reads", "early_C_reads", "late_P_reads", "late_C_reads",
  "any_P_or_C_zero_at_early_or_late", "any_P_plus_C_double_zero_at_early_or_late",
  "main_balance_exclusion_reason"
), drop = FALSE]
exclusion_ledger$main_delta_B_calculable <- exclusion_ledger$main_balance_exclusion_reason == ""

# Two independent hand calculations from raw counts for sign/formula verification.
handcheck_candidates <- pairs[
  ((pairs$project %in% natural_projects & pairs$in_frozen_primary_pair_population) |
     (pairs$project == "PRJNA1125274" & pairs$A_stage_eligible)) &
    is.finite(pairs$delta_B_pc0_5),
  , drop = FALSE
]
handcheck_candidates <- handcheck_candidates[order(handcheck_candidates$project, handcheck_candidates$patient_id), , drop = FALSE]
handcheck <- handcheck_candidates[c(1L, nrow(handcheck_candidates)), c(
  "project", "patient_id", "early_P_reads", "early_C_reads", "late_P_reads", "late_C_reads",
  "early_B_pc0_5", "late_B_pc0_5", "delta_B_pc0_5"
), drop = FALSE]
handcheck$manual_early_B <- log((handcheck$early_P_reads + 0.5) / (handcheck$early_C_reads + 0.5))
handcheck$manual_late_B <- log((handcheck$late_P_reads + 0.5) / (handcheck$late_C_reads + 0.5))
handcheck$manual_delta_B <- handcheck$manual_late_B - handcheck$manual_early_B
handcheck$max_absolute_error <- pmax(
  abs(handcheck$manual_early_B - handcheck$early_B_pc0_5),
  abs(handcheck$manual_late_B - handcheck$late_B_pc0_5),
  abs(handcheck$manual_delta_B - handcheck$delta_B_pc0_5)
)
handcheck$status <- ifelse(handcheck$max_absolute_error <= 1e-12, "PASS", "FAIL")

# Output tables before figures.
write_csv_once(samples, "B01_sample_direction_metrics.csv")
write_csv_once(pairs, "B02_patient_paired_direction_metrics.csv")
write_csv_once(exclusion_ledger, "B03_pair_exclusion_ledger.csv")
write_csv_once(cohort_effects, "B04_primary_cohort_effects.csv")
write_csv_once(meta_summary, "B05_natural_history_meta_analysis.csv")
write_csv_once(meta_weight_table, "B05B_meta_analysis_weights.csv")
write_csv_once(primary_family, "B06_primary_holm_family.csv")
write_csv_once(sensitivity_effects, "B07_cohort_and_external_sensitivity_results.csv")
write_csv_once(healthy_sample_metrics, "B08_PRJNA851469_healthy_reference_sample_metrics.csv")
write_csv_once(healthy_effect, "B09_PRJNA851469_healthy_reference_paired_result.csv")
write_csv_once(associations, "B10_D_B_H_associations.csv")
write_csv_once(transition_patient, "B11_transition_quadrant_patient_table.csv")
write_csv_once(transition_summary, "B11B_transition_quadrant_summary.csv")
write_csv_once(handcheck, "B12_formula_handcheck.csv")
write_csv_once(metafor_check, "B13_meta_implementation_crosscheck.csv")

# Figures use only displayed primary/secondary results; no threshold optimization.
figure_rows <- rbind(
  data.frame(
    label = cohort_effects$project,
    estimate = cohort_effects$mean,
    low = cohort_effects$ci_low,
    high = cohort_effects$ci_high,
    type = ifelse(cohort_effects$project == "PRJNA1125274", "External primary", "Natural cohort"),
    stringsAsFactors = FALSE
  ),
  data.frame(
    label = "Natural-history REML-HK",
    estimate = primary_meta$pooled_mean,
    low = primary_meta$ci_low,
    high = primary_meta$ci_high,
    type = "Pooled",
    stringsAsFactors = FALSE
  )
)
write_csv_once(figure_rows, "B14_forest_plot_source.csv")

draw_forest <- function() {
  table <- figure_rows[rev(seq_len(nrow(figure_rows))), , drop = FALSE]
  finite <- is.finite(table$estimate) & is.finite(table$low) & is.finite(table$high)
  xlim <- range(c(table$low[finite], table$high[finite], 0), finite = TRUE)
  pad <- diff(xlim) * 0.12
  plot(table$estimate, seq_len(nrow(table)), xlim = xlim + c(-pad, pad),
       ylim = c(0.5, nrow(table) + 0.5), yaxt = "n", ylab = "", xlab = "Mean early-to-late delta_B",
       pch = ifelse(table$type == "Pooled", 18, 19), cex = ifelse(table$type == "Pooled", 1.4, 1.0),
       main = "Direction balance change across cohorts")
  abline(v = 0, lty = 2, col = "grey50")
  segments(table$low, seq_len(nrow(table)), table$high, seq_len(nrow(table)), lwd = 2)
  axis(2, at = seq_len(nrow(table)), labels = table$label, las = 1, cex.axis = 0.85)
  legend("bottomright", legend = c("Natural cohort", "External primary", "Pooled"),
         pch = c(19, 19, 18), col = c("black", "black", "black"), bty = "n")
}
png(file.path(out_dir, "figures", "Figure_B1_direction_balance_forest.png"), width = 2400, height = 1600, res = 300)
par(mar = c(5, 12, 3, 2)); draw_forest(); dev.off()
pdf(file.path(out_dir, "figures", "Figure_B1_direction_balance_forest.pdf"), width = 10, height = 5.5)
par(mar = c(5, 12, 3, 2)); draw_forest(); dev.off()

trajectory_data <- pairs[main_pair & is.finite(pairs$delta_B_pc0_5),
                         c("project", "patient_id", "early_B_pc0_5", "late_B_pc0_5"), drop = FALSE]
draw_trajectories <- function() {
  par(mfrow = c(2, 2), mar = c(4, 4, 3, 1))
  y_range <- range(c(trajectory_data$early_B_pc0_5, trajectory_data$late_B_pc0_5), finite = TRUE)
  for (project in c(natural_projects, "PRJNA1125274")) {
    d <- trajectory_data[trajectory_data$project == project, , drop = FALSE]
    plot(c(1, 2), y_range, type = "n", xaxt = "n", xlab = "", ylab = "B = log[(P+0.5)/(C+0.5)]", main = project)
    axis(1, at = c(1, 2), labels = c("Early", "Late"))
    abline(h = 0, lty = 3, col = "grey70")
    for (i in seq_len(nrow(d))) lines(c(1, 2), c(d$early_B_pc0_5[i], d$late_B_pc0_5[i]), col = rgb(0, 0, 0, 0.35))
    points(rep(1, nrow(d)), d$early_B_pc0_5, pch = 16, cex = 0.55)
    points(rep(2, nrow(d)), d$late_B_pc0_5, pch = 16, cex = 0.55)
  }
}
png(file.path(out_dir, "figures", "Figure_B2_patient_balance_trajectories.png"), width = 2400, height = 1800, res = 300)
draw_trajectories(); dev.off()
pdf(file.path(out_dir, "figures", "Figure_B2_patient_balance_trajectories.pdf"), width = 9, height = 7)
draw_trajectories(); dev.off()

draw_db <- function() {
  par(mfrow = c(2, 2), mar = c(4, 4, 3, 1))
  for (project in c(natural_projects, "PRJNA1125274")) {
    d <- transition_patient[transition_patient$project == project, , drop = FALSE]
    plot(d$delta_D, d$delta_B_pc0_5, pch = 19, xlab = "delta_D (late - early Bray displacement)",
         ylab = "delta_B (late - early)", main = project)
    abline(h = 0, v = 0, lty = 3, col = "grey60")
  }
}
png(file.path(out_dir, "figures", "Figure_B3_displacement_vs_direction.png"), width = 2400, height = 1800, res = 300)
draw_db(); dev.off()
pdf(file.path(out_dir, "figures", "Figure_B3_displacement_vs_direction.pdf"), width = 9, height = 7)
draw_db(); dev.off()

draw_h <- function() {
  d <- pairs[h_selected, , drop = FALSE]
  y_range <- range(c(d$early_H, d$late_H), finite = TRUE)
  plot(c(1, 2), y_range, type = "n", xaxt = "n", xlab = "", ylab = "Mean Bray-Curtis distance to 13 healthy references",
       main = "PRJNA851469 distance to healthy reference")
  axis(1, at = c(1, 2), labels = c("Day-3", "Day-7"))
  for (i in seq_len(nrow(d))) lines(c(1, 2), c(d$early_H[i], d$late_H[i]), col = rgb(0, 0, 0, 0.4))
  points(rep(1, nrow(d)), d$early_H, pch = 16); points(rep(2, nrow(d)), d$late_H, pch = 16)
}
png(file.path(out_dir, "figures", "Figure_B4_healthy_reference_trajectory.png"), width = 1800, height = 1600, res = 300)
par(mar = c(5, 5, 3, 1)); draw_h(); dev.off()
pdf(file.path(out_dir, "figures", "Figure_B4_healthy_reference_trajectory.pdf"), width = 6.5, height = 5.5)
par(mar = c(5, 5, 3, 1)); draw_h(); dev.off()

# Machine-readable run decision and manuscript-safe candidate wording.
decision_table <- data.frame(
  item = c(
    "primary_conclusion_rule", "natural_meta_direction", "external_direction",
    "both_holm_p_below_0_05", "healthy_reference_available", "patient_level_host_model_available"
  ),
  value = c(
    primary_conclusion,
    primary_meta$direction,
    external_effect$direction,
    as.character(both_holm),
    "TRUE_13_FROZEN_HEALTHY_SAMPLES",
    "FALSE_PUBLIC_JOINABLE_HOST_VALUES_ABSENT"
  ),
  stringsAsFactors = FALSE
)
write_csv_once(decision_table, "B15_result_decision_table.csv")

manuscript_lines <- c(
  "# Stage-B candidate Results text",
  "",
  "## Material Passport",
  "",
  "- Origin Skill: academic-research-suite / experiment-agent",
  "- Origin Mode: run",
  "- Verification Status: ANALYSIS_COMPLETE_PENDING_VALIDATION_AND_HUMAN_EDITORIAL_REVIEW",
  "- Version: direction_stage_b_results_candidate_v1",
  "",
  "## Primary direction result",
  "",
  sprintf(
    "Across the three prespecified natural-history cohorts, the random-effects REML estimate for early-to-late change in the fixed P/C log-ratio was %s (95%% CI %s to %s; Hartung-Knapp p=%s; tau-squared=%s; I-squared=%s%%).",
    safe_value(primary_meta$pooled_mean, 3), safe_value(primary_meta$ci_low, 3), safe_value(primary_meta$ci_high, 3),
    safe_value(primary_meta$p_value, 4), safe_value(primary_meta$tau2_REML, 3), safe_value(primary_meta$I2_percent, 1)
  ),
  sprintf(
    "In the frozen 24-patient PRJNA1125274 primary population, the mean change was %s (95%% CI %s to %s; two-sided t p=%s).",
    safe_value(external_effect$mean, 3), safe_value(external_effect$ci_low, 3), safe_value(external_effect$ci_high, 3), safe_value(external_effect$t_p, 4)
  ),
  sprintf(
    "The prespecified two-test Holm-adjusted p values were %s and %s, yielding the rule-based conclusion %s.",
    safe_value(primary_family$holm_p[1], 4), safe_value(primary_family$holm_p[2], 4), primary_conclusion
  ),
  "Positive delta_B denotes an increase in P relative to C and does not establish absolute expansion or host causation.",
  "",
  "## Healthy-reference secondary result",
  "",
  sprintf(
    "Using 13 co-processed healthy PRJNA851469 samples as a same-cohort reference, the mean Day-7 minus Day-3 change in distance to the healthy reference was %s (ordinary paired 95%% CI %s to %s; two-level bootstrap 95%% CI %s to %s; p=%s; n=%d).",
    safe_value(healthy_effect$mean, 3), safe_value(healthy_effect$ci_low, 3), safe_value(healthy_effect$ci_high, 3),
    safe_value(healthy_effect$two_level_bootstrap_ci_low, 3), safe_value(healthy_effect$two_level_bootstrap_ci_high, 3),
    safe_value(healthy_effect$t_p, 4), healthy_effect$n_analyzed
  ),
  "This healthy reference is not an independent validation cohort.",
  "",
  "## Required limitations",
  "",
  "- PRJNA516701 loses three frozen pairs because P+C was zero at an early or late sample; these values were not converted to a neutral balance.",
  "- PRJNA1125274 has substantially lower Family-level resolution than the natural-history objects; unresolved reads remain in OTHER and in total-read accounting.",
  "- All D-B associations use related summaries derived from the same 16S data and are exploratory consistency analyses.",
  "- No joinable patient-level host measurement was available, so no host model was fitted.",
  "- The previously negative PRJNA851469 clinical outcome analysis remains unchanged and is not replaced by this direction analysis."
)
write_lines_once(manuscript_lines, "B16_candidate_results_text.md")

session_path <- file.path(out_dir, "session_info.txt")
if (file.exists(session_path)) stop("Refusing to overwrite: ", session_path)
writeLines(c(
  paste0("timestamp=", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
  "script=E:/sepsis_project/code/03_data_processing/99_DIRECTION_UPGRADE/99B1_direction_analysis.R",
  "scope=stage B fixed secondary analysis; no sequence preprocessing",
  capture.output(sessionInfo())
), session_path, useBytes = TRUE)

message("Stage-B analysis completed. Rule-based conclusion: ", primary_conclusion)
