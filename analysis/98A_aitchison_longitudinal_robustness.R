# ============================================================
# Sepsis V2 Upgrade - Step 98A
# Aitchison (CLR / compositional) longitudinal robustness
# ------------------------------------------------------------
# Purpose:
#   Re-run the frozen Step88A2 within-patient longitudinal
#   displacement framework under the Aitchison compositional
#   geometry, per cohort, WITHOUT altering any frozen input or
#   output.
#
# Frozen inputs used (read-only):
#   - Step87B cohort-specific RDS analysis objects
#     (obj$counts, obj$metadata, obj$provenance)
#   - Step88A2 sample-level Bray-Curtis displacement CSV
#     (framework of samples / patients / baseline / time labels)
#   - Step88A2 model global time tests (for QC mirroring)
#   - Step88B paired early-vs-late contrast registry (for the
#     pre-specified contrast structure; BC values recomputed here
#     and cross-checked against frozen Step88B numbers)
#
# Rules (mirroring Step88A2/88B):
#   - cohort-specific analysis only; no cross-cohort ASV merge
#   - patients with >=2 timepoints (653-sample GE2 framework)
#   - each patient's earliest longitudinal sample = personal
#     baseline; baseline distance = 0 (QC hard guard)
#   - Bray-Curtis on relative abundance (recomputed -> QC)
#   - Aitchison distance = Euclidean distance between CLR vectors
#   - zero handling fixed BEFORE results: primary = CZM
#     (zCompositions::cmultRepl), sensitivity 1 = pseudocount 0.5,
#     sensitivity 2 = pseudocount 1. No result-driven selection.
#   - follow-up-only mixed models (baseline rows excluded from
#     inferential models)
#   - early-vs-late paired contrasts identical to Step88B
#     (ALL_PAIRED + COMMON_ANCHOR_SENSITIVITY)
#
# All outputs -> E:/sepsis_project/results/V2_UPGRADE_20260907/
#                    98A_AITCHISON_LONGITUDINAL_ROBUSTNESS/
# ============================================================

options(stringsAsFactors = FALSE)

required_pkgs <- c(
  "readr", "dplyr", "tidyr", "purrr", "tibble", "stringr",
  "vegan", "lme4", "lmerTest", "ggplot2", "zCompositions"
)

missing_pkgs <- required_pkgs[
  !vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_pkgs)) {
  install.packages(missing_pkgs, repos = "https://cloud.r-project.org")
}

suppressPackageStartupMessages({
  library(readr)
  library(tibble)
  library(purrr)
  library(stringr)
  library(vegan)
  library(lme4)
  library(lmerTest)
  library(zCompositions)
  # attach dplyr/tidyr last so dplyr::select masks MASS::select (MASS is
  # attached as a dependency of vegan)
  library(dplyr)
  library(tidyr)
})

set.seed(20260907)

# ------------------------------------------------------------
# Paths
# ------------------------------------------------------------

ROOT <- "E:/sepsis_project"

STEP87B_FREEZE <- file.path(
  ROOT, "results",
  "V2_27B_ANALYSIS_OBJECTS_AND_REPLICATE_FREEZE"
)

STEP88A2_RES <- file.path(
  ROOT, "results",
  "V2_28A2_STEP88A_LONGITUDINAL_DIVERSITY_AND_DISPLACEMENT_FIXED"
)

STEP88B_RES <- file.path(
  ROOT, "results",
  "V2_28B_STEP88B_ANCHOR_ROBUSTNESS_AND_PAIRED_CONTRASTS"
)

OUT <- file.path(
  ROOT, "results", "V2_UPGRADE_20260907",
  "98A_AITCHISON_LONGITUDINAL_ROBUSTNESS"
)

dirs <- c(
  OUT,
  file.path(OUT, "01_SAMPLE_DISPLACEMENT"),
  file.path(OUT, "02_TIMEPOINT_SUMMARIES"),
  file.path(OUT, "03_PAIRED_CONTRASTS"),
  file.path(OUT, "03_PAIRED_CONTRASTS", "PATIENT_LEVEL_DIFFERENCES"),
  file.path(OUT, "04_MIXED_MODELS"),
  file.path(OUT, "05_QC_AND_CORRELATION"),
  file.path(OUT, "06_ZERO_HANDLING_SENSITIVITY"),
  file.path(OUT, "07_EXPLORATORY_PLOTS")
)

for (d in dirs) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

LOG <- file.path(OUT, "_STEP98A_runtime_checkpoints.txt")
ERR <- file.path(OUT, "_STEP98A_FATAL_ERROR.txt")

if (file.exists(LOG)) unlink(LOG)
if (file.exists(ERR)) unlink(ERR)

ck <- function(x) {
  cat(paste0(x, ": ", Sys.time(), "\n"), file = LOG, append = TRUE)
}

write_csv_safe <- function(x, p) {
  write_excel_csv(x, p, na = "")
}

safe_csv <- function(p) {
  suppressMessages(
    read_csv(p, show_col_types = FALSE, progress = FALSE, name_repair = "unique")
  )
}

bh_adjust <- function(p) {
  out <- rep(NA_real_, length(p))
  ok <- !is.na(p)
  if (any(ok)) out[ok] <- p.adjust(p[ok], method = "BH")
  out
}

# ------------------------------------------------------------
# Copied verbatim from Step88A2/88B so the mirroring is exact
# ------------------------------------------------------------

paired_boot_ci <- function(diff, stat = c("mean", "median"),
                           B = 2000, seed = 20260907) {
  stat <- match.arg(stat)
  x <- as.numeric(diff[!is.na(diff)])
  n <- length(x)
  if (n < 3) return(c(low = NA_real_, high = NA_real_))
  set.seed(seed)
  vals <- replicate(B, {
    z <- sample(x, size = n, replace = TRUE)
    if (stat == "mean") mean(z) else median(z)
  })
  as.numeric(quantile(vals, probs = c(0.025, 0.975), na.rm = TRUE,
                      names = FALSE)) |>
    setNames(c("low", "high"))
}

safe_wilcox <- function(x) {
  x <- as.numeric(x[!is.na(x)])
  if (!length(x)) return(NA_real_)
  if (all(abs(x) < 1e-15)) return(1)
  tryCatch(
    suppressWarnings(
      wilcox.test(x, mu = 0, paired = FALSE, exact = FALSE,
                  correct = FALSE)$p.value
    ),
    error = function(e) NA_real_
  )
}

paired_effect <- function(d, project, early_label, late_label, metric,
                          subset_type = "ALL_PAIRED",
                          required_anchor = NA_character_) {
  dd <- d |> filter(.data$project == .env$project)

  if (subset_type == "COMMON_ANCHOR_SENSITIVITY") {
    if (is.na(required_anchor) || required_anchor == "") {
      stop(project, ": common-anchor sensitivity requested without an anchor.")
    }
    anchored_ids <- dd |>
      filter(is_patient_reference, time_factor == required_anchor) |>
      distinct(patient_id) |>
      pull(patient_id)
    dd <- dd |> filter(patient_id %in% anchored_ids)
  }

  pair <- dd |>
    filter(time_factor %in% c(early_label, late_label)) |>
    select(patient_id, time_factor, all_of(metric)) |>
    distinct(patient_id, time_factor, .keep_all = TRUE) |>
    pivot_wider(names_from = time_factor, values_from = all_of(metric))

  if (!(early_label %in% names(pair)) || !(late_label %in% names(pair))) {
    return(tibble(
      project = project, metric = metric, subset_type = subset_type,
      required_anchor = required_anchor, early_label = early_label,
      late_label = late_label, n_pairs = 0L,
      early_median = NA_real_, late_median = NA_real_,
      median_difference = NA_real_, mean_difference = NA_real_,
      mean_diff_ci_low = NA_real_, mean_diff_ci_high = NA_real_,
      median_diff_ci_low = NA_real_, median_diff_ci_high = NA_real_,
      paired_sd_difference = NA_real_, paired_effect_dz = NA_real_,
      wilcoxon_p = NA_real_, direction = "NO_DATA"
    ))
  }

  pair <- pair |> filter(!is.na(.data[[early_label]]),
                         !is.na(.data[[late_label]]))
  diff <- pair[[late_label]] - pair[[early_label]]
  n <- length(diff)

  if (!n) {
    return(tibble(
      project = project, metric = metric, subset_type = subset_type,
      required_anchor = required_anchor, early_label = early_label,
      late_label = late_label, n_pairs = 0L,
      early_median = NA_real_, late_median = NA_real_,
      median_difference = NA_real_, mean_difference = NA_real_,
      mean_diff_ci_low = NA_real_, mean_diff_ci_high = NA_real_,
      median_diff_ci_low = NA_real_, median_diff_ci_high = NA_real_,
      paired_sd_difference = NA_real_, paired_effect_dz = NA_real_,
      wilcoxon_p = NA_real_, direction = "NO_DATA"
    ))
  }

  sd_diff <- if (n >= 2) sd(diff, na.rm = TRUE) else NA_real_
  dz <- if (!is.na(sd_diff) && sd_diff > 0) {
    mean(diff, na.rm = TRUE) / sd_diff
  } else NA_real_

  ci_mean <- paired_boot_ci(diff, stat = "mean")
  ci_med <- paired_boot_ci(diff, stat = "median")
  mean_diff <- mean(diff, na.rm = TRUE)

  direction <- case_when(
    mean_diff > 0 ~ "INCREASE",
    mean_diff < 0 ~ "DECREASE",
    TRUE ~ "NO_CHANGE"
  )

  tibble(
    project = project, metric = metric, subset_type = subset_type,
    required_anchor = required_anchor, early_label = early_label,
    late_label = late_label, n_pairs = n,
    early_median = median(pair[[early_label]], na.rm = TRUE),
    late_median = median(pair[[late_label]], na.rm = TRUE),
    median_difference = median(diff, na.rm = TRUE),
    mean_difference = mean_diff,
    mean_diff_ci_low = ci_mean["low"],
    mean_diff_ci_high = ci_mean["high"],
    median_diff_ci_low = ci_med["low"],
    median_diff_ci_high = ci_med["high"],
    paired_sd_difference = sd_diff,
    paired_effect_dz = dz,
    wilcoxon_p = safe_wilcox(diff),
    direction = direction
  )
}

# safe_lmer verbatim from Step88A2
safe_lmer <- function(d, outcome, time_mode = c("factor", "numeric"),
                      include_depth = TRUE, adjustment = NA_character_) {
  time_mode <- match.arg(time_mode)
  result_empty <- list(status = "NOT_RUN", global = tibble(),
                       coefficients = tibble())
  needed <- c(outcome, "patient_id")
  if (time_mode == "factor") needed <- c(needed, "time_factor")
  else needed <- c(needed, "time_day")
  if (include_depth) needed <- c(needed, "log10_library_size")
  if (!is.na(adjustment)) needed <- c(needed, adjustment)
  needed <- unique(needed)
  if (!all(needed %in% names(d))) {
    result_empty$status <- paste0(
      "MISSING_COLUMNS: ",
      paste(setdiff(needed, names(d)), collapse = ", "))
    return(result_empty)
  }
  dd <- d |> select(all_of(needed)) |> drop_na()
  if (nrow(dd) < 8 || n_distinct(dd$patient_id) < 4) {
    result_empty$status <- "TOO_FEW_COMPLETE_OBSERVATIONS"
    return(result_empty)
  }
  if (time_mode == "factor") {
    dd$time_factor <- droplevels(factor(dd$time_factor))
    if (nlevels(dd$time_factor) < 2) {
      result_empty$status <- "TIME_FACTOR_HAS_LT2_LEVELS"
      return(result_empty)
    }
    time_term <- "time_factor"
  } else {
    if (n_distinct(dd$time_day) < 2) {
      result_empty$status <- "TIME_NUMERIC_HAS_LT2_VALUES"
      return(result_empty)
    }
    time_term <- "time_day"
  }
  fixed_base <- character(0)
  if (include_depth) fixed_base <- c(fixed_base, "log10_library_size")
  if (!is.na(adjustment)) {
    dd[[adjustment]] <- factor(dd[[adjustment]])
    if (nlevels(dd[[adjustment]]) >= 2) fixed_base <- c(fixed_base, adjustment)
  }
  full_terms <- c(time_term, fixed_base)
  reduced_terms <- fixed_base
  rhs_full <- paste(c(full_terms, "(1 | patient_id)"), collapse = " + ")
  rhs_reduced <- paste(c(
    if (length(reduced_terms)) reduced_terms else "1",
    "(1 | patient_id)"), collapse = " + ")
  f_full <- as.formula(paste(outcome, "~", rhs_full))
  f_reduced <- as.formula(paste(outcome, "~", rhs_reduced))
  fit_full <- tryCatch(
    lmerTest::lmer(f_full, data = dd, REML = FALSE,
                   control = lme4::lmerControl(optimizer = "bobyqa",
                                               optCtrl = list(maxfun = 2e5))),
    error = function(e) e)
  if (inherits(fit_full, "error")) {
    result_empty$status <- paste0("FULL_MODEL_ERROR: ", conditionMessage(fit_full))
    return(result_empty)
  }
  fit_reduced <- tryCatch(
    lmerTest::lmer(f_reduced, data = dd, REML = FALSE,
                   control = lme4::lmerControl(optimizer = "bobyqa",
                                               optCtrl = list(maxfun = 2e5))),
    error = function(e) e)
  if (inherits(fit_reduced, "error")) {
    result_empty$status <- paste0("REDUCED_MODEL_ERROR: ", conditionMessage(fit_reduced))
    return(result_empty)
  }
  lrt <- tryCatch(anova(fit_reduced, fit_full), error = function(e) e)
  global_p <- NA_real_; chisq <- NA_real_; df_diff <- NA_real_
  if (!inherits(lrt, "error")) {
    lrt_df <- as.data.frame(lrt)
    if (nrow(lrt_df) >= 2) {
      if ("Pr(>Chisq)" %in% names(lrt_df)) global_p <- lrt_df[2, "Pr(>Chisq)"]
      if ("Chisq" %in% names(lrt_df)) chisq <- lrt_df[2, "Chisq"]
      if ("Chi Df" %in% names(lrt_df)) df_diff <- lrt_df[2, "Chi Df"]
    }
  }
  cc <- as.data.frame(coef(summary(fit_full)))
  cc <- rownames_to_column(cc, "term")
  names(cc) <- make.names(names(cc))
  p_col <- grep("^Pr", names(cc), value = TRUE)
  coef_out <- tibble(
    term = cc$term,
    estimate = cc$Estimate,
    std_error = cc$Std..Error,
    df = if ("df" %in% names(cc)) cc$df else NA_real_,
    statistic = if ("t.value" %in% names(cc)) cc$t.value else NA_real_,
    p_value = if (length(p_col)) cc[[p_col[1]]] else NA_real_)
  singular <- lme4::isSingular(fit_full, tol = 1e-4)
  conv_messages <- fit_full@optinfo$conv$lme4$messages
  if (is.null(conv_messages)) conv_messages <- "" else {
    conv_messages <- paste(conv_messages, collapse = "; ")
  }
  global <- tibble(
    n_observations = nrow(dd),
    n_patients = n_distinct(dd$patient_id),
    time_mode = time_mode,
    outcome = outcome,
    adjustment = ifelse(is.na(adjustment), "", adjustment),
    include_depth = include_depth,
    likelihood_ratio_chisq = chisq,
    df_difference = df_diff,
    global_time_p = global_p,
    singular_fit = singular,
    convergence_message = conv_messages)
  list(status = "OK", global = global, coefficients = coef_out)
}

resolve_adjustment <- function(project, d) {
  candidate <- switch(project,
    "PRJEB82425" = "phenotype",
    "PRJNA1166732" = "intervention_arm",
    "PRJNA430161" = "intervention_arm",
    NA_character_)
  if (is.na(candidate)) return(NA_character_)
  if (!(candidate %in% names(d))) return(NA_character_)
  vals <- d[[candidate]]
  if (is.factor(vals)) vals <- as.character(vals)
  vals <- vals[!is.na(vals)]
  if (length(unique(trimws(vals))) < 2) return(NA_character_)
  candidate
}

# ------------------------------------------------------------
# CLR / Aitchison helper
# ------------------------------------------------------------

clr_matrix <- function(x) {
  # x: matrix with rows = samples, cols = parts, all > 0
  logx <- log(x)
  logx - rowMeans(logx)
}

zero_impute <- function(counts, method = c("CZM", "PC0.5", "PC1")) {
  method <- match.arg(method)
  if (method == "CZM") {
    out <- tryCatch(
      zCompositions::cmultRepl(
        counts, label = 0, method = "CZM",
        output = "p-counts", z.warning = 1,
        z.delete = FALSE, suppress.print = TRUE
      ),
      error = function(e) e
    )
    if (inherits(out, "error")) {
      stop("CZM failed: ", conditionMessage(out))
    }
    out
  } else if (method == "PC0.5") {
    counts + 0.5
  } else {
    counts + 1
  }
}

compute_aitchison_displacement <- function(counts, framework) {
  # counts: rows = run_id (framework rows only), cols = ASVs
  # framework: data.frame with run_id, patient_id, baseline_run_id
  # returns named vector aitchison_from_patient_baseline for each framework row
  feat_total <- colSums(counts)
  if (any(feat_total == 0)) {
    counts <- counts[, feat_total > 0, drop = FALSE]
  }
  imp <- zero_impute(counts, method = "CZM")
  clr <- clr_matrix(imp)
  dst <- as.matrix(dist(clr, method = "euclidean"))
  out <- numeric(nrow(framework))
  for (i in seq_len(nrow(framework))) {
    rid <- framework$run_id[i]
    bid <- framework$baseline_run_id[i]
    out[i] <- dst[rid, bid]
  }
  out
}

# ------------------------------------------------------------
# Load frozen framework
# ------------------------------------------------------------

ck("STEP98A STARTED")

registry <- safe_csv(file.path(
  STEP87B_FREEZE, "V2_STEP87B_analysis_object_registry.csv"))

disp_all <- safe_csv(file.path(
  STEP88A2_RES, "03_BETA_DISPLACEMENT",
  "V2_STEP88A2_ALL_within_patient_bray_displacement.csv"))

# Step88A2 frozen model global tests (Bray displacement rows) for QC
frozen_global <- safe_csv(file.path(
  STEP88A2_RES, "05_MODELS", "V2_STEP88A2_model_global_time_tests.csv"))
frozen_model_status <- safe_csv(file.path(
  STEP88A2_RES, "05_MODELS", "V2_STEP88A2_model_status.csv"))

# Step88B pre-specified contrasts (Bray), used as the structural registry
step88b_contrasts <- safe_csv(file.path(
  STEP88B_RES, "02_PAIRED_CONTRASTS",
  "V2_STEP88B_paired_early_vs_late_followup_contrasts.csv"))

contrast_registry <- step88b_contrasts |>
  filter(metric == "bray_from_patient_baseline", subset_type == "ALL_PAIRED") |>
  select(project, required_anchor, early_label, late_label)

projects <- registry$project
longitudinal_projects <- unique(disp_all$project)

cat("Framework rows:", nrow(disp_all), "\n")
cat("Cohorts in framework:", paste(longitudinal_projects, collapse = ", "), "\n")

stopifnot(all(contrast_registry$project %in% longitudinal_projects))

# ------------------------------------------------------------
# Output collectors
# ------------------------------------------------------------

sample_rows <- list()
bray_qc_rows <- list()
summary_rows <- list()
paired_rows <- list()
patient_diff_rows <- list()
model_status_rows <- list()
model_global_rows <- list()
model_coef_rows <- list()
model_qc_rows <- list()
spearman_rows <- list()

for (proj in longitudinal_projects) {

  ck(paste0("PROCESSING ", proj))

  rr <- registry |> filter(project == proj)
  obj <- readRDS(rr$analysis_object_path[1])
  counts_all <- obj$counts
  md <- obj$metadata

  framework <- disp_all |>
    filter(project == proj) |>
    mutate(
      time_order_numeric = as.numeric(analysis_time_order),
      time_factor = as.character(time_factor),
      time_day = suppressWarnings(as.numeric(time_day))
    ) |>
    arrange(patient_id, time_order_numeric, run_id)

  # counts rows needed
  need_rows <- unique(c(framework$run_id, framework$baseline_run_id))
  missing_runs <- setdiff(need_rows, rownames(counts_all))
  if (length(missing_runs)) {
    stop(proj, ": run_ids missing from RDS counts: ",
         paste(head(missing_runs, 10), collapse = ", "))
  }

  cnt <- counts_all[framework$run_id, , drop = FALSE]

  # ---- Bray-Curtis recomputation (QC vs frozen) ----
  rel <- cnt / rowSums(cnt)
  bray_mat <- as.matrix(vegan::vegdist(rel, method = "bray"))
  framework$bray_recomputed <- vapply(
    seq_len(nrow(framework)),
    function(i) bray_mat[framework$run_id[i], framework$baseline_run_id[i]],
    numeric(1))

  bray_qc <- framework |>
    transmute(
      project = proj,
      run_id,
      bray_frozen = as.numeric(bray_from_patient_baseline),
      bray_recomputed,
      abs_diff = abs(as.numeric(bray_from_patient_baseline) - bray_recomputed),
      max_abs_diff = max(abs(as.numeric(bray_from_patient_baseline) - bray_recomputed)),
      n_mismatch = sum(abs(as.numeric(bray_from_patient_baseline) - bray_recomputed) > 1e-9)
    ) |>
    slice(1) |>
    select(-run_id)

  if (bray_qc$n_mismatch[1] > 0) {
    stop(proj, ": Bray-Curtis recomputation mismatch vs Step88A2 frozen values.")
  }
  bray_qc_rows[[length(bray_qc_rows) + 1]] <- bray_qc

  framework$bray_from_patient_baseline <- as.numeric(framework$bray_from_patient_baseline)
  framework$is_patient_reference <- framework$run_id == framework$baseline_run_id

  if (any(abs(framework$bray_from_patient_baseline[framework$is_patient_reference]) > 1e-12)) {
    stop(proj, ": frozen baseline Bray distances are not zero.")
  }

  # ---- Aitchison displacement (primary CZM) ----
  framework$aitchison_from_patient_baseline_CZM <- compute_aitchison_displacement(cnt, framework)

  # Zero-handling sensitivities (fixed before looking at results)
  imp05 <- zero_impute(cnt, method = "PC0.5")
  clr05 <- clr_matrix(imp05)
  dst05 <- as.matrix(dist(clr05, method = "euclidean"))
  framework$aitchison_from_patient_baseline_PC0.5 <- vapply(
    seq_len(nrow(framework)),
    function(i) dst05[framework$run_id[i], framework$baseline_run_id[i]],
    numeric(1))

  imp1 <- zero_impute(cnt, method = "PC1")
  clr1 <- clr_matrix(imp1)
  dst1 <- as.matrix(dist(clr1, method = "euclidean"))
  framework$aitchison_from_patient_baseline_PC1 <- vapply(
    seq_len(nrow(framework)),
    function(i) dst1[framework$run_id[i], framework$baseline_run_id[i]],
    numeric(1))

  # baseline = 0 QC (all schemes)
  for (metric in c("aitchison_from_patient_baseline_CZM",
                   "aitchison_from_patient_baseline_PC0.5",
                   "aitchison_from_patient_baseline_PC1")) {
    if (any(abs(framework[[metric]][framework$is_patient_reference]) > 1e-9)) {
      stop(proj, ": baseline Aitchison distance not zero for ", metric)
    }
  }

  write_csv_safe(
    framework |>
      select(project, run_id, patient_id, time_label, time_order_numeric,
             time_factor, is_patient_reference, baseline_run_id,
             bray_from_patient_baseline,
             aitchison_from_patient_baseline_CZM,
             aitchison_from_patient_baseline_PC0.5,
             aitchison_from_patient_baseline_PC1),
    file.path(OUT, "01_SAMPLE_DISPLACEMENT",
              paste0(proj, "_within_patient_displacement_bray_aitchison.csv")))

  sample_rows[[length(sample_rows) + 1]] <- framework |>
    select(project, run_id, patient_id, time_label, time_order_numeric,
           time_factor, is_patient_reference, baseline_run_id,
           bray_from_patient_baseline,
           aitchison_from_patient_baseline_CZM,
           aitchison_from_patient_baseline_PC0.5,
           aitchison_from_patient_baseline_PC1)

  # ---- Timepoint summaries (Aitchison CZM primary + bray for reference) ----
  summ <- framework |>
    group_by(project, analysis_role = rr$analysis_role[1],
             time_order_numeric, time_label) |>
    summarise(
      n_samples = n(),
      n_patients = n_distinct(patient_id),
      bray_median = median(bray_from_patient_baseline, na.rm = TRUE),
      bray_mean = mean(bray_from_patient_baseline, na.rm = TRUE),
      bray_sd = sd(bray_from_patient_baseline, na.rm = TRUE),
      aitchison_median_CZM = median(aitchison_from_patient_baseline_CZM, na.rm = TRUE),
      aitchison_q1_CZM = quantile(aitchison_from_patient_baseline_CZM, 0.25, na.rm = TRUE),
      aitchison_q3_CZM = quantile(aitchison_from_patient_baseline_CZM, 0.75, na.rm = TRUE),
      aitchison_mean_CZM = mean(aitchison_from_patient_baseline_CZM, na.rm = TRUE),
      aitchison_sd_CZM = sd(aitchison_from_patient_baseline_CZM, na.rm = TRUE),
      aitchison_median_PC0.5 = median(aitchison_from_patient_baseline_PC0.5, na.rm = TRUE),
      aitchison_median_PC1 = median(aitchison_from_patient_baseline_PC1, na.rm = TRUE),
      .groups = "drop")

  write_csv_safe(summ, file.path(
    OUT, "02_TIMEPOINT_SUMMARIES",
    paste0(proj, "_displacement_timepoint_summary.csv")))

  summary_rows[[length(summary_rows) + 1]] <- summ

  # ---- Paired contrasts (Bray recomputed & Aitchison) ----
  cr <- contrast_registry |> filter(project == proj)

  # patient-level paired differences are also exported for meta-analysis
  pair_df <- function(metric_col, scheme) {
    framework |>
      filter(time_factor %in% c(cr$early_label, cr$late_label)) |>
      select(patient_id, time_factor, val = all_of(metric_col)) |>
      distinct(patient_id, time_factor, .keep_all = TRUE) |>
      pivot_wider(names_from = time_factor, values_from = val) |>
      filter(!is.na(.data[[cr$early_label]]), !is.na(.data[[cr$late_label]])) |>
      transmute(
        project = proj,
        patient_id,
        scheme = scheme,
        early_time_label = cr$early_label,
        late_time_label = cr$late_label,
        early_disp = .data[[cr$early_label]],
        late_disp = .data[[cr$late_label]],
        paired_diff = late_disp - early_disp
      )
  }

  # all-paired patient diffs per metric/scheme
  diff_specs <- tribble(
    ~metric_col, ~scheme,
    "bray_from_patient_baseline", "bray",
    "aitchison_from_patient_baseline_CZM", "aitchison_CZM",
    "aitchison_from_patient_baseline_PC0.5", "aitchison_PC0.5",
    "aitchison_from_patient_baseline_PC1", "aitchison_PC1"
  )

  for (k in seq_len(nrow(diff_specs))) {
    pdf <- pair_df(diff_specs$metric_col[k], diff_specs$scheme[k])
    if (nrow(pdf)) {
      patient_diff_rows[[length(patient_diff_rows) + 1]] <- pdf
      write_csv_safe(pdf, file.path(
        OUT, "03_PAIRED_CONTRASTS", "PATIENT_LEVEL_DIFFERENCES",
        paste0(proj, "_patient_paired_differences_",
               diff_specs$scheme[k], ".csv")))
    }
  }

  for (k in seq_len(nrow(diff_specs))) {
    metric <- if (diff_specs$scheme[k] == "bray") {
      "bray_from_patient_baseline"
    } else {
      diff_specs$metric_col[k]
    }
    paired_rows[[length(paired_rows) + 1]] <- paired_effect(
      framework, proj, cr$early_label, cr$late_label, metric,
      subset_type = "ALL_PAIRED", required_anchor = cr$required_anchor)
    paired_rows[[length(paired_rows) + 1]] <- paired_effect(
      framework, proj, cr$early_label, cr$late_label, metric,
      subset_type = "COMMON_ANCHOR_SENSITIVITY",
      required_anchor = cr$required_anchor)
  }

  # ---- Spearman correlation Bray vs Aitchison (follow-up rows) ----
  fu <- framework |> filter(!is_patient_reference)
  sp <- tibble(
    project = proj,
    n_followup_samples = nrow(fu),
    spearman_bray_CZM = cor(fu$bray_from_patient_baseline,
                            fu$aitchison_from_patient_baseline_CZM,
                            method = "spearman"),
    spearman_bray_PC0.5 = cor(fu$bray_from_patient_baseline,
                              fu$aitchison_from_patient_baseline_PC0.5,
                              method = "spearman"),
    spearman_bray_PC1 = cor(fu$bray_from_patient_baseline,
                            fu$aitchison_from_patient_baseline_PC1,
                            method = "spearman"),
    spearman_CZM_PC1 = cor(fu$aitchison_from_patient_baseline_CZM,
                           fu$aitchison_from_patient_baseline_PC1,
                           method = "spearman")
  )
  spearman_rows[[length(spearman_rows) + 1]] <- sp

  # ---- Mixed models on follow-up displacement ----
  adj <- resolve_adjustment(proj, md)
  followup <- framework |>
    filter(!is_patient_reference)

  for (metric in c("bray_from_patient_baseline",
                   "aitchison_from_patient_baseline_CZM")) {
    m_f <- safe_lmer(followup, outcome = metric, time_mode = "factor",
                     include_depth = FALSE, adjustment = adj)
    model_status_rows[[length(model_status_rows) + 1]] <- tibble(
      project = proj, analysis_family = "WITHIN_PATIENT_DISPLACEMENT",
      outcome = metric, model_type = "FOLLOWUP_CATEGORICAL_TIME_MIXED_MODEL",
      status = m_f$status)
    if (m_f$status == "OK") {
      model_global_rows[[length(model_global_rows) + 1]] <- m_f$global |>
        mutate(project = proj, analysis_family = "WITHIN_PATIENT_DISPLACEMENT",
               model_type = "FOLLOWUP_CATEGORICAL_TIME_MIXED_MODEL", .before = 1)
      model_coef_rows[[length(model_coef_rows) + 1]] <- m_f$coefficients |>
        mutate(project = proj, analysis_family = "WITHIN_PATIENT_DISPLACEMENT",
               outcome = metric,
               model_type = "FOLLOWUP_CATEGORICAL_TIME_MIXED_MODEL", .before = 1)
    }

    if (nrow(followup) > 0 && all(!is.na(followup$time_day)) &&
        n_distinct(followup$time_day) >= 2) {
      m_n <- safe_lmer(followup, outcome = metric, time_mode = "numeric",
                       include_depth = FALSE, adjustment = adj)
      model_status_rows[[length(model_status_rows) + 1]] <- tibble(
        project = proj, analysis_family = "WITHIN_PATIENT_DISPLACEMENT",
        outcome = metric, model_type = "FOLLOWUP_CONTINUOUS_DAY_MIXED_MODEL",
        status = m_n$status)
      if (m_n$status == "OK") {
        model_global_rows[[length(model_global_rows) + 1]] <- m_n$global |>
          mutate(project = proj, analysis_family = "WITHIN_PATIENT_DISPLACEMENT",
                 model_type = "FOLLOWUP_CONTINUOUS_DAY_MIXED_MODEL", .before = 1)
        model_coef_rows[[length(model_coef_rows) + 1]] <- m_n$coefficients |>
          mutate(project = proj, analysis_family = "WITHIN_PATIENT_DISPLACEMENT",
                 outcome = metric,
                 model_type = "FOLLOWUP_CONTINUOUS_DAY_MIXED_MODEL", .before = 1)
      }
    }
  }

  ck(paste0("FINISHED ", proj))
}

# ------------------------------------------------------------
# Bind and write combined tables
# ------------------------------------------------------------

ck("COMBINING OUTPUTS")

sample_all <- bind_rows(sample_rows)
write_csv_safe(sample_all, file.path(
  OUT, "01_SAMPLE_DISPLACEMENT",
  "V2_98A_ALL_within_patient_displacement_bray_aitchison.csv"))

bray_qc_all <- bind_rows(bray_qc_rows)
write_csv_safe(bray_qc_all, file.path(
  OUT, "05_QC_AND_CORRELATION",
  "V2_98A_bray_recompute_vs_frozen_QC.csv"))

summary_all <- bind_rows(summary_rows)
write_csv_safe(summary_all, file.path(
  OUT, "02_TIMEPOINT_SUMMARIES",
  "V2_98A_ALL_displacement_timepoint_summary.csv"))

# Paired contrasts with FDR within (metric, subset_type)
paired_all <- bind_rows(paired_rows)
paired_all <- paired_all |>
  group_by(metric, subset_type) |>
  mutate(wilcoxon_fdr_within_metric = bh_adjust(wilcoxon_p)) |>
  ungroup() |>
  arrange(metric, subset_type, project)

write_csv_safe(paired_all, file.path(
  OUT, "03_PAIRED_CONTRASTS",
  "V2_98A_paired_early_vs_late_contrasts.csv"))

patient_diff_all <- bind_rows(patient_diff_rows)
write_csv_safe(patient_diff_all, file.path(
  OUT, "03_PAIRED_CONTRASTS", "PATIENT_LEVEL_DIFFERENCES",
  "V2_98A_ALL_patient_paired_differences.csv"))

# Mixed model tables
model_status_all <- bind_rows(model_status_rows)
write_csv_safe(model_status_all, file.path(
  OUT, "04_MIXED_MODELS", "V2_98A_displacement_model_status.csv"))

model_global_all <- bind_rows(model_global_rows)
write_csv_safe(model_global_all, file.path(
  OUT, "04_MIXED_MODELS", "V2_98A_displacement_model_global_time_tests.csv"))

model_coef_all <- bind_rows(model_coef_rows)
write_csv_safe(model_coef_all, file.path(
  OUT, "04_MIXED_MODELS", "V2_98A_displacement_model_coefficients.csv"))

# QC mirror of Bray models vs frozen Step88A2 global tests
bray_model_global <- model_global_all |>
  filter(outcome == "bray_from_patient_baseline")

frozen_bray_global <- frozen_global |>
  filter(analysis_family == "WITHIN_PATIENT_BRAY_DISPLACEMENT") |>
  mutate(outcome = "bray_from_patient_baseline",
         project = as.character(project)) |>
  select(project, model_type, n_observations, n_patients, time_mode,
         adjustment, include_depth, likelihood_ratio_chisq, global_time_p)

model_qc <- bray_model_global |>
  select(project, model_type, n_observations, n_patients,
         likelihood_ratio_chisq, global_time_p) |>
  left_join(frozen_bray_global,
            by = c("project", "model_type"),
            suffix = c("_recomputed", "_frozen")) |>
  mutate(
    chisq_diff = abs(likelihood_ratio_chisq_recomputed - likelihood_ratio_chisq_frozen),
    p_diff = abs(global_time_p_recomputed - global_time_p_frozen),
    model_matches = chisq_diff < 1e-6 & p_diff < 1e-8
  )

write_csv_safe(model_qc, file.path(
  OUT, "05_QC_AND_CORRELATION",
  "V2_98A_bray_mixed_model_recompute_vs_frozen_QC.csv"))

if (any(!model_qc$model_matches, na.rm = TRUE)) {
  stop("Bray mixed-model recomputation does not match frozen Step88A2 results.")
}

spearman_all <- bind_rows(spearman_rows)
write_csv_safe(spearman_all, file.path(
  OUT, "05_QC_AND_CORRELATION",
  "V2_98A_bray_aitchison_spearman_correlations.csv"))

# ------------------------------------------------------------
# Direction consistency registry (Bray vs Aitchison)
# ------------------------------------------------------------

dir_reg <- paired_all |>
  filter(subset_type == "ALL_PAIRED") |>
  select(project, metric, n_pairs, mean_difference, paired_effect_dz,
         wilcoxon_p, wilcoxon_fdr_within_metric, direction) |>
  pivot_wider(id_cols = project,
              names_from = metric,
              values_from = c(n_pairs, mean_difference, paired_effect_dz,
                              wilcoxon_p, wilcoxon_fdr_within_metric,
                              direction)) |>
  mutate(
    direction_consistent_bray_vs_aitchison_CZM = case_when(
      direction_bray_from_patient_baseline ==
        direction_aitchison_from_patient_baseline_CZM ~ TRUE,
      TRUE ~ FALSE)
  )

write_csv_safe(dir_reg, file.path(
  OUT, "05_QC_AND_CORRELATION",
  "V2_98A_bray_aitchison_direction_consistency_registry.csv"))

# ------------------------------------------------------------
# Zero-handling sensitivity table
# ------------------------------------------------------------

sens <- paired_all |>
  filter(subset_type == "ALL_PAIRED",
         metric %in% c("bray_from_patient_baseline",
                       "aitchison_from_patient_baseline_CZM",
                       "aitchison_from_patient_baseline_PC0.5",
                       "aitchison_from_patient_baseline_PC1")) |>
  select(project, metric, n_pairs, mean_difference, paired_effect_dz,
         wilcoxon_p, direction)

model_sens <- model_global_all |>
  filter(model_type == "FOLLOWUP_CATEGORICAL_TIME_MIXED_MODEL",
         outcome %in% c("bray_from_patient_baseline",
                        "aitchison_from_patient_baseline_CZM")) |>
  select(project, outcome, global_time_p) |>
  pivot_wider(id_cols = project, names_from = outcome, values_from = global_time_p)

sens <- sens |>
  left_join(model_sens, by = "project")

write_csv_safe(sens, file.path(
  OUT, "06_ZERO_HANDLING_SENSITIVITY",
  "V2_98A_zero_handling_sensitivity_summary.csv"))

# ------------------------------------------------------------
# Done marker / README
# ------------------------------------------------------------

ck("STEP98A COMPLETE")

cat("Step98A complete.\n",
    file = file.path(OUT, "_STEP98A_COMPLETE.ok"))

readme <- c(
  "SEPSIS V2 UPGRADE - STEP98A AITCHISON LONGITUDINAL ROBUSTNESS",
  paste0("Created: ", Sys.time()),
  "",
  "Framework: Step88A2 GE2 repeated-measures subset (653 samples / 258 patients),",
  "cohort-specific, no cross-cohort ASV merge.",
  "Baseline: each patient's earliest longitudinal sample (identical to Step88A2).",
  "Aitchison distance = Euclidean distance between CLR vectors.",
  "Zero handling (fixed before results): CZM multiplicative replacement (primary),",
  "pseudocount 0.5 and pseudocount 1 (sensitivity).",
  "Bray-Curtis recomputed from frozen RDS and matched to Step88A2 values (QC).",
  "Paired contrasts follow Step88B structure (ALL_PAIRED + COMMON_ANCHOR_SENSITIVITY).",
  "",
  "Output directories:",
  "  01_SAMPLE_DISPLACEMENT   per-sample bray + aitchison displacement",
  "  02_TIMEPOINT_SUMMARIES   per-timepoint summary",
  "  03_PAIRED_CONTRASTS      early-vs-late contrasts + patient-level differences",
  "  04_MIXED_MODELS          follow-up categorical/continuous day mixed models",
  "  05_QC_AND_CORRELATION    Bray recompute QC, model QC, Spearman, direction registry",
  "  06_ZERO_HANDLING_SENSITIVITY zero-handling comparison",
  "",
  "Frozen Step87B/88A2/88B outputs were only read; nothing was overwritten."
)

writeLines(readme, file.path(OUT, "README_STEP98A.txt"))

cat("DONE\n")
