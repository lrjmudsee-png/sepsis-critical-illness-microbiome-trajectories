# ============================================================
# Sepsis V2 Upgrade - Step 98B
# Formal random-effects meta-analysis of the within-patient
# longitudinal early-vs-late ecological displacement.
# ------------------------------------------------------------
# Design (frozen BEFORE inspecting any Step98A significance):
#   - eligibility: META_ELIGIBILITY_FREEZE.csv (design-based)
#   - primary set (natural-history): PRJNA691455, PRJNA851469,
#     PRJNA516701
#   - sensitivity set: primary + PRJEB82425 (supportive, labelled)
#   - excluded from pooling: PRJNA578267 (nonsepsis control),
#     PRJNA1166732/PRJNA430161 (intervention support), PRJNA978257
#     (static) - reported as separate cohort evidence only
#   - Bray-Curtis and Aitchison are NEVER pooled together
#   - patient-level early-vs-late paired differences (from Step98A)
#   - effect: Cohen dz and small-sample corrected Hedges g_z
#   - pooling: metafor::rma REML with Hartung-Knapp (HKSJ) CIs
#   - sensitivity: leave-one-cohort-out + patient-level bootstrap
#   - prediction interval explicitly labelled exploratory
# ============================================================

options(stringsAsFactors = FALSE)

required_pkgs <- c(
  "readr", "dplyr", "tidyr", "purrr", "tibble", "stringr", "metafor"
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
  library(dplyr)
  library(tidyr)
  library(metafor)
})

set.seed(20260907)

ROOT <- "E:/sepsis_project"
OUT <- file.path(
  ROOT, "results", "V2_UPGRADE_20260907",
  "98B_RANDOM_EFFECTS_META_ANALYSIS"
)

STEP98A_OUT <- file.path(
  ROOT, "results", "V2_UPGRADE_20260907",
  "98A_AITCHISON_LONGITUDINAL_ROBUSTNESS"
)

dirs <- c(
  OUT,
  file.path(OUT, "01_COHORT_EFFECT_TABLES"),
  file.path(OUT, "02_POOLED_RESULTS"),
  file.path(OUT, "03_FOREST"),
  file.path(OUT, "04_LEAVE_ONE_OUT"),
  file.path(OUT, "05_BOOTSTRAP")
)
for (d in dirs) dir.create(d, recursive = TRUE, showWarnings = FALSE)

LOG <- file.path(OUT, "_STEP98B_runtime_checkpoints.txt")
ERR <- file.path(OUT, "_STEP98B_FATAL_ERROR.txt")
if (file.exists(LOG)) unlink(LOG)
if (file.exists(ERR)) unlink(ERR)

ck <- function(x) cat(paste0(x, ": ", Sys.time(), "\n"), file = LOG, append = TRUE)
write_csv_safe <- function(x, p) write_excel_csv(x, p, na = "")
safe_csv <- function(p) {
  suppressMessages(read_csv(p, show_col_types = FALSE, progress = FALSE,
                            name_repair = "unique"))
}

# hedges small-sample correction for paired d with df = n-1
J_paired <- function(n) 1 - 3 / (4 * (n - 1) - 1)

cohort_effect <- function(pd, project, scheme) {
  d <- pd |>
    filter(project == .env$project, scheme == .env$scheme) |>
    pull(paired_diff)
  n <- length(d)
  if (n < 2) return(NULL)
  m <- mean(d)
  s <- sd(d)
  dz <- if (s > 0) m / s else NA_real_
  gz <- if (!is.na(dz)) J_paired(n) * dz else NA_real_
  se_dz <- sqrt(1 / n + dz^2 / (2 * n))
  se_g <- J_paired(n) * se_dz
  tstat <- if (s > 0) m / (s / sqrt(n)) else NA_real_
  p_t <- if (!is.na(tstat)) 2 * pt(-abs(tstat), df = n - 1) else NA_real_
  tcrit <- qt(0.975, df = n - 1)
  tibble(
    project = project, scheme = scheme, n_pairs = n,
    mean_paired_diff = m, sd_paired_diff = s,
    median_paired_diff = median(d),
    cohen_dz = dz, hedges_gz = gz,
    se_gz = se_g,
    gz_ci_low = gz - tcrit * se_g,
    gz_ci_high = gz + tcrit * se_g,
    paired_t_p = p_t,
    wilcoxon_p = tryCatch(suppressWarnings(
      wilcox.test(d, mu = 0, paired = FALSE, exact = FALSE,
                  correct = FALSE)$p.value),
      error = function(e) NA_real_),
    direction = case_when(m > 0 ~ "INCREASE", m < 0 ~ "DECREASE",
                          TRUE ~ "NO_CHANGE")
  )
}

ck("STEP98B STARTED")

stopifnot(file.exists(file.path(OUT, "META_ELIGIBILITY_FREEZE.csv")))

elig <- safe_csv(file.path(OUT, "META_ELIGIBILITY_FREEZE.csv"))

# ---- read Step98A patient-level paired differences ----
pd_all_path <- file.path(
  STEP98A_OUT, "03_PAIRED_CONTRASTS", "PATIENT_LEVEL_DIFFERENCES",
  "V2_98A_ALL_patient_paired_differences.csv")
stopifnot(file.exists(pd_all_path))
pd_all <- safe_csv(pd_all_path)

schemes <- c("bray", "aitchison_CZM", "aitchison_PC0.5", "aitchison_PC1")
primary_set <- c("PRJNA691455", "PRJNA851469", "PRJNA516701")
sensitivity_extra <- c("PRJEB82425")
support_cohorts <- c("PRJNA578267", "PRJNA1166732", "PRJNA430161")

pooled_rows <- list()
cohort_rows <- list()
loo_rows <- list()
boot_rows <- list()

pool_one <- function(cohorts, label, scheme, pd) {
  ce <- map_dfr(cohorts, function(proj) cohort_effect(pd, proj, scheme)) |>
    filter(!is.na(hedges_gz))
  if (nrow(ce) < 2) {
    return(list(tab = NULL, res = NULL, ce = ce))
  }
  res <- tryCatch(
    metafor::rma(
      yi = hedges_gz, sei = se_gz, data = ce,
      method = "REML", test = "knha", weighted = TRUE
    ),
    error = function(e) e
  )
  if (inherits(res, "error")) {
    res <- tryCatch(
      metafor::rma(yi = hedges_gz, sei = se_gz, data = ce,
                   method = "REML", weighted = TRUE),
      error = function(e) e
    )
  }
  if (inherits(res, "error")) {
    return(list(tab = NULL, res = NULL, ce = ce))
  }
  pred <- tryCatch(predict(res), error = function(e) NULL)
  tab <- tibble(
    pool_label = label, scheme = scheme, k_studies = res$k,
    pooled_hedges_gz = as.numeric(res$b),
    se = res$se,
    ci_low = res$ci.lb, ci_high = res$ci.ub,
    tau2 = res$tau2, tau = res$tau,
    I2 = ifelse(is.null(res$I2), NA_real_, res$I2),
    H2 = ifelse(is.null(res$H2), NA_real_, res$H2),
    Q = res$QE, Q_df = res$k - 1, Q_p = res$QEp,
    test = res$test,
    pred_low = if (!is.null(pred)) pred$pi.lb else NA_real_,
    pred_high = if (!is.null(pred)) pred$pi.ub else NA_real_,
    pred_labelled = "EXPLORATORY"
  )
  list(tab = tab, res = res, ce = ce)
}

run_scheme <- function(scheme) {
  pd <- pd_all |> filter(scheme == .env$scheme)

  # Cohort-level evidence table (ALL cohorts; no significance filter)
  all_proj <- elig$project
  ce_all <- map_dfr(all_proj, function(proj) cohort_effect(pd, proj, scheme)) |>
    left_join(
      elig |> select(project, analysis_role, primary_meta_eligible,
                     sensitivity_meta_eligible, exclusion_reason),
      by = "project")
  write_csv_safe(ce_all, file.path(
    OUT, "01_COHORT_EFFECT_TABLES",
    paste0("V2_98B_cohort_effect_table_", scheme, ".csv")))

  # ---- primary pool ----
  p1 <- pool_one(primary_set, "PRIMARY_META", scheme, pd)
  if (!is.null(p1$tab)) pooled_rows[[length(pooled_rows) + 1]] <<- p1$tab
  if (!is.null(p1$res)) {
    pdf(file.path(OUT, "03_FOREST",
                  paste0("V2_98B_forest_PRMARY_", scheme, ".pdf")))
    metafor::forest(p1$res, main = paste("Primary meta-analysis:", scheme),
                    xlab = "Hedges g_z (early-to-late paired displacement)",
                    transf = identity, cex = 1)
    dev.off()
  }

  # ---- sensitivity pool (primary + PRJEB82425) ----
  s1 <- pool_one(c(primary_set, sensitivity_extra),
                 "SENSITIVITY_INCLUDING_PRJEB82425", scheme, pd)
  if (!is.null(s1$tab)) pooled_rows[[length(pooled_rows) + 1]] <<- s1$tab
  if (!is.null(s1$res)) {
    pdf(file.path(OUT, "03_FOREST",
                  paste0("V2_98B_forest_SENSITIVITY_", scheme, ".pdf")))
    metafor::forest(s1$res,
                    main = paste("Sensitivity (incl. PRJEB82425):", scheme),
                    xlab = "Hedges g_z", cex = 1)
    dev.off()
  }

  # ---- leave-one-cohort-out (primary set) ----
  loo <- map_dfr(primary_set, function(drop) {
    kp <- setdiff(primary_set, drop)
    z <- pool_one(kp, paste0("LEAVE_OUT_", drop), scheme, pd)
    if (!is.null(z$tab)) {
      z$tab |> mutate(removed_cohort = drop, remaining = paste(kp, collapse = "+"))
    } else {
      tibble()
    }
  })
  if (nrow(loo)) {
    write_csv_safe(loo, file.path(
      OUT, "04_LEAVE_ONE_OUT",
      paste0("V2_98B_leave_one_out_", scheme, ".csv")))
    loo_rows[[length(loo_rows) + 1]] <<- loo
  }

  # ---- patient-level bootstrap ----
  B <- 1000
  boot_est <- numeric(B)
  boot_se <- numeric(B)
  boot_p <- numeric(B)
  pd_sub <- pd |> filter(project %in% primary_set)
  proj_vec <- pd_sub$project
  pat_vec <- pd_sub$patient_id
  diff_vec <- pd_sub$paired_diff
  # resample patients WITHIN cohort
  for (b in seq_len(B)) {
    bd <- map_dfr(primary_set, function(proj) {
      dsub <- pd_sub |> filter(project == proj)
      idx <- sample(seq_len(nrow(dsub)), nrow(dsub), replace = TRUE)
      dsub[idx, ]
    })
    ce_b <- map_dfr(primary_set, function(proj) cohort_effect(bd, proj, scheme)) |>
      filter(!is.na(hedges_gz))
    if (nrow(ce_b) < 2) {
      boot_est[b] <- NA_real_
      boot_se[b] <- NA_real_
      boot_p[b] <- NA_real_
      next
    }
    rb <- tryCatch(
      metafor::rma(yi = hedges_gz, sei = se_gz, data = ce_b,
                   method = "REML", test = "knha"),
      error = function(e) NULL)
    if (is.null(rb)) {
      rb <- tryCatch(
        metafor::rma(yi = hedges_gz, sei = se_gz, data = ce_b,
                     method = "REML"),
        error = function(e) NULL)
    }
    if (is.null(rb)) {
      boot_est[b] <- NA_real_; boot_se[b] <- NA_real_; boot_p[b] <- NA_real_
    } else {
      boot_est[b] <- as.numeric(rb$b)
      boot_se[b] <- rb$se
      boot_p[b] <- rb$pval
    }
  }
  ok <- !is.na(boot_est)
  boot_tab <- tibble(
    scheme = scheme,
    n_boot = B,
    n_valid = sum(ok),
    boot_median = median(boot_est[ok]),
    boot_mean = mean(boot_est[ok]),
    boot_ci_low = quantile(boot_est[ok], 0.025),
    boot_ci_high = quantile(boot_est[ok], 0.975),
    frac_positive = mean(boot_est[ok] > 0),
    frac_p_lt_0.05 = mean(boot_p[ok] < 0.05, na.rm = TRUE),
    min_boot_est = min(boot_est[ok]),
    max_boot_est = max(boot_est[ok])
  )
  write_csv_safe(boot_tab, file.path(
    OUT, "05_BOOTSTRAP",
    paste0("V2_98B_patient_bootstrap_", scheme, ".csv")))
  boot_rows[[length(boot_rows) + 1]] <<- boot_tab
  ck(paste0("SCHEME DONE ", scheme))
}

for (s in schemes) run_scheme(s)

# Combined pooled table
pooled_all <- bind_rows(pooled_rows)
write_csv_safe(pooled_all, file.path(
  OUT, "02_POOLED_RESULTS", "V2_98B_pooled_results_all_schemes.csv"))

# Human-readable quick view of the primary pools
view <- pooled_all |>
  filter(grepl("PRIMARY_META", pool_label)) |>
  select(pool_label, scheme, k_studies, pooled_hedges_gz, ci_low, ci_high,
         tau2, I2, Q, Q_p)
write_csv_safe(view, file.path(
  OUT, "02_POOLED_RESULTS", "V2_98B_PRIMARY_pooled_summary.csv"))

# ---- evidence summary (direction registry incl. meta pool) ----
cat("Step98B complete.\n", file = file.path(OUT, "_STEP98B_COMPLETE.ok"))

readme <- c(
  "SEPSIS V2 UPGRADE - STEP98B RANDOM-EFFECTS META-ANALYSIS",
  paste0("Created: ", Sys.time()),
  "",
  "Eligibility: META_ELIGIBILITY_FREEZE.csv (design-based, frozen before",
  "inspecting Step98A significance).",
  "PRIMARY pool: PRJNA691455 + PRJNA851469 + PRJNA516701 (natural history).",
  "SENSITIVITY pool: primary + PRJEB82425 (supportive only, labelled).",
  "Never pooled: PRJNA578267 (nonsepsis control), PRJNA1166732/PRJNA430161",
  "(intervention support) - see cohort effect tables.",
  "Metrics: Bray-Curtis and Aitchison run separately (never mixed).",
  "Effect: cohort Hedges g_z (small-sample corrected paired d);",
  "pooling via metafor::rma REML + Hartung-Knapp (test='knha');",
  "prediction intervals marked EXPLORATORY.",
  "Sensitivity: leave-one-cohort-out and patient-level block bootstrap.",
  "Inputs: Step98A patient-level paired differences + Step87B registry roles."
)
writeLines(readme, file.path(OUT, "README_STEP98B.txt"))

cat("DONE\n")
