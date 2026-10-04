# ============================================================
# Sepsis V2 - Step94B2B
# CRA002354 CORRECTED FORMAL SOURCE x LONGITUDINAL ECOLOGY
#
# Requires Step94B1D:
# - true non-chimeric OTU97 counts
# - corrected OTU relative abundance
# - primary alpha rarefied to 4000 reads
# - sensitivity alpha rarefied to 2000 reads
#
# This analysis intentionally does NOT require genus taxonomy.
# Taxonomy-ID repair is handled separately in Step94B1E.
# ============================================================

options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(stringr)
  library(tibble)
  library(lme4)
})

if (!requireNamespace("vegan", quietly = TRUE)) {
  stop("R package 'vegan' is required.")
}

ROOT <- "E:/sepsis_project"

DATA_DIR <- file.path(
  ROOT, "data", "CRA002354", "03_vsearch97_silva1382"
)

B1D_OUT <- file.path(
  ROOT, "results",
  "V2_34B1D_CRA002354_NONCHIMERIC_OTU_TABLE_AND_RAREFACTION_FIX"
)

SPEC_DIR <- file.path(
  ROOT, "results",
  "V2_34B2A3_SOURCE_MODEL_SPEC_QC_FIX"
)

OUT <- file.path(
  ROOT, "results",
  "V2_34B2B_CRA002354_CORRECTED_FORMAL_SOURCE_LONGITUDINAL_ANALYSIS"
)

dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------------------
# 0. Hard gates
# ------------------------------------------------------------

ready_file <- file.path(
  B1D_OUT,
  "06_STEP94B1D_CORRECTED_ANALYSIS_READINESS.csv"
)

if (!file.exists(ready_file)) {
  stop("Step94B1D corrected readiness file is missing.")
}

ready <- read_csv(ready_file, show_col_types = FALSE)

if (
  !"ready_for_corrected_step94B2" %in% names(ready) ||
  !isTRUE(ready$ready_for_corrected_step94B2[1])
) {
  stop("Step94B1D is not ready for corrected Step94B2B.")
}

if (!file.exists(file.path(SPEC_DIR, "_STEP94B2A3_COMPLETE.ok"))) {
  stop("Step94B2A3 specification freeze/QC is incomplete.")
}

REL_FILE <- file.path(
  DATA_DIR,
  "CRA002354_OTU97_relative_abundance_NONCHIMERIC_CORRECTED.csv"
)

ALPHA4000_FILE <- file.path(
  DATA_DIR,
  "CRA002354_OTU97_alpha_PRIMARY_rarefied4000.csv"
)

META4000_FILE <- file.path(
  DATA_DIR,
  "CRA002354_analysis_metadata_PRIMARY_NONCHIMERIC_min4000.csv"
)

RARE4000_FILE <- file.path(
  DATA_DIR,
  "CRA002354_OTU97_rarefied_PRIMARY_4000.csv"
)

ALPHA2000_FILE <- file.path(
  DATA_DIR,
  "CRA002354_OTU97_alpha_SENSITIVITY_rarefied2000.csv"
)

META2000_FILE <- file.path(
  DATA_DIR,
  "CRA002354_analysis_metadata_SENSITIVITY_NONCHIMERIC_min2000.csv"
)

needed <- c(
  REL_FILE, ALPHA4000_FILE, META4000_FILE, RARE4000_FILE,
  ALPHA2000_FILE, META2000_FILE
)

missing <- needed[!file.exists(needed)]
if (length(missing) > 0) {
  stop(paste("Missing corrected B1D data:", paste(missing, collapse = "; ")))
}

# ------------------------------------------------------------
# Helpers
# ------------------------------------------------------------

bray_two <- function(x, y) {
  den <- sum(x + y)
  if (den <= 0) return(NA_real_)
  sum(abs(x - y)) / den
}

rank_biserial <- function(x, g) {
  keep <- !is.na(x) & !is.na(g)
  x <- x[keep]
  g <- g[keep]

  a <- x[g == "PULMONARY"]
  b <- x[g == "NONPULMONARY_RECORDED"]

  if (length(a) == 0 || length(b) == 0) {
    return(tibble(
      n_pulmonary = length(a),
      n_nonpulmonary = length(b),
      median_pulmonary = NA_real_,
      median_nonpulmonary = NA_real_,
      median_difference = NA_real_,
      rank_biserial = NA_real_,
      p_value = NA_real_
    ))
  }

  wt <- suppressWarnings(wilcox.test(a, b, exact = FALSE))
  U <- unname(wt$statistic)

  tibble(
    n_pulmonary = length(a),
    n_nonpulmonary = length(b),
    median_pulmonary = median(a),
    median_nonpulmonary = median(b),
    median_difference = median(a) - median(b),
    rank_biserial = 2 * U / (length(a) * length(b)) - 1,
    p_value = wt$p.value
  )
}

fit_interaction <- function(
  df,
  outcome,
  timevar,
  extra_covars = character(),
  analysis_label = "PRIMARY"
) {

  vars <- unique(c(
    outcome, timevar, "pulmonary_binary", "patient_id", extra_covars
  ))

  d <- df %>%
    select(all_of(vars)) %>%
    filter(if_all(everything(), ~ !is.na(.)))

  if (nrow(d) < 20 || n_distinct(d$patient_id) < 10) {
    return(tibble(
      analysis = analysis_label,
      outcome = outcome,
      timevar = timevar,
      nobs = nrow(d),
      npatients = n_distinct(d$patient_id),
      beta_interaction = NA_real_,
      se_interaction = NA_real_,
      ci95_low = NA_real_,
      ci95_high = NA_real_,
      p_LRT = NA_real_,
      singular = NA,
      status = "INSUFFICIENT_DATA"
    ))
  }

  d$pulmonary_binary <- factor(
    d$pulmonary_binary,
    levels = c("NONPULMONARY_RECORDED", "PULMONARY")
  )

  # Ensure categorical covariates are factors.
  for (vv in extra_covars) {
    if (is.character(d[[vv]])) d[[vv]] <- factor(d[[vv]])
  }

  rhs_full <- paste0(timevar, " * pulmonary_binary")
  rhs_red  <- paste0(timevar, " + pulmonary_binary")

  if (length(extra_covars) > 0) {
    ec <- paste(extra_covars, collapse = " + ")
    rhs_full <- paste(rhs_full, ec, sep = " + ")
    rhs_red  <- paste(rhs_red, ec, sep = " + ")
  }

  f_full <- as.formula(
    paste(outcome, "~", rhs_full, "+ (1|patient_id)")
  )

  f_red <- as.formula(
    paste(outcome, "~", rhs_red, "+ (1|patient_id)")
  )

  full <- tryCatch(
    lmer(
      f_full,
      data = d,
      REML = FALSE,
      control = lmerControl(optimizer = "bobyqa")
    ),
    error = function(e) NULL
  )

  red <- tryCatch(
    lmer(
      f_red,
      data = d,
      REML = FALSE,
      control = lmerControl(optimizer = "bobyqa")
    ),
    error = function(e) NULL
  )

  if (is.null(full) || is.null(red)) {
    return(tibble(
      analysis = analysis_label,
      outcome = outcome,
      timevar = timevar,
      nobs = nrow(d),
      npatients = n_distinct(d$patient_id),
      beta_interaction = NA_real_,
      se_interaction = NA_real_,
      ci95_low = NA_real_,
      ci95_high = NA_real_,
      p_LRT = NA_real_,
      singular = NA,
      status = "MODEL_ERROR"
    ))
  }

  cn <- names(fixef(full))
  target1 <- paste0(timevar, ":pulmonary_binaryPULMONARY")
  target2 <- paste0("pulmonary_binaryPULMONARY:", timevar)

  int_name <- cn[cn %in% c(target1, target2)]

  if (length(int_name) != 1) {
    int_name <- cn[
      grepl(timevar, cn, fixed = TRUE) &
      grepl("pulmonary_binaryPULMONARY", cn, fixed = TRUE)
    ]
  }

  beta <- if (length(int_name)) unname(fixef(full)[int_name[1]]) else NA_real_
  se <- if (length(int_name)) unname(sqrt(diag(vcov(full)))[int_name[1]]) else NA_real_

  lrt <- anova(red, full)
  p <- lrt$`Pr(>Chisq)`[2]

  tibble(
    analysis = analysis_label,
    outcome = outcome,
    timevar = timevar,
    nobs = nrow(d),
    npatients = n_distinct(d$patient_id),
    beta_interaction = beta,
    se_interaction = se,
    ci95_low = beta - 1.96 * se,
    ci95_high = beta + 1.96 * se,
    p_LRT = p,
    singular = isSingular(full, tol = 1e-4),
    status = "OK"
  )
}

prepare_rel_matrix <- function(file, keep_ids) {
  x <- read_csv(file, show_col_types = FALSE, name_repair = "minimal")
  feats <- setdiff(names(x), "Run_ID")
  M <- as.matrix(x[, feats, drop = FALSE])
  storage.mode(M) <- "numeric"
  rownames(M) <- x$Run_ID
  M <- M[intersect(keep_ids, rownames(M)), , drop = FALSE]
  rs <- rowSums(M)
  if (any(rs <= 0)) stop("Zero-sum OTU relative-abundance sample found.")
  M / rs
}

prepare_rarefied_matrix <- function(file, keep_ids) {
  x <- read_csv(file, show_col_types = FALSE, name_repair = "minimal")
  feats <- setdiff(names(x), "Run_ID")
  M <- as.matrix(x[, feats, drop = FALSE])
  storage.mode(M) <- "numeric"
  rownames(M) <- x$Run_ID
  M <- M[intersect(keep_ids, rownames(M)), , drop = FALSE]
  M / rowSums(M)
}

add_personal_time_and_bray <- function(meta, M) {

  meta <- meta %>%
    filter(Run_ID %in% rownames(M)) %>%
    arrange(match(Run_ID, rownames(M))) %>%
    group_by(patient_id) %>%
    mutate(
      personal_baseline_day = min(time_day, na.rm = TRUE),
      days_since_personal_baseline = time_day - personal_baseline_day,
      is_personal_baseline = time_day == personal_baseline_day
    ) %>%
    ungroup()

  M <- M[meta$Run_ID, , drop = FALSE]

  meta$Bray_from_personal_baseline <- NA_real_

  for (pid in unique(meta$patient_id)) {
    ii <- which(meta$patient_id == pid)
    base_i <- ii[which.min(meta$time_day[ii])]
    x0 <- M[meta$Run_ID[base_i], , drop = TRUE]

    for (jj in ii) {
      meta$Bray_from_personal_baseline[jj] <- bray_two(
        x0,
        M[meta$Run_ID[jj], , drop = TRUE]
      )
    }
  }

  list(meta = meta, M = M)
}

# ------------------------------------------------------------
# 1. Primary corrected dataset
# ------------------------------------------------------------

meta4000 <- read_csv(META4000_FILE, show_col_types = FALSE)
alpha4000 <- read_csv(ALPHA4000_FILE, show_col_types = FALSE)

meta4000 <- meta4000 %>%
  select(-any_of(c("Observed_OTU97","Shannon_OTU97","Simpson_OTU97","rarefied_depth"))) %>%
  left_join(
    alpha4000 %>%
      select(
        Run_ID,
        rarefied_depth,
        Observed_OTU97,
        Shannon_OTU97,
        Simpson_OTU97
      ),
    by = "Run_ID"
  )

Mrel4000 <- prepare_rel_matrix(
  REL_FILE,
  meta4000$Run_ID
)

prim <- add_personal_time_and_bray(meta4000, Mrel4000)
dat <- prim$meta
Mrel4000 <- prim$M

write_csv(
  dat,
  file.path(OUT, "01_PRIMARY_CORRECTED_ANALYSIS_DATA.csv")
)

# ------------------------------------------------------------
# 2. Sample/patient QC
# ------------------------------------------------------------

qc <- dat %>%
  group_by(pulmonary_binary) %>%
  summarise(
    samples = n(),
    patients = n_distinct(patient_id),
    patients_GE2 = n_distinct(patient_id[duplicated(patient_id) | duplicated(patient_id, fromLast = TRUE)]),
    baseline_samples = sum(is_personal_baseline),
    followup_samples = sum(days_since_personal_baseline > 0),
    .groups = "drop"
  )

write_csv(
  qc,
  file.path(OUT, "02_PRIMARY_SAMPLE_PATIENT_QC.csv")
)

# ------------------------------------------------------------
# 3. Baseline alpha
# ------------------------------------------------------------

base <- dat %>% filter(is_personal_baseline)

alpha_metrics <- c(
  "Observed_OTU97",
  "Shannon_OTU97",
  "Simpson_OTU97"
)

baseline_alpha <- bind_rows(
  lapply(alpha_metrics, function(mm) {
    rr <- rank_biserial(base[[mm]], base$pulmonary_binary)
    rr$metric <- mm
    rr
  })
) %>%
  mutate(FDR = p.adjust(p_value, method = "BH")) %>%
  select(metric, everything())

write_csv(
  baseline_alpha,
  file.path(OUT, "03_PRIMARY_BASELINE_ALPHA_SOURCE_CONTRAST.csv")
)

# ------------------------------------------------------------
# 4. Baseline Bray PERMANOVA + dispersion
# ------------------------------------------------------------

Mbase <- Mrel4000[base$Run_ID, , drop = FALSE]
d_base <- vegan::vegdist(Mbase, method = "bray")

set.seed(20260824)
perm <- vegan::adonis2(
  d_base ~ pulmonary_binary,
  data = base,
  permutations = 9999
)

bd <- vegan::betadisper(d_base, base$pulmonary_binary)

set.seed(20260824)
bdp <- vegan::permutest(bd, permutations = 9999)

write_csv(
  tibble(
    n = nrow(base),
    pseudo_F = perm$F[1],
    R2 = perm$R2[1],
    p_value = perm$`Pr(>F)`[1]
  ),
  file.path(OUT, "04_PRIMARY_BASELINE_BRAY_PERMANOVA.csv")
)

write_csv(
  tibble(
    n = nrow(base),
    dispersion_F = bdp$tab$F[1],
    dispersion_p = bdp$tab$`Pr(>F)`[1]
  ),
  file.path(OUT, "05_PRIMARY_BASELINE_BRAY_DISPERSION.csv")
)

# ------------------------------------------------------------
# 5. Primary longitudinal alpha interactions
# ------------------------------------------------------------

alpha_models <- bind_rows(
  lapply(alpha_metrics, function(mm) {
    fit_interaction(
      dat,
      outcome = mm,
      timevar = "days_since_personal_baseline",
      analysis_label = "PRIMARY_RAREFIED4000"
    )
  })
) %>%
  mutate(FDR = p.adjust(p_LRT, method = "BH"))

write_csv(
  alpha_models,
  file.path(OUT, "06_PRIMARY_ALPHA_LONGITUDINAL_INTERACTIONS.csv")
)

# ------------------------------------------------------------
# 6. Primary Bray displacement interaction
# Exclude structural baseline zeros from formal inference.
# ------------------------------------------------------------

follow <- dat %>%
  filter(days_since_personal_baseline > 0)

bray_primary <- fit_interaction(
  follow,
  outcome = "Bray_from_personal_baseline",
  timevar = "days_since_personal_baseline",
  analysis_label = "PRIMARY_CORRECTED_RELATIVE_ABUNDANCE"
)

write_csv(
  bray_primary,
  file.path(OUT, "07_PRIMARY_BRAY_DISPLACEMENT_INTERACTION.csv")
)

# ------------------------------------------------------------
# 7. Baseline-to-latest supportive contrast
# ------------------------------------------------------------

latest <- dat %>%
  group_by(patient_id) %>%
  filter(days_since_personal_baseline == max(days_since_personal_baseline)) %>%
  slice_tail(n = 1) %>%
  ungroup() %>%
  filter(days_since_personal_baseline > 0)

latest_rr <- rank_biserial(
  latest$Bray_from_personal_baseline,
  latest$pulmonary_binary
) %>%
  mutate(
    estimand = "Personal-baseline to latest observed Bray displacement"
  )

write_csv(
  latest_rr,
  file.path(OUT, "08_SUPPORTIVE_BASELINE_TO_LATEST_BRAY.csv")
)

# ------------------------------------------------------------
# 8. Absolute ICU day sensitivity
# ------------------------------------------------------------

abs_models <- bind_rows(
  lapply(alpha_metrics, function(mm) {
    fit_interaction(
      dat,
      outcome = mm,
      timevar = "time_day",
      analysis_label = "ABSOLUTE_ICU_DAY"
    )
  }),
  fit_interaction(
    follow,
    outcome = "Bray_from_personal_baseline",
    timevar = "time_day",
    analysis_label = "ABSOLUTE_ICU_DAY"
  )
)

write_csv(
  abs_models,
  file.path(OUT, "09_SENSITIVITY_ABSOLUTE_ICU_DAY.csv")
)

# ------------------------------------------------------------
# 9. Separate severity-adjusted Bray sensitivities
# ------------------------------------------------------------

sev_specs <- list(
  SOFA = c("age","sex","baseline_sofa"),
  APACHE = c("age","sex","baseline_apache_ii"),
  LACTATE = c("age","sex","baseline_lactate")
)

sev <- bind_rows(
  lapply(names(sev_specs), function(nm) {
    fit_interaction(
      follow,
      outcome = "Bray_from_personal_baseline",
      timevar = "days_since_personal_baseline",
      extra_covars = sev_specs[[nm]],
      analysis_label = paste0("ADJUSTED_", nm)
    )
  })
)

write_csv(
  sev,
  file.path(OUT, "10_SENSITIVITY_SEVERITY_ADJUSTED_BRAY.csv")
)

# ------------------------------------------------------------
# 10. Exclude OTHER_UNKNOWN sensitivity
# ------------------------------------------------------------

known_follow <- follow %>%
  filter(infection_source_group != "OTHER_UNKNOWN")

known_model <- fit_interaction(
  known_follow,
  outcome = "Bray_from_personal_baseline",
  timevar = "days_since_personal_baseline",
  analysis_label = "EXCLUDE_OTHER_UNKNOWN"
)

write_csv(
  known_model,
  file.path(OUT, "11_SENSITIVITY_EXCLUDE_OTHER_UNKNOWN_BRAY.csv")
)

# ------------------------------------------------------------
# 11. Baseline first-sample <= ICU Day3 sensitivity
# ------------------------------------------------------------

base3 <- base %>%
  filter(personal_baseline_day <= 3)

base3_alpha <- bind_rows(
  lapply(alpha_metrics, function(mm) {
    rr <- rank_biserial(base3[[mm]], base3$pulmonary_binary)
    rr$metric <- mm
    rr
  })
) %>%
  mutate(FDR = p.adjust(p_value, method = "BH")) %>%
  select(metric, everything())

write_csv(
  base3_alpha,
  file.path(OUT, "12_SENSITIVITY_BASELINE_DAY3_ALPHA.csv")
)

if (nrow(base3) >= 20 && n_distinct(base3$pulmonary_binary) == 2) {
  M3 <- Mrel4000[base3$Run_ID, , drop = FALSE]
  d3 <- vegan::vegdist(M3, method = "bray")

  set.seed(20260824)
  p3 <- vegan::adonis2(
    d3 ~ pulmonary_binary,
    data = base3,
    permutations = 9999
  )

  bd3 <- vegan::betadisper(d3, base3$pulmonary_binary)

  set.seed(20260824)
  bd3p <- vegan::permutest(bd3, permutations = 9999)

  base3_beta <- tibble(
    n = nrow(base3),
    PERMANOVA_R2 = p3$R2[1],
    PERMANOVA_p = p3$`Pr(>F)`[1],
    dispersion_F = bd3p$tab$F[1],
    dispersion_p = bd3p$tab$`Pr(>F)`[1]
  )
} else {
  base3_beta <- tibble(
    n = nrow(base3),
    PERMANOVA_R2 = NA_real_,
    PERMANOVA_p = NA_real_,
    dispersion_F = NA_real_,
    dispersion_p = NA_real_
  )
}

write_csv(
  base3_beta,
  file.path(OUT, "13_SENSITIVITY_BASELINE_DAY3_BETA.csv")
)

# ------------------------------------------------------------
# 12. Rarefied-4000 Bray sensitivity
# ------------------------------------------------------------

Mrare4000 <- prepare_rarefied_matrix(
  RARE4000_FILE,
  dat$Run_ID
)

rare_dat <- add_personal_time_and_bray(
  dat %>% select(-Bray_from_personal_baseline),
  Mrare4000
)$meta

rare_follow <- rare_dat %>%
  filter(days_since_personal_baseline > 0)

rare_bray <- fit_interaction(
  rare_follow,
  outcome = "Bray_from_personal_baseline",
  timevar = "days_since_personal_baseline",
  analysis_label = "RAREFIED4000_BRAY"
)

write_csv(
  rare_bray,
  file.path(OUT, "14_SENSITIVITY_RAREFIED4000_BRAY.csv")
)

# ------------------------------------------------------------
# 13. >=2000 depth sensitivity
# ------------------------------------------------------------

meta2000 <- read_csv(META2000_FILE, show_col_types = FALSE)
alpha2000 <- read_csv(ALPHA2000_FILE, show_col_types = FALSE)

meta2000 <- meta2000 %>%
  select(-any_of(c("Observed_OTU97","Shannon_OTU97","Simpson_OTU97","rarefied_depth"))) %>%
  left_join(
    alpha2000 %>%
      select(
        Run_ID,
        rarefied_depth,
        Observed_OTU97,
        Shannon_OTU97,
        Simpson_OTU97
      ),
    by = "Run_ID"
  )

Mrel2000 <- prepare_rel_matrix(
  REL_FILE,
  meta2000$Run_ID
)

sens2 <- add_personal_time_and_bray(meta2000, Mrel2000)
dat2 <- sens2$meta

alpha2 <- bind_rows(
  lapply(alpha_metrics, function(mm) {
    fit_interaction(
      dat2,
      outcome = mm,
      timevar = "days_since_personal_baseline",
      analysis_label = "SENSITIVITY_RAREFIED2000"
    )
  })
) %>%
  mutate(FDR = p.adjust(p_LRT, method = "BH"))

bray2 <- fit_interaction(
  dat2 %>% filter(days_since_personal_baseline > 0),
  outcome = "Bray_from_personal_baseline",
  timevar = "days_since_personal_baseline",
  analysis_label = "SENSITIVITY_MIN2000"
)

write_csv(
  alpha2,
  file.path(OUT, "15_SENSITIVITY_MIN2000_ALPHA.csv")
)

write_csv(
  bray2,
  file.path(OUT, "16_SENSITIVITY_MIN2000_BRAY.csv")
)

# ------------------------------------------------------------
# 14. Trajectory summary
# ------------------------------------------------------------

trajectory <- dat %>%
  mutate(
    time_bin = case_when(
      days_since_personal_baseline == 0 ~ "0",
      days_since_personal_baseline == 1 ~ "1",
      days_since_personal_baseline == 2 ~ "2",
      days_since_personal_baseline == 3 ~ "3",
      days_since_personal_baseline >= 4 ~ "4+"
    )
  ) %>%
  group_by(pulmonary_binary, time_bin) %>%
  summarise(
    n_samples = n(),
    n_patients = n_distinct(patient_id),
    median_Bray = median(Bray_from_personal_baseline),
    q1_Bray = unname(quantile(Bray_from_personal_baseline, 0.25)),
    q3_Bray = unname(quantile(Bray_from_personal_baseline, 0.75)),
    median_Observed_OTU97 = median(Observed_OTU97),
    median_Shannon_OTU97 = median(Shannon_OTU97),
    median_Simpson_OTU97 = median(Simpson_OTU97),
    .groups = "drop"
  )

write_csv(
  trajectory,
  file.path(OUT, "17_PRIMARY_TRAJECTORY_SUMMARY.csv")
)

# ------------------------------------------------------------
# 15. Evidence classification
# ------------------------------------------------------------

p_bray <- bray_primary$p_LRT[1]
alpha_sig <- sum(alpha_models$FDR < 0.05, na.rm = TRUE)

sens_ps <- c(
  known_model$p_LRT[1],
  rare_bray$p_LRT[1],
  bray2$p_LRT[1]
)

sens_direction <- c(
  known_model$beta_interaction[1],
  rare_bray$beta_interaction[1],
  bray2$beta_interaction[1]
)

primary_dir <- sign(bray_primary$beta_interaction[1])

direction_consistent <- all(
  sign(sens_direction[!is.na(sens_direction)]) == primary_dir
)

all_key_sens_sig <- all(
  !is.na(sens_ps) & sens_ps < 0.05
)

tier <- if (
  !is.na(p_bray) &&
  p_bray < 0.05 &&
  all_key_sens_sig &&
  direction_consistent
) {
  "ROBUST_SOURCE_ASSOCIATED_TRAJECTORY_DIFFERENCE"
} else if (
  !is.na(p_bray) &&
  p_bray < 0.05 &&
  direction_consistent
) {
  "PRIMARY_SOURCE_ASSOCIATED_TRAJECTORY_DIFFERENCE_WITH_PARTIAL_SENSITIVITY_SUPPORT"
} else if (
  !is.na(p_bray) &&
  p_bray >= 0.05 &&
  alpha_sig == 0
) {
  "NO_CLEAR_SOURCE_MODIFICATION_OF_LONGITUDINAL_ECOLOGICAL_DISPLACEMENT"
} else {
  "MIXED_OR_METRIC_SPECIFIC_SOURCE_ASSOCIATION"
}

evidence <- tibble(
  primary_Bray_interaction_beta = bray_primary$beta_interaction[1],
  primary_Bray_interaction_p = p_bray,
  primary_alpha_FDR_significant = alpha_sig,
  exclude_OTHER_UNKNOWN_p = known_model$p_LRT[1],
  rarefied4000_Bray_p = rare_bray$p_LRT[1],
  min2000_Bray_p = bray2$p_LRT[1],
  key_sensitivity_direction_consistent = direction_consistent,
  evidence_tier = tier
)

write_csv(
  evidence,
  file.path(OUT, "18_SOURCE_TRAJECTORY_EVIDENCE_TIER.csv")
)

# ------------------------------------------------------------
# 16. Simple figure PDF
# ------------------------------------------------------------

pdf(
  file.path(OUT, "Figure_STEP94B2B_source_longitudinal_ecology.pdf"),
  width = 10,
  height = 8
)

par(mfrow = c(2,2), mar = c(4.5,4.5,3,1))

plot(
  dat$days_since_personal_baseline[
    dat$pulmonary_binary == "NONPULMONARY_RECORDED"
  ],
  dat$Bray_from_personal_baseline[
    dat$pulmonary_binary == "NONPULMONARY_RECORDED"
  ],
  xlab = "Days since personal baseline",
  ylab = "Bray displacement",
  main = "A. Recorded non-pulmonary",
  pch = 16
)

plot(
  dat$days_since_personal_baseline[
    dat$pulmonary_binary == "PULMONARY"
  ],
  dat$Bray_from_personal_baseline[
    dat$pulmonary_binary == "PULMONARY"
  ],
  xlab = "Days since personal baseline",
  ylab = "Bray displacement",
  main = "B. Pulmonary",
  pch = 16
)

boxplot(
  Bray_from_personal_baseline ~ pulmonary_binary,
  data = latest,
  ylab = "Baseline-to-latest Bray",
  xlab = "",
  main = "C. Latest displacement"
)

boxplot(
  Shannon_OTU97 ~ pulmonary_binary,
  data = base,
  ylab = "Rarefied Shannon",
  xlab = "",
  main = "D. Baseline OTU diversity"
)

dev.off()

# ------------------------------------------------------------
# 17. Interpretation text
# ------------------------------------------------------------

direction_text <- ifelse(
  is.na(bray_primary$beta_interaction[1]),
  "not estimable",
  ifelse(
    bray_primary$beta_interaction[1] > 0,
    "pulmonary patients show a steeper increase in displacement over personal time",
    "pulmonary patients show a shallower increase in displacement over personal time"
  )
)

writeLines(
  c(
    "STEP94B2B CORRECTED FORMAL SOURCE LONGITUDINAL ANALYSIS",
    "",
    paste0("Primary Bray interaction beta: ", signif(bray_primary$beta_interaction[1], 5)),
    paste0("Primary Bray interaction p: ", signif(p_bray, 5)),
    paste0("Primary alpha interactions significant after FDR: ", alpha_sig, "/3"),
    paste0("Direction: ", direction_text),
    paste0("Evidence tier: ", tier),
    "",
    "Inference guardrail:",
    "This is an observational infection-source association, not a causal effect.",
    "A non-significant source-by-time interaction means no clear source modification was detected; it does not prove identical trajectories."
  ),
  file.path(OUT, "19_STEP94B2B_INTERPRETATION.txt")
)

writeLines(
  c(
    paste0("Completed: ", Sys.time()),
    paste0("Evidence tier: ", tier),
    "STEP94B2B COMPLETE"
  ),
  file.path(OUT, "_STEP94B2B_COMPLETE.ok")
)

cat("STEP94B2B COMPLETE\n")
