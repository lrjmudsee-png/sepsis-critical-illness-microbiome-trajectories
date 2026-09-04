# ============================================================
# Sepsis V2 - Step94B2C
# CRA002354 FINAL SOURCE x LONGITUDINAL ECOLOGY + EVIDENCE FREEZE
#
# Supersedes incomplete Step94B2B.
#
# Fix:
# B2B stopped at severity-adjusted sensitivity because B1D metadata
# did not retain age/sex. B2C restores age/sex from the frozen B2A2
# metadata (with frozen-master fallback), then reruns the full analysis.
#
# Primary ecology:
# - OTU97 alpha rarefied to 4000 reads
# - corrected true-nonchimeric OTU97 relative abundance for Bray
# - personal-baseline time axis
#
# Primary inferential target:
# days_since_personal_baseline x pulmonary_binary
#
# Important:
# genus taxonomy is NOT required here.
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

B2A2_DIR <- file.path(
  ROOT, "results",
  "V2_34B2A2_CRA002354_SOURCE_MODEL_SPECIFICATION_FREEZE"
)

B2A3_DIR <- file.path(
  ROOT, "results",
  "V2_34B2A3_SOURCE_MODEL_SPEC_QC_FIX"
)

MASTER <- file.path(
  ROOT,
  "data",
  "_V2_ANALYSIS_READY",
  "00_FREEZE",
  "V2_FINAL_FROZEN_master_metadata_CLINICAL_REPAIRED_20260814_131620.csv"
)

OUT <- file.path(
  ROOT, "results",
  "V2_34B2C_CRA002354_FINAL_SOURCE_LONGITUDINAL_ANALYSIS_AND_FREEZE"
)

dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------------------
# 0. Hard readiness gates
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
  stop("Step94B1D is not ready for formal source analysis.")
}

if (!file.exists(file.path(B2A3_DIR, "_STEP94B2A3_COMPLETE.ok"))) {
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
  REL_FILE,
  ALPHA4000_FILE,
  META4000_FILE,
  RARE4000_FILE,
  ALPHA2000_FILE,
  META2000_FILE
)

missing <- needed[!file.exists(needed)]

if (length(missing) > 0) {
  stop(
    paste0(
      "Missing corrected B1D data: ",
      paste(missing, collapse = "; ")
    )
  )
}

# ------------------------------------------------------------
# 1. Clinical covariate bridge: restore age/sex
# ------------------------------------------------------------

B2A2_META <- file.path(
  B2A2_DIR,
  "01_METADATA_WITH_PERSONAL_BASELINE_TIME_AXIS.csv"
)

if (file.exists(B2A2_META)) {

  bridge_raw <- read_csv(
    B2A2_META,
    show_col_types = FALSE
  )

  bridge <- bridge_raw %>%
    select(
      Run_ID,
      patient_id,
      any_of(c("age", "sex"))
    ) %>%
    distinct(Run_ID, .keep_all = TRUE)

  bridge_source <- "B2A2_FROZEN_METADATA"

} else {

  if (!file.exists(MASTER)) {
    stop("Neither B2A2 frozen metadata nor frozen master is available.")
  }

  mm <- read_csv(
    MASTER,
    show_col_types = FALSE,
    guess_max = 50000,
    name_repair = "unique"
  ) %>%
    filter(toupper(project) == "CRA002354")

  if (!all(c("run_id","patient_id","age","sex") %in% names(mm))) {
    stop("Frozen master fallback lacks required age/sex columns.")
  }

  bridge <- mm %>%
    transmute(
      Run_ID = as.character(run_id),
      patient_id = as.character(patient_id),
      age = suppressWarnings(as.numeric(age)),
      sex = as.character(sex)
    ) %>%
    distinct(Run_ID, .keep_all = TRUE)

  bridge_source <- "FROZEN_MASTER_FALLBACK"
}

write_csv(
  bridge,
  file.path(OUT, "00_CLINICAL_COVARIATE_BRIDGE.csv")
)

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
    return(
      tibble(
        n_pulmonary = length(a),
        n_nonpulmonary = length(b),
        median_pulmonary = NA_real_,
        median_nonpulmonary = NA_real_,
        median_difference = NA_real_,
        rank_biserial = NA_real_,
        p_value = NA_real_
      )
    )
  }

  wt <- suppressWarnings(
    wilcox.test(a, b, exact = FALSE)
  )

  U <- unname(wt$statistic)

  tibble(
    n_pulmonary = length(a),
    n_nonpulmonary = length(b),
    median_pulmonary = median(a),
    median_nonpulmonary = median(b),
    median_difference = median(a) - median(b),
    rank_biserial =
      2 * U / (length(a) * length(b)) - 1,
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

  vars <- unique(
    c(
      outcome,
      timevar,
      "pulmonary_binary",
      "patient_id",
      extra_covars
    )
  )

  missing_vars <- setdiff(vars, names(df))

  if (length(missing_vars) > 0) {
    return(
      tibble(
        analysis = analysis_label,
        outcome = outcome,
        timevar = timevar,
        nobs = NA_integer_,
        npatients = NA_integer_,
        beta_interaction = NA_real_,
        se_interaction = NA_real_,
        ci95_low = NA_real_,
        ci95_high = NA_real_,
        p_LRT = NA_real_,
        singular = NA,
        status = paste0(
          "MISSING_VARIABLES:",
          paste(missing_vars, collapse = ";")
        )
      )
    )
  }

  d <- df %>%
    select(all_of(vars)) %>%
    filter(if_all(everything(), ~ !is.na(.)))

  if (
    nrow(d) < 20 ||
    n_distinct(d$patient_id) < 10
  ) {
    return(
      tibble(
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
      )
    )
  }

  d$pulmonary_binary <- factor(
    d$pulmonary_binary,
    levels = c(
      "NONPULMONARY_RECORDED",
      "PULMONARY"
    )
  )

  for (vv in extra_covars) {
    if (is.character(d[[vv]])) {
      d[[vv]] <- factor(d[[vv]])
    }
  }

  rhs_full <- paste0(
    timevar,
    " * pulmonary_binary"
  )

  rhs_red <- paste0(
    timevar,
    " + pulmonary_binary"
  )

  if (length(extra_covars) > 0) {

    ec <- paste(
      extra_covars,
      collapse = " + "
    )

    rhs_full <- paste(
      rhs_full,
      ec,
      sep = " + "
    )

    rhs_red <- paste(
      rhs_red,
      ec,
      sep = " + "
    )
  }

  f_full <- as.formula(
    paste(
      outcome,
      "~",
      rhs_full,
      "+ (1|patient_id)"
    )
  )

  f_red <- as.formula(
    paste(
      outcome,
      "~",
      rhs_red,
      "+ (1|patient_id)"
    )
  )

  full <- tryCatch(
    lmer(
      f_full,
      data = d,
      REML = FALSE,
      control = lmerControl(
        optimizer = "bobyqa"
      )
    ),
    error = function(e) NULL
  )

  red <- tryCatch(
    lmer(
      f_red,
      data = d,
      REML = FALSE,
      control = lmerControl(
        optimizer = "bobyqa"
      )
    ),
    error = function(e) NULL
  )

  if (is.null(full) || is.null(red)) {
    return(
      tibble(
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
      )
    )
  }

  cn <- names(fixef(full))

  target1 <- paste0(
    timevar,
    ":pulmonary_binaryPULMONARY"
  )

  target2 <- paste0(
    "pulmonary_binaryPULMONARY:",
    timevar
  )

  int_name <- cn[
    cn %in% c(target1, target2)
  ]

  if (length(int_name) != 1) {
    int_name <- cn[
      grepl(
        timevar,
        cn,
        fixed = TRUE
      ) &
      grepl(
        "pulmonary_binaryPULMONARY",
        cn,
        fixed = TRUE
      )
    ]
  }

  beta <- if (length(int_name)) {
    unname(fixef(full)[int_name[1]])
  } else {
    NA_real_
  }

  se <- if (length(int_name)) {
    unname(
      sqrt(diag(vcov(full)))[int_name[1]]
    )
  } else {
    NA_real_
  }

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

fit_group_slope <- function(
  df,
  group_name
) {

  d <- df %>%
    filter(
      pulmonary_binary == group_name,
      days_since_personal_baseline > 0,
      !is.na(Bray_from_personal_baseline)
    )

  if (
    nrow(d) < 10 ||
    n_distinct(d$patient_id) < 5
  ) {
    return(
      tibble(
        pulmonary_binary = group_name,
        nobs = nrow(d),
        npatients = n_distinct(d$patient_id),
        slope_per_day = NA_real_,
        se = NA_real_,
        ci95_low = NA_real_,
        ci95_high = NA_real_,
        p_LRT = NA_real_,
        status = "INSUFFICIENT_DATA"
      )
    )
  }

  full <- tryCatch(
    lmer(
      Bray_from_personal_baseline ~
        days_since_personal_baseline +
        (1|patient_id),
      data = d,
      REML = FALSE
    ),
    error = function(e) NULL
  )

  red <- tryCatch(
    lmer(
      Bray_from_personal_baseline ~
        1 +
        (1|patient_id),
      data = d,
      REML = FALSE
    ),
    error = function(e) NULL
  )

  if (
    is.null(full) ||
    is.null(red)
  ) {
    return(
      tibble(
        pulmonary_binary = group_name,
        nobs = nrow(d),
        npatients = n_distinct(d$patient_id),
        slope_per_day = NA_real_,
        se = NA_real_,
        ci95_low = NA_real_,
        ci95_high = NA_real_,
        p_LRT = NA_real_,
        status = "MODEL_ERROR"
      )
    )
  }

  beta <- unname(
    fixef(full)[
      "days_since_personal_baseline"
    ]
  )

  se <- unname(
    sqrt(diag(vcov(full)))[
      "days_since_personal_baseline"
    ]
  )

  p <- anova(red, full)$`Pr(>Chisq)`[2]

  tibble(
    pulmonary_binary = group_name,
    nobs = nrow(d),
    npatients = n_distinct(d$patient_id),
    slope_per_day = beta,
    se = se,
    ci95_low = beta - 1.96 * se,
    ci95_high = beta + 1.96 * se,
    p_LRT = p,
    status = "OK"
  )
}

prepare_rel_matrix <- function(
  file,
  keep_ids
) {

  x <- read_csv(
    file,
    show_col_types = FALSE,
    name_repair = "minimal"
  )

  feats <- setdiff(
    names(x),
    "Run_ID"
  )

  M <- as.matrix(
    x[, feats, drop = FALSE]
  )

  storage.mode(M) <- "numeric"
  rownames(M) <- x$Run_ID

  M <- M[
    intersect(
      keep_ids,
      rownames(M)
    ),
    ,
    drop = FALSE
  ]

  rs <- rowSums(M)

  if (any(rs <= 0)) {
    stop("Zero-sum relative-abundance sample found.")
  }

  M / rs
}

prepare_rarefied_matrix <- function(
  file,
  keep_ids
) {

  x <- read_csv(
    file,
    show_col_types = FALSE,
    name_repair = "minimal"
  )

  feats <- setdiff(
    names(x),
    "Run_ID"
  )

  M <- as.matrix(
    x[, feats, drop = FALSE]
  )

  storage.mode(M) <- "numeric"
  rownames(M) <- x$Run_ID

  M <- M[
    intersect(
      keep_ids,
      rownames(M)
    ),
    ,
    drop = FALSE
  ]

  M / rowSums(M)
}

add_personal_time_and_bray <- function(
  meta,
  M
) {

  meta <- meta %>%
    filter(
      Run_ID %in% rownames(M)
    ) %>%
    arrange(
      match(
        Run_ID,
        rownames(M)
      )
    ) %>%
    group_by(patient_id) %>%
    mutate(
      personal_baseline_day =
        min(time_day, na.rm = TRUE),
      days_since_personal_baseline =
        time_day -
        personal_baseline_day,
      is_personal_baseline =
        time_day ==
        personal_baseline_day
    ) %>%
    ungroup()

  M <- M[
    meta$Run_ID,
    ,
    drop = FALSE
  ]

  meta$Bray_from_personal_baseline <-
    NA_real_

  for (pid in unique(meta$patient_id)) {

    ii <- which(
      meta$patient_id == pid
    )

    base_i <- ii[
      which.min(
        meta$time_day[ii]
      )
    ]

    x0 <- M[
      meta$Run_ID[base_i],
      ,
      drop = TRUE
    ]

    for (jj in ii) {

      meta$Bray_from_personal_baseline[jj] <-
        bray_two(
          x0,
          M[
            meta$Run_ID[jj],
            ,
            drop = TRUE
          ]
        )
    }
  }

  list(
    meta = meta,
    M = M
  )
}

# ------------------------------------------------------------
# 2. Primary corrected dataset + clinical bridge
# ------------------------------------------------------------

meta4000 <- read_csv(
  META4000_FILE,
  show_col_types = FALSE
)

alpha4000 <- read_csv(
  ALPHA4000_FILE,
  show_col_types = FALSE
)

meta4000 <- meta4000 %>%
  select(
    -any_of(
      c(
        "Observed_OTU97",
        "Shannon_OTU97",
        "Simpson_OTU97",
        "rarefied_depth",
        "age",
        "sex"
      )
    )
  ) %>%
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
  ) %>%
  left_join(
    bridge %>%
      select(
        Run_ID,
        age,
        sex
      ),
    by = "Run_ID"
  )

cov_bridge_qc <- tibble(
  bridge_source = bridge_source,
  primary_samples = nrow(meta4000),
  age_nonmissing = sum(!is.na(meta4000$age)),
  sex_nonmissing = sum(!is.na(meta4000$sex)),
  age_coverage = mean(!is.na(meta4000$age)),
  sex_coverage = mean(!is.na(meta4000$sex))
)

write_csv(
  cov_bridge_qc,
  file.path(
    OUT,
    "01_CLINICAL_COVARIATE_BRIDGE_QC.csv"
  )
)

if (
  cov_bridge_qc$age_coverage < 0.95 ||
  cov_bridge_qc$sex_coverage < 0.95
) {
  stop("Age/sex clinical bridge coverage <95%.")
}

Mrel4000 <- prepare_rel_matrix(
  REL_FILE,
  meta4000$Run_ID
)

prim <- add_personal_time_and_bray(
  meta4000,
  Mrel4000
)

dat <- prim$meta
Mrel4000 <- prim$M

write_csv(
  dat,
  file.path(
    OUT,
    "02_PRIMARY_FINAL_ANALYSIS_DATA.csv"
  )
)

# ------------------------------------------------------------
# 3. Sample/patient QC
# ------------------------------------------------------------

patient_counts <- dat %>%
  count(
    patient_id,
    pulmonary_binary,
    name = "n_samples"
  )

sample_qc <- patient_counts %>%
  group_by(pulmonary_binary) %>%
  summarise(
    patients = n(),
    patients_GE2 =
      sum(n_samples >= 2),
    patients_GE3 =
      sum(n_samples >= 3),
    .groups = "drop"
  ) %>%
  left_join(
    dat %>%
      group_by(pulmonary_binary) %>%
      summarise(
        samples = n(),
        baseline_samples =
          sum(is_personal_baseline),
        followup_samples =
          sum(
            days_since_personal_baseline > 0
          ),
        .groups = "drop"
      ),
    by = "pulmonary_binary"
  )

write_csv(
  sample_qc,
  file.path(
    OUT,
    "03_PRIMARY_SAMPLE_PATIENT_QC.csv"
  )
)

# ------------------------------------------------------------
# 4. Baseline alpha
# ------------------------------------------------------------

base <- dat %>%
  filter(is_personal_baseline)

alpha_metrics <- c(
  "Observed_OTU97",
  "Shannon_OTU97",
  "Simpson_OTU97"
)

baseline_alpha <- bind_rows(
  lapply(
    alpha_metrics,
    function(mm) {

      rr <- rank_biserial(
        base[[mm]],
        base$pulmonary_binary
      )

      rr$metric <- mm

      rr
    }
  )
) %>%
  mutate(
    FDR =
      p.adjust(
        p_value,
        method = "BH"
      )
  ) %>%
  select(
    metric,
    everything()
  )

write_csv(
  baseline_alpha,
  file.path(
    OUT,
    "04_PRIMARY_BASELINE_ALPHA_SOURCE_CONTRAST.csv"
  )
)

# ------------------------------------------------------------
# 5. Baseline Bray PERMANOVA + dispersion
# ------------------------------------------------------------

Mbase <- Mrel4000[
  base$Run_ID,
  ,
  drop = FALSE
]

d_base <- vegan::vegdist(
  Mbase,
  method = "bray"
)

set.seed(20260824)

perm <- vegan::adonis2(
  d_base ~ pulmonary_binary,
  data = base,
  permutations = 9999
)

bd <- vegan::betadisper(
  d_base,
  base$pulmonary_binary
)

set.seed(20260824)

bdp <- vegan::permutest(
  bd,
  permutations = 9999
)

baseline_beta <- tibble(
  n = nrow(base),
  pseudo_F = perm$F[1],
  R2 = perm$R2[1],
  PERMANOVA_p =
    perm$`Pr(>F)`[1],
  dispersion_F =
    bdp$tab$F[1],
  dispersion_p =
    bdp$tab$`Pr(>F)`[1],
  centroid_interpretation = ifelse(
    perm$`Pr(>F)`[1] < 0.05 &
    bdp$tab$`Pr(>F)`[1] >= 0.05,
    "CLEAN_CENTROID_DIFFERENTIATION",
    ifelse(
      perm$`Pr(>F)`[1] < 0.05 &
      bdp$tab$`Pr(>F)`[1] < 0.05,
      "PERMANOVA_SIGNIFICANT_BUT_DISPERSION_CONFOUNDED",
      "NO_CLEAR_BASELINE_BETA_DIFFERENCE"
    )
  )
)

write_csv(
  baseline_beta,
  file.path(
    OUT,
    "05_PRIMARY_BASELINE_BRAY_PERMANOVA_AND_DISPERSION.csv"
  )
)

# ------------------------------------------------------------
# 6. Primary longitudinal alpha interactions
# ------------------------------------------------------------

alpha_models <- bind_rows(
  lapply(
    alpha_metrics,
    function(mm) {

      fit_interaction(
        dat,
        outcome = mm,
        timevar =
          "days_since_personal_baseline",
        analysis_label =
          "PRIMARY_RAREFIED4000"
      )
    }
  )
) %>%
  mutate(
    FDR =
      p.adjust(
        p_LRT,
        method = "BH"
      )
  )

write_csv(
  alpha_models,
  file.path(
    OUT,
    "06_PRIMARY_ALPHA_LONGITUDINAL_INTERACTIONS.csv"
  )
)

# ------------------------------------------------------------
# 7. Primary Bray source x time interaction
# ------------------------------------------------------------

follow <- dat %>%
  filter(
    days_since_personal_baseline > 0
  )

bray_primary <- fit_interaction(
  follow,
  outcome =
    "Bray_from_personal_baseline",
  timevar =
    "days_since_personal_baseline",
  analysis_label =
    "PRIMARY_CORRECTED_RELATIVE_ABUNDANCE"
)

write_csv(
  bray_primary,
  file.path(
    OUT,
    "07_PRIMARY_BRAY_DISPLACEMENT_INTERACTION.csv"
  )
)

# ------------------------------------------------------------
# 8. Source-specific Bray slopes
# ------------------------------------------------------------

source_slopes <- bind_rows(
  fit_group_slope(
    dat,
    "PULMONARY"
  ),
  fit_group_slope(
    dat,
    "NONPULMONARY_RECORDED"
  )
)

write_csv(
  source_slopes,
  file.path(
    OUT,
    "08_PRIMARY_SOURCE_SPECIFIC_BRAY_SLOPES.csv"
  )
)

# ------------------------------------------------------------
# 9. Baseline-to-latest supportive contrast
# ------------------------------------------------------------

latest <- dat %>%
  group_by(patient_id) %>%
  filter(
    days_since_personal_baseline ==
      max(
        days_since_personal_baseline
      )
  ) %>%
  slice_tail(n = 1) %>%
  ungroup() %>%
  filter(
    days_since_personal_baseline > 0
  )

latest_rr <- rank_biserial(
  latest$Bray_from_personal_baseline,
  latest$pulmonary_binary
) %>%
  mutate(
    estimand =
      "Personal-baseline to latest observed Bray displacement"
  )

write_csv(
  latest_rr,
  file.path(
    OUT,
    "09_SUPPORTIVE_BASELINE_TO_LATEST_BRAY.csv"
  )
)

# ------------------------------------------------------------
# 10. Absolute ICU-day sensitivity
# ------------------------------------------------------------

abs_alpha <- bind_rows(
  lapply(
    alpha_metrics,
    function(mm) {

      fit_interaction(
        dat,
        outcome = mm,
        timevar = "time_day",
        analysis_label =
          "ABSOLUTE_ICU_DAY"
      )
    }
  )
) %>%
  mutate(
    FDR =
      p.adjust(
        p_LRT,
        method = "BH"
      )
  )

abs_bray <- fit_interaction(
  follow,
  outcome =
    "Bray_from_personal_baseline",
  timevar =
    "time_day",
  analysis_label =
    "ABSOLUTE_ICU_DAY"
)

write_csv(
  abs_alpha,
  file.path(
    OUT,
    "10_SENSITIVITY_ABSOLUTE_ICU_DAY_ALPHA.csv"
  )
)

write_csv(
  abs_bray,
  file.path(
    OUT,
    "11_SENSITIVITY_ABSOLUTE_ICU_DAY_BRAY.csv"
  )
)

# ------------------------------------------------------------
# 11. Severity-adjusted Bray sensitivities
# ------------------------------------------------------------

sev_specs <- list(
  SOFA =
    c(
      "age",
      "sex",
      "baseline_sofa"
    ),
  APACHE =
    c(
      "age",
      "sex",
      "baseline_apache_ii"
    ),
  LACTATE =
    c(
      "age",
      "sex",
      "baseline_lactate"
    )
)

sev <- bind_rows(
  lapply(
    names(sev_specs),
    function(nm) {

      fit_interaction(
        follow,
        outcome =
          "Bray_from_personal_baseline",
        timevar =
          "days_since_personal_baseline",
        extra_covars =
          sev_specs[[nm]],
        analysis_label =
          paste0(
            "ADJUSTED_",
            nm
          )
      )
    }
  )
)

write_csv(
  sev,
  file.path(
    OUT,
    "12_SENSITIVITY_SEVERITY_ADJUSTED_BRAY.csv"
  )
)

# ------------------------------------------------------------
# 12. Exclude OTHER_UNKNOWN
# ------------------------------------------------------------

known_follow <- follow %>%
  filter(
    infection_source_group !=
      "OTHER_UNKNOWN"
  )

known_model <- fit_interaction(
  known_follow,
  outcome =
    "Bray_from_personal_baseline",
  timevar =
    "days_since_personal_baseline",
  analysis_label =
    "EXCLUDE_OTHER_UNKNOWN"
)

write_csv(
  known_model,
  file.path(
    OUT,
    "13_SENSITIVITY_EXCLUDE_OTHER_UNKNOWN_BRAY.csv"
  )
)

# ------------------------------------------------------------
# 13. Baseline <= ICU Day3 sensitivity
# ------------------------------------------------------------

base3 <- base %>%
  filter(
    personal_baseline_day <= 3
  )

base3_alpha <- bind_rows(
  lapply(
    alpha_metrics,
    function(mm) {

      rr <- rank_biserial(
        base3[[mm]],
        base3$pulmonary_binary
      )

      rr$metric <- mm

      rr
    }
  )
) %>%
  mutate(
    FDR =
      p.adjust(
        p_value,
        method = "BH"
      )
  ) %>%
  select(
    metric,
    everything()
  )

write_csv(
  base3_alpha,
  file.path(
    OUT,
    "14_SENSITIVITY_BASELINE_DAY3_ALPHA.csv"
  )
)

if (
  nrow(base3) >= 20 &&
  n_distinct(
    base3$pulmonary_binary
  ) == 2
) {

  M3 <- Mrel4000[
    base3$Run_ID,
    ,
    drop = FALSE
  ]

  d3 <- vegan::vegdist(
    M3,
    method = "bray"
  )

  set.seed(20260824)

  p3 <- vegan::adonis2(
    d3 ~ pulmonary_binary,
    data = base3,
    permutations = 9999
  )

  bd3 <- vegan::betadisper(
    d3,
    base3$pulmonary_binary
  )

  set.seed(20260824)

  bd3p <- vegan::permutest(
    bd3,
    permutations = 9999
  )

  base3_beta <- tibble(
    n = nrow(base3),
    PERMANOVA_R2 =
      p3$R2[1],
    PERMANOVA_p =
      p3$`Pr(>F)`[1],
    dispersion_F =
      bd3p$tab$F[1],
    dispersion_p =
      bd3p$tab$`Pr(>F)`[1]
  )

} else {

  base3_beta <- tibble(
    n = nrow(base3),
    PERMANOVA_R2 =
      NA_real_,
    PERMANOVA_p =
      NA_real_,
    dispersion_F =
      NA_real_,
    dispersion_p =
      NA_real_
  )
}

write_csv(
  base3_beta,
  file.path(
    OUT,
    "15_SENSITIVITY_BASELINE_DAY3_BETA.csv"
  )
)

# ------------------------------------------------------------
# 14. Rarefied-4000 Bray sensitivity
# ------------------------------------------------------------

Mrare4000 <- prepare_rarefied_matrix(
  RARE4000_FILE,
  dat$Run_ID
)

rare_obj <- add_personal_time_and_bray(
  dat %>%
    select(
      -Bray_from_personal_baseline
    ),
  Mrare4000
)

rare_follow <- rare_obj$meta %>%
  filter(
    days_since_personal_baseline > 0
  )

rare_bray <- fit_interaction(
  rare_follow,
  outcome =
    "Bray_from_personal_baseline",
  timevar =
    "days_since_personal_baseline",
  analysis_label =
    "RAREFIED4000_BRAY"
)

write_csv(
  rare_bray,
  file.path(
    OUT,
    "16_SENSITIVITY_RAREFIED4000_BRAY.csv"
  )
)

# ------------------------------------------------------------
# 15. >=2000 / rarefied-2000 sensitivity
# ------------------------------------------------------------

meta2000 <- read_csv(
  META2000_FILE,
  show_col_types = FALSE
)

alpha2000 <- read_csv(
  ALPHA2000_FILE,
  show_col_types = FALSE
)

meta2000 <- meta2000 %>%
  select(
    -any_of(
      c(
        "Observed_OTU97",
        "Shannon_OTU97",
        "Simpson_OTU97",
        "rarefied_depth",
        "age",
        "sex"
      )
    )
  ) %>%
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
  ) %>%
  left_join(
    bridge %>%
      select(
        Run_ID,
        age,
        sex
      ),
    by = "Run_ID"
  )

Mrel2000 <- prepare_rel_matrix(
  REL_FILE,
  meta2000$Run_ID
)

sens2 <- add_personal_time_and_bray(
  meta2000,
  Mrel2000
)

dat2 <- sens2$meta

alpha2 <- bind_rows(
  lapply(
    alpha_metrics,
    function(mm) {

      fit_interaction(
        dat2,
        outcome = mm,
        timevar =
          "days_since_personal_baseline",
        analysis_label =
          "SENSITIVITY_RAREFIED2000"
      )
    }
  )
) %>%
  mutate(
    FDR =
      p.adjust(
        p_LRT,
        method = "BH"
      )
  )

bray2 <- fit_interaction(
  dat2 %>%
    filter(
      days_since_personal_baseline > 0
    ),
  outcome =
    "Bray_from_personal_baseline",
  timevar =
    "days_since_personal_baseline",
  analysis_label =
    "SENSITIVITY_MIN2000"
)

write_csv(
  alpha2,
  file.path(
    OUT,
    "17_SENSITIVITY_MIN2000_ALPHA.csv"
  )
)

write_csv(
  bray2,
  file.path(
    OUT,
    "18_SENSITIVITY_MIN2000_BRAY.csv"
  )
)

# ------------------------------------------------------------
# 16. Trajectory summaries
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
  group_by(
    pulmonary_binary,
    time_bin
  ) %>%
  summarise(
    n_samples = n(),
    n_patients =
      n_distinct(patient_id),
    median_Bray =
      median(
        Bray_from_personal_baseline
      ),
    q1_Bray =
      unname(
        quantile(
          Bray_from_personal_baseline,
          0.25
        )
      ),
    q3_Bray =
      unname(
        quantile(
          Bray_from_personal_baseline,
          0.75
        )
      ),
    median_Observed_OTU97 =
      median(
        Observed_OTU97
      ),
    median_Shannon_OTU97 =
      median(
        Shannon_OTU97
      ),
    median_Simpson_OTU97 =
      median(
        Simpson_OTU97
      ),
    .groups = "drop"
  )

write_csv(
  trajectory,
  file.path(
    OUT,
    "19_PRIMARY_TRAJECTORY_SUMMARY.csv"
  )
)

# ------------------------------------------------------------
# 17. Evidence classification
# ------------------------------------------------------------

p_bray <-
  bray_primary$p_LRT[1]

alpha_sig <-
  sum(
    alpha_models$FDR < 0.05,
    na.rm = TRUE
  )

abs_alpha_sig <-
  sum(
    abs_alpha$FDR < 0.05,
    na.rm = TRUE
  )

key_sens <- tibble(
  sensitivity = c(
    "exclude_OTHER_UNKNOWN",
    "rarefied4000_Bray",
    "min2000_Bray",
    "SOFA_adjusted",
    "APACHE_adjusted",
    "lactate_adjusted"
  ),
  beta = c(
    known_model$beta_interaction[1],
    rare_bray$beta_interaction[1],
    bray2$beta_interaction[1],
    sev$beta_interaction[
      sev$analysis ==
        "ADJUSTED_SOFA"
    ][1],
    sev$beta_interaction[
      sev$analysis ==
        "ADJUSTED_APACHE"
    ][1],
    sev$beta_interaction[
      sev$analysis ==
        "ADJUSTED_LACTATE"
    ][1]
  ),
  p = c(
    known_model$p_LRT[1],
    rare_bray$p_LRT[1],
    bray2$p_LRT[1],
    sev$p_LRT[
      sev$analysis ==
        "ADJUSTED_SOFA"
    ][1],
    sev$p_LRT[
      sev$analysis ==
        "ADJUSTED_APACHE"
    ][1],
    sev$p_LRT[
      sev$analysis ==
        "ADJUSTED_LACTATE"
    ][1]
  )
)

write_csv(
  key_sens,
  file.path(
    OUT,
    "20_KEY_BRAY_SENSITIVITY_SUMMARY.csv"
  )
)

primary_dir <-
  sign(
    bray_primary$beta_interaction[1]
  )

available_dir <-
  sign(
    key_sens$beta[
      !is.na(
        key_sens$beta
      )
    ]
  )

direction_consistent <-
  if (
    length(available_dir) == 0 ||
    is.na(primary_dir)
  ) {
    NA
  } else {
    all(
      available_dir ==
        primary_dir
    )
  }

source_slopes_positive <-
  all(
    source_slopes$slope_per_day > 0,
    na.rm = TRUE
  )

tier <- case_when(

  !is.na(p_bray) &&
    p_bray < 0.05 &&
    all(
      key_sens$p[
        !is.na(key_sens$p)
      ] < 0.05
    ) &&
    isTRUE(direction_consistent) ~
    "ROBUST_SOURCE_ASSOCIATED_TRAJECTORY_DIFFERENCE",

  !is.na(p_bray) &&
    p_bray < 0.05 ~
    "PRIMARY_SOURCE_ASSOCIATED_TRAJECTORY_DIFFERENCE_WITH_PARTIAL_SENSITIVITY_SUPPORT",

  !is.na(p_bray) &&
    p_bray >= 0.05 &&
    alpha_sig == 0 &&
    isTRUE(source_slopes_positive) ~
    "NO_CLEAR_SOURCE_MODIFICATION_WITH_SHARED_POSITIVE_DISPLACEMENT_DIRECTION",

  !is.na(p_bray) &&
    p_bray >= 0.05 &&
    alpha_sig == 0 ~
    "NO_CLEAR_SOURCE_MODIFICATION_OF_LONGITUDINAL_ECOLOGICAL_DISPLACEMENT",

  TRUE ~
    "MIXED_OR_METRIC_SPECIFIC_SOURCE_ASSOCIATION"
)

evidence <- tibble(
  primary_Bray_interaction_beta =
    bray_primary$beta_interaction[1],
  primary_Bray_interaction_CI_low =
    bray_primary$ci95_low[1],
  primary_Bray_interaction_CI_high =
    bray_primary$ci95_high[1],
  primary_Bray_interaction_p =
    p_bray,
  primary_alpha_FDR_significant_tests =
    alpha_sig,
  absolute_ICU_day_alpha_FDR_significant_tests =
    abs_alpha_sig,
  baseline_PERMANOVA_p =
    baseline_beta$PERMANOVA_p,
  baseline_dispersion_p =
    baseline_beta$dispersion_p,
  baseline_beta_interpretation =
    baseline_beta$centroid_interpretation,
  pulmonary_Bray_slope =
    source_slopes$slope_per_day[
      source_slopes$pulmonary_binary ==
        "PULMONARY"
    ][1],
  pulmonary_Bray_slope_p =
    source_slopes$p_LRT[
      source_slopes$pulmonary_binary ==
        "PULMONARY"
    ][1],
  nonpulmonary_Bray_slope =
    source_slopes$slope_per_day[
      source_slopes$pulmonary_binary ==
        "NONPULMONARY_RECORDED"
    ][1],
  nonpulmonary_Bray_slope_p =
    source_slopes$p_LRT[
      source_slopes$pulmonary_binary ==
        "NONPULMONARY_RECORDED"
    ][1],
  key_Bray_sensitivity_direction_consistent =
    direction_consistent,
  evidence_tier =
    tier
)

write_csv(
  evidence,
  file.path(
    OUT,
    "21_FINAL_SOURCE_TRAJECTORY_EVIDENCE_FREEZE.csv"
  )
)

# ------------------------------------------------------------
# 18. Manuscript wording suggestions
# ------------------------------------------------------------

if (
  grepl(
    "^NO_CLEAR_SOURCE_MODIFICATION",
    tier
  )
) {

  primary_wording <- paste0(
    "In CRA002354, pulmonary and recorded non-pulmonary sepsis did not show a clear difference in the rate of within-patient ecological displacement from each patient's first available microbiome sample (source-by-time interaction beta ",
    sprintf("%.3f", bray_primary$beta_interaction[1]),
    ", 95% CI ",
    sprintf("%.3f", bray_primary$ci95_low[1]),
    " to ",
    sprintf("%.3f", bray_primary$ci95_high[1]),
    "; likelihood-ratio p=",
    signif(p_bray, 3),
    "). None of the three rarefied OTU-level alpha-diversity trajectories showed a significant source-by-time interaction after FDR correction."
  )

} else {

  primary_wording <- paste0(
    "In CRA002354, infection source was associated with longitudinal ecological displacement (source-by-time interaction beta ",
    sprintf("%.3f", bray_primary$beta_interaction[1]),
    ", 95% CI ",
    sprintf("%.3f", bray_primary$ci95_low[1]),
    " to ",
    sprintf("%.3f", bray_primary$ci95_high[1]),
    "; likelihood-ratio p=",
    signif(p_bray, 3),
    ")."
  )
}

baseline_wording <- if (
  baseline_beta$PERMANOVA_p < 0.05 &&
  baseline_beta$dispersion_p < 0.05
) {
  paste0(
    "Baseline Bray-Curtis composition differed by source in PERMANOVA (R2=",
    sprintf("%.3f", baseline_beta$R2),
    ", p=",
    signif(baseline_beta$PERMANOVA_p, 3),
    "), but dispersion also differed (p=",
    signif(baseline_beta$dispersion_p, 3),
    "), precluding interpretation as an unconfounded centroid separation."
  )
} else {
  "Baseline beta-diversity results should be interpreted according to the accompanying PERMANOVA and dispersion tests."
}

writeLines(
  c(
    "STEP94B2C FINAL MANUSCRIPT WORDING",
    "",
    primary_wording,
    "",
    baseline_wording,
    "",
    paste0(
      "Evidence tier: ",
      tier
    ),
    "",
    "Guardrail: infection source was observational; do not use causal language."
  ),
  file.path(
    OUT,
    "22_FROZEN_MANUSCRIPT_WORDING.txt"
  )
)

# ------------------------------------------------------------
# 19. Figure PDF
# ------------------------------------------------------------

pdf(
  file.path(
    OUT,
    "Figure_STEP94B2C_source_longitudinal_ecology.pdf"
  ),
  width = 10,
  height = 8
)

par(
  mfrow = c(2,2),
  mar = c(4.5,4.5,3,1)
)

plot(
  dat$days_since_personal_baseline[
    dat$pulmonary_binary ==
      "NONPULMONARY_RECORDED"
  ],
  dat$Bray_from_personal_baseline[
    dat$pulmonary_binary ==
      "NONPULMONARY_RECORDED"
  ],
  xlab =
    "Days since personal baseline",
  ylab =
    "Bray displacement",
  main =
    "A. Recorded non-pulmonary",
  pch = 16
)

plot(
  dat$days_since_personal_baseline[
    dat$pulmonary_binary ==
      "PULMONARY"
  ],
  dat$Bray_from_personal_baseline[
    dat$pulmonary_binary ==
      "PULMONARY"
  ],
  xlab =
    "Days since personal baseline",
  ylab =
    "Bray displacement",
  main =
    "B. Pulmonary",
  pch = 16
)

boxplot(
  Bray_from_personal_baseline ~
    pulmonary_binary,
  data = latest,
  ylab =
    "Baseline-to-latest Bray",
  xlab = "",
  main =
    "C. Latest displacement"
)

boxplot(
  Shannon_OTU97 ~
    pulmonary_binary,
  data = base,
  ylab =
    "Rarefied OTU Shannon",
  xlab = "",
  main =
    "D. Baseline diversity"
)

dev.off()

# ------------------------------------------------------------
# 20. Completion
# ------------------------------------------------------------

writeLines(
  c(
    "STEP94B2C FINAL SOURCE LONGITUDINAL ANALYSIS",
    "",
    paste0(
      "Clinical bridge source: ",
      bridge_source
    ),
    paste0(
      "Primary Bray interaction beta: ",
      signif(
        bray_primary$beta_interaction[1],
        5
      )
    ),
    paste0(
      "Primary Bray interaction p: ",
      signif(
        p_bray,
        5
      )
    ),
    paste0(
      "Primary alpha FDR significant: ",
      alpha_sig,
      "/3"
    ),
    paste0(
      "Evidence tier: ",
      tier
    ),
    "",
    "STEP94B2C supersedes incomplete Step94B2B."
  ),
  file.path(
    OUT,
    "23_STEP94B2C_INTERPRETATION.txt"
  )
)

writeLines(
  c(
    paste0(
      "Completed: ",
      Sys.time()
    ),
    paste0(
      "Evidence tier: ",
      tier
    ),
    "STEP94B2C COMPLETE"
  ),
  file.path(
    OUT,
    "_STEP94B2C_COMPLETE.ok"
  )
)

cat("STEP94B2C COMPLETE\n")
