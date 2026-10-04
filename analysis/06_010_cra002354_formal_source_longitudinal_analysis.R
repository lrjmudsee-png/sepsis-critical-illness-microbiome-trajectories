# ============================================================
# Sepsis V2 - Step94B2
# CRA002354 FORMAL PULMONARY vs NON-PULMONARY
# LONGITUDINAL ECOLOGICAL ANALYSIS
#
# DO NOT RUN until Step94B1B completes successfully and
# 10_STEP94B1B_ANALYSIS_READINESS.csv reports ready_for_step94B2=TRUE.
#
# Frozen design source:
# Step94B2A2 + Step94B2A3
#
# Primary abundance threshold:
# >=4000 assigned OTU97 reads
#
# Sensitivity:
# >=2000 assigned OTU97 reads
#
# Primary time axis:
# days_since_personal_baseline
#
# Primary source contrast:
# PULMONARY vs NONPULMONARY_RECORDED
#
# Formal analyses:
# A) baseline alpha
# B) baseline Bray PERMANOVA + dispersion
# C) longitudinal alpha interaction
# D) Bray-from-personal-baseline interaction
# E) baseline-to-latest Bray by source
#
# Sensitivities:
# - absolute ICU day
# - age + sex + SOFA
# - age + sex + APACHE II
# - age + sex + lactate
# - exclude OTHER_UNKNOWN
# - baseline first sample <= ICU Day3
# - >=2000 reads
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
  ROOT,
  "data",
  "CRA002354",
  "03_vsearch97_silva1382"
)

B1B_OUT <- file.path(
  ROOT,
  "results",
  "V2_34B1B_CRA002354_FASTA_WRITE_FIX_AND_OTU97_RESUME"
)

SPEC_DIR <- file.path(
  ROOT,
  "results",
  "V2_34B2A3_SOURCE_MODEL_SPEC_QC_FIX"
)

OUT <- file.path(
  ROOT,
  "results",
  "V2_34B2_CRA002354_FORMAL_SOURCE_LONGITUDINAL_ANALYSIS"
)

dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------------------
# 0. Readiness gate
# ------------------------------------------------------------

readiness_file <- file.path(
  B1B_OUT,
  "10_STEP94B1B_ANALYSIS_READINESS.csv"
)

if (!file.exists(readiness_file)) {
  stop("Step94B1B readiness file is missing. Do not run Step94B2 yet.")
}

ready <- read_csv(
  readiness_file,
  show_col_types = FALSE
)

if (
  !"ready_for_step94B2" %in% names(ready) ||
  !isTRUE(ready$ready_for_step94B2[1])
) {
  stop("Step94B1B is not ready for Step94B2.")
}

spec_complete <- file.path(
  SPEC_DIR,
  "_STEP94B2A3_COMPLETE.ok"
)

if (!file.exists(spec_complete)) {
  stop("Step94B2A3 specification QC fix is missing.")
}

# ------------------------------------------------------------
# Helpers
# ------------------------------------------------------------

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

bray_two <- function(x, y) {
  den <- sum(x + y)
  if (den == 0) return(NA_real_)
  sum(abs(x - y)) / den
}

calc_alpha <- function(M) {
  observed <- rowSums(M > 0)

  shannon <- apply(M, 1, function(p) {
    p <- p[p > 0]
    -sum(p * log(p))
  })

  simpson <- apply(M, 1, function(p) {
    1 - sum(p^2)
  })

  tibble(
    Run_ID = rownames(M),
    Observed_Genera = observed,
    Shannon = shannon,
    Simpson = simpson
  )
}

fit_lmer_interaction <- function(
  df,
  outcome,
  timevar,
  extra_covars = character()
) {

  vars <- c(
    outcome,
    timevar,
    "pulmonary_binary",
    "patient_id",
    extra_covars
  )

  d <- df %>%
    select(all_of(vars)) %>%
    filter(if_all(everything(), ~ !is.na(.)))

  if (nrow(d) < 20 || n_distinct(d$patient_id) < 10) {
    return(
      tibble(
        outcome = outcome,
        timevar = timevar,
        nobs = nrow(d),
        npatients = n_distinct(d$patient_id),
        interaction_beta = NA_real_,
        interaction_se = NA_real_,
        interaction_t = NA_real_,
        interaction_p_LRT = NA_real_,
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

  rhs <- paste0(
    timevar,
    " * pulmonary_binary"
  )

  if (length(extra_covars) > 0) {
    rhs <- paste(
      rhs,
      paste(extra_covars, collapse = " + "),
      sep = " + "
    )
  }

  f_full <- as.formula(
    paste0(
      outcome,
      " ~ ",
      rhs,
      " + (1|patient_id)"
    )
  )

  rhs_red <- paste0(
    timevar,
    " + pulmonary_binary"
  )

  if (length(extra_covars) > 0) {
    rhs_red <- paste(
      rhs_red,
      paste(extra_covars, collapse = " + "),
      sep = " + "
    )
  }

  f_red <- as.formula(
    paste0(
      outcome,
      " ~ ",
      rhs_red,
      " + (1|patient_id)"
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
        outcome = outcome,
        timevar = timevar,
        nobs = nrow(d),
        npatients = n_distinct(d$patient_id),
        interaction_beta = NA_real_,
        interaction_se = NA_real_,
        interaction_t = NA_real_,
        interaction_p_LRT = NA_real_,
        singular = NA,
        status = "MODEL_ERROR"
      )
    )
  }

  cn <- names(fixef(full))

  int_name <- cn[
    str_detect(
      cn,
      paste0(
        "^",
        fixed(timevar),
        ":pulmonary_binaryPULMONARY$|",
        "^pulmonary_binaryPULMONARY:",
        fixed(timevar),
        "$"
      )
    )
  ]

  if (length(int_name) == 0) {
    int_name <- cn[
      str_detect(
        cn,
        "pulmonary_binaryPULMONARY"
      ) &
        str_detect(
          cn,
          fixed(timevar)
        )
    ]
  }

  beta <- if (length(int_name)) {
    fixef(full)[int_name[1]]
  } else NA_real_

  se <- if (length(int_name)) {
    sqrt(diag(vcov(full)))[int_name[1]]
  } else NA_real_

  tt <- beta / se

  lrt <- anova(red, full)

  p_lrt <- lrt$`Pr(>Chisq)`[2]

  tibble(
    outcome = outcome,
    timevar = timevar,
    nobs = nrow(d),
    npatients = n_distinct(d$patient_id),
    interaction_beta = beta,
    interaction_se = se,
    interaction_t = tt,
    interaction_p_LRT = p_lrt,
    singular = isSingular(full, tol = 1e-4),
    status = "OK"
  )
}

build_dataset <- function(
  abundance_file,
  metadata_file,
  label
) {

  abund <- read_csv(
    abundance_file,
    show_col_types = FALSE,
    name_repair = "minimal"
  )

  meta <- read_csv(
    metadata_file,
    show_col_types = FALSE
  )

  if (!"Run_ID" %in% names(abund)) {
    stop(paste(label, "abundance missing Run_ID"))
  }

  genus_cols <- setdiff(names(abund), "Run_ID")

  M <- as.matrix(
    abund[, genus_cols, drop = FALSE]
  )

  storage.mode(M) <- "numeric"
  rownames(M) <- abund$Run_ID

  # Normalize defensively.
  rs <- rowSums(M)
  M <- M / rs

  common <- intersect(
    rownames(M),
    meta$Run_ID
  )

  M <- M[common, , drop = FALSE]

  meta <- meta %>%
    filter(Run_ID %in% common) %>%
    arrange(match(Run_ID, common))

  M <- M[meta$Run_ID, , drop = FALSE]

  meta <- meta %>%
    group_by(patient_id) %>%
    mutate(
      personal_baseline_day =
        min(time_day, na.rm = TRUE),
      days_since_personal_baseline =
        time_day - personal_baseline_day,
      is_personal_baseline =
        time_day == personal_baseline_day
    ) %>%
    ungroup()

  alpha <- calc_alpha(M)

  dat <- meta %>%
    left_join(
      alpha,
      by = "Run_ID"
    )

  # Bray from each patient's first available sample.
  dat$Bray_from_personal_baseline <- NA_real_

  for (pid in unique(dat$patient_id)) {

    ii <- which(dat$patient_id == pid)

    dd <- dat[ii, , drop = FALSE]

    base_row <- ii[
      which.min(dd$time_day)
    ]

    x0 <- M[
      dat$Run_ID[base_row],
      ,
      drop = TRUE
    ]

    for (jj in ii) {
      dat$Bray_from_personal_baseline[jj] <-
        bray_two(
          x0,
          M[
            dat$Run_ID[jj],
            ,
            drop = TRUE
          ]
        )
    }
  }

  list(
    M = M,
    meta = dat,
    label = label
  )
}

# ------------------------------------------------------------
# 1. Load primary and sensitivity datasets
# ------------------------------------------------------------

primary <- build_dataset(
  file.path(
    DATA_DIR,
    "CRA002354_genus_relative_abundance_PRIMARY_min4000.csv"
  ),
  file.path(
    DATA_DIR,
    "CRA002354_analysis_metadata_PRIMARY_min4000.csv"
  ),
  "PRIMARY_MIN4000"
)

sens2000 <- build_dataset(
  file.path(
    DATA_DIR,
    "CRA002354_genus_relative_abundance_SENSITIVITY_min2000.csv"
  ),
  file.path(
    DATA_DIR,
    "CRA002354_analysis_metadata_SENSITIVITY_min2000.csv"
  ),
  "SENSITIVITY_MIN2000"
)

write_csv(
  primary$meta,
  file.path(
    OUT,
    "01_PRIMARY_ANALYSIS_DATA_WITH_ECOLOGICAL_METRICS.csv"
  )
)

# ------------------------------------------------------------
# 2. Baseline alpha
# ------------------------------------------------------------

base <- primary$meta %>%
  filter(is_personal_baseline)

alpha_metrics <- c(
  "Observed_Genera",
  "Shannon",
  "Simpson"
)

baseline_alpha <- bind_rows(
  lapply(alpha_metrics, function(mm) {
    rr <- rank_biserial(
      base[[mm]],
      base$pulmonary_binary
    )
    rr$metric <- mm
    rr
  })
) %>%
  mutate(
    FDR = p.adjust(
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
    "02_PRIMARY_BASELINE_ALPHA_SOURCE_CONTRAST.csv"
  )
)

# ------------------------------------------------------------
# 3. Baseline Bray PERMANOVA + dispersion
# ------------------------------------------------------------

base_ids <- base$Run_ID

Mbase <- primary$M[
  base_ids,
  ,
  drop = FALSE
]

bray_base <- vegan::vegdist(
  Mbase,
  method = "bray"
)

perm <- vegan::adonis2(
  bray_base ~ pulmonary_binary,
  data = base,
  permutations = 9999
)

perm_out <- tibble(
  n = nrow(base),
  pseudo_F = perm$F[1],
  R2 = perm$R2[1],
  p_value = perm$`Pr(>F)`[1]
)

write_csv(
  perm_out,
  file.path(
    OUT,
    "03_PRIMARY_BASELINE_BRAY_PERMANOVA.csv"
  )
)

bd <- vegan::betadisper(
  bray_base,
  group = base$pulmonary_binary
)

bdp <- vegan::permutest(
  bd,
  permutations = 9999
)

bd_out <- tibble(
  n = nrow(base),
  F = bdp$tab$F[1],
  p_value = bdp$tab$`Pr(>F)`[1]
)

write_csv(
  bd_out,
  file.path(
    OUT,
    "04_PRIMARY_BASELINE_BRAY_DISPERSION.csv"
  )
)

# ------------------------------------------------------------
# 4. Primary longitudinal alpha interactions
# ------------------------------------------------------------

alpha_models <- bind_rows(
  lapply(
    alpha_metrics,
    function(mm) {
      fit_lmer_interaction(
        primary$meta,
        outcome = mm,
        timevar =
          "days_since_personal_baseline"
      )
    }
  )
) %>%
  mutate(
    interaction_FDR =
      p.adjust(
        interaction_p_LRT,
        method = "BH"
      )
  )

write_csv(
  alpha_models,
  file.path(
    OUT,
    "05_PRIMARY_ALPHA_LONGITUDINAL_INTERACTIONS.csv"
  )
)

# ------------------------------------------------------------
# 5. Primary Bray displacement interaction
# Baseline structural zeros excluded from inference.
# ------------------------------------------------------------

bray_follow <- primary$meta %>%
  filter(
    days_since_personal_baseline > 0
  )

bray_model <- fit_lmer_interaction(
  bray_follow,
  outcome =
    "Bray_from_personal_baseline",
  timevar =
    "days_since_personal_baseline"
)

write_csv(
  bray_model,
  file.path(
    OUT,
    "06_PRIMARY_BRAY_DISPLACEMENT_INTERACTION.csv"
  )
)

# ------------------------------------------------------------
# 6. Early-to-latest intuitive contrast
# ------------------------------------------------------------

latest <- primary$meta %>%
  group_by(patient_id) %>%
  filter(
    days_since_personal_baseline ==
      max(days_since_personal_baseline)
  ) %>%
  slice_tail(n = 1) %>%
  ungroup() %>%
  filter(
    days_since_personal_baseline > 0
  )

latest_contrast <- rank_biserial(
  latest$Bray_from_personal_baseline,
  latest$pulmonary_binary
) %>%
  mutate(
    estimand =
      "Bray displacement from personal baseline to latest observed sample"
  )

write_csv(
  latest_contrast,
  file.path(
    OUT,
    "07_SUPPORTIVE_BASELINE_TO_LATEST_BRAY_BY_SOURCE.csv"
  )
)

# ------------------------------------------------------------
# 7. Absolute ICU-day sensitivity
# ------------------------------------------------------------

absolute_models <- bind_rows(
  lapply(
    c(
      alpha_metrics,
      "Bray_from_personal_baseline"
    ),
    function(mm) {

      dd <- primary$meta

      if (
        mm ==
          "Bray_from_personal_baseline"
      ) {
        dd <- dd %>%
          filter(
            days_since_personal_baseline > 0
          )
      }

      fit_lmer_interaction(
        dd,
        outcome = mm,
        timevar = "time_day"
      )
    }
  )
)

write_csv(
  absolute_models,
  file.path(
    OUT,
    "08_SENSITIVITY_ABSOLUTE_ICU_DAY_INTERACTIONS.csv"
  )
)

# ------------------------------------------------------------
# 8. Severity-adjusted Bray sensitivities
# ------------------------------------------------------------

severity_specs <- list(
  SOFA = c(
    "age",
    "sex",
    "baseline_sofa"
  ),
  APACHE = c(
    "age",
    "sex",
    "baseline_apache_ii"
  ),
  LACTATE = c(
    "age",
    "sex",
    "baseline_lactate"
  )
)

severity_bray <- bind_rows(
  lapply(
    names(severity_specs),
    function(nm) {

      rr <- fit_lmer_interaction(
        bray_follow,
        outcome =
          "Bray_from_personal_baseline",
        timevar =
          "days_since_personal_baseline",
        extra_covars =
          severity_specs[[nm]]
      )

      rr$adjustment <- nm
      rr
    }
  )
)

write_csv(
  severity_bray,
  file.path(
    OUT,
    "09_SENSITIVITY_SEVERITY_ADJUSTED_BRAY_INTERACTIONS.csv"
  )
)

# ------------------------------------------------------------
# 9. Source-definition sensitivity:
# exclude OTHER_UNKNOWN
# ------------------------------------------------------------

known_source <- primary$meta %>%
  filter(
    infection_source_group !=
      "OTHER_UNKNOWN"
  )

known_bray <- known_source %>%
  filter(
    days_since_personal_baseline > 0
  )

known_model <- fit_lmer_interaction(
  known_bray,
  outcome =
    "Bray_from_personal_baseline",
  timevar =
    "days_since_personal_baseline"
)

write_csv(
  known_model,
  file.path(
    OUT,
    "10_SENSITIVITY_EXCLUDE_OTHER_UNKNOWN_BRAY_INTERACTION.csv"
  )
)

# ------------------------------------------------------------
# 10. Baseline-window sensitivity: first sample <= ICU Day3
# ------------------------------------------------------------

base_day3 <- base %>%
  filter(
    personal_baseline_day <= 3
  )

baseline_alpha_day3 <- bind_rows(
  lapply(
    alpha_metrics,
    function(mm) {
      rr <- rank_biserial(
        base_day3[[mm]],
        base_day3$pulmonary_binary
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
  baseline_alpha_day3,
  file.path(
    OUT,
    "11_SENSITIVITY_BASELINE_DAY3_ALPHA.csv"
  )
)

if (
  nrow(base_day3) >= 20 &&
  n_distinct(
    base_day3$pulmonary_binary
  ) == 2
) {

  Mbase3 <- primary$M[
    base_day3$Run_ID,
    ,
    drop = FALSE
  ]

  d3 <- vegan::vegdist(
    Mbase3,
    method = "bray"
  )

  p3 <- vegan::adonis2(
    d3 ~ pulmonary_binary,
    data = base_day3,
    permutations = 9999
  )

  bd3 <- vegan::betadisper(
    d3,
    base_day3$pulmonary_binary
  )

  bd3p <- vegan::permutest(
    bd3,
    permutations = 9999
  )

  baseline_beta_day3 <- tibble(
    n = nrow(base_day3),
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

  baseline_beta_day3 <- tibble(
    n = nrow(base_day3),
    PERMANOVA_R2 = NA_real_,
    PERMANOVA_p = NA_real_,
    dispersion_F = NA_real_,
    dispersion_p = NA_real_
  )
}

write_csv(
  baseline_beta_day3,
  file.path(
    OUT,
    "12_SENSITIVITY_BASELINE_DAY3_BETA.csv"
  )
)

# ------------------------------------------------------------
# 11. >=2000-read depth sensitivity
# ------------------------------------------------------------

sens_bray_follow <- sens2000$meta %>%
  filter(
    days_since_personal_baseline > 0
  )

sens_bray_model <- fit_lmer_interaction(
  sens_bray_follow,
  outcome =
    "Bray_from_personal_baseline",
  timevar =
    "days_since_personal_baseline"
)

sens_alpha <- bind_rows(
  lapply(
    alpha_metrics,
    function(mm) {
      fit_lmer_interaction(
        sens2000$meta,
        outcome = mm,
        timevar =
          "days_since_personal_baseline"
      )
    }
  )
) %>%
  mutate(
    interaction_FDR =
      p.adjust(
        interaction_p_LRT,
        method = "BH"
      )
  )

write_csv(
  sens_bray_model,
  file.path(
    OUT,
    "13_SENSITIVITY_MIN2000_BRAY_INTERACTION.csv"
  )
)

write_csv(
  sens_alpha,
  file.path(
    OUT,
    "14_SENSITIVITY_MIN2000_ALPHA_INTERACTIONS.csv"
  )
)

# ------------------------------------------------------------
# 12. Trajectory summaries for figure/manuscript
# ------------------------------------------------------------

trajectory_summary <- primary$meta %>%
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
        Bray_from_personal_baseline,
        na.rm = TRUE
      ),
    q1_Bray =
      quantile(
        Bray_from_personal_baseline,
        0.25,
        na.rm = TRUE
      ),
    q3_Bray =
      quantile(
        Bray_from_personal_baseline,
        0.75,
        na.rm = TRUE
      ),
    median_Shannon =
      median(
        Shannon,
        na.rm = TRUE
      ),
    median_Observed =
      median(
        Observed_Genera,
        na.rm = TRUE
      ),
    .groups = "drop"
  )

write_csv(
  trajectory_summary,
  file.path(
    OUT,
    "15_PRIMARY_TRAJECTORY_SUMMARY_BY_SOURCE.csv"
  )
)

# ------------------------------------------------------------
# 13. Evidence tier
# ------------------------------------------------------------

bray_p <- bray_model$interaction_p_LRT[1]
alpha_sig_n <- sum(
  alpha_models$interaction_FDR < 0.05,
  na.rm = TRUE
)

severity_consistent <- all(
  severity_bray$interaction_p_LRT < 0.05,
  na.rm = TRUE
)

source_sens_sig <-
  known_model$interaction_p_LRT[1] < 0.05

depth_sens_sig <-
  sens_bray_model$interaction_p_LRT[1] < 0.05

tier <- case_when(

  !is.na(bray_p) &&
    bray_p < 0.05 &&
    source_sens_sig &&
    depth_sens_sig &&
    severity_consistent ~
    "ROBUST_SOURCE_ASSOCIATED_TRAJECTORY_DIFFERENCE",

  !is.na(bray_p) &&
    bray_p < 0.05 &&
    (source_sens_sig || depth_sens_sig) ~
    "SUPPORTED_SOURCE_ASSOCIATED_TRAJECTORY_DIFFERENCE",

  (!is.na(bray_p) && bray_p >= 0.05) &&
    alpha_sig_n == 0 ~
    "NO_CLEAR_SOURCE_MODIFICATION_OF_ECOLOGICAL_TRAJECTORY",

  TRUE ~
    "MIXED_OR_METRIC_SPECIFIC_SOURCE_ASSOCIATION"
)

evidence <- tibble(
  primary_Bray_interaction_p = bray_p,
  primary_alpha_FDR_significant_tests =
    alpha_sig_n,
  severity_adjusted_Bray_all_p_lt_0_05 =
    severity_consistent,
  exclude_OTHER_UNKNOWN_Bray_p =
    known_model$interaction_p_LRT[1],
  min2000_Bray_p =
    sens_bray_model$interaction_p_LRT[1],
  evidence_tier = tier
)

write_csv(
  evidence,
  file.path(
    OUT,
    "16_SOURCE_TRAJECTORY_EVIDENCE_TIER.csv"
  )
)

# ------------------------------------------------------------
# 14. Figure
# ------------------------------------------------------------

pdf(
  file.path(
    OUT,
    "Figure_STEP94B2_source_longitudinal_ecology.pdf"
  ),
  width = 10,
  height = 8
)

par(
  mfrow = c(2,2),
  mar = c(4.5,4.5,3,1)
)

# A. Bray by source over personal time
plot(
  primary$meta$days_since_personal_baseline[
    primary$meta$pulmonary_binary ==
      "NONPULMONARY_RECORDED"
  ],
  primary$meta$Bray_from_personal_baseline[
    primary$meta$pulmonary_binary ==
      "NONPULMONARY_RECORDED"
  ],
  xlab = "Days since personal baseline",
  ylab = "Bray displacement",
  main = "A. Recorded non-pulmonary",
  pch = 16
)

plot(
  primary$meta$days_since_personal_baseline[
    primary$meta$pulmonary_binary ==
      "PULMONARY"
  ],
  primary$meta$Bray_from_personal_baseline[
    primary$meta$pulmonary_binary ==
      "PULMONARY"
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
  Shannon ~ pulmonary_binary,
  data = base,
  ylab = "Shannon",
  xlab = "",
  main = "D. Baseline diversity"
)

dev.off()

# ------------------------------------------------------------
# 15. Interpretation
# ------------------------------------------------------------

txt <- c(
  "STEP94B2 FORMAL SOURCE LONGITUDINAL ANALYSIS",
  "",
  paste0(
    "Primary Bray source-by-time interaction p = ",
    signif(bray_p, 4)
  ),
  paste0(
    "Primary alpha interaction tests significant after FDR = ",
    alpha_sig_n,
    "/3"
  ),
  paste0(
    "Evidence tier: ",
    tier
  ),
  "",
  "Interpretation rule:",
  "If the primary source-by-time interaction is non-significant, do not conclude that pulmonary and non-pulmonary groups are identical; conclude that no clear source modification of longitudinal ecological displacement was detected.",
  "If the interaction is significant, retain causal guardrails because infection source was not randomized.",
  "",
  "Detailed anatomical source categories remain descriptive/exploratory."
)

writeLines(
  txt,
  file.path(
    OUT,
    "17_STEP94B2_INTERPRETATION.txt"
  )
)

writeLines(
  c(
    paste0("Completed: ", Sys.time()),
    paste0("Evidence tier: ", tier),
    "STEP94B2 COMPLETE"
  ),
  file.path(
    OUT,
    "_STEP94B2_COMPLETE.ok"
  )
)

cat("STEP94B2 COMPLETE\n")
