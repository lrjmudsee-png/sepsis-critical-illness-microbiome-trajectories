# ============================================================
# Sepsis V2 - Step94B2A2
# CRA002354 SOURCE MODEL SPECIFICATION FREEZE
#
# SAFE TO RUN WHILE VSEARCH STEP94B1B IS STILL RUNNING.
#
# Purpose:
# Correct/freeze the model specification before microbiome abundance
# is available.
#
# Key decisions:
# 1) Primary within-patient time axis:
#       days_since_personal_baseline
#    not absolute ICU day.
#
# 2) Absolute ICU day remains a secondary time-axis sensitivity.
#
# 3) Primary source model is minimally adjusted.
#
# 4) Do NOT adjust for 28-day outcome.
#
# 5) Do NOT put SOFA + APACHE II + lactate simultaneously into the
#    same adjusted sensitivity model by default.
#
# 6) Use separate severity-adjusted sensitivities:
#       age + sex + SOFA
#       age + sex + APACHE II
#
# 7) Exclude OTHER_UNKNOWN in a source-definition sensitivity analysis.
# ============================================================

options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(stringr)
  library(tibble)
})

ROOT <- "E:/sepsis_project"

B2A <- file.path(
  ROOT,
  "results",
  "V2_34B2A_CRA002354_SOURCE_MODEL_DESIGN_PREP"
)

OUT <- file.path(
  ROOT,
  "results",
  "V2_34B2A2_CRA002354_SOURCE_MODEL_SPECIFICATION_FREEZE"
)

dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

META_FILE <- file.path(
  B2A,
  "01_CRA002354_SOURCE_ANALYSIS_METADATA.csv"
)

PAT_FILE <- file.path(
  B2A,
  "02_PATIENT_LEVEL_SOURCE_DESIGN.csv"
)

if (!file.exists(META_FILE) || !file.exists(PAT_FILE)) {
  stop("Step94B2A outputs are missing.")
}

meta <- read_csv(
  META_FILE,
  show_col_types = FALSE
)

pat <- read_csv(
  PAT_FILE,
  show_col_types = FALSE
)

# ------------------------------------------------------------
# 1. Build personal-baseline time axis
# ------------------------------------------------------------

meta_time <- meta %>%
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

write_csv(
  meta_time,
  file.path(
    OUT,
    "01_METADATA_WITH_PERSONAL_BASELINE_TIME_AXIS.csv"
  )
)

# ------------------------------------------------------------
# 2. First-sample timing audit
# ------------------------------------------------------------

first_day_audit <- pat %>%
  count(
    pulmonary_binary,
    first_day,
    name = "patients"
  ) %>%
  group_by(pulmonary_binary) %>%
  mutate(
    fraction =
      patients / sum(patients)
  ) %>%
  ungroup()

write_csv(
  first_day_audit,
  file.path(
    OUT,
    "02_FIRST_SAMPLE_DAY_DISTRIBUTION_BY_SOURCE.csv"
  )
)

baseline_window <- pat %>%
  group_by(pulmonary_binary) %>%
  summarise(
    patients = n(),
    first_day_le3 =
      sum(first_day <= 3),
    fraction_first_day_le3 =
      mean(first_day <= 3),
    first_day_le5 =
      sum(first_day <= 5),
    fraction_first_day_le5 =
      mean(first_day <= 5),
    median_first_day =
      median(first_day),
    max_first_day =
      max(first_day),
    .groups = "drop"
  )

write_csv(
  baseline_window,
  file.path(
    OUT,
    "03_BASELINE_WINDOW_FEASIBILITY.csv"
  )
)

# ------------------------------------------------------------
# 3. Baseline group-balance audit
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

  w <- suppressWarnings(
    wilcox.test(
      a,
      b,
      exact = FALSE
    )
  )

  # Mann-Whitney U from Wilcoxon rank-sum W
  n1 <- length(a)
  n2 <- length(b)
  W <- unname(w$statistic)
  U <- W - n1 * (n1 + 1) / 2

  rb <- 2 * U / (n1 * n2) - 1

  tibble(
    n_pulmonary = n1,
    n_nonpulmonary = n2,
    median_pulmonary = median(a),
    median_nonpulmonary = median(b),
    median_difference =
      median(a) - median(b),
    rank_biserial = rb,
    p_value = w$p.value
  )
}

num_vars <- c(
  "baseline_sofa",
  "baseline_apache_ii",
  "baseline_lactate",
  "age",
  "first_day",
  "followup_span_days"
)

balance_numeric <- bind_rows(
  lapply(num_vars, function(v) {

    out <- rank_biserial(
      pat[[v]],
      pat$pulmonary_binary
    )

    out$variable <- v
    out
  })
) %>%
  select(
    variable,
    everything()
  )

write_csv(
  balance_numeric,
  file.path(
    OUT,
    "04_BASELINE_NUMERIC_GROUP_BALANCE.csv"
  )
)

categorical_test <- function(v) {

  dd <- pat %>%
    filter(
      !is.na(.data[[v]]),
      !is.na(pulmonary_binary)
    )

  tab <- table(
    dd$pulmonary_binary,
    dd[[v]]
  )

  p <- if (
    nrow(tab) == 2 &&
    ncol(tab) == 2
  ) {
    fisher.test(tab)$p.value
  } else {
    suppressWarnings(
      chisq.test(tab)$p.value
    )
  }

  tibble(
    variable = v,
    pulmonary_distribution =
      paste(
        names(table(
          dd[[v]][
            dd$pulmonary_binary == "PULMONARY"
          ]
        )),
        as.integer(
          table(
            dd[[v]][
              dd$pulmonary_binary == "PULMONARY"
            ]
          )
        ),
        sep = "=",
        collapse = ";"
      ),
    nonpulmonary_distribution =
      paste(
        names(table(
          dd[[v]][
            dd$pulmonary_binary ==
              "NONPULMONARY_RECORDED"
          ]
        )),
        as.integer(
          table(
            dd[[v]][
              dd$pulmonary_binary ==
                "NONPULMONARY_RECORDED"
            ]
          )
        ),
        sep = "=",
        collapse = ";"
      ),
    p_value = p
  )
}

balance_cat <- bind_rows(
  lapply(
    c("sex", "outcome_28d"),
    categorical_test
  )
)

write_csv(
  balance_cat,
  file.path(
    OUT,
    "05_BASELINE_CATEGORICAL_GROUP_BALANCE.csv"
  )
)

# ------------------------------------------------------------
# 4. Severity-covariate redundancy
# ------------------------------------------------------------

sev <- pat %>%
  select(
    baseline_sofa,
    baseline_apache_ii,
    baseline_lactate,
    age
  )

cor_mat <- suppressWarnings(
  cor(
    sev,
    use = "pairwise.complete.obs",
    method = "spearman"
  )
)

cor_df <- as.data.frame(
  as.table(cor_mat),
  stringsAsFactors = FALSE
) %>%
  rename(
    variable_1 = Var1,
    variable_2 = Var2,
    spearman_rho = Freq
  )

write_csv(
  cor_df,
  file.path(
    OUT,
    "06_BASELINE_COVARIATE_SPEARMAN_CORRELATIONS.csv"
  )
)

# ------------------------------------------------------------
# 5. Freeze formal estimands and formulas
# ------------------------------------------------------------

estimands <- tibble(
  analysis_id = c(
    "A_BASELINE_ALPHA",
    "B_BASELINE_BETA",
    "C_ALPHA_LONGITUDINAL",
    "D_PERSONAL_BASELINE_BRAY",
    "E_EARLY_LATE_PAIRED",
    "F_ABSOLUTE_ICU_DAY_SENSITIVITY"
  ),
  estimand = c(
    "Difference in ecological diversity at each patient's first available sample",
    "Difference in baseline community composition at each patient's first available sample",
    "Difference in within-patient alpha-diversity trajectory after personal baseline",
    "Difference in rate/magnitude of ecological displacement from personal baseline",
    "Difference in within-patient early-to-late displacement between source groups",
    "Difference in ecological trajectory as a function of absolute ICU day"
  ),
  primary_time_axis = c(
    "personal_baseline_sample",
    "personal_baseline_sample",
    "days_since_personal_baseline",
    "days_since_personal_baseline",
    "patient-specific early_to_late",
    "absolute time_day"
  ),
  primary_model = c(
    "Wilcoxon source contrast; report effect size",
    "PERMANOVA source contrast + betadisper",
    "alpha ~ days_since_personal_baseline * pulmonary_binary + (1|patient_id)",
    "Bray_from_personal_baseline ~ days_since_personal_baseline * pulmonary_binary + (1|patient_id)",
    "paired within-patient change; compare change distributions by source group",
    "metric ~ time_day * pulmonary_binary + (1|patient_id)"
  ),
  tier = c(
    "PRIMARY",
    "PRIMARY",
    "PRIMARY",
    "PRIMARY",
    "SUPPORTIVE_INTUITIVE_CONTRAST",
    "TIME_AXIS_SENSITIVITY"
  )
)

write_csv(
  estimands,
  file.path(
    OUT,
    "07_FROZEN_SOURCE_ANALYSIS_ESTIMANDS.csv"
  )
)

# ------------------------------------------------------------
# 6. Adjustment strategy
# ------------------------------------------------------------

adjustment <- tibble(
  model_tier = c(
    "PRIMARY",
    "SENSITIVITY_SEVERITY_SOFA",
    "SENSITIVITY_SEVERITY_APACHE",
    "SENSITIVITY_LACTATE",
    "SENSITIVITY_SOURCE_DEFINITION",
    "SENSITIVITY_BASELINE_WINDOW",
    "PROHIBITED_ADJUSTMENT"
  ),
  covariates_or_rule = c(
    "No clinical covariates; source-by-time association with patient random intercept",
    "age + sex + baseline_sofa",
    "age + sex + baseline_apache_ii",
    "age + sex + baseline_lactate",
    "Exclude OTHER_UNKNOWN; compare known pulmonary vs known non-pulmonary sources",
    "For baseline cross-sectional analyses, repeat after restricting first available sample to ICU day <=3",
    "Do not adjust for outcome_28d"
  ),
  rationale = c(
    "Primary estimand is descriptive/associational and avoids overfitting 64 patients.",
    "SOFA is a prespecified baseline severity sensitivity covariate.",
    "APACHE II is an alternative severity adjustment; do not combine with SOFA by default.",
    "Lactate is a separate physiologic-severity sensitivity, not added on top of all severity scores.",
    "Tests whether OTHER_UNKNOWN drives the binary source contrast.",
    "Protects baseline comparison against late first sampling.",
    "28-day outcome occurs downstream and is not a baseline confounder."
  )
)

write_csv(
  adjustment,
  file.path(
    OUT,
    "08_FROZEN_ADJUSTMENT_STRATEGY.csv"
  )
)

# ------------------------------------------------------------
# 7. Multiple-testing families
# ------------------------------------------------------------

multiplicity <- tibble(
  family = c(
    "Primary alpha longitudinal",
    "Primary beta/displacement longitudinal",
    "Baseline ecological state",
    "Sensitivity analyses",
    "Exploratory taxa"
  ),
  tests = c(
    "Observed genera; Shannon; Simpson source-by-time interactions",
    "Bray displacement source-by-time interaction",
    "Three alpha metrics + baseline PERMANOVA",
    "Report separately; do not mix with primary FDR family",
    "Genus-level tests use BH-FDR and remain exploratory"
  ),
  correction = c(
    "BH-FDR across 3 alpha interaction tests",
    "Single prespecified primary Bray interaction; nominal p plus CI/effect size",
    "BH-FDR across alpha tests; PERMANOVA reported separately",
    "No claim upgrade from sensitivity-only significance",
    "BH-FDR"
  )
)

write_csv(
  multiplicity,
  file.path(
    OUT,
    "09_FROZEN_MULTIPLICITY_PLAN.csv"
  )
)

# ------------------------------------------------------------
# 8. Guardrails
# ------------------------------------------------------------

guardrails <- tibble(
  issue = c(
    "Personal baseline is not always ICU Day1",
    "Non-pulmonary group contains OTHER_UNKNOWN",
    "Severity differs somewhat between source groups",
    "SOFA and APACHE II may be redundant",
    "28-day mortality is downstream",
    "Detailed source categories are sparse"
  ),
  decision = c(
    "Use days_since_personal_baseline for primary displacement trajectory",
    "Keep in primary recorded-nonpulmonary group; exclude in sensitivity",
    "Use separate severity-adjusted sensitivity models",
    "Do not automatically include both together",
    "Do not use as adjustment covariate",
    "No formal multicategory inferential model"
  )
)

write_csv(
  guardrails,
  file.path(
    OUT,
    "10_FROZEN_SOURCE_MODEL_GUARDRAILS.csv"
  )
)

# ------------------------------------------------------------
# 9. Final model-specification text
# ------------------------------------------------------------

txt <- c(
  "STEP94B2A2 SOURCE MODEL SPECIFICATION FREEZE",
  "",
  "PRIMARY TIME AXIS",
  "Use days_since_personal_baseline for within-patient displacement and alpha-diversity trajectories.",
  "Each patient's first available microbiome sample defines Bray displacement = 0.",
  "",
  "SECONDARY TIME AXIS",
  "Use absolute ICU time_day as a sensitivity analysis because first sampling occurs on different ICU days across patients.",
  "",
  "PRIMARY SOURCE CONTRAST",
  "PULMONARY vs NONPULMONARY_RECORDED.",
  "",
  "SOURCE-DEFINITION SENSITIVITY",
  "Repeat after excluding OTHER_UNKNOWN from NONPULMONARY_RECORDED.",
  "",
  "PRIMARY MODELS",
  "Alpha: alpha ~ days_since_personal_baseline * pulmonary_binary + (1|patient_id)",
  "Displacement: Bray_from_personal_baseline ~ days_since_personal_baseline * pulmonary_binary + (1|patient_id)",
  "",
  "ADJUSTED SENSITIVITIES",
  "A: primary model + age + sex + baseline_sofa",
  "B: primary model + age + sex + baseline_apache_ii",
  "C: primary model + age + sex + baseline_lactate",
  "",
  "DO NOT",
  "- adjust for outcome_28d",
  "- automatically include SOFA + APACHE II + lactate simultaneously",
  "- interpret source effects causally",
  "- run formal detailed-source multicategory models",
  "",
  "BASELINE CROSS-SECTIONAL SENSITIVITY",
  "Repeat baseline ecological-state comparisons restricted to patients whose first sample was obtained by ICU Day3."
)

writeLines(
  txt,
  file.path(
    OUT,
    "11_FROZEN_MODEL_SPECIFICATION.txt"
  )
)

summary <- tibble(
  patients = nrow(pat),
  pulmonary_patients =
    sum(pat$pulmonary_binary == "PULMONARY"),
  nonpulmonary_patients =
    sum(
      pat$pulmonary_binary ==
        "NONPULMONARY_RECORDED"
    ),
  pulmonary_GE2 =
    sum(
      pat$pulmonary_binary == "PULMONARY" &
        pat$n_samples >= 2
    ),
  nonpulmonary_GE2 =
    sum(
      pat$pulmonary_binary ==
        "NONPULMONARY_RECORDED" &
        pat$n_samples >= 2
    ),
  primary_time_axis =
    "days_since_personal_baseline",
  secondary_time_axis =
    "absolute_time_day",
  outcome_28d_adjustment = FALSE,
  formal_multicategory_source_model = FALSE,
  specification_frozen = TRUE
)

write_csv(
  summary,
  file.path(
    OUT,
    "12_STEP94B2A2_SPECIFICATION_SUMMARY.csv"
  )
)

writeLines(
  c(
    paste0("Completed: ", Sys.time()),
    "STEP94B2A2 COMPLETE"
  ),
  file.path(
    OUT,
    "_STEP94B2A2_COMPLETE.ok"
  )
)

cat("STEP94B2A2 COMPLETE\n")
