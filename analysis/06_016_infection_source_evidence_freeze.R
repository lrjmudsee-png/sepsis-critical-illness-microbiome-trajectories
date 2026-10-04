# ============================================================
# Sepsis V2 - Step94B2D
# CRA002354 INFECTION-SOURCE EVIDENCE FREEZE
#
# Purpose:
# Freeze the final, corrected infection-source branch after Step94B2C.
# This step performs no new inferential analysis.
# It copies the authoritative outputs, creates a compact evidence
# registry, and writes manuscript-ready wording/guardrails.
# ============================================================

options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tibble)
})

ROOT <- "E:/sepsis_project"

SRC <- file.path(
  ROOT,
  "results",
  "V2_34B2C_CRA002354_FINAL_SOURCE_LONGITUDINAL_ANALYSIS_AND_FREEZE"
)

OUT <- file.path(
  ROOT,
  "results",
  "V2_34B2D_INFECTION_SOURCE_EVIDENCE_FREEZE"
)

dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

required <- c(
  "_STEP94B2C_COMPLETE.ok",
  "04_PRIMARY_BASELINE_ALPHA_SOURCE_CONTRAST.csv",
  "05_PRIMARY_BASELINE_BRAY_PERMANOVA_AND_DISPERSION.csv",
  "06_PRIMARY_ALPHA_LONGITUDINAL_INTERACTIONS.csv",
  "07_PRIMARY_BRAY_DISPLACEMENT_INTERACTION.csv",
  "08_PRIMARY_SOURCE_SPECIFIC_BRAY_SLOPES.csv",
  "09_SUPPORTIVE_BASELINE_TO_LATEST_BRAY.csv",
  "12_SENSITIVITY_SEVERITY_ADJUSTED_BRAY.csv",
  "13_SENSITIVITY_EXCLUDE_OTHER_UNKNOWN_BRAY.csv",
  "14_SENSITIVITY_BASELINE_DAY3_ALPHA.csv",
  "15_SENSITIVITY_BASELINE_DAY3_BETA.csv",
  "16_SENSITIVITY_RAREFIED4000_BRAY.csv",
  "17_SENSITIVITY_MIN2000_ALPHA.csv",
  "18_SENSITIVITY_MIN2000_BRAY.csv",
  "20_KEY_BRAY_SENSITIVITY_SUMMARY.csv",
  "21_FINAL_SOURCE_TRAJECTORY_EVIDENCE_FREEZE.csv",
  "22_FROZEN_MANUSCRIPT_WORDING.txt",
  "23_STEP94B2C_INTERPRETATION.txt"
)

missing <- required[
  !file.exists(
    file.path(SRC, required)
  )
]

if (length(missing) > 0) {
  stop(
    paste0(
      "Missing Step94B2C outputs: ",
      paste(missing, collapse = "; ")
    )
  )
}

# ------------------------------------------------------------
# 1. Read authoritative results
# ------------------------------------------------------------

alpha_base <- read_csv(
  file.path(
    SRC,
    "04_PRIMARY_BASELINE_ALPHA_SOURCE_CONTRAST.csv"
  ),
  show_col_types = FALSE
)

beta_base <- read_csv(
  file.path(
    SRC,
    "05_PRIMARY_BASELINE_BRAY_PERMANOVA_AND_DISPERSION.csv"
  ),
  show_col_types = FALSE
)

alpha_long <- read_csv(
  file.path(
    SRC,
    "06_PRIMARY_ALPHA_LONGITUDINAL_INTERACTIONS.csv"
  ),
  show_col_types = FALSE
)

bray_primary <- read_csv(
  file.path(
    SRC,
    "07_PRIMARY_BRAY_DISPLACEMENT_INTERACTION.csv"
  ),
  show_col_types = FALSE
)

slopes <- read_csv(
  file.path(
    SRC,
    "08_PRIMARY_SOURCE_SPECIFIC_BRAY_SLOPES.csv"
  ),
  show_col_types = FALSE
)

latest <- read_csv(
  file.path(
    SRC,
    "09_SUPPORTIVE_BASELINE_TO_LATEST_BRAY.csv"
  ),
  show_col_types = FALSE
)

sev <- read_csv(
  file.path(
    SRC,
    "12_SENSITIVITY_SEVERITY_ADJUSTED_BRAY.csv"
  ),
  show_col_types = FALSE
)

known <- read_csv(
  file.path(
    SRC,
    "13_SENSITIVITY_EXCLUDE_OTHER_UNKNOWN_BRAY.csv"
  ),
  show_col_types = FALSE
)

base3_beta <- read_csv(
  file.path(
    SRC,
    "15_SENSITIVITY_BASELINE_DAY3_BETA.csv"
  ),
  show_col_types = FALSE
)

rare4 <- read_csv(
  file.path(
    SRC,
    "16_SENSITIVITY_RAREFIED4000_BRAY.csv"
  ),
  show_col_types = FALSE
)

min2 <- read_csv(
  file.path(
    SRC,
    "18_SENSITIVITY_MIN2000_BRAY.csv"
  ),
  show_col_types = FALSE
)

freeze <- read_csv(
  file.path(
    SRC,
    "21_FINAL_SOURCE_TRAJECTORY_EVIDENCE_FREEZE.csv"
  ),
  show_col_types = FALSE
)

# ------------------------------------------------------------
# 2. Compact evidence registry
# ------------------------------------------------------------

registry <- bind_rows(
  tibble(
    domain = "PRIMARY_LONGITUDINAL_BRAY",
    estimand = "source_by_personal_time_interaction",
    effect = bray_primary$beta_interaction[1],
    ci_low = bray_primary$ci95_low[1],
    ci_high = bray_primary$ci95_high[1],
    p_value = bray_primary$p_LRT[1],
    interpretation = "NO_CLEAR_SOURCE_MODIFICATION"
  ),
  tibble(
    domain = "PRIMARY_LONGITUDINAL_ALPHA",
    estimand = alpha_long$outcome,
    effect = alpha_long$beta_interaction,
    ci_low = alpha_long$ci95_low,
    ci_high = alpha_long$ci95_high,
    p_value = alpha_long$FDR,
    interpretation = ifelse(
      alpha_long$FDR < 0.05,
      "SOURCE_MODIFICATION_DETECTED",
      "NO_CLEAR_SOURCE_MODIFICATION"
    )
  ),
  tibble(
    domain = "BASELINE_BETA",
    estimand = "PERMANOVA_source",
    effect = beta_base$R2[1],
    ci_low = NA_real_,
    ci_high = NA_real_,
    p_value = beta_base$PERMANOVA_p[1],
    interpretation =
      "SIGNIFICANT_BUT_DISPERSION_CONFOUNDED"
  ),
  tibble(
    domain = "BASELINE_BETA",
    estimand = "dispersion_source",
    effect = beta_base$dispersion_F[1],
    ci_low = NA_real_,
    ci_high = NA_real_,
    p_value = beta_base$dispersion_p[1],
    interpretation =
      "DISPERSION_HETEROGENEITY_PRESENT"
  ),
  tibble(
    domain = "SUPPORTIVE_LATEST_BRAY",
    estimand =
      "baseline_to_latest_source_contrast",
    effect = latest$rank_biserial[1],
    ci_low = NA_real_,
    ci_high = NA_real_,
    p_value = latest$p_value[1],
    interpretation =
      "NO_CLEAR_SOURCE_DIFFERENCE"
  )
)

write_csv(
  registry,
  file.path(
    OUT,
    "01_FINAL_INFECTION_SOURCE_EVIDENCE_REGISTRY.csv"
  )
)

# ------------------------------------------------------------
# 3. Sensitivity robustness registry
# ------------------------------------------------------------

sens <- bind_rows(
  sev %>%
    transmute(
      sensitivity = analysis,
      beta = beta_interaction,
      ci_low = ci95_low,
      ci_high = ci95_high,
      p_value = p_LRT
    ),
  known %>%
    transmute(
      sensitivity = analysis,
      beta = beta_interaction,
      ci_low = ci95_low,
      ci_high = ci95_high,
      p_value = p_LRT
    ),
  rare4 %>%
    transmute(
      sensitivity = analysis,
      beta = beta_interaction,
      ci_low = ci95_low,
      ci_high = ci95_high,
      p_value = p_LRT
    ),
  min2 %>%
    transmute(
      sensitivity = analysis,
      beta = beta_interaction,
      ci_low = ci95_low,
      ci_high = ci95_high,
      p_value = p_LRT
    )
) %>%
  mutate(
    same_direction_as_primary =
      sign(beta) ==
      sign(
        bray_primary$beta_interaction[1]
      ),
    significant_0_05 =
      p_value < 0.05
  )

write_csv(
  sens,
  file.path(
    OUT,
    "02_FINAL_BRAY_SENSITIVITY_ROBUSTNESS_REGISTRY.csv"
  )
)

# ------------------------------------------------------------
# 4. Source-specific slope registry
# ------------------------------------------------------------

slope_registry <- slopes %>%
  mutate(
    direction = case_when(
      slope_per_day > 0 ~ "POSITIVE",
      slope_per_day < 0 ~ "NEGATIVE",
      TRUE ~ "ZERO"
    ),
    significant_0_05 =
      p_LRT < 0.05
  )

write_csv(
  slope_registry,
  file.path(
    OUT,
    "03_SOURCE_SPECIFIC_SLOPE_REGISTRY.csv"
  )
)

# ------------------------------------------------------------
# 5. Baseline-window robustness
# ------------------------------------------------------------

baseline_window <- tibble(
  primary_PERMANOVA_R2 =
    beta_base$R2[1],
  primary_PERMANOVA_p =
    beta_base$PERMANOVA_p[1],
  primary_dispersion_p =
    beta_base$dispersion_p[1],
  day3_restricted_PERMANOVA_R2 =
    base3_beta$PERMANOVA_R2[1],
  day3_restricted_PERMANOVA_p =
    base3_beta$PERMANOVA_p[1],
  day3_restricted_dispersion_p =
    base3_beta$dispersion_p[1],
  conclusion =
    "BASELINE_SOURCE_DIFFERENCE_NOT_ROBUST_TO_DAY3_RESTRICTION_AND_PRIMARY_IS_DISPERSION_CONFOUNDED"
)

write_csv(
  baseline_window,
  file.path(
    OUT,
    "04_BASELINE_BETA_ROBUSTNESS_FREEZE.csv"
  )
)

# ------------------------------------------------------------
# 6. Final branch decision
# ------------------------------------------------------------

branch <- tibble(
  cohort = "CRA002354",
  branch = "INFECTION_SOURCE_LONGITUDINAL",
  primary_exposure =
    "PULMONARY_vs_NONPULMONARY_RECORDED",
  primary_time_axis =
    "DAYS_SINCE_PERSONAL_BASELINE",
  primary_Bray_interaction_beta =
    bray_primary$beta_interaction[1],
  primary_Bray_interaction_p =
    bray_primary$p_LRT[1],
  primary_alpha_FDR_significant =
    sum(
      alpha_long$FDR < 0.05,
      na.rm = TRUE
    ),
  all_key_Bray_sensitivities_same_direction =
    all(
      sens$same_direction_as_primary,
      na.rm = TRUE
    ),
  any_key_Bray_sensitivity_significant =
    any(
      sens$significant_0_05,
      na.rm = TRUE
    ),
  pulmonary_slope_positive =
    slopes$slope_per_day[
      slopes$pulmonary_binary ==
        "PULMONARY"
    ][1] > 0,
  pulmonary_slope_p =
    slopes$p_LRT[
      slopes$pulmonary_binary ==
        "PULMONARY"
    ][1],
  nonpulmonary_slope_positive =
    slopes$slope_per_day[
      slopes$pulmonary_binary ==
        "NONPULMONARY_RECORDED"
    ][1] > 0,
  nonpulmonary_slope_p =
    slopes$p_LRT[
      slopes$pulmonary_binary ==
        "NONPULMONARY_RECORDED"
    ][1],
  baseline_beta_status =
    "DISPERSION_CONFOUNDED_AND_NOT_ROBUST_TO_DAY3_RESTRICTION",
  evidence_tier =
    freeze$evidence_tier[1],
  branch_frozen = TRUE
)

write_csv(
  branch,
  file.path(
    OUT,
    "05_FINAL_INFECTION_SOURCE_BRANCH_FREEZE.csv"
  )
)

# ------------------------------------------------------------
# 7. Frozen manuscript wording
# ------------------------------------------------------------

main_wording <- paste0(
  "In CRA002354, pulmonary and recorded non-pulmonary sepsis did not show a clear difference in the rate of within-patient ecological displacement from each patient's first available microbiome sample (source-by-time interaction beta ",
  sprintf("%.3f", bray_primary$beta_interaction[1]),
  ", 95% CI ",
  sprintf("%.3f", bray_primary$ci95_low[1]),
  " to ",
  sprintf("%.3f", bray_primary$ci95_high[1]),
  "; likelihood-ratio p=",
  signif(bray_primary$p_LRT[1], 3),
  "). None of the three rarefied OTU-level alpha-diversity trajectories showed a significant source-by-time interaction after FDR correction."
)

direction_wording <- paste0(
  "Both source groups showed positive estimated Bray displacement slopes, although neither source-specific slope reached conventional statistical significance (pulmonary beta=",
  sprintf(
    "%.3f",
    slopes$slope_per_day[
      slopes$pulmonary_binary ==
        "PULMONARY"
    ][1]
  ),
  ", p=",
  signif(
    slopes$p_LRT[
      slopes$pulmonary_binary ==
        "PULMONARY"
    ][1],
    3
  ),
  "; recorded non-pulmonary beta=",
  sprintf(
    "%.3f",
    slopes$slope_per_day[
      slopes$pulmonary_binary ==
        "NONPULMONARY_RECORDED"
    ][1]
  ),
  ", p=",
  signif(
    slopes$p_LRT[
      slopes$pulmonary_binary ==
        "NONPULMONARY_RECORDED"
    ][1],
    3
  ),
  ")."
)

baseline_wording <- paste0(
  "Baseline Bray-Curtis composition differed by source in the unrestricted PERMANOVA (R2=",
  sprintf("%.3f", beta_base$R2[1]),
  ", p=",
  signif(beta_base$PERMANOVA_p[1], 3),
  "), but dispersion also differed (p=",
  signif(beta_base$dispersion_p[1], 3),
  "), and the association was not reproduced when baseline sampling was restricted to ICU Day 3 or earlier (R2=",
  sprintf("%.3f", base3_beta$PERMANOVA_R2[1]),
  ", p=",
  signif(base3_beta$PERMANOVA_p[1], 3),
  ")."
)

writeLines(
  c(
    "STEP94B2D FROZEN MANUSCRIPT WORDING",
    "",
    main_wording,
    "",
    direction_wording,
    "",
    baseline_wording,
    "",
    "Recommended synthesis:",
    "These findings do not support a strong infection-source modification of longitudinal ecological displacement. The positive slope direction in both source groups is compatible with a shared direction of ecological change, but the source-specific slopes themselves were not statistically significant.",
    "",
    "Do not write:",
    "- pulmonary and non-pulmonary trajectories are identical",
    "- infection source has no effect",
    "- both groups showed statistically significant progressive displacement",
    "- baseline community composition clearly differed by infection source",
    "- infection source causally determines microbiome trajectory"
  ),
  file.path(
    OUT,
    "06_FROZEN_MANUSCRIPT_WORDING_AND_GUARDRAILS.txt"
  )
)

# ------------------------------------------------------------
# 8. Copy authoritative source outputs
# ------------------------------------------------------------

copy_files <- c(
  "04_PRIMARY_BASELINE_ALPHA_SOURCE_CONTRAST.csv",
  "05_PRIMARY_BASELINE_BRAY_PERMANOVA_AND_DISPERSION.csv",
  "06_PRIMARY_ALPHA_LONGITUDINAL_INTERACTIONS.csv",
  "07_PRIMARY_BRAY_DISPLACEMENT_INTERACTION.csv",
  "08_PRIMARY_SOURCE_SPECIFIC_BRAY_SLOPES.csv",
  "09_SUPPORTIVE_BASELINE_TO_LATEST_BRAY.csv",
  "12_SENSITIVITY_SEVERITY_ADJUSTED_BRAY.csv",
  "13_SENSITIVITY_EXCLUDE_OTHER_UNKNOWN_BRAY.csv",
  "15_SENSITIVITY_BASELINE_DAY3_BETA.csv",
  "16_SENSITIVITY_RAREFIED4000_BRAY.csv",
  "17_SENSITIVITY_MIN2000_ALPHA.csv",
  "18_SENSITIVITY_MIN2000_BRAY.csv",
  "20_KEY_BRAY_SENSITIVITY_SUMMARY.csv",
  "21_FINAL_SOURCE_TRAJECTORY_EVIDENCE_FREEZE.csv",
  "Figure_STEP94B2C_source_longitudinal_ecology.pdf"
)

audit <- tibble(
  file = copy_files,
  copied = FALSE
)

for (i in seq_along(copy_files)) {

  src <- file.path(
    SRC,
    copy_files[i]
  )

  if (file.exists(src)) {

    file.copy(
      src,
      file.path(
        OUT,
        paste0(
          "SOURCE_",
          copy_files[i]
        )
      ),
      overwrite = TRUE
    )

    audit$copied[i] <- TRUE
  }
}

write_csv(
  audit,
  file.path(
    OUT,
    "07_FREEZE_COPY_AUDIT.csv"
  )
)

writeLines(
  c(
    "STEP94B2D INFECTION-SOURCE EVIDENCE FREEZE",
    "",
    paste0(
      "Evidence tier: ",
      freeze$evidence_tier[1]
    ),
    "",
    "Final interpretation:",
    "No clear infection-source modification of longitudinal ecological displacement was detected.",
    "Both source groups had positive estimated displacement slopes, but neither source-specific slope was individually significant.",
    "The unrestricted baseline PERMANOVA signal was dispersion-confounded and disappeared in the early-baseline sensitivity.",
    "",
    "This branch is frozen."
  ),
  file.path(
    OUT,
    "08_STEP94B2D_INTERPRETATION.txt"
  )
)

writeLines(
  c(
    paste0(
      "Completed: ",
      Sys.time()
    ),
    "Branch frozen: TRUE",
    "STEP94B2D COMPLETE"
  ),
  file.path(
    OUT,
    "_STEP94B2D_COMPLETE.ok"
  )
)

cat("STEP94B2D COMPLETE\n")
