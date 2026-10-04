# ============================================================
# Sepsis V2 - Step94B3C
# CRA002354 EXPLORATORY GENUS EVIDENCE FREEZE
#
# Purpose:
# Freeze the final genus-level infection-source branch using:
# - Step94B3A primary exploratory genus analysis
# - Step94B3B genus-coverage robustness sensitivity
#
# No new inferential analysis is performed.
# ============================================================

options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tibble)
})

ROOT <- "E:/sepsis_project"

B3A <- file.path(
  ROOT,
  "results",
  "V2_34B3A_CRA002354_EXPLORATORY_GENUS_SOURCE_ANALYSIS"
)

B3B <- file.path(
  ROOT,
  "results",
  "V2_34B3B_CRA002354_GENUS_COVERAGE_ROBUSTNESS"
)

OUT <- file.path(
  ROOT,
  "results",
  "V2_34B3C_CRA002354_GENUS_EVIDENCE_FREEZE"
)

dir.create(
  OUT,
  recursive = TRUE,
  showWarnings = FALSE
)

if (
  !file.exists(
    file.path(
      B3A,
      "_STEP94B3A_COMPLETE.ok"
    )
  )
) {
  stop("Step94B3A is not complete.")
}

if (
  !file.exists(
    file.path(
      B3B,
      "_STEP94B3B_COMPLETE.ok"
    )
  )
) {
  stop("Step94B3B is not complete.")
}

base <- read_csv(
  file.path(
    B3A,
    "04_EXPLORATORY_BASELINE_GENUS_SOURCE_CONTRAST.csv"
  ),
  show_col_types = FALSE
)

long <- read_csv(
  file.path(
    B3A,
    "06_EXPLORATORY_LONGITUDINAL_GENUS_SOURCE_TIME_INTERACTIONS.csv"
  ),
  show_col_types = FALSE
)

abd <- read_csv(
  file.path(
    B3A,
    "09_DESCRIPTIVE_PULMONARY_VS_ABDOMINAL_GI_BASELINE_GENUS.csv"
  ),
  show_col_types = FALSE
)

b3a_sum <- read_csv(
  file.path(
    B3A,
    "12_GENUS_SOURCE_EVIDENCE_SUMMARY.csv"
  ),
  show_col_types = FALSE
)

rob <- read_csv(
  file.path(
    B3B,
    "04_PRIMARY_BASELINE_HIT_COVERAGE_ROBUSTNESS.csv"
  ),
  show_col_types = FALSE
)

b3b_sum <- read_csv(
  file.path(
    B3B,
    "05_GENUS_COVERAGE_ROBUSTNESS_SUMMARY.csv"
  ),
  show_col_types = FALSE
)

long_cov <- read_csv(
  file.path(
    B3B,
    "03_LONGITUDINAL_GENUS_COVERAGE_SENSITIVITY.csv"
  ),
  show_col_types = FALSE
)

# ------------------------------------------------------------
# 1. Final baseline hit registry
# ------------------------------------------------------------

baseline_hits <- rob %>%
  mutate(
    final_interpretation = case_when(
      robust_FDR_70 &
        robust_FDR_80 ~
        "ROBUST_ACROSS_70_AND_80_PERCENT_COVERAGE",

      robust_FDR_70 &
        !robust_FDR_80 ~
        "PARTIAL_SUPPORT_AT_70_PERCENT_ONLY",

      !robust_FDR_70 &
        !robust_FDR_80 ~
        "PRIMARY_ONLY_NOT_COVERAGE_ROBUST",

      TRUE ~
        "MIXED"
    ),
    reporting_priority = case_when(
      genus == "Enterococcus" ~
        "BIOLOGICALLY_INTERPRETABLE_BUT_NOT_COVERAGE_ROBUST",
      robust_FDR_70 &
        !low_abundance_both_groups_lt_0_1pct ~
        "SUPPLEMENTARY_PRIORITY",
      robust_FDR_70 &
        low_abundance_both_groups_lt_0_1pct ~
        "LOW_ABUNDANCE_SUPPLEMENTARY",
      TRUE ~
        "DEEMPHASIZE"
    )
  )

write_csv(
  baseline_hits,
  file.path(
    OUT,
    "01_FINAL_BASELINE_GENUS_HIT_REGISTRY.csv"
  )
)

# ------------------------------------------------------------
# 2. Longitudinal primary/sensitivity registry
# ------------------------------------------------------------

primary_long_hits <- long %>%
  filter(
    !is.na(FDR),
    FDR < 0.05
  )

cov70_hits <- long_cov %>%
  filter(
    threshold == 0.70,
    !is.na(FDR),
    FDR < 0.05
  ) %>%
  transmute(
    genus,
    beta_70 = beta,
    p_70 = p_LRT,
    FDR_70 = FDR
  )

cov80_hits <- long_cov %>%
  filter(
    threshold == 0.80,
    !is.na(FDR),
    FDR < 0.05
  ) %>%
  transmute(
    genus,
    beta_80 = beta,
    p_80 = p_LRT,
    FDR_80 = FDR
  )

sensitivity_only <- full_join(
  cov70_hits,
  cov80_hits,
  by = "genus"
) %>%
  mutate(
    primary_FDR_significant =
      genus %in%
      primary_long_hits$genus,
    sensitivity_class = case_when(
      !is.na(FDR_70) &
        !is.na(FDR_80) ~
        "SIGNIFICANT_AT_BOTH_COVERAGE_THRESHOLDS",
      !is.na(FDR_70) ~
        "SIGNIFICANT_AT_70_PERCENT_ONLY",
      !is.na(FDR_80) ~
        "SIGNIFICANT_AT_80_PERCENT_ONLY",
      TRUE ~
        "NONE"
    ),
    interpretation =
      "SENSITIVITY_ONLY_DO_NOT_UPGRADE_PRIMARY_LONGITUDINAL_CONCLUSION"
  )

write_csv(
  sensitivity_only,
  file.path(
    OUT,
    "02_LONGITUDINAL_COVERAGE_SENSITIVITY_ONLY_SIGNALS.csv"
  )
)

# ------------------------------------------------------------
# 3. Final evidence summary
# ------------------------------------------------------------

n_primary_base <-
  sum(
    base$FDR < 0.05,
    na.rm = TRUE
  )

n_primary_long <-
  sum(
    long$FDR < 0.05,
    na.rm = TRUE
  )

n_rob70 <-
  sum(
    baseline_hits$robust_FDR_70,
    na.rm = TRUE
  )

n_rob80 <-
  sum(
    baseline_hits$robust_FDR_80,
    na.rm = TRUE
  )

n_robboth <-
  sum(
    baseline_hits$robust_FDR_70 &
    baseline_hits$robust_FDR_80,
    na.rm = TRUE
  )

n_long70 <-
  sum(
    long_cov$threshold == 0.70 &
    long_cov$FDR < 0.05,
    na.rm = TRUE
  )

n_long80 <-
  sum(
    long_cov$threshold == 0.80 &
    long_cov$FDR < 0.05,
    na.rm = TRUE
  )

n_abd <-
  sum(
    abd$FDR < 0.05,
    na.rm = TRUE
  )

tier <-
  "EXPLORATORY_BASELINE_GENUS_SIGNALS_WITH_LIMITED_COVERAGE_ROBUSTNESS_AND_NO_PRIMARY_LONGITUDINAL_SIGNAL"

summary <- tibble(
  cohort = "CRA002354",
  branch =
    "INFECTION_SOURCE_EXPLORATORY_GENUS",
  abundance_weighted_genus_assignment =
    b3a_sum$abundance_weighted_genus_assignment[1],
  baseline_genera_tested =
    b3a_sum$baseline_genera_tested[1],
  primary_baseline_FDR_hits =
    n_primary_base,
  baseline_hits_FDR_robust_at_70pct =
    n_rob70,
  baseline_hits_FDR_robust_at_80pct =
    n_rob80,
  baseline_hits_FDR_robust_at_both =
    n_robboth,
  longitudinal_genera_tested =
    b3a_sum$longitudinal_genera_tested[1],
  primary_longitudinal_FDR_hits =
    n_primary_long,
  longitudinal_FDR_hits_at_70pct_coverage =
    n_long70,
  longitudinal_FDR_hits_at_80pct_coverage =
    n_long80,
  pulmonary_vs_abdominal_GI_FDR_hits =
    n_abd,
  evidence_tier =
    tier,
  branch_role =
    "EXPLORATORY_SUPPORTIVE_ONLY",
  branch_frozen =
    TRUE
)

write_csv(
  summary,
  file.path(
    OUT,
    "03_FINAL_GENUS_BRANCH_EVIDENCE_FREEZE.csv"
  )
)

# ------------------------------------------------------------
# 4. Reporting recommendation
# ------------------------------------------------------------

report <- tibble(
  genus =
    baseline_hits$genus,
  primary_FDR =
    baseline_hits$primary_FDR,
  robust_FDR_70 =
    baseline_hits$robust_FDR_70,
  robust_FDR_80 =
    baseline_hits$robust_FDR_80,
  low_abundance_both_groups_lt_0_1pct =
    baseline_hits$low_abundance_both_groups_lt_0_1pct,
  final_interpretation =
    baseline_hits$final_interpretation,
  recommended_location = case_when(
    baseline_hits$robust_FDR_70 &
      !baseline_hits$low_abundance_both_groups_lt_0_1pct ~
      "SUPPLEMENTARY_RESULTS",

    baseline_hits$robust_FDR_70 &
      baseline_hits$low_abundance_both_groups_lt_0_1pct ~
      "SUPPLEMENTARY_TABLE_ONLY",

    baseline_hits$genus ==
      "Enterococcus" ~
      "OPTIONAL_DISCUSSION_ONLY_WITH_ROBUSTNESS_CAVEAT",

    TRUE ~
      "SUPPLEMENTARY_TABLE_ONLY"
  )
)

write_csv(
  report,
  file.path(
    OUT,
    "04_GENUS_REPORTING_RECOMMENDATION.csv"
  )
)

# ------------------------------------------------------------
# 5. Frozen wording + guardrails
# ------------------------------------------------------------

writeLines(
  c(
    "STEP94B3C FROZEN GENUS-BRANCH WORDING",
    "",
    "Primary exploratory genus analysis:",
    paste0(
      n_primary_base,
      " of ",
      b3a_sum$baseline_genera_tested[1],
      " eligible baseline genera were significant after BH-FDR correction, whereas none of ",
      b3a_sum$longitudinal_genera_tested[1],
      " longitudinal genus-level source-by-time interactions survived FDR correction."
    ),
    "",
    "Coverage robustness:",
    paste0(
      "Of the ",
      n_primary_base,
      " baseline FDR-significant genera, ",
      n_rob70,
      " remained significant after restricting to samples with at least 70% genus-level taxonomic coverage, and none remained significant at the 80% threshold."
    ),
    "",
    "Recommended interpretation:",
    "The genus-level findings are therefore best interpreted as exploratory baseline taxonomic differences with limited robustness to sample-level taxonomy coverage, rather than a stable infection-source-specific taxonomic signature.",
    "",
    "Longitudinal sensitivity:",
    paste0(
      n_long70,
      " genus-level source-by-time signal(s) reached FDR significance only in the >=70% coverage sensitivity and ",
      n_long80,
      " at the >=80% sensitivity. Because the prespecified primary B3A longitudinal family contained no FDR-significant genera, sensitivity-only signals must not upgrade the primary longitudinal conclusion."
    ),
    "",
    "Specific taxa:",
    "Roseateles and Enhydrobacter retained FDR significance at the >=70% coverage threshold but not at >=80%.",
    "Caldimonas and Enterococcus were not FDR-robust after coverage filtering.",
    "Enterococcus had the largest biologically interpretable abundance contrast among the primary hits, but its statistical robustness to coverage filtering was limited.",
    "Bacteroides emerged only in the >=70% longitudinal coverage sensitivity and was not FDR-significant at >=80%; it should be treated as a coverage-sensitive exploratory signal only.",
    "",
    "Do not write:",
    "- a reproducible pulmonary-sepsis genus signature was identified",
    "- Enterococcus is a robust infection-source biomarker",
    "- Bacteroides shows a validated source-specific longitudinal trajectory",
    "- infection source determines genus-level microbiome evolution",
    "- baseline genus differences establish causal source effects",
    "",
    paste0(
      "Final evidence tier: ",
      tier
    )
  ),
  file.path(
    OUT,
    "05_FROZEN_GENUS_MANUSCRIPT_WORDING_AND_GUARDRAILS.txt"
  )
)

# ------------------------------------------------------------
# 6. Copy authoritative source outputs
# ------------------------------------------------------------

copy_map <- c(
  "B3A_04_BASELINE.csv" =
    file.path(
      B3A,
      "04_EXPLORATORY_BASELINE_GENUS_SOURCE_CONTRAST.csv"
    ),
  "B3A_06_LONGITUDINAL.csv" =
    file.path(
      B3A,
      "06_EXPLORATORY_LONGITUDINAL_GENUS_SOURCE_TIME_INTERACTIONS.csv"
    ),
  "B3A_09_ABDOMINAL_GI.csv" =
    file.path(
      B3A,
      "09_DESCRIPTIVE_PULMONARY_VS_ABDOMINAL_GI_BASELINE_GENUS.csv"
    ),
  "B3A_12_SUMMARY.csv" =
    file.path(
      B3A,
      "12_GENUS_SOURCE_EVIDENCE_SUMMARY.csv"
    ),
  "B3B_04_HIT_ROBUSTNESS.csv" =
    file.path(
      B3B,
      "04_PRIMARY_BASELINE_HIT_COVERAGE_ROBUSTNESS.csv"
    ),
  "B3B_05_SUMMARY.csv" =
    file.path(
      B3B,
      "05_GENUS_COVERAGE_ROBUSTNESS_SUMMARY.csv"
    )
)

audit <- tibble(
  frozen_name =
    names(copy_map),
  source =
    unname(copy_map),
  copied =
    FALSE
)

for (
  i in seq_along(copy_map)
) {

  src <- unname(
    copy_map[i]
  )

  dst <- file.path(
    OUT,
    names(copy_map)[i]
  )

  if (
    file.exists(src)
  ) {
    file.copy(
      src,
      dst,
      overwrite = TRUE
    )
    audit$copied[i] <- TRUE
  }
}

write_csv(
  audit,
  file.path(
    OUT,
    "06_FREEZE_COPY_AUDIT.csv"
  )
)

writeLines(
  c(
    "STEP94B3C GENUS EVIDENCE FREEZE",
    "",
    paste0(
      "Final evidence tier: ",
      tier
    ),
    "",
    "This branch is exploratory/supportive only.",
    "No further genus-level infection-source model is required for the current V2 manuscript.",
    "STEP94B3C is the authoritative final freeze for CRA002354 genus-level source findings."
  ),
  file.path(
    OUT,
    "07_STEP94B3C_INTERPRETATION.txt"
  )
)

writeLines(
  c(
    paste0(
      "Completed: ",
      Sys.time()
    ),
    "Branch frozen: TRUE",
    "STEP94B3C COMPLETE"
  ),
  file.path(
    OUT,
    "_STEP94B3C_COMPLETE.ok"
  )
)

cat(
  "STEP94B3C COMPLETE\n"
)
