# ============================================================
# Sepsis V2 - Step93W
# LONGITUDINAL EVIDENCE FREEZE
#
# Purpose:
# Freeze the longitudinal ecological evidence after Step93V3 reconciliation.
#
# Combines, WITHOUT mixing estimands:
#   1) primary all-patient early->late paired Bray contrast
#   2) common-anchor paired sensitivity
#   3) categorical global time model
#   4) continuous-day model
#   5) anchor quality
#   6) formal interpretation tier
#   7) Step93U core-sepsis taxonomic trajectory evidence
#
# This becomes the manuscript-facing evidence matrix.
# ============================================================

options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(stringr)
  library(tibble)
  library(tidyr)
})

ROOT <- "E:/sepsis_project"

IN_V3 <- file.path(
  ROOT,
  "results",
  "V2_33V3_STEP88B_ESTIMAND_RECONCILIATION"
)

IN_U <- file.path(
  ROOT,
  "results",
  "V2_33U_TRAJECTORY_ECOLOGY_ROBUSTNESS"
)

OUT <- file.path(
  ROOT,
  "results",
  "V2_33W_LONGITUDINAL_EVIDENCE_FREEZE"
)

dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

TARGETS <- c(
  "PRJNA691455",
  "PRJEB82425",
  "PRJNA516701",
  "PRJNA851469",
  "PRJNA578267",
  "PRJNA1166732",
  "PRJNA430161"
)

role_map <- tibble(
  project = TARGETS,
  cohort_role = c(
    "CORE_SEPSIS_LONGITUDINAL",
    "ICU_INFECTION_LONGITUDINAL_EXTERNAL",
    "ICU_BACKGROUND_LONGITUDINAL",
    "ICU_BACKGROUND_LONGITUDINAL",
    "NONSEPSIS_LONGITUDINAL_CONTROL",
    "INTERVENTION_LONGITUDINAL_SUPPORT",
    "INTERVENTION_LONGITUDINAL_SUPPORT"
  )
)

required_v3 <- c(
  "02_ACTUAL_early_vs_late_Bray_contrasts.csv",
  "03_ACTUAL_Bray_direction_robustness_registry.csv",
  "05_ACTUAL_anchor_audit.csv",
  "06_ACTUAL_categorical_global_models_Bray.csv",
  "07_ACTUAL_continuous_day_models_Bray.csv"
)

for (f in required_v3) {
  if (!file.exists(file.path(IN_V3, f))) {
    stop(paste0("Missing Step93V3 file: ", f))
  }
}

primary <- read_csv(
  file.path(IN_V3, required_v3[1]),
  show_col_types = FALSE
)

registry <- read_csv(
  file.path(IN_V3, required_v3[2]),
  show_col_types = FALSE
)

anchor <- read_csv(
  file.path(IN_V3, required_v3[3]),
  show_col_types = FALSE
)

catmod <- read_csv(
  file.path(IN_V3, required_v3[4]),
  show_col_types = FALSE
)

contmod <- read_csv(
  file.path(IN_V3, required_v3[5]),
  show_col_types = FALSE
)

# ------------------------------------------------------------
# 1. Standardize all evidence sources
# ------------------------------------------------------------

primary2 <- primary %>%
  select(
    project,
    primary_n_pairs = n_pairs,
    primary_dz = paired_effect_dz,
    primary_p = wilcoxon_p,
    primary_fdr = fdr,
    primary_direction = early_late_direction,
    primary_evidence = early_late_evidence
  )

registry2 <- registry %>%
  transmute(
    project,
    analysis_role_group,
    required_anchor,
    required_anchor_fraction,
    required_anchor_flag,
    all_n_pairs = n_pairs_all,
    all_mean_diff = bray_mean_diff_all,
    all_dz = bray_dz_all,
    all_p = bray_p_all,
    all_fdr = bray_fdr_all,
    all_direction = direction_all,
    common_anchor_n_pairs = n_pairs_common_anchor,
    common_anchor_mean_diff = bray_mean_diff_common_anchor,
    common_anchor_dz = bray_dz_common_anchor,
    common_anchor_p = bray_p_common_anchor,
    common_anchor_fdr = bray_fdr_common_anchor,
    common_anchor_direction = direction_common_anchor,
    direction_robust,
    formal_interpretation_tier
  )

anchor2 <- anchor %>%
  select(
    project,
    early_followup,
    late_followup,
    anchor_structure,
    n_reference_patients_total,
    n_required_anchor_patients
  )

cat2 <- catmod %>%
  transmute(
    project,
    categorical_global_time_p = global_time_p,
    categorical_global_time_fdr = global_time_fdr_within_outcome
  )

cont2 <- contmod %>%
  transmute(
    project,
    continuous_day_global_time_p = global_time_p,
    continuous_day_global_time_fdr = global_time_fdr_within_outcome
  )

evidence <- role_map %>%
  left_join(primary2, by = "project") %>%
  left_join(registry2, by = "project") %>%
  left_join(anchor2, by = "project") %>%
  left_join(cat2, by = "project") %>%
  left_join(cont2, by = "project")

# ------------------------------------------------------------
# 2. Evidence classification
# ------------------------------------------------------------

evidence <- evidence %>%
  mutate(
    primary_significant = !is.na(primary_fdr) & primary_fdr < 0.05,
    common_anchor_significant = !is.na(common_anchor_fdr) & common_anchor_fdr < 0.05,
    categorical_time_significant =
      !is.na(categorical_global_time_fdr) &
      categorical_global_time_fdr < 0.05,
    continuous_day_significant =
      !is.na(continuous_day_global_time_fdr) &
      continuous_day_global_time_fdr < 0.05,

    frozen_evidence_class = case_when(

      project %in% c("PRJNA691455", "PRJNA851469", "PRJNA516701") &
        primary_significant &
        primary_dz > 0 &
        direction_robust ~
        "ROBUST_PROGRESSIVE_DISPLACEMENT",

      project == "PRJEB82425" &
        primary_significant &
        primary_dz > 0 ~
        "SUPPORTIVE_PROGRESSIVE_DISPLACEMENT_HETEROGENEOUS_ANCHOR",

      project == "PRJNA578267" &
        primary_dz < 0 &
        !primary_significant &
        common_anchor_significant &
        common_anchor_dz < 0 &
        categorical_time_significant ~
        "PRIMARY_RECOVERY_TREND_SENSITIVITY_SUPPORTED_RECONVERGENCE",

      project %in% c("PRJNA1166732", "PRJNA430161") &
        !primary_significant &
        !common_anchor_significant ~
        "NO_CLEAR_PAIRED_CHANGE_INTERVENTION_SUPPORT",

      TRUE ~
        "MIXED_OR_UNRESOLVED"
    ),

    frozen_manuscript_wording = case_when(

      frozen_evidence_class == "ROBUST_PROGRESSIVE_DISPLACEMENT" ~
        "Progressive within-patient ecological displacement was supported by the primary early-to-late contrast and was directionally robust to common-anchor sensitivity analysis.",

      frozen_evidence_class ==
        "SUPPORTIVE_PROGRESSIVE_DISPLACEMENT_HETEROGENEOUS_ANCHOR" ~
        "Progressive displacement was observed, but heterogeneous baseline anchoring limits this cohort to supportive external evidence.",

      frozen_evidence_class ==
        "PRIMARY_RECOVERY_TREND_SENSITIVITY_SUPPORTED_RECONVERGENCE" ~
        "The primary early-to-late contrast showed a non-significant recovery-direction trend, while common-anchor sensitivity analysis and the global time model supported re-convergence toward baseline.",

      frozen_evidence_class ==
        "NO_CLEAR_PAIRED_CHANGE_INTERVENTION_SUPPORT" ~
        "No clear early-to-late paired ecological displacement was supported in this intervention-support cohort.",

      TRUE ~
        "Evidence remains mixed or unresolved."
    )
  )

write_csv(
  evidence,
  file.path(OUT, "01_FROZEN_longitudinal_evidence_matrix.csv")
)

# ------------------------------------------------------------
# 3. Main manuscript table
# ------------------------------------------------------------

manuscript_table <- evidence %>%
  transmute(
    project,
    cohort_role,
    anchor_quality = required_anchor_flag,
    early_followup,
    late_followup,
    primary_n_pairs,
    primary_dz,
    primary_p,
    primary_fdr,
    primary_direction,
    common_anchor_n_pairs,
    common_anchor_dz,
    common_anchor_p,
    common_anchor_fdr,
    common_anchor_direction,
    categorical_global_time_fdr,
    continuous_day_global_time_fdr,
    formal_interpretation_tier,
    frozen_evidence_class
  )

write_csv(
  manuscript_table,
  file.path(OUT, "02_MANUSCRIPT_longitudinal_evidence_table.csv")
)

# ------------------------------------------------------------
# 4. Primary vs common-anchor comparison table
# ------------------------------------------------------------

comparison <- evidence %>%
  transmute(
    project,
    cohort_role,
    required_anchor_flag,
    primary_dz,
    primary_fdr,
    common_anchor_dz,
    common_anchor_fdr,
    effect_shift_common_minus_primary =
      common_anchor_dz - primary_dz,
    direction_same =
      sign(primary_dz) == sign(common_anchor_dz),
    direction_robust,
    formal_interpretation_tier
  )

write_csv(
  comparison,
  file.path(OUT, "03_PRIMARY_vs_COMMON_ANCHOR_comparison.csv")
)

# ------------------------------------------------------------
# 5. Main + sensitivity figure
# ------------------------------------------------------------

plot_df <- evidence %>%
  filter(!is.na(primary_dz)) %>%
  arrange(primary_dz)

make_plot <- function() {

  par(mar = c(5, 14, 3, 2))

  y <- seq_len(nrow(plot_df))

  lim <- range(
    c(
      plot_df$primary_dz,
      plot_df$common_anchor_dz,
      0
    ),
    na.rm = TRUE
  )

  plot(
    plot_df$primary_dz,
    y,
    xlim = lim,
    ylim = c(0.5, nrow(plot_df) + 0.5),
    yaxt = "n",
    ylab = "",
    xlab = "Paired standardized Bray-Curtis change (dz)",
    pch = 19,
    main = "Primary and common-anchor sensitivity contrasts"
  )

  abline(v = 0, lty = 2)

  for (i in seq_len(nrow(plot_df))) {
    if (!is.na(plot_df$common_anchor_dz[i])) {
      segments(
        plot_df$primary_dz[i],
        y[i],
        plot_df$common_anchor_dz[i],
        y[i],
        lty = 3
      )
    }
  }

  points(
    plot_df$common_anchor_dz,
    y,
    pch = 1,
    cex = 1.2
  )

  axis(
    2,
    at = y,
    labels = paste0(
      plot_df$project,
      "  [",
      plot_df$required_anchor_flag,
      "]"
    ),
    las = 1,
    cex.axis = 0.78
  )

  legend(
    "bottomright",
    legend = c(
      "Primary early-to-late",
      "Common-anchor sensitivity"
    ),
    pch = c(19, 1),
    bty = "n"
  )
}

pdf(
  file.path(OUT, "Figure_STEP93W_primary_vs_common_anchor.pdf"),
  width = 9.5,
  height = 6.5
)
make_plot()
dev.off()

png(
  file.path(OUT, "Figure_STEP93W_primary_vs_common_anchor.png"),
  width = 1900,
  height = 1300,
  res = 180
)
make_plot()
dev.off()

# ------------------------------------------------------------
# 6. Integrate Step93U core taxonomic trajectory evidence
# ------------------------------------------------------------

u_files <- c(
  heterogeneity = file.path(
    IN_U,
    "06_trajectory_heterogeneity_summary.csv"
  ),
  robustness = file.path(
    IN_U,
    "09_trajectory_ecology_robustness.csv"
  ),
  paired_bray = file.path(
    IN_U,
    "04_paired_Day3_Day7_Bray_test.csv"
  )
)

for (f in u_files) {
  if (!file.exists(f)) {
    stop(paste0("Missing Step93U source: ", f))
  }
}

u_hetero <- read_csv(u_files[["heterogeneity"]], show_col_types = FALSE)
u_robust <- read_csv(u_files[["robustness"]], show_col_types = FALSE)
u_bray <- read_csv(u_files[["paired_bray"]], show_col_types = FALSE)

taxonomic_core <- tibble(
  metric = c(
    "PRJNA691455 Day3 Bray median",
    "PRJNA691455 Day7 Bray median",
    "PRJNA691455 Day3-Day7 paired p",
    "Signed genus trajectory features",
    "Median pairwise cosine similarity",
    "Proportion pairwise cosine <= 0"
  ),
  value = c(
    u_bray$median_Day3[1],
    u_bray$median_Day7[1],
    u_bray$wilcoxon_p[1],
    u_hetero$n_features[1],
    u_hetero$median_cosine_similarity[1],
    u_hetero$proportion_cosine_le_0[1]
  )
)

write_csv(
  taxonomic_core,
  file.path(OUT, "04_CORE_SEPSIS_taxonomic_trajectory_frozen_summary.csv")
)

# ------------------------------------------------------------
# 7. Frozen claim guardrails
# ------------------------------------------------------------

guardrails <- tibble(
  claim = c(
    "Progressive ecological displacement is supported in the core sepsis cohort",
    "Progressive ecological displacement is supported in both ICU-background cohorts",
    "PRJEB82425 provides supportive external ICU-infection evidence",
    "The non-sepsis control has established recovery in the primary analysis",
    "The non-sepsis control shows a recovery-direction pattern strengthened by sensitivity analysis",
    "Sepsis-specific ecological instability is established",
    "Two taxonomic trajectory subtypes are established",
    "EII independently validates taxonomic trajectories"
  ),
  status = c(
    "YES",
    "YES",
    "YES_WITH_ANCHOR_LIMITATION",
    "NO",
    "YES",
    "NO",
    "NO",
    "NO"
  ),
  reason = c(
    "Primary and common-anchor paired Bray effects are positive and statistically supported.",
    "Both cohorts have positive primary paired effects; common-anchor direction is concordant.",
    "Large positive displacement is present, but only 48.6% share the required Inclusion anchor.",
    "Primary early-to-late contrast is negative but not FDR-significant.",
    "Common-anchor sensitivity is significantly negative and the categorical time model is significant.",
    "ICU-background cohorts show similar progressive displacement.",
    "k=2 solution was 9 vs 1 and is not a validated subtype structure.",
    "EII includes Bray/alpha-instability components and is not an independent endpoint."
  )
)

write_csv(
  guardrails,
  file.path(OUT, "05_FROZEN_claim_guardrails.csv")
)

# ------------------------------------------------------------
# 8. Manuscript-ready longitudinal Results paragraph
# ------------------------------------------------------------

getrow <- function(p) evidence %>% filter(project == p)

core <- getrow("PRJNA691455")
bg1 <- getrow("PRJNA851469")
bg2 <- getrow("PRJNA516701")
ctrl <- getrow("PRJNA578267")
ext <- getrow("PRJEB82425")

results_text <- c(
  "FROZEN LONGITUDINAL RESULTS TEXT",
  "",
  paste0(
    "In the core sepsis cohort (PRJNA691455), within-patient Bray-Curtis displacement increased from the early to late follow-up interval (n=",
    core$primary_n_pairs,
    ", dz=",
    formatC(core$primary_dz, digits=3, format="f"),
    ", FDR=",
    formatC(core$primary_fdr, digits=3, format="f"),
    "), with concordant direction in the common-anchor analysis."
  ),
  paste0(
    "The same progressive-displacement direction was independently observed in both ICU-background cohorts (PRJNA851469: dz=",
    formatC(bg1$primary_dz, digits=3, format="f"),
    ", FDR=",
    formatC(bg1$primary_fdr, digits=3, format="f"),
    "; PRJNA516701: dz=",
    formatC(bg2$primary_dz, digits=3, format="f"),
    ", FDR=",
    formatC(bg2$primary_fdr, digits=3, format="f"),
    ")."
  ),
  paste0(
    "The external ICU-infection cohort PRJEB82425 also showed strong progressive displacement (dz=",
    formatC(ext$primary_dz, digits=3, format="f"),
    ", FDR=",
    formatC(ext$primary_fdr, digits=4, format="f"),
    "), although its heterogeneous baseline anchoring limits this cohort to supportive sensitivity evidence."
  ),
  paste0(
    "In contrast, the non-sepsis longitudinal control PRJNA578267 showed a negative early-to-late direction in the primary analysis (dz=",
    formatC(ctrl$primary_dz, digits=3, format="f"),
    ", FDR=",
    formatC(ctrl$primary_fdr, digits=3, format="f"),
    "), which was not statistically significant. However, restriction to patients sharing the common baseline anchor strengthened the negative effect (dz=",
    formatC(ctrl$common_anchor_dz, digits=3, format="f"),
    ", FDR=",
    formatC(ctrl$common_anchor_fdr, digits=3, format="f"),
    "), supporting re-convergence toward baseline as a sensitivity finding."
  ),
  "",
  "Interpretation:",
  "The longitudinal evidence supports progressive ecological displacement during critical illness rather than a sepsis-specific phenomenon. The non-sepsis control did not show the same progressive pattern and showed sensitivity-supported re-convergence toward baseline."
)

writeLines(
  results_text,
  file.path(OUT, "06_FROZEN_longitudinal_results_paragraph.txt")
)

# ------------------------------------------------------------
# 9. Methods wording
# ------------------------------------------------------------

methods_text <- c(
  "FROZEN LONGITUDINAL METHODS WORDING",
  "",
  "The prespecified primary paired estimand compared the designated early and late follow-up time points among patients with both observations.",
  "A separate common-anchor sensitivity analysis restricted each cohort to patients sharing the required baseline/reference stage.",
  "Categorical global-time mixed models and, where meaningful, continuous-day mixed models were retained as complementary longitudinal tests.",
  "Because cohort definitions, reference stages, and follow-up schedules were heterogeneous, effect estimates were not meta-pooled across cohorts.",
  "Results from distinct estimands were reported separately and were not combined within a single effect estimate."
)

writeLines(
  methods_text,
  file.path(OUT, "07_FROZEN_longitudinal_methods_wording.txt")
)

# ------------------------------------------------------------
# 10. QC
# ------------------------------------------------------------

qc <- tibble(
  cohorts = nrow(evidence),
  robust_progressive = sum(
    evidence$frozen_evidence_class == "ROBUST_PROGRESSIVE_DISPLACEMENT"
  ),
  supportive_external = sum(
    evidence$frozen_evidence_class ==
      "SUPPORTIVE_PROGRESSIVE_DISPLACEMENT_HETEROGENEOUS_ANCHOR"
  ),
  sensitivity_supported_reconvergence = sum(
    evidence$frozen_evidence_class ==
      "PRIMARY_RECOVERY_TREND_SENSITIVITY_SUPPORTED_RECONVERGENCE"
  ),
  intervention_no_clear_change = sum(
    evidence$frozen_evidence_class ==
      "NO_CLEAR_PAIRED_CHANGE_INTERVENTION_SUPPORT"
  )
)

write_csv(
  qc,
  file.path(OUT, "08_STEP93W_QC_SUMMARY.csv")
)

writeLines(
  c(
    paste0("Completed: ", Sys.time()),
    "Primary and sensitivity estimands retained separately.",
    "LONGITUDINAL EVIDENCE FROZEN.",
    "STEP93W COMPLETE"
  ),
  file.path(OUT, "_STEP93W_COMPLETE.ok")
)

cat("STEP93W COMPLETE\n")
