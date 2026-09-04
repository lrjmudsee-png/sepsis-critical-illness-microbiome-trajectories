# ============================================================
# Sepsis V2 - Step93V3
# STEP88B ESTIMAND RECONCILIATION
#
# Purpose:
# Reconcile the actual Step88B outputs WITHOUT using remembered/fallback
# numbers. This prevents mixing different estimands:
#
#   A) early -> late paired follow-up contrast
#   B) anchor robustness / Bray direction registry
#   C) categorical global time model
#   D) continuous-day model
#
# Primary goal:
# Decide which evidence supports:
# - progressive displacement
# - recovery/re-convergence
# - direction-only / no-clear-change
#
# No cherry-picking and no cross-estimand substitution.
# ============================================================

options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(stringr)
  library(tibble)
})

ROOT <- "E:/sepsis_project"

BASE <- file.path(
  ROOT,
  "results",
  "V2_28B_STEP88B_ANCHOR_ROBUSTNESS_AND_PAIRED_CONTRASTS"
)

OUT <- file.path(
  ROOT,
  "results",
  "V2_33V3_STEP88B_ESTIMAND_RECONCILIATION"
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

FILES <- c(
  anchor_audit = file.path(
    BASE,
    "01_ANCHOR_AUDIT",
    "V2_STEP88B_anchor_audit.csv"
  ),
  paired_early_late = file.path(
    BASE,
    "02_PAIRED_CONTRASTS",
    "V2_STEP88B_paired_early_vs_late_followup_contrasts.csv"
  ),
  categorical_models = file.path(
    BASE,
    "03_MODEL_FDR",
    "V2_STEP88B_PRIMARY_categorical_global_models_with_FDR.csv"
  ),
  continuous_models = file.path(
    BASE,
    "03_MODEL_FDR",
    "V2_STEP88B_SECONDARY_continuous_day_models_with_FDR.csv"
  ),
  bray_registry = file.path(
    BASE,
    "04_SYNTHESIS",
    "V2_STEP88B_Bray_direction_and_robustness_registry.csv"
  )
)

safe_read <- function(f) {
  if (!file.exists(f)) return(NULL)
  tryCatch(
    read_csv(f, show_col_types = FALSE, guess_max = 20000, name_repair = "unique"),
    error = function(e) NULL
  )
}

norm_project <- function(x) {
  toupper(str_trim(as.character(x)))
}

find_col <- function(nms, patterns) {
  for (p in patterns) {
    hit <- nms[str_detect(nms, regex(p, ignore_case = TRUE))]
    if (length(hit) > 0) return(hit[1])
  }
  NA_character_
}

# ------------------------------------------------------------
# 1. Source file audit
# ------------------------------------------------------------

source_audit <- lapply(names(FILES), function(nm) {
  f <- FILES[[nm]]
  x <- safe_read(f)

  tibble(
    source = nm,
    file = f,
    exists = file.exists(f),
    readable = !is.null(x),
    rows = if (is.null(x)) NA_integer_ else nrow(x),
    cols = if (is.null(x)) NA_integer_ else ncol(x),
    columns = if (is.null(x)) NA_character_ else paste(names(x), collapse = ";")
  )
}) %>% bind_rows()

write_csv(
  source_audit,
  file.path(OUT, "01_STEP88B_source_file_audit.csv")
)

if (!all(source_audit$readable)) {
  warning("One or more Step88B source files are missing/unreadable. Audit written.")
}

# ------------------------------------------------------------
# 2. Read all available source tables
# ------------------------------------------------------------

anchor <- safe_read(FILES[["anchor_audit"]])
paired <- safe_read(FILES[["paired_early_late"]])
catmod <- safe_read(FILES[["categorical_models"]])
contmod <- safe_read(FILES[["continuous_models"]])
registry <- safe_read(FILES[["bray_registry"]])

if (is.null(paired)) {
  stop("Actual Step88B paired early-vs-late table is missing.")
}

# ------------------------------------------------------------
# 3. Exact early -> late paired Bray table
# ------------------------------------------------------------

p_project <- find_col(
  names(paired),
  c("^project$", "bioproject", "cohort")
)

p_effect <- find_col(
  names(paired),
  c("^paired_effect_dz$", "^dz$", "effect.*dz", "cohen")
)

p_n <- find_col(
  names(paired),
  c("^n_pairs$", "paired.*n", "n_pair")
)

p_p <- find_col(
  names(paired),
  c("^wilcoxon_p$", "wilcoxon.*p", "^p_value$", "^p$")
)

p_fdr <- find_col(
  names(paired),
  c("^wilcoxon_fdr_within_metric$", "fdr", "padj", "q_value")
)

p_metric <- find_col(
  names(paired),
  c("^metric$", "outcome", "measure")
)

if (is.na(p_project) || is.na(p_effect)) {
  stop("Cannot identify project/effect columns in actual paired table.")
}

paired2 <- paired

if (!is.na(p_metric)) {
  bray_rows <- str_detect(
    tolower(as.character(paired2[[p_metric]])),
    "bray"
  )
  if (any(bray_rows, na.rm = TRUE)) {
    paired2 <- paired2[bray_rows %in% TRUE, , drop = FALSE]
  }
}

early_late <- tibble(
  project = norm_project(paired2[[p_project]]),
  n_pairs = if (!is.na(p_n)) suppressWarnings(as.numeric(paired2[[p_n]])) else NA_real_,
  paired_effect_dz = suppressWarnings(as.numeric(paired2[[p_effect]])),
  wilcoxon_p = if (!is.na(p_p)) suppressWarnings(as.numeric(paired2[[p_p]])) else NA_real_,
  fdr = if (!is.na(p_fdr)) suppressWarnings(as.numeric(paired2[[p_fdr]])) else NA_real_
) %>%
  filter(project %in% TARGETS) %>%
  distinct(project, .keep_all = TRUE) %>%
  left_join(role_map, by = "project") %>%
  mutate(
    early_late_direction = case_when(
      paired_effect_dz > 0 ~ "increasing_displacement",
      paired_effect_dz < 0 ~ "recovery_direction",
      TRUE ~ "neutral"
    ),
    early_late_evidence = case_when(
      !is.na(fdr) & fdr < 0.05 & paired_effect_dz > 0 ~ "FDR_SIGNIFICANT_PROGRESSIVE",
      !is.na(fdr) & fdr < 0.05 & paired_effect_dz < 0 ~ "FDR_SIGNIFICANT_RECOVERY",
      !is.na(fdr) & fdr >= 0.05 & paired_effect_dz > 0 ~ "POSITIVE_DIRECTION_NOT_FDR_SIGNIFICANT",
      !is.na(fdr) & fdr >= 0.05 & paired_effect_dz < 0 ~ "NEGATIVE_DIRECTION_NOT_FDR_SIGNIFICANT",
      !is.na(wilcoxon_p) & wilcoxon_p >= 0.05 & paired_effect_dz > 0 ~ "POSITIVE_DIRECTION_NO_CLEAR_CHANGE",
      !is.na(wilcoxon_p) & wilcoxon_p >= 0.05 & paired_effect_dz < 0 ~ "NEGATIVE_DIRECTION_NO_CLEAR_CHANGE",
      TRUE ~ "UNCERTAIN"
    )
  )

write_csv(
  early_late,
  file.path(OUT, "02_ACTUAL_early_vs_late_Bray_contrasts.csv")
)

# ------------------------------------------------------------
# 4. Preserve the exact Bray robustness registry
# ------------------------------------------------------------

if (!is.null(registry)) {
  write_csv(
    registry,
    file.path(OUT, "03_ACTUAL_Bray_direction_robustness_registry.csv")
  )

  registry_profile <- tibble(
    column = names(registry),
    class = vapply(registry, function(x) class(x)[1], character(1)),
    nonmissing = vapply(registry, function(x) sum(!is.na(x)), integer(1)),
    unique_values = vapply(registry, function(x) length(unique(x[!is.na(x)])), integer(1)),
    example_values = vapply(
      registry,
      function(x) paste(head(unique(as.character(x[!is.na(x)])), 5), collapse = " | "),
      character(1)
    )
  )

  write_csv(
    registry_profile,
    file.path(OUT, "04_Bray_registry_column_profile.csv")
  )
}

# ------------------------------------------------------------
# 5. Preserve exact anchor audit
# ------------------------------------------------------------

if (!is.null(anchor)) {
  write_csv(
    anchor,
    file.path(OUT, "05_ACTUAL_anchor_audit.csv")
  )
}

# ------------------------------------------------------------
# 6. Preserve Bray-relevant global models
# ------------------------------------------------------------

extract_bray_or_all <- function(x) {
  if (is.null(x)) return(NULL)

  metric_col <- find_col(
    names(x),
    c("^metric$", "outcome", "measure", "endpoint")
  )

  if (!is.na(metric_col)) {
    br <- str_detect(
      tolower(as.character(x[[metric_col]])),
      "bray"
    )
    if (any(br, na.rm = TRUE)) {
      return(x[br %in% TRUE, , drop = FALSE])
    }
  }

  x
}

cat_bray <- extract_bray_or_all(catmod)
cont_bray <- extract_bray_or_all(contmod)

if (!is.null(cat_bray)) {
  write_csv(
    cat_bray,
    file.path(OUT, "06_ACTUAL_categorical_global_models_Bray.csv")
  )
}

if (!is.null(cont_bray)) {
  write_csv(
    cont_bray,
    file.path(OUT, "07_ACTUAL_continuous_day_models_Bray.csv")
  )
}

# ------------------------------------------------------------
# 7. Cohort-level primary claim strength based ONLY on actual
#    early -> late paired estimand
# ------------------------------------------------------------

claim_strength <- early_late %>%
  transmute(
    project,
    cohort_role,
    n_pairs,
    paired_effect_dz,
    wilcoxon_p,
    fdr,
    early_late_direction,
    early_late_evidence,
    manuscript_primary_wording = case_when(
      early_late_evidence == "FDR_SIGNIFICANT_PROGRESSIVE" ~
        "Progressive within-patient ecological displacement was supported.",
      early_late_evidence == "FDR_SIGNIFICANT_RECOVERY" ~
        "Recovery/re-convergence toward baseline was supported.",
      early_late_evidence == "NEGATIVE_DIRECTION_NOT_FDR_SIGNIFICANT" ~
        "The cohort showed a recovery-direction trend, but the early-to-late paired contrast was not FDR-significant.",
      early_late_evidence == "POSITIVE_DIRECTION_NOT_FDR_SIGNIFICANT" ~
        "The cohort showed a positive displacement direction, but the early-to-late paired contrast was not FDR-significant.",
      early_late_evidence == "POSITIVE_DIRECTION_NO_CLEAR_CHANGE" ~
        "No clear early-to-late paired change was supported.",
      early_late_evidence == "NEGATIVE_DIRECTION_NO_CLEAR_CHANGE" ~
        "No clear early-to-late paired change was supported; direction was negative.",
      TRUE ~
        "Evidence was uncertain under this estimand."
    )
  )

write_csv(
  claim_strength,
  file.path(OUT, "08_PRIMARY_claim_strength_by_cohort.csv")
)

# ------------------------------------------------------------
# 8. Reconciliation checklist:
#    identify where older/common-anchor claims require sensitivity
#    rather than primary wording
# ------------------------------------------------------------

reconciliation <- claim_strength %>%
  mutate(
    primary_role = "EARLY_VS_LATE_PAIRED_ESTIMAND",
    sensitivity_required = TRUE,
    sensitivity_source = "V2_STEP88B_Bray_direction_and_robustness_registry.csv + anchor audit + global models",
    rule = case_when(
      early_late_evidence == "FDR_SIGNIFICANT_PROGRESSIVE" ~
        "Can state progressive displacement as primary paired result; use registry/global models for robustness.",
      early_late_evidence == "FDR_SIGNIFICANT_RECOVERY" ~
        "Can state recovery as primary paired result; use registry/global models for robustness.",
      str_detect(early_late_evidence, "NEGATIVE_DIRECTION") ~
        "Do NOT state established recovery from early-vs-late alone. Any stronger recovery result must be explicitly labeled sensitivity/common-anchor evidence.",
      str_detect(early_late_evidence, "POSITIVE_DIRECTION") ~
        "Do NOT state significant progression from early-vs-late alone. Any stronger result from another estimand must be explicitly labeled sensitivity evidence.",
      TRUE ~
        "Retain descriptive/supportive wording."
    )
  )

write_csv(
  reconciliation,
  file.path(OUT, "09_ESTIMAND_RECONCILIATION_RULES.csv")
)

# ------------------------------------------------------------
# 9. Figure using only actual early -> late effect estimates
# ------------------------------------------------------------

plot_df <- early_late %>%
  arrange(paired_effect_dz)

make_plot <- function() {
  par(mar = c(5, 14, 3, 2))
  y <- seq_len(nrow(plot_df))

  plot(
    plot_df$paired_effect_dz,
    y,
    xlim = range(c(plot_df$paired_effect_dz, 0), na.rm = TRUE),
    ylim = c(0.5, nrow(plot_df) + 0.5),
    yaxt = "n",
    ylab = "",
    xlab = "Early-to-late paired standardized change (dz)",
    pch = 19,
    main = "Actual Step88B early-to-late Bray-Curtis contrasts"
  )

  abline(v = 0, lty = 2)

  axis(
    2,
    at = y,
    labels = paste0(plot_df$project, "  [", plot_df$cohort_role, "]"),
    las = 1,
    cex.axis = 0.75
  )

  text(
    plot_df$paired_effect_dz,
    y,
    labels = paste0(
      "  n=", plot_df$n_pairs,
      "; FDR=",
      ifelse(
        is.na(plot_df$fdr),
        "NA",
        formatC(plot_df$fdr, digits = 3, format = "f")
      )
    ),
    pos = 4,
    cex = 0.68
  )
}

pdf(
  file.path(OUT, "Figure_STEP93V3_actual_early_late_estimand.pdf"),
  width = 9.5,
  height = 6.5
)
make_plot()
dev.off()

png(
  file.path(OUT, "Figure_STEP93V3_actual_early_late_estimand.png"),
  width = 1900,
  height = 1300,
  res = 180
)
make_plot()
dev.off()

# ------------------------------------------------------------
# 10. Current manuscript guardrails
# ------------------------------------------------------------

control <- claim_strength %>%
  filter(project == "PRJNA578267")

icu_crit <- claim_strength %>%
  filter(
    cohort_role %in% c(
      "CORE_SEPSIS_LONGITUDINAL",
      "ICU_INFECTION_LONGITUDINAL_EXTERNAL",
      "ICU_BACKGROUND_LONGITUDINAL"
    )
  )

lines <- c(
  "STEP93V3 ESTIMAND RECONCILIATION",
  "",
  "Primary paired estimand:",
  "Use the ACTUAL Step88B early-vs-late follow-up contrast as the primary paired cross-cohort effect table.",
  "",
  "Critical-illness cohorts:",
  paste0(
    sum(icu_crit$early_late_evidence == "FDR_SIGNIFICANT_PROGRESSIVE"),
    " of ",
    nrow(icu_crit),
    " critical-illness cohorts are FDR-significant in the progressive-displacement direction under this estimand."
  ),
  "",
  "Non-sepsis control:",
  if (nrow(control) == 1) {
    paste0(
      "PRJNA578267: dz=",
      formatC(control$paired_effect_dz, digits = 3, format = "f"),
      ", p=",
      formatC(control$wilcoxon_p, digits = 4, format = "f"),
      ", FDR=",
      formatC(control$fdr, digits = 3, format = "f"),
      ". Therefore describe this as a negative/re-convergence direction or trend under the primary early-vs-late estimand, NOT established recovery."
    )
  } else {
    "PRJNA578267 not found."
  },
  "",
  "Sensitivity analyses:",
  "Any common-anchor, anchor-robustness-registry, categorical-global-model, or continuous-day result must be reported explicitly as a separate estimand/sensitivity analysis.",
  "",
  "Do not mix n, dz, p, or FDR values across estimands in the same sentence/table row.",
  "",
  "Preferred conservative cross-cohort wording before sensitivity reconciliation:",
  "Progressive ecological displacement was observed across the core sepsis, external ICU-infection, and ICU-background cohorts, whereas the non-sepsis longitudinal control did not show the same progressive pattern and instead showed a negative direction toward baseline."
)

writeLines(
  lines,
  file.path(OUT, "10_STEP93V3_CURRENT_MANUSCRIPT_RULES.txt")
)

writeLines(
  c(
    paste0("Completed: ", Sys.time()),
    "All primary values derived directly from actual Step88B source files.",
    "STEP93V3 COMPLETE"
  ),
  file.path(OUT, "_STEP93V3_COMPLETE.ok")
)

cat("STEP93V3 COMPLETE\n")
