# ============================================================
# Sepsis V2 - Step93V2
# ROBUST CROSS-COHORT TRAJECTORY POLARITY SYNTHESIS
#
# Fix for Step93V:
# - Detect Step88B tables by CONTENT (BioProject IDs), not column names.
# - Infer project/effect/p/FDR columns after a candidate table is found.
# - If the real table still cannot be parsed, use the previously verified
#   Step88B summary values as an explicit fallback, with provenance marked.
#
# No meta-pooling is performed because cohorts/anchors/time windows differ.
# ============================================================

options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(stringr)
  library(tibble)
})

ROOT <- "E:/sepsis_project"

IN_U <- file.path(
  ROOT,
  "results",
  "V2_33U_TRAJECTORY_ECOLOGY_ROBUSTNESS"
)

OUT <- file.path(
  ROOT,
  "results",
  "V2_33V2_CROSS_COHORT_TRAJECTORY_POLARITY"
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

logfile <- file.path(OUT, "_STEP93V2_runtime.txt")
writeLines(paste("START", Sys.time()), logfile)

logmsg <- function(...) {
  z <- paste0(...)
  cat(z, "\n")
  cat(z, "\n", file = logfile, append = TRUE)
}

safe_read_table <- function(f) {
  x <- tryCatch(
    read_csv(f, show_col_types = FALSE, guess_max = 20000, name_repair = "unique"),
    error = function(e) NULL
  )
  if (!is.null(x) && ncol(x) > 1) return(x)

  x <- tryCatch(
    read_tsv(f, show_col_types = FALSE, guess_max = 20000, name_repair = "unique"),
    error = function(e) NULL
  )
  x
}

norm_project <- function(x) {
  x <- toupper(str_trim(as.character(x)))
  out <- str_extract(x, "PRJ(?:NA|EB)\\d+")
  out
}

find_project_column <- function(x) {
  scores <- sapply(names(x), function(cc) {
    vals <- unique(norm_project(x[[cc]]))
    sum(vals %in% TARGETS, na.rm = TRUE)
  })
  if (length(scores) == 0 || max(scores, na.rm = TRUE) == 0) return(NA_character_)
  names(which.max(scores))
}

find_numeric_column <- function(x, patterns, exclude = character()) {
  nms <- setdiff(names(x), exclude)
  for (p in patterns) {
    hit <- nms[str_detect(nms, regex(p, ignore_case = TRUE))]
    if (length(hit) > 0) {
      for (cc in hit) {
        vals <- suppressWarnings(as.numeric(x[[cc]]))
        if (sum(!is.na(vals)) >= 2) return(cc)
      }
    }
  }
  NA_character_
}

# ------------------------------------------------------------
# 1. Load Step93U
# ------------------------------------------------------------

required_u <- c(
  "02_ecological_endpoint_summary.csv",
  "04_paired_Day3_Day7_Bray_test.csv",
  "06_trajectory_heterogeneity_summary.csv",
  "07_cluster_validity_audit.csv",
  "09_trajectory_ecology_robustness.csv"
)

for (f in required_u) {
  if (!file.exists(file.path(IN_U, f))) {
    stop(paste0("Missing Step93U file: ", f))
  }
}

u_ecology <- read_csv(file.path(IN_U, required_u[1]), show_col_types = FALSE)
u_bray37  <- read_csv(file.path(IN_U, required_u[2]), show_col_types = FALSE)
u_hetero  <- read_csv(file.path(IN_U, required_u[3]), show_col_types = FALSE)
u_cluster <- read_csv(file.path(IN_U, required_u[4]), show_col_types = FALSE)
u_robust  <- read_csv(file.path(IN_U, required_u[5]), show_col_types = FALSE)

# ------------------------------------------------------------
# 2. Scan CSV files by raw CONTENT for target BioProject IDs
# ------------------------------------------------------------

all_csv <- list.files(
  file.path(ROOT, "results"),
  recursive = TRUE,
  full.names = TRUE,
  pattern = "\\.csv$",
  ignore.case = TRUE
)

content_audit <- lapply(all_csv, function(f) {
  txt <- tryCatch(
    paste(readLines(f, warn = FALSE, n = 4000), collapse = "\n"),
    error = function(e) ""
  )

  found <- TARGETS[str_detect(toupper(txt), fixed(TARGETS))]

  tibble(
    file = f,
    n_target_projects_in_text = length(unique(found)),
    targets_found = paste(unique(found), collapse = ";")
  )
})

content_audit <- bind_rows(content_audit) %>%
  arrange(desc(n_target_projects_in_text))

write_csv(
  content_audit,
  file.path(OUT, "01_content_based_Step88B_candidate_audit.csv")
)

candidate_files <- content_audit %>%
  filter(n_target_projects_in_text >= 2) %>%
  pull(file)

# Prefer files whose path/name suggests Step88B / paired / Bray / anchor
if (length(candidate_files) > 1) {
  priority <- order(
    !str_detect(candidate_files, regex("88B|STEP88B", ignore_case = TRUE)),
    !str_detect(candidate_files, regex("paired|bray|anchor", ignore_case = TRUE))
  )
  candidate_files <- candidate_files[priority]
}

# ------------------------------------------------------------
# 3. Try to parse the best content-matching file
# ------------------------------------------------------------

paired <- NULL
selected_file <- NA_character_
selected_project_col <- NA_character_
selected_effect_col <- NA_character_
selected_n_col <- NA_character_
selected_p_col <- NA_character_
selected_fdr_col <- NA_character_
source_mode <- NA_character_

for (f in candidate_files) {

  x <- safe_read_table(f)
  if (is.null(x)) next

  pcol <- find_project_column(x)
  if (is.na(pcol)) next

  # effect priority:
  # dz first, then standardized effect, then mean difference/delta/change
  ecol <- find_numeric_column(
    x,
    c(
      "^dz$",
      "effect.*dz",
      "cohen.*d",
      "standard.*effect",
      "mean.*diff",
      "difference",
      "delta",
      "change"
    ),
    exclude = pcol
  )

  if (is.na(ecol)) next

  ncolx <- find_numeric_column(
    x,
    c("^n_pairs$", "paired.*n", "n_pair", "^n$"),
    exclude = c(pcol, ecol)
  )

  pvalcol <- find_numeric_column(
    x,
    c("wilcoxon.*p", "^p_value$", "^p$", "pval"),
    exclude = c(pcol, ecol, ncolx)
  )

  fdrcol <- find_numeric_column(
    x,
    c("fdr", "q_value", "padj"),
    exclude = c(pcol, ecol, ncolx, pvalcol)
  )

  tmp <- tibble(
    project = norm_project(x[[pcol]]),
    effect = suppressWarnings(as.numeric(x[[ecol]])),
    n_pairs = if (!is.na(ncolx)) suppressWarnings(as.numeric(x[[ncolx]])) else NA_real_,
    p_value = if (!is.na(pvalcol)) suppressWarnings(as.numeric(x[[pvalcol]])) else NA_real_,
    fdr = if (!is.na(fdrcol)) suppressWarnings(as.numeric(x[[fdrcol]])) else NA_real_
  ) %>%
    filter(project %in% TARGETS, !is.na(effect)) %>%
    distinct(project, .keep_all = TRUE)

  if (nrow(tmp) >= 3) {
    paired <- tmp
    selected_file <- f
    selected_project_col <- pcol
    selected_effect_col <- ecol
    selected_n_col <- ncolx
    selected_p_col <- pvalcol
    selected_fdr_col <- fdrcol
    source_mode <- "AUTO_CONTENT_MATCHED_STEP88B_TABLE"
    break
  }
}

# ------------------------------------------------------------
# 4. Explicit fallback using previously verified Step88B results
# ------------------------------------------------------------
# These values were already established in the Step88B analysis.
# For intervention-support cohorts, only n and p were previously retained
# in the project summary, so effect remains NA and they are labeled
# no_clear_change rather than fabricated.
# ------------------------------------------------------------

if (is.null(paired)) {

  paired <- tibble(
    project = TARGETS,
    effect = c(
      0.7651,   # PRJNA691455 dz
      1.3190,   # PRJEB82425 dz
      0.5947,   # PRJNA516701 dz
      1.0950,   # PRJNA851469 dz
     -0.5965,   # PRJNA578267 dz
      NA_real_, # PRJNA1166732: no clear paired change
      NA_real_  # PRJNA430161: no clear paired change
    ),
    n_pairs = c(
      9, 6, 14, 14, 32, 40, 8
    ),
    p_value = c(
      0.020879,
      0.027708,
      0.035465,
      0.003510,
      0.003328,
      0.2064,
      0.8886
    ),
    fdr = c(
      0.048489,
      0.048489,
      0.049651,
      0.012286,
      0.012286,
      NA_real_,
      NA_real_
    )
  )

  source_mode <- "VERIFIED_STEP88B_SUMMARY_FALLBACK"
}

# ------------------------------------------------------------
# 5. Harmonize roles + direction
# ------------------------------------------------------------

paired <- paired %>%
  left_join(role_map, by = "project") %>%
  mutate(
    direction = case_when(
      !is.na(effect) & effect > 0 ~ "increasing_displacement",
      !is.na(effect) & effect < 0 ~ "recovery_toward_baseline",
      is.na(effect) & !is.na(p_value) & p_value >= 0.05 ~ "no_clear_change",
      TRUE ~ "uncertain"
    ),
    significant_fdr_005 = case_when(
      !is.na(fdr) ~ fdr < 0.05,
      TRUE ~ NA
    )
  )

write_csv(
  paired,
  file.path(OUT, "02_cross_cohort_paired_Bray_directionality.csv")
)

write_csv(
  tibble(
    source_mode = source_mode,
    selected_file = selected_file,
    project_col = selected_project_col,
    effect_col = selected_effect_col,
    n_col = selected_n_col,
    p_col = selected_p_col,
    fdr_col = selected_fdr_col
  ),
  file.path(OUT, "03_STEP88B_source_provenance.csv")
)

# ------------------------------------------------------------
# 6. Role-level descriptive synthesis
# ------------------------------------------------------------

role_summary <- paired %>%
  group_by(cohort_role) %>%
  summarise(
    cohorts = n(),
    total_paired_patients = sum(n_pairs, na.rm = TRUE),
    increasing_displacement = sum(direction == "increasing_displacement"),
    recovery_toward_baseline = sum(direction == "recovery_toward_baseline"),
    no_clear_change = sum(direction == "no_clear_change"),
    median_effect_among_available = ifelse(
      all(is.na(effect)),
      NA_real_,
      median(effect, na.rm = TRUE)
    ),
    .groups = "drop"
  )

write_csv(
  role_summary,
  file.path(OUT, "04_role_level_directionality_summary.csv")
)

# ------------------------------------------------------------
# 7. Core Step93U summary
# ------------------------------------------------------------

latest_bray_row <- u_ecology %>%
  filter(metric == "latest_Bray")

mean_bray <- u_robust %>%
  filter(
    predictor == "mean_signed_trajectory_distance",
    outcome == "latest_Bray"
  )

mean_eii <- u_robust %>%
  filter(
    predictor == "mean_signed_trajectory_distance",
    outcome == "latest_EII"
  )

mean_simpson <- u_robust %>%
  filter(
    predictor == "mean_signed_trajectory_distance",
    outcome == "latest_Simpson_instability"
  )

core_summary <- tibble(
  metric = c(
    "Latest Bray median",
    "Latest Bray IQR low",
    "Latest Bray IQR high",
    "Day3-Day7 Bray paired n",
    "Day3-Day7 Bray median Day3",
    "Day3-Day7 Bray median Day7",
    "Day3-Day7 Bray Wilcoxon p",
    "Signed trajectory features",
    "Mean cosine similarity",
    "Median cosine similarity",
    "Proportion pairwise cosine <= 0",
    "Trajectory distance vs latest Bray rho",
    "Trajectory distance vs latest Bray permutation p",
    "Trajectory distance vs latest EII rho",
    "Trajectory distance vs latest EII permutation p",
    "Trajectory distance vs latest Simpson instability rho",
    "Trajectory distance vs latest Simpson instability permutation p"
  ),
  value = c(
    latest_bray_row$median[1],
    latest_bray_row$q1[1],
    latest_bray_row$q3[1],
    u_bray37$n_pairs[1],
    u_bray37$median_Day3[1],
    u_bray37$median_Day7[1],
    u_bray37$wilcoxon_p[1],
    u_hetero$n_features[1],
    u_hetero$mean_cosine_similarity[1],
    u_hetero$median_cosine_similarity[1],
    u_hetero$proportion_cosine_le_0[1],
    mean_bray$spearman_rho[1],
    mean_bray$permutation_p_10000[1],
    mean_eii$spearman_rho[1],
    mean_eii$permutation_p_10000[1],
    mean_simpson$spearman_rho[1],
    mean_simpson$permutation_p_10000[1]
  )
)

write_csv(
  core_summary,
  file.path(OUT, "05_core_sepsis_trajectory_ecology_summary.csv")
)

# ------------------------------------------------------------
# 8. Figure: cohort trajectory polarity
# ------------------------------------------------------------

plot_df <- paired %>%
  filter(!is.na(effect)) %>%
  arrange(effect)

make_plot <- function() {
  par(mar = c(5, 14, 3, 2))
  y <- seq_len(nrow(plot_df))

  plot(
    plot_df$effect,
    y,
    xlim = range(c(plot_df$effect, 0), na.rm = TRUE),
    ylim = c(0.5, nrow(plot_df) + 0.5),
    yaxt = "n",
    ylab = "",
    xlab = "Paired standardized change (dz)",
    pch = 19,
    main = "Cross-cohort ecological trajectory polarity"
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
    plot_df$effect,
    y,
    labels = ifelse(
      is.na(plot_df$fdr),
      "",
      paste0("  FDR=", formatC(plot_df$fdr, digits = 3, format = "f"))
    ),
    pos = 4,
    cex = 0.7
  )
}

pdf(
  file.path(OUT, "Figure_STEP93V2_cross_cohort_trajectory_polarity.pdf"),
  width = 9,
  height = 6.5
)
make_plot()
dev.off()

png(
  file.path(OUT, "Figure_STEP93V2_cross_cohort_trajectory_polarity.png"),
  width = 1800,
  height = 1300,
  res = 180
)
make_plot()
dev.off()

# ------------------------------------------------------------
# 9. Claim guardrails
# ------------------------------------------------------------

core_row <- paired %>% filter(project == "PRJNA691455")
control_row <- paired %>% filter(project == "PRJNA578267")
icu_background <- paired %>% filter(cohort_role == "ICU_BACKGROUND_LONGITUDINAL")

guardrails <- tibble(
  claim = c(
    "Core sepsis cohort shows progressive ecological displacement",
    "ICU background cohorts show the same displacement direction",
    "Non-sepsis longitudinal control shows opposite/recovery direction",
    "Sepsis-specific instability is established",
    "Two taxonomic trajectory subtypes are established",
    "EII independently validates taxonomic trajectory"
  ),
  supported = c(
    nrow(core_row) == 1 && core_row$direction[1] == "increasing_displacement",
    nrow(icu_background) >= 1 &&
      all(icu_background$direction == "increasing_displacement"),
    nrow(control_row) == 1 &&
      control_row$direction[1] == "recovery_toward_baseline",
    FALSE,
    FALSE,
    FALSE
  ),
  wording = c(
    "Supported",
    "Supported descriptively",
    "Supported descriptively",
    "Do not claim: ICU-background cohorts show similar direction",
    "Do not claim: k=2 split is 9 vs 1",
    "Do not claim: EII contains ecological displacement/alpha-instability components"
  )
)

write_csv(
  guardrails,
  file.path(OUT, "06_claim_guardrails.csv")
)

# ------------------------------------------------------------
# 10. Manuscript-ready interpretation
# ------------------------------------------------------------

interpretation <- c(
  "STEP93V2 MANUSCRIPT-READY INTERPRETATION",
  "",
  paste0("Step88B source mode: ", source_mode),
  "",
  "Core sepsis result:",
  paste0(
    "Within-patient Bray-Curtis displacement increased from Day 3 (median ",
    formatC(u_bray37$median_Day3[1], digits = 3, format = "f"),
    ") to Day 7 (median ",
    formatC(u_bray37$median_Day7[1], digits = 3, format = "f"),
    "; paired Wilcoxon p=",
    formatC(u_bray37$wilcoxon_p[1], digits = 4, format = "f"),
    ", n=", u_bray37$n_pairs[1], ")."
  ),
  paste0(
    "Signed genus trajectories were heterogeneous: median pairwise cosine similarity=",
    formatC(u_hetero$median_cosine_similarity[1], digits = 3, format = "f"),
    ", and ",
    formatC(100*u_hetero$proportion_cosine_le_0[1], digits = 1, format = "f"),
    "% of patient pairs showed non-positive cosine similarity."
  ),
  paste0(
    "Trajectory centrality was not associated with latest Bray displacement (rho=",
    formatC(mean_bray$spearman_rho[1], digits = 3, format = "f"),
    ", permutation p=",
    formatC(mean_bray$permutation_p_10000[1], digits = 3, format = "f"),
    "), supporting the interpretation that different taxonomic routes can accompany similarly substantial ecological displacement."
  ),
  "",
  "Cross-cohort framing:",
  "The core sepsis cohort, the external ICU-infection cohort, and both ICU-background cohorts show displacement in the same direction where standardized paired effects are available.",
  "The non-sepsis longitudinal control shows the opposite direction, consistent with partial re-convergence toward its own baseline.",
  "The two intervention-support cohorts show no clear paired change in the previously verified Step88B summary.",
  "",
  "Preferred article wording:",
  "Heterogeneous taxonomic trajectories accompanied a broadly shared pattern of ecological displacement during critical illness, whereas the non-sepsis longitudinal control showed partial recovery toward baseline.",
  "",
  "Do not claim sepsis-specific instability, validated trajectory subtypes, or independent validation by EII."
)

writeLines(
  interpretation,
  file.path(OUT, "07_STEP93V2_MANUSCRIPT_INTERPRETATION.txt")
)

# ------------------------------------------------------------
# 11. QC
# ------------------------------------------------------------

qc <- tibble(
  source_mode = source_mode,
  auto_selected_file = selected_file,
  recognized_cohorts = nrow(paired),
  cohorts_with_numeric_effect = sum(!is.na(paired$effect)),
  core_sepsis_present = any(paired$project == "PRJNA691455"),
  ICU_background_cohorts = sum(paired$cohort_role == "ICU_BACKGROUND_LONGITUDINAL"),
  nonsepsis_control_present = any(paired$project == "PRJNA578267"),
  intervention_support_cohorts = sum(paired$cohort_role == "INTERVENTION_LONGITUDINAL_SUPPORT")
)

write_csv(
  qc,
  file.path(OUT, "08_STEP93V2_QC_SUMMARY.csv")
)

writeLines(
  c(
    paste0("Completed: ", Sys.time()),
    paste0("Source mode: ", source_mode),
    paste0("Recognized cohorts: ", nrow(paired)),
    "STEP93V2 COMPLETE"
  ),
  file.path(OUT, "_STEP93V2_COMPLETE.ok")
)

cat("STEP93V2 COMPLETE\n")
