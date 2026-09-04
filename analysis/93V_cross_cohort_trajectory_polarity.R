# ============================================================
# Sepsis V2 - Step93V
# CROSS-COHORT ECOLOGICAL TRAJECTORY POLARITY SYNTHESIS
#
# Goal
# ----
# Integrate the already-completed Step88B paired Bray analyses with the
# Step93U taxonomic-trajectory findings.
#
# This step deliberately DOES NOT meta-pool all cohorts, because anchors,
# populations, and time windows are not fully exchangeable.
#
# Main outputs:
# 1) locate the real Step88B paired-Bray table automatically
# 2) harmonize cohort roles
# 3) summarize direction of ecological displacement/recovery
# 4) generate a forest-style directionality figure
# 5) combine the cross-cohort result with Step93U core-sepsis
#    taxonomic heterogeneity results
# 6) write manuscript-ready guarded interpretation
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
  "V2_33V_CROSS_COHORT_TRAJECTORY_POLARITY"
)

dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

logfile <- file.path(OUT, "_STEP93V_runtime.txt")
writeLines(paste("START", Sys.time()), logfile)

logmsg <- function(...) {
  x <- paste0(...)
  cat(x, "\n")
  cat(x, "\n", file = logfile, append = TRUE)
}

safe_read <- function(f) {
  tryCatch(
    read_csv(f, show_col_types = FALSE, guess_max = 10000),
    error = function(e) NULL
  )
}

first_match <- function(nms, patterns) {
  for (p in patterns) {
    hit <- nms[str_detect(nms, regex(p, ignore_case = TRUE))]
    if (length(hit) > 0) return(hit[1])
  }
  NA_character_
}

# ------------------------------------------------------------
# 1. Load Step93U key results
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

u_ecology <- read_csv(
  file.path(IN_U, "02_ecological_endpoint_summary.csv"),
  show_col_types = FALSE
)

u_bray37 <- read_csv(
  file.path(IN_U, "04_paired_Day3_Day7_Bray_test.csv"),
  show_col_types = FALSE
)

u_hetero <- read_csv(
  file.path(IN_U, "06_trajectory_heterogeneity_summary.csv"),
  show_col_types = FALSE
)

u_cluster <- read_csv(
  file.path(IN_U, "07_cluster_validity_audit.csv"),
  show_col_types = FALSE
)

u_robust <- read_csv(
  file.path(IN_U, "09_trajectory_ecology_robustness.csv"),
  show_col_types = FALSE
)

# ------------------------------------------------------------
# 2. Find Step88B candidate CSVs
# ------------------------------------------------------------

all_csv <- list.files(
  file.path(ROOT, "results"),
  recursive = TRUE,
  full.names = TRUE,
  pattern = "\\.csv$",
  ignore.case = TRUE
)

cand_paths <- all_csv[
  str_detect(
    all_csv,
    regex("88B|STEP88B", ignore_case = TRUE)
  )
]

if (length(cand_paths) == 0) {
  # broad fallback search: filenames likely to hold paired/common-anchor Bray
  cand_paths <- all_csv[
    str_detect(
      basename(all_csv),
      regex("paired|common.*anchor|anchor.*paired|bray", ignore_case = TRUE)
    )
  ]
}

audit <- list()
loaded <- list()

for (f in cand_paths) {

  x <- safe_read(f)
  if (is.null(x)) next

  nms <- names(x)

  project_col <- first_match(
    nms,
    c("^project$", "bioproject", "^cohort$", "project")
  )

  n_col <- first_match(
    nms,
    c("^n_pairs$", "paired_n", "n_pair", "^n$")
  )

  dz_col <- first_match(
    nms,
    c("^dz$", "effect.*dz", "paired.*dz", "cohen.*d")
  )

  p_col <- first_match(
    nms,
    c("wilcoxon.*p", "^p_value$", "^p$", "pval")
  )

  fdr_col <- first_match(
    nms,
    c("fdr", "q_value", "padj")
  )

  metric_col <- first_match(
    nms,
    c("^metric$", "outcome", "measure", "endpoint")
  )

  score <-
    5 * !is.na(project_col) +
    5 * !is.na(dz_col) +
    2 * !is.na(n_col) +
    2 * !is.na(p_col) +
    1 * !is.na(fdr_col) +
    1 * !is.na(metric_col)

  audit[[length(audit) + 1]] <- tibble(
    file = f,
    rows = nrow(x),
    cols = ncol(x),
    score = score,
    project_col = project_col,
    n_col = n_col,
    dz_col = dz_col,
    p_col = p_col,
    fdr_col = fdr_col,
    metric_col = metric_col,
    columns = paste(nms, collapse = ";")
  )

  loaded[[f]] <- x
}

audit <- bind_rows(audit) %>%
  arrange(desc(score), desc(rows))

write_csv(
  audit,
  file.path(OUT, "01_STEP88B_candidate_table_audit.csv")
)

if (nrow(audit) == 0) {
  stop("No readable Step88B candidate CSV found.")
}

best <- audit[1, ]

if (best$score[[1]] < 10) {
  writeLines(
    c(
      "No high-confidence Step88B paired-effect table was found.",
      "Inspect 01_STEP88B_candidate_table_audit.csv."
    ),
    file.path(OUT, "_STEP93V_NEEDS_STEP88B_TABLE_REVIEW.txt")
  )
  stop("STEP93V could not confidently identify the Step88B paired-effect table.")
}

selected_file <- best$file[[1]]
x <- loaded[[selected_file]]

write_csv(
  best,
  file.path(OUT, "02_SELECTED_STEP88B_TABLE.csv")
)

logmsg("Selected Step88B table: ", selected_file)

# ------------------------------------------------------------
# 3. Normalize the paired Bray table
# ------------------------------------------------------------

project_col <- best$project_col[[1]]
n_col <- best$n_col[[1]]
dz_col <- best$dz_col[[1]]
p_col <- best$p_col[[1]]
fdr_col <- best$fdr_col[[1]]
metric_col <- best$metric_col[[1]]

if (!is.na(metric_col)) {
  metric_text <- tolower(as.character(x[[metric_col]]))
  keep_bray <- str_detect(metric_text, "bray")

  # only filter if Bray rows actually exist
  if (any(keep_bray, na.rm = TRUE)) {
    x <- x[keep_bray %in% TRUE, , drop = FALSE]
  }
}

paired <- tibble(
  project = as.character(x[[project_col]]),
  n_pairs = if (!is.na(n_col)) suppressWarnings(as.numeric(x[[n_col]])) else NA_real_,
  dz = suppressWarnings(as.numeric(x[[dz_col]])),
  p_value = if (!is.na(p_col)) suppressWarnings(as.numeric(x[[p_col]])) else NA_real_,
  fdr = if (!is.na(fdr_col)) suppressWarnings(as.numeric(x[[fdr_col]])) else NA_real_
)

# Keep recognized cohorts only
role_map <- tibble(
  project = c(
    "PRJNA691455",
    "PRJEB82425",
    "PRJNA516701",
    "PRJNA851469",
    "PRJNA578267",
    "PRJNA1166732",
    "PRJNA430161"
  ),
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

paired <- paired %>%
  mutate(
    project = toupper(str_trim(project))
  ) %>%
  inner_join(role_map, by = "project") %>%
  filter(!is.na(dz)) %>%
  distinct(project, .keep_all = TRUE) %>%
  mutate(
    direction = case_when(
      dz > 0 ~ "increasing_displacement",
      dz < 0 ~ "recovery_toward_baseline",
      TRUE ~ "neutral"
    ),
    significant_fdr_005 = ifelse(!is.na(fdr), fdr < 0.05, NA)
  )

write_csv(
  paired,
  file.path(OUT, "03_cross_cohort_paired_Bray_directionality.csv")
)

# ------------------------------------------------------------
# 4. Role-level synthesis (descriptive only; no meta-pooling)
# ------------------------------------------------------------

role_summary <- paired %>%
  group_by(cohort_role) %>%
  summarise(
    cohorts = n(),
    total_paired_patients = sum(n_pairs, na.rm = TRUE),
    positive_dz = sum(dz > 0, na.rm = TRUE),
    negative_dz = sum(dz < 0, na.rm = TRUE),
    median_dz = median(dz, na.rm = TRUE),
    min_dz = min(dz, na.rm = TRUE),
    max_dz = max(dz, na.rm = TRUE),
    .groups = "drop"
  )

write_csv(
  role_summary,
  file.path(OUT, "04_role_level_directionality_summary.csv")
)

# ------------------------------------------------------------
# 5. Extract key Step93U results into one table
# ------------------------------------------------------------

latest_bray_row <- u_ecology %>%
  filter(metric == "latest_Bray")

mean_eii <- u_robust %>%
  filter(
    predictor == "mean_signed_trajectory_distance",
    outcome == "latest_EII"
  )

mean_bray <- u_robust %>%
  filter(
    predictor == "mean_signed_trajectory_distance",
    outcome == "latest_Bray"
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
# 6. Cross-cohort figure
# ------------------------------------------------------------

plot_df <- paired %>%
  arrange(dz)

pdf(
  file.path(OUT, "Figure_STEP93V_cross_cohort_trajectory_polarity.pdf"),
  width = 9,
  height = 6.5
)

par(mar = c(5, 14, 3, 2))

y <- seq_len(nrow(plot_df))

plot(
  plot_df$dz,
  y,
  xlim = range(c(plot_df$dz, 0), na.rm = TRUE),
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
  plot_df$dz,
  y,
  labels = ifelse(
    is.na(plot_df$fdr),
    "",
    paste0("  FDR=", formatC(plot_df$fdr, digits = 3, format = "f"))
  ),
  pos = 4,
  cex = 0.7
)

dev.off()

png(
  file.path(OUT, "Figure_STEP93V_cross_cohort_trajectory_polarity.png"),
  width = 1800,
  height = 1300,
  res = 180
)

par(mar = c(5, 14, 3, 2))

plot(
  plot_df$dz,
  y,
  xlim = range(c(plot_df$dz, 0), na.rm = TRUE),
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
  plot_df$dz,
  y,
  labels = ifelse(
    is.na(plot_df$fdr),
    "",
    paste0("  FDR=", formatC(plot_df$fdr, digits = 3, format = "f"))
  ),
  pos = 4,
  cex = 0.7
)

dev.off()

# ------------------------------------------------------------
# 7. Guardrails + manuscript-ready interpretation
# ------------------------------------------------------------

core_row <- paired %>%
  filter(project == "PRJNA691455")

control_row <- paired %>%
  filter(project == "PRJNA578267")

icu_background <- paired %>%
  filter(
    cohort_role == "ICU_BACKGROUND_LONGITUDINAL"
  )

external_icu <- paired %>%
  filter(
    cohort_role == "ICU_INFECTION_LONGITUDINAL_EXTERNAL"
  )

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
    nrow(core_row) == 1 && core_row$dz[1] > 0,
    nrow(icu_background) >= 1 && all(icu_background$dz > 0),
    nrow(control_row) == 1 && control_row$dz[1] < 0,
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

interpretation <- c(
  "STEP93V MANUSCRIPT-READY INTERPRETATION",
  "",
  "Core result:",
  paste0(
    "In PRJNA691455, within-patient Bray-Curtis displacement increased from Day 3 (median ",
    formatC(u_bray37$median_Day3[1], digits = 3, format = "f"),
    ") to Day 7 (median ",
    formatC(u_bray37$median_Day7[1], digits = 3, format = "f"),
    "; paired Wilcoxon p=",
    formatC(u_bray37$wilcoxon_p[1], digits = 4, format = "f"),
    ", n=",
    u_bray37$n_pairs[1],
    ")."
  ),
  paste0(
    "At the genus level, signed trajectories were heterogeneous across patients: median pairwise cosine similarity=",
    formatC(u_hetero$median_cosine_similarity[1], digits = 3, format = "f"),
    ", with ",
    formatC(100*u_hetero$proportion_cosine_le_0[1], digits = 1, format = "f"),
    "% of patient pairs showing non-positive cosine similarity."
  ),
  paste0(
    "Trajectory centrality was not associated with the magnitude of latest Bray displacement (rho=",
    formatC(mean_bray$spearman_rho[1], digits = 3, format = "f"),
    ", permutation p=",
    formatC(mean_bray$permutation_p_10000[1], digits = 3, format = "f"),
    "), indicating that markedly different taxonomic routes can accompany similarly large ecological displacement."
  ),
  "",
  "Cross-cohort framing:",
  "Do not describe this as sepsis-specific instability.",
  "ICU-background cohorts also show increasing displacement, whereas the non-sepsis longitudinal control shows the opposite direction toward its own baseline.",
  "The preferred framing is therefore critical-illness-associated ecological destabilization versus recovery, with the sepsis cohort providing the deepest taxonomic trajectory analysis.",
  "",
  "EII interpretation:",
  "EII associations are supportive internal ecological consistency analyses only; EII is not an independent prognostic or external validation endpoint.",
  "",
  "Cluster interpretation:",
  "Do not use the k=2 9-vs-1 split as evidence for taxonomic trajectory subtypes."
)

writeLines(
  interpretation,
  file.path(OUT, "07_STEP93V_MANUSCRIPT_INTERPRETATION.txt")
)

# ------------------------------------------------------------
# 8. QC
# ------------------------------------------------------------

qc <- tibble(
  selected_step88B_file = selected_file,
  recognized_cohorts = nrow(paired),
  core_sepsis_present = any(paired$project == "PRJNA691455"),
  ICU_background_cohorts = sum(paired$cohort_role == "ICU_BACKGROUND_LONGITUDINAL"),
  nonsepsis_control_present = any(paired$project == "PRJNA578267"),
  intervention_support_cohorts = sum(paired$cohort_role == "INTERVENTION_LONGITUDINAL_SUPPORT")
)

write_csv(
  qc,
  file.path(OUT, "08_STEP93V_QC_SUMMARY.csv")
)

writeLines(
  c(
    paste0("Completed: ", Sys.time()),
    paste0("Selected Step88B table: ", selected_file),
    paste0("Recognized cohorts: ", nrow(paired)),
    "STEP93V COMPLETE"
  ),
  file.path(OUT, "_STEP93V_COMPLETE.ok")
)

cat("STEP93V COMPLETE\n")
