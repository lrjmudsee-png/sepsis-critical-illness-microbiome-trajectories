
# ============================================================
# V2_96B RESULTS INTEGRATION
#
# Purpose:
# Integrate the frozen Step36J cross-cohort taxonomic
# reproducibility result into the frozen Results v5 manuscript.
#
# NO statistical model is rerun.
# NO frozen numeric result is changed.
# ============================================================

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(stringr)
})

ROOT <- "E:/sepsis_project"
RESULTS <- file.path(ROOT, "results")

RESULTS_V5 <- file.path(
  RESULTS,
  "V2_35F3_RESULTS_AND_TABLE1_FINAL_FREEZE_SENTENCE_AUDIT",
  "02_FINAL_RESULTS",
  "Results_Draft_v5_FINAL_MANUSCRIPT_FREEZE.txt"
)

FREEZE_DIR <- file.path(
  RESULTS,
  "V2_36J_MICROBIOME_ENHANCEMENT_EVIDENCE_FREEZE"
)

METRICS_FILE <- file.path(
  FREEZE_DIR,
  "02_FROZEN_ENHANCEMENT_METRICS.csv"
)

BRANCH_FILE <- file.path(
  FREEZE_DIR,
  "01_ENHANCEMENT_BRANCH_FREEZE.csv"
)

OUT <- file.path(
  RESULTS,
  "V2_96B_RESULTS_INTEGRATION"
)

dir.create(OUT, recursive=TRUE, showWarnings=FALSE)

if (!file.exists(RESULTS_V5)) {
  stop("Frozen Results v5 not found: ", RESULTS_V5)
}

if (!file.exists(METRICS_FILE)) {
  stop("Step36J frozen metrics not found: ", METRICS_FILE)
}

if (!file.exists(BRANCH_FILE)) {
  stop("Step36J branch freeze not found: ", BRANCH_FILE)
}

# ------------------------------------------------------------
# 1. Validate frozen branch
# ------------------------------------------------------------
branch <- read_csv(BRANCH_FILE, show_col_types=FALSE)

tax_branch <- branch %>%
  filter(branch == "CROSS_COHORT_TAXONOMIC_REPRODUCIBILITY")

if (nrow(tax_branch) != 1) {
  stop("Could not uniquely identify frozen taxonomic reproducibility branch.")
}

if (
  tax_branch$evidence_tier[[1]] !=
  "ECOLOGICAL_DIRECTION_MORE_REPRODUCIBLE_THAN_TAXONOMIC_EFFECTS"
) {
  stop("Unexpected Step36J taxonomic evidence tier.")
}

# ------------------------------------------------------------
# 2. Read exact frozen metrics
# ------------------------------------------------------------
m <- read_csv(METRICS_FILE, show_col_types=FALSE)

get_metric <- function(name) {
  z <- m %>% filter(metric == name)
  if (nrow(z) != 1) {
    stop("Metric missing or duplicated: ", name)
  }
  as.numeric(z$value[[1]])
}

eco <- get_metric("primary_ecological_direction_concordance")

g_rho <- get_metric("genus_median_pairwise_spearman_rho")
g_sign <- get_metric("genus_median_sign_agreement")
g_jac <- get_metric("genus_median_top_jaccard")
g_three <- get_metric("genus_threeway_same_direction_proportion")

f_rho <- get_metric("family_median_pairwise_spearman_rho")
f_sign <- get_metric("family_median_sign_agreement")
f_jac <- get_metric("family_median_top_jaccard")
f_three <- get_metric("family_threeway_same_direction_proportion")

source_metrics <- tibble(
  metric = c(
    "primary_ecological_direction_concordance",
    "genus_median_pairwise_spearman_rho",
    "genus_median_sign_agreement",
    "genus_median_top_jaccard",
    "genus_threeway_same_direction_proportion",
    "family_median_pairwise_spearman_rho",
    "family_median_sign_agreement",
    "family_median_top_jaccard",
    "family_threeway_same_direction_proportion"
  ),
  value = c(
    eco, g_rho, g_sign, g_jac, g_three,
    f_rho, f_sign, f_jac, f_three
  )
)

write_csv(
  source_metrics,
  file.path(OUT, "02_SOURCE_METRICS_USED.csv")
)

# ------------------------------------------------------------
# 3. Build manuscript paragraph from frozen metrics
# ------------------------------------------------------------
new_heading <- "Cross-cohort ecological concordance accompanied taxonomic heterogeneity"

new_paragraph <- paste0(
  "Cross-cohort taxonomic effects were compared after estimating longitudinal paired CLR changes independently within each cohort, without pooling ASV or OTU abundance matrices across studies. ",
  "Among the three cohorts with robust progressive ecological displacement (PRJNA691455, PRJNA851469, and PRJNA516701), the direction of ecosystem-level change was concordant across all three cohorts, whereas the corresponding taxonomic effects were substantially less reproducible. ",
  "At the genus level, the median pairwise Spearman correlation between within-cohort CLR effect vectors was ",
  sprintf("%.2f", g_rho),
  ", with a median directional agreement of ",
  sprintf("%.1f", 100 * g_sign),
  "% and a median Jaccard overlap of ",
  sprintf("%.1f", 100 * g_jac),
  "% among the strongest effects. ",
  "At the family level, the median pairwise effect-vector correlation was ",
  sprintf("%.2f", f_rho),
  ", with a median directional agreement of ",
  sprintf("%.1f", 100 * f_sign),
  "%. ",
  "Among taxa shared across all three cohorts, only ",
  sprintf("%.1f", 100 * g_three),
  "% of genera and ",
  sprintf("%.1f", 100 * f_three),
  "% of families changed in the same direction in all three cohorts. ",
  "The common-anchor sensitivity analysis showed the same overall pattern. ",
  "These findings indicate that concordant ecosystem-level displacement across independent cohorts was accompanied by only low-to-moderate concordance of longitudinal taxonomic effects, supporting heterogeneous and cohort-dependent taxonomic routes to a shared direction of ecological disturbance."
)

# ------------------------------------------------------------
# 4. Read Results v5 and locate insertion point
# ------------------------------------------------------------
lines <- readLines(
  RESULTS_V5,
  warn=FALSE,
  encoding="UTF-8"
)

# Prevent accidental duplicate insertion
if (any(str_detect(
  tolower(lines),
  "cross-cohort ecological concordance accompanied taxonomic heterogeneity"
))) {
  stop("Results v5 already appears to contain the Step36J insertion.")
}

# Preferred insertion:
# after Section 3 and before the external validation section.
external_candidates <- which(
  str_detect(
    tolower(lines),
    "external ecological-state differentiation beyond severe trauma"
  )
)

if (length(external_candidates) == 0) {
  external_candidates <- which(
    str_detect(
      tolower(lines),
      "external.*trauma"
    )
  )
}

if (length(external_candidates) != 1) {

  heading_candidates <- tibble(
    line_number = seq_along(lines),
    text = lines
  ) %>%
    filter(
      str_detect(
        tolower(text),
        "trajectory|external|trauma|infection source|ecological"
      )
    )

  write_csv(
    heading_candidates,
    file.path(OUT, "00_HEADING_CANDIDATES_FOR_MANUAL_REVIEW.csv")
  )

  stop(
    "Could not uniquely locate the external-validation Results section. ",
    "See 00_HEADING_CANDIDATES_FOR_MANUAL_REVIEW.csv"
  )
}

insert_before <- external_candidates[[1]]

before <- if (insert_before > 1) {
  lines[1:(insert_before - 1)]
} else {
  character()
}

after <- lines[insert_before:length(lines)]

insert_block <- c(
  "",
  new_heading,
  "",
  new_paragraph,
  ""
)

new_lines <- c(
  before,
  insert_block,
  after
)

OUT_RESULTS <- file.path(
  OUT,
  "Results_Draft_v6_WITH_CROSS_COHORT_REPRODUCIBILITY.txt"
)

writeLines(
  new_lines,
  OUT_RESULTS,
  useBytes=TRUE
)

# ------------------------------------------------------------
# 5. Audit
# ------------------------------------------------------------
audit <- tibble(
  source_results = RESULTS_V5,
  frozen_metrics_source = METRICS_FILE,
  insertion_before_original_line = insert_before,
  inserted_heading = new_heading,
  original_line_count = length(lines),
  new_line_count = length(new_lines),
  numeric_inference_rerun = FALSE,
  frozen_results_modified = FALSE,
  manuscript_copy_created = TRUE
)

write_csv(
  audit,
  file.path(OUT, "01_INSERTION_AUDIT.csv")
)

writeLines(
  c(
    "V2 STEP96B RESULTS INTEGRATION",
    "",
    paste0("Source: ", RESULTS_V5),
    paste0("Output: ", OUT_RESULTS),
    "",
    "Inserted frozen cross-cohort taxonomic reproducibility result before the external trauma-validation Results section.",
    "",
    "No statistical model was rerun.",
    "No frozen result was overwritten.",
    "The source Results v5 file remains unchanged.",
    "",
    "STEP96B COMPLETE"
  ),
  file.path(OUT, "03_STEP96B_SUMMARY.txt")
)

writeLines(
  c(
    paste0("Completed: ", Sys.time()),
    "STEP96B COMPLETE",
    "Results v6 created.",
    "No new statistical inference."
  ),
  file.path(OUT, "_STEP96B_COMPLETE.txt")
)

cat("\nSTEP96B COMPLETE\n")
cat("Output:", OUT_RESULTS, "\n")
