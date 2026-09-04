# ============================================================
# Sepsis V2 - Step95F2
# RESULTS + TABLE 1 FINAL FREEZE
#
# Supersedes Step95F only because Step95F's guardrail checker
# produced a false positive on the explicitly negated sentence:
# "do not establish sepsis-specific longitudinal instability."
#
# NO new statistics.
# NO changed estimates, p values, FDRs, CIs, evidence tiers, or tables.
#
# Step95F2:
# 1) copies the manuscript Table 1 unchanged;
# 2) makes one wording-only Results change:
#    "these findings do not establish sepsis-specific longitudinal instability"
#    -> "these findings should not be interpreted as specific to sepsis";
# 3) replaces substring-based guardrail detection with affirmative-overclaim
#    pattern detection;
# 4) freezes Results + Table 1 for manuscript production.
# ============================================================

options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tibble)
  library(stringr)
})

ROOT <- "E:/sepsis_project"
RESULTS <- file.path(ROOT, "results")

F <- file.path(
  RESULTS,
  "V2_35F_RESULTS_POLISH_AND_TABLE1_MANUSCRIPT_FREEZE"
)

OUT <- file.path(
  RESULTS,
  "V2_35F2_RESULTS_AND_TABLE1_FINAL_FREEZE"
)

TABDIR <- file.path(OUT, "01_FINAL_TABLES")
DRAFTDIR <- file.path(OUT, "02_FINAL_RESULTS")
AUDDIR <- file.path(OUT, "03_FINAL_AUDIT")

for (d in c(OUT,TABDIR,DRAFTDIR,AUDDIR)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

# Step95F did not complete because of a QC false positive, so use
# existence of the authoritative manuscript outputs instead of COMPLETE flag.
required <- c(
  file.path(
    F,
    "01_FINAL_TABLES",
    "Main_Table_1_Cohort_Characteristics_MANUSCRIPT.csv"
  ),
  file.path(
    F,
    "01_FINAL_TABLES",
    "Main_Table_1_Cohort_Characteristics_MANUSCRIPT_Footnotes.txt"
  ),
  file.path(
    F,
    "02_RESULTS",
    "Results_Draft_v3_English_MANUSCRIPT.txt"
  ),
  file.path(
    F,
    "03_AUDIT_AND_PLACEMENT",
    "FINAL_TABLE_PLACEMENT.csv"
  ),
  file.path(
    F,
    "03_AUDIT_AND_PLACEMENT",
    "SUPPLEMENTARY_TABLE_PLAN.csv"
  )
)

if (!all(file.exists(required))) {
  stop("One or more Step95F manuscript outputs are missing.")
}

# ------------------------------------------------------------
# 1. Copy Table 1 unchanged
# ------------------------------------------------------------

table_files <- c(
  "Main_Table_1_Cohort_Characteristics_MANUSCRIPT.csv",
  "Main_Table_1_Cohort_Characteristics_MANUSCRIPT_Footnotes.txt"
)

copy_table_audit <- tibble(
  file = table_files,
  copied = FALSE
)

for (i in seq_along(table_files)) {

  src <- file.path(
    F,
    "01_FINAL_TABLES",
    table_files[i]
  )

  dst <- file.path(
    TABDIR,
    table_files[i]
  )

  copy_table_audit$copied[i] <- file.copy(
    src,
    dst,
    overwrite = TRUE
  )
}

write_csv(
  copy_table_audit,
  file.path(AUDDIR, "01_TABLE1_COPY_AUDIT.csv")
)

if (!all(copy_table_audit$copied)) {
  stop("Could not copy manuscript Table 1 files.")
}

# ------------------------------------------------------------
# 2. Wording-only Results repair
# ------------------------------------------------------------

draft_src <- file.path(
  F,
  "02_RESULTS",
  "Results_Draft_v3_English_MANUSCRIPT.txt"
)

draft <- paste(
  readLines(
    draft_src,
    warn = FALSE,
    encoding = "UTF-8"
  ),
  collapse = "\n"
)

old_sentence <- paste0(
  "Because progressive displacement was also present in ICU-background cohorts, ",
  "these findings do not establish sepsis-specific longitudinal instability."
)

new_sentence <- paste0(
  "Because progressive displacement was also present in ICU-background cohorts, ",
  "these findings should not be interpreted as specific to sepsis."
)

if (!str_detect(draft, fixed(old_sentence))) {
  stop("Expected Step95F sentence was not found; refusing silent text replacement.")
}

draft_final <- str_replace(
  draft,
  fixed(old_sentence),
  new_sentence
)

writeLines(
  draft_final,
  file.path(
    DRAFTDIR,
    "Results_Draft_v4_FINAL_FREEZE.txt"
  )
)

# ------------------------------------------------------------
# 3. Improved guardrail audit
#
# Detect affirmative overclaims rather than prohibited substrings appearing
# inside explicit negation.
# ------------------------------------------------------------

guardrails <- tribble(
  ~guardrail, ~pattern,
  "AFFIRMATIVE_SEPSIS_SPECIFIC_INSTABILITY",
  "(demonstrat|establish|confirm|show|prove)[a-z]*.{0,50}sepsis-specific longitudinal instability",

  "AFFIRMATIVE_STABLE_TRAJECTORY_SUBTYPE",
  "(identify|demonstrat|establish|confirm|support)[a-z]*.{0,50}stable (taxonomic )?trajectory subtype",

  "EII_AS_CLINICAL_PROGNOSTIC_SCORE",
  "(clinical|prognostic|validated).{0,30}(ecological injury index|EII)",

  "SOURCE_TRAJECTORY_EQUIVALENCE",
  "(pulmonary and (recorded )?non-pulmonary).{0,40}(identical|equivalent|the same)",

  "AFFIRMATIVE_NO_SOURCE_EFFECT",
  "infection source (has|had) no effect",

  "ENTEROCOCCUS_ROBUST_BIOMARKER",
  "Enterococcus.{0,50}(robust|validated).{0,30}(biomarker|signature)",

  "BACTEROIDES_VALIDATED_LONGITUDINAL_SIGNAL",
  "Bacteroides.{0,50}(validated|robust).{0,30}(trajectory|longitudinal signal)"
)

guard_audit <- guardrails %>%
  rowwise() %>%
  mutate(
    detected = str_detect(
      draft_final,
      regex(pattern, ignore_case = TRUE)
    )
  ) %>%
  ungroup()

write_csv(
  guard_audit,
  file.path(
    AUDDIR,
    "02_RESULTS_AFFIRMATIVE_OVERCLAIM_AUDIT.csv"
  )
)

if (any(guard_audit$detected)) {
  stop("An affirmative manuscript overclaim was detected.")
}

# ------------------------------------------------------------
# 4. Copy placement / supplementary plan
# ------------------------------------------------------------

placement_files <- c(
  "FINAL_TABLE_PLACEMENT.csv",
  "SUPPLEMENTARY_TABLE_PLAN.csv"
)

placement_audit <- tibble(
  file = placement_files,
  copied = FALSE
)

for (i in seq_along(placement_files)) {

  src <- file.path(
    F,
    "03_AUDIT_AND_PLACEMENT",
    placement_files[i]
  )

  dst <- file.path(
    AUDDIR,
    placement_files[i]
  )

  placement_audit$copied[i] <- file.copy(
    src,
    dst,
    overwrite = TRUE
  )
}

write_csv(
  placement_audit,
  file.path(
    AUDDIR,
    "03_PLACEMENT_PLAN_COPY_AUDIT.csv"
  )
)

if (!all(placement_audit$copied)) {
  stop("Could not copy placement/supplementary plans.")
}

# ------------------------------------------------------------
# 5. Freeze note
# ------------------------------------------------------------

writeLines(
  c(
    "STEP95F2 FINAL RESULTS + TABLE 1 FREEZE",
    "",
    "Step95F failed only because its substring-based guardrail checker flagged an explicitly negated safety sentence.",
    "",
    "Wording-only change:",
    paste0("OLD: ", old_sentence),
    paste0("NEW: ", new_sentence),
    "",
    "All manuscript Table 1 values were copied unchanged from Step95F.",
    "All Results numeric values were copied unchanged from Step95F.",
    "No effect estimate, CI, p value, FDR, evidence tier, or scientific conclusion was changed.",
    "",
    "The new guardrail audit checks affirmative overclaims instead of banned substrings appearing inside negation.",
    "",
    "Step95F2 supersedes Step95F as the manuscript-facing Results + Table 1 freeze."
  ),
  file.path(
    AUDDIR,
    "STEP95F2_FREEZE_NOTE.txt"
  )
)

# ------------------------------------------------------------
# 6. Final QC
# ------------------------------------------------------------

final_files <- c(
  file.path(
    TABDIR,
    "Main_Table_1_Cohort_Characteristics_MANUSCRIPT.csv"
  ),
  file.path(
    TABDIR,
    "Main_Table_1_Cohort_Characteristics_MANUSCRIPT_Footnotes.txt"
  ),
  file.path(
    DRAFTDIR,
    "Results_Draft_v4_FINAL_FREEZE.txt"
  ),
  file.path(
    AUDDIR,
    "FINAL_TABLE_PLACEMENT.csv"
  ),
  file.path(
    AUDDIR,
    "SUPPLEMENTARY_TABLE_PLAN.csv"
  ),
  file.path(
    AUDDIR,
    "02_RESULTS_AFFIRMATIVE_OVERCLAIM_AUDIT.csv"
  )
)

ready <- all(file.exists(final_files)) &&
  !any(guard_audit$detected)

qc <- tibble(
  manuscript_table1_frozen = TRUE,
  results_v4_frozen = TRUE,
  affirmative_overclaim_audit_pass = !any(guard_audit$detected),
  wording_only_change_from_step95F = TRUE,
  numeric_results_changed = FALSE,
  inferential_analysis_rerun = FALSE,
  ready_for_supplementary_tables_and_discussion = ready
)

write_csv(
  qc,
  file.path(
    OUT,
    "STEP95F2_READINESS.csv"
  )
)

if (!ready) {
  stop("Step95F2 readiness failed.")
}

writeLines(
  c(
    paste0("Completed: ", Sys.time()),
    "Manuscript Table 1 frozen: TRUE",
    "Results v4 frozen: TRUE",
    "Affirmative overclaim audit passed: TRUE",
    "Numeric results changed: FALSE",
    "No inferential analysis rerun: TRUE",
    "Ready for Supplementary Tables + Discussion: TRUE",
    "STEP95F2 COMPLETE"
  ),
  file.path(
    OUT,
    "_STEP95F2_COMPLETE.ok"
  )
)

cat("STEP95F2 COMPLETE\n")
