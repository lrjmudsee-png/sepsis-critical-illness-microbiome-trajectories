# ============================================================
# Sepsis V2 - Step95F3
# FINAL RESULTS + TABLE 1 FREEZE WITH SENTENCE-LEVEL GUARDRAIL AUDIT
#
# Supersedes Step95F2.
#
# Step95F2 failed only because its regex still matched:
# "Exploratory clustering did not support stable trajectory subtypes..."
#
# Step95F3:
# - copies Step95F2 Table 1 unchanged;
# - uses Step95F2 Results v4 as the basis;
# - makes one additional style-only cleanup:
#   "harmonized primary analysis object" -> "harmonized primary analysis set";
# - performs sentence-level overclaim checking with explicit negation handling;
# - freezes Table 1 + Results if all checks pass.
#
# NO new inferential analysis.
# NO changed numerical results.
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

F2 <- file.path(
  RESULTS,
  "V2_35F2_RESULTS_AND_TABLE1_FINAL_FREEZE"
)

F <- file.path(
  RESULTS,
  "V2_35F_RESULTS_POLISH_AND_TABLE1_MANUSCRIPT_FREEZE"
)

OUT <- file.path(
  RESULTS,
  "V2_35F3_RESULTS_AND_TABLE1_FINAL_FREEZE_SENTENCE_AUDIT"
)

TABDIR <- file.path(OUT, "01_FINAL_TABLES")
DRAFTDIR <- file.path(OUT, "02_FINAL_RESULTS")
AUDDIR <- file.path(OUT, "03_FINAL_AUDIT")

for (d in c(OUT,TABDIR,DRAFTDIR,AUDDIR)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

required_f2 <- c(
  file.path(
    F2,
    "01_FINAL_TABLES",
    "Main_Table_1_Cohort_Characteristics_MANUSCRIPT.csv"
  ),
  file.path(
    F2,
    "01_FINAL_TABLES",
    "Main_Table_1_Cohort_Characteristics_MANUSCRIPT_Footnotes.txt"
  ),
  file.path(
    F2,
    "02_FINAL_RESULTS",
    "Results_Draft_v4_FINAL_FREEZE.txt"
  )
)

if (!all(file.exists(required_f2))) {
  stop("Required Step95F2 manuscript outputs are missing.")
}

# ------------------------------------------------------------
# 1. Copy Table 1 unchanged
# ------------------------------------------------------------

table_files <- c(
  "Main_Table_1_Cohort_Characteristics_MANUSCRIPT.csv",
  "Main_Table_1_Cohort_Characteristics_MANUSCRIPT_Footnotes.txt"
)

copy_audit <- tibble(
  file = table_files,
  copied = FALSE
)

for (i in seq_along(table_files)) {

  src <- file.path(
    F2,
    "01_FINAL_TABLES",
    table_files[i]
  )

  dst <- file.path(
    TABDIR,
    table_files[i]
  )

  copy_audit$copied[i] <- file.copy(
    src,
    dst,
    overwrite = TRUE
  )
}

write_csv(
  copy_audit,
  file.path(AUDDIR, "01_TABLE1_COPY_AUDIT.csv")
)

if (!all(copy_audit$copied)) {
  stop("Table 1 copy audit failed.")
}

# ------------------------------------------------------------
# 2. Results wording-only cleanup
# ------------------------------------------------------------

draft <- paste(
  readLines(
    file.path(
      F2,
      "02_FINAL_RESULTS",
      "Results_Draft_v4_FINAL_FREEZE.txt"
    ),
    warn = FALSE,
    encoding = "UTF-8"
  ),
  collapse = "\n"
)

old_phrase <- "harmonized primary analysis object"
new_phrase <- "harmonized primary analysis set"

if (!str_detect(draft, fixed(old_phrase))) {
  stop("Expected internal phrase not found; refusing silent replacement.")
}

draft_final <- str_replace(
  draft,
  fixed(old_phrase),
  new_phrase
)

writeLines(
  draft_final,
  file.path(
    DRAFTDIR,
    "Results_Draft_v5_FINAL_MANUSCRIPT_FREEZE.txt"
  )
)

# ------------------------------------------------------------
# 3. Sentence-level affirmative-overclaim audit
# ------------------------------------------------------------

# Collapse section breaks for sentence splitting, preserving content.
audit_text <- str_replace_all(draft_final, "\n+", " ")

sentences <- unlist(
  str_split(
    audit_text,
    "(?<=[.!?])\\s+"
  )
)

sentences <- str_trim(sentences)
sentences <- sentences[sentences != ""]

is_negated <- function(sentence) {

  str_detect(
    sentence,
    regex(
      paste0(
        "\\b(",
        "did not|does not|do not|not support|not establish|not demonstrate|",
        "not confirm|should not|cannot|could not|no evidence|no clear|",
        "failed to|neither|nor|without establishing|rather than",
        ")\\b"
      ),
      ignore_case = TRUE
    )
  )
}

rules <- tribble(
  ~guardrail, ~concept_pattern, ~affirmative_pattern,

  "SEPSIS_SPECIFIC_LONGITUDINAL_INSTABILITY",
  "sepsis-specific longitudinal instability|specific to sepsis",
  "demonstrat|establish|confirm|prove|show",

  "STABLE_TRAJECTORY_SUBTYPE",
  "stable (taxonomic )?trajectory subtype|stable trajectory subtypes",
  "identify|demonstrat|establish|confirm|support",

  "EII_AS_VALIDATED_CLINICAL_SCORE",
  "ecological injury index|\\bEII\\b",
  "clinical|prognostic|validated|predict",

  "SOURCE_TRAJECTORY_EQUIVALENCE",
  "pulmonary.*non-pulmonary|non-pulmonary.*pulmonary",
  "identical|equivalent|the same",

  "ABSOLUTE_NO_SOURCE_EFFECT",
  "infection source",
  "has no effect|had no effect|does not affect",

  "ENTEROCOCCUS_ROBUST_BIOMARKER",
  "Enterococcus",
  "robust|validated|biomarker|signature",

  "BACTEROIDES_VALIDATED_LONGITUDINAL_SIGNAL",
  "Bacteroides",
  "validated|robust|trajectory|longitudinal signal"
)

audit_rows <- list()

for (i in seq_len(nrow(rules))) {

  concept_hits <- sentences[
    str_detect(
      sentences,
      regex(
        rules$concept_pattern[i],
        ignore_case = TRUE
      )
    )
  ]

  if (length(concept_hits) == 0) {

    audit_rows[[length(audit_rows)+1]] <- tibble(
      guardrail = rules$guardrail[i],
      sentence = NA_character_,
      concept_present = FALSE,
      affirmative_language_present = FALSE,
      negation_present = FALSE,
      overclaim_detected = FALSE
    )

    next
  }

  for (s in concept_hits) {

    aff <- str_detect(
      s,
      regex(
        rules$affirmative_pattern[i],
        ignore_case = TRUE
      )
    )

    neg <- is_negated(s)

    audit_rows[[length(audit_rows)+1]] <- tibble(
      guardrail = rules$guardrail[i],
      sentence = s,
      concept_present = TRUE,
      affirmative_language_present = aff,
      negation_present = neg,
      overclaim_detected = aff && !neg
    )
  }
}

guard_audit <- bind_rows(audit_rows)

write_csv(
  guard_audit,
  file.path(
    AUDDIR,
    "02_RESULTS_SENTENCE_LEVEL_OVERCLAIM_AUDIT.csv"
  )
)

if (any(guard_audit$overclaim_detected)) {
  stop("Sentence-level audit detected an affirmative manuscript overclaim.")
}

# ------------------------------------------------------------
# 4. Copy table-placement and supplementary-table plans
# ------------------------------------------------------------

plans <- c(
  "FINAL_TABLE_PLACEMENT.csv",
  "SUPPLEMENTARY_TABLE_PLAN.csv"
)

plan_audit <- tibble(
  file = plans,
  copied = FALSE
)

for (i in seq_along(plans)) {

  # Step95F produced these before its QC false-positive stop.
  src <- file.path(
    F,
    "03_AUDIT_AND_PLACEMENT",
    plans[i]
  )

  if (!file.exists(src)) {
    stop(
      paste0(
        "Required Step95F plan is missing: ",
        plans[i]
      )
    )
  }

  plan_audit$copied[i] <- file.copy(
    src,
    file.path(AUDDIR, plans[i]),
    overwrite = TRUE
  )
}

write_csv(
  plan_audit,
  file.path(
    AUDDIR,
    "03_PLAN_COPY_AUDIT.csv"
  )
)

if (!all(plan_audit$copied)) {
  stop("Placement/supplementary plan copy failed.")
}

# ------------------------------------------------------------
# 5. Change / freeze note
# ------------------------------------------------------------

writeLines(
  c(
    "STEP95F3 FINAL FREEZE NOTE",
    "",
    "Why Step95F2 stopped:",
    "Its regex matched 'support stable trajectory subtypes' inside the explicitly negated sentence 'did not support stable trajectory subtypes'.",
    "",
    "Step95F3 audit behavior:",
    "Overclaim rules are evaluated at the sentence level.",
    "A candidate statement is only flagged when affirmative language is present without an explicit negation/limitation context.",
    "",
    "Additional style-only edit:",
    paste0("OLD: ", old_phrase),
    paste0("NEW: ", new_phrase),
    "",
    "Table 1 was copied unchanged.",
    "All numeric Results values were copied unchanged.",
    "No effect estimate, 95% CI, p value, FDR, evidence tier, or scientific conclusion was altered.",
    "",
    "Step95F3 is the authoritative manuscript-facing Table 1 + Results freeze."
  ),
  file.path(
    AUDDIR,
    "STEP95F3_FINAL_FREEZE_NOTE.txt"
  )
)

# ------------------------------------------------------------
# 6. Final readiness
# ------------------------------------------------------------

required_final <- c(
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
    "Results_Draft_v5_FINAL_MANUSCRIPT_FREEZE.txt"
  ),
  file.path(
    AUDDIR,
    "02_RESULTS_SENTENCE_LEVEL_OVERCLAIM_AUDIT.csv"
  ),
  file.path(
    AUDDIR,
    "FINAL_TABLE_PLACEMENT.csv"
  ),
  file.path(
    AUDDIR,
    "SUPPLEMENTARY_TABLE_PLAN.csv"
  )
)

ready <- all(file.exists(required_final)) &&
  !any(guard_audit$overclaim_detected)

qc <- tibble(
  manuscript_table1_frozen = TRUE,
  results_v5_frozen = TRUE,
  sentence_level_overclaim_audit_pass =
    !any(guard_audit$overclaim_detected),
  table1_numeric_values_changed = FALSE,
  results_numeric_values_changed = FALSE,
  inferential_analysis_rerun = FALSE,
  ready_for_supplementary_tables_and_discussion = ready
)

write_csv(
  qc,
  file.path(
    OUT,
    "STEP95F3_READINESS.csv"
  )
)

if (!ready) {
  stop("Step95F3 final readiness failed.")
}

writeLines(
  c(
    paste0("Completed: ", Sys.time()),
    "Manuscript Table 1 frozen: TRUE",
    "Results v5 frozen: TRUE",
    "Sentence-level overclaim audit passed: TRUE",
    "Table 1 numeric values changed: FALSE",
    "Results numeric values changed: FALSE",
    "No inferential analysis rerun: TRUE",
    "Ready for Supplementary Tables + Discussion: TRUE",
    "STEP95F3 COMPLETE"
  ),
  file.path(
    OUT,
    "_STEP95F3_COMPLETE.ok"
  )
)

cat("STEP95F3 COMPLETE\n")
