
# ============================================================
# V2_96D4 METHODS FINAL MANUSCRIPT POLISH
#
# Purpose:
# Convert source-verified Methods v2 into manuscript-facing v3.
#
# This step ONLY:
# - removes internal workflow labels;
# - fixes draft heading/version language;
# - removes redundant software wording;
# - improves transparency of cohort selection wording.
#
# NO analysis rerun.
# NO numerical/statistical changes.
# ============================================================

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(stringr)
})

ROOT <- "E:/sepsis_project"
RESULTS <- file.path(ROOT,"results")

INFILE <- file.path(
  RESULTS,
  "V2_96D3_METHODS_MANUSCRIPTIZATION_VERIFIED",
  "Methods_Draft_v2_SOURCE_VERIFIED.txt"
)

OUT <- file.path(
  RESULTS,
  "V2_96D4_METHODS_FINAL_MANUSCRIPT_POLISH"
)

dir.create(OUT, recursive=TRUE, showWarnings=FALSE)

if(!file.exists(INFILE)){
  stop("Methods v2 not found: ",INFILE)
}

lines <- readLines(
  INFILE,
  warn=FALSE,
  encoding="UTF-8"
)

txt <- paste(lines,collapse="\n")

# ------------------------------------------------------------
# 1. Remove internal draft/workflow language
# ------------------------------------------------------------
txt <- str_replace(
  txt,
  "^METHODS DRAFT V1 — MANUSCRIPTIZATION FROM FROZEN V2 ANALYSES",
  "Methods"
)

txt <- str_replace_all(
  txt,
  fixed(
    "paired taxonomic effects were first estimated independently within each cohort using the frozen Step89A paired CLR analysis."
  ),
  "paired taxonomic effects were first estimated independently within each cohort using within-patient paired centered-log-ratio (CLR) contrasts."
)

txt <- str_replace_all(
  txt,
  fixed(
    "The primary cross-cohort analysis used the ALL_PAIRED definition, while COMMON_ANCHOR_SENSITIVITY was retained as a prespecified sensitivity analysis."
  ),
  "The primary cross-cohort analysis used all eligible within-patient early-to-late pairs. A common-anchor analysis, restricted to comparisons sharing the prespecified reference time point, was retained as a sensitivity analysis."
)

# ------------------------------------------------------------
# 2. Make cohort-selection wording transparent
# ------------------------------------------------------------
txt <- str_replace_all(
  txt,
  fixed(
    "For each pair of the three primary progressive cohorts (PRJNA691455, PRJNA851469, and PRJNA516701), reproducibility was summarized using the Spearman correlation between shared-taxon CLR effect vectors, the proportion of shared taxa changing in the same direction, and the Jaccard overlap among taxa with the largest absolute effects."
  ),
  "Cross-cohort taxonomic reproducibility was evaluated among the three longitudinal cohorts classified in the frozen evidence synthesis as showing robust progressive ecological displacement (PRJNA691455, PRJNA851469, and PRJNA516701). For each cohort pair, reproducibility was summarized using the Spearman correlation between shared-taxon CLR effect vectors, the proportion of shared taxa changing in the same direction, and the Jaccard overlap among taxa with the largest absolute effects."
)

# ------------------------------------------------------------
# 3. Remove duplicated generic software paragraph
# ------------------------------------------------------------
dup_para <- paste0(
  "Data processing and statistical analyses were performed in R (version 4.4.0) with external command-line tools used where required for sequence processing. ",
  "VSEARCH was used for the CRA002354 97% OTU reconstruction. ",
  "All manuscript-level results were generated from frozen cohort-specific analysis objects and evidence-freeze outputs. ",
  "No ASV or OTU feature matrix was directly merged across heterogeneous cohorts.\n\n"
)

txt <- str_replace(
  txt,
  fixed(dup_para),
  ""
)

# Strengthen remaining software paragraph with reproducibility guardrail.
old_soft <- paste0(
  "The source-verification environment used R version 4.4.0 (2024-04-24 ucrt), with DADA2 1.34.0, vegan 2.7.3, lme4 2.0.1, and lmerTest 3.2.1. ",
  "VSEARCH 2.31.0 was used for the CRA002354 OTU rescue workflow. ",
  "Analysis scripts and frozen evidence outputs were retained under the project workflow to preserve provenance."
)

new_soft <- paste0(
  "Analyses were performed in R version 4.4.0 (2024-04-24 ucrt), with DADA2 1.34.0, vegan 2.7.3, lme4 2.0.1, and lmerTest 3.2.1. ",
  "VSEARCH 2.31.0 was used for the CRA002354 97% OTU rescue workflow. ",
  "Cohort-specific feature matrices remained separate throughout the analysis, and manuscript-level synthesis was performed using frozen cohort-specific ecological or taxonomic effect estimates."
)

txt <- str_replace(
  txt,
  fixed(old_soft),
  new_soft
)

# ------------------------------------------------------------
# 4. Replace remaining internal "frozen" language where it is
#    provenance-related rather than scientifically necessary
# ------------------------------------------------------------
txt <- str_replace_all(
  txt,
  fixed(
    "Replicate handling was frozen before longitudinal inference."
  ),
  "Replicate handling rules were finalized before longitudinal inference."
)

txt <- str_replace_all(
  txt,
  fixed(
    "platform- and cohort-specific filtering parameters were frozen after pilot quality audits"
  ),
  "platform- and cohort-specific filtering parameters were finalized after pilot quality audits"
)

txt <- str_replace_all(
  txt,
  fixed(
    "a minimum of 2,000 retained prokaryotic reads per sample was frozen for the eight taxonomy-ready cohorts"
  ),
  "a minimum of 2,000 retained prokaryotic reads per sample was required for the eight taxonomy-ready cohorts"
)

# ------------------------------------------------------------
# 5. Manuscript-facing terminology cleanup
# ------------------------------------------------------------
txt <- str_replace_all(
  txt,
  "OTU rescue workflow",
  "OTU reconstruction workflow"
)

txt <- str_replace_all(
  txt,
  "taxonomy-ready ASV cohorts",
  "ASV-based cohorts"
)

# ------------------------------------------------------------
# 6. Hard QC
# ------------------------------------------------------------
internal_patterns <- c(
  "\\[VERIFY",
  "Step89A",
  "Step88",
  "Step94",
  "Step96",
  "ALL_PAIRED",
  "COMMON_ANCHOR_SENSITIVITY",
  "METHODS DRAFT V1"
)

internal_hits <- tibble(
  pattern=internal_patterns,
  hit=vapply(
    internal_patterns,
    function(p) str_detect(txt,regex(p,ignore_case=TRUE)),
    logical(1)
  )
)

write_csv(
  internal_hits,
  file.path(
    OUT,
    "01_INTERNAL_WORKFLOW_LANGUAGE_QC.csv"
  )
)

if(any(internal_hits$hit)){
  stop(
    "Internal workflow language remains: ",
    paste(
      internal_hits$pattern[internal_hits$hit],
      collapse=", "
    )
  )
}

required_strings <- c(
  "SILVA 138.2",
  "2,000",
  "REML=FALSE",
  "patient-specific random intercept",
  "5,000 permutations",
  "VSEARCH 2.31.0",
  "DADA2 1.34.0",
  "PRJNA691455",
  "PRJNA851469",
  "PRJNA516701"
)

required_qc <- tibble(
  item=required_strings,
  present=vapply(
    required_strings,
    function(s) str_detect(txt,fixed(s)),
    logical(1)
  )
)

write_csv(
  required_qc,
  file.path(
    OUT,
    "02_METHODS_CONTENT_PRESERVATION_QC.csv"
  )
)

if(!all(required_qc$present)){
  stop(
    "Required Methods content missing: ",
    paste(
      required_qc$item[!required_qc$present],
      collapse=", "
    )
  )
}

# ------------------------------------------------------------
# 7. Write final manuscript-facing Methods
# ------------------------------------------------------------
OUTFILE <- file.path(
  OUT,
  "Methods_v3_MANUSCRIPT_READY.txt"
)

writeLines(
  str_split(txt,"\n",simplify=FALSE)[[1]],
  OUTFILE,
  useBytes=TRUE
)

audit <- tibble(
  source_file=INFILE,
  output_file=OUTFILE,
  analysis_rerun=FALSE,
  numeric_result_changed=FALSE,
  verify_markers_remaining=FALSE,
  internal_step_labels_remaining=FALSE,
  manuscript_ready=TRUE
)

write_csv(
  audit,
  file.path(
    OUT,
    "03_METHODS_V3_FINAL_AUDIT.csv"
  )
)

writeLines(
  c(
    "V2 STEP96D4 METHODS FINAL MANUSCRIPT POLISH",
    "",
    "Methods v2 source verification: retained.",
    "Internal Step labels: removed.",
    "VERIFY markers: 0.",
    "Software duplication: removed.",
    "Cross-cohort cohort-selection wording: made transparent.",
    "Statistical/microbiome analysis rerun: NO.",
    "Numerical results changed: NO.",
    "",
    paste0("Final Methods: ",OUTFILE),
    "",
    "STEP96D4 COMPLETE"
  ),
  file.path(
    OUT,
    "04_STEP96D4_SUMMARY.txt"
  )
)

writeLines(
  c(
    paste0("Completed: ",Sys.time()),
    "STEP96D4 COMPLETE",
    "Methods v3 manuscript-ready."
  ),
  file.path(
    OUT,
    "_STEP96D4_COMPLETE.txt"
  )
)

cat("STEP96D4 COMPLETE\n")
