
# ============================================================
# V2_96D3 VERIFIED METHODS MANUSCRIPTIZATION
#
# Replaces the four [VERIFY] markers in Methods v1 using the
# source-grounded Step96D2 audit.
#
# NO microbiome processing or statistical model is rerun.
# ============================================================

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(stringr)
})

ROOT <- "E:/sepsis_project"
RESULTS <- file.path(ROOT,"results")

METHODS_V1 <- file.path(
  RESULTS,
  "V2_96D_METHODS_MANUSCRIPTIZATION",
  "Methods_Draft_v1_FROM_FROZEN_V2_ANALYSES.txt"
)

AUDIT_DIR <- file.path(
  RESULTS,
  "V2_96D2_METHODS_SOURCE_VERIFICATION_AUDIT"
)

VERSIONS <- file.path(
  AUDIT_DIR,
  "04_VERIFY4_SOFTWARE_VERSIONS.csv"
)

OUT <- file.path(
  RESULTS,
  "V2_96D3_METHODS_MANUSCRIPTIZATION_VERIFIED"
)
dir.create(OUT,recursive=TRUE,showWarnings=FALSE)

if(!file.exists(METHODS_V1)){
  stop("Methods v1 missing: ",METHODS_V1)
}

required_audit <- c(
  "01B_VERIFY1_UPSTREAM_PROCESSING_EVIDENCE.csv",
  "01C_VERIFY1_CRA002354_PROCESSING_EVIDENCE.csv",
  "02B_STEP88A2_NUMBERED_SOURCE.txt",
  "03B_NUMBERED_SOURCE_94B2C_CRA002354_final_source_longitudinal_analysis_and_freeze.R.txt",
  "04_VERIFY4_SOFTWARE_VERSIONS.csv"
)

for(f in required_audit){
  p <- file.path(AUDIT_DIR,f)
  if(!file.exists(p)){
    stop("Missing Step96D2 source-audit file: ",p)
  }
}

# ------------------------------------------------------------
# Version helper
# ------------------------------------------------------------
ver <- read_csv(
  VERSIONS,
  show_col_types=FALSE
)

get_version <- function(name){
  x <- ver %>% filter(software==name)
  if(nrow(x)<1) return(NA_character_)
  as.character(x$version[[1]])
}

r_version <- get_version("R")
dada2_version <- get_version("dada2")
vegan_version <- get_version("vegan")
lme4_version <- get_version("lme4")
lmertest_version <- get_version("lmerTest")

# The source audit can contain a failed false-positive VSEARCH path.
# Select the row containing the actual vsearch version string.
vsearch_rows <- ver %>%
  filter(
    software=="VSEARCH",
    str_detect(
      version,
      regex("vsearch v[0-9]",ignore_case=TRUE)
    )
  )

vsearch_version <- NA_character_

if(nrow(vsearch_rows)>=1){
  mm <- str_match(
    vsearch_rows$version[[1]],
    regex("vsearch v([0-9.]+)",ignore_case=TRUE)
  )
  if(!is.na(mm[1,2])){
    vsearch_version <- mm[1,2]
  }
}

if(any(is.na(c(
  r_version,
  dada2_version,
  vegan_version,
  lme4_version,
  lmertest_version,
  vsearch_version
)))){
  stop("Could not resolve all core software versions from Step96D2.")
}

# ------------------------------------------------------------
# Source-grounded replacement paragraphs
# ------------------------------------------------------------

verify1_text <- paste0(
  "For the eight taxonomy-ready ASV cohorts, feature tables remained cohort-specific throughout processing and were never merged across different 16S target regions. ",
  "Three cohorts (PRJNA430161, PRJNA691455, and PRJNA978257) entered the canonical workspace with previously validated ASV tables, whereas PRJEB82425, PRJNA516701, PRJNA578267, PRJNA851469, and PRJNA1166732 underwent the project workflow for cohort-specific paired-end DADA2 processing. ",
  "For DADA2-processed cohorts, platform- and cohort-specific filtering parameters were frozen after pilot quality audits, followed by filtering and trimming, cohort-specific error learning from the full cohort, dereplication, sample-wise DADA inference without pooling, paired-read merging, sequence-table construction, and consensus chimera removal. ",
  "Taxonomy was assigned independently within each cohort using the SILVA 138.2 SSU NR99 DADA2 training set. ",
  "The taxonomic feature filter retained Bacteria and Archaea and removed chloroplast-, mitochondria-, other-kingdom-, and kingdom-unassigned features. ",
  "After taxonomic filtering, a minimum of 2,000 retained prokaryotic reads per sample was frozen for the eight taxonomy-ready cohorts; special-route datasets were subjected to separate platform-aware depth audits."
)

verify2_text <- paste0(
  "Longitudinal alpha-diversity and Bray-displacement models were fitted separately within each cohort using linear mixed-effects models with a patient-specific random intercept. ",
  "Time was evaluated both as a categorical factor and, when numeric sampling days were available, as a continuous day variable. ",
  "Full models contained the time term and prespecified cohort-level adjustments, whereas reduced models omitted the time term; global time effects were tested by likelihood-ratio comparison of maximum-likelihood fits (REML=FALSE). ",
  "Alpha-diversity models additionally adjusted for log10 library size. ",
  "Where applicable, PRJEB82425 models adjusted for phenotype and PRJNA430161 and PRJNA1166732 models adjusted for intervention arm. ",
  "For Bray–Curtis displacement, baseline/reference observations were retained for descriptive trajectories but excluded from inferential mixed models because their displacement from themselves was structurally zero; Bray models did not include sequencing-depth adjustment."
)

verify3_text <- paste0(
  "For the CRA002354 primary longitudinal source analysis, linear mixed-effects models included days since each patient's personal baseline, pulmonary versus recorded non-pulmonary source, their interaction, and a patient-specific random intercept. ",
  "The primary source-by-time test compared a full model containing the interaction with a reduced model containing only the corresponding main effects, using maximum-likelihood estimation (REML=FALSE) and a likelihood-ratio test. ",
  "The primary Bray model was fitted to follow-up observations only (days since personal baseline >0), with Bray–Curtis displacement from personal baseline as the outcome. ",
  "The three rarefied-4,000-read alpha-diversity outcomes (Observed OTU97, Shannon OTU97, and Simpson OTU97) were analyzed using the same source-by-time interaction framework, with Benjamini–Hochberg correction across the three alpha-diversity interaction tests. ",
  "Prespecified sensitivity models added available covariates or altered the analysis population/depth definition without changing the primary interaction estimand."
)

verify4_text <- paste0(
  "The source-verification environment used ",
  r_version,
  ", with DADA2 ",
  dada2_version,
  ", vegan ",
  vegan_version,
  ", lme4 ",
  lme4_version,
  ", and lmerTest ",
  lmertest_version,
  ". ",
  "VSEARCH ",
  vsearch_version,
  " was used for the CRA002354 OTU rescue workflow. ",
  "Analysis scripts and frozen evidence outputs were retained under the project workflow to preserve provenance."
)

# ------------------------------------------------------------
# Replace only the exact [VERIFY] paragraphs
# ------------------------------------------------------------

lines <- readLines(
  METHODS_V1,
  warn=FALSE,
  encoding="UTF-8"
)

replace_marker <- function(lines,marker,replacement){

  idx <- which(str_detect(lines,fixed(marker)))

  if(length(idx)!=1){
    stop(
      "Expected exactly one ",
      marker,
      " marker; found ",
      length(idx)
    )
  }

  lines[idx] <- replacement
  lines
}

lines <- replace_marker(
  lines,
  "[VERIFY-1]",
  verify1_text
)

lines <- replace_marker(
  lines,
  "[VERIFY-2]",
  verify2_text
)

lines <- replace_marker(
  lines,
  "[VERIFY-3]",
  verify3_text
)

lines <- replace_marker(
  lines,
  "[VERIFY-4]",
  verify4_text
)

# Hard guard: no unresolved marker remains.
remaining <- lines[
  str_detect(
    lines,
    regex("\\[VERIFY",ignore_case=TRUE)
  )
]

if(length(remaining)>0){
  stop(
    "Unresolved VERIFY marker(s) remain: ",
    paste(remaining,collapse=" | ")
  )
}

METHODS_V2 <- file.path(
  OUT,
  "Methods_Draft_v2_SOURCE_VERIFIED.txt"
)

writeLines(
  lines,
  METHODS_V2,
  useBytes=TRUE
)

# ------------------------------------------------------------
# Source-resolution registry
# ------------------------------------------------------------

resolution <- tibble(
  verify_id=c(
    "VERIFY-1",
    "VERIFY-2",
    "VERIFY-3",
    "VERIFY-4"
  ),
  status="RESOLVED_FROM_LOCAL_SOURCE_AUDIT",
  principal_evidence=c(
    "Step82C/83C/84B/85A3/85C-86B and Step94B1 source evidence",
    "88A2_longitudinal_diversity_and_displacement_FIXED.R",
    "94B2C_CRA002354_final_source_longitudinal_analysis_and_freeze.R",
    "Step96D2 local packageVersion()/VSEARCH --version audit"
  ),
  manuscript_action="VERIFY_MARKER_REPLACED"
)

write_csv(
  resolution,
  file.path(
    OUT,
    "01_VERIFY_RESOLUTION_REGISTRY.csv"
  )
)

versions_used <- tibble(
  software=c(
    "R",
    "dada2",
    "vegan",
    "lme4",
    "lmerTest",
    "VSEARCH"
  ),
  version=c(
    r_version,
    dada2_version,
    vegan_version,
    lme4_version,
    lmertest_version,
    vsearch_version
  )
)

write_csv(
  versions_used,
  file.path(
    OUT,
    "02_CORE_SOFTWARE_VERSIONS_USED_IN_METHODS.csv"
  )
)

# ------------------------------------------------------------
# QC report
# ------------------------------------------------------------

methods_text <- paste(lines,collapse="\n")

qc <- tibble(
  check=c(
    "unresolved_VERIFY_markers",
    "mentions_no_cross_cohort_feature_pooling",
    "mentions_SILVA_138_2",
    "mentions_depth_2000",
    "mentions_REML_FALSE_logic",
    "mentions_patient_random_intercept",
    "mentions_CRA_OTU97",
    "mentions_VSEARCH_version"
  ),
  passed=c(
    !str_detect(methods_text,"\\[VERIFY"),
    str_detect(
      tolower(methods_text),
      "never merged|not pooled|were not pooled"
    ),
    str_detect(methods_text,"SILVA 138.2"),
    str_detect(methods_text,"2,000"),
    str_detect(
      methods_text,
      "REML=FALSE"
    ),
    str_detect(
      tolower(methods_text),
      "patient-specific random intercept"
    ),
    str_detect(
      methods_text,
      "OTU97"
    ),
    str_detect(
      methods_text,
      paste0("VSEARCH ",vsearch_version)
    )
  )
)

write_csv(
  qc,
  file.path(
    OUT,
    "03_METHODS_V2_QC.csv"
  )
)

if(!all(qc$passed)){
  stop(
    "Methods v2 QC failed: ",
    paste(
      qc$check[!qc$passed],
      collapse=", "
    )
  )
}

writeLines(
  c(
    "V2 STEP96D3 VERIFIED METHODS MANUSCRIPTIZATION",
    "",
    "All four Step96D [VERIFY] markers were resolved from the Step96D2 local source audit.",
    "",
    paste0("Output: ",METHODS_V2),
    "",
    "No microbiome processing was rerun.",
    "No statistical model was rerun.",
    "No frozen result was overwritten.",
    "",
    "STEP96D3 COMPLETE"
  ),
  file.path(
    OUT,
    "04_STEP96D3_SUMMARY.txt"
  )
)

writeLines(
  c(
    paste0("Completed: ",Sys.time()),
    "STEP96D3 COMPLETE",
    "Methods v2 source-verified.",
    "Unresolved VERIFY markers: 0"
  ),
  file.path(
    OUT,
    "_STEP96D3_COMPLETE.txt"
  )
)

cat("STEP96D3 COMPLETE\n")
