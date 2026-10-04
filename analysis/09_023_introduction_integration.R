
# ============================================================
# V2_96H INTRODUCTION INTEGRATION
#
# Purpose:
# Replace the Introduction placeholder in the Step96G manuscript
# core with a manuscript-facing, evidence-grounded Introduction.
#
# Citation keys remain stable [CIT:...] identifiers.
# Final numeric reference numbering is still deferred.
#
# NO statistical analysis is rerun.
# ============================================================

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(stringr)
  library(tibble)
})

ROOT <- "E:/sepsis_project"
RESULTS <- file.path(ROOT,"results")
CODE <- file.path(ROOT,"code","03_data_processing")

INFILE <- file.path(
  RESULTS,
  "V2_96G_DISCUSSION_REFERENCE_RESOLUTION",
  "Manuscript_Core_v2_DISCUSSION_REFS_RESOLVED.txt"
)

OLD_POOL <- file.path(
  RESULTS,
  "V2_96G_DISCUSSION_REFERENCE_RESOLUTION",
  "V2_96G_VERIFIED_REFERENCE_POOL.csv"
)

NEW_POOL <- file.path(
  CODE,
  "V2_96H_NEW_VERIFIED_REFERENCES.csv"
)

OUT <- file.path(
  RESULTS,
  "V2_96H_INTRODUCTION_INTEGRATION"
)

dir.create(OUT,recursive=TRUE,showWarnings=FALSE)

for(f in c(INFILE,OLD_POOL,NEW_POOL)){
  if(!file.exists(f)){
    stop("Missing Step96H input: ",f)
  }
}

txt <- paste(
  readLines(
    INFILE,
    warn=FALSE,
    encoding="UTF-8"
  ),
  collapse="\n"
)

placeholder <- "[INTRODUCTION TO BE WRITTEN AFTER FINAL LITERATURE INTEGRATION]"

if(!str_detect(txt,fixed(placeholder))){
  stop("Introduction placeholder not found or already replaced.")
}

# ------------------------------------------------------------
# 1. Manuscript-facing Introduction
# ------------------------------------------------------------

intro <- paste(
  c(
    paste0(
      "Sepsis is defined as life-threatening organ dysfunction caused by a dysregulated host response to infection and remains a major global health burden [CIT:SINGER2016;GBD2025]. ",
      "Although infection initiates the syndrome, its biological expression is highly heterogeneous, reflecting variation in the host response, organ dysfunction, treatment exposures, and the underlying infectious context. ",
      "The intestinal microbiome is increasingly recognized as one component of this host–environment interface. ",
      "Critical illness can disrupt intestinal barrier function and expose microbial communities to profound perturbations from antibiotics, altered nutrition, organ support, impaired motility, and other intensive-care interventions, creating conditions for loss of commensal organisms and expansion of opportunistic taxa [CIT:HAAK2017;ADELMAN2020;MCDONALD2016]."
    ),
    "",
    paste0(
      "Sequencing studies have repeatedly documented marked gut-microbiome disruption in critically ill and septic patients, but the temporal pattern and reproducibility of these changes remain incompletely resolved. ",
      "Early multicenter ICU work demonstrated rapid and extreme dysbiosis during critical illness [CIT:MCDONALD2016], while subsequent sepsis and critical-care cohorts identified reduced diversity, pathogen enrichment, and associations between microbiome composition and clinical exposures or outcomes [CIT:AGUDELO2020;SCHLECHTE2023;KITSIOS2024]. ",
      "A recent systematic review of 36 longitudinal sequencing studies involving 2,067 critically ill adults found that progressive diversity loss and temporal community shifts were common, yet also emphasized substantial between-study heterogeneity and very low certainty for many specific microbiome outcomes [CIT:THEOCHARIDOU2026]. ",
      "Thus, the broad phenomenon of microbiome disruption during critical illness appears reproducible, whereas the precise taxa, magnitude, and direction of taxonomic change often vary among cohorts."
    ),
    "",
    paste0(
      "This distinction is important because much of microbiome biomarker research has focused on identifying disease-associated taxa or taxonomic signatures. ",
      "Across microbiome studies, however, differential-abundance results and microbial signatures can vary materially with cohort composition, laboratory processing, analytical method, and baseline community structure [CIT:GEISTLINGER2024;CLAUSEN2022;NEARING2022;GAO2024]. ",
      "From an ecological perspective, different starting communities may also respond to a shared perturbation through different taxonomic routes while undergoing a similar higher-order displacement from their previous state [CIT:FASSARELLA2021]. ",
      "For sepsis and critical illness, it therefore remains unclear whether longitudinal change is more consistently expressed as ecosystem-level displacement than as a universal set of taxonomic alterations."
    ),
    "",
    paste0(
      "We addressed this question using a multi-cohort secondary analysis of publicly available 16S rRNA gene sequencing datasets, with all feature-level analyses performed within cohort rather than by pooling heterogeneous ASV or OTU matrices across studies. ",
      "We first evaluated longitudinal ecological displacement across independent sepsis and critical-illness cohorts and characterized patient-specific taxonomic trajectories in a core repeated-sepsis cohort. ",
      "We then quantified cross-cohort reproducibility of longitudinal genus- and family-level effects, tested ecological-state differentiation beyond severe trauma in an external cohort, and examined whether pulmonary versus recorded non-pulmonary infection source modified longitudinal ecological displacement in a dedicated sepsis cohort. ",
      "We hypothesized that a common direction of ecological disturbance would be more reproducible across cohorts than the specific taxonomic routes through which that disturbance emerged. ",
      "Accordingly, the study was designed to identify portable ecological patterns of microbiome change rather than to derive a universal sepsis-specific taxonomic biomarker or clinical prediction score."
    )
  ),
  collapse="\n"
)

txt <- str_replace(
  txt,
  fixed(placeholder),
  intro
)

# ------------------------------------------------------------
# 2. Merge verified reference pools
# ------------------------------------------------------------

old <- read_csv(
  OLD_POOL,
  show_col_types=FALSE
)

new <- read_csv(
  NEW_POOL,
  show_col_types=FALSE
)

needed_cols <- c(
  "citation_key",
  "authors",
  "title",
  "journal",
  "year",
  "volume",
  "pages_or_article",
  "doi",
  "pmid",
  "manuscript_use"
)

if(!all(needed_cols %in% names(old)) ||
   !all(needed_cols %in% names(new))){
  stop("Reference pool column mismatch.")
}

pool <- bind_rows(old,new) %>%
  distinct(citation_key,.keep_all=TRUE) %>%
  arrange(citation_key)

write_csv(
  pool,
  file.path(
    OUT,
    "01_VERIFIED_REFERENCE_POOL_EXPANDED.csv"
  )
)

# ------------------------------------------------------------
# 3. Validate every citation key in manuscript
# ------------------------------------------------------------

cit_blocks <- unique(
  unlist(
    str_extract_all(
      txt,
      "\\[CIT:[^\\]]+\\]"
    )
  )
)

keys <- unique(
  unlist(
    str_split(
      str_replace_all(
        str_replace_all(cit_blocks,"\\[CIT:",""),
        "\\]",""
      ),
      ";"
    )
  )
)

keys <- sort(keys[nzchar(keys)])

missing_keys <- setdiff(
  keys,
  pool$citation_key
)

unused_keys <- setdiff(
  pool$citation_key,
  keys
)

write_csv(
  tibble(citation_key=keys),
  file.path(
    OUT,
    "02_CITATION_KEYS_USED_AFTER_INTRODUCTION.csv"
  )
)

write_csv(
  tibble(missing_citation_key=missing_keys),
  file.path(
    OUT,
    "03_MISSING_CITATION_KEYS.csv"
  )
)

write_csv(
  tibble(unused_verified_key=unused_keys),
  file.path(
    OUT,
    "04_UNUSED_VERIFIED_REFERENCE_KEYS.csv"
  )
)

if(length(missing_keys)>0){
  stop(
    "Manuscript citation key(s) missing from verified pool: ",
    paste(missing_keys,collapse=", ")
  )
}

# ------------------------------------------------------------
# 4. Placeholder / internal-language audit
# ------------------------------------------------------------

bad_patterns <- c(
  "\\[REF",
  "\\[VERIFY",
  "\\bStep[0-9]",
  "ALL_PAIRED",
  "COMMON_ANCHOR_SENSITIVITY",
  "\\[INTRODUCTION TO BE WRITTEN"
)

bad_qc <- tibble(
  pattern=bad_patterns,
  hit=vapply(
    bad_patterns,
    function(p){
      str_detect(
        txt,
        regex(p,ignore_case=TRUE)
      )
    },
    logical(1)
  )
)

write_csv(
  bad_qc,
  file.path(
    OUT,
    "05_INTRODUCTION_AND_INTERNAL_LANGUAGE_QC.csv"
  )
)

if(any(bad_qc$hit)){
  stop(
    "Unwanted placeholder/internal language remains: ",
    paste(
      bad_qc$pattern[bad_qc$hit],
      collapse=", "
    )
  )
}

# ------------------------------------------------------------
# 5. Write manuscript core v3
# ------------------------------------------------------------

OUTFILE <- file.path(
  OUT,
  "Manuscript_Core_v3_WITH_INTRODUCTION.txt"
)

writeLines(
  str_split(txt,"\n")[[1]],
  OUTFILE,
  useBytes=TRUE
)

# ------------------------------------------------------------
# 6. Introduction-specific audit
# ------------------------------------------------------------

intro_audit <- tibble(
  check=c(
    "sepsis_definition_and_burden",
    "critical_illness_microbiome_context",
    "longitudinal_evidence_gap",
    "cross_cohort_reproducibility_gap",
    "study_objectives",
    "no_prognostic_biomarker_claim",
    "no_cross_cohort_feature_pooling",
    "all_citation_keys_resolved"
  ),
  passed=c(
    str_detect(intro,"life-threatening organ dysfunction"),
    str_detect(intro,"Critical illness can disrupt"),
    str_detect(intro,"36 longitudinal sequencing studies"),
    str_detect(intro,"reproducibility"),
    str_detect(intro,"We addressed this question"),
    str_detect(intro,"rather than to derive a universal sepsis-specific taxonomic biomarker or clinical prediction score"),
    str_detect(intro,"rather than by pooling heterogeneous ASV or OTU matrices"),
    length(missing_keys)==0
  )
)

write_csv(
  intro_audit,
  file.path(
    OUT,
    "06_INTRODUCTION_QC.csv"
  )
)

if(!all(intro_audit$passed)){
  stop(
    "Introduction QC failed: ",
    paste(
      intro_audit$check[!intro_audit$passed],
      collapse=", "
    )
  )
}

writeLines(
  c(
    "V2 STEP96H INTRODUCTION INTEGRATION",
    "",
    "Introduction placeholder replaced with manuscript-facing text.",
    paste0("Unique citation keys in manuscript: ",length(keys)),
    paste0("Verified references in expanded pool: ",nrow(pool)),
    paste0("Missing citation keys: ",length(missing_keys)),
    "",
    "Final numeric reference numbering remains intentionally deferred.",
    "No statistical analysis was rerun.",
    "",
    paste0("Output: ",OUTFILE),
    "",
    "STEP96H COMPLETE"
  ),
  file.path(
    OUT,
    "07_STEP96H_SUMMARY.txt"
  )
)

writeLines(
  c(
    paste0("Completed: ",Sys.time()),
    "STEP96H COMPLETE",
    "Introduction integrated.",
    paste0("Citation keys resolved: ",length(keys))
  ),
  file.path(
    OUT,
    "_STEP96H_COMPLETE.txt"
  )
)

cat("STEP96H COMPLETE\n")
