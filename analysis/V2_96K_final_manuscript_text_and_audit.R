
# ============================================================
# V2_96K FINAL MANUSCRIPT TEXT + SUBMISSION-FACING AUDIT
#
# Purpose:
# 1. Insert final Title, structured Abstract, and Keywords.
# 2. Remove residual internal workflow language.
# 3. Standardize Results subheadings.
# 4. Preserve the already-finalized 19-reference numbering.
# 5. Stage final main figures / Table 1 / new supplementary figure.
# 6. Run a final scientific-content and manuscript-language audit.
#
# NO statistical analysis is rerun.
# NO inferential result is changed.
# ============================================================

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(stringr)
  library(tibble)
})

ROOT <- "E:/sepsis_project"
RESULTS <- file.path(ROOT,"results")

MANUSCRIPT_IN <- file.path(
  RESULTS,
  "V2_96I2_REFERENCE_FORMAT_POLISH",
  "Manuscript_Core_v4_1_REFERENCES_POLISHED.txt"
)

REF_MAP <- file.path(
  RESULTS,
  "V2_96I2_REFERENCE_FORMAT_POLISH",
  "02_FINAL_REFERENCE_NUMBER_MAP_UNCHANGED.csv"
)

FIG96F <- file.path(
  RESULTS,
  "V2_96F_MANUSCRIPT_CORE_ASSEMBLY_AND_GAP_AUDIT",
  "02_FIGURE_ASSETS"
)

FIG3_FINAL_DIR <- file.path(
  RESULTS,
  "V2_96J3_FINAL_FIGURE3_SELF_CONTAINED"
)

TABLE1_IN <- file.path(
  RESULTS,
  "V2_96F_MANUSCRIPT_CORE_ASSEMBLY_AND_GAP_AUDIT",
  "Main_Table_1_Cohort_Characteristics_MANUSCRIPT.csv"
)

OUT <- file.path(
  RESULTS,
  "V2_96K_FINAL_MANUSCRIPT_TEXT_AND_AUDIT"
)

ASSET_OUT <- file.path(OUT,"02_FINAL_MAIN_FIGURES")
SUPP_OUT <- file.path(OUT,"03_NEW_SUPPLEMENTARY_FIGURE")

dir.create(OUT,recursive=TRUE,showWarnings=FALSE)
dir.create(ASSET_OUT,recursive=TRUE,showWarnings=FALSE)
dir.create(SUPP_OUT,recursive=TRUE,showWarnings=FALSE)

for(f in c(MANUSCRIPT_IN,REF_MAP,TABLE1_IN)){
  if(!file.exists(f)){
    stop("Missing Step96K input: ",f)
  }
}

txt <- paste(
  readLines(
    MANUSCRIPT_IN,
    warn=FALSE,
    encoding="UTF-8"
  ),
  collapse="\n"
)

# ------------------------------------------------------------
# 1. Final title / structured abstract / keywords
# ------------------------------------------------------------

TITLE <- paste(
  "Concordant longitudinal ecological displacement accompanies",
  "heterogeneous taxonomic trajectories across gut microbiome",
  "cohorts in sepsis and critical illness"
)

ABSTRACT <- paste(
  c(
    "Abstract",
    "",
    paste0(
      "Background: Gut-microbiome disruption is repeatedly observed during sepsis and critical illness, but the longitudinal reproducibility of ecosystem-level change relative to taxon-specific change remains unclear. ",
      "We tested whether within-patient ecological displacement shows a concordant cross-cohort direction despite heterogeneous taxonomic trajectories."
    ),
    "",
    paste0(
      "Results: We analyzed ten public 16S rRNA gene amplicon sequencing cohorts with prespecified analytical roles while retaining cohort-specific feature matrices. ",
      "Three cohorts showed FDR-supported progressive early-to-late Bray-Curtis displacement (paired n=9, 14, and 15; standardized paired effects dz=0.765, 1.095, and 0.623), with an additional ICU-infection cohort providing supportive evidence. ",
      "A non-sepsis surgical longitudinal control instead showed recovery-direction change and significant re-convergence in the common-anchor sensitivity analysis. ",
      "In the core repeated-sepsis cohort, median displacement increased from 0.725 at Day 3 to 0.904 at Day 7, while patient-specific signed genus trajectories were heterogeneous (mean cosine similarity=0.162; 35.6% of patient-pair comparisons <=0). ",
      "Across the three primary progressive cohorts, ecosystem-level direction was concordant, whereas taxonomic effect-vector concordance was lower (median pairwise Spearman rho=0.36 at genus level and 0.08 at family level); only 36.0% of shared genera and 31.8% of shared families changed in the same direction across all three cohorts. ",
      "In an external Control-Trauma-Sepsis cohort, sepsis samples were farther from the healthy-control centroid than trauma samples (0.860 vs 0.657; FDR=0.00132), with compositional sensitivity analyses supporting the separation. ",
      "Pulmonary versus recorded non-pulmonary infection source did not clearly modify longitudinal displacement (source-by-time beta=0.017, 95% CI -0.024 to 0.058; p=0.426)."
    ),
    "",
    paste0(
      "Conclusions: Across heterogeneous public cohorts, a concordant direction of longitudinal ecological displacement was accompanied by substantially more cohort-dependent taxonomic change. ",
      "These findings support ecological-state and trajectory-based approaches for studying microbiome disruption in sepsis and critical illness, while arguing against a universal sepsis-specific taxonomic signature. ",
      "Because progressive displacement also occurred in ICU-background cohorts, the ecological pattern should not be interpreted as uniquely sepsis-specific and requires prospective validation."
    )
  ),
  collapse="\n"
)

KEYWORDS <- paste0(
  "Keywords: Sepsis; Critical illness; Gut microbiome; Longitudinal analysis; ",
  "Ecological displacement; Microbiome trajectories; Taxonomic reproducibility; ",
  "16S rRNA gene amplicon sequencing"
)

for(ph in c(
  "[TITLE TO BE FINALIZED]",
  "[ABSTRACT TO BE WRITTEN AFTER FULL-TEXT FREEZE]",
  "[KEYWORDS TO BE FINALIZED]"
)){
  if(!str_detect(txt,fixed(ph))){
    stop("Expected manuscript placeholder missing: ",ph)
  }
}

txt <- str_replace(
  txt,
  fixed("[TITLE TO BE FINALIZED]"),
  TITLE
)

txt <- str_replace(
  txt,
  fixed("[ABSTRACT TO BE WRITTEN AFTER FULL-TEXT FREEZE]"),
  ABSTRACT
)

txt <- str_replace(
  txt,
  fixed("[KEYWORDS TO BE FINALIZED]"),
  KEYWORDS
)

# ------------------------------------------------------------
# 2. Remove manuscript-internal workflow labels
# ------------------------------------------------------------

txt <- str_replace(
  txt,
  "^MANUSCRIPT CORE DRAFT\\n\\n",
  ""
)

txt <- str_replace_all(
  txt,
  fixed("DISCUSSION DRAFT V1 — EVIDENCE-GROUNDED WITH LITERATURE PLACEHOLDERS\n\n"),
  ""
)

txt <- str_replace_all(
  txt,
  fixed("The V2 analysis integrated"),
  "The analysis integrated"
)

txt <- str_replace_all(
  txt,
  fixed("Taken together, the V2 results support"),
  "Taken together, these results support"
)

txt <- str_replace_all(
  txt,
  fixed("classified in the frozen evidence synthesis as showing robust progressive ecological displacement"),
  "showing robust progressive ecological displacement in the prespecified cross-cohort evidence synthesis"
)

txt <- str_replace_all(
  txt,
  fixed("manuscript-level synthesis was performed using frozen cohort-specific ecological or taxonomic effect estimates"),
  "manuscript-level synthesis was performed using finalized cohort-specific ecological or taxonomic effect estimates"
)

txt <- str_replace_all(
  txt,
  fixed("entered the canonical workspace with previously validated ASV tables"),
  "entered the harmonized analysis workflow with previously validated ASV tables"
)

txt <- str_replace_all(
  txt,
  fixed("underwent the project workflow for cohort-specific paired-end DADA2 processing"),
  "underwent cohort-specific paired-end DADA2 processing within the analysis workflow"
)

txt <- str_replace_all(
  txt,
  fixed("across the eight primary/supporting longitudinal datasets"),
  "across the eight longitudinal datasets represented in the harmonized analysis set"
)

# ------------------------------------------------------------
# 3. Standardize Results subheadings: journal-facing, unnumbered
# ------------------------------------------------------------

heading_map <- c(
  "1. Cohort architecture and analytical framework" =
    "Cohort architecture and analytical framework",
  "2. Progressive ecological displacement across longitudinal critical-illness cohorts" =
    "Progressive ecological displacement across longitudinal critical-illness cohorts",
  "3. Heterogeneous taxonomic routes accompany ecological displacement in core sepsis" =
    "Heterogeneous taxonomic routes accompany ecological displacement in core sepsis",
  "4. External ecological-state differentiation beyond severe trauma" =
    "External ecological-state differentiation beyond severe trauma",
  "5. Infection source does not clearly modify longitudinal ecological displacement" =
    "Infection source does not clearly modify longitudinal ecological displacement"
)

for(old in names(heading_map)){
  txt <- str_replace_all(
    txt,
    fixed(old),
    heading_map[[old]]
  )
}

# ------------------------------------------------------------
# 4. Reference numbering must remain exactly unchanged
# ------------------------------------------------------------

ref_map <- read_csv(
  REF_MAP,
  show_col_types=FALSE
)

expected_ref_n <- seq_len(nrow(ref_map))

if(nrow(ref_map)!=19){
  stop("Expected 19 finalized references, found ",nrow(ref_map))
}

if(!identical(as.integer(ref_map$reference_number),expected_ref_n)){
  stop("Reference map is not contiguous 1-19.")
}

numeric_citations <- str_extract_all(
  txt,
  "\\[[0-9][0-9,\\- ]*\\]"
)[[1]]

extract_nums <- function(x){
  raw <- str_remove_all(x,"\\[|\\]")
  pieces <- unlist(str_split(raw,",\\s*"))
  nums <- integer()

  for(p in pieces){
    p <- str_trim(p)

    if(str_detect(p,"^[0-9]+-[0-9]+$")){
      z <- as.integer(str_split(p,"-")[[1]])
      nums <- c(nums,seq(z[1],z[2]))
    }else if(str_detect(p,"^[0-9]+$")){
      nums <- c(nums,as.integer(p))
    }
  }

  nums
}

cited_nums <- sort(unique(unlist(lapply(numeric_citations,extract_nums))))

if(!identical(cited_nums,expected_ref_n)){
  stop(
    "Numeric citation coverage changed after final text insertion. Found: ",
    paste(cited_nums,collapse=",")
  )
}

# ------------------------------------------------------------
# 5. Abstract numerical lock
# ------------------------------------------------------------

abstract_checks <- c(
  "dz=0.765",
  "1.095",
  "0.623",
  "0.725",
  "0.904",
  "0.162",
  "35.6%",
  "rho=0.36",
  "0.08",
  "36.0%",
  "31.8%",
  "0.860",
  "0.657",
  "FDR=0.00132",
  "beta=0.017",
  "95% CI -0.024 to 0.058",
  "p=0.426"
)

abstract_numeric_qc <- tibble(
  expected_string=abstract_checks,
  present=vapply(
    abstract_checks,
    function(x) str_detect(ABSTRACT,fixed(x)),
    logical(1)
  )
)

write_csv(
  abstract_numeric_qc,
  file.path(
    OUT,
    "01_ABSTRACT_NUMERIC_LOCK_QC.csv"
  )
)

if(!all(abstract_numeric_qc$present)){
  stop("Abstract numeric lock QC failed.")
}

# ------------------------------------------------------------
# 6. Internal-language / placeholder / terminology audit
# ------------------------------------------------------------

bad_patterns <- c(
  "MANUSCRIPT CORE DRAFT",
  "DISCUSSION DRAFT",
  "\\[TITLE TO BE",
  "\\[ABSTRACT TO BE",
  "\\[KEYWORDS TO BE",
  "\\[REF",
  "\\[CIT:",
  "\\[VERIFY",
  "\\bStep[0-9]",
  "\\bV2 analysis\\b",
  "\\bV2 results\\b",
  "frozen evidence synthesis",
  "frozen cohort-specific",
  "canonical workspace",
  "project workflow",
  "metagenomic sequencing",
  "\\b16S sequencing\\b"
)

language_qc <- tibble(
  pattern=bad_patterns,
  hit=vapply(
    bad_patterns,
    function(p) str_detect(txt,regex(p,ignore_case=TRUE)),
    logical(1)
  )
)

write_csv(
  language_qc,
  file.path(
    OUT,
    "02_FINAL_INTERNAL_LANGUAGE_QC.csv"
  )
)

if(any(language_qc$hit)){
  stop(
    "Residual internal/undesired manuscript language: ",
    paste(language_qc$pattern[language_qc$hit],collapse=", ")
  )
}

# ------------------------------------------------------------
# 7. Scientific guardrails audit
# ------------------------------------------------------------

guardrails <- tibble(
  claim=c(
    "not_sepsis_specific",
    "no_universal_taxonomic_signature",
    "no_cross_cohort_ASV_OTU_pooling",
    "EII_not_clinical_prediction",
    "infection_source_not_claimed_equivalent",
    "common_anchor_not_independent_replication",
    "taxonomic_vs_ecological_not_formal_commensurable_test"
  ),
  evidence_present=c(
    str_detect(txt,"should not be interpreted as uniquely sepsis-specific"),
    str_detect(txt,"universal sepsis-specific taxonomic signature"),
    str_detect(txt,"ASV or OTU abundance matrices were not pooled across studies"),
    str_detect(txt,"not developed or validated as a clinical prediction or prognostic score"),
    str_detect(txt,"cannot establish mechanistic equivalence between source groups"),
    str_detect(txt,"should not be interpreted as independent replication"),
    str_detect(txt,"rather than as a formal statistical comparison of directly commensurable metrics")
  )
)

write_csv(
  guardrails,
  file.path(
    OUT,
    "03_SCIENTIFIC_GUARDRAIL_QC.csv"
  )
)

if(!all(guardrails$evidence_present)){
  stop(
    "Scientific guardrail missing: ",
    paste(guardrails$claim[!guardrails$evidence_present],collapse=", ")
  )
}

# ------------------------------------------------------------
# 8. Stage final manuscript text
# ------------------------------------------------------------

MANUSCRIPT_OUT <- file.path(
  OUT,
  "Manuscript_v5_SUBMISSION_DRAFT.txt"
)

writeLines(
  str_split(txt,"\n")[[1]],
  MANUSCRIPT_OUT,
  useBytes=TRUE
)

write_csv(
  ref_map,
  file.path(
    OUT,
    "04_FINAL_REFERENCE_NUMBER_MAP.csv"
  )
)

file.copy(
  TABLE1_IN,
  file.path(
    OUT,
    "Main_Table_1_Cohort_Characteristics_MANUSCRIPT.csv"
  ),
  overwrite=TRUE
)

# ------------------------------------------------------------
# 9. Stage Figure 1,2,4,5 from Step96F and final Figure3 from J3
# ------------------------------------------------------------

copy_first_match <- function(
  dir,
  pattern,
  dest_name,
  target_dir
){
  cand <- list.files(
    dir,
    recursive=TRUE,
    full.names=TRUE,
    pattern=pattern,
    ignore.case=TRUE
  )

  if(length(cand)<1){
    stop(
      "Could not find figure asset: ",
      pattern,
      " under ",
      dir
    )
  }

  # Prefer PDF for PDF destination and PNG for PNG destination.
  ext <- tools::file_ext(dest_name)
  same_ext <- cand[
    tolower(tools::file_ext(cand))==tolower(ext)
  ]

  if(length(same_ext)>=1){
    src <- same_ext[[1]]
  }else{
    src <- cand[[1]]
  }

  dest <- file.path(
    target_dir,
    dest_name
  )

  file.copy(
    src,
    dest,
    overwrite=TRUE
  )

  src
}

asset_rows <- list()

for(fig in c(1,2,4,5)){

  pdf_src <- copy_first_match(
    FIG96F,
    paste0(
      "FROZEN_Figure_",fig,
      ".*\\.pdf$"
    ),
    paste0(
      "Figure_",fig,"_FINAL.pdf"
    ),
    ASSET_OUT
  )

  png_src <- copy_first_match(
    FIG96F,
    paste0(
      "FROZEN_Figure_",fig,
      ".*\\.png$"
    ),
    paste0(
      "Figure_",fig,"_FINAL.png"
    ),
    ASSET_OUT
  )

  asset_rows[[length(asset_rows)+1]] <- tibble(
    figure=paste0("Figure ",fig),
    pdf_source=pdf_src,
    png_source=png_src,
    role="MAIN"
  )
}

fig3_pdf <- file.path(
  FIG3_FINAL_DIR,
  "V2_96J2_Figure_3_FINAL_ABCD.pdf"
)

fig3_png <- file.path(
  FIG3_FINAL_DIR,
  "V2_96J2_Figure_3_FINAL_ABCD_600dpi.png"
)

fig3_legend <- file.path(
  FIG3_FINAL_DIR,
  "V2_96J2_Figure_3_FINAL_LEGEND.txt"
)

for(f in c(fig3_pdf,fig3_png,fig3_legend)){
  if(!file.exists(f)){
    stop("Missing final Step96J3 Figure3 asset: ",f)
  }
}

file.copy(
  fig3_pdf,
  file.path(
    ASSET_OUT,
    "Figure_3_FINAL_ABCD.pdf"
  ),
  overwrite=TRUE
)

file.copy(
  fig3_png,
  file.path(
    ASSET_OUT,
    "Figure_3_FINAL_ABCD_600dpi.png"
  ),
  overwrite=TRUE
)

file.copy(
  fig3_legend,
  file.path(
    ASSET_OUT,
    "Figure_3_FINAL_LEGEND.txt"
  ),
  overwrite=TRUE
)

asset_rows[[length(asset_rows)+1]] <- tibble(
  figure="Figure 3",
  pdf_source=fig3_pdf,
  png_source=fig3_png,
  role="MAIN_FINAL_ABCD"
)

# New supplementary concordance sensitivity figure
supp_pdf <- file.path(
  FIG3_FINAL_DIR,
  "V2_96J2_Supplementary_Taxonomic_Concordance_Sensitivity.pdf"
)

supp_png <- file.path(
  FIG3_FINAL_DIR,
  "V2_96J2_Supplementary_Taxonomic_Concordance_Sensitivity.png"
)

supp_legend <- file.path(
  FIG3_FINAL_DIR,
  "V2_96J2_Supplementary_Taxonomic_Concordance_LEGEND.txt"
)

for(f in c(supp_pdf,supp_png,supp_legend)){
  if(!file.exists(f)){
    stop("Missing Step96J3 supplementary asset: ",f)
  }
}

file.copy(
  supp_pdf,
  file.path(
    SUPP_OUT,
    "Supplementary_Figure_Taxonomic_Concordance_Sensitivity.pdf"
  ),
  overwrite=TRUE
)

file.copy(
  supp_png,
  file.path(
    SUPP_OUT,
    "Supplementary_Figure_Taxonomic_Concordance_Sensitivity.png"
  ),
  overwrite=TRUE
)

file.copy(
  supp_legend,
  file.path(
    SUPP_OUT,
    "Supplementary_Figure_Taxonomic_Concordance_Sensitivity_LEGEND.txt"
  ),
  overwrite=TRUE
)

asset_registry <- bind_rows(asset_rows)

write_csv(
  asset_registry,
  file.path(
    OUT,
    "05_FINAL_MAIN_FIGURE_ASSET_REGISTRY.csv"
  )
)

# ------------------------------------------------------------
# 10. Text metrics / abstract length
# ------------------------------------------------------------

word_count <- function(x){
  x <- str_replace_all(x,"\\s+"," ")
  length(
    unlist(
      str_split(
        str_trim(x),
        "\\s+"
      )
    )
  )
}

text_metrics <- tibble(
  component=c(
    "Title",
    "Abstract",
    "Full manuscript including references"
  ),
  word_count=c(
    word_count(TITLE),
    word_count(ABSTRACT),
    word_count(txt)
  ),
  character_count=c(
    nchar(TITLE),
    nchar(ABSTRACT),
    nchar(txt)
  )
)

write_csv(
  text_metrics,
  file.path(
    OUT,
    "06_TEXT_LENGTH_METRICS.csv"
  )
)

# ------------------------------------------------------------
# 11. Final submission-readiness audit
# ------------------------------------------------------------

readiness <- tibble(
  item=c(
    "Title",
    "Structured abstract",
    "Keywords",
    "Introduction",
    "Methods",
    "Results",
    "Discussion",
    "References 1-19",
    "Table 1",
    "Main Figure 1",
    "Main Figure 2",
    "Main Figure 3 A-D",
    "Main Figure 4",
    "Main Figure 5",
    "New supplementary taxonomic-concordance figure",
    "Scientific/internal-language QC",
    "Author names and affiliations",
    "Author contributions",
    "Funding",
    "Acknowledgements",
    "Competing interests",
    "Data/code availability statement",
    "Ethics/public-data statement",
    "Journal submission metadata"
  ),
  status=c(
    rep("READY",16),
    rep("AUTHOR_LEVEL_INFORMATION_STILL_REQUIRED",8)
  )
)

write_csv(
  readiness,
  file.path(
    OUT,
    "07_SUBMISSION_READINESS_AUDIT.csv"
  )
)

# ------------------------------------------------------------
# 12. Summary
# ------------------------------------------------------------

writeLines(
  c(
    "V2 STEP96K FINAL MANUSCRIPT TEXT AND AUDIT",
    "",
    paste0("Title: ",TITLE),
    "",
    paste0(
      "Abstract word count: ",
      word_count(ABSTRACT)
    ),
    "Final scientific references: 19, numbering preserved.",
    "Main Figures 1-5 staged; Figure 3 is final A-D composite.",
    "New supplementary taxonomic-concordance sensitivity figure staged.",
    "",
    "Residual REF/CIT/VERIFY/Step labels: 0.",
    "Residual V2/DRAFT/frozen-workflow manuscript language: 0.",
    "Scientific guardrails: PASS.",
    "",
    "Science-facing manuscript text is now submission-draft ready.",
    "Only author-level/declaration/submission metadata remain to be supplied separately.",
    "",
    "NO statistical analysis rerun.",
    "NO inferential result changed.",
    "",
    "STEP96K COMPLETE"
  ),
  file.path(
    OUT,
    "08_STEP96K_SUMMARY.txt"
  )
)

writeLines(
  c(
    paste0("Completed: ",Sys.time()),
    "STEP96K COMPLETE",
    "Final manuscript text assembled and audited.",
    "Scientific manuscript content READY."
  ),
  file.path(
    OUT,
    "_STEP96K_COMPLETE.txt"
  )
)

cat("STEP96K COMPLETE\n")
cat("Abstract words:",word_count(ABSTRACT),"\n")
