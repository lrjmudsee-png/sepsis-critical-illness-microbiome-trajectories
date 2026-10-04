
# ============================================================
# V2_96F MANUSCRIPT CORE ASSEMBLY + GAP AUDIT
#
# Purpose:
# Assemble the manuscript-facing Methods, Results, Discussion,
# Table 1, frozen main figures and new Figure 3D assets into a
# single reproducible manuscript-core package.
#
# NO statistical model is rerun.
# NO frozen result is overwritten.
# ============================================================

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(stringr)
  library(tibble)
})

ROOT <- "E:/sepsis_project"
RESULTS <- file.path(ROOT,"results")

METHODS <- file.path(
  RESULTS,
  "V2_96D4_METHODS_FINAL_MANUSCRIPT_POLISH",
  "Methods_v3_MANUSCRIPT_READY.txt"
)

RESULTS_TXT <- file.path(
  RESULTS,
  "V2_96B_RESULTS_INTEGRATION",
  "Results_Draft_v6_WITH_CROSS_COHORT_REPRODUCIBILITY.txt"
)

DISCUSSION <- file.path(
  RESULTS,
  "V2_96C_DISCUSSION_INTEGRATION",
  "Discussion_Draft_v2_WITH_CROSS_COHORT_REPRODUCIBILITY.txt"
)

TABLE1 <- file.path(
  RESULTS,
  "V2_35F3_RESULTS_AND_TABLE1_FINAL_FREEZE_SENTENCE_AUDIT",
  "01_FINAL_TABLES",
  "Main_Table_1_Cohort_Characteristics_MANUSCRIPT.csv"
)

TABLE1_FOOT <- file.path(
  RESULTS,
  "V2_35F3_RESULTS_AND_TABLE1_FINAL_FREEZE_SENTENCE_AUDIT",
  "01_FINAL_TABLES"
)

FIG_DIR <- file.path(
  RESULTS,
  "V2_35D3_FINAL_FIGURE_FORMATTING_FIX"
)

FIG3D_DIR <- file.path(
  RESULTS,
  "V2_96E_FIGURE3D_CROSS_COHORT_REPRODUCIBILITY_FINAL"
)

DISC_LIT_MAP <- file.path(
  RESULTS,
  "V2_96C_DISCUSSION_INTEGRATION",
  "02_VERIFIED_LITERATURE_MAP.csv"
)

OUT <- file.path(
  RESULTS,
  "V2_96F_MANUSCRIPT_CORE_ASSEMBLY_AND_GAP_AUDIT"
)

ASSET_DIR <- file.path(OUT,"02_FIGURE_ASSETS")

dir.create(OUT,recursive=TRUE,showWarnings=FALSE)
dir.create(ASSET_DIR,recursive=TRUE,showWarnings=FALSE)

critical <- c(METHODS,RESULTS_TXT,DISCUSSION,TABLE1)

if(any(!file.exists(critical))){
  stop(
    "Missing manuscript input(s): ",
    paste(critical[!file.exists(critical)],collapse=" ; ")
  )
}

# ------------------------------------------------------------
# 1. Read manuscript sections
# ------------------------------------------------------------
methods_lines <- readLines(METHODS,warn=FALSE,encoding="UTF-8")
results_lines <- readLines(RESULTS_TXT,warn=FALSE,encoding="UTF-8")
discussion_lines <- readLines(DISCUSSION,warn=FALSE,encoding="UTF-8")

# Remove redundant first-line "Methods" from section content if present.
if(length(methods_lines)>0 && trimws(tolower(methods_lines[1]))=="methods"){
  methods_lines <- methods_lines[-1]
}

assembled <- c(
  "MANUSCRIPT CORE DRAFT",
  "",
  "[TITLE TO BE FINALIZED]",
  "",
  "[ABSTRACT TO BE WRITTEN AFTER FULL-TEXT FREEZE]",
  "",
  "[KEYWORDS TO BE FINALIZED]",
  "",
  "Introduction",
  "",
  "[INTRODUCTION TO BE WRITTEN AFTER FINAL LITERATURE INTEGRATION]",
  "",
  "Methods",
  "",
  methods_lines,
  "",
  "Results",
  "",
  results_lines,
  "",
  "Discussion",
  "",
  discussion_lines,
  "",
  "References",
  "",
  "[REFERENCE LIST TO BE ASSEMBLED AFTER ALL REF PLACEHOLDERS ARE RESOLVED]"
)

MANUSCRIPT <- file.path(
  OUT,
  "01_MANUSCRIPT_CORE_DRAFT_v1.txt"
)

writeLines(
  assembled,
  MANUSCRIPT,
  useBytes=TRUE
)

# ------------------------------------------------------------
# 2. Table 1
# ------------------------------------------------------------
file.copy(
  TABLE1,
  file.path(
    OUT,
    "Main_Table_1_Cohort_Characteristics_MANUSCRIPT.csv"
  ),
  overwrite=TRUE
)

# Copy likely Table1 footnote companion(s).
if(dir.exists(TABLE1_FOOT)){

  foot_candidates <- list.files(
    TABLE1_FOOT,
    pattern="foot|note",
    full.names=TRUE,
    ignore.case=TRUE
  )

  if(length(foot_candidates)>0){
    for(f in foot_candidates){
      file.copy(
        f,
        file.path(OUT,basename(f)),
        overwrite=TRUE
      )
    }
  }
}

# ------------------------------------------------------------
# 3. Main frozen figure assets
# ------------------------------------------------------------
figure_registry <- tibble()

if(dir.exists(FIG_DIR)){

  old_figs <- list.files(
    FIG_DIR,
    pattern="\\.(pdf|png)$",
    recursive=TRUE,
    full.names=TRUE,
    ignore.case=TRUE
  )

  # Keep only files that look like manuscript Figure 1-5.
  old_figs <- old_figs[
    str_detect(
      basename(old_figs),
      regex("figure[_ -]?[1-5]|^fig[_ -]?[1-5]",ignore_case=TRUE)
    )
  ]

  if(length(old_figs)>0){

    for(f in old_figs){

      dest <- file.path(
        ASSET_DIR,
        paste0("FROZEN_",basename(f))
      )

      file.copy(f,dest,overwrite=TRUE)

      figure_registry <- bind_rows(
        figure_registry,
        tibble(
          asset=basename(dest),
          source=f,
          role="FROZEN_MAIN_FIGURE_ASSET",
          final_composite_status=ifelse(
            str_detect(
              basename(f),
              regex("figure[_ -]?3|^fig[_ -]?3",ignore_case=TRUE)
            ),
            "FIGURE3_BASE_REQUIRES_NEW_PANEL_D_INTEGRATION",
            "FROZEN"
          )
        )
      )
    }
  }
}

# ------------------------------------------------------------
# 4. New Figure 3D / supplementary assets
# ------------------------------------------------------------
new_assets <- c(
  "Figure3D_Cross_Cohort_Taxonomic_Concordance_MAIN.pdf",
  "Figure3D_Cross_Cohort_Taxonomic_Concordance_MAIN.png",
  "Supplementary_Figure_Cross_Cohort_Taxonomic_Concordance_Sensitivity.pdf",
  "Supplementary_Figure_Cross_Cohort_Taxonomic_Concordance_Sensitivity.png",
  "03_FIGURE3D_LEGEND.txt",
  "04_SUPPLEMENTARY_FIGURE_LEGEND.txt"
)

for(a in new_assets){

  src <- file.path(FIG3D_DIR,a)

  if(!file.exists(src)){
    stop("Missing Step96E figure asset: ",src)
  }

  dest <- file.path(ASSET_DIR,a)
  file.copy(src,dest,overwrite=TRUE)

  figure_registry <- bind_rows(
    figure_registry,
    tibble(
      asset=a,
      source=src,
      role=case_when(
        str_detect(a,"Figure3D_") ~ "NEW_MAIN_FIGURE3_PANEL_D",
        str_detect(a,"Supplementary_Figure") ~ "NEW_SUPPLEMENTARY_FIGURE",
        TRUE ~ "FIGURE_LEGEND"
      ),
      final_composite_status=case_when(
        str_detect(a,"Figure3D_") ~
          "READY_AS_PANEL_D; NEEDS_COMPOSITION_WITH_EXISTING_FIGURE3",
        str_detect(a,"Supplementary_Figure") ~
          "READY",
        TRUE ~
          "TEXT_READY"
      )
    )
  )
}

write_csv(
  figure_registry,
  file.path(
    OUT,
    "03_FIGURE_ASSET_REGISTRY.csv"
  )
)

# ------------------------------------------------------------
# 5. Literature map
# ------------------------------------------------------------
if(file.exists(DISC_LIT_MAP)){
  file.copy(
    DISC_LIT_MAP,
    file.path(
      OUT,
      "04_VERIFIED_DISCUSSION_LITERATURE_MAP.csv"
    ),
    overwrite=TRUE
  )
}

# ------------------------------------------------------------
# 6. Gap audit
# ------------------------------------------------------------
full_text <- paste(assembled,collapse="\n")

# Extract bracketed REF placeholders.
ref_tags <- unique(
  unlist(
    str_extract_all(
      full_text,
      "\\[(?:REF[^\\]]*)\\]"
    )
  )
)

ref_tags <- ref_tags[
  nzchar(ref_tags)
]

if(length(ref_tags)>0){
  write_csv(
    tibble(reference_placeholder=sort(ref_tags)),
    file.path(
      OUT,
      "05_UNRESOLVED_REFERENCE_PLACEHOLDERS.csv"
    )
  )
}else{
  write_csv(
    tibble(reference_placeholder=character()),
    file.path(
      OUT,
      "05_UNRESOLVED_REFERENCE_PLACEHOLDERS.csv"
    )
  )
}

# Detect internal workflow wording.
internal_patterns <- c(
  "\\[VERIFY",
  "\\bStep[0-9]",
  "ALL_PAIRED",
  "COMMON_ANCHOR_SENSITIVITY"
)

internal_hits <- tibble(
  pattern=internal_patterns,
  hit=vapply(
    internal_patterns,
    function(p) str_detect(full_text,regex(p,ignore_case=TRUE)),
    logical(1)
  )
)

write_csv(
  internal_hits,
  file.path(
    OUT,
    "06_INTERNAL_LANGUAGE_AUDIT.csv"
  )
)

gap_audit <- tibble(
  item=c(
    "Methods",
    "Results",
    "Discussion",
    "Table 1",
    "Main Figures 1-5 assets",
    "New Figure 3D panel",
    "Supplementary taxonomic-concordance sensitivity figure",
    "Introduction",
    "Abstract",
    "Title",
    "Keywords",
    "Reference placeholders",
    "Final numbered reference list",
    "Final Figure 3 composite",
    "Final full-manuscript consistency audit"
  ),
  status=c(
    "READY",
    "READY",
    "DRAFT_READY_BUT_REFERENCE_PLACEHOLDERS_REMAIN",
    "READY",
    ifelse(nrow(figure_registry)>0,"AVAILABLE","CHECK_REQUIRED"),
    "READY",
    "READY",
    "NOT_YET_WRITTEN",
    "NOT_YET_WRITTEN",
    "NOT_FINALIZED",
    "NOT_FINALIZED",
    ifelse(length(ref_tags)>0,paste0(length(ref_tags)," UNIQUE_PLACEHOLDERS"),"NONE"),
    "NOT_YET_ASSEMBLED",
    "PANEL_D_READY; COMPOSITE_PENDING",
    "PENDING"
  ),
  next_action=c(
    "Use Methods v3 in full manuscript.",
    "Use Results v6 in full manuscript.",
    "Resolve literature placeholders and polish complete Discussion.",
    "Use frozen Table 1.",
    "Retain frozen inferential content.",
    "Integrate with existing Figure 3 base.",
    "Include in Supplement.",
    "Write after literature/reference integration.",
    "Write last after main text stabilizes.",
    "Finalize after full manuscript.",
    "Finalize after abstract/title.",
    "Resolve using verified literature + targeted web review.",
    "Renumber references in order of first appearance.",
    "Compose existing Figure 3 with new panel D.",
    "Run sentence/numeric/figure/table/reference consistency audit."
  )
)

write_csv(
  gap_audit,
  file.path(
    OUT,
    "07_MANUSCRIPT_GAP_AUDIT.csv"
  )
)

# ------------------------------------------------------------
# 7. Next-stage decision
# ------------------------------------------------------------
decision <- c(
  "V2 STEP96F MANUSCRIPT CORE ASSEMBLY",
  "",
  "READY:",
  "- Methods v3",
  "- Results v6",
  "- Discussion v2 evidence structure",
  "- Table 1",
  "- frozen Figure 1-5 assets (when located)",
  "- Figure 3D panel",
  "- supplementary taxonomic-concordance sensitivity figure",
  "",
  "REMAINING BEFORE A COMPLETE SUBMISSION DRAFT:",
  paste0("- Resolve ",length(ref_tags)," unique reference placeholder(s)."),
  "- Write/finalize Introduction.",
  "- Assemble final numbered reference list.",
  "- Compose final Figure 3 with new panel D.",
  "- Write Abstract and finalize title/keywords.",
  "- Run final manuscript consistency audit.",
  "",
  "NO FURTHER STATISTICAL ANALYSIS IS REQUIRED.",
  "",
  "STEP96F COMPLETE"
)

writeLines(
  decision,
  file.path(
    OUT,
    "08_NEXT_STAGE_DECISION.txt"
  )
)

writeLines(
  c(
    paste0("Completed: ",Sys.time()),
    "STEP96F COMPLETE",
    "Manuscript core assembled.",
    paste0("Unique unresolved REF placeholders: ",length(ref_tags)),
    "No statistical inference rerun."
  ),
  file.path(
    OUT,
    "_STEP96F_COMPLETE.txt"
  )
)

cat("STEP96F COMPLETE\n")
cat("Unresolved reference placeholders:",length(ref_tags),"\n")
