
# ============================================================
# V2_36J MICROBIOME ENHANCEMENT EVIDENCE FREEZE
#
# Freeze two enhancement branches:
# A. PRJNA851469 Day-3 clinical outcome branch:
#    NO_CLEAR_CLINICAL_OUTCOME_ASSOCIATION
#    -> retain as negative exploratory sensitivity only.
#
# B. Cross-cohort taxonomic reproducibility branch:
#    ECOLOGICAL_DIRECTION_MORE_REPRODUCIBLE_THAN_TAXONOMIC_EFFECTS
#    -> integrate into main manuscript narrative.
#
# No new statistical inference.
# ============================================================

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
})

ROOT <- "E:/sepsis_project"
RESULTS <- file.path(ROOT,"results")

TAX_DIR <- file.path(
  RESULTS,
  "V2_36I3_CROSS_COHORT_TAXONOMIC_REPRODUCIBILITY_FINAL"
)

CLIN_DIR <- file.path(
  RESULTS,
  "V2_36G2_PRJNA851469_DAY3_CLINICAL_OUTCOME_MODELS_FIXED"
)

OUT <- file.path(
  RESULTS,
  "V2_36J_MICROBIOME_ENHANCEMENT_EVIDENCE_FREEZE"
)

dir.create(OUT,recursive=TRUE,showWarnings=FALSE)

required_tax <- c(
  "03_PRIMARY3_PAIRWISE_CONCORDANCE.csv",
  "06_THREEWAY_DIRECTION_SUMMARY.csv",
  "08_ECOLOGICAL_VS_TAXONOMIC_SYNTHESIS.csv",
  "10_EVIDENCE_DECISION.csv",
  "11_Figure_candidate_taxonomic_concordance.pdf",
  "11_Figure_candidate_taxonomic_concordance.png",
  "12_MANUSCRIPT_SAFE_SUMMARY.txt"
)

for(f in required_tax){
  p <- file.path(TAX_DIR,f)
  if(!file.exists(p)) stop("Missing Step36I3 source: ",p)
}

tax_dec <- read_csv(
  file.path(TAX_DIR,"10_EVIDENCE_DECISION.csv"),
  show_col_types=FALSE
)

tax_syn <- read_csv(
  file.path(TAX_DIR,"08_ECOLOGICAL_VS_TAXONOMIC_SYNTHESIS.csv"),
  show_col_types=FALSE
)

three <- read_csv(
  file.path(TAX_DIR,"06_THREEWAY_DIRECTION_SUMMARY.csv"),
  show_col_types=FALSE
)

if(
  nrow(tax_dec)!=1 ||
  tax_dec$evidence_tier[[1]] !=
    "ECOLOGICAL_DIRECTION_MORE_REPRODUCIBLE_THAN_TAXONOMIC_EFFECTS"
){
  stop("Step36I3 evidence tier is not the expected freeze tier.")
}

# Clinical branch is optional if local folder has been preserved.
clinical_status <- "NOT_FOUND"
clinical_primary_hr <- NA_real_
clinical_primary_p <- NA_real_
clinical_adjusted_hr <- NA_real_
clinical_adjusted_p <- NA_real_

clin_file <- file.path(CLIN_DIR,"08_EVIDENCE_DECISION.csv")

if(file.exists(clin_file)){

  clin <- read_csv(
    clin_file,
    show_col_types=FALSE
  )

  if(nrow(clin)>=1){

    clinical_status <- clin$evidence_tier[[1]]
    clinical_primary_hr <- clin$primary_HR[[1]]
    clinical_primary_p <- clin$primary_p[[1]]
    clinical_adjusted_hr <- clin$antibiotic_adjusted_HR[[1]]
    clinical_adjusted_p <- clin$antibiotic_adjusted_p[[1]]
  }
}

get_rank <- function(rk){
  tax_syn %>% filter(rank==rk)
}

g <- get_rank("GENUS")
f <- get_rank("FAMILY")

three_primary_g <- three %>%
  filter(
    analysis_role=="PRIMARY",
    rank=="GENUS"
  )

three_primary_f <- three %>%
  filter(
    analysis_role=="PRIMARY",
    rank=="FAMILY"
  )

freeze <- tibble(
  branch=c(
    "CROSS_COHORT_TAXONOMIC_REPRODUCIBILITY",
    "PRJNA851469_DAY3_CLINICAL_OUTCOME"
  ),

  role=c(
    "MAIN_MANUSCRIPT_ENHANCEMENT",
    "NEGATIVE_EXPLORATORY_SENSITIVITY_ONLY"
  ),

  evidence_tier=c(
    tax_dec$evidence_tier[[1]],
    clinical_status
  ),

  manuscript_action=c(
    "INTEGRATE_MAIN_TEXT_AND_MAIN_FIGURE_OR_MAIN_FIGURE_PANEL",
    "DO_NOT_PROMOTE_TO_MAIN_CLAIM"
  )
)

write_csv(
  freeze,
  file.path(
    OUT,
    "01_ENHANCEMENT_BRANCH_FREEZE.csv"
  )
)

metrics <- tibble(
  metric=c(
    "primary_ecological_direction_concordance",
    "genus_median_pairwise_spearman_rho",
    "genus_median_sign_agreement",
    "genus_median_top_jaccard",
    "genus_threeway_same_direction_proportion",
    "family_median_pairwise_spearman_rho",
    "family_median_sign_agreement",
    "family_median_top_jaccard",
    "family_threeway_same_direction_proportion",
    "clinical_primary_HR",
    "clinical_primary_p",
    "clinical_antibiotic_adjusted_HR",
    "clinical_antibiotic_adjusted_p"
  ),
  value=c(
    tax_dec$primary_ecological_direction_concordance[[1]],
    g$median_spearman_rho[[1]],
    g$median_sign_agreement[[1]],
    g$median_top_jaccard[[1]],
    three_primary_g$proportion_same_direction_all3[[1]],
    f$median_spearman_rho[[1]],
    f$median_sign_agreement[[1]],
    f$median_top_jaccard[[1]],
    three_primary_f$proportion_same_direction_all3[[1]],
    clinical_primary_hr,
    clinical_primary_p,
    clinical_adjusted_hr,
    clinical_adjusted_p
  )
)

write_csv(
  metrics,
  file.path(
    OUT,
    "02_FROZEN_ENHANCEMENT_METRICS.csv"
  )
)

# Frozen manuscript wording
wording <- c(
  "FROZEN MAIN-TEXT WORDING — CROSS-COHORT TAXONOMIC REPRODUCIBILITY",
  "",
  "Across the three cohorts showing robust progressive ecological displacement, the direction of ecosystem-level change was concordant, whereas the underlying taxonomic effects were substantially less reproducible across cohorts.",
  "",
  paste0(
    "At the genus level, the median pairwise Spearman correlation between within-cohort paired CLR effect vectors was ",
    sprintf("%.2f",g$median_spearman_rho[[1]]),
    ", with a median directional agreement of ",
    sprintf("%.1f",100*g$median_sign_agreement[[1]]),
    "% and a median overlap among the strongest taxonomic effects of ",
    sprintf("%.1f",100*g$median_top_jaccard[[1]]),
    "%."
  ),
  "",
  paste0(
    "At the family level, the corresponding median effect-vector correlation was ",
    sprintf("%.2f",f$median_spearman_rho[[1]]),
    ", with ",
    sprintf("%.1f",100*f$median_sign_agreement[[1]]),
    "% directional agreement."
  ),
  "",
  paste0(
    "Among taxa shared across all three cohorts, only ",
    sprintf("%.1f",100*three_primary_g$proportion_same_direction_all3[[1]]),
    "% of genera and ",
    sprintf("%.1f",100*three_primary_f$proportion_same_direction_all3[[1]]),
    "% of families changed in the same direction in all three cohorts."
  ),
  "",
  "These findings support a model in which a more reproducible ecosystem-level displacement can arise through heterogeneous taxonomic routes, rather than through a universal longitudinal taxonomic signature.",
  "",
  "GUARDRAILS:",
  "Do not state that ecological and taxonomic reproducibility are directly numerically comparable.",
  "Do not claim universal ecological trajectories at the individual-patient level.",
  "Do not claim a universal sepsis-specific taxonomic signature.",
  "Do not interpret modest concordance between PRJNA851469 and PRJNA516701 as evidence against taxonomic heterogeneity; the overall three-cohort pattern remains low/modest and heterogeneous."
)

writeLines(
  wording,
  file.path(
    OUT,
    "03_FROZEN_MAIN_TEXT_WORDING.txt"
  )
)

discussion <- c(
  "DISCUSSION INTEGRATION POINTS",
  "",
  "1. The new analysis strengthens the central distinction between ecological-level reproducibility and taxonomic-level heterogeneity.",
  "2. It should be discussed as a difference in biological resolution, not as a direct statistical contest between Bray displacement and CLR-effect correlations.",
  "3. The strongest pairwise taxonomic concordance occurred between PRJNA851469 and PRJNA516701, indicating that taxonomic overlap is not absent; rather, it is incomplete and cohort-dependent.",
  "4. The non-sepsis surgical control showed weak or inverse genus-level concordance with the progressive ICU cohorts, which is directionally compatible with its recovery/re-convergence ecological behavior.",
  "5. The PRJNA851469 Day-3 clinical outcome analysis was null and should not be used to claim prognostic value of ecological displacement.",
  "6. The null clinical analysis can be mentioned in limitations or supplementary results if useful, but should not dilute the main ecological reproducibility story."
)

writeLines(
  discussion,
  file.path(
    OUT,
    "04_DISCUSSION_INTEGRATION_POINTS.txt"
  )
)

figure_plan <- c(
  "FIGURE ARCHITECTURE RECOMMENDATION",
  "",
  "Preferred option:",
  "Integrate cross-cohort taxonomic reproducibility as a new panel in Figure 3 rather than creating an unrelated standalone Figure 6.",
  "",
  "Suggested Figure 3 structure:",
  "A. Core-sepsis ecological displacement trajectory.",
  "B. Individual taxonomic-route heterogeneity.",
  "C. Relationship between trajectory heterogeneity and ecological consistency.",
  "D. NEW: cross-cohort GENUS/FAMILY effect-vector concordance across PRJNA691455, PRJNA851469, and PRJNA516701.",
  "",
  "If journal layout becomes crowded, move the detailed FAMILY/GENUS pairwise concordance plot to Supplement and retain a concise synthesis panel in the main figure.",
  "",
  "Do not promote the null PRJNA851469 clinical-outcome plot to a main figure."
)

writeLines(
  figure_plan,
  file.path(
    OUT,
    "05_FIGURE_ARCHITECTURE_RECOMMENDATION.txt"
  )
)

# Copy key provenance source files.
file.copy(
  file.path(TAX_DIR,"03_PRIMARY3_PAIRWISE_CONCORDANCE.csv"),
  file.path(OUT,"SOURCE_primary3_pairwise_concordance.csv"),
  overwrite=TRUE
)

file.copy(
  file.path(TAX_DIR,"06_THREEWAY_DIRECTION_SUMMARY.csv"),
  file.path(OUT,"SOURCE_threeway_direction_summary.csv"),
  overwrite=TRUE
)

file.copy(
  file.path(TAX_DIR,"11_Figure_candidate_taxonomic_concordance.pdf"),
  file.path(OUT,"SOURCE_taxonomic_concordance_candidate.pdf"),
  overwrite=TRUE
)

file.copy(
  file.path(TAX_DIR,"11_Figure_candidate_taxonomic_concordance.png"),
  file.path(OUT,"SOURCE_taxonomic_concordance_candidate.png"),
  overwrite=TRUE
)

writeLines(
  c(
    paste0("Completed: ",Sys.time()),
    "STEP36J COMPLETE",
    "Cross-cohort taxonomic reproducibility frozen for manuscript integration.",
    paste0("Clinical outcome branch status: ",clinical_status),
    "No new statistical inference performed."
  ),
  file.path(
    OUT,
    "_STEP36J_COMPLETE.txt"
  )
)

cat("STEP36J COMPLETE\n")
