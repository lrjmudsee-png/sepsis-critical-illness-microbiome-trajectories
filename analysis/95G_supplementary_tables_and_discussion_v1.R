# ============================================================
# Sepsis V2 - Step95G
# SUPPLEMENTARY TABLES S1-S8 + DISCUSSION DRAFT V1
#
# Inputs:
# - Step95F3 final Table 1 + Results freeze
# - Step95B frozen manuscript source pack
# - Step95D2 frozen constants/provenance
# - Step95E2 corrected longitudinal analysis table
#
# Purpose:
# 1) assemble manuscript-facing Supplementary Tables S1-S8;
# 2) retain provenance for every supplementary table;
# 3) draft an evidence-grounded Discussion v1;
# 4) use [REF] placeholders wherever prior literature must later be cited.
#
# NO new inferential analysis.
# ============================================================

options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(stringr)
})

ROOT <- "E:/sepsis_project"
RESULTS <- file.path(ROOT, "results")

F3 <- file.path(
  RESULTS,
  "V2_35F3_RESULTS_AND_TABLE1_FINAL_FREEZE_SENTENCE_AUDIT"
)

SRC <- file.path(
  RESULTS,
  "V2_35B_MANUSCRIPT_SOURCE_PACK_ASSEMBLY"
)

D2 <- file.path(
  RESULTS,
  "V2_35D2_PUBLICATION_FIGURES_SOURCEPACK_ONLY"
)

E2 <- file.path(
  RESULTS,
  "V2_35E2_TABLE1_ANALYSIS_POPULATION_FIX_AND_RESULTS_V2"
)

OUT <- file.path(
  RESULTS,
  "V2_35G_SUPPLEMENTARY_TABLES_AND_DISCUSSION_V1"
)

SUPPDIR <- file.path(OUT, "01_SUPPLEMENTARY_TABLES")
DISCDIR <- file.path(OUT, "02_DISCUSSION_DRAFT")
AUDDIR <- file.path(OUT, "03_PROVENANCE_AND_QC")

for (d in c(OUT,SUPPDIR,DISCDIR,AUDDIR)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

if (!file.exists(file.path(F3, "_STEP95F3_COMPLETE.ok"))) {
  stop("Step95F3 final Table 1 + Results freeze is incomplete.")
}

# ------------------------------------------------------------
# Helpers
# ------------------------------------------------------------

fmt <- function(x, digits=3) {
  ifelse(
    is.na(x),
    "NA",
    formatC(x, format="f", digits=digits)
  )
}

fmtp <- function(x) {
  ifelse(
    is.na(x),
    "NA",
    formatC(x, format="g", digits=3)
  )
}

first_present <- function(nms, choices) {
  x <- choices[choices %in% nms]
  if (length(x) == 0) return(NA_character_)
  x[1]
}

# ------------------------------------------------------------
# Frozen sources
# ------------------------------------------------------------

W <- read_csv(
  file.path(
    E2,
    "01_TABLES",
    "Main_or_Supp_Table_Cross_Cohort_Longitudinal_Results.csv"
  ),
  show_col_types=FALSE
)

CORE_RAW <- read_csv(
  file.path(
    SRC,
    "01_MAIN_TABLE_SOURCES",
    "MAIN_T3_CORE_TRAJECTORY_ENDPOINTS.csv"
  ),
  show_col_types=FALSE
)

EXT_RAW <- read_csv(
  file.path(
    SRC,
    "01_MAIN_TABLE_SOURCES",
    "MAIN_T4_EXTERNAL_VALIDATION.csv"
  ),
  show_col_types=FALSE
)

ECO <- read_csv(
  file.path(
    SRC,
    "02_SUPPLEMENTARY_TABLE_SOURCES",
    "SUPP_S4_SOURCE_ECOLOGY_REGISTRY.csv"
  ),
  show_col_types=FALSE
)

SENS <- read_csv(
  file.path(
    SRC,
    "02_SUPPLEMENTARY_TABLE_SOURCES",
    "SUPP_S5_SOURCE_BRAY_SENSITIVITY.csv"
  ),
  show_col_types=FALSE
)

SLOPES <- read_csv(
  file.path(
    SRC,
    "02_SUPPLEMENTARY_TABLE_SOURCES",
    "SUPP_S6_SOURCE_SLOPES.csv"
  ),
  show_col_types=FALSE
)

BASEB <- read_csv(
  file.path(
    SRC,
    "02_SUPPLEMENTARY_TABLE_SOURCES",
    "SUPP_S7_SOURCE_BASELINE_BETA.csv"
  ),
  show_col_types=FALSE
)

GENUS_BASE <- read_csv(
  file.path(
    SRC,
    "02_SUPPLEMENTARY_TABLE_SOURCES",
    "SUPP_S8_GENUS_BASELINE_HITS.csv"
  ),
  show_col_types=FALSE
)

GENUS_LONG <- read_csv(
  file.path(
    SRC,
    "02_SUPPLEMENTARY_TABLE_SOURCES",
    "SUPP_S9_GENUS_LONGITUDINAL_SENSITIVITY_ONLY.csv"
  ),
  show_col_types=FALSE
)

GENUS_FREEZE <- read_csv(
  file.path(
    SRC,
    "02_SUPPLEMENTARY_TABLE_SOURCES",
    "SUPP_S10_GENUS_BRANCH_FREEZE.csv"
  ),
  show_col_types=FALSE
)

GENUS_REC <- read_csv(
  file.path(
    SRC,
    "02_SUPPLEMENTARY_TABLE_SOURCES",
    "SUPP_S11_GENUS_REPORTING_RECOMMENDATION.csv"
  ),
  show_col_types=FALSE
)

CONST <- read_csv(
  file.path(
    D2,
    "03_PROVENANCE",
    "01_FROZEN_CONSTANTS_USED.csv"
  ),
  show_col_types=FALSE
)

cval <- function(metric) {
  CONST$value[CONST$metric == metric][1]
}

# ------------------------------------------------------------
# S1 - Full cross-cohort longitudinal evidence
# ------------------------------------------------------------

S1 <- W %>%
  transmute(
    Cohort = Cohort,
    `Cohort role` = Cohort_role,
    `Primary paired participants` = Primary_pairs,
    `Primary dz` = Primary_effect_dz,
    `Primary p` = Primary_p,
    `Primary BH-FDR` = Primary_FDR,
    `Primary direction` = Primary_direction,
    `Common-anchor paired participants` = Common_anchor_pairs,
    `Common-anchor dz` = Common_anchor_effect_dz,
    `Common-anchor p` = Common_anchor_p,
    `Common-anchor BH-FDR` = Common_anchor_FDR,
    `Direction robust` = Direction_robust,
    `Frozen evidence class` = Evidence_class
  )

write_csv(
  S1,
  file.path(
    SUPPDIR,
    "Supplementary_Table_S1_Cross_Cohort_Longitudinal_Evidence.csv"
  )
)

# ------------------------------------------------------------
# S2 - Core-sepsis trajectory ecology endpoints
# ------------------------------------------------------------

S2_summary <- tribble(
  ~Metric, ~Value, ~Interpretation,
  "Median Bray displacement at Day 3",
    cval("Day3_Bray_median"),
    "Ecological displacement from patient-specific baseline",
  "Median Bray displacement at Day 7",
    cval("Day7_Bray_median"),
    "Ecological displacement from patient-specific baseline",
  "Paired Day 3-Day 7 participants",
    cval("Day3_to_Day7_paired_n"),
    "Participants contributing to paired Day 3-Day 7 change",
  "Mean paired Day 3-Day 7 Bray change",
    cval("Day3_to_Day7_mean_change"),
    "Positive values indicate progressive displacement",
  "Mean pairwise signed-trajectory Euclidean distance",
    cval("pairwise_signed_Euclidean_mean"),
    "Between-patient taxonomic-route heterogeneity",
  "Median pairwise signed-trajectory Euclidean distance",
    cval("pairwise_signed_Euclidean_median"),
    "Between-patient taxonomic-route heterogeneity",
  "Mean pairwise cosine similarity",
    cval("mean_cosine_similarity"),
    "Directional similarity of signed taxonomic trajectories",
  "Fraction of patient pairs with cosine <= 0",
    cval("fraction_cosine_le_0"),
    "Opposed or non-aligned taxonomic trajectory directions",
  "Spearman rho: trajectory distance vs latest Bray",
    cval("rho_latest_Bray"),
    "No clear association with latest Bray displacement",
  "p: trajectory distance vs latest Bray",
    cval("p_latest_Bray"),
    "Nominal correlation p value",
  "Spearman rho: trajectory distance vs latest EII",
    cval("rho_latest_EII"),
    "Exploratory ecological-consistency association",
  "Permutation p: trajectory distance vs latest EII",
    cval("perm_p_latest_EII"),
    "10,000-permutation support",
  "Spearman rho: trajectory distance vs Simpson instability",
    cval("rho_Simpson_instability"),
    "Ecological-consistency association",
  "Permutation p: trajectory distance vs Simpson instability",
    cval("perm_p_Simpson_instability"),
    "10,000-permutation support"
)

write_csv(
  S2_summary,
  file.path(
    SUPPDIR,
    "Supplementary_Table_S2_Core_Sepsis_Trajectory_Ecology.csv"
  )
)

# Preserve source endpoint table as an audit companion.
write_csv(
  CORE_RAW,
  file.path(
    AUDDIR,
    "S2_SOURCE_Core_Trajectory_Endpoint_Table.csv"
  )
)

# ------------------------------------------------------------
# S3 - External Control-Trauma-Sepsis validation
# ------------------------------------------------------------

S3 <- tribble(
  ~Evidence_family, ~Metric, ~Value_or_range, ~Reporting_role,
  "Analysis population",
    "Control / Trauma / Sepsis",
    "17 / 18 / 18",
    "External supportive",
  "Healthy-centroid displacement",
    "Trauma median distance",
    fmt(cval("Trauma_healthy_centroid_distance"),3),
    "Primary ecological state",
  "Healthy-centroid displacement",
    "Sepsis median distance",
    fmt(cval("Sepsis_healthy_centroid_distance"),3),
    "Primary ecological state",
  "Healthy-centroid displacement",
    "Sepsis minus Trauma difference",
    fmt(cval("Sepsis_minus_Trauma_centroid_difference"),3),
    "Primary ecological state",
  "Healthy-centroid displacement",
    "BH-FDR",
    fmtp(cval("centroid_difference_FDR")),
    "Primary ecological state",
  "CLR/Aitchison robustness",
    "PERMANOVA R2 range",
    paste0(
      fmt(cval("CLR_Aitchison_R2_min"),3),
      "-",
      fmt(cval("CLR_Aitchison_R2_max"),3)
    ),
    "Robust location evidence",
  "CLR/Aitchison robustness",
    "PERMANOVA FDR range",
    paste0(
      fmtp(cval("CLR_Aitchison_PERMANOVA_FDR_min")),
      "-",
      fmtp(cval("CLR_Aitchison_PERMANOVA_FDR_max"))
    ),
    "Robust location evidence",
  "CLR/Aitchison robustness",
    "Dispersion FDR range",
    paste0(
      fmt(cval("CLR_Aitchison_dispersion_FDR_min"),3),
      "-",
      fmt(cval("CLR_Aitchison_dispersion_FDR_max"),3)
    ),
    "No significant dispersion heterogeneity",
  "Taxonomic alignment",
    "External genus-level result",
    "Exploratory only",
    "Does not establish universal sepsis-specific taxonomic signature"
)

write_csv(
  S3,
  file.path(
    SUPPDIR,
    "Supplementary_Table_S3_External_Control_Trauma_Sepsis_Validation.csv"
  )
)

write_csv(
  EXT_RAW,
  file.path(
    AUDDIR,
    "S3_SOURCE_External_Taxonomic_Summary.csv"
  )
)

# ------------------------------------------------------------
# S4 - CRA002354 primary ecology + source-specific slopes
# ------------------------------------------------------------

source_col <- first_present(
  names(SLOPES),
  c("pulmonary_binary","source_group","group","source")
)

if (is.na(source_col)) {
  stop("Could not identify source-group column in source-specific slope table.")
}

S4_slopes <- SLOPES %>%
  mutate(
    source_group = as.character(.data[[source_col]])
  ) %>%
  select(
    source_group,
    everything(),
    -all_of(source_col)
  )

write_csv(
  ECO,
  file.path(
    SUPPDIR,
    "Supplementary_Table_S4A_CRA002354_Primary_Source_Ecology_Registry.csv"
  )
)

write_csv(
  S4_slopes,
  file.path(
    SUPPDIR,
    "Supplementary_Table_S4B_CRA002354_Source_Specific_Slopes.csv"
  )
)

# ------------------------------------------------------------
# S5 - CRA002354 sensitivities + baseline beta robustness
# ------------------------------------------------------------

write_csv(
  SENS,
  file.path(
    SUPPDIR,
    "Supplementary_Table_S5A_CRA002354_Bray_Sensitivity_Analyses.csv"
  )
)

write_csv(
  BASEB,
  file.path(
    SUPPDIR,
    "Supplementary_Table_S5B_CRA002354_Baseline_Beta_Robustness.csv"
  )
)

# ------------------------------------------------------------
# S6 - Baseline genus hits + coverage robustness
# ------------------------------------------------------------

write_csv(
  GENUS_BASE,
  file.path(
    SUPPDIR,
    "Supplementary_Table_S6_CRA002354_Baseline_Genus_Hits_and_Coverage_Robustness.csv"
  )
)

# ------------------------------------------------------------
# S7 - Longitudinal genus sensitivity-only signals + reporting
# ------------------------------------------------------------

write_csv(
  GENUS_LONG,
  file.path(
    SUPPDIR,
    "Supplementary_Table_S7A_CRA002354_Longitudinal_Genus_Sensitivity_Only.csv"
  )
)

write_csv(
  GENUS_REC,
  file.path(
    SUPPDIR,
    "Supplementary_Table_S7B_CRA002354_Genus_Reporting_Recommendations.csv"
  )
)

write_csv(
  GENUS_FREEZE,
  file.path(
    AUDDIR,
    "S7_SOURCE_Final_Genus_Branch_Freeze.csv"
  )
)

# ------------------------------------------------------------
# S8 - Exploratory EII construction / ecological consistency
# ------------------------------------------------------------

S8 <- tribble(
  ~Component, ~Definition_or_result, ~Role,
  "EII construction",
    "Row-mean of standardized Bray-Curtis displacement and absolute diversity-change components, rescaled to 0-100",
    "Exploratory internal ecological index",
  "EII clinical status",
    "Not a clinical prognostic score and not independently validated",
    "Guardrail",
  "Trajectory distance vs latest EII",
    paste0(
      "Spearman rho=",
      fmt(cval("rho_latest_EII"),3),
      "; permutation p=",
      fmtp(cval("perm_p_latest_EII"))
    ),
    "Exploratory ecological-consistency association",
  "Trajectory distance vs Simpson instability",
    paste0(
      "Spearman rho=",
      fmt(cval("rho_Simpson_instability"),3),
      "; permutation p=",
      fmtp(cval("perm_p_Simpson_instability"))
    ),
    "Ecological-consistency association",
  "Trajectory distance vs latest Bray",
    paste0(
      "Spearman rho=",
      fmt(cval("rho_latest_Bray"),3),
      "; p=",
      fmtp(cval("p_latest_Bray"))
    ),
    "No clear association"
)

write_csv(
  S8,
  file.path(
    SUPPDIR,
    "Supplementary_Table_S8_Exploratory_EII_and_Ecological_Consistency.csv"
  )
)

# ------------------------------------------------------------
# Supplementary-table provenance registry
# ------------------------------------------------------------

supp_registry <- tribble(
  ~Supplementary_item, ~Authoritative_source, ~Inference_status,
  "Table S1",
    "Step93W frozen longitudinal evidence via Step95E2",
    "No new inference",
  "Table S2",
    "Step93U frozen trajectory/ecology evidence via Step95D2 constants",
    "No new inference",
  "Table S3",
    "Step93X2/X3/X4 frozen external-validation evidence",
    "No new inference",
  "Table S4A-S4B",
    "Step94B2D frozen infection-source ecology + source-specific slopes",
    "No new inference",
  "Table S5A-S5B",
    "Step94B2D frozen Bray sensitivities + baseline beta robustness",
    "No new inference",
  "Table S6",
    "Step94B3C baseline genus registry after Step94B3B coverage sensitivity",
    "Exploratory only; no new inference",
  "Table S7A-S7B",
    "Step94B3C longitudinal sensitivity-only genus registry + reporting recommendations",
    "Exploratory only; no new inference",
  "Table S8",
    "Step92B2 / Step93U frozen exploratory EII evidence",
    "Exploratory only; no new inference"
)

write_csv(
  supp_registry,
  file.path(
    AUDDIR,
    "SUPPLEMENTARY_TABLE_PROVENANCE_REGISTRY.csv"
  )
)

# ------------------------------------------------------------
# Discussion Draft v1
# Evidence-grounded only. Literature comparison placeholders are explicit.
# ------------------------------------------------------------

discussion <- paste0(
"DISCUSSION DRAFT V1 — EVIDENCE-GROUNDED WITH LITERATURE PLACEHOLDERS\n\n",

"Principal findings\n\n",
"Across multiple longitudinal critical-illness cohorts, we observed a recurring pattern of increasing within-patient gut-microbiome ecological displacement over time. This signal was evident in the core repeated-sepsis cohort and in two ICU-background cohorts, with an additional large positive contrast in an external ICU-infection cohort. Importantly, the non-sepsis surgical longitudinal control showed a different pattern, with recovery-direction change in the primary comparison and significant re-convergence in the common-anchor sensitivity analysis. These observations support progressive ecological displacement as a reproducible feature of several critical-illness settings, but they also argue against interpreting this pattern as uniquely sepsis-specific. [REF: longitudinal microbiome disruption during critical illness and sepsis]\n\n",

"Ecological displacement may be more reproducible than individual taxonomic trajectories\n\n",
"In the core repeated-sepsis cohort, increasing Bray-Curtis displacement was accompanied by marked heterogeneity in the signed genus-level trajectories followed by individual patients. Pairwise trajectory distances were substantial, mean cosine similarity was low, and more than one-third of patient-pair comparisons showed non-aligned or opposed directions. At the same time, greater trajectory heterogeneity was associated with the exploratory ecological injury index and with Simpson-diversity instability. Together, these findings suggest that patients can move toward a comparably disrupted ecological state through different taxonomic routes. This distinction may help explain why taxon-specific signatures reported across sepsis microbiome studies have often been difficult to reproduce across cohorts, sequencing strategies, and clinical settings. [REF: cross-cohort inconsistency of sepsis-associated taxa; ecological-state approaches in microbiome research]\n\n",

"External validation supports ecological-state differentiation rather than a universal taxonomic signature\n\n",
"The external Control-Trauma-Sepsis cohort provided an important test of whether the observed ecological changes simply reflected severe physiological stress. Sepsis samples were farther from the healthy-control centroid than severe-trauma samples, and the Sepsis-versus-Trauma separation remained robust in CLR/Aitchison analyses across pseudocount choices without significant dispersion heterogeneity. This supports a sepsis-associated ecological-state difference beyond severe trauma in that cohort. However, the external analysis did not establish a reproducible universal genus-level signature. The combined evidence therefore favors an ecological-state interpretation over a fixed biomarker panel: broad community displacement may be more transportable across cohorts than the identities or directions of individual taxa. [REF: trauma microbiome studies; compositional microbiome methods; sepsis biomarker studies]\n\n",

"Infection source showed limited modification of the longitudinal ecological trajectory\n\n",
"In CRA002354, pulmonary and recorded non-pulmonary sepsis did not show a clear difference in the rate of longitudinal Bray-Curtis displacement. Estimated slopes were positive in both source groups, and the direction of the source-by-time interaction remained positive across prespecified sensitivity analyses, but neither the primary interaction nor either source-specific slope was statistically significant. A nominal baseline source association was accompanied by dispersion heterogeneity and disappeared when baseline sampling was restricted to ICU Day 3 or earlier. These results do not support a strong infection-source modification of the longitudinal ecological-displacement process. They are instead compatible with the possibility that critical illness exerts shared ecological pressures that partly transcend the anatomical origin of infection, although the current analysis cannot establish mechanistic equivalence between source groups. [REF: gut microbiome across sepsis infection sources; gut-lung axis / systemic critical-illness ecology]\n\n",

"Genus-level source associations were exploratory and sensitive to taxonomic coverage\n\n",
"Several baseline genera differed between pulmonary and recorded non-pulmonary sepsis in the primary exploratory analysis, but robustness decreased after restricting to samples with higher genus-level taxonomic assignment coverage. None of the primary baseline hits remained FDR-significant at both the >=70% and >=80% coverage thresholds, and no genus-level source-by-time interaction survived FDR correction in the primary longitudinal analysis. A Bacteroides longitudinal signal appeared only in the >=70% coverage sensitivity analysis and did not remain significant at >=80% coverage. These results reinforce the need to distinguish exploratory taxonomic signals from reproducible ecological effects, particularly when public 16S datasets differ in sequence quality, taxonomic resolution, and preprocessing constraints. [REF: robustness and reproducibility of 16S genus-level differential abundance]\n\n",

"Implications\n\n",
"Taken together, the V2 results support a framework in which longitudinal ecological displacement is the more consistent cross-cohort signal, whereas the taxonomic routes producing that displacement are heterogeneous. This framework has two practical implications. First, future microbiome studies in sepsis and critical illness may benefit from prioritizing within-patient ecological change, trajectory geometry, and compositional state measures rather than relying exclusively on single-time-point taxonomic biomarkers. Second, external validation should test whether ecological patterns remain distinguishable from relevant critically ill comparators, not only from healthy controls. These ideas require prospective testing before they can inform clinical prediction or intervention. [REF: longitudinal precision microbiome / ecological trajectory methods]\n\n",

"Strengths and limitations\n\n",
"This study has several strengths. It used prespecified cohort roles, separated core, supportive, control, external-validation, and infection-source analyses, and preserved within-cohort taxonomic features rather than directly merging ASVs or OTUs generated from heterogeneous 16S regions. The principal ecological conclusions were evaluated across multiple cohorts and supported by sensitivity analyses, while exploratory taxonomic findings were explicitly downgraded when they lacked robustness.\n\n",
"Several limitations should also be considered. First, the study relied on public datasets with heterogeneous patient populations, sampling schedules, sequencing regions, laboratory protocols, and available clinical metadata. These differences limit direct cross-cohort taxonomic harmonization and preclude causal interpretation. Second, the core repeated-sepsis trajectory analysis was small, including nine patients in the paired Day 3-to-Day 7 comparison, and therefore the continuous trajectory correlations should be regarded as exploratory. Third, the ecological injury index was constructed as an internal exploratory summary of ecological displacement and diversity instability; it was not linked to sufficiently complete clinical outcomes and should not be interpreted as a validated prognostic score. Fourth, PRJNA1010969 provided cross-sectional external ecological-state validation rather than independent longitudinal replication. Fifth, CRA002354 was distributed as FASTA sequences without per-base quality scores. We therefore used recoverable sequence filters followed by 97% OTU reconstruction and SILVA-based taxonomy rather than claiming an ASV-level reanalysis equivalent to the original raw-read pipeline. Finally, infection-source analyses used relatively broad source categories, and smaller non-pulmonary subgroups limited more granular source-specific inference.\n\n",

"Conclusion\n\n",
"Across heterogeneous public cohorts, progressive gut-microbiome ecological displacement emerged as a recurring longitudinal feature of critical illness, while the taxonomic routes underlying that displacement varied substantially between patients. External comparison supported sepsis-associated ecological-state differentiation beyond severe trauma without identifying a universal sepsis-specific taxonomic signature, and infection source did not clearly modify the longitudinal ecological-displacement trajectory. These findings support a shift from fixed taxonomic signatures toward longitudinal ecological-state and trajectory-based approaches for studying microbiome disruption in sepsis and critical illness. [REF: concluding contextual citations]\n"
)

writeLines(
  discussion,
  file.path(
    DISCDIR,
    "Discussion_Draft_v1_Evidence_Grounded.txt"
  )
)

# ------------------------------------------------------------
# Discussion claim guardrails
# ------------------------------------------------------------

discussion_guardrails <- c(
  "Do not describe ecological displacement as sepsis-specific because ICU-background cohorts also progressed.",
  "Do not claim stable taxonomic trajectory subtypes from the 9-vs-1 exploratory clustering.",
  "Do not present EII as a clinical prognostic or validated score.",
  "Do not claim the external cohort validates a universal sepsis-specific genus signature.",
  "Do not claim pulmonary and non-pulmonary trajectories are identical or equivalent.",
  "Do not state that infection source has no effect; use 'no clear modification detected'.",
  "Do not present Enterococcus as a robust infection-source biomarker.",
  "Do not present Bacteroides as a validated longitudinal source signal.",
  "Do not imply that CRA002354 was reprocessed with DADA2/ASV quality filtering; FASTA lacked per-base qualities.",
  "Do not merge ASV/OTU identities across heterogeneous cohorts/16S regions."
)

writeLines(
  c(
    "DISCUSSION GUARDRAILS",
    "",
    paste0("- ",discussion_guardrails)
  ),
  file.path(
    DISCDIR,
    "Discussion_Guardrails.txt"
  )
)

# ------------------------------------------------------------
# Literature-integration checklist
# ------------------------------------------------------------

lit_plan <- tribble(
  ~Discussion_section, ~Literature_needed,
  "Principal findings",
    "Longitudinal gut microbiome disruption in sepsis/ICU critical illness; comparator studies with non-sepsis ICU patients",
  "Ecological displacement vs taxonomic routes",
    "Cross-cohort inconsistency of sepsis taxa; ecological-state / trajectory-based microbiome methods",
  "External trauma comparison",
    "Trauma-associated microbiome disruption; compositional CLR/Aitchison methodology; severe-illness comparator studies",
  "Infection source",
    "Microbiome differences by infection source; gut-lung axis; systemic effects of critical illness",
  "Taxonomic robustness",
    "16S taxonomic-resolution limitations; reproducibility of differential-abundance biomarkers",
  "Implications",
    "Longitudinal personalized microbiome trajectories; ecological-state biomarkers",
  "Limitations",
    "Public-dataset heterogeneity, batch effects, 16S region heterogeneity, FASTA-only reanalysis constraints"
)

write_csv(
  lit_plan,
  file.path(
    DISCDIR,
    "Discussion_Literature_Integration_Plan.csv"
  )
)

# ------------------------------------------------------------
# Final QC
# ------------------------------------------------------------

supp_expected <- c(
  "Supplementary_Table_S1_Cross_Cohort_Longitudinal_Evidence.csv",
  "Supplementary_Table_S2_Core_Sepsis_Trajectory_Ecology.csv",
  "Supplementary_Table_S3_External_Control_Trauma_Sepsis_Validation.csv",
  "Supplementary_Table_S4A_CRA002354_Primary_Source_Ecology_Registry.csv",
  "Supplementary_Table_S4B_CRA002354_Source_Specific_Slopes.csv",
  "Supplementary_Table_S5A_CRA002354_Bray_Sensitivity_Analyses.csv",
  "Supplementary_Table_S5B_CRA002354_Baseline_Beta_Robustness.csv",
  "Supplementary_Table_S6_CRA002354_Baseline_Genus_Hits_and_Coverage_Robustness.csv",
  "Supplementary_Table_S7A_CRA002354_Longitudinal_Genus_Sensitivity_Only.csv",
  "Supplementary_Table_S7B_CRA002354_Genus_Reporting_Recommendations.csv",
  "Supplementary_Table_S8_Exploratory_EII_and_Ecological_Consistency.csv"
)

supp_audit <- tibble(
  file=supp_expected,
  exists=file.exists(file.path(SUPPDIR,supp_expected))
)

write_csv(
  supp_audit,
  file.path(
    AUDDIR,
    "SUPPLEMENTARY_TABLE_FILE_AUDIT.csv"
  )
)

discussion_file <- file.path(
  DISCDIR,
  "Discussion_Draft_v1_Evidence_Grounded.txt"
)

ready <- all(supp_audit$exists) &&
  file.exists(discussion_file)

qc <- tibble(
  step95F3_complete=TRUE,
  supplementary_items_complete=all(supp_audit$exists),
  discussion_draft_v1_created=file.exists(discussion_file),
  literature_claims_left_as_REF_placeholders=TRUE,
  no_new_inferential_analysis=TRUE,
  ready_for_supplementary_review_and_literature_integration=ready
)

write_csv(
  qc,
  file.path(
    OUT,
    "STEP95G_READINESS.csv"
  )
)

if (!ready) {
  stop("Step95G readiness failed.")
}

writeLines(
  c(
    paste0("Completed: ",Sys.time()),
    "Supplementary Tables S1-S8 assembled: TRUE",
    "Discussion Draft v1 created: TRUE",
    "Literature comparisons use explicit [REF] placeholders: TRUE",
    "No new inferential analysis: TRUE",
    "Ready for Supplementary review + literature integration: TRUE",
    "STEP95G COMPLETE"
  ),
  file.path(
    OUT,
    "_STEP95G_COMPLETE.ok"
  )
)

cat("STEP95G COMPLETE\n")
