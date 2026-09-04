# ============================================================
# Sepsis V2 - Step95F
# RESULTS POLISH + MANUSCRIPT TABLE 1 FREEZE
#
# Inputs:
# - Step95E2 corrected analysis-population Table 1
# - Step95B frozen result registries
# - Step95D2 frozen constants
#
# Purpose:
# 1) convert internal Table 1 into a journal-facing cohort table;
# 2) polish Results Draft v2 into manuscript-facing Results Draft v3;
# 3) add exact baseline-beta robustness numbers;
# 4) explicitly label EII as exploratory;
# 5) freeze table placement and supplementary-table plan.
#
# NO new inferential analysis.
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

E2 <- file.path(
  RESULTS,
  "V2_35E2_TABLE1_ANALYSIS_POPULATION_FIX_AND_RESULTS_V2"
)

SRC <- file.path(
  RESULTS,
  "V2_35B_MANUSCRIPT_SOURCE_PACK_ASSEMBLY"
)

D2 <- file.path(
  RESULTS,
  "V2_35D2_PUBLICATION_FIGURES_SOURCEPACK_ONLY"
)

OUT <- file.path(
  RESULTS,
  "V2_35F_RESULTS_POLISH_AND_TABLE1_MANUSCRIPT_FREEZE"
)

TABDIR <- file.path(OUT, "01_FINAL_TABLES")
DRAFTDIR <- file.path(OUT, "02_RESULTS")
AUDDIR <- file.path(OUT, "03_AUDIT_AND_PLACEMENT")

for (d in c(OUT,TABDIR,DRAFTDIR,AUDDIR)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

if (!file.exists(file.path(E2, "_STEP95E2_COMPLETE.ok"))) {
  stop("Step95E2 is incomplete.")
}

# ------------------------------------------------------------
# 1. Read corrected Table 1 and frozen sources
# ------------------------------------------------------------

T1 <- read_csv(
  file.path(
    E2,
    "01_TABLES",
    "Main_Table_1_Cohort_and_Study_Characteristics_ANALYSIS_POPULATION.csv"
  ),
  show_col_types = FALSE
)

W <- read_csv(
  file.path(
    E2,
    "01_TABLES",
    "Main_or_Supp_Table_Cross_Cohort_Longitudinal_Results.csv"
  ),
  show_col_types = FALSE
)

ECO <- read_csv(
  file.path(
    SRC,
    "02_SUPPLEMENTARY_TABLE_SOURCES",
    "SUPP_S4_SOURCE_ECOLOGY_REGISTRY.csv"
  ),
  show_col_types = FALSE
)

SLOPES <- read_csv(
  file.path(
    SRC,
    "02_SUPPLEMENTARY_TABLE_SOURCES",
    "SUPP_S6_SOURCE_SLOPES.csv"
  ),
  show_col_types = FALSE
)

BASEB <- read_csv(
  file.path(
    SRC,
    "02_SUPPLEMENTARY_TABLE_SOURCES",
    "SUPP_S7_SOURCE_BASELINE_BETA.csv"
  ),
  show_col_types = FALSE
)

GENUS <- read_csv(
  file.path(
    SRC,
    "02_SUPPLEMENTARY_TABLE_SOURCES",
    "SUPP_S10_GENUS_BRANCH_FREEZE.csv"
  ),
  show_col_types = FALSE
)

CONST <- read_csv(
  file.path(
    D2,
    "03_PROVENANCE",
    "01_FROZEN_CONSTANTS_USED.csv"
  ),
  show_col_types = FALSE
)

cval <- function(metric) {
  CONST$value[CONST$metric == metric][1]
}

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

# ------------------------------------------------------------
# 2. Journal-facing Table 1
# ------------------------------------------------------------

T1_clean <- T1 %>%
  transmute(
    Cohort = Cohort,
    `Clinical context` = Clinical_context,
    `Samples analyzed` = Samples_analyzed,
    `Participants analyzed` = Participants_analyzed,
    `Participants with repeated samples` =
      Participants_with_repeated_samples,
    `Observed sampling window` =
      ifelse(
        is.na(Observed_time_window) |
        Observed_time_window == "NR",
        "NR",
        Observed_time_window
      ),
    `Primary analytical use` = case_when(
      Cohort == "PRJNA691455" ~
        "Core longitudinal sepsis trajectory analysis",
      Cohort == "PRJEB82425" ~
        "External longitudinal ICU-infection support",
      Cohort %in% c("PRJNA516701","PRJNA851469") ~
        "ICU-background longitudinal context",
      Cohort == "PRJNA578267" ~
        "Non-sepsis longitudinal control",
      Cohort %in% c("PRJNA430161","PRJNA1166732") ~
        "Intervention-support longitudinal analysis",
      Cohort == "PRJNA978257" ~
        "Static supportive analysis",
      Cohort == "PRJNA1010969" ~
        "External Control-Trauma-Sepsis validation",
      Cohort == "CRA002354" ~
        "Infection-source longitudinal extension",
      TRUE ~ Analysis_role
    )
  )

write_csv(
  T1_clean,
  file.path(
    TABDIR,
    "Main_Table_1_Cohort_Characteristics_MANUSCRIPT.csv"
  )
)

writeLines(
  c(
    "Table 1. Cohort and analysis-population characteristics.",
    "",
    "Footnotes:",
    "Counts refer to the samples and participants entering the corresponding V2 analytical branch, not to all samples originally deposited in the public study.",
    "For the eight harmonized Step87B cohorts, counts were derived from the primary 785-observation analysis metadata after replicate resolution.",
    "PRJNA1010969 comprised 53 analyzed samples: 17 Control, 18 Trauma, and 18 Sepsis.",
    "CRA002354 comprised 130 primary-analysis samples from 64 participants after true-nonchimeric OTU reconstruction and the >=4000-read depth threshold.",
    "NR, a harmonized longitudinal day scale was not represented for that cohort in the final analysis metadata."
  ),
  file.path(
    TABDIR,
    "Main_Table_1_Cohort_Characteristics_MANUSCRIPT_Footnotes.txt"
  )
)

# ------------------------------------------------------------
# 3. Normalize slope registry
# ------------------------------------------------------------

source_col <- if ("pulmonary_binary" %in% names(SLOPES)) {
  "pulmonary_binary"
} else if ("source_group" %in% names(SLOPES)) {
  "source_group"
} else {
  stop("Could not identify source group in slope registry.")
}

S <- SLOPES %>%
  mutate(
    .source = as.character(.data[[source_col]])
  )

pulm <- S %>%
  filter(
    str_detect(toupper(.source), "PULMONARY") &
    !str_detect(toupper(.source), "NON")
  ) %>%
  slice(1)

nonp <- S %>%
  filter(
    str_detect(toupper(.source), "NON")
  ) %>%
  slice(1)

if (nrow(pulm)==0 || nrow(nonp)==0) {
  stop("Could not resolve pulmonary and non-pulmonary slope rows.")
}

primary_cra <- ECO %>%
  filter(
    domain == "PRIMARY_LONGITUDINAL_BRAY",
    estimand == "source_by_personal_time_interaction"
  ) %>%
  slice(1)

# ------------------------------------------------------------
# 4. Results Draft v3
# ------------------------------------------------------------

getW <- function(cohort) {
  W %>% filter(Cohort == cohort) %>% slice(1)
}

r_core <- getW("PRJNA691455")
r_851 <- getW("PRJNA851469")
r_516 <- getW("PRJNA516701")
r_ext <- getW("PRJEB82425")
r_ctrl <- getW("PRJNA578267")

sec1 <- paste0(
  "1. Cohort architecture and analytical framework\n\n",
  "The V2 analysis integrated ten public gut-microbiome cohorts with prespecified analytical roles (Table 1; Figure 1). ",
  "Eight cohorts were represented in the harmonized primary analysis object, including seven longitudinal cohorts and one static supportive cohort. ",
  "PRJNA691455 served as the core repeated-sepsis cohort for trajectory-level interpretation, PRJNA578267 as the non-sepsis longitudinal control, and PRJNA430161 and PRJNA1166732 as intervention-support cohorts. ",
  "PRJNA1010969 provided an external Control-Trauma-Sepsis ecological-state comparison using 53 analyzed samples (17 Control, 18 Trauma, and 18 Sepsis), whereas CRA002354 provided an independent infection-source extension using 130 primary-analysis samples from 64 participants. ",
  "Table 1 reports the populations entering each analytical branch rather than the total number of samples deposited in the original public studies."
)

sec2 <- paste0(
  "2. Progressive ecological displacement across longitudinal critical-illness cohorts\n\n",
  "Three cohorts showed FDR-supported progressive displacement in the primary early-to-late analysis: PRJNA691455 (n=",
  r_core$Primary_pairs,
  " paired participants, dz=",
  fmt(r_core$Primary_effect_dz,3),
  ", FDR=",
  fmtp(r_core$Primary_FDR),
  "), PRJNA851469 (n=",
  r_851$Primary_pairs,
  ", dz=",
  fmt(r_851$Primary_effect_dz,3),
  ", FDR=",
  fmtp(r_851$Primary_FDR),
  "), and PRJNA516701 (n=",
  r_516$Primary_pairs,
  ", dz=",
  fmt(r_516$Primary_effect_dz,3),
  ", FDR=",
  fmtp(r_516$Primary_FDR),
  "). ",
  "PRJEB82425 also showed a large positive primary contrast (n=",
  r_ext$Primary_pairs,
  ", dz=",
  fmt(r_ext$Primary_effect_dz,3),
  ", FDR=",
  fmtp(r_ext$Primary_FDR),
  "), but heterogeneous anchor structure supported its classification as external supportive evidence rather than core replication. ",
  "In contrast, the non-sepsis surgical control PRJNA578267 showed a negative recovery-direction primary contrast (n=",
  r_ctrl$Primary_pairs,
  ", dz=",
  fmt(r_ctrl$Primary_effect_dz,3),
  ", FDR=",
  fmtp(r_ctrl$Primary_FDR),
  ") and significant re-convergence in the common-anchor sensitivity analysis (n=",
  r_ctrl$Common_anchor_pairs,
  ", dz=",
  fmt(r_ctrl$Common_anchor_effect_dz,3),
  ", FDR=",
  fmtp(r_ctrl$Common_anchor_FDR),
  "). ",
  "Neither intervention-support cohort showed a clear paired longitudinal change after multiplicity correction. ",
  "Thus, progressive within-patient ecological displacement was observed across several critical-illness cohorts, whereas the non-sepsis longitudinal control showed a contrasting recovery/re-convergence pattern in sensitivity analysis (Figure 2). ",
  "Because progressive displacement was also present in ICU-background cohorts, these findings do not establish sepsis-specific longitudinal instability."
)

sec3 <- paste0(
  "3. Heterogeneous taxonomic routes accompany ecological displacement in core sepsis\n\n",
  "Within PRJNA691455, median Bray-Curtis displacement from the patient-specific baseline increased from ",
  fmt(cval("Day3_Bray_median"),3),
  " at Day 3 to ",
  fmt(cval("Day7_Bray_median"),3),
  " at Day 7. ",
  "Among ",
  as.integer(cval("Day3_to_Day7_paired_n")),
  " patients with paired Day 3 and Day 7 observations, the mean increase was ",
  fmt(cval("Day3_to_Day7_mean_change"),3),
  ". ",
  "Despite this shared ecological movement, patient-specific taxonomic trajectories were heterogeneous: the mean pairwise signed-trajectory Euclidean distance was ",
  fmt(cval("pairwise_signed_Euclidean_mean"),3),
  " (median ",
  fmt(cval("pairwise_signed_Euclidean_median"),3),
  "), the mean cosine similarity was ",
  fmt(cval("mean_cosine_similarity"),3),
  ", and ",
  fmt(100*cval("fraction_cosine_le_0"),1),
  "% of patient-pair comparisons had cosine similarity <=0. ",
  "Mean trajectory distance was positively associated with the exploratory ecological injury index (EII; Spearman rho=",
  fmt(cval("rho_latest_EII"),3),
  ", permutation p=",
  fmtp(cval("perm_p_latest_EII")),
  ") and with Simpson-diversity instability (rho=",
  fmt(cval("rho_Simpson_instability"),3),
  ", permutation p=",
  fmtp(cval("perm_p_Simpson_instability")),
  "), but not clearly with latest Bray-Curtis displacement (rho=",
  fmt(cval("rho_latest_Bray"),3),
  ", p=",
  fmtp(cval("p_latest_Bray")),
  "). ",
  "Exploratory clustering did not support stable trajectory subtypes because the k=2 solution separated nine patients from a single patient. ",
  "These findings support heterogeneous taxonomic routes accompanying a broadly shared ecological displacement process rather than a single conserved genus-level trajectory (Figure 3)."
)

sec4 <- paste0(
  "4. External ecological-state differentiation beyond severe trauma\n\n",
  "In PRJNA1010969, sepsis samples were farther from the healthy-control centroid than severe-trauma samples (",
  fmt(cval("Sepsis_healthy_centroid_distance"),3),
  " vs ",
  fmt(cval("Trauma_healthy_centroid_distance"),3),
  "; difference +",
  fmt(cval("Sepsis_minus_Trauma_centroid_difference"),3),
  "; FDR=",
  fmtp(cval("centroid_difference_FDR")),
  "). ",
  "CLR/Aitchison sensitivity analyses showed stable Sepsis-versus-Trauma separation across pseudocount choices, with PERMANOVA R2 ranging from ",
  fmt(cval("CLR_Aitchison_R2_min"),3),
  " to ",
  fmt(cval("CLR_Aitchison_R2_max"),3),
  " and PERMANOVA FDR from ",
  fmtp(cval("CLR_Aitchison_PERMANOVA_FDR_min")),
  " to ",
  fmtp(cval("CLR_Aitchison_PERMANOVA_FDR_max")),
  ", while dispersion FDR remained non-significant (",
  fmt(cval("CLR_Aitchison_dispersion_FDR_min"),3),
  "-",
  fmt(cval("CLR_Aitchison_dispersion_FDR_max"),3),
  "). ",
  "The external cohort therefore supports sepsis-associated ecological-state differentiation beyond severe trauma, but does not establish a reproducible sepsis-specific taxonomic signature (Figure 4)."
)

sec5 <- paste0(
  "5. Infection source does not clearly modify longitudinal ecological displacement\n\n",
  "In CRA002354, pulmonary versus recorded non-pulmonary infection source did not show a clear difference in the rate of Bray-Curtis displacement from each participant's first available microbiome sample (source-by-time beta=",
  fmt(primary_cra$effect,3),
  ", 95% CI ",
  fmt(primary_cra$ci_low,3),
  " to ",
  fmt(primary_cra$ci_high,3),
  "; p=",
  fmtp(primary_cra$p_value),
  "). ",
  "Estimated source-specific slopes were positive in both groups (pulmonary beta=",
  fmt(pulm$slope_per_day,3),
  ", p=",
  fmtp(pulm$p_LRT),
  "; recorded non-pulmonary beta=",
  fmt(nonp$slope_per_day,3),
  ", p=",
  fmtp(nonp$p_LRT),
  "), although neither was individually significant. ",
  "All prespecified key Bray-Curtis sensitivity analyses retained the same positive interaction direction without statistical significance. ",
  "At baseline, the unrestricted source association was small (PERMANOVA R2=",
  fmt(BASEB$primary_PERMANOVA_R2[1],3),
  ", p=",
  fmtp(BASEB$primary_PERMANOVA_p[1]),
  ") and accompanied by dispersion heterogeneity (p=",
  fmtp(BASEB$primary_dispersion_p[1]),
  "). ",
  "When baseline sampling was restricted to ICU Day 3 or earlier, the association was not reproduced (R2=",
  fmt(BASEB$day3_restricted_PERMANOVA_R2[1],3),
  ", p=",
  fmtp(BASEB$day3_restricted_PERMANOVA_p[1]),
  "; dispersion p=",
  fmtp(BASEB$day3_restricted_dispersion_p[1]),
  "). ",
  "Exploratory genus-level analysis yielded ",
  GENUS$primary_baseline_FDR_hits[1],
  " baseline FDR-significant genera in the primary analysis, but none remained significant at both the >=70% and >=80% sample-level genus-coverage thresholds, and no primary longitudinal genus source-by-time interaction survived FDR correction. ",
  "These findings do not support a strong infection-source modification of longitudinal ecological displacement and do not establish a robust infection-source-specific taxonomic signature (Figure 5)."
)

draft_v3 <- paste(sec1,sec2,sec3,sec4,sec5,sep="\n\n")

writeLines(
  draft_v3,
  file.path(
    DRAFTDIR,
    "Results_Draft_v3_English_MANUSCRIPT.txt"
  )
)

# ------------------------------------------------------------
# 5. Change log
# ------------------------------------------------------------

writeLines(
  c(
    "Results v2 -> v3 manuscript-polish changes",
    "",
    "1. Replaced internal wording such as 'Step87B analysis object' with manuscript-facing language.",
    "2. Changed 'strong positive contrast' to 'large positive contrast' for PRJEB82425 to avoid qualitative overstatement.",
    "3. Replaced 'reproducible' with the more conservative 'observed across several cohorts' in the synthesis sentence.",
    "4. Defined EII at first Results use as an exploratory ecological injury index.",
    "5. Kept the k=2 9-vs-1 clustering result explicitly negative for stable subtype interpretation.",
    "6. Added exact CRA002354 unrestricted and Day-3-restricted PERMANOVA and dispersion statistics.",
    "7. Retained genus-level findings as exploratory and robustness-limited.",
    "8. No frozen numeric result or evidence tier was changed."
  ),
  file.path(
    DRAFTDIR,
    "Results_v2_to_v3_Change_Log.txt"
  )
)

# ------------------------------------------------------------
# 6. Table placement freeze
# ------------------------------------------------------------

table_placement <- tribble(
  ~item, ~recommended_location, ~status,
  "Main Table 1 - Cohort characteristics",
    "MAIN MANUSCRIPT",
    "KEEP",
  "Cross-cohort longitudinal result table",
    "MAIN MANUSCRIPT if journal allows >=2 tables; otherwise Supplementary",
    "KEEP",
  "Evidence architecture table",
    "INTERNAL / optional Supplementary",
    "DO NOT USE AS MAIN TABLE",
  "Compact mixed supporting-result table",
    "INTERNAL DRAFTING AID",
    "DO NOT SUBMIT IN CURRENT FORM"
)

write_csv(
  table_placement,
  file.path(
    AUDDIR,
    "FINAL_TABLE_PLACEMENT.csv"
  )
)

# ------------------------------------------------------------
# 7. Supplementary table plan
# ------------------------------------------------------------

supp_plan <- tribble(
  ~supp_table, ~content, ~frozen_source,
  "Table S1",
    "Full cross-cohort longitudinal primary/common-anchor estimates and evidence classes",
    "Step93W / Step95E2 longitudinal result table",
  "Table S2",
    "Core-sepsis trajectory ecology endpoints and robustness metrics",
    "Step93U",
  "Table S3",
    "External Control-Trauma-Sepsis ecological-state and CLR/Aitchison robustness results",
    "Step93X4",
  "Table S4",
    "CRA002354 primary infection-source ecological models and source-specific slopes",
    "Step94B2D",
  "Table S5",
    "CRA002354 key Bray sensitivity analyses and baseline beta-diversity robustness",
    "Step94B2D",
  "Table S6",
    "Exploratory baseline genus hits and taxonomy-coverage robustness",
    "Step94B3C",
  "Table S7",
    "Sensitivity-only longitudinal genus signals and reporting guardrails",
    "Step94B3C",
  "Table S8",
    "Exploratory EII construction / ecological-consistency analyses",
    "Step92B2 / Step93U"
)

write_csv(
  supp_plan,
  file.path(
    AUDDIR,
    "SUPPLEMENTARY_TABLE_PLAN.csv"
  )
)

# ------------------------------------------------------------
# 8. Guardrail audit
# ------------------------------------------------------------

prohibited <- tribble(
  ~phrase_family, ~should_not_appear,
  "sepsis specificity", "sepsis-specific longitudinal instability",
  "stable trajectory subtype", "stable taxonomic trajectory subtype",
  "clinical EII", "prognostic EII",
  "source equivalence", "pulmonary and non-pulmonary trajectories were identical",
  "source null", "infection source has no effect",
  "Enterococcus biomarker", "Enterococcus is a robust infection-source biomarker",
  "Bacteroides validation", "validated source-specific Bacteroides trajectory"
)

guard_audit <- prohibited %>%
  mutate(
    detected =
      str_detect(
        tolower(draft_v3),
        fixed(tolower(should_not_appear))
      )
  )

write_csv(
  guard_audit,
  file.path(
    AUDDIR,
    "RESULTS_GUARDRAIL_AUDIT.csv"
  )
)

if (any(guard_audit$detected)) {
  stop("A prohibited overclaim phrase was detected in Results Draft v3.")
}

# ------------------------------------------------------------
# 9. Final freeze
# ------------------------------------------------------------

qc <- tibble(
  step95E2_complete = TRUE,
  manuscript_table1_created =
    file.exists(
      file.path(
        TABDIR,
        "Main_Table_1_Cohort_Characteristics_MANUSCRIPT.csv"
      )
    ),
  results_v3_created =
    file.exists(
      file.path(
        DRAFTDIR,
        "Results_Draft_v3_English_MANUSCRIPT.txt"
      )
    ),
  guardrail_audit_pass =
    !any(guard_audit$detected),
  no_new_inferential_analysis = TRUE
) %>%
  mutate(
    ready_for_supplementary_tables_and_discussion =
      manuscript_table1_created &
      results_v3_created &
      guardrail_audit_pass
  )

write_csv(
  qc,
  file.path(
    OUT,
    "STEP95F_READINESS.csv"
  )
)

if (!isTRUE(qc$ready_for_supplementary_tables_and_discussion)) {
  stop("Step95F readiness failed.")
}

writeLines(
  c(
    paste0("Completed: ",Sys.time()),
    "Manuscript Table 1 created: TRUE",
    "Results Draft v3 created: TRUE",
    "Guardrail audit passed: TRUE",
    "No new inferential analysis: TRUE",
    "Ready for Supplementary Tables + Discussion: TRUE",
    "STEP95F COMPLETE"
  ),
  file.path(
    OUT,
    "_STEP95F_COMPLETE.ok"
  )
)

cat("STEP95F COMPLETE\n")
