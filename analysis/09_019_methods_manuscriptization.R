
# ============================================================
# V2_96D METHODS MANUSCRIPTIZATION + PROVENANCE AUDIT
#
# Purpose:
# 1. Audit the local frozen analysis provenance.
# 2. Assemble a manuscript-style Methods draft from verified,
#    frozen analysis decisions.
# 3. Mark only details that still require source-level checking.
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
CODE <- file.path(ROOT,"code","03_data_processing")
RESULTS <- file.path(ROOT,"results")

OUT <- file.path(
  RESULTS,
  "V2_96D_METHODS_MANUSCRIPTIZATION"
)
dir.create(OUT,recursive=TRUE,showWarnings=FALSE)

# ------------------------------------------------------------
# 1. Frozen analysis directories expected from the V2 workflow
# ------------------------------------------------------------
dirs <- tribble(
  ~step, ~role, ~path,
  "27B", "primary analysis metadata / replicate freeze",
  file.path(RESULTS,"V2_27B_ANALYSIS_OBJECTS_AND_REPLICATE_FREEZE"),
  "28A2", "longitudinal diversity and Bray displacement",
  file.path(RESULTS,"V2_28A2_STEP88A_LONGITUDINAL_DIVERSITY_AND_DISPLACEMENT_FIXED"),
  "29A2", "paired taxonomic CLR trajectories",
  file.path(RESULTS,"V2_29A2_STEP89A_TAXONOMIC_PAIRED_TRAJECTORIES_FIXED"),
  "32B2", "exploratory EII harmonization",
  file.path(RESULTS,"V2_32B2_STEP92B_BASELINE_HARMONIZATION_FIXED"),
  "33U", "core trajectory ecology robustness",
  file.path(RESULTS,"V2_33U_TRAJECTORY_ECOLOGY_ROBUSTNESS"),
  "33W", "cross-cohort longitudinal evidence freeze",
  file.path(RESULTS,"V2_33W_LONGITUDINAL_EVIDENCE_FREEZE"),
  "33X4", "external trauma comparison evidence freeze",
  file.path(RESULTS,"V2_33X4_EXTERNAL_VALIDATION_EVIDENCE_FREEZE"),
  "34B2D", "infection-source ecological evidence freeze",
  file.path(RESULTS,"V2_34B2D_INFECTION_SOURCE_EVIDENCE_FREEZE"),
  "34B3C", "infection-source genus evidence freeze",
  file.path(RESULTS,"V2_34B3C_CRA002354_GENUS_EVIDENCE_FREEZE"),
  "35F3", "final Results and Table 1 freeze",
  file.path(RESULTS,"V2_35F3_RESULTS_AND_TABLE1_FINAL_FREEZE_SENTENCE_AUDIT"),
  "36I3", "cross-cohort taxonomic reproducibility",
  file.path(RESULTS,"V2_36I3_CROSS_COHORT_TAXONOMIC_REPRODUCIBILITY_FINAL"),
  "36J", "Microbiome enhancement evidence freeze",
  file.path(RESULTS,"V2_36J_MICROBIOME_ENHANCEMENT_EVIDENCE_FREEZE")
) %>%
  mutate(exists=dir.exists(path))

write_csv(
  dirs,
  file.path(OUT,"01_METHODS_PROVENANCE_DIRECTORY_REGISTRY.csv")
)

critical_missing <- dirs %>%
  filter(
    step %in% c("27B","28A2","29A2","33U","33W","33X4","34B2D","34B3C","36I3","36J"),
    !exists
  )

if(nrow(critical_missing)>0){
  stop(
    "Missing critical frozen analysis directory/directories: ",
    paste(critical_missing$path,collapse=" ; ")
  )
}

# ------------------------------------------------------------
# 2. Build code provenance registry
# ------------------------------------------------------------
code_files <- list.files(
  CODE,
  recursive=TRUE,
  full.names=TRUE,
  pattern="\\.(R|r|ps1|PS1)$"
)

step_patterns <- c(
  "27B","87B",
  "28A","88A",
  "29A","89A",
  "92B",
  "93U","93W","93X",
  "94B",
  "36I","36J"
)

is_relevant <- vapply(
  code_files,
  function(f){
    any(str_detect(
      basename(f),
      regex(paste(step_patterns,collapse="|"),ignore_case=TRUE)
    ))
  },
  logical(1)
)

relevant_code <- code_files[is_relevant]

code_registry <- tibble(
  file=relevant_code,
  filename=basename(relevant_code),
  size_bytes=file.info(relevant_code)$size,
  modified=file.info(relevant_code)$mtime
) %>%
  arrange(filename)

write_csv(
  code_registry,
  file.path(OUT,"02_RELEVANT_ANALYSIS_CODE_REGISTRY.csv")
)

# ------------------------------------------------------------
# 3. Key frozen inputs/outputs
# ------------------------------------------------------------
key_files <- tribble(
  ~item, ~path,
  "primary_analysis_metadata_785",
  file.path(
    RESULTS,
    "V2_27B_ANALYSIS_OBJECTS_AND_REPLICATE_FREEZE",
    "V2_STEP87B_PRIMARY_analysis_metadata_785_patient_time_observations.csv"
  ),
  "step88A2_bray_displacement",
  file.path(
    RESULTS,
    "V2_28A2_STEP88A_LONGITUDINAL_DIVERSITY_AND_DISPLACEMENT_FIXED",
    "03_BETA_DISPLACEMENT",
    "V2_STEP88A2_ALL_within_patient_bray_displacement.csv"
  ),
  "step36I3_pairwise_concordance",
  file.path(
    RESULTS,
    "V2_36I3_CROSS_COHORT_TAXONOMIC_REPRODUCIBILITY_FINAL",
    "03_PRIMARY3_PAIRWISE_CONCORDANCE.csv"
  ),
  "step36I3_threeway_summary",
  file.path(
    RESULTS,
    "V2_36I3_CROSS_COHORT_TAXONOMIC_REPRODUCIBILITY_FINAL",
    "06_THREEWAY_DIRECTION_SUMMARY.csv"
  ),
  "step36J_frozen_metrics",
  file.path(
    RESULTS,
    "V2_36J_MICROBIOME_ENHANCEMENT_EVIDENCE_FREEZE",
    "02_FROZEN_ENHANCEMENT_METRICS.csv"
  )
) %>%
  mutate(exists=file.exists(path))

write_csv(
  key_files,
  file.path(OUT,"03_KEY_FROZEN_FILE_REGISTRY.csv")
)

if(any(!key_files$exists)){
  stop(
    "Missing key frozen file(s): ",
    paste(key_files$path[!key_files$exists],collapse=" ; ")
  )
}

# ------------------------------------------------------------
# 4. Manuscript-style Methods draft
# ------------------------------------------------------------
methods <- c(
"METHODS DRAFT V1 — MANUSCRIPTIZATION FROM FROZEN V2 ANALYSES",
"",
"Study design and analytical framework",
"",
"This study used a multi-cohort secondary analysis of publicly available 16S rRNA gene sequencing datasets to characterize longitudinal gut-microbiome change during sepsis and critical illness. Cohorts were assigned prespecified analytical roles before synthesis, including a core repeated-sepsis cohort, longitudinal ICU-background cohorts, a non-sepsis longitudinal control, intervention-support cohorts, an external Control–Trauma–Sepsis comparison, and an infection-source extension. Because cohorts differed in sequencing region, laboratory workflow, sampling schedule, and upstream feature definition, ASV or OTU abundance matrices were not pooled across studies. Analyses were conducted within cohort and synthesized at the level of ecological effects or harmonized taxonomic effect estimates.",
"",
"Cohort-level data processing and analysis populations",
"",
"Sample-level metadata and microbiome feature tables were harmonized within cohort. The final longitudinal analysis metadata contained 785 patient–time observations across the eight primary/supporting longitudinal datasets. Replicate handling was frozen before longitudinal inference. Where duplicate observations represented the same patient and time point, a single observation was retained according to the prespecified sequencing-depth rule. The external Control–Trauma–Sepsis cohort and the CRA002354 infection-source cohort were analyzed in their dedicated branches rather than merged into the 785-observation longitudinal metadata object.",
"",
"[VERIFY-1] Insert the exact upstream denoising/QC and taxonomic-assignment workflow for each public cohort from the original processing scripts. The manuscript should distinguish ASV-based cohorts from CRA002354, which required a separate 97% OTU rescue workflow.",
"",
"Alpha diversity and longitudinal ecological displacement",
"",
"Within each longitudinal cohort, alpha-diversity metrics included observed richness, Shannon diversity, and Simpson diversity. Longitudinal alpha-diversity effects were evaluated within cohort using the prespecified repeated-measures framework. Ecological displacement was quantified as the Bray–Curtis distance between each longitudinal sample and that patient's baseline or first available microbiome sample. Primary early-to-late contrasts and a common-anchor sensitivity analysis were retained as distinct, prespecified summaries of longitudinal change. Cohort-specific estimates were interpreted within cohort rather than by directly pooling microbiome feature matrices across studies.",
"",
"[VERIFY-2] Confirm the exact fixed/random effects and covariate specification used in the Step88A2 longitudinal mixed models from the frozen Step88A2 script before final Methods submission.",
"",
"Core-sepsis taxonomic trajectory analysis",
"",
"In the core repeated-sepsis cohort (PRJNA691455), longitudinal taxonomic change was represented as patient-specific signed genus-level trajectories. For each patient with the required repeated samples, taxonomic changes across the prespecified time points were used to characterize the direction and magnitude of compositional change. Pairwise Euclidean distances and cosine similarities between patient trajectory vectors were used to quantify between-patient taxonomic heterogeneity. These trajectory-level measures were evaluated alongside ecological displacement and diversity instability. Unsupervised clustering was explored only as a descriptive sensitivity analysis; because the resulting two-cluster solution was highly imbalanced, no stable trajectory subtype interpretation was retained.",
"",
"Exploratory ecological injury index",
"",
"An exploratory ecological injury index (EII) was constructed as an internal summary of microbiome disturbance by combining standardized ecological displacement and absolute diversity-change components and rescaling the resulting score to a 0–100 range. The EII was used only as a supportive ecological-state descriptor. It was not developed or validated as a clinical prediction or prognostic score.",
"",
"External ecological-state validation against severe trauma",
"",
"PRJNA1010969 was analyzed as an external cross-sectional ecological-state comparison comprising Control, Trauma, and Sepsis groups. Alpha-diversity and community-level separation were evaluated using the cohort's harmonized analysis population. To assess whether Sepsis–Trauma separation was robust to compositional analysis and dispersion effects, CLR/Aitchison analyses were repeated across prespecified pseudocount choices and accompanied by dispersion testing. This branch was interpreted as external ecological-state differentiation beyond severe trauma, not as independent longitudinal replication or validation of a universal taxonomic signature.",
"",
"CRA002354 infection-source extension",
"",
"CRA002354 was analyzed as a dedicated infection-source extension. The public dataset was distributed as FASTA sequences without per-base quality scores and therefore could not be reprocessed using an ASV-denoising workflow equivalent to datasets with raw quality information. Recoverable sequence-level quality filters were applied, followed by chimera removal, 97% OTU reconstruction with VSEARCH, and SILVA 138.2-based taxonomic assignment. The primary analysis population required at least 4,000 non-chimeric reads and included 130 samples from 64 participants.",
"",
"Pulmonary and recorded non-pulmonary sepsis were compared using longitudinal within-patient Bray–Curtis displacement from each patient's first available microbiome sample. The primary model tested the source-by-time interaction, with prespecified sensitivity analyses incorporating available clinical variables, exclusion of OTHER_UNKNOWN source assignments, rarefaction, and an alternative minimum-depth threshold. Alpha-diversity source-by-time interactions were evaluated in parallel. Genus-level source analyses were exploratory and were subjected to FDR correction and additional taxonomic-coverage sensitivity analyses.",
"",
"[VERIFY-3] Confirm the exact CRA002354 model family, random-effects specification, and likelihood-ratio-test implementation from the frozen Step94B2 scripts before final Methods submission.",
"",
"Cross-cohort taxonomic reproducibility",
"",
"To quantify the reproducibility of longitudinal taxonomic effects across cohorts without pooling heterogeneous feature matrices, paired taxonomic effects were first estimated independently within each cohort using the frozen Step89A paired CLR analysis. Genus- and family-level effects were analyzed separately. The primary cross-cohort analysis used the ALL_PAIRED definition, while COMMON_ANCHOR_SENSITIVITY was retained as a prespecified sensitivity analysis.",
"",
"Taxon labels were harmonized only after within-cohort effect estimation. For each pair of the three primary progressive cohorts (PRJNA691455, PRJNA851469, and PRJNA516701), reproducibility was summarized using the Spearman correlation between shared-taxon CLR effect vectors, the proportion of shared taxa changing in the same direction, and the Jaccard overlap among taxa with the largest absolute effects. Null distributions for the correlation, directional agreement, and top-effect overlap were generated using 5,000 permutations of one cohort's effect vector. Analyses were repeated at genus and family levels and were supplemented by a three-cohort analysis of taxa shared across all three datasets. The common-anchor analysis served as a sensitivity analysis rather than an independent replication.",
"",
"Statistical interpretation and multiplicity",
"",
"Statistical inference was conducted within each prespecified analytical branch. Where multiple alpha-diversity metrics or taxonomic features were evaluated within a defined family of tests, false-discovery-rate correction was applied. Permutation testing was used for the cross-cohort taxonomic reproducibility analysis. Statistical significance was not used as the sole criterion for evidence classification; direction, sensitivity analyses, dispersion diagnostics, taxonomic-coverage robustness, and prespecified cohort roles were considered when determining the manuscript interpretation.",
"",
"Software and reproducibility",
"",
"Data processing and statistical analyses were performed in R (version 4.4.0) with external command-line tools used where required for sequence processing. VSEARCH was used for the CRA002354 97% OTU reconstruction. All manuscript-level results were generated from frozen cohort-specific analysis objects and evidence-freeze outputs. No ASV or OTU feature matrix was directly merged across heterogeneous cohorts.",
"",
"[VERIFY-4] Populate the exact R package names/versions and VSEARCH version from the local analysis environment or frozen session information before submission."
)

METHODS_FILE <- file.path(
  OUT,
  "Methods_Draft_v1_FROM_FROZEN_V2_ANALYSES.txt"
)

writeLines(
  methods,
  METHODS_FILE,
  useBytes=TRUE
)

# ------------------------------------------------------------
# 5. Verification checklist
# ------------------------------------------------------------
verify <- tibble(
  id=c("VERIFY-1","VERIFY-2","VERIFY-3","VERIFY-4"),
  priority=c("HIGH","HIGH","HIGH","MEDIUM"),
  item=c(
    "Exact upstream cohort-specific QC/denoising/taxonomy workflow; distinguish ASV cohorts from CRA002354 OTU rescue.",
    "Exact Step88A2 mixed-model fixed/random effects and any covariates.",
    "Exact Step94B2 CRA002354 longitudinal model family/random effect/LRT implementation.",
    "Exact R package versions and VSEARCH version."
  ),
  why_needed=c(
    "Required for reproducible microbiome Methods and to avoid claiming identical preprocessing across heterogeneous public datasets.",
    "Needed for precise statistical Methods wording.",
    "Needed for precise infection-source Methods wording.",
    "Required for final software/reproducibility statement."
  )
)

write_csv(
  verify,
  file.path(OUT,"04_METHODS_ITEMS_REQUIRING_VERIFICATION.csv")
)

# ------------------------------------------------------------
# 6. Summary
# ------------------------------------------------------------
writeLines(
  c(
    "V2 STEP96D METHODS MANUSCRIPTIZATION",
    "",
    paste0("Methods draft: ",METHODS_FILE),
    "",
    "The draft was assembled from frozen V2 analytical decisions.",
    "No statistical model was rerun.",
    "No existing frozen result was overwritten.",
    "",
    "Four verification items remain and are explicitly listed in 04_METHODS_ITEMS_REQUIRING_VERIFICATION.csv.",
    "",
    "STEP96D COMPLETE"
  ),
  file.path(OUT,"05_STEP96D_SUMMARY.txt")
)

writeLines(
  c(
    paste0("Completed: ",Sys.time()),
    "STEP96D COMPLETE",
    "Methods draft v1 created.",
    "No new statistical inference."
  ),
  file.path(OUT,"_STEP96D_COMPLETE.txt")
)

cat("STEP96D COMPLETE\n")
cat("Methods:",METHODS_FILE,"\n")
