
# ============================================================
# V2_96C DISCUSSION INTEGRATION
#
# Purpose:
# Integrate the frozen Step36I3/Step36J result into the existing
# evidence-grounded Discussion v1.
#
# NO statistical model is rerun.
# NO frozen result is overwritten.
# Literature placeholders inserted here have a verified map
# supplied in this package.
# ============================================================

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(stringr)
})

ROOT <- "E:/sepsis_project"
RESULTS <- file.path(ROOT,"results")

DISCUSSION_DIR <- file.path(
  RESULTS,
  "V2_35G_SUPPLEMENTARY_TABLES_AND_DISCUSSION_V1"
)

# Locate discussion v1 robustly.
candidate_files <- list.files(
  DISCUSSION_DIR,
  pattern="Discussion.*Draft.*v1.*\\.txt$",
  full.names=TRUE,
  recursive=TRUE,
  ignore.case=TRUE
)

if(length(candidate_files)==0){
  candidate_files <- list.files(
    DISCUSSION_DIR,
    pattern="Discussion.*\\.txt$",
    full.names=TRUE,
    recursive=TRUE,
    ignore.case=TRUE
  )
}

if(length(candidate_files)==0){
  stop("Could not locate Discussion v1 under: ",DISCUSSION_DIR)
}

# Prefer exact expected name if present.
exact_expected <- file.path(
  DISCUSSION_DIR,
  "Discussion_Draft_v1_Evidence_Grounded.txt"
)

if(file.exists(exact_expected)){
  DISCUSSION_V1 <- exact_expected
}else{
  DISCUSSION_V1 <- candidate_files[[1]]
}

FREEZE_DIR <- file.path(
  RESULTS,
  "V2_36J_MICROBIOME_ENHANCEMENT_EVIDENCE_FREEZE"
)

METRICS_FILE <- file.path(
  FREEZE_DIR,
  "02_FROZEN_ENHANCEMENT_METRICS.csv"
)

if(!file.exists(METRICS_FILE)){
  stop("Missing Step36J frozen metrics: ",METRICS_FILE)
}

OUT <- file.path(
  RESULTS,
  "V2_96C_DISCUSSION_INTEGRATION"
)
dir.create(OUT,recursive=TRUE,showWarnings=FALSE)

m <- read_csv(METRICS_FILE,show_col_types=FALSE)

gm <- function(name){
  z <- m %>% filter(metric==name)
  if(nrow(z)!=1) stop("Metric missing/duplicated: ",name)
  as.numeric(z$value[[1]])
}

g_rho <- gm("genus_median_pairwise_spearman_rho")
f_rho <- gm("family_median_pairwise_spearman_rho")
g_same <- gm("genus_threeway_same_direction_proportion")
f_same <- gm("family_threeway_same_direction_proportion")

lines <- readLines(
  DISCUSSION_V1,
  warn=FALSE,
  encoding="UTF-8"
)

# Prevent duplicate insertion.
if(any(str_detect(
  tolower(lines),
  "cross-cohort taxonomic effect vectors"
))){
  stop("Discussion source already appears to contain the Step96C insertion.")
}

# Preferred location:
# before the existing external-validation discussion section.
idx <- which(
  str_detect(
    tolower(lines),
    "external validation supports ecological-state differentiation"
  )
)

if(length(idx)==0){
  idx <- which(
    str_detect(
      tolower(lines),
      "external validation"
    )
  )
}

if(length(idx)!=1){

  heading_candidates <- tibble(
    line_number=seq_along(lines),
    text=lines
  ) %>%
    filter(
      str_detect(
        tolower(text),
        "ecological|taxonomic|external|infection source|implications|limitations"
      )
    )

  write_csv(
    heading_candidates,
    file.path(
      OUT,
      "00_DISCUSSION_HEADING_CANDIDATES.csv"
    )
  )

  stop(
    "Could not uniquely locate external-validation Discussion section. ",
    "See 00_DISCUSSION_HEADING_CANDIDATES.csv"
  )
}

insert_before <- idx[[1]]

new_block <- c(
  "",
  "Cross-cohort ecological concordance and taxonomic heterogeneity",
  "",
  paste0(
    "The cross-cohort analysis further separated reproducibility at the ecosystem level from reproducibility at the level of individual taxa. ",
    "Across the three primary cohorts with progressive ecological displacement, the direction of ecosystem-level change was concordant, whereas longitudinal paired CLR taxonomic effect vectors showed only low-to-moderate concordance (median pairwise Spearman rho=",
    sprintf("%.2f",g_rho),
    " at genus level and ",
    sprintf("%.2f",f_rho),
    " at family level). ",
    "Among taxa shared across all three cohorts, only ",
    sprintf("%.1f",100*g_same),
    "% of genera and ",
    sprintf("%.1f",100*f_same),
    "% of families changed in the same direction in all three cohorts. ",
    "This pattern is compatible with a model in which critically ill patients undergo a broadly shared direction of ecological disturbance through heterogeneous and cohort-dependent taxonomic routes."
  ),
  "",
  paste0(
    "This interpretation is consistent with longitudinal studies showing progressive microbiome disruption during critical illness. ",
    "Schlechte et al. reported dynamic and progressive intestinal dysbiosis during the first week of critical illness, including Enterobacteriaceae enrichment and links to nosocomial infection [REF_SCHLECHTE_2023]. ",
    "In a larger multicompartment study of mechanically ventilated patients with acute respiratory failure, Kitsios et al. likewise observed progressive dysbiosis across gut, oral, and lung communities and associations between microbial trajectories and clinical exposures [REF_KITSIOS_2024]. ",
    "A recent systematic review of 36 longitudinal sequencing studies involving 2,067 critically ill adults found that progressive diversity loss and temporal beta-diversity shifts were common, but also emphasized substantial between-study heterogeneity and very low certainty for many specific microbiome outcomes [REF_THEOCHARIDOU_2026]."
  ),
  "",
  paste0(
    "The limited reproducibility of individual taxonomic effects should therefore not be interpreted as an absence of shared biology. ",
    "Rather, it highlights the difference between a community-level ecological response and the specific taxa through which that response is realized. ",
    "Cross-study microbiome research has shown that the consistency of differential-abundance signatures varies markedly across conditions and cohorts [REF_GEISTLINGER_2024], while experimental and analytical differences in 16S workflows can materially reduce cross-laboratory replicability of taxonomic signals [REF_CLAUSEN_2022]. ",
    "For this reason, the present study estimated taxonomic changes within each cohort and harmonized taxon labels only at the effect level, instead of pooling ASV or OTU abundance matrices across heterogeneous datasets."
  ),
  "",
  "Importantly, the contrast between ecological and taxonomic reproducibility should be interpreted across biological resolutions rather than as a formal statistical comparison of directly commensurable metrics. The data support concordant ecosystem-level direction accompanied by incomplete, cohort-dependent taxonomic concordance; they do not establish a universal sepsis-specific taxonomic signature.",
  ""
)

before <- if(insert_before>1) lines[1:(insert_before-1)] else character()
after <- lines[insert_before:length(lines)]

new_lines <- c(
  before,
  new_block,
  after
)

OUT_DISC <- file.path(
  OUT,
  "Discussion_Draft_v2_WITH_CROSS_COHORT_REPRODUCIBILITY.txt"
)

writeLines(
  new_lines,
  OUT_DISC,
  useBytes=TRUE
)

audit <- tibble(
  source_discussion=DISCUSSION_V1,
  insertion_before_original_line=insert_before,
  original_line_count=length(lines),
  new_line_count=length(new_lines),
  statistical_inference_rerun=FALSE,
  source_discussion_overwritten=FALSE,
  literature_placeholders_verified_map_supplied=TRUE
)

write_csv(
  audit,
  file.path(OUT,"01_DISCUSSION_INTEGRATION_AUDIT.csv")
)

writeLines(
  c(
    "V2 STEP96C DISCUSSION INTEGRATION",
    "",
    paste0("Source Discussion: ",DISCUSSION_V1),
    paste0("Output Discussion: ",OUT_DISC),
    "",
    "Inserted the cross-cohort ecological-vs-taxonomic interpretation before the external-validation section.",
    "Literature tags use the verified literature map supplied with this package.",
    "",
    "No statistical model was rerun.",
    "The original Discussion v1 remains unchanged.",
    "",
    "STEP96C COMPLETE"
  ),
  file.path(OUT,"03_STEP96C_SUMMARY.txt")
)

writeLines(
  c(
    paste0("Completed: ",Sys.time()),
    "STEP96C COMPLETE",
    "Discussion v2 created.",
    "No new statistical inference."
  ),
  file.path(OUT,"_STEP96C_COMPLETE.txt")
)

cat("STEP96C COMPLETE\n")
cat("Output:",OUT_DISC,"\n")
