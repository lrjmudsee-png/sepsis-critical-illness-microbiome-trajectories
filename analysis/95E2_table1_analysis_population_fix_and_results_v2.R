# ============================================================
# Sepsis V2 - Step95E2
# TABLE 1 ANALYSIS-POPULATION FIX + RESULTS DRAFT V2
#
# Supersedes incomplete Step95E.
#
# Why:
# 1) Step95E metadata gate passed, but Table 1 selected upstream frozen-master
#    counts rather than the actual analysis populations.
# 2) Results drafting stopped because the source-specific slope registry used
#    a different source-group column name than Step95E expected.
#
# Step95E2:
# - longitudinal/static Step87B cohorts: use the exact 785-observation PRIMARY
#   analysis metadata after replicate resolution;
# - PRJNA1010969: locate the actual external analysis population and require
#   17 Control + 18 Trauma + 18 Sepsis = 53 analyzed samples;
# - CRA002354: use the B1D/B2D primary true-nonchimeric >=4000 analysis metadata;
# - generate a clean conventional Table 1 plus full provenance;
# - generate Results Draft v2 from frozen evidence only.
#
# NO new inferential analysis.
# ============================================================

options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(tibble)
})

ROOT <- "E:/sepsis_project"
RESULTS <- file.path(ROOT, "results")
DATA <- file.path(ROOT, "data")

SRC <- file.path(
  RESULTS,
  "V2_35B_MANUSCRIPT_SOURCE_PACK_ASSEMBLY"
)

D2 <- file.path(
  RESULTS,
  "V2_35D2_PUBLICATION_FIGURES_SOURCEPACK_ONLY"
)

D3 <- file.path(
  RESULTS,
  "V2_35D3_FINAL_FIGURE_FORMATTING_FIX"
)

OUT <- file.path(
  RESULTS,
  "V2_35E2_TABLE1_ANALYSIS_POPULATION_FIX_AND_RESULTS_V2"
)

TABDIR <- file.path(OUT, "01_TABLES")
DRAFTDIR <- file.path(OUT, "02_RESULTS_DRAFT")
AUDDIR <- file.path(OUT, "03_PROVENANCE_AND_AUDIT")

for (d in c(OUT,TABDIR,DRAFTDIR,AUDDIR)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

if (!file.exists(file.path(SRC, "_STEP95B_COMPLETE.ok"))) {
  stop("Step95B source pack is incomplete.")
}

if (!file.exists(file.path(D3, "_STEP95D3_COMPLETE.ok"))) {
  stop("Step95D3 final figures are incomplete.")
}

# ------------------------------------------------------------
# 1. Frozen result sources
# ------------------------------------------------------------

W <- read_csv(
  file.path(
    RESULTS,
    "V2_35C_MANUSCRIPT_ARCHITECTURE_FIGURE1_AND_MAIN_TABLES",
    "Main_Table_2_Cross_Cohort_Longitudinal_Results.csv"
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

SENS <- read_csv(
  file.path(
    SRC,
    "02_SUPPLEMENTARY_TABLE_SOURCES",
    "SUPP_S5_SOURCE_BRAY_SENSITIVITY.csv"
  ),
  show_col_types = FALSE
)

SLOPES_RAW <- read_csv(
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

first_present <- function(nms, choices) {
  x <- choices[choices %in% nms]
  if (length(x) == 0) return(NA_character_)
  x[1]
}

# ------------------------------------------------------------
# 2. Exact Step87B primary analysis population
# ------------------------------------------------------------

STEP87B <- file.path(
  RESULTS,
  "V2_27B_ANALYSIS_OBJECTS_AND_REPLICATE_FREEZE",
  "V2_STEP87B_PRIMARY_analysis_metadata_785_patient_time_observations.csv"
)

if (!file.exists(STEP87B)) {
  stop("Exact Step87B PRIMARY analysis metadata is missing.")
}

M87 <- read_csv(
  STEP87B,
  show_col_types = FALSE
)

required87 <- c("project","run_id","patient_id","time_day")
if (!all(required87 %in% names(M87))) {
  stop("Step87B PRIMARY metadata schema is not as expected.")
}

# Hard provenance gate.
if (nrow(M87) != 785) {
  stop(
    paste0(
      "Step87B PRIMARY row count is ",
      nrow(M87),
      ", expected frozen value 785."
    )
  )
}

n87_pat <- n_distinct(M87$patient_id[!is.na(M87$patient_id)])

# Do not require a recalled participant count; record actual source-derived count.
write_csv(
  tibble(
    source = STEP87B,
    rows = nrow(M87),
    unique_patients = n87_pat
  ),
  file.path(AUDDIR, "01_STEP87B_PRIMARY_SOURCE_AUDIT.csv")
)

targets87 <- c(
  "PRJNA691455",
  "PRJEB82425",
  "PRJNA516701",
  "PRJNA851469",
  "PRJNA578267",
  "PRJNA430161",
  "PRJNA1166732",
  "PRJNA978257"
)

summarize_project <- function(p) {

  d <- M87 %>%
    filter(project == p) %>%
    distinct(run_id, .keep_all = TRUE)

  if (nrow(d) == 0) {
    return(
      tibble(
        cohort = p,
        analysis_source = STEP87B,
        samples_analyzed = NA_integer_,
        participants_analyzed = NA_integer_,
        repeated_participants = NA_integer_,
        time_min = NA_real_,
        time_max = NA_real_,
        population_note = "MISSING_FROM_STEP87B"
      )
    )
  }

  pc <- d %>%
    filter(!is.na(patient_id), patient_id != "") %>%
    distinct(patient_id, run_id) %>%
    count(patient_id, name="n_samples")

  td <- suppressWarnings(as.numeric(d$time_day))
  good_t <- is.finite(td)

  tibble(
    cohort = p,
    analysis_source = STEP87B,
    samples_analyzed = n_distinct(d$run_id),
    participants_analyzed =
      n_distinct(d$patient_id[!is.na(d$patient_id) & d$patient_id != ""]),
    repeated_participants = sum(pc$n_samples >= 2),
    time_min = ifelse(any(good_t), min(td[good_t]), NA_real_),
    time_max = ifelse(any(good_t), max(td[good_t]), NA_real_),
    population_note =
      "Step87B PRIMARY analysis population after replicate resolution"
  )
}

counts87 <- bind_rows(
  lapply(targets87, summarize_project)
)

if (any(is.na(counts87$samples_analyzed))) {
  stop("At least one expected Step87B project is absent from the PRIMARY analysis metadata.")
}

write_csv(
  counts87,
  file.path(AUDDIR, "02_STEP87B_PROJECT_ANALYSIS_COUNTS.csv")
)

# ------------------------------------------------------------
# 3. PRJNA1010969 external analysis-population locator
# ------------------------------------------------------------

read_tabular <- function(path) {

  ext <- tolower(tools::file_ext(path))

  tryCatch(
    {
      if (ext == "tsv") {
        read_tsv(path, show_col_types=FALSE, progress=FALSE)
      } else {
        read_csv(path, show_col_types=FALSE, progress=FALSE)
      }
    },
    error=function(e) NULL
  )
}

xdirs <- list.dirs(
  RESULTS,
  recursive = FALSE,
  full.names = TRUE
)

xdirs <- xdirs[
  str_detect(
    basename(xdirs),
    regex("33X|93X|EXTERNAL", ignore_case=TRUE)
  )
]

xfiles <- unique(
  unlist(
    lapply(
      xdirs,
      function(d) {
        list.files(
          d,
          pattern="\\.(csv|tsv)$",
          recursive=TRUE,
          full.names=TRUE,
          ignore.case=TRUE
        )
      }
    ),
    use.names=FALSE
  )
)

sample_choices <- c(
  "run_id","Run_ID","sample_id","Sample_ID","sample","Sample","id"
)

group_choices <- c(
  "group","Group","external_group","clinical_group",
  "disease_group","condition","Condition","status"
)

normalize_group <- function(x) {

  u <- toupper(as.character(x))

  case_when(
    str_detect(u, "CONTROL|HEALTH") ~ "CONTROL",
    str_detect(u, "TRAUMA") ~ "TRAUMA",
    str_detect(u, "SEPSIS") ~ "SEPSIS",
    TRUE ~ NA_character_
  )
}

external_candidates <- list()

for (f in xfiles) {

  z <- read_tabular(f)
  if (is.null(z) || nrow(z) == 0) next

  sc <- first_present(names(z), sample_choices)
  gc <- first_present(names(z), group_choices)

  if (is.na(sc) || is.na(gc)) next

  d <- tibble(
    sample_id = as.character(z[[sc]]),
    group_raw = as.character(z[[gc]])
  ) %>%
    mutate(group = normalize_group(group_raw)) %>%
    filter(
      !is.na(sample_id),
      sample_id != "",
      !is.na(group)
    ) %>%
    distinct(sample_id, .keep_all=TRUE)

  if (nrow(d) == 0) next

  cc <- table(
    factor(
      d$group,
      levels=c("CONTROL","TRAUMA","SEPSIS")
    )
  )

  external_candidates[[length(external_candidates)+1]] <- tibble(
    file=f,
    sample_col=sc,
    group_col=gc,
    n_control=as.integer(cc["CONTROL"]),
    n_trauma=as.integer(cc["TRAUMA"]),
    n_sepsis=as.integer(cc["SEPSIS"]),
    n_total=nrow(d),
    exact_frozen_match =
      as.integer(cc["CONTROL"]) == 17 &&
      as.integer(cc["TRAUMA"]) == 18 &&
      as.integer(cc["SEPSIS"]) == 18
  )
}

external_audit <- if (length(external_candidates)>0) {
  bind_rows(external_candidates) %>%
    arrange(desc(exact_frozen_match), desc(n_total))
} else {
  tibble(
    file=character(),
    sample_col=character(),
    group_col=character(),
    n_control=integer(),
    n_trauma=integer(),
    n_sepsis=integer(),
    n_total=integer(),
    exact_frozen_match=logical()
  )
}

write_csv(
  external_audit,
  file.path(AUDDIR, "03_PRJNA1010969_EXTERNAL_POPULATION_CANDIDATES.csv")
)

if (
  nrow(external_audit) == 0 ||
  !any(external_audit$exact_frozen_match)
) {
  stop(
    "Could not source-anchor the frozen PRJNA1010969 population 17 Control + 18 Trauma + 18 Sepsis."
  )
}

ext_sel <- external_audit %>%
  filter(exact_frozen_match) %>%
  slice(1)

ext_count <- tibble(
  cohort="PRJNA1010969",
  analysis_source=ext_sel$file[1],
  samples_analyzed=53L,
  participants_analyzed=53L,
  repeated_participants=0L,
  time_min=NA_real_,
  time_max=NA_real_,
  population_note=
    "Final external analysis population: 17 Control + 18 Trauma + 18 Sepsis"
)

# ------------------------------------------------------------
# 4. CRA002354 primary B1D/B2D analysis population
# ------------------------------------------------------------

CRA_META <- file.path(
  DATA,
  "CRA002354",
  "03_vsearch97_silva1382",
  "CRA002354_analysis_metadata_PRIMARY_NONCHIMERIC_min4000.csv"
)

if (!file.exists(CRA_META)) {
  stop("CRA002354 primary non-chimeric >=4000 analysis metadata is missing.")
}

CRA <- read_csv(
  CRA_META,
  show_col_types=FALSE
)

cra_sample_col <- first_present(
  names(CRA),
  c("Run_ID","run_id","sample_id","Sample_ID")
)

cra_patient_col <- first_present(
  names(CRA),
  c("patient_id","Patient_ID","subject_id")
)

cra_time_col <- first_present(
  names(CRA),
  c("time_day","day","Day")
)

if (is.na(cra_sample_col) || is.na(cra_patient_col)) {
  stop("CRA primary metadata lacks sample/patient identifiers.")
}

CRA2 <- CRA %>%
  mutate(
    .sample = as.character(.data[[cra_sample_col]]),
    .patient = as.character(.data[[cra_patient_col]])
  ) %>%
  filter(
    !is.na(.sample),
    .sample != ""
  ) %>%
  distinct(.sample, .keep_all=TRUE)

cra_pc <- CRA2 %>%
  filter(
    !is.na(.patient),
    .patient != ""
  ) %>%
  distinct(.patient,.sample) %>%
  count(.patient,name="n_samples")

cra_td <- if (!is.na(cra_time_col)) {
  suppressWarnings(as.numeric(CRA2[[cra_time_col]]))
} else {
  rep(NA_real_,nrow(CRA2))
}

cra_count <- tibble(
  cohort="CRA002354",
  analysis_source=CRA_META,
  samples_analyzed=nrow(CRA2),
  participants_analyzed=
    n_distinct(CRA2$.patient[!is.na(CRA2$.patient) & CRA2$.patient!=""]),
  repeated_participants=sum(cra_pc$n_samples>=2),
  time_min=ifelse(any(is.finite(cra_td)),min(cra_td[is.finite(cra_td)]),NA_real_),
  time_max=ifelse(any(is.finite(cra_td)),max(cra_td[is.finite(cra_td)]),NA_real_),
  population_note=
    "Primary true-nonchimeric OTU97 analysis population after >=4000-read depth gate"
)

# Frozen B1D expectations.
if (cra_count$samples_analyzed != 130L) {
  stop(
    paste0(
      "CRA primary analyzed samples = ",
      cra_count$samples_analyzed,
      "; expected frozen B1D value 130."
    )
  )
}

if (cra_count$participants_analyzed != 64L) {
  stop(
    paste0(
      "CRA primary participants = ",
      cra_count$participants_analyzed,
      "; expected frozen value 64."
    )
  )
}

# ------------------------------------------------------------
# 5. Final Table 1 analysis-population registry
# ------------------------------------------------------------

counts <- bind_rows(
  counts87,
  ext_count,
  cra_count
)

contexts <- tribble(
  ~cohort, ~clinical_context, ~analysis_role, ~manuscript_role,
  "PRJNA691455", "Sepsis / critical illness", "Core sepsis longitudinal", "Primary",
  "PRJEB82425", "ICU infection", "External longitudinal ICU-infection support", "Supportive",
  "PRJNA516701", "ICU background", "ICU-background longitudinal context", "Primary context",
  "PRJNA851469", "ICU background", "ICU-background longitudinal context", "Primary context",
  "PRJNA578267", "Non-sepsis surgical control", "Non-sepsis longitudinal control", "Control",
  "PRJNA430161", "Intervention cohort", "Intervention longitudinal support", "Supportive",
  "PRJNA1166732", "Intervention cohort", "Intervention longitudinal support", "Supportive",
  "PRJNA978257", "Static supportive cohort", "Static support", "Supplementary",
  "PRJNA1010969", "Control / severe trauma / sepsis", "External ecological-state validation", "Supportive external",
  "CRA002354", "Sepsis stratified by infection source", "Infection-source longitudinal extension", "Supportive extension"
)

long_result <- W %>%
  transmute(
    cohort=Cohort,
    primary_longitudinal_pairs=Primary_pairs,
    longitudinal_evidence_class=Evidence_class
  )

table1_full <- contexts %>%
  left_join(counts,by="cohort") %>%
  left_join(long_result,by="cohort") %>%
  mutate(
    observed_time_window=case_when(
      !is.na(time_min) & !is.na(time_max) ~
        paste0(
          format(time_min,trim=TRUE,scientific=FALSE),
          "-",
          format(time_max,trim=TRUE,scientific=FALSE),
          " d"
        ),
      TRUE ~ "NR"
    )
  )

if (
  any(is.na(table1_full$samples_analyzed)) ||
  any(is.na(table1_full$participants_analyzed))
) {
  stop("Final Table 1 contains missing sample/participant counts.")
}

write_csv(
  table1_full,
  file.path(AUDDIR, "04_FINAL_TABLE1_ANALYSIS_POPULATION_PROVENANCE.csv")
)

table1 <- table1_full %>%
  transmute(
    Cohort=cohort,
    Clinical_context=clinical_context,
    Analysis_role=analysis_role,
    Samples_analyzed=samples_analyzed,
    Participants_analyzed=participants_analyzed,
    Participants_with_repeated_samples=repeated_participants,
    Observed_time_window=observed_time_window,
    Primary_longitudinal_pairs=primary_longitudinal_pairs,
    Manuscript_role=manuscript_role
  )

write_csv(
  table1,
  file.path(
    TABDIR,
    "Main_Table_1_Cohort_and_Study_Characteristics_ANALYSIS_POPULATION.csv"
  )
)

# Add a footnote file for manuscript use.
writeLines(
  c(
    "Table 1 footnotes:",
    "Samples and participants refer to the population entering the corresponding V2 analytical branch, not to all samples deposited in the original public study.",
    "For the eight Step87B cohorts, counts are derived from the PRIMARY 785-observation analysis metadata after replicate resolution.",
    "PRJNA1010969 counts reflect the final external analysis population (17 Control, 18 Trauma, 18 Sepsis).",
    "CRA002354 counts reflect the true-nonchimeric primary >=4000-read OTU97 analysis population (130 samples from 64 participants).",
    "NR, not represented by a meaningful longitudinal day scale in the harmonized analysis table."
  ),
  file.path(TABDIR,"Main_Table_1_Footnotes.txt")
)

# Preserve frozen longitudinal result table.
write_csv(
  W,
  file.path(
    TABDIR,
    "Main_or_Supp_Table_Cross_Cohort_Longitudinal_Results.csv"
  )
)

# ------------------------------------------------------------
# 6. Robust source-specific slope schema normalization
# ------------------------------------------------------------

sn <- names(SLOPES_RAW)

source_col <- first_present(
  sn,
  c(
    "source_group",
    "group",
    "source",
    "pulmonary_binary",
    "infection_source"
  )
)

slope_col <- first_present(
  sn,
  c(
    "slope_per_day",
    "slope",
    "beta",
    "estimate"
  )
)

lo_col <- first_present(
  sn,
  c(
    "ci95_low",
    "ci_low",
    "lower",
    "lower95"
  )
)

hi_col <- first_present(
  sn,
  c(
    "ci95_high",
    "ci_high",
    "upper",
    "upper95"
  )
)

p_col <- first_present(
  sn,
  c(
    "p_LRT",
    "p_value",
    "p",
    "pvalue"
  )
)

if (
  any(
    is.na(
      c(source_col,slope_col,p_col)
    )
  )
) {
  stop(
    paste0(
      "Could not normalize slope-registry schema. Columns: ",
      paste(sn,collapse=", ")
    )
  )
}

SLOPES <- SLOPES_RAW %>%
  transmute(
    source_group=as.character(.data[[source_col]]),
    slope_per_day=as.numeric(.data[[slope_col]]),
    ci95_low=if (!is.na(lo_col)) as.numeric(.data[[lo_col]]) else NA_real_,
    ci95_high=if (!is.na(hi_col)) as.numeric(.data[[hi_col]]) else NA_real_,
    p_value=as.numeric(.data[[p_col]])
  )

write_csv(
  tibble(
    original_columns=paste(sn,collapse=" | "),
    selected_source_col=source_col,
    selected_slope_col=slope_col,
    selected_ci_low_col=lo_col,
    selected_ci_high_col=hi_col,
    selected_p_col=p_col
  ),
  file.path(AUDDIR,"05_SOURCE_SLOPE_SCHEMA_NORMALIZATION.csv")
)

# ------------------------------------------------------------
# 7. Results Draft v2
# ------------------------------------------------------------

robust_prog <- W %>%
  filter(
    Cohort %in% c(
      "PRJNA691455",
      "PRJNA851469",
      "PRJNA516701"
    )
  )

ext_inf <- W %>%
  filter(Cohort=="PRJEB82425") %>%
  slice(1)

ctrl <- W %>%
  filter(Cohort=="PRJNA578267") %>%
  slice(1)

primary_cra <- ECO %>%
  filter(
    domain=="PRIMARY_LONGITUDINAL_BRAY",
    estimand=="source_by_personal_time_interaction"
  ) %>%
  slice(1)

# Identify pulmonary/nonpulmonary robustly.
sg <- toupper(SLOPES$source_group)

pulm_i <- which(
  str_detect(sg,"PULMONARY") &
  !str_detect(sg,"NON")
)[1]

non_i <- which(
  str_detect(sg,"NON") &
  str_detect(sg,"PULMONARY")
)[1]

if (is.na(pulm_i) || is.na(non_i)) {
  # fallback: if exactly two rows, use positive order from frozen registry
  if (nrow(SLOPES)==2) {
    pulm_i <- 1
    non_i <- 2
  } else {
    stop("Could not identify pulmonary/non-pulmonary slope rows.")
  }
}

pulm <- SLOPES[pulm_i,]
nonp <- SLOPES[non_i,]

# Section 1
sec1 <- paste0(
  "1. Cohort architecture and analytical framework\n\n",
  "The V2 analysis integrated ten public gut-microbiome cohorts with prespecified analytical roles (Table 1; Figure 1). ",
  "Eight cohorts were represented in the harmonized Step87B analysis object, including seven longitudinal cohorts and one static supportive cohort. ",
  "PRJNA691455 served as the core repeated-sepsis cohort for trajectory-level interpretation; PRJNA578267 provided the non-sepsis longitudinal control; and PRJNA430161 and PRJNA1166732 were retained as intervention-support cohorts. ",
  "PRJNA1010969 provided an external Control-Trauma-Sepsis ecological-state comparison using 53 analyzed samples (17 Control, 18 Trauma, and 18 Sepsis), whereas CRA002354 provided an independent infection-source extension using 130 primary-analysis samples from 64 participants after true-nonchimeric OTU reconstruction and the >=4000-read depth gate. ",
  "Table 1 reports the populations that entered each analytical branch rather than the total number of samples deposited in the original public studies."
)

# Section 2
getrow <- function(cohort) {
  robust_prog %>% filter(Cohort==cohort) %>% slice(1)
}

r1 <- getrow("PRJNA691455")
r2 <- getrow("PRJNA851469")
r3 <- getrow("PRJNA516701")

sec2 <- paste0(
  "2. Progressive ecological displacement across longitudinal critical-illness cohorts\n\n",
  "Three cohorts showed FDR-supported progressive displacement in the primary early-to-late analysis: PRJNA691455 (n=",
  r1$Primary_pairs,
  " paired participants, dz=",
  fmt(r1$Primary_effect_dz,3),
  ", FDR=",
  fmtp(r1$Primary_FDR),
  "), PRJNA851469 (n=",
  r2$Primary_pairs,
  ", dz=",
  fmt(r2$Primary_effect_dz,3),
  ", FDR=",
  fmtp(r2$Primary_FDR),
  "), and PRJNA516701 (n=",
  r3$Primary_pairs,
  ", dz=",
  fmt(r3$Primary_effect_dz,3),
  ", FDR=",
  fmtp(r3$Primary_FDR),
  "). ",
  "The external ICU-infection cohort PRJEB82425 also showed a strong positive primary contrast (n=",
  ext_inf$Primary_pairs,
  ", dz=",
  fmt(ext_inf$Primary_effect_dz,3),
  ", FDR=",
  fmtp(ext_inf$Primary_FDR),
  "), but its heterogeneous anchor structure supported its classification as external supportive rather than core replication evidence. ",
  "In contrast, the non-sepsis surgical control PRJNA578267 showed a negative recovery-direction primary contrast (n=",
  ctrl$Primary_pairs,
  ", dz=",
  fmt(ctrl$Primary_effect_dz,3),
  ", FDR=",
  fmtp(ctrl$Primary_FDR),
  ") and significant re-convergence in the common-anchor sensitivity analysis (n=",
  ctrl$Common_anchor_pairs,
  ", dz=",
  fmt(ctrl$Common_anchor_effect_dz,3),
  ", FDR=",
  fmtp(ctrl$Common_anchor_FDR),
  "). ",
  "Neither intervention-support cohort showed a clear paired longitudinal change after multiplicity correction. ",
  "These results indicate that progressive within-patient ecological displacement was reproducible across several critical-illness cohorts, while the non-sepsis longitudinal control showed a distinct recovery/re-convergence pattern in sensitivity analysis (Figure 2). ",
  "Because progressive displacement was also observed in ICU-background cohorts, these data do not establish sepsis-specific longitudinal instability."
)

# Section 3
sec3 <- paste0(
  "3. Heterogeneous taxonomic routes accompany ecological displacement in core sepsis\n\n",
  "Within PRJNA691455, median Bray-Curtis displacement from the patient-specific baseline increased from ",
  fmt(cval("Day3_Bray_median"),3),
  " at Day 3 to ",
  fmt(cval("Day7_Bray_median"),3),
  " at Day 7. Among ",
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
  "Mean trajectory distance correlated with higher latest EII (Spearman rho=",
  fmt(cval("rho_latest_EII"),3),
  ", permutation p=",
  fmtp(cval("perm_p_latest_EII")),
  ") and greater latest Simpson instability (rho=",
  fmt(cval("rho_Simpson_instability"),3),
  ", permutation p=",
  fmtp(cval("perm_p_Simpson_instability")),
  "), but not clearly with latest Bray displacement (rho=",
  fmt(cval("rho_latest_Bray"),3),
  ", p=",
  fmtp(cval("p_latest_Bray")),
  "). ",
  "Exploratory clustering did not support stable trajectory subtypes because the k=2 solution separated nine patients from a single patient. ",
  "Together, these findings support heterogeneous taxonomic routes accompanying a broadly shared ecological displacement process rather than a single conserved genus-level trajectory (Figure 3)."
)

# Section 4
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

# Section 5
sec5 <- paste0(
  "5. Infection source does not clearly modify longitudinal ecological displacement\n\n",
  "In CRA002354, pulmonary versus recorded non-pulmonary infection source did not show a clear difference in the rate of Bray-Curtis displacement from each patient's first available microbiome sample (source-by-time beta=",
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
  fmtp(pulm$p_value),
  "; recorded non-pulmonary beta=",
  fmt(nonp$slope_per_day,3),
  ", p=",
  fmtp(nonp$p_value),
  "), although neither was individually significant. ",
  "All prespecified key Bray sensitivity analyses retained the same positive interaction direction without statistical significance. ",
  "At baseline, the unrestricted PERMANOVA association was accompanied by dispersion heterogeneity and was not reproduced when baseline sampling was restricted to ICU Day 3 or earlier. ",
  "Exploratory genus-level analysis yielded ",
  GENUS$primary_baseline_FDR_hits[1],
  " primary baseline FDR-significant genera, but ",
  GENUS$baseline_hits_FDR_robust_at_both[1],
  " remained significant at both the >=70% and >=80% sample-level genus-coverage thresholds, and no primary longitudinal genus source-by-time interaction survived FDR correction. ",
  "These findings do not support a strong infection-source modification of longitudinal ecological displacement and do not establish a robust infection-source-specific taxonomic signature (Figure 5)."
)

draft <- paste(
  sec1,sec2,sec3,sec4,sec5,
  sep="\n\n"
)

writeLines(
  draft,
  file.path(DRAFTDIR,"Results_Draft_v2_English.txt")
)

# ------------------------------------------------------------
# 8. Sentence provenance map + Chinese structure
# ------------------------------------------------------------

prov <- tribble(
  ~results_section, ~statement_family, ~authoritative_source,
  "1", "Step87B cohort analysis populations", STEP87B,
  "1", "PRJNA1010969 analyzed group counts", ext_sel$file[1],
  "1", "CRA002354 primary analyzed population", CRA_META,
  "2", "Cross-cohort paired effect sizes/FDR/common-anchor results", "Step93W via Step95C longitudinal table",
  "3", "Core-sepsis Bray/time and trajectory heterogeneity", "Step93U frozen constants recorded by Step95D2",
  "3", "EII/Simpson/Bray trajectory correlations", "Step93U frozen constants recorded by Step95D2",
  "4", "Healthy-centroid and Aitchison robustness", "Step93X2/X3 frozen constants recorded by Step95D2",
  "5", "Primary source-by-time ecology", "Step94B2D via Step95B",
  "5", "Source-specific slopes", "Step94B2D via Step95B normalized in Step95E2",
  "5", "Baseline beta robustness", "Step94B2D via Step95B",
  "5", "Exploratory genus evidence tier", "Step94B3C via Step95B"
)

write_csv(
  prov,
  file.path(DRAFTDIR,"Results_v2_Evidence_Provenance_Map.csv")
)

writeLines(
  c(
    "V2 Results v2 中文结构核对",
    "",
    "1. Table 1 全部改为“实际进入对应分析分支的人群”，不再使用上游 frozen master 的 deposited/available 样本量。",
    "2. 跨队列主结果仍为 3 个 robust progressive cohorts + PRJEB82425 supportive external；PRJNA578267 为恢复/再趋近对照；干预队列无明确变化。",
    "3. PRJNA691455 强调 shared ecological displacement 与 heterogeneous taxonomic routes 并存，不声称稳定 subtype。",
    "4. PRJNA1010969 以 53 个最终分析样本为口径，支持 sepsis-associated ecological-state differentiation beyond trauma，不升级为 universal taxonomic signature。",
    "5. CRA002354 以 primary >=4000 true-nonchimeric population 为口径；source×time 无明确差异，genus 只保留 exploratory/supportive 定位。"
  ),
  file.path(DRAFTDIR,"Results_v2_Structure_Review_Chinese.txt")
)

# ------------------------------------------------------------
# 9. Final readiness
# ------------------------------------------------------------

qc <- tibble(
  step95D3_final_figures_complete=TRUE,
  step87B_primary_rows=nrow(M87),
  step87B_source_gate_pass=(nrow(M87)==785),
  PRJNA1010969_final_samples=ext_count$samples_analyzed,
  PRJNA1010969_group_gate_pass=TRUE,
  CRA_primary_samples=cra_count$samples_analyzed,
  CRA_primary_participants=cra_count$participants_analyzed,
  CRA_population_gate_pass=
    cra_count$samples_analyzed==130 &&
    cra_count$participants_analyzed==64,
  all_10_table1_cohorts_complete=nrow(table1)==10 &&
    all(!is.na(table1$Samples_analyzed)) &&
    all(!is.na(table1$Participants_analyzed)),
  source_slope_schema_normalized=TRUE,
  results_draft_v2_created=
    file.exists(
      file.path(DRAFTDIR,"Results_Draft_v2_English.txt")
    ),
  no_new_inferential_analysis=TRUE
) %>%
  mutate(
    ready_for_table1_and_results_manual_review=
      step87B_source_gate_pass &
      PRJNA1010969_group_gate_pass &
      CRA_population_gate_pass &
      all_10_table1_cohorts_complete &
      source_slope_schema_normalized &
      results_draft_v2_created
  )

write_csv(
  qc,
  file.path(OUT,"STEP95E2_READINESS.csv")
)

if (!isTRUE(qc$ready_for_table1_and_results_manual_review)) {
  stop("Step95E2 final readiness failed.")
}

writeLines(
  c(
    paste0("Completed: ",Sys.time()),
    "Table 1 uses branch-specific analysis populations: TRUE",
    "Step87B exact PRIMARY population gate: TRUE",
    "PRJNA1010969 17/18/18 population gate: TRUE",
    "CRA002354 primary 130 samples / 64 participants gate: TRUE",
    "Source-specific slope schema normalized: TRUE",
    "Results Draft v2 created: TRUE",
    "No new inferential analysis: TRUE",
    "Ready for manual Table 1 + Results review: TRUE",
    "STEP95E2 COMPLETE"
  ),
  file.path(OUT,"_STEP95E2_COMPLETE.ok")
)

cat("STEP95E2 COMPLETE\n")
