# ============================================================
# Sepsis V2 - Step95E
# CONVENTIONAL COHORT TABLE 1 + RESULTS DRAFT V1
#
# Inputs:
# - Step95B manuscript source pack
# - Step95D2 frozen constants/provenance
# - Step95D3 final figures
# - analysis metadata discovered from the project itself
#
# Purpose:
# 1) build a conventional study/cohort characteristics Table 1;
# 2) source-anchor sample / participant counts;
# 3) generate a first manuscript-style English Results draft;
# 4) generate an evidence-to-sentence provenance map.
#
# NO new inferential analysis.
# ============================================================

options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(stringr)
  library(tidyr)
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
  "V2_35E_TABLE1_AND_RESULTS_DRAFT"
)

dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

TABDIR <- file.path(OUT, "01_TABLES")
DRAFTDIR <- file.path(OUT, "02_RESULTS_DRAFT")
AUDDIR <- file.path(OUT, "03_PROVENANCE_AND_AUDIT")

for (d in c(TABDIR,DRAFTDIR,AUDDIR)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

if (!file.exists(file.path(SRC, "_STEP95B_COMPLETE.ok"))) {
  stop("Step95B source pack is incomplete.")
}

if (!file.exists(file.path(D2, "_STEP95D2_COMPLETE.ok"))) {
  stop("Step95D2 is incomplete.")
}

if (!file.exists(file.path(D3, "_STEP95D3_COMPLETE.ok"))) {
  stop("Step95D3 final figures are incomplete.")
}

# ------------------------------------------------------------
# 1. Frozen analysis sources
# ------------------------------------------------------------

W <- read_csv(
  file.path(
    RESULTS,
    "V2_35C_MANUSCRIPT_ARCHITECTURE_FIGURE1_AND_MAIN_TABLES",
    "Main_Table_2_Cross_Cohort_Longitudinal_Results.csv"
  ),
  show_col_types = FALSE
)

EVID <- read_csv(
  file.path(
    SRC,
    "05_PROVENANCE",
    "STEP95A2_06_FINAL_GLOBAL_V2_EVIDENCE_REGISTRY.csv"
  ),
  show_col_types = FALSE
)

EXT_TAX <- read_csv(
  file.path(
    SRC,
    "01_MAIN_TABLE_SOURCES",
    "MAIN_T4_EXTERNAL_VALIDATION.csv"
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

# ------------------------------------------------------------
# 2. Target cohort registry
# ------------------------------------------------------------

targets <- tribble(
  ~cohort, ~clinical_context, ~analysis_role, ~main_manuscript_role,
  "PRJNA691455", "Sepsis / critical illness", "CORE_SEPSIS_LONGITUDINAL", "PRIMARY",
  "PRJEB82425", "ICU infection", "ICU_INFECTION_LONGITUDINAL_EXTERNAL", "SUPPORTIVE",
  "PRJNA516701", "ICU background", "ICU_BACKGROUND_LONGITUDINAL", "PRIMARY_CONTEXT",
  "PRJNA851469", "ICU background", "ICU_BACKGROUND_LONGITUDINAL", "PRIMARY_CONTEXT",
  "PRJNA578267", "Non-sepsis surgical longitudinal control", "NONSEPSIS_LONGITUDINAL_CONTROL", "CONTROL",
  "PRJNA430161", "Intervention longitudinal cohort", "INTERVENTION_LONGITUDINAL_SUPPORT", "SUPPORTIVE",
  "PRJNA1166732", "Intervention longitudinal cohort", "INTERVENTION_LONGITUDINAL_SUPPORT", "SUPPORTIVE",
  "PRJNA978257", "Static supportive cohort", "STATIC_SUPPORT", "SUPPLEMENTARY",
  "PRJNA1010969", "Control / severe trauma / sepsis", "EXTERNAL_ECOLOGICAL_STATE_VALIDATION", "SUPPORTIVE_EXTERNAL",
  "CRA002354", "Sepsis stratified by infection source", "INFECTION_SOURCE_LONGITUDINAL", "SUPPORTIVE_EXTENSION"
)

write_csv(
  targets,
  file.path(AUDDIR, "01_TARGET_COHORT_REGISTRY.csv")
)

# ------------------------------------------------------------
# 3. Metadata-source discovery
# ------------------------------------------------------------

search_roots <- c(
  file.path(DATA, "_V2_ANALYSIS_READY"),
  file.path(DATA, "CRA002354"),
  RESULTS
)

search_roots <- search_roots[dir.exists(search_roots)]

all_files <- unlist(
  lapply(
    search_roots,
    function(d) {
      list.files(
        d,
        recursive = TRUE,
        full.names = TRUE
      )
    }
  ),
  use.names = FALSE
)

candidate_files <- all_files[
  str_detect(
    tolower(all_files),
    "\\.(csv|tsv|txt)$"
  ) &
  str_detect(
    tolower(basename(all_files)),
    "metadata|master|analysis|sample|87a|87b|patient"
  )
]

# Exclude manuscript-production copies and very large files.
candidate_files <- candidate_files[
  !str_detect(
    candidate_files,
    "V2_35[A-Z0-9_]*"
  )
]

candidate_files <- candidate_files[
  file.info(candidate_files)$size <= 100 * 1024^2
]

read_candidate <- function(path, n_max = Inf) {

  ext <- tolower(tools::file_ext(path))

  tryCatch(
    {
      if (ext == "tsv") {
        read_tsv(
          path,
          show_col_types = FALSE,
          n_max = n_max,
          progress = FALSE
        )
      } else {
        read_csv(
          path,
          show_col_types = FALSE,
          n_max = n_max,
          progress = FALSE
        )
      }
    },
    error = function(e) NULL
  )
}

first_col <- function(nms, choices) {
  hit <- choices[choices %in% nms]
  if (length(hit) == 0) return(NA_character_)
  hit[1]
}

project_choices <- c(
  "project",
  "project_raw",
  "Project",
  "project_id",
  "project_accession",
  "BioProject",
  "bioproject",
  "study_accession"
)

sample_choices <- c(
  "Run_ID",
  "run_id",
  "run",
  "sample_id",
  "Sample_ID",
  "sample",
  "accession"
)

patient_choices <- c(
  "patient_id",
  "Patient_ID",
  "subject_id",
  "Subject_ID",
  "participant_id",
  "host_subject_id",
  "individual_id"
)

time_choices <- c(
  "time_day",
  "day",
  "Day",
  "days",
  "sampling_day",
  "days_since_baseline",
  "days_since_personal_baseline"
)

candidate_header <- bind_rows(
  lapply(
    candidate_files,
    function(f) {

      z <- read_candidate(f, n_max = 5)

      if (is.null(z)) {
        return(
          tibble(
            file = f,
            readable = FALSE,
            project_col = NA_character_,
            sample_col = NA_character_,
            patient_col = NA_character_,
            time_col = NA_character_
          )
        )
      }

      nms <- names(z)

      tibble(
        file = f,
        readable = TRUE,
        project_col = first_col(nms, project_choices),
        sample_col = first_col(nms, sample_choices),
        patient_col = first_col(nms, patient_choices),
        time_col = first_col(nms, time_choices)
      )
    }
  )
)

candidate_header <- candidate_header %>%
  mutate(
    path_contains_target =
      vapply(
        file,
        function(x) {
          any(
            str_detect(
              toupper(x),
              fixed(targets$cohort)
            )
          )
        },
        logical(1)
      ),
    structural_score =
      5 * !is.na(sample_col) +
      4 * !is.na(patient_col) +
      3 * !is.na(project_col) +
      1 * !is.na(time_col) +
      2 * path_contains_target +
      2 * str_detect(tolower(basename(file)), "master|87a|87b") +
      1 * str_detect(tolower(basename(file)), "analysis_metadata")
  ) %>%
  arrange(desc(structural_score))

write_csv(
  candidate_header,
  file.path(AUDDIR, "02_METADATA_CANDIDATE_HEADER_AUDIT.csv")
)

# ------------------------------------------------------------
# 4. Build cohort counts from source-anchored metadata
# ------------------------------------------------------------

extract_for_target <- function(target) {

  rows <- list()

  # Try candidates in priority order.
  for (i in seq_len(nrow(candidate_header))) {

    info <- candidate_header[i,]

    if (
      !isTRUE(info$readable) ||
      is.na(info$sample_col)
    ) next

    # Candidate must either have a project column or target in its path.
    if (
      is.na(info$project_col) &&
      !str_detect(
        toupper(info$file),
        fixed(target)
      )
    ) next

    z <- read_candidate(info$file)

    if (is.null(z) || nrow(z) == 0) next

    # Project filtering.
    if (!is.na(info$project_col)) {

      pv <- toupper(
        as.character(
          z[[info$project_col]]
        )
      )

      hit <- str_detect(
        pv,
        fixed(toupper(target))
      )

      if (!any(hit, na.rm = TRUE)) {
        # Some rows may use aliases; only CRA currently needs explicit alias support.
        if (target == "CRA002354") {
          hit <- str_detect(
            pv,
            "CRA002354|PRJCA002241"
          )
        }
      }

      if (!any(hit, na.rm = TRUE)) next

      z <- z[which(hit),,drop=FALSE]
    }

    if (nrow(z) == 0) next

    sid <- as.character(
      z[[info$sample_col]]
    )

    valid_sid <- !is.na(sid) & sid != ""

    z <- z[valid_sid,,drop=FALSE]
    sid <- sid[valid_sid]

    if (length(sid) == 0) next

    # Deduplicate samples.
    keep <- !duplicated(sid)
    z <- z[keep,,drop=FALSE]
    sid <- sid[keep]

    n_samples <- length(unique(sid))

    n_patients <- NA_integer_
    n_repeated <- NA_integer_

    if (!is.na(info$patient_col)) {

      pid <- as.character(
        z[[info$patient_col]]
      )

      pid[pid == ""] <- NA_character_

      n_patients <- n_distinct(
        pid[!is.na(pid)]
      )

      if (n_patients > 0) {

        patient_counts <- tibble(
          patient_id = pid,
          sample_id = sid
        ) %>%
          filter(
            !is.na(patient_id)
          ) %>%
          distinct(
            patient_id,
            sample_id
          ) %>%
          count(
            patient_id,
            name = "n"
          )

        n_repeated <- sum(
          patient_counts$n >= 2
        )
      }
    }

    day_min <- NA_real_
    day_max <- NA_real_

    if (!is.na(info$time_col)) {

      td <- suppressWarnings(
        as.numeric(
          z[[info$time_col]]
        )
      )

      if (any(is.finite(td))) {
        day_min <- min(td[is.finite(td)])
        day_max <- max(td[is.finite(td)])
      }
    }

    rows[[length(rows)+1]] <- tibble(
      cohort = target,
      metadata_source = info$file,
      sample_col = info$sample_col,
      patient_col = info$patient_col,
      time_col = info$time_col,
      n_samples = n_samples,
      n_participants = n_patients,
      n_repeated_participants = n_repeated,
      observed_day_min = day_min,
      observed_day_max = day_max,
      source_priority_score = info$structural_score
    )
  }

  if (length(rows) == 0) {
    return(
      tibble(
        cohort = target,
        metadata_source = NA_character_,
        sample_col = NA_character_,
        patient_col = NA_character_,
        time_col = NA_character_,
        n_samples = NA_integer_,
        n_participants = NA_integer_,
        n_repeated_participants = NA_integer_,
        observed_day_min = NA_real_,
        observed_day_max = NA_real_,
        source_priority_score = NA_real_
      )
    )
  }

  rr <- bind_rows(rows) %>%
    mutate(
      has_patient_count = !is.na(n_participants) & n_participants > 0,
      count_plausibility =
        !is.na(n_samples) &
        n_samples > 0 &
        (
          is.na(n_participants) |
          n_samples >= n_participants
        ),
      rank_score =
        source_priority_score +
        4 * has_patient_count +
        2 * count_plausibility
    ) %>%
    arrange(
      desc(rank_score),
      desc(n_samples)
    )

  # Keep selected source.
  rr[1,]
}

metadata_counts <- bind_rows(
  lapply(
    targets$cohort,
    extract_for_target
  )
)

write_csv(
  metadata_counts,
  file.path(AUDDIR, "03_SELECTED_COHORT_METADATA_SOURCES.csv")
)

# Hard gate:
# - every cohort needs a source-anchored sample count;
# - longitudinal/core/source cohorts need participant counts;
# - static/external cross-sectional cohorts may report NR if participant
#   identifier is not present in the selected metadata.
main_longitudinal <- c(
  "PRJNA691455",
  "PRJEB82425",
  "PRJNA516701",
  "PRJNA851469",
  "PRJNA578267",
  "PRJNA430161",
  "PRJNA1166732",
  "CRA002354"
)

metadata_gate <- metadata_counts %>%
  mutate(
    sample_count_pass =
      !is.na(n_samples) &
      n_samples > 0,
    participant_count_required =
      cohort %in% main_longitudinal,
    participant_count_pass =
      ifelse(
        participant_count_required,
        !is.na(n_participants) &
          n_participants > 0,
        TRUE
      ),
    gate_pass =
      sample_count_pass &
      participant_count_pass
  )

write_csv(
  metadata_gate,
  file.path(AUDDIR, "04_TABLE1_METADATA_GATE.csv")
)

if (!all(metadata_gate$gate_pass)) {

  writeLines(
    c(
      "STEP95E STOPPED BEFORE TABLE 1 / RESULTS DRAFT",
      "",
      "At least one cohort lacks a source-anchored sample count or a required longitudinal participant count.",
      "No missing count was filled from memory.",
      "",
      "Inspect:",
      "02_METADATA_CANDIDATE_HEADER_AUDIT.csv",
      "03_SELECTED_COHORT_METADATA_SOURCES.csv",
      "04_TABLE1_METADATA_GATE.csv"
    ),
    file.path(OUT, "_STEP95E_METADATA_GATE_FAILED.txt")
  )

  stop("Table 1 metadata source gate failed.")
}

# ------------------------------------------------------------
# 5. Build conventional manuscript Table 1
# ------------------------------------------------------------

long_evidence <- W %>%
  transmute(
    cohort = Cohort,
    primary_pairs = Primary_pairs,
    common_anchor_pairs = Common_anchor_pairs,
    longitudinal_evidence_class = Evidence_class
  )

table1 <- targets %>%
  left_join(
    metadata_counts %>%
      select(
        cohort,
        n_samples,
        n_participants,
        n_repeated_participants,
        observed_day_min,
        observed_day_max,
        metadata_source
      ),
    by = "cohort"
  ) %>%
  left_join(
    long_evidence,
    by = "cohort"
  ) %>%
  mutate(
    participant_display =
      ifelse(
        is.na(n_participants),
        "NR",
        as.character(n_participants)
      ),
    repeated_display =
      ifelse(
        is.na(n_repeated_participants),
        "NR",
        as.character(n_repeated_participants)
      ),
    time_window =
      case_when(
        !is.na(observed_day_min) &
          !is.na(observed_day_max) ~
          paste0(
            observed_day_min,
            "-",
            observed_day_max,
            " d"
          ),
        TRUE ~
          "NR"
      )
  ) %>%
  transmute(
    Cohort = cohort,
    Clinical_context = clinical_context,
    Analysis_role = analysis_role,
    Samples = n_samples,
    Participants = participant_display,
    Repeated_participants = repeated_display,
    Observed_time_window = time_window,
    Primary_longitudinal_pairs =
      ifelse(
        is.na(primary_pairs),
        NA,
        primary_pairs
      ),
    Longitudinal_evidence_class =
      longitudinal_evidence_class,
    Manuscript_role =
      main_manuscript_role
  )

write_csv(
  table1,
  file.path(
    TABDIR,
    "Main_Table_1_Cohort_and_Study_Characteristics.csv"
  )
)

# Full provenance version.
table1_prov <- targets %>%
  left_join(
    metadata_counts,
    by = "cohort"
  ) %>%
  left_join(
    long_evidence,
    by = "cohort"
  )

write_csv(
  table1_prov,
  file.path(
    AUDDIR,
    "05_TABLE1_FULL_PROVENANCE.csv"
  )
)

# ------------------------------------------------------------
# 6. Keep the frozen longitudinal result table as result table
# ------------------------------------------------------------

write_csv(
  W,
  file.path(
    TABDIR,
    "Main_or_Supp_Table_Cross_Cohort_Longitudinal_Results.csv"
  )
)

# ------------------------------------------------------------
# 7. Results evidence registry
# ------------------------------------------------------------

primary_cra <- ECO %>%
  filter(
    domain == "PRIMARY_LONGITUDINAL_BRAY",
    estimand == "source_by_personal_time_interaction"
  ) %>%
  slice(1)

source_p <- SLOPES %>%
  transmute(
    source_group,
    slope_per_day,
    ci95_low,
    ci95_high,
    p_LRT
  )

robust_prog <- W %>%
  filter(
    Cohort %in%
      c(
        "PRJNA691455",
        "PRJNA851469",
        "PRJNA516701"
      )
  )

external_infect <- W %>%
  filter(
    Cohort == "PRJEB82425"
  )

control <- W %>%
  filter(
    Cohort == "PRJNA578267"
  )

intervention <- W %>%
  filter(
    Cohort %in%
      c(
        "PRJNA1166732",
        "PRJNA430161"
      )
  )

results_evidence <- tribble(
  ~section, ~claim_key, ~source,
  "2", "Three robust progressive longitudinal cohorts", "Step93W / Step95C longitudinal table",
  "2", "External ICU-infection cohort supportive but anchor-heterogeneous", "Step93W",
  "2", "Non-sepsis surgical control recovery/re-convergence sensitivity", "Step93W",
  "2", "Intervention cohorts no clear longitudinal pattern", "Step93W",
  "3", "Core sepsis Day3-Day7 ecological displacement", "Step93U frozen constants",
  "3", "Heterogeneous taxonomic routes", "Step93U frozen constants",
  "3", "Trajectory heterogeneity associated with EII and Simpson instability", "Step93U frozen constants",
  "4", "Sepsis farther from healthy centroid than severe trauma", "Step93X2 frozen constants",
  "4", "CLR/Aitchison separation robust without dispersion signal", "Step93X3 frozen constants",
  "4", "No reproducible universal genus signature established", "Step93X4 frozen wording",
  "5", "No clear infection-source modification of Bray displacement", "Step94B2D frozen registry",
  "5", "Positive estimated source-specific slopes but neither individually significant", "Step94B2D frozen slopes",
  "5", "Baseline PERMANOVA not robust and dispersion-confounded", "Step94B2D baseline-beta freeze",
  "5", "Genus results limited to exploratory baseline signals", "Step94B3C freeze"
)

write_csv(
  results_evidence,
  file.path(
    DRAFTDIR,
    "Results_Evidence_to_Sentence_Map.csv"
  )
)

# ------------------------------------------------------------
# 8. Manuscript-style Results draft
# ------------------------------------------------------------

fmt <- function(x, digits=3) {
  ifelse(
    is.na(x),
    "NA",
    formatC(
      x,
      format = "f",
      digits = digits
    )
  )
}

fmtp <- function(x) {
  ifelse(
    is.na(x),
    "NA",
    formatC(
      x,
      format = "g",
      digits = 3
    )
  )
}

rp_sentence <- paste0(
  "Three longitudinal cohorts showed FDR-supported progressive displacement in the primary early-to-late analysis: PRJNA691455 (n=",
  robust_prog$Primary_pairs[robust_prog$Cohort=="PRJNA691455"],
  " paired participants, dz=",
  fmt(robust_prog$Primary_effect_dz[robust_prog$Cohort=="PRJNA691455"],3),
  ", FDR=",
  fmtp(robust_prog$Primary_FDR[robust_prog$Cohort=="PRJNA691455"]),
  "), PRJNA851469 (n=",
  robust_prog$Primary_pairs[robust_prog$Cohort=="PRJNA851469"],
  ", dz=",
  fmt(robust_prog$Primary_effect_dz[robust_prog$Cohort=="PRJNA851469"],3),
  ", FDR=",
  fmtp(robust_prog$Primary_FDR[robust_prog$Cohort=="PRJNA851469"]),
  "), and PRJNA516701 (n=",
  robust_prog$Primary_pairs[robust_prog$Cohort=="PRJNA516701"],
  ", dz=",
  fmt(robust_prog$Primary_effect_dz[robust_prog$Cohort=="PRJNA516701"],3),
  ", FDR=",
  fmtp(robust_prog$Primary_FDR[robust_prog$Cohort=="PRJNA516701"]),
  ")."
)

ext_sentence <- paste0(
  "The external ICU-infection cohort PRJEB82425 showed a larger positive primary effect (n=",
  external_infect$Primary_pairs[1],
  ", dz=",
  fmt(external_infect$Primary_effect_dz[1],3),
  ", FDR=",
  fmtp(external_infect$Primary_FDR[1]),
  "), but heterogeneous anchor structure supported its use as external supportive rather than core replication evidence."
)

control_sentence <- paste0(
  "By contrast, the non-sepsis surgical longitudinal control PRJNA578267 showed a negative recovery-direction primary contrast (n=",
  control$Primary_pairs[1],
  ", dz=",
  fmt(control$Primary_effect_dz[1],3),
  ", FDR=",
  fmtp(control$Primary_FDR[1]),
  ") and significant re-convergence in the common-anchor sensitivity analysis (n=",
  control$Common_anchor_pairs[1],
  ", dz=",
  fmt(control$Common_anchor_effect_dz[1],3),
  ", FDR=",
  fmtp(control$Common_anchor_FDR[1]),
  ")."
)

intervention_sentence <- paste0(
  "Neither intervention-support cohort showed a clear primary longitudinal displacement pattern after multiplicity correction."
)

section1 <- paste0(
  "1. Cohort architecture and analytical framework\n\n",
  "The V2 analysis integrated ten public microbiome cohorts with distinct prespecified analytical roles (Table 1; Figure 1). Seven cohorts contributed to the longitudinal critical-illness framework, with PRJNA691455 designated as the core repeated-sepsis cohort for trajectory-level interpretation. PRJNA578267 served as the non-sepsis longitudinal control, whereas PRJNA430161 and PRJNA1166732 were retained as intervention-support cohorts. PRJNA1010969 provided an external Control-Trauma-Sepsis ecological-state comparison, and CRA002354 provided an independent infection-source longitudinal extension. PRJNA978257 was retained as static supportive evidence rather than a longitudinal replication cohort. Cohort-level sample and participant counts in Table 1 were derived directly from source-anchored project metadata rather than reconstructed from inferential result files."
)

section2 <- paste0(
  "2. Progressive ecological displacement across longitudinal critical-illness cohorts\n\n",
  rp_sentence, " ",
  ext_sentence, " ",
  control_sentence, " ",
  intervention_sentence, " ",
  "Together, these results indicate that progressive within-patient ecological displacement was reproducible across several critical-illness cohorts, while the longitudinal non-sepsis control showed a distinct recovery/re-convergence pattern in sensitivity analysis (Figure 2). Because progressive displacement was also observed in ICU-background cohorts, these data do not establish sepsis-specific longitudinal instability."
)

section3 <- paste0(
  "3. Heterogeneous taxonomic routes accompany ecological displacement in core sepsis\n\n",
  "Within the core repeated-sepsis cohort PRJNA691455, median Bray-Curtis displacement from the patient-specific baseline increased from ",
  fmt(cval("Day3_Bray_median"),3),
  " at Day 3 to ",
  fmt(cval("Day7_Bray_median"),3),
  " at Day 7. Among ",
  as.integer(cval("Day3_to_Day7_paired_n")),
  " patients with paired Day 3 and Day 7 observations, the mean increase was ",
  fmt(cval("Day3_to_Day7_mean_change"),3),
  ". Despite this shared ecological movement, patient-specific taxonomic trajectories were heterogeneous: the mean pairwise signed-trajectory Euclidean distance was ",
  fmt(cval("pairwise_signed_Euclidean_mean"),3),
  " (median ",
  fmt(cval("pairwise_signed_Euclidean_median"),3),
  "), the mean cosine similarity was only ",
  fmt(cval("mean_cosine_similarity"),3),
  ", and ",
  fmt(100*cval("fraction_cosine_le_0"),1),
  "% of patient-pair trajectory comparisons had cosine similarity <=0. Mean trajectory distance was associated with greater latest EII (Spearman rho=",
  fmt(cval("rho_latest_EII"),3),
  ", permutation p=",
  fmtp(cval("perm_p_latest_EII")),
  ") and greater latest Simpson instability (rho=",
  fmt(cval("rho_Simpson_instability"),3),
  ", permutation p=",
  fmtp(cval("perm_p_Simpson_instability")),
  "), but not clearly with latest Bray displacement itself (rho=",
  fmt(cval("rho_latest_Bray"),3),
  ", p=",
  fmtp(cval("p_latest_Bray")),
  "). Exploratory clustering did not support stable trajectory subtypes because the k=2 solution separated nine patients from a single patient. These findings therefore support heterogeneous taxonomic routes accompanying a broadly shared ecological displacement process rather than a single conserved genus-level trajectory (Figure 3)."
)

section4 <- paste0(
  "4. External ecological-state differentiation beyond severe trauma\n\n",
  "In PRJNA1010969, sepsis samples were farther from the healthy-control centroid than severe-trauma samples (",
  fmt(cval("Sepsis_healthy_centroid_distance"),3),
  " vs ",
  fmt(cval("Trauma_healthy_centroid_distance"),3),
  "; difference +",
  fmt(cval("Sepsis_minus_Trauma_centroid_difference"),3),
  "; FDR=",
  fmtp(cval("centroid_difference_FDR")),
  "). CLR/Aitchison sensitivity analyses showed stable Sepsis-versus-Trauma separation across pseudocount choices, with PERMANOVA R2 ranging from ",
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
  "). The external cohort therefore supports sepsis-associated ecological-state differentiation beyond severe trauma. However, taxonomic alignment analyses were exploratory and did not establish a reproducible sepsis-specific genus signature (Figure 4)."
)

pulm <- source_p %>%
  filter(
    str_detect(
      toupper(source_group),
      "PULMONARY"
    ) &
    !str_detect(
      toupper(source_group),
      "NON"
    )
  ) %>%
  slice(1)

nonpulm <- source_p %>%
  filter(
    str_detect(
      toupper(source_group),
      "NON"
    )
  ) %>%
  slice(1)

section5 <- paste0(
  "5. Infection source does not clearly modify longitudinal ecological displacement\n\n",
  "In CRA002354, pulmonary versus recorded non-pulmonary infection source did not show a clear difference in the rate of Bray-Curtis displacement from each patient's first available microbiome sample (source-by-time beta=",
  fmt(primary_cra$effect,3),
  ", 95% CI ",
  fmt(primary_cra$ci_low,3),
  " to ",
  fmt(primary_cra$ci_high,3),
  "; p=",
  fmtp(primary_cra$p_value),
  "). Estimated source-specific slopes were positive in both groups (pulmonary beta=",
  fmt(pulm$slope_per_day,3),
  ", p=",
  fmtp(pulm$p_LRT),
  "; recorded non-pulmonary beta=",
  fmt(nonpulm$slope_per_day,3),
  ", p=",
  fmtp(nonpulm$p_LRT),
  "), although neither slope was individually significant. All prespecified key Bray sensitivity analyses retained the same positive interaction direction without statistical significance. At baseline, the unrestricted PERMANOVA association was accompanied by significant dispersion heterogeneity and was not reproduced when baseline sampling was restricted to ICU Day 3 or earlier. Exploratory genus-level analysis identified four baseline FDR-significant genera in the primary analysis, but none remained significant at both >=70% and >=80% sample-level genus-coverage thresholds, and no primary longitudinal genus source-by-time interaction survived FDR correction. These findings do not support a strong infection-source modification of longitudinal ecological displacement and do not establish a robust infection-source-specific taxonomic signature (Figure 5)."
)

draft <- paste(
  section1,
  section2,
  section3,
  section4,
  section5,
  sep = "\n\n"
)

writeLines(
  draft,
  file.path(
    DRAFTDIR,
    "Results_Draft_v1_English.txt"
  )
)

# Chinese structural outline for manual review; not a translation of every sentence.
outline_cn <- c(
  "V2 Results 结构核对版",
  "",
  "1. 队列架构：7个纵向队列 + PRJNA1010969外部Control-Trauma-Sepsis + CRA002354感染来源扩展 + PRJNA978257静态支持。",
  "2. 跨队列纵向结果：3个队列primary FDR支持progressive displacement；PRJEB82425作为外部支持；非脓毒症手术对照呈恢复/再趋近方向；干预队列无明确模式。",
  "3. Core sepsis：生态位移增强，但患者之间taxonomic trajectory高度异质；trajectory heterogeneity与EII、Simpson instability相关；不支持稳定trajectory subtype。",
  "4. 外部验证：sepsis相较severe trauma离healthy centroid更远；CLR/Aitchison结果稳健且无dispersion混杂；不能升级为universal taxonomic signature。",
  "5. Infection source：无明确source×time modification；两组斜率方向均为正但均不显著；baseline beta结果不稳健；genus仅有稳健性有限的exploratory baseline signals。"
)

writeLines(
  outline_cn,
  file.path(
    DRAFTDIR,
    "Results_Structure_Review_Chinese.txt"
  )
)

# ------------------------------------------------------------
# 9. Final readiness
# ------------------------------------------------------------

qc <- tibble(
  step95D3_final_figures_complete = TRUE,
  all_target_cohorts_sample_count_source_anchored =
    all(metadata_gate$sample_count_pass),
  all_required_longitudinal_participant_counts_source_anchored =
    all(
      metadata_gate$participant_count_pass
    ),
  main_table1_created =
    file.exists(
      file.path(
        TABDIR,
        "Main_Table_1_Cohort_and_Study_Characteristics.csv"
      )
    ),
  results_draft_v1_created =
    file.exists(
      file.path(
        DRAFTDIR,
        "Results_Draft_v1_English.txt"
      )
    ),
  no_new_inferential_analysis = TRUE
) %>%
  mutate(
    ready_for_manual_table1_and_results_review =
      all_target_cohorts_sample_count_source_anchored &
      all_required_longitudinal_participant_counts_source_anchored &
      main_table1_created &
      results_draft_v1_created
  )

write_csv(
  qc,
  file.path(
    OUT,
    "STEP95E_READINESS.csv"
  )
)

if (!isTRUE(qc$ready_for_manual_table1_and_results_review)) {
  stop("Step95E final readiness failed.")
}

writeLines(
  c(
    paste0("Completed: ", Sys.time()),
    "Table 1 metadata source gate: TRUE",
    "Conventional cohort Table 1 created: TRUE",
    "Results Draft v1 created: TRUE",
    "No new inferential analysis: TRUE",
    "Ready for manual Table 1 + Results review: TRUE",
    "STEP95E COMPLETE"
  ),
  file.path(
    OUT,
    "_STEP95E_COMPLETE.ok"
  )
)

cat("STEP95E COMPLETE\n")
