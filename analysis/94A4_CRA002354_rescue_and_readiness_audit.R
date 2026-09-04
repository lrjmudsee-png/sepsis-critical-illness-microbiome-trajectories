# ============================================================
# Sepsis V2 - Step94A4
# CRA002354 / PRJCA002241 RESCUE + ANALYSIS-READINESS AUDIT
#
# Background
# ----------
# Step94A3 did NOT truly show that infection-source analysis was impossible.
# The 131 infection_source_standard rows all belong to project_raw CRA002354.
# The prior regex recognized PRJNA/PRJEB only, so CRA002354 was falsely
# classified as unresolved.
#
# This step:
# 1) repairs accession parsing for CRA / PRJCA;
# 2) freezes patient-level infection-source structure;
# 3) audits longitudinal depth (Day 1-9);
# 4) evaluates pulmonary vs non-pulmonary feasibility;
# 5) scans the local project for already-processed CRA002354/PRJCA002241
#    microbiome tables, including tables discoverable only by CRR run overlap;
# 6) decides whether this cohort can be rescued into V2 without reprocessing.
#
# AUDIT ONLY. No source-stratified microbiome modeling yet.
# ============================================================

options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(stringr)
  library(tibble)
  library(purrr)
})

ROOT <- "E:/sepsis_project"
DATA_ROOT <- file.path(ROOT, "data")

MASTER <- file.path(
  DATA_ROOT,
  "_V2_ANALYSIS_READY",
  "00_FREEZE",
  "V2_FINAL_FROZEN_master_metadata_CLINICAL_REPAIRED_20260814_131620.csv"
)

OUT <- file.path(
  ROOT,
  "results",
  "V2_34A4_CRA002354_RESCUE_AND_READINESS_AUDIT"
)

dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

logfile <- file.path(OUT, "_STEP94A4_runtime.txt")
writeLines(paste("START", Sys.time()), logfile)

logmsg <- function(...) {
  z <- paste0(...)
  cat(z, "\n")
  cat(z, "\n", file = logfile, append = TRUE)
}

norm_chr <- function(x) {
  z <- str_squish(as.character(x))
  z[z %in% c("", "NA", "N/A", "NULL", "None", "none", ".", "-")] <- NA_character_
  z
}

extract_project_accession <- function(x) {
  z <- toupper(norm_chr(x))

  # Support NCBI/ENA and CNCB/GSA project accessions.
  out <- str_extract(
    z,
    "PRJ(?:NA|EB|CA)\\d+|CRA\\d+"
  )

  out
}

safe_num <- function(x) {
  suppressWarnings(as.numeric(as.character(x)))
}

harmonize_source <- function(x) {
  z <- tolower(norm_chr(x))

  case_when(
    is.na(z) ~ NA_character_,

    z %in% c(
      "lung",
      "pulmonary",
      "respiratory",
      "pneumonia",
      "vap"
    ) ~ "RESPIRATORY",

    z %in% c(
      "abdominal",
      "abdomen",
      "intestinal",
      "intra_abdominal",
      "intra-abdominal",
      "gastrointestinal",
      "gi"
    ) ~ "ABDOMINAL_GI",

    z %in% c(
      "blood",
      "bloodstream",
      "bacteremia",
      "bacteraemia"
    ) ~ "BLOODSTREAM",

    z %in% c(
      "urinary",
      "urine",
      "uti",
      "urosepsis"
    ) ~ "URINARY",

    z %in% c(
      "surgical_site",
      "surgical site",
      "wound",
      "skin",
      "soft_tissue"
    ) ~ "SKIN_SOFT_TISSUE",

    z %in% c(
      "other",
      "unknown",
      "undetermined"
    ) ~ "OTHER_UNKNOWN",

    TRUE ~ "UNMAPPED"
  )
}

safe_read <- function(f) {
  ext <- tolower(tools::file_ext(f))

  if (ext == "csv") {
    return(
      tryCatch(
        read_csv(
          f,
          show_col_types = FALSE,
          guess_max = 20000,
          name_repair = "unique"
        ),
        error = function(e) NULL
      )
    )
  }

  if (ext %in% c("tsv", "txt")) {
    return(
      tryCatch(
        read_tsv(
          f,
          show_col_types = FALSE,
          guess_max = 20000,
          name_repair = "unique"
        ),
        error = function(e) NULL
      )
    )
  }

  NULL
}

classify_candidate_file <- function(f, x = NULL) {

  z <- tolower(basename(f))

  if (str_detect(z, "genus") && str_detect(z, "abundance|relative|table|profile")) {
    return("GENUS_ABUNDANCE")
  }

  if (str_detect(z, "asv|otu|feature") && str_detect(z, "table|abundance|count")) {
    return("FEATURE_TABLE")
  }

  if (str_detect(z, "taxonomy")) {
    return("TAXONOMY")
  }

  if (str_detect(z, "metadata|manifest|mapping|run.?table|sample.?sheet")) {
    return("METADATA")
  }

  if (str_detect(z, "phyloseq|\\.rds$")) {
    return("R_OBJECT")
  }

  if (str_detect(z, "alpha|bray|distance|ordination|pcoa")) {
    return("DOWNSTREAM_RESULT")
  }

  if (!is.null(x)) {
    numeric_cols <- sum(vapply(x, is.numeric, logical(1)))
    if (numeric_cols >= 10) return("NUMERIC_MICROBIOME_TABLE")
  }

  "OTHER_TABULAR"
}

# ------------------------------------------------------------
# 1. Read frozen master and repair project accession
# ------------------------------------------------------------

if (!file.exists(MASTER)) {
  stop(paste0("Missing frozen master metadata: ", MASTER))
}

m <- read_csv(
  MASTER,
  show_col_types = FALSE,
  guess_max = 50000,
  name_repair = "unique"
)

required_cols <- c(
  "project",
  "patient_id",
  "run_id",
  "sample_id",
  "time_day",
  "infection_source_standard"
)

missing_required <- setdiff(required_cols, names(m))

if (length(missing_required) > 0) {
  stop(
    paste0(
      "Missing required master columns: ",
      paste(missing_required, collapse = "; ")
    )
  )
}

resolved_master <- m %>%
  transmute(
    project_raw = norm_chr(project),
    project_accession = extract_project_accession(project),
    patient_id = norm_chr(patient_id),
    run_id = norm_chr(run_id),
    sample_id = norm_chr(sample_id),
    time_day = safe_num(time_day),
    infection_source_raw =
      norm_chr(infection_source_standard)
  )

source_rows <- resolved_master %>%
  filter(!is.na(infection_source_raw)) %>%
  mutate(
    infection_source_harmonized =
      harmonize_source(infection_source_raw),
    binary_pulmonary_source = case_when(
      infection_source_harmonized == "RESPIRATORY" ~ "PULMONARY",
      !is.na(infection_source_harmonized) ~ "NONPULMONARY_RECORDED",
      TRUE ~ NA_character_
    )
  )

write_csv(
  source_rows,
  file.path(OUT, "01_REPAIRED_master_source_rows.csv")
)

project_counts <- source_rows %>%
  count(
    project_raw,
    project_accession,
    name = "source_rows"
  )

write_csv(
  project_counts,
  file.path(OUT, "02_REPAIRED_source_rows_by_project.csv")
)

# ------------------------------------------------------------
# 2. CRA002354 alias validation
# ------------------------------------------------------------

cra <- source_rows %>%
  filter(project_accession == "CRA002354")

if (nrow(cra) == 0) {
  stop("No CRA002354 source rows found after repaired accession parsing.")
}

alias_table <- tibble(
  GSA_accession = "CRA002354",
  BioProject_alias = "PRJCA002241",
  cohort_description =
    "Chinese ICU sepsis/septic-shock gut microbiome cohort",
  expected_samples_publication = 131,
  expected_patients_publication = 64,
  observed_source_rows = nrow(cra),
  observed_patients =
    n_distinct(cra$patient_id)
)

write_csv(
  alias_table,
  file.path(OUT, "03_CRA002354_ACCESSION_ALIAS.csv")
)

# ------------------------------------------------------------
# 3. Patient-level source consistency and category counts
# ------------------------------------------------------------

patient_source <- cra %>%
  filter(!is.na(patient_id)) %>%
  group_by(patient_id) %>%
  summarise(
    n_samples = n(),
    min_day = min(time_day, na.rm = TRUE),
    max_day = max(time_day, na.rm = TRUE),
    n_distinct_raw_source =
      n_distinct(infection_source_raw),
    n_distinct_harmonized_source =
      n_distinct(infection_source_harmonized),
    raw_sources =
      paste(
        sort(unique(infection_source_raw)),
        collapse = ";"
      ),
    harmonized_sources =
      paste(
        sort(unique(infection_source_harmonized)),
        collapse = ";"
      ),
    binary_pulmonary_source =
      paste(
        sort(unique(binary_pulmonary_source)),
        collapse = ";"
      ),
    source_consistent =
      n_distinct_harmonized_source == 1,
    .groups = "drop"
  )

write_csv(
  patient_source,
  file.path(OUT, "04_CRA002354_PATIENT_SOURCE_CONSISTENCY.csv")
)

patient_category_counts <- patient_source %>%
  filter(source_consistent) %>%
  count(
    harmonized_sources,
    name = "patients"
  ) %>%
  arrange(desc(patients))

write_csv(
  patient_category_counts,
  file.path(OUT, "05_CRA002354_PATIENT_SOURCE_COUNTS.csv")
)

binary_counts <- patient_source %>%
  filter(source_consistent) %>%
  count(
    binary_pulmonary_source,
    name = "patients"
  ) %>%
  arrange(desc(patients))

write_csv(
  binary_counts,
  file.path(OUT, "06_CRA002354_PULMONARY_VS_NONPULMONARY_COUNTS.csv")
)

# ------------------------------------------------------------
# 4. Longitudinal depth audit
# ------------------------------------------------------------

longitudinal_summary <- tibble(
  total_patients = nrow(patient_source),
  total_samples = nrow(cra),
  patients_GE2 = sum(patient_source$n_samples >= 2),
  patients_GE3 = sum(patient_source$n_samples >= 3),
  patients_GE4 = sum(patient_source$n_samples >= 4),
  median_samples_per_patient =
    median(patient_source$n_samples),
  max_samples_per_patient =
    max(patient_source$n_samples),
  earliest_day = min(cra$time_day, na.rm = TRUE),
  latest_day = max(cra$time_day, na.rm = TRUE),
  source_consistent_patients =
    sum(patient_source$source_consistent),
  source_inconsistent_patients =
    sum(!patient_source$source_consistent)
)

write_csv(
  longitudinal_summary,
  file.path(OUT, "07_CRA002354_LONGITUDINAL_DEPTH_SUMMARY.csv")
)

sample_count_distribution <- patient_source %>%
  count(
    n_samples,
    name = "patients"
  ) %>%
  arrange(n_samples)

write_csv(
  sample_count_distribution,
  file.path(OUT, "08_PATIENT_SAMPLE_COUNT_DISTRIBUTION.csv")
)

day_distribution <- cra %>%
  count(
    time_day,
    name = "samples"
  ) %>%
  arrange(time_day)

write_csv(
  day_distribution,
  file.path(OUT, "09_DAY_SAMPLE_DISTRIBUTION.csv")
)

source_longitudinal_depth <- patient_source %>%
  group_by(
    binary_pulmonary_source
  ) %>%
  summarise(
    patients = n(),
    GE2 = sum(n_samples >= 2),
    GE3 = sum(n_samples >= 3),
    median_samples = median(n_samples),
    .groups = "drop"
  )

write_csv(
  source_longitudinal_depth,
  file.path(OUT, "10_SOURCE_STRATIFIED_LONGITUDINAL_DEPTH.csv")
)

# ------------------------------------------------------------
# 5. Local filesystem path/name scan
# ------------------------------------------------------------

all_files <- list.files(
  DATA_ROOT,
  recursive = TRUE,
  full.names = TRUE
)

direct_hits <- all_files[
  str_detect(
    toupper(all_files),
    "CRA002354|PRJCA002241|CRR117"
  )
]

direct_inventory <- tibble(
  file = direct_hits,
  size_bytes = ifelse(
    file.exists(direct_hits),
    file.info(direct_hits)$size,
    NA_real_
  ),
  extension = tolower(
    tools::file_ext(direct_hits)
  ),
  direct_match_type = case_when(
    str_detect(toupper(direct_hits), "CRA002354") ~ "CRA002354_PATH",
    str_detect(toupper(direct_hits), "PRJCA002241") ~ "PRJCA002241_PATH",
    str_detect(toupper(direct_hits), "CRR117") ~ "CRR_RUN_PATH",
    TRUE ~ "OTHER"
  )
)

write_csv(
  direct_inventory,
  file.path(OUT, "11_LOCAL_DIRECT_FILE_INVENTORY.csv")
)

# ------------------------------------------------------------
# 6. Tabular run-overlap scan
# ------------------------------------------------------------
# This catches already-processed tables that do not carry CRA002354
# in the filename/path but contain the CRR run IDs internally.
# ------------------------------------------------------------

tab_files <- list.files(
  DATA_ROOT,
  recursive = TRUE,
  full.names = TRUE,
  pattern = "\\.(csv|tsv|txt)$",
  ignore.case = TRUE
)

# Avoid giant raw report files and restrict to manageable tables.
fi <- file.info(tab_files)

tab_files <- tab_files[
  !is.na(fi$size) &
    fi$size > 0 &
    fi$size <= 200 * 1024^2
]

cra_runs <- unique(na.omit(cra$run_id))

overlap_audit <- list()

for (f in tab_files) {

  x <- safe_read(f)
  if (is.null(x) || nrow(x) == 0 || ncol(x) == 0) next

  # Search all low-cardinality/string-like columns for exact CRR overlap.
  overlaps <- sapply(names(x), function(cc) {
    vals <- toupper(norm_chr(x[[cc]]))
    length(
      intersect(
        unique(vals),
        toupper(cra_runs)
      )
    )
  })

  max_overlap <- max(overlaps, na.rm = TRUE)

  if (!is.finite(max_overlap) || max_overlap == 0) next

  best_col <- names(overlaps)[which.max(overlaps)]

  numeric_cols <- sum(
    vapply(
      x,
      is.numeric,
      logical(1)
    )
  )

  overlap_audit[[length(overlap_audit) + 1]] <- tibble(
    file = f,
    rows = nrow(x),
    cols = ncol(x),
    best_run_column = best_col,
    exact_CRR_run_overlap = max_overlap,
    overlap_fraction_of_131 =
      max_overlap / length(cra_runs),
    numeric_columns = numeric_cols,
    candidate_type =
      classify_candidate_file(f, x),
    columns =
      paste(names(x), collapse = ";")
  )
}

overlap_audit <- bind_rows(overlap_audit)

if (nrow(overlap_audit) > 0) {
  overlap_audit <- overlap_audit %>%
    arrange(
      desc(exact_CRR_run_overlap),
      desc(numeric_columns)
    )
}

write_csv(
  overlap_audit,
  file.path(OUT, "12_CRR_RUN_OVERLAP_TABULAR_AUDIT.csv")
)

# ------------------------------------------------------------
# 7. Processed microbiome readiness
# ------------------------------------------------------------

micro_candidates <- overlap_audit %>%
  filter(
    candidate_type %in% c(
      "GENUS_ABUNDANCE",
      "FEATURE_TABLE",
      "NUMERIC_MICROBIOME_TABLE",
      "R_OBJECT",
      "DOWNSTREAM_RESULT"
    ) |
      numeric_columns >= 10
  )

write_csv(
  micro_candidates,
  file.path(OUT, "13_PROCESSED_MICROBIOME_CANDIDATES.csv")
)

best_micro_overlap <- if (nrow(micro_candidates) > 0) {
  max(
    micro_candidates$exact_CRR_run_overlap,
    na.rm = TRUE
  )
} else 0

best_micro_fraction <- if (length(cra_runs) > 0) {
  best_micro_overlap / length(cra_runs)
} else 0

# Direct raw sequence evidence
raw_sequence_hits <- direct_inventory %>%
  filter(
    str_detect(
      tolower(file),
      "\\.(fastq|fq|fasta|fa)(\\.gz)?$"
    )
  )

# ------------------------------------------------------------
# 8. Feasibility decision
# ------------------------------------------------------------

meaningful_counts <- patient_category_counts %>%
  filter(
    !harmonized_sources %in%
      c(
        "OTHER_UNKNOWN",
        "UNMAPPED"
      )
  ) %>%
  arrange(desc(patients))

multi_source_formal <- (
  nrow(meaningful_counts) >= 3 &&
    meaningful_counts$patients[3] >= 5
)

binary_formal <- (
  nrow(binary_counts) == 2 &&
    min(binary_counts$patients) >= 15 &&
    all(patient_source$source_consistent)
)

longitudinal_formal <- (
  longitudinal_summary$patients_GE2 >= 30 &&
    longitudinal_summary$patients_GE3 >= 15
)

processed_ready <- best_micro_fraction >= 0.80

raw_available <- nrow(raw_sequence_hits) >= 100

rescue_status <- case_when(

  binary_formal &&
    longitudinal_formal &&
    processed_ready ~
    "READY_FOR_PULMONARY_VS_NONPULMONARY_LONGITUDINAL_ANALYSIS",

  binary_formal &&
    longitudinal_formal &&
    raw_available ~
    "SOURCE_ANALYSIS_FEASIBLE_BUT_REPROCESSING_REQUIRED",

  binary_formal &&
    longitudinal_formal ~
    "SOURCE_METADATA_STRONG_BUT_MICROBIOME_TABLE_NOT_LOCATED",

  TRUE ~
    "SOURCE_ANALYSIS_NOT_YET_FORMALLY_FEASIBLE"
)

decision <- tibble(
  cohort = "CRA002354",
  alias = "PRJCA002241",
  samples = nrow(cra),
  patients = n_distinct(cra$patient_id),
  source_consistent_patients =
    sum(patient_source$source_consistent),
  pulmonary_patients =
    binary_counts$patients[
      binary_counts$binary_pulmonary_source ==
        "PULMONARY"
    ][1],
  nonpulmonary_recorded_patients =
    binary_counts$patients[
      binary_counts$binary_pulmonary_source ==
        "NONPULMONARY_RECORDED"
    ][1],
  patients_GE2 =
    longitudinal_summary$patients_GE2,
  patients_GE3 =
    longitudinal_summary$patients_GE3,
  formal_multicategory_source_analysis =
    multi_source_formal,
  formal_binary_pulmonary_nonpulmonary =
    binary_formal,
  formal_longitudinal_depth =
    longitudinal_formal,
  best_processed_run_overlap =
    best_micro_overlap,
  best_processed_run_overlap_fraction =
    best_micro_fraction,
  raw_sequence_files_detected =
    nrow(raw_sequence_hits),
  rescue_status = rescue_status,
  recommended_next_action = case_when(
    rescue_status ==
      "READY_FOR_PULMONARY_VS_NONPULMONARY_LONGITUDINAL_ANALYSIS" ~
      "Integrate the best processed microbiome table with frozen master metadata, QC depth, then run pulmonary-vs-nonpulmonary longitudinal ecological models.",

    rescue_status ==
      "SOURCE_ANALYSIS_FEASIBLE_BUT_REPROCESSING_REQUIRED" ~
      "Process CRA002354 raw 16S data with the current V2 taxonomy pipeline before source-stratified modeling.",

    rescue_status ==
      "SOURCE_METADATA_STRONG_BUT_MICROBIOME_TABLE_NOT_LOCATED" ~
      "Locate or reconstruct CRA002354 abundance/ASV tables before modeling.",

    TRUE ~
      "Retain infection-site data descriptively and do not run formal source-stratified models yet."
  )
)

write_csv(
  decision,
  file.path(OUT, "14_FINAL_CRA002354_RESCUE_DECISION.csv")
)

# ------------------------------------------------------------
# 9. Interpretation
# ------------------------------------------------------------

txt <- c(
  "STEP94A4 CRA002354 RESCUE + READINESS AUDIT",
  "",
  paste0(
    "Observed samples: ",
    nrow(cra)
  ),
  paste0(
    "Observed patients: ",
    n_distinct(cra$patient_id)
  ),
  paste0(
    "Patients with >=2 samples: ",
    longitudinal_summary$patients_GE2
  ),
  paste0(
    "Patients with >=3 samples: ",
    longitudinal_summary$patients_GE3
  ),
  paste0(
    "Source-consistent patients: ",
    sum(patient_source$source_consistent),
    "/",
    nrow(patient_source)
  ),
  "",
  paste0(
    "Pulmonary patients: ",
    binary_counts$patients[
      binary_counts$binary_pulmonary_source ==
        "PULMONARY"
    ][1]
  ),
  paste0(
    "Non-pulmonary recorded patients: ",
    binary_counts$patients[
      binary_counts$binary_pulmonary_source ==
        "NONPULMONARY_RECORDED"
    ][1]
  ),
  "",
  "Important inference rule:",
  "The detailed non-pulmonary source categories are sparse, so respiratory/abdominal/urinary/bloodstream/skin should NOT be modeled as a many-level inferential factor.",
  "A pulmonary-vs-nonpulmonary analysis is potentially feasible because the original cohort/publication explicitly discussed pulmonary versus nonpulmonary infection-site patterns.",
  "",
  paste0(
    "Processed microbiome readiness status: ",
    rescue_status
  ),
  "",
  "Do not run source-stratified microbiome models until the processed-table readiness audit is reviewed."
)

writeLines(
  txt,
  file.path(OUT, "15_STEP94A4_INTERPRETATION.txt")
)

# ------------------------------------------------------------
# 10. QC
# ------------------------------------------------------------

qc <- tibble(
  repaired_project_accession = "CRA002354",
  alias = "PRJCA002241",
  source_rows = nrow(cra),
  patients = n_distinct(cra$patient_id),
  all_patient_sources_consistent =
    all(patient_source$source_consistent),
  GE2_patients =
    longitudinal_summary$patients_GE2,
  GE3_patients =
    longitudinal_summary$patients_GE3,
  day_min =
    longitudinal_summary$earliest_day,
  day_max =
    longitudinal_summary$latest_day,
  direct_files_detected =
    nrow(direct_inventory),
  tabular_files_with_CRR_overlap =
    nrow(overlap_audit),
  processed_microbiome_candidates =
    nrow(micro_candidates),
  best_processed_overlap =
    best_micro_overlap,
  rescue_status =
    rescue_status
)

write_csv(
  qc,
  file.path(OUT, "16_STEP94A4_QC_SUMMARY.csv")
)

writeLines(
  c(
    paste0("Completed: ", Sys.time()),
    paste0("Rescue status: ", rescue_status),
    "STEP94A4 COMPLETE"
  ),
  file.path(OUT, "_STEP94A4_COMPLETE.ok")
)

cat("STEP94A4 COMPLETE\n")
