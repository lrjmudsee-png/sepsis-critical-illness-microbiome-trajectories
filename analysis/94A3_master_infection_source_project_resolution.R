# ============================================================
# Sepsis V2 - Step94A3
# MASTER infection_source_standard PROJECT/PATIENT RESOLUTION
#
# Why:
# Step94A2 found a genuine standardized clinical field:
#   infection_source_standard
# in the frozen master metadata, with 131 non-missing rows.
#
# It was labeled project=UNKNOWN only because project was inferred from
# the FILE PATH rather than from columns inside the master table.
#
# This step resolves those 131 rows back to:
# - BioProject/cohort
# - patient
# - run/sample
# and determines whether infection-source analysis is actually feasible.
#
# AUDIT ONLY. No metadata are overwritten.
# ============================================================

options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(stringr)
  library(tibble)
})

ROOT <- "E:/sepsis_project"

MASTER <- file.path(
  ROOT,
  "data",
  "_V2_ANALYSIS_READY",
  "00_FREEZE",
  "V2_FINAL_FROZEN_master_metadata_CLINICAL_REPAIRED_20260814_131620.csv"
)

OUT <- file.path(
  ROOT,
  "results",
  "V2_34A3_MASTER_INFECTION_SOURCE_PROJECT_RESOLUTION"
)

dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

logfile <- file.path(OUT, "_STEP94A3_runtime.txt")
writeLines(paste("START", Sys.time()), logfile)

norm_chr <- function(x) {
  z <- str_squish(as.character(x))
  z[z %in% c("", "NA", "N/A", "NULL", "None", "none", ".", "-")] <- NA_character_
  z
}

pick_col <- function(nms, candidates, regex_fallback = NULL) {
  exact <- candidates[candidates %in% nms]
  if (length(exact) > 0) return(exact[1])

  if (!is.null(regex_fallback)) {
    hit <- nms[
      str_detect(
        nms,
        regex(regex_fallback, ignore_case = TRUE)
      )
    ]
    if (length(hit) > 0) return(hit[1])
  }

  NA_character_
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
      "urinary",
      "urine",
      "uti",
      "urosepsis"
    ) ~ "URINARY",

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
      "surgical_site",
      "surgical site",
      "wound",
      "skin",
      "soft_tissue"
    ) ~ "SKIN_SOFT_TISSUE",

    z %in% c(
      "cns",
      "meningitis",
      "encephalitis"
    ) ~ "CNS",

    z %in% c(
      "catheter",
      "device",
      "line"
    ) ~ "DEVICE_LINE",

    z %in% c(
      "other",
      "unknown",
      "undetermined"
    ) ~ "OTHER_UNKNOWN",

    TRUE ~ "UNMAPPED"
  )
}

if (!file.exists(MASTER)) {
  stop(paste0("Missing frozen master metadata: ", MASTER))
}

m <- read_csv(
  MASTER,
  show_col_types = FALSE,
  guess_max = 50000,
  name_repair = "unique"
)

if (!"infection_source_standard" %in% names(m)) {
  stop("infection_source_standard is absent from frozen master metadata.")
}

# ------------------------------------------------------------
# 1. Detect identifiers
# ------------------------------------------------------------

project_col <- pick_col(
  names(m),
  c(
    "project",
    "bioproject",
    "BioProject",
    "project_id",
    "project_accession",
    "study_accession",
    "cohort",
    "dataset"
  ),
  "bioproject|project|study_accession|cohort|dataset"
)

patient_col <- pick_col(
  names(m),
  c(
    "patient_id",
    "master_patient_id",
    "Patient_ID",
    "patient_uid",
    "Patient_true",
    "subject_id",
    "Subject_ID",
    "participant_id",
    "patient"
  ),
  "patient|subject|participant"
)

run_col <- pick_col(
  names(m),
  c(
    "run_id",
    "Run",
    "Run_ID",
    "run",
    "RunID",
    "run_accession"
  ),
  "^run$|run.*id|run.*access"
)

sample_col <- pick_col(
  names(m),
  c(
    "sample_id",
    "Sample_ID",
    "sample",
    "sample_name",
    "Sample_Name",
    "biosample"
  ),
  "sample.?id|sample.?name|biosample"
)

time_col <- pick_col(
  names(m),
  c(
    "time_day",
    "time_raw",
    "timepoint",
    "visit",
    "day"
  ),
  "time|visit|day"
)

id_detection <- tibble(
  field_type = c(
    "project",
    "patient",
    "run",
    "sample",
    "time"
  ),
  detected_column = c(
    project_col,
    patient_col,
    run_col,
    sample_col,
    time_col
  )
)

write_csv(
  id_detection,
  file.path(OUT, "01_IDENTIFIER_COLUMN_DETECTION.csv")
)

if (is.na(project_col)) {
  stop(
    "No project/cohort column detected in master metadata. Inspect 01_IDENTIFIER_COLUMN_DETECTION.csv."
  )
}

# ------------------------------------------------------------
# 2. Resolve the 131 infection-source rows
# ------------------------------------------------------------

src <- tibble(
  row_index = seq_len(nrow(m)),
  project_raw = norm_chr(m[[project_col]]),
  patient_id = if (!is.na(patient_col)) norm_chr(m[[patient_col]]) else NA_character_,
  run_id = if (!is.na(run_col)) norm_chr(m[[run_col]]) else NA_character_,
  sample_id = if (!is.na(sample_col)) norm_chr(m[[sample_col]]) else NA_character_,
  time_value = if (!is.na(time_col)) norm_chr(m[[time_col]]) else NA_character_,
  infection_source_raw = norm_chr(m$infection_source_standard)
) %>%
  mutate(
    project = toupper(
      str_extract(
        toupper(project_raw),
        "PRJ(?:NA|EB)\\d+"
      )
    ),
    infection_source_harmonized =
      harmonize_source(infection_source_raw)
  ) %>%
  filter(!is.na(infection_source_raw))

write_csv(
  src,
  file.path(OUT, "02_RESOLVED_infection_source_rows.csv")
)

# ------------------------------------------------------------
# 3. Raw and harmonized category dictionary
# ------------------------------------------------------------

dictionary <- src %>%
  count(
    infection_source_raw,
    infection_source_harmonized,
    name = "n"
  ) %>%
  arrange(desc(n))

write_csv(
  dictionary,
  file.path(OUT, "03_INFECTION_SOURCE_DICTIONARY.csv")
)

# ------------------------------------------------------------
# 4. Resolve project coverage
# ------------------------------------------------------------

project_row_counts <- src %>%
  count(
    project,
    name = "source_rows"
  ) %>%
  arrange(desc(source_rows))

write_csv(
  project_row_counts,
  file.path(OUT, "04_SOURCE_ROWS_BY_PROJECT.csv")
)

# Full master denominators by project
master_project <- tibble(
  project_raw = norm_chr(m[[project_col]])
) %>%
  mutate(
    project = toupper(
      str_extract(
        toupper(project_raw),
        "PRJ(?:NA|EB)\\d+"
      )
    )
  ) %>%
  count(project, name = "master_rows")

project_coverage <- master_project %>%
  full_join(project_row_counts, by = "project") %>%
  mutate(
    master_rows = coalesce(master_rows, 0L),
    source_rows = coalesce(source_rows, 0L),
    row_coverage = ifelse(
      master_rows > 0,
      source_rows / master_rows,
      NA_real_
    )
  ) %>%
  arrange(desc(source_rows))

write_csv(
  project_coverage,
  file.path(OUT, "05_PROJECT_ROW_COVERAGE.csv")
)

# ------------------------------------------------------------
# 5. Patient-level source consistency
# ------------------------------------------------------------

if (!is.na(patient_col)) {

  patient_src <- src %>%
    filter(
      !is.na(project),
      !is.na(patient_id)
    ) %>%
    group_by(
      project,
      patient_id
    ) %>%
    summarise(
      n_source_rows = n(),
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
      patient_source_consistent =
        n_distinct_harmonized_source == 1,
      .groups = "drop"
    )

  write_csv(
    patient_src,
    file.path(OUT, "06_PATIENT_SOURCE_CONSISTENCY.csv")
  )

  patient_summary <- patient_src %>%
    group_by(project) %>%
    summarise(
      patients_with_source = n(),
      consistent_patients =
        sum(patient_source_consistent),
      inconsistent_patients =
        sum(!patient_source_consistent),
      consistency_fraction =
        mean(patient_source_consistent),
      .groups = "drop"
    )

  write_csv(
    patient_summary,
    file.path(OUT, "07_PATIENT_SOURCE_COVERAGE_AND_CONSISTENCY.csv")
  )

} else {

  patient_src <- tibble()
  patient_summary <- tibble(
    project = character(),
    patients_with_source = integer(),
    consistent_patients = integer(),
    inconsistent_patients = integer(),
    consistency_fraction = numeric()
  )

  write_csv(
    patient_summary,
    file.path(OUT, "07_PATIENT_SOURCE_COVERAGE_AND_CONSISTENCY.csv")
  )
}

# ------------------------------------------------------------
# 6. Patient-level category counts by project
# ------------------------------------------------------------

if (nrow(patient_src) > 0) {

  patient_category_counts <- patient_src %>%
    filter(patient_source_consistent) %>%
    separate_rows(
      harmonized_sources,
      sep = ";"
    ) %>%
    count(
      project,
      harmonized_sources,
      name = "patients"
    ) %>%
    rename(
      infection_source_harmonized =
        harmonized_sources
    ) %>%
    arrange(
      project,
      desc(patients)
    )

} else {

  # fallback sample-row counts if no patient identifier exists
  patient_category_counts <- src %>%
    count(
      project,
      infection_source_harmonized,
      name = "patients"
    )
}

write_csv(
  patient_category_counts,
  file.path(OUT, "08_SOURCE_CATEGORY_COUNTS_BY_PROJECT.csv")
)

# ------------------------------------------------------------
# 7. Determine actual source-analysis feasibility
# ------------------------------------------------------------

projects <- sort(
  unique(
    na.omit(
      c(
        project_coverage$project,
        patient_category_counts$project
      )
    )
  )
)

feas <- bind_rows(
  lapply(projects, function(p) {

    cov <- project_coverage %>%
      filter(project == p)

    pcat <- patient_category_counts %>%
      filter(project == p) %>%
      filter(
        infection_source_harmonized != "OTHER_UNKNOWN",
        infection_source_harmonized != "UNMAPPED"
      ) %>%
      arrange(desc(patients))

    category_n <- nrow(pcat)

    largest <- ifelse(
      category_n >= 1,
      pcat$patients[1],
      0
    )

    second <- ifelse(
      category_n >= 2,
      pcat$patients[2],
      0
    )

    third <- ifelse(
      category_n >= 3,
      pcat$patients[3],
      0
    )

    ps <- patient_summary %>%
      filter(project == p)

    n_patients <- ifelse(
      nrow(ps) == 1,
      ps$patients_with_source,
      NA_integer_
    )

    consistency <- ifelse(
      nrow(ps) == 1,
      ps$consistency_fraction,
      NA_real_
    )

    # Conservative thresholds:
    # Formal: >=3 meaningful source categories; >=5 patients in top 3;
    # >=20 patients with source; >=90% patient-level consistency.
    formal <- (
      category_n >= 3 &&
      third >= 5 &&
      (is.na(n_patients) || n_patients >= 20) &&
      (is.na(consistency) || consistency >= 0.90)
    )

    # Supportive: >=2 meaningful categories; >=5 in second category;
    # >=12 patients; >=90% consistency.
    supportive <- (
      category_n >= 2 &&
      second >= 5 &&
      (is.na(n_patients) || n_patients >= 12) &&
      (is.na(consistency) || consistency >= 0.90)
    )

    class <- case_when(
      formal ~ "FORMAL_INFECTION_SOURCE_HETEROGENEITY_FEASIBLE",
      supportive ~ "SUPPORTIVE_INFECTION_SOURCE_ANALYSIS_FEASIBLE",
      category_n >= 2 ~ "DESCRIPTIVE_ONLY_INFECTION_SOURCE",
      TRUE ~ "INSUFFICIENT_SOURCE_VARIATION"
    )

    tibble(
      project = p,
      master_rows = ifelse(
        nrow(cov) == 1,
        cov$master_rows,
        NA_integer_
      ),
      source_rows = ifelse(
        nrow(cov) == 1,
        cov$source_rows,
        NA_integer_
      ),
      row_coverage = ifelse(
        nrow(cov) == 1,
        cov$row_coverage,
        NA_real_
      ),
      patients_with_source = n_patients,
      patient_consistency_fraction = consistency,
      meaningful_source_categories = category_n,
      largest_category_patients = largest,
      second_category_patients = second,
      third_category_patients = third,
      feasibility_class = class
    )
  })
)

write_csv(
  feas,
  file.path(OUT, "09_PROJECT_LEVEL_SOURCE_FEASIBILITY.csv")
)

# ------------------------------------------------------------
# 8. Final branch decision
# ------------------------------------------------------------

formal_projects <- feas %>%
  filter(
    feasibility_class ==
      "FORMAL_INFECTION_SOURCE_HETEROGENEITY_FEASIBLE"
  ) %>%
  pull(project)

supportive_projects <- feas %>%
  filter(
    feasibility_class ==
      "SUPPORTIVE_INFECTION_SOURCE_ANALYSIS_FEASIBLE"
  ) %>%
  pull(project)

descriptive_projects <- feas %>%
  filter(
    feasibility_class ==
      "DESCRIPTIVE_ONLY_INFECTION_SOURCE"
  ) %>%
  pull(project)

branch_status <- case_when(
  length(formal_projects) >= 2 ~
    "MULTI_COHORT_INFECTION_SOURCE_BRANCH_FEASIBLE",

  length(formal_projects) == 1 ~
    "SINGLE_COHORT_FORMAL_SOURCE_ANALYSIS_FEASIBLE",

  length(supportive_projects) >= 1 ~
    "SUPPORTIVE_SOURCE_ANALYSIS_ONLY",

  length(descriptive_projects) >= 1 ~
    "DESCRIPTIVE_SOURCE_CONTEXT_ONLY",

  TRUE ~
    "INFECTION_SOURCE_BRANCH_NOT_FEASIBLE"
)

decision <- tibble(
  branch_status = branch_status,
  formal_projects =
    paste(formal_projects, collapse = ";"),
  supportive_projects =
    paste(supportive_projects, collapse = ";"),
  descriptive_projects =
    paste(descriptive_projects, collapse = ";"),
  total_source_rows = nrow(src),
  projects_resolved =
    n_distinct(src$project, na.rm = TRUE),
  unresolved_project_rows =
    sum(is.na(src$project)),
  unmapped_source_rows =
    sum(src$infection_source_harmonized == "UNMAPPED"),
  recommended_next_action = case_when(
    branch_status ==
      "MULTI_COHORT_INFECTION_SOURCE_BRANCH_FEASIBLE" ~
      "Proceed to within-cohort source-stratified ecological analyses, then synthesize directions across cohorts.",

    branch_status ==
      "SINGLE_COHORT_FORMAL_SOURCE_ANALYSIS_FEASIBLE" ~
      "Proceed with one formal within-cohort source analysis; do not claim multi-cohort replication.",

    branch_status ==
      "SUPPORTIVE_SOURCE_ANALYSIS_ONLY" ~
      "Run supportive source-stratified analysis only; keep outside the main causal narrative.",

    branch_status ==
      "DESCRIPTIVE_SOURCE_CONTEXT_ONLY" ~
      "Retain source distributions descriptively; do not model inferentially.",

    TRUE ~
      "Close infection-source branch and proceed to final statistical synthesis."
  )
)

write_csv(
  decision,
  file.path(OUT, "10_FINAL_MASTER_SOURCE_BRANCH_DECISION.csv")
)

# ------------------------------------------------------------
# 9. Interpretation
# ------------------------------------------------------------

txt <- c(
  "STEP94A3 MASTER INFECTION-SOURCE PROJECT RESOLUTION",
  "",
  paste0(
    "Frozen master rows: ",
    nrow(m)
  ),
  paste0(
    "Rows with infection_source_standard: ",
    nrow(src)
  ),
  paste0(
    "Resolved BioProjects among source rows: ",
    n_distinct(src$project, na.rm = TRUE)
  ),
  paste0(
    "Rows with unresolved project: ",
    sum(is.na(src$project))
  ),
  "",
  paste0(
    "Final branch status: ",
    branch_status
  ),
  "",
  "Important:",
  "Step94A2's project=UNKNOWN result for infection_source_standard was caused by path-based project inference. Step94A3 resolves project identity from the master metadata itself.",
  "",
  "Do not close the infection-source branch until 10_FINAL_MASTER_SOURCE_BRANCH_DECISION.csv has been reviewed."
)

writeLines(
  txt,
  file.path(OUT, "11_STEP94A3_INTERPRETATION.txt")
)

# ------------------------------------------------------------
# 10. QC
# ------------------------------------------------------------

qc <- tibble(
  master_rows = nrow(m),
  source_rows = nrow(src),
  resolved_project_rows = sum(!is.na(src$project)),
  unresolved_project_rows = sum(is.na(src$project)),
  projects_with_source = n_distinct(src$project, na.rm = TRUE),
  patient_identifier_available = !is.na(patient_col),
  run_identifier_available = !is.na(run_col),
  unmapped_source_rows =
    sum(src$infection_source_harmonized == "UNMAPPED"),
  branch_status = branch_status
)

write_csv(
  qc,
  file.path(OUT, "12_STEP94A3_QC_SUMMARY.csv")
)

writeLines(
  c(
    paste0("Completed: ", Sys.time()),
    paste0("Branch status: ", branch_status),
    "STEP94A3 COMPLETE"
  ),
  file.path(OUT, "_STEP94A3_COMPLETE.ok")
)

cat("STEP94A3 COMPLETE\n")
