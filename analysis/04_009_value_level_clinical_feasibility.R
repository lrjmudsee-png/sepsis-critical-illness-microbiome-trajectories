
# ============================================================
# V2_36D VALUE-LEVEL CLINICAL FEASIBILITY AUDIT
#
# PURPOSE
#   Confirm whether clinical variables are ACTUALLY populated
#   and analysis-ready in the frozen V2 cohorts.
#
# IMPORTANT
#   - No microbiome inferential model is run here.
#   - Step36C column-presence counts are NOT treated as evidence
#     of usable clinical data.
#   - Uses the frozen master metadata as the authoritative layer.
# ============================================================

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(stringr)
  library(purrr)
  library(tidyr)
})

ROOT <- "E:/sepsis_project"
DATA <- file.path(ROOT, "data")
OUT  <- file.path(ROOT, "results", "V2_36D_VALUE_LEVEL_CLINICAL_FEASIBILITY")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------------------
# 1. Locate authoritative frozen master metadata
# ------------------------------------------------------------
master_candidates <- list.files(
  file.path(DATA, "_V2_ANALYSIS_READY", "00_FREEZE"),
  pattern = "V2_FINAL_FROZEN_master_metadata_CLINICAL_REPAIRED.*\\.csv$",
  full.names = TRUE
)

if (length(master_candidates) == 0) {
  stop("Authoritative frozen clinical master metadata not found.")
}

master_file <- master_candidates[which.max(file.info(master_candidates)$mtime)]
cat("MASTER:", master_file, "\n")

dat <- read_csv(master_file, show_col_types = FALSE, progress = FALSE)

required_id <- c("project", "patient_uid", "sample_uid")
missing_id <- setdiff(required_id, names(dat))
if (length(missing_id) > 0) {
  stop("Missing required ID fields: ", paste(missing_id, collapse = ", "))
}

# Keep the known frozen V2 cohorts.
target_projects <- c(
  "PRJNA691455",
  "PRJEB82425",
  "PRJNA516701",
  "PRJNA851469",
  "PRJNA578267",
  "PRJNA430161",
  "PRJNA1166732",
  "PRJNA978257"
)

dat <- dat %>%
  filter(project %in% target_projects)

# Prefer rows that are in the frozen analysis layer when fields exist.
if ("freeze_metadata_include" %in% names(dat)) {
  keep <- dat$freeze_metadata_include
  if (is.logical(keep)) {
    dat <- dat %>% filter(is.na(freeze_metadata_include) | freeze_metadata_include)
  }
}

# ------------------------------------------------------------
# 2. Utility functions
# ------------------------------------------------------------
is_informative <- function(x) {
  y <- str_trim(as.character(x))
  !is.na(x) &
    y != "" &
    !tolower(y) %in% c(
      "na", "n/a", "nan", "null",
      "unknown", "not available", "not_available",
      "not reported", "not_reported", "missing"
    )
}

clean_value <- function(x) {
  y <- str_trim(as.character(x))
  y[!is_informative(y)] <- NA_character_
  y
}

safe_examples <- function(x, n = 8) {
  y <- unique(clean_value(x))
  y <- y[!is.na(y)]
  if (length(y) == 0) return("")
  paste(head(y, n), collapse = " || ")
}

numeric_parse <- function(x) {
  suppressWarnings(parse_number(as.character(x)))
}

# ------------------------------------------------------------
# 3. Value-level inventory
# ------------------------------------------------------------
key_vars <- c(
  "antibiotics",
  "severity",
  "outcome",
  "outcome_28d",
  "baseline_sofa",
  "baseline_apache_ii",
  "baseline_lactate",
  "organ_dysfunction_status",
  "infection_source_standard",
  "intervention",
  "intervention_arm",
  "clinical_metadata_source",
  "clinical_mapping_status"
)
key_vars <- intersect(key_vars, names(dat))

summarize_var <- function(df, project_name, v) {
  x <- df[[v]]
  inf <- is_informative(x)

  p_all <- unique(df$patient_uid[!is.na(df$patient_uid)])
  p_inf <- unique(df$patient_uid[inf & !is.na(df$patient_uid)])

  val <- clean_value(x)
  uniq <- unique(val[!is.na(val)])

  num <- numeric_parse(val)
  numeric_fraction <- if (sum(inf) == 0) NA_real_ else sum(!is.na(num) & inf) / sum(inf)

  tibble(
    project = project_name,
    variable = v,
    samples_total = nrow(df),
    samples_informative = sum(inf),
    sample_coverage_pct = round(100 * mean(inf), 2),
    patients_total = length(p_all),
    patients_informative = length(p_inf),
    patient_coverage_pct = ifelse(length(p_all) == 0, NA_real_,
                                  round(100 * length(p_inf) / length(p_all), 2)),
    n_unique_nonmissing = length(uniq),
    has_variation = length(uniq) >= 2,
    numeric_fraction = round(numeric_fraction, 3),
    numeric_n = sum(!is.na(num)),
    numeric_median = ifelse(sum(!is.na(num)) > 0, median(num, na.rm = TRUE), NA_real_),
    numeric_q1 = ifelse(sum(!is.na(num)) > 0, quantile(num, .25, na.rm = TRUE), NA_real_),
    numeric_q3 = ifelse(sum(!is.na(num)) > 0, quantile(num, .75, na.rm = TRUE), NA_real_),
    example_values = safe_examples(x)
  )
}

inventory <- map_dfr(target_projects, function(prj) {
  d <- dat %>% filter(project == prj)
  if (nrow(d) == 0) return(tibble())
  map_dfr(key_vars, ~summarize_var(d, prj, .x))
})

write_csv(
  inventory,
  file.path(OUT, "01_value_level_clinical_variable_inventory.csv")
)

# ------------------------------------------------------------
# 4. Patient-level longitudinal eligibility
# ------------------------------------------------------------
# Calculate how many patients actually have >=2 distinct timepoints.
time_col <- if ("analysis_time_order" %in% names(dat)) {
  "analysis_time_order"
} else if ("time_day" %in% names(dat)) {
  "time_day"
} else if ("time_order" %in% names(dat)) {
  "time_order"
} else {
  NA_character_
}

if (!is.na(time_col)) {
  long_status <- dat %>%
    mutate(.time = .data[[time_col]]) %>%
    group_by(project, patient_uid) %>%
    summarise(
      n_time = n_distinct(.time[!is.na(.time)]),
      repeated = n_time >= 2,
      .groups = "drop"
    )
} else {
  long_status <- dat %>%
    distinct(project, patient_uid) %>%
    mutate(n_time = NA_integer_, repeated = NA)
}

# ------------------------------------------------------------
# 5. Analysis-specific feasibility
# ------------------------------------------------------------
patient_var_summary <- function(df, v) {
  if (!v %in% names(df)) {
    return(tibble(
      patients_informative = 0,
      repeated_patients_informative = 0,
      n_unique = 0,
      numeric_patients = 0
    ))
  }

  tmp <- df %>%
    mutate(.value = clean_value(.data[[v]])) %>%
    group_by(project, patient_uid) %>%
    summarise(
      value = first(.value[!is.na(.value)], default = NA_character_),
      .groups = "drop"
    ) %>%
    left_join(long_status, by = c("project", "patient_uid"))

  tibble(
    patients_informative = sum(!is.na(tmp$value)),
    repeated_patients_informative = sum(!is.na(tmp$value) & tmp$repeated %in% TRUE),
    n_unique = n_distinct(tmp$value[!is.na(tmp$value)]),
    numeric_patients = sum(!is.na(numeric_parse(tmp$value)))
  )
}

feas_list <- list()

for (prj in target_projects) {
  d <- dat %>% filter(project == prj)
  if (nrow(d) == 0) next

  # Antibiotic adjustment
  a <- patient_var_summary(d, "antibiotics")
  antibiotic_ready <- (
    a$patients_informative >= 10 &&
    a$repeated_patients_informative >= 8 &&
    a$n_unique >= 2
  )

  feas_list[[length(feas_list)+1]] <- tibble(
    project = prj,
    analysis = "ANTIBIOTIC_ADJUSTED_LONGITUDINAL_DISPLACEMENT",
    variable = "antibiotics",
    patients_informative = a$patients_informative,
    repeated_patients_informative = a$repeated_patients_informative,
    n_unique_nonmissing = a$n_unique,
    preliminary_ready = antibiotic_ready,
    reason = ifelse(
      antibiotic_ready,
      "Patient-linked antibiotic field has coverage, repeated patients, and variation.",
      "Insufficient populated patient-linked antibiotic coverage, repeated patients, or variation."
    )
  )

  # Baseline SOFA
  s <- patient_var_summary(d, "baseline_sofa")
  sofa_ready <- (
    s$numeric_patients >= 10 &&
    s$repeated_patients_informative >= 8 &&
    s$n_unique >= 3
  )

  feas_list[[length(feas_list)+1]] <- tibble(
    project = prj,
    analysis = "BASELINE_SOFA_VS_LONGITUDINAL_ECOLOGICAL_DISPLACEMENT",
    variable = "baseline_sofa",
    patients_informative = s$patients_informative,
    repeated_patients_informative = s$repeated_patients_informative,
    n_unique_nonmissing = s$n_unique,
    preliminary_ready = sofa_ready,
    reason = ifelse(
      sofa_ready,
      "Numeric baseline SOFA is populated in enough longitudinally sampled patients.",
      "Insufficient numeric baseline SOFA among longitudinally sampled patients."
    )
  )

  # Baseline APACHE II
  ap <- patient_var_summary(d, "baseline_apache_ii")
  apache_ready <- (
    ap$numeric_patients >= 10 &&
    ap$repeated_patients_informative >= 8 &&
    ap$n_unique >= 3
  )

  feas_list[[length(feas_list)+1]] <- tibble(
    project = prj,
    analysis = "BASELINE_APACHEII_VS_LONGITUDINAL_ECOLOGICAL_DISPLACEMENT",
    variable = "baseline_apache_ii",
    patients_informative = ap$patients_informative,
    repeated_patients_informative = ap$repeated_patients_informative,
    n_unique_nonmissing = ap$n_unique,
    preliminary_ready = apache_ready,
    reason = ifelse(
      apache_ready,
      "Numeric baseline APACHE II is populated in enough longitudinally sampled patients.",
      "Insufficient numeric baseline APACHE II among longitudinally sampled patients."
    )
  )

  # Baseline lactate
  l <- patient_var_summary(d, "baseline_lactate")
  lactate_ready <- (
    l$numeric_patients >= 10 &&
    l$repeated_patients_informative >= 8 &&
    l$n_unique >= 3
  )

  feas_list[[length(feas_list)+1]] <- tibble(
    project = prj,
    analysis = "BASELINE_LACTATE_VS_LONGITUDINAL_ECOLOGICAL_DISPLACEMENT",
    variable = "baseline_lactate",
    patients_informative = l$patients_informative,
    repeated_patients_informative = l$repeated_patients_informative,
    n_unique_nonmissing = l$n_unique,
    preliminary_ready = lactate_ready,
    reason = ifelse(
      lactate_ready,
      "Numeric baseline lactate is populated in enough longitudinally sampled patients.",
      "Insufficient numeric baseline lactate among longitudinally sampled patients."
    )
  )

  # 28-day outcome
  o <- patient_var_summary(d, "outcome_28d")
  out_ready <- (
    o$patients_informative >= 15 &&
    o$repeated_patients_informative >= 10 &&
    o$n_unique >= 2
  )

  feas_list[[length(feas_list)+1]] <- tibble(
    project = prj,
    analysis = "OUTCOME28D_VS_LONGITUDINAL_ECOLOGICAL_DISPLACEMENT",
    variable = "outcome_28d",
    patients_informative = o$patients_informative,
    repeated_patients_informative = o$repeated_patients_informative,
    n_unique_nonmissing = o$n_unique,
    preliminary_ready = out_ready,
    reason = ifelse(
      out_ready,
      "28-day outcome has enough linked longitudinal patients and outcome variation.",
      "Insufficient linked 28-day outcome coverage, longitudinal patients, or outcome variation."
    )
  )
}

feasibility <- bind_rows(feas_list)

write_csv(
  feasibility,
  file.path(OUT, "02_analysis_specific_feasibility.csv")
)

# ------------------------------------------------------------
# 6. Inspect antibiotic detail files anywhere in data tree
# ------------------------------------------------------------
all_csv <- list.files(DATA, pattern = "\\.csv$", recursive = TRUE, full.names = TRUE)

detail_terms <- c(
  "what type of antimicrobial",
  "dol_start",
  "dol_stop",
  "date started",
  "date stopped",
  "ampicillin",
  "gentamicin",
  "vancomycin",
  "meropenem",
  "cefazolin",
  "cefotaxime",
  "cefepime",
  "clindamycin"
)

scan_detail_file <- function(f) {
  hdr <- tryCatch(
    names(read_csv(f, n_max = 0, show_col_types = FALSE, progress = FALSE)),
    error = function(e) character()
  )
  low <- tolower(hdr)
  hits <- hdr[str_detect(low, paste(detail_terms, collapse = "|"))]
  if (length(hits) == 0) return(tibble())

  prj_match <- str_extract(
    f,
    "PRJNA[0-9]+|PRJEB[0-9]+|CRA[0-9]+"
  )

  tibble(
    project_from_path = prj_match,
    file = f,
    n_antimicrobial_detail_columns = length(hits),
    antimicrobial_detail_columns = paste(hits, collapse = " | ")
  )
}

detail_inventory <- map_dfr(all_csv, scan_detail_file)

write_csv(
  detail_inventory,
  file.path(OUT, "03_antimicrobial_detail_file_inventory.csv")
)

# ------------------------------------------------------------
# 7. Human-readable decision summary
# ------------------------------------------------------------
ready <- feasibility %>% filter(preliminary_ready)

summary_lines <- c(
  "V2 STEP36D VALUE-LEVEL CLINICAL FEASIBILITY",
  paste0("Authoritative master: ", master_file),
  "",
  paste0("Frozen target cohorts checked: ", length(unique(dat$project))),
  paste0("Candidate analyses passing preliminary value-level criteria: ", nrow(ready)),
  ""
)

if (nrow(ready) == 0) {
  summary_lines <- c(
    summary_lines,
    "DECISION: No clinical enhancement analysis should be launched yet from the frozen master metadata.",
    "Do not infer usability from column presence alone."
  )
} else {
  summary_lines <- c(
    summary_lines,
    "PRELIMINARY READY CANDIDATES:",
    paste0(
      ready$project, " :: ", ready$analysis,
      " :: n patients=", ready$patients_informative,
      " :: repeated=", ready$repeated_patients_informative,
      " :: unique=", ready$n_unique_nonmissing
    ),
    "",
    "These candidates still require model-specific QC before manuscript use."
  )
}

writeLines(
  summary_lines,
  file.path(OUT, "04_DECISION_SUMMARY.txt")
)

writeLines(
  c(
    paste0("Completed: ", Sys.time()),
    "No microbiome inferential model was run.",
    "STEP36D COMPLETE"
  ),
  file.path(OUT, "_STEP36D_COMPLETE.txt")
)

cat("\nSTEP36D COMPLETE\n")
cat("Output:", OUT, "\n")
