# ============================================================
# Sepsis V2 - Step94A
# INFECTION-SOURCE METADATA FEASIBILITY AUDIT
#
# Goal
# ----
# Before any infection-source modeling, determine:
# 1) which cohorts actually contain infection-source information;
# 2) which columns are genuine clinical infection-source fields rather
#    than specimen/sample-source fields;
# 3) coverage at sample/patient level;
# 4) whether source categories can be harmonized;
# 5) which cohorts are eligible for formal, supportive, or no analysis.
#
# IMPORTANT
# ---------
# This step is an AUDIT ONLY.
# It does not overwrite metadata and does not force uncertain values into
# infection-source categories.
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

CANONICAL_ROOT <- file.path(
  ROOT,
  "data",
  "_V2_ANALYSIS_READY",
  "01_PROJECTS"
)

DATA_ROOT <- file.path(
  ROOT,
  "data"
)

OUT <- file.path(
  ROOT,
  "results",
  "V2_34A_INFECTION_SOURCE_FEASIBILITY_AUDIT"
)

dir.create(
  OUT,
  recursive = TRUE,
  showWarnings = FALSE
)

TARGET_PROJECTS <- c(
  "PRJNA691455",
  "PRJEB82425",
  "PRJNA516701",
  "PRJNA851469",
  "PRJNA578267",
  "PRJNA1166732",
  "PRJNA430161",
  "PRJNA978257",
  "PRJNA1010969"
)

logfile <- file.path(OUT, "_STEP94A_runtime.txt")
writeLines(paste("START", Sys.time()), logfile)

logmsg <- function(...) {
  z <- paste0(...)
  cat(z, "\n")
  cat(z, "\n", file = logfile, append = TRUE)
}

norm_chr <- function(x) {
  x <- str_squish(as.character(x))
  x[x %in% c("", "NA", "N/A", "NULL", "None", "none", ".", "-")] <- NA_character_
  x
}

detect_project <- function(path) {
  p <- toupper(path)
  hit <- str_extract(
    p,
    "PRJ(?:NA|EB)\\d+"
  )
  ifelse(is.na(hit), "UNKNOWN", hit)
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

candidate_score <- function(colname) {

  z <- tolower(colname)

  score <- 0
  reason <- character()

  # Strong clinical-source patterns
  if (str_detect(
    z,
    "infection.?source|source.?of.?infection|infection.?site|site.?of.?infection|sepsis.?source|primary.?infection|focus.?of.?infection|infectious.?focus"
  )) {
    score <- score + 12
    reason <- c(reason, "strong infection-source field name")
  }

  # Moderate clinical patterns
  if (str_detect(z, "infection")) {
    score <- score + 4
    reason <- c(reason, "contains infection")
  }

  if (str_detect(z, "diagnos|syndrome|admission.?reason|primary.?dx|clinical.?source")) {
    score <- score + 3
    reason <- c(reason, "clinical diagnosis/syndrome field")
  }

  if (str_detect(z, "(^|_)source($|_)|(^|_)site($|_)|focus")) {
    score <- score + 2
    reason <- c(reason, "generic source/site field")
  }

  # Negative evidence: likely biological specimen/sample source,
  # not infection source.
  if (str_detect(
    z,
    "sample.?source|specimen|source.?material|tissue|body.?site|body.?habitat|biome|isolation.?source|sample.?type|sample.?site|collection.?site"
  )) {
    score <- score - 10
    reason <- c(reason, "likely specimen/sample-source field")
  }

  if (str_detect(
    z,
    "run|accession|biosample|bioproject|library|platform|assay|taxonomy|genus|species"
  )) {
    score <- score - 5
    reason <- c(reason, "technical accession/taxonomy field")
  }

  list(
    score = score,
    reason = paste(reason, collapse = "; ")
  )
}

detect_id_col <- function(nms, type = c("patient", "run")) {

  type <- match.arg(type)

  if (type == "patient") {
    priority <- c(
      "patient_id",
      "master_patient_id",
      "Patient_ID",
      "patient_uid",
      "Patient_true",
      "subject_id",
      "Subject_ID",
      "participant_id",
      "patient"
    )
  } else {
    priority <- c(
      "run_id",
      "Run",
      "Run_ID",
      "run",
      "RunID",
      "run_accession"
    )
  }

  exact <- priority[priority %in% nms]

  if (length(exact) > 0) {
    return(exact[1])
  }

  if (type == "patient") {
    hit <- nms[
      str_detect(
        nms,
        regex("patient|subject|participant", ignore_case = TRUE)
      )
    ]
  } else {
    hit <- nms[
      str_detect(
        nms,
        regex("^run$|run.*id|run.*access", ignore_case = TRUE)
      )
    ]
  }

  if (length(hit) == 0) return(NA_character_)
  hit[1]
}

source_category_audit_only <- function(x) {

  z <- tolower(norm_chr(x))

  case_when(
    is.na(z) ~ NA_character_,

    str_detect(
      z,
      "pneum|respirat|lung|pulmon|lower respiratory|upper respiratory|ventilator"
    ) ~ "RESPIRATORY",

    str_detect(
      z,
      "urinary|uti|urine|urosepsis|pyelo|kidney"
    ) ~ "URINARY",

    str_detect(
      z,
      "intra.?abdom|abdom|gastro|intestinal|bowel|periton|biliary|cholang|append|pancrea|liver abscess"
    ) ~ "ABDOMINAL_GI",

    str_detect(
      z,
      "bloodstream|bacteremia|bacteraemia|septicemia|septicaemia|blood culture"
    ) ~ "BLOODSTREAM",

    str_detect(
      z,
      "skin|soft tissue|cellul|wound|fasci|surgical site"
    ) ~ "SKIN_SOFT_TISSUE",

    str_detect(
      z,
      "mening|encephal|central nervous|cns"
    ) ~ "CNS",

    str_detect(
      z,
      "catheter|central line|device|line infection|clabsi"
    ) ~ "DEVICE_LINE",

    str_detect(
      z,
      "unknown|unclear|undetermined|no source|not identified|other"
    ) ~ "UNKNOWN_OTHER",

    TRUE ~ "UNMAPPED_TEXT"
  )
}

# ------------------------------------------------------------
# 1. Build candidate metadata file inventory
# ------------------------------------------------------------

canonical_files <- character()

if (dir.exists(CANONICAL_ROOT)) {
  canonical_files <- list.files(
    CANONICAL_ROOT,
    recursive = TRUE,
    full.names = TRUE,
    pattern = "\\.(csv|tsv|txt)$",
    ignore.case = TRUE
  )
}

# Raw/manual metadata fallback, but avoid scanning every abundance/count file.
raw_files <- list.files(
  DATA_ROOT,
  recursive = TRUE,
  full.names = TRUE,
  pattern = "\\.(csv|tsv|txt)$",
  ignore.case = TRUE
)

raw_files <- raw_files[
  str_detect(
    basename(raw_files),
    regex(
      "metadata|manifest|run.?table|srarun|clinical|phenotype|mapping|sample.?sheet|subject|patient",
      ignore_case = TRUE
    )
  )
]

all_files <- unique(c(canonical_files, raw_files))

# Keep target projects or canonical shared metadata.
file_inventory <- tibble(
  file = all_files,
  project = vapply(all_files, detect_project, character(1)),
  source_layer = case_when(
    str_detect(
      all_files,
      fixed(CANONICAL_ROOT, ignore_case = TRUE)
    ) ~ "CANONICAL_ANALYSIS_READY",
    TRUE ~ "RAW_OR_MANUAL_METADATA"
  )
) %>%
  filter(
    project %in% TARGET_PROJECTS |
      source_layer == "CANONICAL_ANALYSIS_READY"
  )

write_csv(
  file_inventory,
  file.path(OUT, "01_metadata_file_inventory.csv")
)

logmsg("Metadata files inventoried: ", nrow(file_inventory))

# ------------------------------------------------------------
# 2. Scan columns
# ------------------------------------------------------------

column_audit <- list()
value_rows <- list()
sample_rows <- list()

for (i in seq_len(nrow(file_inventory))) {

  f <- file_inventory$file[i]
  project <- file_inventory$project[i]
  layer <- file_inventory$source_layer[i]

  x <- safe_read(f)
  if (is.null(x)) next

  patient_col <- detect_id_col(names(x), "patient")
  run_col <- detect_id_col(names(x), "run")

  for (cc in names(x)) {

    sc <- candidate_score(cc)

    # retain moderate+ candidates, plus weak generic candidates
    # only if their values look clinically informative later.
    if (sc$score < 2) next

    vals <- norm_chr(x[[cc]])
    nonmissing <- sum(!is.na(vals))
    unique_nonmissing <- length(unique(vals[!is.na(vals)]))

    if (nonmissing == 0) next

    mapped_cat <- source_category_audit_only(vals)

    mapped_fraction <- mean(
      mapped_cat %in% c(
        "RESPIRATORY",
        "URINARY",
        "ABDOMINAL_GI",
        "BLOODSTREAM",
        "SKIN_SOFT_TISSUE",
        "CNS",
        "DEVICE_LINE",
        "UNKNOWN_OTHER"
      ),
      na.rm = TRUE
    )

    column_audit[[length(column_audit) + 1]] <- tibble(
      project = project,
      source_layer = layer,
      file = f,
      column = cc,
      candidate_score = sc$score,
      candidate_reason = sc$reason,
      rows = nrow(x),
      nonmissing = nonmissing,
      row_coverage = nonmissing / nrow(x),
      unique_nonmissing = unique_nonmissing,
      auto_mapped_fraction_audit_only = mapped_fraction,
      detected_patient_col = patient_col,
      detected_run_col = run_col
    )

    # value dictionary
    value_tab <- tibble(
      raw_value = vals,
      audit_category = mapped_cat
    ) %>%
      filter(!is.na(raw_value)) %>%
      count(raw_value, audit_category, name = "n") %>%
      arrange(desc(n))

    if (nrow(value_tab) > 0) {
      value_rows[[length(value_rows) + 1]] <- value_tab %>%
        mutate(
          project = project,
          file = f,
          column = cc,
          .before = 1
        )
    }

    # row-level extraction
    row_df <- tibble(
      project = project,
      file = f,
      source_column = cc,
      row_index = seq_len(nrow(x)),
      raw_source_value = vals,
      audit_category = mapped_cat
    )

    if (!is.na(patient_col)) {
      row_df$patient_id_candidate <- norm_chr(x[[patient_col]])
    } else {
      row_df$patient_id_candidate <- NA_character_
    }

    if (!is.na(run_col)) {
      row_df$run_id_candidate <- norm_chr(x[[run_col]])
    } else {
      row_df$run_id_candidate <- NA_character_
    }

    row_df <- row_df %>%
      filter(!is.na(raw_source_value))

    if (nrow(row_df) > 0) {
      sample_rows[[length(sample_rows) + 1]] <- row_df
    }
  }
}

column_audit <- bind_rows(column_audit) %>%
  arrange(
    project,
    desc(candidate_score),
    desc(row_coverage),
    desc(nonmissing)
  )

value_dictionary <- bind_rows(value_rows)

row_candidates <- bind_rows(sample_rows)

write_csv(
  column_audit,
  file.path(OUT, "02_candidate_infection_source_columns.csv")
)

write_csv(
  value_dictionary,
  file.path(OUT, "03_candidate_value_dictionary.csv")
)

write_csv(
  row_candidates,
  file.path(OUT, "04_row_level_infection_source_candidates.csv")
)

# ------------------------------------------------------------
# 3. Identify best candidate field per project
# ------------------------------------------------------------

best_field <- column_audit %>%
  mutate(
    clinical_priority =
      candidate_score +
      5 * (auto_mapped_fraction_audit_only >= 0.50) +
      3 * (row_coverage >= 0.50) +
      2 * (source_layer == "CANONICAL_ANALYSIS_READY")
  ) %>%
  group_by(project) %>%
  arrange(
    desc(clinical_priority),
    desc(candidate_score),
    desc(row_coverage),
    desc(nonmissing)
  ) %>%
  slice_head(n = 1) %>%
  ungroup()

write_csv(
  best_field,
  file.path(OUT, "05_best_candidate_field_by_project.csv")
)

# ------------------------------------------------------------
# 4. Coverage / category feasibility for the best field
# ------------------------------------------------------------

project_feasibility <- list()

for (p in TARGET_PROJECTS) {

  bf <- best_field %>%
    filter(project == p)

  if (nrow(bf) == 0) {

    project_feasibility[[length(project_feasibility) + 1]] <- tibble(
      project = p,
      candidate_found = FALSE,
      best_file = NA_character_,
      best_column = NA_character_,
      nonmissing_rows = 0,
      row_coverage = 0,
      unique_raw_values = 0,
      mapped_categories = 0,
      largest_category_n = 0,
      second_largest_category_n = 0,
      patient_ids_available = FALSE,
      unique_patients_with_source = NA_integer_,
      feasibility_class = "NO_USABLE_INFECTION_SOURCE_FIELD",
      recommended_role = "DO_NOT_MODEL"
    )

    next
  }

  rows <- row_candidates %>%
    filter(
      project == p,
      file == bf$file[1],
      source_column == bf$column[1]
    )

  cats <- rows %>%
    filter(
      !is.na(audit_category),
      audit_category != "UNMAPPED_TEXT"
    ) %>%
    count(audit_category, name = "n") %>%
    arrange(desc(n))

  category_sizes <- cats$n
  largest <- ifelse(length(category_sizes) >= 1, category_sizes[1], 0)
  second <- ifelse(length(category_sizes) >= 2, category_sizes[2], 0)

  patient_available <- any(
    !is.na(rows$patient_id_candidate)
  )

  n_patients <- if (patient_available) {
    n_distinct(
      rows$patient_id_candidate[
        !is.na(rows$patient_id_candidate)
      ]
    )
  } else {
    NA_integer_
  }

  # Conservative feasibility classification.
  formal_ok <- (
    bf$row_coverage[1] >= 0.70 &&
      nrow(cats) >= 2 &&
      second >= 5 &&
      (
        is.na(n_patients) ||
          n_patients >= 20
      )
  )

  supportive_ok <- (
    bf$row_coverage[1] >= 0.40 &&
      nrow(cats) >= 2 &&
      second >= 3 &&
      (
        is.na(n_patients) ||
          n_patients >= 10
      )
  )

  class <- case_when(
    formal_ok ~ "FORMAL_SOURCE_HETEROGENEITY_CANDIDATE",
    supportive_ok ~ "SUPPORTIVE_SOURCE_ANALYSIS_CANDIDATE",
    TRUE ~ "INSUFFICIENT_FOR_INFERENCE"
  )

  role <- case_when(
    class == "FORMAL_SOURCE_HETEROGENEITY_CANDIDATE" ~
      "ELIGIBLE_FOR_FORMAL_WITHIN_COHORT_ANALYSIS",
    class == "SUPPORTIVE_SOURCE_ANALYSIS_CANDIDATE" ~
      "SUPPORTIVE_OR_SENSITIVITY_ONLY",
    TRUE ~
      "DO_NOT_MODEL"
  )

  project_feasibility[[length(project_feasibility) + 1]] <- tibble(
    project = p,
    candidate_found = TRUE,
    best_file = bf$file[1],
    best_column = bf$column[1],
    nonmissing_rows = bf$nonmissing[1],
    row_coverage = bf$row_coverage[1],
    unique_raw_values = bf$unique_nonmissing[1],
    mapped_categories = nrow(cats),
    largest_category_n = largest,
    second_largest_category_n = second,
    patient_ids_available = patient_available,
    unique_patients_with_source = n_patients,
    feasibility_class = class,
    recommended_role = role
  )
}

project_feasibility <- bind_rows(project_feasibility)

write_csv(
  project_feasibility,
  file.path(OUT, "06_project_level_infection_source_feasibility.csv")
)

# ------------------------------------------------------------
# 5. Harmonization suggestions for every observed raw value
# ------------------------------------------------------------

harmonization <- value_dictionary %>%
  transmute(
    project,
    file,
    column,
    raw_value,
    n,
    suggested_category_audit_only = audit_category,
    harmonization_status = case_when(
      audit_category == "UNMAPPED_TEXT" ~ "MANUAL_REVIEW_REQUIRED",
      audit_category == "UNKNOWN_OTHER" ~ "RETAIN_AS_UNKNOWN_OTHER",
      TRUE ~ "AUTO_SUGGESTION_REQUIRES_MANUAL_CONFIRMATION"
    )
  )

write_csv(
  harmonization,
  file.path(OUT, "07_infection_source_harmonization_suggestions.csv")
)

# ------------------------------------------------------------
# 6. Summary counts across categories from best fields
# ------------------------------------------------------------

best_category_summary <- bind_rows(
  lapply(TARGET_PROJECTS, function(p) {

    bf <- best_field %>%
      filter(project == p)

    if (nrow(bf) == 0) return(NULL)

    row_candidates %>%
      filter(
        project == p,
        file == bf$file[1],
        source_column == bf$column[1]
      ) %>%
      count(audit_category, name = "n") %>%
      mutate(project = p, .before = 1)
  })
)

write_csv(
  best_category_summary,
  file.path(OUT, "08_best_field_category_counts.csv")
)

# ------------------------------------------------------------
# 7. Decision matrix
# ------------------------------------------------------------

formal_n <- sum(
  project_feasibility$feasibility_class ==
    "FORMAL_SOURCE_HETEROGENEITY_CANDIDATE"
)

supportive_n <- sum(
  project_feasibility$feasibility_class ==
    "SUPPORTIVE_SOURCE_ANALYSIS_CANDIDATE"
)

decision <- tibble(
  criterion = c(
    "Projects with any candidate infection-source field",
    "Projects eligible for formal source heterogeneity analysis",
    "Projects eligible for supportive source analysis",
    "Can infection source be a major V2 analysis branch now?",
    "Next action"
  ),
  value = c(
    sum(project_feasibility$candidate_found),
    formal_n,
    supportive_n,
    ifelse(
      formal_n >= 2,
      "YES: multi-cohort within-cohort source analyses are potentially feasible",
      ifelse(
        formal_n >= 1 && supportive_n >= 1,
        "LIMITED: one formal cohort plus supportive replication may be feasible",
        "NO/UNCERTAIN: infection-source metadata is currently too sparse for a major inferential branch"
      )
    ),
    "Manually confirm the best-field raw-value dictionary before any modeling"
  )
)

write_csv(
  decision,
  file.path(OUT, "09_STEP94A_DECISION_MATRIX.csv")
)

# ------------------------------------------------------------
# 8. Human-readable interpretation
# ------------------------------------------------------------

formal_projects <- project_feasibility %>%
  filter(
    feasibility_class ==
      "FORMAL_SOURCE_HETEROGENEITY_CANDIDATE"
  ) %>%
  pull(project)

supportive_projects <- project_feasibility %>%
  filter(
    feasibility_class ==
      "SUPPORTIVE_SOURCE_ANALYSIS_CANDIDATE"
  ) %>%
  pull(project)

txt <- c(
  "STEP94A INFECTION-SOURCE FEASIBILITY AUDIT",
  "",
  paste0(
    "Metadata files inventoried: ",
    nrow(file_inventory)
  ),
  paste0(
    "Candidate infection-source columns: ",
    nrow(column_audit)
  ),
  "",
  paste0(
    "Formal candidate cohorts: ",
    ifelse(
      length(formal_projects) == 0,
      "NONE",
      paste(formal_projects, collapse = "; ")
    )
  ),
  paste0(
    "Supportive candidate cohorts: ",
    ifelse(
      length(supportive_projects) == 0,
      "NONE",
      paste(supportive_projects, collapse = "; ")
    )
  ),
  "",
  "CRITICAL RULE:",
  "Do not model infection source yet. First inspect 05_best_candidate_field_by_project.csv and 07_infection_source_harmonization_suggestions.csv.",
  "",
  "The automatic categories are audit suggestions only. They must be manually confirmed against the original cohort metadata/publication before harmonization is frozen.",
  "",
  "Specimen source, body site, sample type, and source_material_id are not automatically treated as clinical infection source."
)

writeLines(
  txt,
  file.path(OUT, "10_STEP94A_INTERPRETATION.txt")
)

# ------------------------------------------------------------
# 9. QC
# ------------------------------------------------------------

qc <- tibble(
  metadata_files = nrow(file_inventory),
  candidate_columns = nrow(column_audit),
  projects_with_candidate = sum(project_feasibility$candidate_found),
  formal_candidate_projects = formal_n,
  supportive_candidate_projects = supportive_n,
  projects_without_usable_field =
    sum(
      project_feasibility$feasibility_class ==
        "NO_USABLE_INFECTION_SOURCE_FIELD"
    ) +
    sum(
      project_feasibility$feasibility_class ==
        "INSUFFICIENT_FOR_INFERENCE"
    )
)

write_csv(
  qc,
  file.path(OUT, "11_STEP94A_QC_SUMMARY.csv")
)

writeLines(
  c(
    paste0("Completed: ", Sys.time()),
    "Infection-source metadata audited; no modeling performed.",
    "STEP94A COMPLETE"
  ),
  file.path(OUT, "_STEP94A_COMPLETE.ok")
)

cat("STEP94A COMPLETE\n")
