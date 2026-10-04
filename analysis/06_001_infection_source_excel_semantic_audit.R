# ============================================================
# Sepsis V2 - Step94A2
# INFECTION-SOURCE EXCEL-INCLUSIVE + SEMANTIC AUDIT
#
# Why:
# Step94A found only one "formal" candidate (PRJEB82425 infection_group),
# but:
#   - that field is VAP / other infection / control phenotype, not a
#     detailed clinical infection-source taxonomy;
#   - Step94A did not scan XLS/XLSX workbooks.
#
# Step94A2 therefore:
#   1) scans CSV/TSV/TXT + XLS/XLSX metadata files;
#   2) scans every spreadsheet sheet;
#   3) scores both column NAMES and VALUES;
#   4) distinguishes:
#      CLINICAL_INFECTION_SOURCE
#      INFECTION_PHENOTYPE
#      CASE_CONTROL_DISEASE_GROUP
#      SAMPLE_SPECIMEN_FIELD
#      TECHNICAL_FIELD
#      UNKNOWN_REVIEW
#   5) produces a final feasibility decision for the infection-source branch.
#
# AUDIT ONLY. No metadata are overwritten.
# ============================================================

options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(stringr)
  library(tibble)
  library(purrr)
})

if (!requireNamespace("readxl", quietly = TRUE)) {
  stop(
    "Step94A2 requires package 'readxl' because XLS/XLSX files must be audited."
  )
}

ROOT <- "E:/sepsis_project"

DATA_ROOT <- file.path(ROOT, "data")

OUT <- file.path(
  ROOT,
  "results",
  "V2_34A2_INFECTION_SOURCE_EXCEL_SEMANTIC_AUDIT"
)

dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

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

logfile <- file.path(OUT, "_STEP94A2_runtime.txt")
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

detect_project <- function(path) {
  z <- toupper(path)
  p <- str_extract(z, "PRJ(?:NA|EB)\\d+")
  ifelse(is.na(p), "UNKNOWN", p)
}

# ------------------------------------------------------------
# Infection-source value classifier
# ------------------------------------------------------------

source_category <- function(x) {

  z <- tolower(norm_chr(x))

  case_when(
    is.na(z) ~ NA_character_,

    str_detect(
      z,
      "pneum|respirat|pulmon|lung infection|lower respiratory|upper respiratory|vap|ventilator.associated"
    ) ~ "RESPIRATORY",

    str_detect(
      z,
      "urinary|uti|urosepsis|pyelo|cystitis|urinary tract"
    ) ~ "URINARY",

    str_detect(
      z,
      "intra.?abdom|abdominal|periton|biliary|cholang|cholecyst|appendic|bowel perfor|gastrointestinal|gi infection|pancrea|liver abscess"
    ) ~ "ABDOMINAL_GI",

    str_detect(
      z,
      "bloodstream|bacteremia|bacteraemia|septicemia|septicaemia|primary bacter|blood stream"
    ) ~ "BLOODSTREAM",

    str_detect(
      z,
      "skin|soft tissue|cellul|necrotizing|fasci|wound infection|surgical site"
    ) ~ "SKIN_SOFT_TISSUE",

    str_detect(
      z,
      "mening|encephal|central nervous|cns infection"
    ) ~ "CNS",

    str_detect(
      z,
      "catheter|central line|clabsi|device infection|line infection"
    ) ~ "DEVICE_LINE",

    str_detect(
      z,
      "bone|joint|osteomy|septic arthritis"
    ) ~ "BONE_JOINT",

    str_detect(
      z,
      "endocard"
    ) ~ "ENDOCARDITIS",

    str_detect(
      z,
      "unknown source|source unknown|undetermined source|no source identified"
    ) ~ "UNKNOWN_SOURCE",

    TRUE ~ "NO_SOURCE_KEYWORD"
  )
}

column_name_score <- function(cc) {

  z <- tolower(cc)
  score <- 0
  reasons <- character()

  if (str_detect(
    z,
    "source.?of.?infection|infection.?source|infection.?site|site.?of.?infection|infectious.?focus|focus.?of.?infection|sepsis.?source|primary.?infection.?site|infection.?origin"
  )) {
    score <- score + 15
    reasons <- c(reasons, "explicit infection-source/site field")
  }

  if (str_detect(
    z,
    "infection|infectious|sepsis"
  )) {
    score <- score + 4
    reasons <- c(reasons, "infection/sepsis field")
  }

  if (str_detect(
    z,
    "diagnos|diagnosis|disease|condition|phenotype|syndrome|admission.?reason|primary.?dx|clinical"
  )) {
    score <- score + 3
    reasons <- c(reasons, "clinical disease/diagnosis field")
  }

  if (str_detect(
    z,
    "(^|_)source($|_)|(^|_)site($|_)|origin|focus"
  )) {
    score <- score + 2
    reasons <- c(reasons, "generic source/site/origin/focus")
  }

  if (str_detect(
    z,
    "sample.?source|specimen|source.?material|body.?site|body.?habitat|sample.?type|sample.?site|collection.?site|tissue|biome|isolation.?source"
  )) {
    score <- score - 12
    reasons <- c(reasons, "specimen/sample-source penalty")
  }

  if (str_detect(
    z,
    "run|accession|biosample|bioproject|library|platform|assay|taxonomy|genus|species|read|sequence"
  )) {
    score <- score - 8
    reasons <- c(reasons, "technical field penalty")
  }

  list(
    score = score,
    reasons = paste(reasons, collapse = "; ")
  )
}

semantic_classify <- function(
  column_name,
  values,
  source_cats,
  name_score
) {

  vals <- tolower(unique(norm_chr(values)))
  vals <- vals[!is.na(vals)]

  cats <- unique(source_cats)
  cats <- cats[
    !is.na(cats) &
      cats != "NO_SOURCE_KEYWORD"
  ]

  n_real_source_categories <- length(
    setdiff(cats, "UNKNOWN_SOURCE")
  )

  joined <- paste(vals, collapse = " | ")

  # Technical/specimen first.
  z <- tolower(column_name)

  if (str_detect(
    z,
    "sample.?source|specimen|source.?material|body.?site|sample.?type|sample.?site|tissue|biome|isolation.?source"
  )) {
    return("SAMPLE_SPECIMEN_FIELD")
  }

  if (str_detect(
    z,
    "run|accession|biosample|bioproject|library|platform|assay|taxonomy|read"
  )) {
    return("TECHNICAL_FIELD")
  }

  # Classic case-control status labels are NOT infection source.
  case_control_tokens <- c(
    "sepsis",
    "healthy",
    "healthy_control",
    "control",
    "non_sepsis",
    "non_sepsis_icu",
    "case",
    "trauma"
  )

  normalized_vals <- str_replace_all(vals, "[^a-z0-9]+", "_")

  if (
    length(vals) <= 6 &&
    all(
      normalized_vals %in%
        c(
          case_control_tokens,
          "healthy_control",
          "non_sepsis_icu"
        )
    )
  ) {
    return("CASE_CONTROL_DISEASE_GROUP")
  }

  # PRJEB-like infection phenotype: pneumonia/VAP vs other infection vs control.
  if (
    str_detect(joined, "pneum|vap") &&
    str_detect(joined, "other.?infection") &&
    str_detect(joined, "control|uninfect")
  ) {
    return("INFECTION_PHENOTYPE")
  }

  # True source candidate requires at least two distinct anatomic/source
  # categories (e.g. respiratory + urinary/abdominal/bloodstream).
  if (
    n_real_source_categories >= 2 &&
    (
      name_score >= 3 ||
      sum(source_cats != "NO_SOURCE_KEYWORD", na.rm = TRUE) >= 3
    )
  ) {
    return("CLINICAL_INFECTION_SOURCE")
  }

  if (
    n_real_source_categories == 1 &&
    str_detect(joined, "other.?infection|uninfect|control")
  ) {
    return("INFECTION_PHENOTYPE")
  }

  if (
    name_score >= 3 ||
    sum(source_cats != "NO_SOURCE_KEYWORD", na.rm = TRUE) >= 3
  ) {
    return("UNKNOWN_REVIEW")
  }

  "NOT_RELEVANT"
}

# ------------------------------------------------------------
# Read one file/sheet
# ------------------------------------------------------------

safe_read_text <- function(f) {

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

# ------------------------------------------------------------
# 1. File inventory including XLS/XLSX
# ------------------------------------------------------------

all_files <- list.files(
  DATA_ROOT,
  recursive = TRUE,
  full.names = TRUE,
  pattern = "\\.(csv|tsv|txt|xlsx|xls)$",
  ignore.case = TRUE
)

# Metadata-like files; deliberately broader than Step94A.
keep_file <- str_detect(
  basename(all_files),
  regex(
    "metadata|manifest|run.?table|srarun|clinical|phenotype|mapping|sample.?sheet|subject|patient|analysis.?ready|master|cohort|characteristic|demographic",
    ignore_case = TRUE
  )
)

all_files <- all_files[keep_file]

# Avoid obvious abundance/count/taxonomy matrices.
all_files <- all_files[
  !str_detect(
    basename(all_files),
    regex(
      "abundance|taxonomy|otu|asv|feature.?table|count.?matrix|bray|distance|alpha|genus_relative",
      ignore_case = TRUE
    )
  )
]

inventory <- tibble(
  file = all_files,
  project = vapply(all_files, detect_project, character(1)),
  extension = tolower(tools::file_ext(all_files))
) %>%
  filter(
    project %in% TARGET_PROJECTS |
      str_detect(
        file,
        regex("_V2_ANALYSIS_READY", ignore_case = TRUE)
      )
  )

write_csv(
  inventory,
  file.path(OUT, "01_EXCEL_INCLUSIVE_metadata_file_inventory.csv")
)

logmsg("Files inventoried: ", nrow(inventory))

# ------------------------------------------------------------
# 2. Scan every text file and every Excel sheet
# ------------------------------------------------------------

sheet_audit <- list()
column_audit <- list()
value_dict <- list()

scan_table <- function(
  x,
  file,
  sheet,
  project
) {

  if (is.null(x) || nrow(x) == 0 || ncol(x) == 0) {
    return(NULL)
  }

  local_cols <- list()
  local_vals <- list()

  for (cc in names(x)) {

    vals <- norm_chr(x[[cc]])
    nonmissing <- sum(!is.na(vals))

    if (nonmissing == 0) next

    uniq <- length(unique(vals[!is.na(vals)]))

    ns <- column_name_score(cc)
    cats <- source_category(vals)

    source_hit_n <- sum(
      cats != "NO_SOURCE_KEYWORD" &
        !is.na(cats)
    )

    source_hit_fraction <- source_hit_n / nonmissing

    sem <- semantic_classify(
      cc,
      vals,
      cats,
      ns$score
    )

    # Save potentially relevant columns.
    if (
      ns$score >= 2 ||
      source_hit_n >= 3 ||
      sem %in% c(
        "CLINICAL_INFECTION_SOURCE",
        "INFECTION_PHENOTYPE",
        "CASE_CONTROL_DISEASE_GROUP",
        "UNKNOWN_REVIEW"
      )
    ) {

      local_cols[[length(local_cols) + 1]] <- tibble(
        project = project,
        file = file,
        sheet = sheet,
        column = cc,
        rows = nrow(x),
        nonmissing = nonmissing,
        row_coverage = nonmissing / nrow(x),
        unique_nonmissing = uniq,
        name_score = ns$score,
        name_score_reason = ns$reasons,
        source_keyword_hits = source_hit_n,
        source_keyword_fraction = source_hit_fraction,
        semantic_class = sem,
        distinct_source_categories = length(
          unique(
            cats[
              !is.na(cats) &
                cats != "NO_SOURCE_KEYWORD"
            ]
          )
        )
      )

      vv <- tibble(
        raw_value = vals,
        source_category = cats
      ) %>%
        filter(!is.na(raw_value)) %>%
        count(
          raw_value,
          source_category,
          name = "n"
        ) %>%
        arrange(desc(n))

      if (nrow(vv) > 0) {
        local_vals[[length(local_vals) + 1]] <- vv %>%
          mutate(
            project = project,
            file = file,
            sheet = sheet,
            column = cc,
            semantic_class = sem,
            .before = 1
          )
      }
    }
  }

  list(
    columns = bind_rows(local_cols),
    values = bind_rows(local_vals)
  )
}

for (i in seq_len(nrow(inventory))) {

  f <- inventory$file[i]
  project <- inventory$project[i]
  ext <- inventory$extension[i]

  if (ext %in% c("xlsx", "xls")) {

    sheets <- tryCatch(
      readxl::excel_sheets(f),
      error = function(e) character()
    )

    if (length(sheets) == 0) {
      sheet_audit[[length(sheet_audit) + 1]] <- tibble(
        project = project,
        file = f,
        sheet = NA_character_,
        readable = FALSE,
        rows = NA_integer_,
        cols = NA_integer_
      )
      next
    }

    for (sh in sheets) {

      x <- tryCatch(
        readxl::read_excel(
          f,
          sheet = sh,
          .name_repair = "unique"
        ),
        error = function(e) NULL
      )

      sheet_audit[[length(sheet_audit) + 1]] <- tibble(
        project = project,
        file = f,
        sheet = sh,
        readable = !is.null(x),
        rows = ifelse(is.null(x), NA_integer_, nrow(x)),
        cols = ifelse(is.null(x), NA_integer_, ncol(x))
      )

      if (is.null(x)) next

      ans <- scan_table(
        x,
        f,
        sh,
        project
      )

      if (!is.null(ans)) {
        if (nrow(ans$columns) > 0) {
          column_audit[[length(column_audit)+1]] <- ans$columns
        }
        if (nrow(ans$values) > 0) {
          value_dict[[length(value_dict)+1]] <- ans$values
        }
      }
    }

  } else {

    x <- safe_read_text(f)

    sheet_audit[[length(sheet_audit) + 1]] <- tibble(
      project = project,
      file = f,
      sheet = "TEXT_TABLE",
      readable = !is.null(x),
      rows = ifelse(is.null(x), NA_integer_, nrow(x)),
      cols = ifelse(is.null(x), NA_integer_, ncol(x))
    )

    if (is.null(x)) next

    ans <- scan_table(
      x,
      f,
      "TEXT_TABLE",
      project
    )

    if (!is.null(ans)) {
      if (nrow(ans$columns) > 0) {
        column_audit[[length(column_audit)+1]] <- ans$columns
      }
      if (nrow(ans$values) > 0) {
        value_dict[[length(value_dict)+1]] <- ans$values
      }
    }
  }
}

sheet_audit <- bind_rows(sheet_audit)
column_audit <- bind_rows(column_audit)
value_dict <- bind_rows(value_dict)

write_csv(
  sheet_audit,
  file.path(OUT, "02_sheet_read_audit.csv")
)

write_csv(
  column_audit,
  file.path(OUT, "03_SEMANTIC_candidate_column_audit.csv")
)

write_csv(
  value_dict,
  file.path(OUT, "04_SEMANTIC_candidate_value_dictionary.csv")
)

# ------------------------------------------------------------
# 3. True clinical infection-source candidates only
# ------------------------------------------------------------

true_source <- column_audit %>%
  filter(
    semantic_class == "CLINICAL_INFECTION_SOURCE"
  ) %>%
  arrange(
    project,
    desc(row_coverage),
    desc(distinct_source_categories),
    desc(source_keyword_fraction)
  )

write_csv(
  true_source,
  file.path(OUT, "05_TRUE_CLINICAL_INFECTION_SOURCE_CANDIDATES.csv")
)

phenotype_fields <- column_audit %>%
  filter(
    semantic_class == "INFECTION_PHENOTYPE"
  ) %>%
  arrange(project, desc(row_coverage))

write_csv(
  phenotype_fields,
  file.path(OUT, "06_INFECTION_PHENOTYPE_NOT_SOURCE_FIELDS.csv")
)

case_control_fields <- column_audit %>%
  filter(
    semantic_class == "CASE_CONTROL_DISEASE_GROUP"
  ) %>%
  arrange(project, desc(row_coverage))

write_csv(
  case_control_fields,
  file.path(OUT, "07_CASE_CONTROL_NOT_SOURCE_FIELDS.csv")
)

manual_review <- column_audit %>%
  filter(
    semantic_class == "UNKNOWN_REVIEW"
  ) %>%
  arrange(
    project,
    desc(name_score),
    desc(source_keyword_hits)
  )

write_csv(
  manual_review,
  file.path(OUT, "08_MANUAL_REVIEW_REQUIRED_FIELDS.csv")
)

# ------------------------------------------------------------
# 4. Project-level true source feasibility
# ------------------------------------------------------------

project_summary <- bind_rows(
  lapply(TARGET_PROJECTS, function(p) {

    ts <- true_source %>%
      filter(project == p)

    ph <- phenotype_fields %>%
      filter(project == p)

    mr <- manual_review %>%
      filter(project == p)

    n_true <- nrow(ts)
    n_ph <- nrow(ph)
    n_review <- nrow(mr)

    best_true_categories <- if (n_true > 0) {
      max(ts$distinct_source_categories, na.rm = TRUE)
    } else 0

    best_true_coverage <- if (n_true > 0) {
      max(ts$row_coverage, na.rm = TRUE)
    } else 0

    feasibility <- case_when(
      n_true > 0 &&
        best_true_categories >= 3 &&
        best_true_coverage >= 0.70 ~
        "FORMAL_INFECTION_SOURCE_ANALYSIS_FEASIBLE",

      n_true > 0 &&
        best_true_categories >= 2 &&
        best_true_coverage >= 0.40 ~
        "SUPPORTIVE_INFECTION_SOURCE_ANALYSIS_FEASIBLE",

      n_ph > 0 ~
        "INFECTION_PHENOTYPE_ONLY_NOT_SOURCE",

      n_review > 0 ~
        "MANUAL_REVIEW_BEFORE_DECISION",

      TRUE ~
        "NO_INFECTION_SOURCE_INFORMATION"
    )

    tibble(
      project = p,
      true_source_candidate_fields = n_true,
      infection_phenotype_fields = n_ph,
      manual_review_fields = n_review,
      best_true_source_categories = best_true_categories,
      best_true_source_coverage = best_true_coverage,
      final_feasibility = feasibility
    )
  })
)

write_csv(
  project_summary,
  file.path(OUT, "09_PROJECT_LEVEL_TRUE_SOURCE_FEASIBILITY.csv")
)

# ------------------------------------------------------------
# 5. PRJEB82425 publication-consistent phenotype validation
# ------------------------------------------------------------

prjeb_values <- value_dict %>%
  filter(
    project == "PRJEB82425",
    column == "infection_group"
  )

prjeb_expected <- tibble(
  raw_value = c(
    "Pneumonia",
    "Other_infection",
    "Control"
  ),
  publication_interpretation = c(
    "VENTILATOR_ASSOCIATED_PNEUMONIA_VAP",
    "OTHER_INFECTIONS",
    "UNINFECTED_CONTROL"
  ),
  expected_patient_n_from_publication = c(
    17,
    12,
    9
  )
)

# Step94A row-level patient counts if available.
STEP94A_ROWS <- file.path(
  ROOT,
  "results",
  "V2_34A_INFECTION_SOURCE_FEASIBILITY_AUDIT",
  "04_row_level_infection_source_candidates.csv"
)

prjeb_validation <- prjeb_expected

if (file.exists(STEP94A_ROWS)) {

  rr <- read_csv(
    STEP94A_ROWS,
    show_col_types = FALSE
  ) %>%
    filter(
      project == "PRJEB82425",
      source_column == "infection_group"
    )

  cnt <- rr %>%
    group_by(raw_source_value) %>%
    summarise(
      observed_samples = n(),
      observed_patients = n_distinct(
        patient_id_candidate[
          !is.na(patient_id_candidate)
        ]
      ),
      .groups = "drop"
    ) %>%
    rename(raw_value = raw_source_value)

  prjeb_validation <- prjeb_validation %>%
    left_join(cnt, by = "raw_value") %>%
    mutate(
      patient_count_matches_publication =
        observed_patients ==
        expected_patient_n_from_publication
    )
}

write_csv(
  prjeb_validation,
  file.path(OUT, "10_PRJEB82425_PHENOTYPE_PUBLICATION_VALIDATION.csv")
)

# ------------------------------------------------------------
# 6. Final branch decision
# ------------------------------------------------------------

formal_source_projects <- project_summary %>%
  filter(
    final_feasibility ==
      "FORMAL_INFECTION_SOURCE_ANALYSIS_FEASIBLE"
  ) %>%
  pull(project)

supportive_source_projects <- project_summary %>%
  filter(
    final_feasibility ==
      "SUPPORTIVE_INFECTION_SOURCE_ANALYSIS_FEASIBLE"
  ) %>%
  pull(project)

phenotype_projects <- project_summary %>%
  filter(
    final_feasibility ==
      "INFECTION_PHENOTYPE_ONLY_NOT_SOURCE"
  ) %>%
  pull(project)

branch_status <- case_when(
  length(formal_source_projects) >= 2 ~
    "MAJOR_MULTI_COHORT_INFECTION_SOURCE_BRANCH_FEASIBLE",

  length(formal_source_projects) == 1 &&
    length(supportive_source_projects) >= 1 ~
    "LIMITED_INFECTION_SOURCE_BRANCH_FEASIBLE",

  length(formal_source_projects) == 1 ~
    "SINGLE_COHORT_INFECTION_SOURCE_ONLY",

  length(supportive_source_projects) >= 1 ~
    "SUPPORTIVE_INFECTION_SOURCE_ONLY",

  TRUE ~
    "INFECTION_SOURCE_BRANCH_NOT_FEASIBLE_WITH_CURRENT_METADATA"
)

decision <- tibble(
  branch_status = branch_status,
  formal_true_source_projects = paste(
    formal_source_projects,
    collapse = ";"
  ),
  supportive_true_source_projects = paste(
    supportive_source_projects,
    collapse = ";"
  ),
  infection_phenotype_only_projects = paste(
    phenotype_projects,
    collapse = ";"
  ),
  recommended_next_action = case_when(
    branch_status ==
      "INFECTION_SOURCE_BRANCH_NOT_FEASIBLE_WITH_CURRENT_METADATA" ~
      "Close infection-source as a major V2 branch; optionally retain PRJEB82425 VAP/other-infection/control as supportive infection-phenotype context, then proceed to final statistical synthesis.",

    TRUE ~
      "Manually verify true infection-source field values against original publications before source-stratified modeling."
  )
)

write_csv(
  decision,
  file.path(OUT, "11_FINAL_INFECTION_SOURCE_BRANCH_DECISION.csv")
)

# ------------------------------------------------------------
# 7. Interpretation
# ------------------------------------------------------------

txt <- c(
  "STEP94A2 EXCEL-INCLUSIVE INFECTION-SOURCE AUDIT",
  "",
  paste0(
    "Final branch status: ",
    branch_status
  ),
  "",
  paste0(
    "Formal true infection-source projects: ",
    ifelse(
      length(formal_source_projects) == 0,
      "NONE",
      paste(formal_source_projects, collapse = "; ")
    )
  ),
  paste0(
    "Supportive true infection-source projects: ",
    ifelse(
      length(supportive_source_projects) == 0,
      "NONE",
      paste(supportive_source_projects, collapse = "; ")
    )
  ),
  paste0(
    "Infection-phenotype-only projects: ",
    ifelse(
      length(phenotype_projects) == 0,
      "NONE",
      paste(phenotype_projects, collapse = "; ")
    )
  ),
  "",
  "Critical semantic distinction:",
  "A field such as Sepsis/Healthy/Non-sepsis ICU is a disease-group field, not infection source.",
  "A field such as Pneumonia/Other infection/Control is an infection-phenotype field, not a detailed source taxonomy.",
  "A true infection-source field should distinguish clinical foci such as respiratory, urinary, abdominal/GI, bloodstream, skin/soft tissue, CNS, device/line, etc.",
  "",
  "PRJEB82425:",
  "The published cohort consists of 17 VAP patients, 12 patients with other infections, and 9 uninfected controls. This supports infection-phenotype analysis, but 'other infection' is not sufficiently granular to establish multi-source heterogeneity.",
  "",
  "No source-stratified modeling should be performed unless file 05 identifies a genuine clinical infection-source field after Excel-inclusive scanning."
)

writeLines(
  txt,
  file.path(OUT, "12_STEP94A2_INTERPRETATION.txt")
)

# ------------------------------------------------------------
# 8. QC
# ------------------------------------------------------------

qc <- tibble(
  files_inventoried = nrow(inventory),
  sheets_or_text_tables_read = sum(sheet_audit$readable, na.rm = TRUE),
  candidate_columns = nrow(column_audit),
  true_source_candidate_columns = nrow(true_source),
  infection_phenotype_columns = nrow(phenotype_fields),
  case_control_columns = nrow(case_control_fields),
  manual_review_columns = nrow(manual_review),
  formal_true_source_projects = length(formal_source_projects),
  supportive_true_source_projects = length(supportive_source_projects),
  branch_status = branch_status
)

write_csv(
  qc,
  file.path(OUT, "13_STEP94A2_QC_SUMMARY.csv")
)

writeLines(
  c(
    paste0("Completed: ", Sys.time()),
    paste0("Branch status: ", branch_status),
    "STEP94A2 COMPLETE"
  ),
  file.path(OUT, "_STEP94A2_COMPLETE.ok")
)

cat("STEP94A2 COMPLETE\n")
