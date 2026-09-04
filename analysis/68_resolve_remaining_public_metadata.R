# ============================================================
# Sepsis V2 - Step 68
# Resolve remaining public metadata gaps
# R 4.4.0 / Windows
#
# Goals:
# 1. PRJNA884103:
#    merge Step67 second-pass Run IDs back into READY map/QC.
# 2. PRJNA595346:
#    inspect actual local NCBI BioSample column names and build
#    mapping ONLY from explicit patient/time fields; no guessing.
#    Also download Nature Source Data for clinical-field inventory.
# 3. PRJEB67798:
#    profile all parsed supplement sheets + human-only repository
#    metadata; identify whether sample/patient/time linkage exists.
# 4. PRJEB82425:
#    reconstruct DOCX supplementary tables from officer summary
#    as far as possible and export table-cell structure.
# 5. Produce final current manual-collection decision table.
#
# Raw sequencing files and source metadata are not modified.
# ============================================================

options(stringsAsFactors = FALSE)

pkgs <- c(
  "readr","readxl","dplyr","tidyr","stringr","purrr",
  "tibble","httr2","officer"
)

missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing) > 0) {
  message("Installing missing packages: ", paste(missing, collapse = ", "))
  install.packages(missing, repos = "https://cloud.r-project.org")
}

suppressPackageStartupMessages({
  library(readr)
  library(readxl)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(purrr)
  library(tibble)
  library(httr2)
  library(officer)
})

DATA_ROOT <- "E:/sepsis_project/data"
STEP66_ROOT <- "E:/sepsis_project/results/V2_06_sample_patient_time_map"
STEP67_ROOT <- "E:/sepsis_project/results/V2_07_mapping_repair"
OUT_ROOT <- "E:/sepsis_project/results/V2_08_remaining_metadata_resolution"

dir.create(OUT_ROOT, recursive = TRUE, showWarnings = FALSE)

UA <- "SepsisV2-R44-Step68/1.0"


# ------------------------------------------------------------
# Helpers
# ------------------------------------------------------------
clean_chr <- function(x) {
  x <- trimws(as.character(x))
  x[
    is.na(x) |
    x == "" |
    tolower(x) %in% c("na","nan","n/a","null","none")
  ] <- NA_character_
  x
}

norm_id <- function(x) {
  x <- clean_chr(x)
  ifelse(
    is.na(x),
    NA_character_,
    toupper(gsub("[^A-Za-z0-9._-]", "", x))
  )
}

find_latest <- function(folder, pattern) {
  if (!dir.exists(folder)) return(NA_character_)
  x <- list.files(folder, pattern = pattern, full.names = TRUE)
  if (length(x) == 0) return(NA_character_)
  x[which.max(file.info(x)$mtime)]
}

find_recursive <- function(folder, pattern) {
  if (!dir.exists(folder)) return(character())
  x <- list.files(folder, recursive = TRUE, full.names = TRUE)
  x[str_detect(basename(x), regex(pattern, ignore_case = TRUE))]
}

safe_csv <- function(path) {
  if (length(path) == 0 || is.na(path) || !file.exists(path)) return(NULL)
  tryCatch(
    suppressMessages(
      read_csv(
        path,
        show_col_types = FALSE,
        progress = FALSE,
        name_repair = "unique"
      )
    ),
    error = function(e) NULL
  )
}

safe_tsv <- function(path) {
  if (length(path) == 0 || is.na(path) || !file.exists(path)) return(NULL)
  tryCatch(
    suppressMessages(
      read_tsv(
        path,
        show_col_types = FALSE,
        progress = FALSE,
        name_repair = "unique"
      )
    ),
    error = function(e) NULL
  )
}

sample_values <- function(x, n = 8) {
  x <- clean_chr(x)
  x <- unique(x[!is.na(x)])
  paste(head(x, n), collapse = " | ")
}

profile_columns <- function(df, project, source_file) {
  if (is.null(df) || ncol(df) == 0) return(tibble())

  map_dfr(names(df), function(nm) {
    x <- df[[nm]]

    candidate_role <- case_when(
      str_detect(nm, regex("patient|subject|participant|host.*subject|host.*id", ignore_case=TRUE)) ~ "PATIENT_CANDIDATE",
      str_detect(nm, regex("sample|biosample|specimen|alias", ignore_case=TRUE)) ~ "SAMPLE_CANDIDATE",
      str_detect(nm, regex("run|srr|err|crr", ignore_case=TRUE)) ~ "RUN_CANDIDATE",
      str_detect(nm, regex("time|day|visit|date|follow|interval", ignore_case=TRUE)) ~ "TIME_CANDIDATE",
      str_detect(nm, regex("sepsis|infection|bacter|pneum|vap|diagnos|group", ignore_case=TRUE)) ~ "INFECTION_CANDIDATE",
      str_detect(nm, regex("outcome|mortality|death|surviv|vfd|length.*stay", ignore_case=TRUE)) ~ "OUTCOME_CANDIDATE",
      str_detect(nm, regex("antibiotic|antimicrobial|abx", ignore_case=TRUE)) ~ "ANTIBIOTIC_CANDIDATE",
      str_detect(nm, regex("sofa|apache|severity|shock|vasopressor", ignore_case=TRUE)) ~ "SEVERITY_CANDIDATE",
      TRUE ~ ""
    )

    tibble(
      Project = project,
      Source_File = source_file,
      Column = nm,
      Candidate_Role = candidate_role,
      NonMissing = sum(!is.na(clean_chr(x))),
      Unique_N = n_distinct(clean_chr(x), na.rm = TRUE),
      Example_Values = sample_values(x)
    )
  })
}

read_ena <- function(project) {
  base <- file.path(DATA_ROOT, project, "00_metadata")

  p <- find_recursive(base, "^01_ENA_read_run_metadata\\.tsv$")
  if (length(p) == 0) {
    p <- find_recursive(base, paste0("^", project, "_ENA_run_manifest\\.tsv$"))
  }
  if (length(p) == 0) return(NULL)

  safe_tsv(p[1])
}

read_biosample <- function(project) {
  base <- file.path(DATA_ROOT, project, "00_metadata")
  p <- find_recursive(base, "^06_NCBI_BioSample_attributes\\.csv$")
  if (length(p) == 0) return(NULL)
  safe_csv(p[1])
}

pick_first <- function(df, patterns) {
  if (is.null(df)) return(NA_character_)
  for (pat in patterns) {
    hit <- names(df)[str_detect(names(df), regex(pat, ignore_case=TRUE))]
    if (length(hit) > 0) return(hit[1])
  }
  NA_character_
}

download_file <- function(url, dest) {
  dir.create(dirname(dest), recursive = TRUE, showWarnings = FALSE)

  if (file.exists(dest) && file.info(dest)$size > 0) {
    return("EXISTS")
  }

  req <- request(url) |>
    req_user_agent(UA) |>
    req_headers(Accept = "*/*") |>
    req_timeout(300)

  resp <- req_perform(req, path = dest)
  paste0("HTTP_", resp_status(resp))
}


# ============================================================
# A. PRJNA884103 - merge second-pass Run repair into READY
# ============================================================
ready66_path <- find_latest(
  STEP66_ROOT,
  "^V2_sample_patient_time_map_READY_.*\\.csv$"
)

ready66 <- safe_csv(ready66_path)

second884_path <- find_latest(
  STEP67_ROOT,
  "^PRJNA884103_secondpass_run_mapping_.*\\.csv$"
)

second884 <- safe_csv(second884_path)

ready68 <- ready66

if (!is.null(ready66) && !is.null(second884)) {

  base884 <- ready66 |>
    filter(Project == "PRJNA884103")

  other_ready <- ready66 |>
    filter(Project != "PRJNA884103")

  repair_lookup <- second884 |>
    select(
      Sample_ID,
      Run_ID_repair = Run_ID,
      Mapping_Method_repair = any_of("Mapping_Method")
    ) |>
    filter(!is.na(Run_ID_repair)) |>
    distinct(Sample_ID, .keep_all = TRUE)

  base884_repaired <- base884 |>
    left_join(repair_lookup, by = "Sample_ID") |>
    mutate(
      Run_ID = coalesce(Run_ID, Run_ID_repair),
      Mapping_Method = ifelse(
        !is.na(Run_ID_repair),
        coalesce(Mapping_Method_repair, Mapping_Method),
        Mapping_Method
      )
    ) |>
    select(-any_of(c("Run_ID_repair","Mapping_Method_repair")))

  ready68 <- bind_rows(
    other_ready,
    base884_repaired
  )
}


# ============================================================
# B. PRJNA595346 - inspect actual BioSample fields and map
#    only explicit patient/time information.
# ============================================================
project <- "PRJNA595346"

bio595 <- read_biosample(project)
ena595 <- read_ena(project)

profile595 <- profile_columns(
  bio595,
  project,
  "06_NCBI_BioSample_attributes.csv"
)

write_excel_csv(
  profile595,
  file.path(OUT_ROOT, "PRJNA595346_BioSample_column_profile.csv"),
  na = ""
)

map595 <- tibble()
status595 <- "NOT_MAPPED"

if (!is.null(bio595) && !is.null(ena595)) {

  bio_acc_col <- pick_first(
    bio595,
    c("^BioSample_accession$", "biosample.*accession")
  )

  patient_col <- pick_first(
    bio595,
    c(
      "^SubjectID$",
      "^Subject$",
      "^Patient$",
      "patient.?id",
      "subject.?id",
      "host.?subject",
      "participant.?id"
    )
  )

  time_col <- pick_first(
    bio595,
    c(
      "^Day$",
      "^Study_Day$",
      "studyday",
      "timepoint",
      "time.?point",
      "followup",
      "collection.?day",
      "collection.?date"
    )
  )

  sample_type_col <- pick_first(
    bio595,
    c(
      "^sample_type$",
      "^SampleType$",
      "body.?site",
      "isolation.?source",
      "^tissue$"
    )
  )

  # ENA fields
  ena_run_col <- pick_first(ena595, c("^run_accession$", "^Run$"))
  ena_sample_col <- pick_first(ena595, c("^sample_accession$"))
  ena_secondary_col <- pick_first(ena595, c("^secondary_sample_accession$"))
  ena_alias_col <- pick_first(ena595, c("^sample_alias$"))

  if (
    !is.na(bio_acc_col) &&
    !is.na(patient_col) &&
    !is.na(time_col) &&
    !is.na(ena_run_col)
  ) {

    bio_std <- bio595 |>
      transmute(
        BioSample = clean_chr(.data[[bio_acc_col]]),
        Patient_ID = clean_chr(.data[[patient_col]]),
        Time_Raw = clean_chr(.data[[time_col]]),
        Body_Site = if (!is.na(sample_type_col)) clean_chr(.data[[sample_type_col]]) else NA_character_
      )

    ena_primary <- tibble(
      BioSample = if (!is.na(ena_sample_col)) clean_chr(ena595[[ena_sample_col]]) else NA_character_,
      Run_ID = clean_chr(ena595[[ena_run_col]]),
      Sample_Alias = if (!is.na(ena_alias_col)) clean_chr(ena595[[ena_alias_col]]) else NA_character_,
      Match_Source = "ENA_sample_accession"
    )

    ena_secondary <- tibble(
      BioSample = if (!is.na(ena_secondary_col)) clean_chr(ena595[[ena_secondary_col]]) else NA_character_,
      Run_ID = clean_chr(ena595[[ena_run_col]]),
      Sample_Alias = if (!is.na(ena_alias_col)) clean_chr(ena595[[ena_alias_col]]) else NA_character_,
      Match_Source = "ENA_secondary_sample_accession"
    )

    ena_lookup <- bind_rows(
      ena_primary,
      ena_secondary
    ) |>
      filter(!is.na(BioSample), !is.na(Run_ID)) |>
      distinct()

    map595 <- bio_std |>
      left_join(
        ena_lookup,
        by = "BioSample"
      ) |>
      filter(
        !is.na(Patient_ID),
        !is.na(Time_Raw),
        !is.na(Run_ID)
      ) |>
      transmute(
        Project = project,
        Patient_ID,
        Sample_ID = coalesce(Sample_Alias, BioSample),
        BioSample,
        Run_ID,
        Time_Raw,
        Time_Day = suppressWarnings(as.numeric(Time_Raw)),
        Collection_Date = if (
          str_detect(time_col, regex("date", ignore_case=TRUE))
        ) Time_Raw else NA_character_,
        Body_Site,
        Source_File = "06_NCBI_BioSample_attributes.csv + ENA",
        Mapping_Method = Match_Source,
        Patient_Source_Column = patient_col,
        Time_Source_Column = time_col
      ) |>
      distinct()

    if (nrow(map595) > 0) {
      status595 <- "EXPLICIT_BIOSAMPLE_PATIENT_TIME_MAPPING_RECOVERED"
    } else {
      status595 <- "EXPLICIT_COLUMNS_FOUND_BUT_NO_RUN_JOIN"
    }
  } else {

    missing_bits <- c(
      if (is.na(bio_acc_col)) "biosample_accession" else NULL,
      if (is.na(patient_col)) "patient_id" else NULL,
      if (is.na(time_col)) "time" else NULL,
      if (is.na(ena_run_col)) "run_id" else NULL
    )

    status595 <- paste0(
      "EXPLICIT_MAPPING_FIELDS_MISSING:",
      paste(missing_bits, collapse = ";")
    )
  }
}

if (nrow(map595) > 0) {

  write_excel_csv(
    map595,
    file.path(OUT_ROOT, "PRJNA595346_explicit_sample_patient_time_map.csv"),
    na = ""
  )

  # Replace/add PRJNA595346 in READY.
  if (!is.null(ready68)) {
    ready68 <- ready68 |>
      filter(Project != "PRJNA595346") |>
      bind_rows(map595)
  }
}


# ------------------------------------------------------------
# PRJNA595346 Nature Source Data:
# useful for clinical variable inventory, but does not establish
# sample->patient linkage unless identifiers are present.
# ------------------------------------------------------------
nature595_url <- paste0(
  "https://media.springernature.com/original/",
  "springer-static/esm/art%3A10.1038%2Fs41467-024-48819-8/",
  "MediaObjects/41467_2024_48819_MOESM4_ESM.xlsx"
)

nature595_file <- file.path(
  DATA_ROOT,
  "PRJNA595346",
  "00_metadata",
  "auto_fetched",
  "12_Nature_source_data",
  "PRJNA595346_Nature_Source_Data.xlsx"
)

nature595_status <- tryCatch(
  download_file(nature595_url, nature595_file),
  error = function(e) paste0("FAILED: ", conditionMessage(e))
)

source595_profile <- tibble()

if (file.exists(nature595_file)) {

  sheets <- tryCatch(
    excel_sheets(nature595_file),
    error = function(e) character()
  )

  source595_profile <- map_dfr(sheets, function(sh) {

    x <- tryCatch(
      read_excel(
        nature595_file,
        sheet = sh,
        n_max = 50,
        .name_repair = "unique"
      ),
      error = function(e) NULL
    )

    if (is.null(x)) return(tibble())

    tibble(
      Sheet = sh,
      Rows_Previewed = nrow(x),
      Columns = ncol(x),
      ColumnNames = paste(names(x), collapse = ";"),
      Has_Patient_ID_Column = any(str_detect(
        names(x),
        regex("patient|subject|participant", ignore_case=TRUE)
      )),
      Has_Sample_ID_Column = any(str_detect(
        names(x),
        regex("sample.*id|biosample|run_accession", ignore_case=TRUE)
      )),
      Has_Time_Column = any(str_detect(
        names(x),
        regex("time|day|follow|interval", ignore_case=TRUE)
      )),
      Has_Outcome_Column = any(str_detect(
        names(x),
        regex("mortality|survival|outcome|vfd", ignore_case=TRUE)
      )),
      Has_Severity_Column = any(str_detect(
        names(x),
        regex("sofa|apache|severity", ignore_case=TRUE)
      )),
      Has_Antibiotic_Column = any(str_detect(
        names(x),
        regex("antibiotic|abx", ignore_case=TRUE)
      ))
    )
  })

  write_excel_csv(
    source595_profile,
    file.path(OUT_ROOT, "PRJNA595346_Nature_Source_Data_profile.csv"),
    na = ""
  )
}


# ============================================================
# C. PRJEB67798 supplement profiling
# ============================================================
project <- "PRJEB67798"

supp677_dir <- file.path(
  DATA_ROOT,
  project,
  "00_metadata",
  "auto_fetched",
  "07_known_public_supplement"
)

supp677_files <- find_recursive(
  supp677_dir,
  "\\.(xlsx|xls|csv)$"
)

supp677_profile <- list()

for (p in supp677_files) {

  ext <- tolower(tools::file_ext(p))

  if (ext %in% c("xlsx","xls")) {

    sheets <- tryCatch(excel_sheets(p), error=function(e) character())

    for (sh in sheets) {

      x <- tryCatch(
        read_excel(p, sheet=sh, .name_repair="unique"),
        error=function(e) NULL
      )

      if (is.null(x)) next

      supp677_profile[[length(supp677_profile)+1]] <- profile_columns(
        x,
        project,
        paste0(basename(p), "::", sh)
      )
    }

  } else {

    x <- safe_csv(p)

    if (!is.null(x)) {
      supp677_profile[[length(supp677_profile)+1]] <- profile_columns(
        x,
        project,
        basename(p)
      )
    }
  }
}

supp677_profile_df <- if (length(supp677_profile)>0) {
  bind_rows(supp677_profile)
} else {
  tibble()
}

write_excel_csv(
  supp677_profile_df,
  file.path(OUT_ROOT, "PRJEB67798_supplement_column_profile.csv"),
  na = ""
)

# Use previously created human-only metadata if present.
human677_path <- find_recursive(
  STEP67_ROOT,
  "^PRJEB67798_HUMAN_ONLY_repository_metadata\\.csv$"
)

human677 <- if (length(human677_path)>0) safe_csv(human677_path[1]) else NULL

if (!is.null(human677)) {
  write_excel_csv(
    profile_columns(
      human677,
      project,
      basename(human677_path[1])
    ),
    file.path(OUT_ROOT, "PRJEB67798_human_repository_column_profile.csv"),
    na = ""
  )
}


# ============================================================
# D. PRJEB82425 DOCX supplement - reconstruct tables
# ============================================================
project <- "PRJEB82425"

docx824 <- find_recursive(
  file.path(DATA_ROOT, project, "00_metadata"),
  "PRJEB82425.*Supplementary.*\\.docx$"
)

docx824_tables <- tibble()
docx824_profile <- tibble()

if (length(docx824) > 0) {

  sum824 <- tryCatch(
    officer::docx_summary(
      officer::read_docx(docx824[1])
    ) |>
      as_tibble(),
    error = function(e) tibble()
  )

  if (nrow(sum824) > 0) {

    # Preserve full summary.
    write_excel_csv(
      sum824,
      file.path(OUT_ROOT, "PRJEB82425_DOCX_full_summary.csv"),
      na = ""
    )

    # Table-cell extraction when officer provides table/row/cell indices.
    needed <- c("table_index","row_id","cell_id","text")

    if (all(needed %in% names(sum824))) {

      docx824_tables <- sum824 |>
        filter(
          !is.na(table_index),
          !is.na(row_id),
          !is.na(cell_id)
        ) |>
        select(
          table_index,
          row_id,
          cell_id,
          text
        ) |>
        arrange(
          table_index,
          row_id,
          cell_id
        )

      write_excel_csv(
        docx824_tables,
        file.path(OUT_ROOT, "PRJEB82425_DOCX_table_cells.csv"),
        na = ""
      )
    }

    # Profile the textual content for useful clinical terms.
    docx824_profile <- tibble(
      File = basename(docx824[1]),
      Total_Content_Rows = nrow(sum824),
      Patient_Term_Hits = sum(str_detect(
        clean_chr(sum824$text),
        regex("patient|subject", ignore_case=TRUE)
      ), na.rm=TRUE),
      Sample_Term_Hits = sum(str_detect(
        clean_chr(sum824$text),
        regex("sample|rectal|swab|trache", ignore_case=TRUE)
      ), na.rm=TRUE),
      Time_Term_Hits = sum(str_detect(
        clean_chr(sum824$text),
        regex("day|time|sampling|baseline|follow", ignore_case=TRUE)
      ), na.rm=TRUE),
      Infection_Term_Hits = sum(str_detect(
        clean_chr(sum824$text),
        regex("infection|sepsis|vap|pneumonia", ignore_case=TRUE)
      ), na.rm=TRUE),
      Antibiotic_Term_Hits = sum(str_detect(
        clean_chr(sum824$text),
        regex("antibiotic|antimicrobial", ignore_case=TRUE)
      ), na.rm=TRUE)
    )

    write_excel_csv(
      docx824_profile,
      file.path(OUT_ROOT, "PRJEB82425_DOCX_content_profile.csv"),
      na = ""
    )
  }
}


# ============================================================
# E. Updated READY QC after Step68
# ============================================================
ready_qc68 <- tibble()

if (!is.null(ready68) && nrow(ready68) > 0) {

  ready_qc68 <- ready68 |>
    group_by(Project) |>
    summarise(
      Rows = n(),
      Unique_Patients = n_distinct(Patient_ID[!is.na(Patient_ID)]),
      Unique_Samples = n_distinct(Sample_ID[!is.na(Sample_ID)]),
      Unique_Runs = n_distinct(Run_ID[!is.na(Run_ID)]),
      Missing_Patient = sum(is.na(Patient_ID)),
      Missing_Run = sum(is.na(Run_ID)),
      Missing_Time = sum(is.na(Time_Raw)),
      .groups = "drop"
    )

  stamp <- format(Sys.time(), "%Y%m%d_%H%M%S")

  write_excel_csv(
    ready68,
    file.path(
      OUT_ROOT,
      paste0("V2_READY_map_after68_", stamp, ".csv")
    ),
    na = ""
  )

  write_excel_csv(
    ready_qc68,
    file.path(
      OUT_ROOT,
      paste0("V2_READY_QC_after68_", stamp, ".csv")
    ),
    na = ""
  )
}


# ============================================================
# F. Current manual collection decision
# ============================================================
# Assess whether existing files contain explicit structured candidates.
has_role <- function(profile, role) {
  if (is.null(profile) || nrow(profile)==0) return(FALSE)
  any(profile$Candidate_Role == role, na.rm=TRUE)
}

manual_decision <- tibble(
  Project = c(
    "PRJNA516701",
    "PRJNA578267",
    "PRJNA884103",
    "PRJNA595346",
    "PRJEB37289",
    "PRJEB67798",
    "PRJEB82425"
  ),

  Current_Status = c(
    "READY",
    "READY",
    if (
      !is.null(ready68) &&
      all(!is.na(ready68$Run_ID[ready68$Project=="PRJNA884103"]))
    ) "READY_SECOND_PASS_MERGED" else "READY_CHECK_RUNS",
    status595,
    "STATIC_BACKGROUND_ONLY",
    if (
      has_role(supp677_profile_df, "PATIENT_CANDIDATE") &&
      has_role(supp677_profile_df, "TIME_CANDIDATE")
    ) "SUPPLEMENT_HAS_PATIENT_TIME_CANDIDATES" else "PATIENT_LINKAGE_NOT_CONFIRMED",
    if (
      nrow(docx824_tables) > 0
    ) "SUPPLEMENT_TABLES_EXTRACTED" else "SUPPLEMENT_TEXT_ONLY"
  ),

  Manual_Collection_Required_Now = c(
    "NO",
    "NO",
    "NO",
    if (nrow(map595)>0) "NO" else "NOT_YET",
    "NO",
    "NOT_YET",
    "NOT_YET"
  ),

  Next_Action = c(
    "Proceed",
    "Proceed",
    "Proceed with repaired Run mapping",
    if (nrow(map595)>0) {
      "Proceed with explicit BioSample mapping"
    } else {
      "Review exported BioSample column profile. If no explicit patient/time field exists, PRJNA595346 cannot be used for patient-level longitudinal analysis without author data."
    },
    "Use only as static ICU/sepsis background validation",
    "Review supplement column profile and human-only metadata. Contact authors only if no sample-patient linkage is present.",
    "Review extracted DOCX table cells. Manual reading is needed only if sample-patient linkage is encoded in figures/text rather than tables."
  )
)

write_excel_csv(
  manual_decision,
  file.path(OUT_ROOT, "V2_manual_collection_decision_after68.csv"),
  na = ""
)


# ============================================================
# Console
# ============================================================
cat("\n============================================================\n")
cat("SEPSIS V2 - STEP 68 COMPLETE\n")
cat("R version:", R.version.string, "\n")
cat("============================================================\n\n")

cat("PRJNA595346 mapping status:\n", status595, "\n\n")
cat("PRJNA595346 Nature source-data download:\n", nature595_status, "\n\n")

cat("Updated READY QC:\n")
print(ready_qc68, n=Inf, width=Inf)

cat("\nManual collection decision:\n")
print(manual_decision, n=Inf, width=Inf)

cat("\nOutputs:\n", OUT_ROOT, "\n")
cat("============================================================\n")
