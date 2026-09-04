# ============================================================
# Sepsis V2 - Step 67
# Repair/adjudicate sample-patient-time mapping
# Designed for R 4.4.0 on Windows
#
# Main tasks:
# 1) Recover PRJNA595346 release/Zenodo de-identified metadata.
# 2) Repair the 13 PRJNA884103 samples without SRA Run IDs.
# 3) Download and parse PRJEB67798 official XLSX supplement.
# 4) Download and parse PRJEB82425 official DOCX supplement.
# 5) Download PRJEB37289 MOFA processed metadata and reclassify
#    it as STATIC_ICU_BACKGROUND rather than core longitudinal.
#
# No source/raw sequencing data are modified.
# ============================================================

options(stringsAsFactors = FALSE)

pkgs <- c(
  "readr", "readxl", "dplyr", "tidyr", "stringr",
  "purrr", "tibble", "httr2", "jsonlite", "data.table",
  "officer"
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
  library(jsonlite)
  library(data.table)
  library(officer)
})

DATA_ROOT <- "E:/sepsis_project/data"
RESULT_ROOT <- "E:/sepsis_project/results/V2_07_mapping_repair"
dir.create(RESULT_ROOT, recursive = TRUE, showWarnings = FALSE)

UA <- "SepsisV2-R44-MetadataRepair/1.0"

# ------------------------------------------------------------
# Helpers
# ------------------------------------------------------------
clean_chr <- function(x) {
  x <- trimws(as.character(x))
  x[x == "" | is.na(x) | tolower(x) %in% c("na","nan","null","none","n/a")] <- NA_character_
  x
}

norm_id <- function(x) {
  x <- clean_chr(x)
  ifelse(
    is.na(x),
    NA_character_,
    toupper(gsub("[[:space:]]+", "", x))
  )
}

find_latest <- function(folder, pattern) {
  x <- list.files(folder, pattern = pattern, full.names = TRUE)
  if (length(x) == 0) return(NA_character_)
  x[which.max(file.info(x)$mtime)]
}

find_recursive <- function(folder, pattern) {
  if (!dir.exists(folder)) return(character())
  x <- list.files(folder, recursive = TRUE, full.names = TRUE)
  x[str_detect(basename(x), regex(pattern, ignore_case = TRUE))]
}

safe_read_csv <- function(path) {
  if (length(path) == 0 || is.na(path) || !file.exists(path)) return(NULL)
  tryCatch(
    suppressMessages(read_csv(path, show_col_types = FALSE, progress = FALSE, name_repair = "unique")),
    error = function(e) NULL
  )
}

safe_read_tsv <- function(path) {
  if (length(path) == 0 || is.na(path) || !file.exists(path)) return(NULL)
  tryCatch(
    suppressMessages(read_tsv(path, show_col_types = FALSE, progress = FALSE, name_repair = "unique")),
    error = function(e) NULL
  )
}

pick_col <- function(df, exact = character(), regex_pat = NULL) {
  if (is.null(df)) return(NA_character_)
  hit <- exact[exact %in% names(df)]
  if (length(hit) > 0) return(hit[1])
  if (!is.null(regex_pat)) {
    hit <- names(df)[str_detect(names(df), regex(regex_pat, ignore_case = TRUE))]
    if (length(hit) > 0) return(hit[1])
  }
  NA_character_
}

col_or_na <- function(df, col) {
  if (is.null(df) || is.na(col) || !col %in% names(df)) return(rep(NA_character_, nrow(df)))
  clean_chr(df[[col]])
}

download_binary <- function(url, dest, force = FALSE) {
  dir.create(dirname(dest), recursive = TRUE, showWarnings = FALSE)
  if (file.exists(dest) && file.info(dest)$size > 0 && !force) return("EXISTS")

  req <- request(url) |>
    req_user_agent(UA) |>
    req_headers(Accept = "*/*") |>
    req_timeout(300)

  resp <- req_perform(req, path = dest)
  paste0("HTTP_", resp_status(resp))
}

extract_archive <- function(path, outdir) {
  if (!file.exists(path)) return(FALSE)
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  low <- tolower(path)

  ok <- tryCatch({
    if (str_detect(low, "\\.zip$")) {
      unzip(path, exdir = outdir)
      TRUE
    } else if (str_detect(low, "\\.(tar\\.gz|tgz|tar)$")) {
      untar(path, exdir = outdir)
      TRUE
    } else {
      FALSE
    }
  }, error = function(e) FALSE)

  ok
}

read_ena <- function(project) {
  base <- file.path(DATA_ROOT, project, "00_metadata")
  p <- find_recursive(base, "^01_ENA_read_run_metadata\\.tsv$")
  if (length(p) == 0) {
    p <- find_recursive(base, paste0("^", project, "_ENA_run_manifest\\.tsv$"))
  }
  if (length(p) == 0) return(NULL)
  safe_read_tsv(p[1])
}

# ------------------------------------------------------------
# A. Load Step 66 outputs
# ------------------------------------------------------------
step66 <- "E:/sepsis_project/results/V2_06_sample_patient_time_map"

ready_path <- find_latest(step66, "^V2_sample_patient_time_map_READY_.*\\.csv$")
unresolved_path <- find_latest(step66, "^V2_patient_mapping_UNRESOLVED_.*\\.csv$")
unmatched_path <- find_latest(step66, "^V2_unmatched_run_QC_.*\\.csv$")

if (is.na(ready_path)) stop("Step 66 READY file not found.")
if (is.na(unresolved_path)) stop("Step 66 unresolved file not found.")

ready66 <- safe_read_csv(ready_path)
unresolved66 <- safe_read_csv(unresolved_path)
unmatched66 <- safe_read_csv(unmatched_path)

# ------------------------------------------------------------
# B. PRJNA595346: fetch release assets + Zenodo fallback
# ------------------------------------------------------------
project <- "PRJNA595346"
base595 <- file.path(DATA_ROOT, project, "00_metadata", "auto_fetched")
release_dir <- file.path(base595, "10_GitHub_release_assets")
zenodo_dir <- file.path(base595, "11_Zenodo_archive")
dir.create(release_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(zenodo_dir, recursive = TRUE, showWarnings = FALSE)

release_log <- list()

# GitHub release API
release_api <- paste0(
  "https://api.github.com/repos/MicrobiomeALIR/",
  "MultiCompartmentMicrobiome/releases/tags/MulticompMicrobiomev1"
)

try({
  resp <- request(release_api) |>
    req_user_agent(UA) |>
    req_headers(Accept = "application/vnd.github+json") |>
    req_timeout(60) |>
    req_perform()

  obj <- resp_body_json(resp, simplifyVector = TRUE)

  if (!is.null(obj$assets) && length(obj$assets) > 0) {
    assets <- as.data.frame(obj$assets)

    for (i in seq_len(nrow(assets))) {
      nm <- assets$name[i]
      url <- assets$browser_download_url[i]
      dest <- file.path(release_dir, nm)

      st <- tryCatch(
        download_binary(url, dest),
        error = function(e) paste0("FAILED: ", conditionMessage(e))
      )

      release_log[[length(release_log) + 1]] <- tibble(
        Source = "GitHub_release",
        File = nm,
        Status = st,
        LocalPath = dest
      )

      if (file.exists(dest)) {
        extract_archive(dest, file.path(release_dir, paste0("extracted_", i)))
      }
    }
  }
}, silent = TRUE)

# Zenodo fallback / complement
zenodo_api <- "https://zenodo.org/api/records/11109543"

try({
  resp <- request(zenodo_api) |>
    req_user_agent(UA) |>
    req_timeout(60) |>
    req_perform()

  obj <- resp_body_json(resp, simplifyVector = TRUE)

  if (!is.null(obj$files) && length(obj$files) > 0) {
    files <- as.data.frame(obj$files)

    for (i in seq_len(nrow(files))) {
      nm <- files$key[i]

      url <- NA_character_
      if ("links.self" %in% names(files)) {
        url <- files[["links.self"]][i]
      }

      if (is.na(url) || is.null(url)) next

      dest <- file.path(zenodo_dir, nm)

      st <- tryCatch(
        download_binary(url, dest),
        error = function(e) paste0("FAILED: ", conditionMessage(e))
      )

      release_log[[length(release_log) + 1]] <- tibble(
        Source = "Zenodo",
        File = nm,
        Status = st,
        LocalPath = dest
      )

      if (file.exists(dest)) {
        extract_archive(dest, file.path(zenodo_dir, paste0("extracted_", i)))
      }
    }
  }
}, silent = TRUE)

release_log_df <- if (length(release_log) > 0) bind_rows(release_log) else tibble()

# Locate actual de-identified metadata
meta595_files <- find_recursive(
  file.path(DATA_ROOT, project, "00_metadata"),
  "^MetaData_all\\.csv$"
)

map595 <- tibble()
map595_status <- "METADATA_NOT_FOUND"

if (length(meta595_files) > 0) {

  meta595 <- safe_read_csv(meta595_files[1])
  ena595 <- read_ena(project)

  bio_path <- find_recursive(
    file.path(DATA_ROOT, project, "00_metadata"),
    "^06_NCBI_BioSample_attributes\\.csv$"
  )
  bio595 <- if (length(bio_path) > 0) safe_read_csv(bio_path[1]) else NULL

  if (!is.null(meta595) && !is.null(ena595)) {

    meta_sample_col <- pick_col(meta595, exact = c("sample_id"), regex_pat = "^sample.?id$")
    patient_col <- pick_col(meta595, exact = c("SubjectID","subject_id"), regex_pat = "subject|patient")
    time_col <- pick_col(meta595, exact = c("studydayfollowup"), regex_pat = "followup|timepoint|studyday")
    type_col <- pick_col(meta595, exact = c("SampleType","compartment"), regex_pat = "sample.?type|compartment")

    meta_std <- meta595 |>
      mutate(
        Meta_Sample_ID = norm_id(col_or_na(meta595, meta_sample_col)),
        Patient_ID = col_or_na(meta595, patient_col),
        Time_Raw = col_or_na(meta595, time_col),
        Body_Site_Meta = col_or_na(meta595, type_col)
      )

    # ENA lookup
    ena_alias_col <- pick_col(ena595, exact = c("sample_alias"))
    ena_title_col <- pick_col(ena595, exact = c("sample_title"))
    ena_sample_col <- pick_col(ena595, exact = c("sample_accession"))
    ena_secondary_col <- pick_col(ena595, exact = c("secondary_sample_accession"))
    ena_run_col <- pick_col(ena595, exact = c("run_accession"))

    ena_lookup <- bind_rows(
      tibble(
        Match_ID = norm_id(col_or_na(ena595, ena_alias_col)),
        Run_ID = col_or_na(ena595, ena_run_col),
        ENA_Sample = col_or_na(ena595, ena_sample_col),
        Match_Method = "ENA_sample_alias"
      ),
      tibble(
        Match_ID = norm_id(col_or_na(ena595, ena_title_col)),
        Run_ID = col_or_na(ena595, ena_run_col),
        ENA_Sample = col_or_na(ena595, ena_sample_col),
        Match_Method = "ENA_sample_title"
      )
    ) |>
      filter(!is.na(Match_ID)) |>
      distinct()

    direct <- meta_std |>
      left_join(ena_lookup, by = c("Meta_Sample_ID" = "Match_ID"))

    # If direct matching is weak, try BioSample attribute crosswalk
    via_bio <- tibble()

    if (!is.null(bio595)) {

      bio_acc_col <- pick_col(
        bio595,
        exact = c("BioSample_accession"),
        regex_pat = "biosample.*accession"
      )

      id_cols <- names(bio595)[
        str_detect(
          names(bio595),
          regex("sample|alias|name|subject|id", ignore_case = TRUE)
        )
      ]

      # Avoid columns that are purely internal numeric IDs if possible.
      id_cols <- unique(c(bio_acc_col, id_cols))
      id_cols <- id_cols[!is.na(id_cols)]

      if (length(id_cols) > 0) {

        bio_long <- bio595 |>
          mutate(.row_id = row_number()) |>
          select(.row_id, all_of(id_cols)) |>
          pivot_longer(
            cols = - .row_id,
            names_to = "Bio_Field",
            values_to = "Bio_Value"
          ) |>
          mutate(
            Match_ID = norm_id(Bio_Value)
          ) |>
          filter(
            !is.na(Match_ID),
            nchar(Match_ID) >= 3
          )

        bio_acc <- bio595 |>
          mutate(.row_id = row_number()) |>
          transmute(
            .row_id,
            BioSample_accession = col_or_na(bio595, bio_acc_col)
          )

        bio_long <- bio_long |>
          left_join(bio_acc, by = ".row_id")

        ena_bio_lookup <- bind_rows(
          tibble(
            BioSample_accession = col_or_na(ena595, ena_sample_col),
            Run_ID = col_or_na(ena595, ena_run_col),
            ENA_Sample = col_or_na(ena595, ena_sample_col)
          ),
          tibble(
            BioSample_accession = col_or_na(ena595, ena_secondary_col),
            Run_ID = col_or_na(ena595, ena_run_col),
            ENA_Sample = col_or_na(ena595, ena_sample_col)
          )
        ) |>
          filter(!is.na(BioSample_accession)) |>
          distinct()

        via_bio <- meta_std |>
          inner_join(
            bio_long,
            by = c("Meta_Sample_ID" = "Match_ID")
          ) |>
          left_join(
            ena_bio_lookup,
            by = "BioSample_accession"
          ) |>
          mutate(
            Match_Method = paste0("BioSample_attribute:", Bio_Field)
          )
      }
    }

    direct_good <- direct |>
      filter(!is.na(Run_ID))

    combined595 <- bind_rows(
      direct_good,
      via_bio
    ) |>
      filter(!is.na(Run_ID)) |>
      distinct(Meta_Sample_ID, Run_ID, .keep_all = TRUE)

    if (nrow(combined595) > 0) {

      map595 <- combined595 |>
        transmute(
          Project = project,
          Patient_ID = clean_chr(Patient_ID),
          Sample_ID = clean_chr(ifelse(
            is.na(Meta_Sample_ID),
            ENA_Sample,
            Meta_Sample_ID
          )),
          Run_ID = clean_chr(Run_ID),
          Time_Raw = clean_chr(Time_Raw),
          Time_Day = NA_real_,
          Collection_Date = NA_character_,
          Body_Site = clean_chr(Body_Site_Meta),
          Source_File = basename(meta595_files[1]),
          Mapping_Method = Match_Method
        )

      map595_status <- "MAPPED_FROM_RELEASE_METADATA"
    } else {
      map595_status <- "METADATA_FOUND_BUT_RUN_LINKAGE_NOT_RESOLVED"
    }
  }
}

# ------------------------------------------------------------
# C. PRJNA884103: second-pass Run lookup
# ------------------------------------------------------------
unmatched884 <- tibble()
resolved884 <- tibble()

if (!is.null(unmatched66)) {
  unmatched884 <- unmatched66 |>
    filter(Project == "PRJNA884103")
}

if (nrow(unmatched884) > 0) {

  repo_dir <- file.path(
    DATA_ROOT, "PRJNA884103", "00_metadata",
    "auto_fetched", "09_public_analysis_repository", "extracted"
  )

  metadata_tables <- list.files(
    repo_dir,
    pattern = "\\.(csv|tsv|txt)$",
    recursive = TRUE,
    full.names = TRUE
  )

  candidates <- list()

  for (p in metadata_tables) {

    ext <- tolower(tools::file_ext(p))
    df <- tryCatch(
      {
        if (ext == "csv") {
          suppressMessages(read_csv(p, show_col_types = FALSE, progress = FALSE))
        } else {
          suppressMessages(read_tsv(p, show_col_types = FALSE, progress = FALSE))
        }
      },
      error = function(e) NULL
    )

    if (is.null(df) || nrow(df) == 0) next

    sample_cols <- names(df)[
      str_detect(
        names(df),
        regex("sample|sequence|specimen", ignore_case = TRUE)
      )
    ]

    if (length(sample_cols) == 0) next

    for (sc in sample_cols) {

      vals <- norm_id(df[[sc]])

      hit_idx <- which(
        vals %in% norm_id(unmatched884$Sample_ID)
      )

      if (length(hit_idx) == 0) next

      for (r in hit_idx) {

        rowvals <- clean_chr(unlist(df[r, , drop = FALSE], use.names = FALSE))
        run_hits <- rowvals[
          !is.na(rowvals) &
          str_detect(rowvals, "^(SRR|ERR|CRR)[0-9]+$")
        ]

        if (length(run_hits) > 0) {

          candidates[[length(candidates) + 1]] <- tibble(
            Sample_ID = clean_chr(df[[sc]][r]),
            Run_ID_SecondPass = run_hits[1],
            Secondary_Source = basename(p),
            Secondary_Match_Column = sc
          )
        }
      }
    }
  }

  if (length(candidates) > 0) {

    lookup884 <- bind_rows(candidates) |>
      distinct(Sample_ID, Run_ID_SecondPass, .keep_all = TRUE)

    resolved884 <- unmatched884 |>
      left_join(
        lookup884,
        by = "Sample_ID"
      ) |>
      mutate(
        Run_ID = coalesce(Run_ID, Run_ID_SecondPass),
        Mapping_Method = ifelse(
          !is.na(Run_ID_SecondPass),
          paste0(
            Mapping_Method,
            " + secondary GitHub exact-ID lookup"
          ),
          Mapping_Method
        )
      ) |>
      select(-Run_ID_SecondPass)
  } else {
    resolved884 <- unmatched884
  }
}

# ------------------------------------------------------------
# D. PRJEB37289: download public MOFA metadata
# ------------------------------------------------------------
project <- "PRJEB37289"
mofa_dir <- file.path(
  DATA_ROOT, project, "00_metadata", "auto_fetched",
  "10_MOFA_processed_metadata"
)
dir.create(mofa_dir, recursive = TRUE, showWarnings = FALSE)

mofa_url <- "https://ftp.ebi.ac.uk/pub/databases/mofa/microbiome/metadata.txt.gz"
mofa_dest <- file.path(mofa_dir, "metadata.txt.gz")

mofa_status <- tryCatch(
  download_binary(mofa_url, mofa_dest),
  error = function(e) paste0("FAILED: ", conditionMessage(e))
)

mofa_meta <- NULL

if (file.exists(mofa_dest)) {
  mofa_meta <- tryCatch(
    data.table::fread(mofa_dest, data.table = FALSE),
    error = function(e) NULL
  )
}

if (!is.null(mofa_meta)) {
  write_excel_csv(
    as_tibble(mofa_meta) |>
      mutate(
        Project = project,
        V2_Role = "STATIC_ICU_BACKGROUND",
        Longitudinal_Core = "NO"
      ),
    file.path(
      RESULT_ROOT,
      "PRJEB37289_public_MOFA_metadata.csv"
    ),
    na = ""
  )
}

# ------------------------------------------------------------
# E. PRJEB67798: official supplement + human-only filtering
# ------------------------------------------------------------
project <- "PRJEB67798"
supp677_dir <- file.path(
  DATA_ROOT, project, "00_metadata", "auto_fetched",
  "07_known_public_supplement"
)
dir.create(supp677_dir, recursive = TRUE, showWarnings = FALSE)

supp677_url <- paste0(
  "https://media.springernature.com/original/",
  "springer-static/esm/art%3A10.1038%2Fs41598-023-49034-z/",
  "MediaObjects/41598_2023_49034_MOESM2_ESM.xlsx"
)

supp677 <- file.path(
  supp677_dir,
  "PRJEB67798_Supplementary_Information_2_REDOWNLOAD.xlsx"
)

supp677_status <- tryCatch(
  download_binary(supp677_url, supp677, force = TRUE),
  error = function(e) paste0("FAILED: ", conditionMessage(e))
)

supp677_profiles <- list()

if (file.exists(supp677)) {

  sheets <- tryCatch(
    excel_sheets(supp677),
    error = function(e) character()
  )

  for (sh in sheets) {

    x <- tryCatch(
      read_excel(supp677, sheet = sh, .name_repair = "unique"),
      error = function(e) NULL
    )

    if (is.null(x)) next

    out <- file.path(
      supp677_dir,
      paste0(
        "PRJEB67798_Supplement_sheet_",
        gsub("[^A-Za-z0-9]+", "_", sh),
        ".csv"
      )
    )

    write_excel_csv(x, out, na = "")

    supp677_profiles[[length(supp677_profiles) + 1]] <- tibble(
      Sheet = sh,
      Rows = nrow(x),
      Columns = ncol(x),
      ColumnNames = paste(names(x), collapse = ";")
    )
  }
}

supp677_profiles_df <- if (length(supp677_profiles) > 0) {
  bind_rows(supp677_profiles)
} else {
  tibble()
}

# Human-only ENA subset based strictly on organism/host fields.
ena677 <- read_ena(project)
human677 <- tibble()
human_filter_status <- "NOT_RESOLVED"

if (!is.null(ena677)) {

  human_cols <- names(ena677)[
    str_detect(
      names(ena677),
      regex(
        "scientific_name|host|organism|tax",
        ignore_case = TRUE
      )
    )
  ]

  if (length(human_cols) > 0) {

    human_mask <- rep(FALSE, nrow(ena677))

    for (cc in human_cols) {
      v <- clean_chr(ena677[[cc]])
      human_mask <- human_mask |
        (!is.na(v) & str_detect(
          v,
          regex("Homo sapiens|human", ignore_case = TRUE)
        ))
    }

    if (sum(human_mask) > 0) {
      human677 <- ena677[human_mask, , drop = FALSE]
      human_filter_status <- "FILTERED_BY_ENA_ORGANISM_HOST"
    }
  }
}

# Fallback: use an already-created human-only manifest if present locally.
if (nrow(human677) == 0) {

  all_manifest_candidates <- list.files(
    DATA_ROOT,
    pattern = "PRJEB67798.*human.*only.*\\.(tsv|csv)$",
    recursive = TRUE,
    full.names = TRUE,
    ignore.case = TRUE
  )

  if (length(all_manifest_candidates) > 0) {

    p <- all_manifest_candidates[1]

    if (tolower(tools::file_ext(p)) == "csv") {
      human677 <- safe_read_csv(p)
    } else {
      human677 <- safe_read_tsv(p)
    }

    if (!is.null(human677) && nrow(human677) > 0) {
      human_filter_status <- "USED_EXISTING_HUMAN_ONLY_MANIFEST"
    }
  }
}

if (!is.null(human677) && nrow(human677) > 0) {
  write_excel_csv(
    human677,
    file.path(
      RESULT_ROOT,
      "PRJEB67798_HUMAN_ONLY_repository_metadata.csv"
    ),
    na = ""
  )
}

# ------------------------------------------------------------
# F. PRJEB82425: official DOCX supplement
# ------------------------------------------------------------
project <- "PRJEB82425"
supp824_dir <- file.path(
  DATA_ROOT, project, "00_metadata", "auto_fetched",
  "07_known_public_supplement"
)
dir.create(supp824_dir, recursive = TRUE, showWarnings = FALSE)

supp824_url <- paste0(
  "https://media.springernature.com/original/",
  "springer-static/esm/art%3A10.1186%2Fs12879-025-10825-6/",
  "MediaObjects/12879_2025_10825_MOESM1_ESM.docx"
)

supp824 <- file.path(
  supp824_dir,
  "PRJEB82425_Supplementary_Material_1_REDOWNLOAD.docx"
)

supp824_status <- tryCatch(
  download_binary(supp824_url, supp824, force = TRUE),
  error = function(e) paste0("FAILED: ", conditionMessage(e))
)

docx824_summary <- tibble()

if (file.exists(supp824)) {

  docx824_summary <- tryCatch(
    officer::docx_summary(
      officer::read_docx(supp824)
    ) |>
      as_tibble(),
    error = function(e) tibble()
  )

  if (nrow(docx824_summary) > 0) {
    write_excel_csv(
      docx824_summary,
      file.path(
        RESULT_ROOT,
        "PRJEB82425_supplement_docx_summary.csv"
      ),
      na = ""
    )
  }
}

# ------------------------------------------------------------
# G. Build final repair/adjudication summary
# ------------------------------------------------------------
# Ready projects from Step 66
ready_counts <- ready66 |>
  count(Project, name = "Rows_Ready66")

# Core QC from ready66
ready_qc <- ready66 |>
  group_by(Project) |>
  summarise(
    Rows = n(),
    Unique_Patients = n_distinct(Patient_ID[!is.na(Patient_ID)]),
    Unique_Samples = n_distinct(Sample_ID[!is.na(Sample_ID)]),
    Unique_Runs = n_distinct(Run_ID[!is.na(Run_ID)]),
    Missing_Run = sum(is.na(Run_ID)),
    Missing_Time = sum(is.na(Time_Raw)),
    .groups = "drop"
  )

if (nrow(map595) > 0) {
  ready_qc <- bind_rows(
    ready_qc,
    map595 |>
      summarise(
        Project = "PRJNA595346",
        Rows = n(),
        Unique_Patients = n_distinct(Patient_ID[!is.na(Patient_ID)]),
        Unique_Samples = n_distinct(Sample_ID[!is.na(Sample_ID)]),
        Unique_Runs = n_distinct(Run_ID[!is.na(Run_ID)]),
        Missing_Run = sum(is.na(Run_ID)),
        Missing_Time = sum(is.na(Time_Raw))
      )
  )
}

adjudication <- tibble(
  Project = c(
    "PRJNA516701",
    "PRJNA578267",
    "PRJNA595346",
    "PRJNA884103",
    "PRJEB37289",
    "PRJEB67798",
    "PRJEB82425"
  ),

  V2_Role = c(
    "CORE_OR_SUPPORT_LONGITUDINAL",
    "NONSEPSIS_LONGITUDINAL_CONTROL",
    "ICU_LONGITUDINAL_EXTERNAL_VALIDATION",
    "NEONATAL_BSI_LONGITUDINAL_EXTERNAL_VALIDATION",
    "STATIC_ICU_SEPSIS_BACKGROUND",
    "SURGICAL_LONGITUDINAL_CONTROL",
    "ICU_INFECTION_LONGITUDINAL_EXTERNAL_VALIDATION"
  ),

  Status_After_67 = c(
    "READY_FROM_66",
    "READY_FROM_66",
    map595_status,
    ifelse(
      nrow(resolved884) > 0 && sum(is.na(resolved884$Run_ID)) == 0,
      "READY_ALL_RUNS_RESOLVED",
      "READY_WITH_NONSEQUENCED_OR_UNRESOLVED_SAMPLES"
    ),
    ifelse(!is.null(mofa_meta), "PUBLIC_METADATA_RECOVERED_NOT_CORE_LONGITUDINAL", "PUBLIC_METADATA_DOWNLOAD_FAILED"),
    paste0("SUPPLEMENT_", ifelse(nrow(supp677_profiles_df) > 0, "PARSED", "NOT_PARSED"),
           "; HUMAN_FILTER_", human_filter_status),
    ifelse(nrow(docx824_summary) > 0, "SUPPLEMENT_DOCX_PARSED", "SUPPLEMENT_DOCX_NOT_PARSED")
  ),

  Manual_Action_Now = c(
    "NONE",
    "NONE",
    ifelse(
      nrow(map595) > 0,
      "NONE",
      "NONE YET - inspect downloaded release metadata first"
    ),
    ifelse(
      nrow(resolved884) > 0 && sum(is.na(resolved884$Run_ID)) > 0,
      "NONE - unmatched rows may simply lack deposited sequencing runs",
      "NONE"
    ),
    "NONE - use only as static/background validation",
    "NONE YET - inspect parsed human supplement before author contact",
    "NONE YET - inspect parsed DOCX; manual figure reading only if no structured patient linkage exists"
  )
)

# ------------------------------------------------------------
# H. Save outputs
# ------------------------------------------------------------
stamp <- format(Sys.time(), "%Y%m%d_%H%M%S")

if (nrow(release_log_df) > 0) {
  write_excel_csv(
    release_log_df,
    file.path(RESULT_ROOT, paste0("PRJNA595346_release_download_log_", stamp, ".csv")),
    na = ""
  )
}

if (nrow(map595) > 0) {
  write_excel_csv(
    map595,
    file.path(RESULT_ROOT, paste0("PRJNA595346_sample_patient_time_map_", stamp, ".csv")),
    na = ""
  )
}

if (nrow(resolved884) > 0) {
  write_excel_csv(
    resolved884,
    file.path(RESULT_ROOT, paste0("PRJNA884103_secondpass_run_mapping_", stamp, ".csv")),
    na = ""
  )
}

if (nrow(supp677_profiles_df) > 0) {
  write_excel_csv(
    supp677_profiles_df,
    file.path(RESULT_ROOT, paste0("PRJEB67798_supplement_profile_", stamp, ".csv")),
    na = ""
  )
}

write_excel_csv(
  ready_qc,
  file.path(RESULT_ROOT, paste0("V2_mapping_QC_after67_", stamp, ".csv")),
  na = ""
)

write_excel_csv(
  adjudication,
  file.path(RESULT_ROOT, paste0("V2_project_adjudication_after67_", stamp, ".csv")),
  na = ""
)

cat("\n============================================================\n")
cat("SEPSIS V2 - STEP 67 COMPLETE\n")
cat("R version:", R.version.string, "\n")
cat("============================================================\n\n")

cat("PRJNA595346:", map595_status, "\n")
cat("PRJNA884103 second-pass unresolved runs:",
    ifelse(nrow(resolved884) > 0, sum(is.na(resolved884$Run_ID)), NA), "\n")
cat("PRJEB37289 MOFA metadata:", mofa_status, "\n")
cat("PRJEB67798 supplement:", supp677_status, "\n")
cat("PRJEB67798 human filter:", human_filter_status, "\n")
cat("PRJEB82425 supplement:", supp824_status, "\n\n")

print(adjudication, n = Inf, width = Inf)

cat("\nOutputs saved to:\n", RESULT_ROOT, "\n")
cat("============================================================\n")
