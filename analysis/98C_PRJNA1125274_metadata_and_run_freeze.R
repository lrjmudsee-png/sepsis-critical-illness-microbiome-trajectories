# ============================================================
# Sepsis V2 Upgrade - Step 98C
# PRJNA1125274 (SURVEIL) metadata and run freeze
# ------------------------------------------------------------
# Purpose:
#   1) Verify the local FASTQ set against the official current
#      289-run record (NCB/ENA).
#   2) Re-derive patient/time mapping from the formal paper +
#      BioSample/SRA metadata + sample aliases (not from the old
#      STRICT_PROVISIONAL map).
#   3) Freeze the final 289-run manifest and patient-time map
#      BEFORE any sequence processing.
#
# Official design (from the SURVEIL paper, MicrobiologyOpen 2026,
# DOI 10.1002/mbo3.70301, PMC13132800):
#   - 132 hospitalized adults; two Italian hospitals:
#       Alessandria "Santi Antonio e Biagio e Cesare Arrigo" = SAL
#       Novara "Maggiore della Carita"                        = SNO
#   - three rectal-swab timepoints: T0 baseline (enrollment),
#     T1 at bacteremia/sepsis diagnosis, T2 at discharge/death
#   - sample alias "<patient#><letter>-<hospital>",
#     letter A/B/C mapped to T0/T1/T2 (design + chronological
#     collection-date validation performed here).
# ============================================================

options(stringsAsFactors = FALSE)
required_pkgs <- c("readr", "dplyr", "tidyr", "purrr", "tibble", "stringr",
                   "tools")
missing_pkgs <- required_pkgs[
  !vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_pkgs)) {
  install.packages(missing_pkgs, repos = "https://cloud.r-project.org")
}
suppressPackageStartupMessages({
  library(readr)
  library(tibble)
  library(purrr)
  library(stringr)
  library(tools)
  library(dplyr)
  library(tidyr)
})

ROOT <- "E:/sepsis_project"
DATA112 <- file.path(ROOT, "data", "PRJNA1125274")
META_DIR <- file.path(DATA112, "00_metadata")
OUT <- file.path(ROOT, "results", "V2_UPGRADE_20260907",
                 "98C_PRJNA1125274_METADATA_FREEZE")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

HISTORY <- file.path(ROOT, "results",
                     "V2_17_finalize_68229_adjudicate_1125274")

write_csv_safe <- function(x, p) write_excel_csv(x, p, na = "")
safe_csv <- function(p) {
  suppressMessages(read_csv(p, show_col_types = FALSE, progress = FALSE,
                            name_repair = "unique"))
}

# ------------------------------------------------------------
# 1. Load official run set (289) from two independent records
# ------------------------------------------------------------

manifest <- suppressMessages(read_delim(
  file.path(META_DIR, "correct_289_download",
            "PRJNA1125274_CORRECT_289_ENA_manifest.tsv"),
  delim = "\t", show_col_types = FALSE, progress = FALSE,
  name_repair = "unique")) |>
  select(run_accession, sample_accession, secondary_sample_accession,
         sample_alias, library_layout, library_strategy,
         fastq_ftp, fastq_md5, fastq_bytes)

official_list <- safe_csv(file.path(
  HISTORY, "PRJNA1125274_CORRECT_current_289_Run_list.csv"))

stopifnot(nrow(manifest) == 289)
cat("ENA manifest runs:", nrow(manifest), "\n")
cat("Frozen official run-list rows:", nrow(official_list), "\n")

setdiff_official <- setdiff(official_list$Run_ID, manifest$run_accession)
setdiff_manifest <- setdiff(manifest$run_accession, official_list$Run_ID)
cat("runs in official list but not ENA manifest:", length(setdiff_official), "\n")
cat("runs in ENA manifest but not official list:", length(setdiff_manifest), "\n")
if (length(setdiff_official) || length(setdiff_manifest)) {
  stop("Run-set mismatch between official list and ENA manifest.")
}

# ------------------------------------------------------------
# 2. Local FASTQ status (size + md5) for all 578 expected files
# ------------------------------------------------------------

raw_dir <- file.path(DATA112, "01_correct_289_raw")

expected <- manifest |>
  rowwise() |>
  mutate(
    fastq_list = list(str_split(fastq_ftp, ";")[[1]]),
    md5_list = list(str_split(fastq_md5, ";")[[1]]),
    bytes_list = list(as.numeric(str_split(fastq_bytes, ";")[[1]]))
  ) |>
  ungroup() |>
  mutate(file_index = map(fastq_list, seq_along)) |>
  unnest(c(fastq_list, md5_list, bytes_list, file_index)) |>
  mutate(
    file_name = basename(fastq_list),
    url = paste0("https://", fastq_list),
    local_path = file.path(raw_dir, file_name)
  ) |>
  select(run_accession, sample_alias, file_index, file_name, local_path,
         expected_bytes = bytes_list, expected_md5 = md5_list)

expected$local_exists <- file.exists(expected$local_path)
expected$local_bytes <- ifelse(
  expected$local_exists,
  file.info(expected$local_path)$size,
  NA_real_)
expected$local_md5 <- ifelse(
  expected$local_exists,
  unname(tools::md5sum(expected$local_path)),
  NA_character_)

expected <- expected |>
  mutate(
    size_match = !is.na(local_bytes) & (local_bytes == expected_bytes),
    md5_match = !is.na(local_md5) & (local_md5 == expected_md5),
    file_status = case_when(
      !local_exists ~ "MISSING",
      !md5_match ~ "MD5_MISMATCH",
      TRUE ~ "OK"
    )
  )

cat("Expected files:", nrow(expected), "\n")
cat("Local file status:\n")
print(table(expected$file_status))

# ------------------------------------------------------------
# 3. Alias parsing and hospital namespace
# ------------------------------------------------------------

manifest_parsed <- manifest |>
  mutate(
    alias_number = as.integer(str_extract(sample_alias, "^\\d+")),
    alias_letter = str_extract(sample_alias, "^\\d+([A-Za-z])", group = 1),
    alias_hospital = str_extract(sample_alias, "-([A-Za-z0-9_]+)$", group = 1)
  )

hospital_code <- c(SAL = "ALESSANDRIA", SNO = "NOVARA")
stopifnot(all(manifest_parsed$alias_hospital %in% names(hospital_code)))

manifest_parsed <- manifest_parsed |>
  mutate(
    hospital = unname(hospital_code[alias_hospital]),
    patient_id = paste0(alias_hospital, "_", alias_number)
  )

letter_map <- c(A = "T0", B = "T1", C = "T2")
manifest_parsed$timepoint <- unname(letter_map[manifest_parsed$alias_letter])

# per-patient expected letters (unique, chronological)
patient_summary <- manifest_parsed |>
  distinct(patient_id, alias_letter) |>
  arrange(patient_id, alias_letter) |>
  group_by(patient_id) |>
  summarise(letter_combo = paste(alias_letter, collapse = ""), .groups = "drop")

cat("Patients total:", n_distinct(manifest_parsed$patient_id), "\n")
cat("Patients by hospital:\n")
print(manifest_parsed |> distinct(patient_id, hospital) |>
        count(hospital))
cat("Per-patient letter patterns:\n")
print(table(patient_summary$letter_combo))

# ------------------------------------------------------------
# 4. Join sample attributes (collection date, geo) by BioSample
# ------------------------------------------------------------

xml_sum <- safe_csv(file.path(
  META_DIR, "step76_sample_xml", "_sample_xml_summary.csv"))
# each unique sample appears twice (SAMN + SRS secondary); keep one row
xml1 <- xml_sum |>
  rename(alias = Alias, collection_date = CollectionDate, geo = Geo) |>
  mutate(alias = str_trim(alias)) |>
  distinct(alias, .keep_all = TRUE)

manifest_parsed <- manifest_parsed |>
  left_join(
    xml1 |> select(alias, collection_date, geo),
    by = c("sample_alias" = "alias"))

# ------------------------------------------------------------
# 5. Chronological validation of letters within patient
# ------------------------------------------------------------

# per-patient, per-letter minimum collection date (NA when letter absent)
patient_letter_dates <- manifest_parsed |>
  filter(!is.na(collection_date)) |>
  mutate(collection_date = as.Date(collection_date)) |>
  group_by(patient_id, alias_letter) |>
  summarise(letter_min_date = min(collection_date, na.rm = TRUE),
            letter_max_date = max(collection_date, na.rm = TRUE),
            .groups = "drop") |>
  pivot_wider(
    id_cols = patient_id,
    names_from = alias_letter,
    values_from = letter_min_date
  )

violations <- patient_letter_dates |>
  mutate(
    has_A = !is.na(A), has_B = !is.na(B), has_C = !is.na(C),
    ok_AB = ifelse(has_A & has_B, A <= B, TRUE),
    ok_AC = ifelse(has_A & has_C, A <= C, TRUE),
    ok_BC = ifelse(has_B & has_C, B <= C, TRUE),
    order_ok = ok_AB & ok_AC & ok_BC,
    order_violation = !order_ok
  )

cat("Patients with >=2 letters:", sum(violations$has_A + violations$has_B +
                                        violations$has_C >= 2), "\n")
cat("Date-order violations (letters not chronological):",
    sum(violations$order_violation, na.rm = TRUE), "\n")
if (any(violations$order_violation, na.rm = TRUE)) {
  print(violations |> filter(order_violation))
}

# diagnostic: letter -> min date per patient (export for review)
patient_letter_diag <- manifest_parsed |>
  filter(!is.na(collection_date)) |>
  mutate(collection_date = as.Date(collection_date)) |>
  group_by(patient_id, alias_letter) |>
  summarise(
    n_runs = n(),
    min_date = min(collection_date),
    max_date = max(collection_date),
    .groups = "drop")
write_csv_safe(patient_letter_diag, file.path(
  OUT, "PRJNA1125274_patient_letter_date_diagnostic.csv"))

# patients without A (no baseline) cannot anchor T0-based trajectory
manifest_parsed <- manifest_parsed |>
  left_join(patient_summary, by = "patient_id") |>
  mutate(
    has_T0 = grepl("A", letter_combo),
    has_T1 = grepl("B", letter_combo),
    has_T2 = grepl("C", letter_combo),
    complete_T0_T1_T2 = has_T0 & has_T1 & has_T2
  )

cat("Patients with complete T0/T1/T2 (letters ABC):",
    n_distinct(manifest_parsed$patient_id[manifest_parsed$complete_T0_T1_T2]), "\n")
cat("Patients with T0+T1 (>= A+B):",
    n_distinct(manifest_parsed$patient_id[
      manifest_parsed$has_T0 & manifest_parsed$has_T1]), "\n")
cat("Patients with T0+T2 (>= A+C):",
    n_distinct(manifest_parsed$patient_id[
      manifest_parsed$has_T0 & manifest_parsed$has_T2]), "\n")

# ---- QC: every patient's letters are a subsequence of A-B-C ----
subseq_ok <- all(
  grepl("^A?B?C?$", patient_summary$letter_combo)
)
cat("All patient letter patterns are increasing subsequences of A,B,C:",
    subseq_ok, "\n")

# ------------------------------------------------------------
# 6. Final patient-time map
# ------------------------------------------------------------

final_map <- manifest_parsed |>
  left_join(
    violations[, c("patient_id", "order_violation"), drop = FALSE],
    by = "patient_id"
  ) |>
  mutate(
    mapping_status = case_when(
      !is.na(order_violation) & order_violation ~
        "LETTER_ORDER_TIME_CODE_DATE_ANOMALY",
      TRUE ~ "LETTER_ORDER_TIME_CODE"
    ),
    timepoint_notes = paste0(
      "letter=", alias_letter, "; hospital=", hospital,
      "; letter codes are increasing subsequences of A/B/C for every patient",
      "; design: A=T0 baseline, B=T1 sepsis/bacteremia diagnosis, ",
      "C=T2 discharge/death (SURVEIL paper)"
    ),
    project = "PRJNA1125274",
    run_id = run_accession,
    biosample = sample_accession
  )

keep_cols <- intersect(
  c("project", "run_id", "biosample", "secondary_sample_accession",
    "sample_alias", "hospital", "patient_id", "alias_number",
    "alias_letter", "timepoint", "collection_date", "geo",
    "letter_combo", "has_T0", "has_T1", "has_T2",
    "complete_T0_T1_T2", "mapping_status", "timepoint_notes"),
  names(final_map))
missing_cols <- setdiff(
  c("project", "run_id", "biosample", "secondary_sample_accession",
    "sample_alias", "hospital", "patient_id", "alias_number",
    "alias_letter", "timepoint", "collection_date", "geo",
    "letter_combo", "has_T0", "has_T1", "has_T2",
    "complete_T0_T1_T2", "mapping_status", "timepoint_notes"),
  names(final_map))
if (length(missing_cols)) stop("final_map missing columns: ",
                               paste(missing_cols, collapse = ", "))
final_map <- final_map[, keep_cols, drop = FALSE]

write_csv_safe(final_map, file.path(
  OUT, "PRJNA1125274_FINAL_patient_time_map.csv"))

# ------------------------------------------------------------
# 7. Final 289-run manifest
# ------------------------------------------------------------

run_status <- expected |>
  group_by(run_accession) |>
  summarise(
    n_files_expected = n(),
    n_files_ok = sum(file_status == "OK"),
    n_missing = sum(file_status == "MISSING"),
    n_md5_mismatch = sum(file_status == "MD5_MISMATCH"),
    total_expected_bytes = sum(expected_bytes),
    total_local_bytes = sum(local_bytes, na.rm = TRUE),
    run_ready = n_files_ok == n_files_expected,
    .groups = "drop"
  )

manifest_ready <- manifest |>
  left_join(run_status, by = "run_accession") |>
  left_join(
    final_map[, c("run_id", "hospital", "patient_id", "alias_number",
                  "alias_letter", "timepoint", "collection_date",
                  "complete_T0_T1_T2", "mapping_status"), drop = FALSE],
    by = c("run_accession" = "run_id")
  ) |>
  rowwise() |>
  mutate(
    r1_path = ifelse(file.exists(file.path(raw_dir, paste0(run_accession, "_1.fastq.gz"))),
                     file.path(raw_dir, paste0(run_accession, "_1.fastq.gz")), NA_character_),
    r2_path = ifelse(file.exists(file.path(raw_dir, paste0(run_accession, "_2.fastq.gz"))),
                     file.path(raw_dir, paste0(run_accession, "_2.fastq.gz")), NA_character_)
  ) |>
  ungroup() |>
  mutate(
    inclusion = case_when(
      run_ready ~ "INCLUDE",
      TRUE ~ "EXCLUDE_INCOMPLETE_DOWNLOAD"
    ),
    exclusion_reason = case_when(
      run_ready ~ "",
      TRUE ~ "FASTQ missing or MD5 mismatch"
    )
  )

manifest_cols <- intersect(
  c("run_accession", "biosample", "secondary_sample_accession",
    "sample_alias", "library_layout", "patient_id", "alias_number",
    "alias_letter", "timepoint", "hospital", "collection_date",
    "complete_T0_T1_T2", "mapping_status", "fastq_r1_path",
    "fastq_r2_path", "n_files_ok", "n_missing", "run_ready",
    "inclusion", "exclusion_reason"),
  names(manifest_ready))

# build expected biosample/r1/r2 columns
manifest_ready <- manifest_ready |>
  mutate(
    biosample = sample_accession,
    fastq_r1_path = r1_path,
    fastq_r2_path = r2_path
  ) |>
  as.data.frame()

manifest_ready <- manifest_ready[, intersect(
  c("run_accession", "biosample", "secondary_sample_accession",
    "sample_alias", "library_layout", "patient_id", "alias_number",
    "alias_letter", "timepoint", "hospital", "collection_date",
    "complete_T0_T1_T2", "mapping_status", "fastq_r1_path",
    "fastq_r2_path", "n_files_ok", "n_missing", "run_ready",
    "inclusion", "exclusion_reason"), names(manifest_ready)), drop = FALSE]

write_csv_safe(manifest_ready, file.path(
  OUT, "PRJNA1125274_FINAL_289_run_manifest.csv"))

# copy frozen outputs into the local data folder (authoritative location)
write_csv_safe(final_map, file.path(DATA112,
                                    "PRJNA1125274_FINAL_patient_time_map.csv"))
write_csv_safe(manifest_ready, file.path(DATA112,
                                         "PRJNA1125274_FINAL_289_run_manifest.csv"))

cat("Runs INCLUDE:", sum(manifest_ready$inclusion == "INCLUDE"), "\n")
cat("Runs EXCLUDE:", sum(manifest_ready$inclusion == "EXCLUDE_INCOMPLETE_DOWNLOAD"), "\n")

# ------------------------------------------------------------
# QC output
# ------------------------------------------------------------

qc_summary <- tibble(
  metric = c("official_current_ENA_runs", "local_runs_complete",
             "local_runs_incomplete", "samples_total",
             "patients_total", "patients_ALESSANDRIA", "patients_NOVARA",
             "samples_T0", "samples_T1", "samples_T2",
             "patients_T0_T1", "patients_T0_T2", "patients_complete_T0_T1_T2",
             "date_order_violations"),
  value = c(
    nrow(manifest),
    sum(run_status$run_ready),
    sum(!run_status$run_ready),
    nrow(manifest),
    n_distinct(manifest_parsed$patient_id),
    n_distinct(manifest_parsed$patient_id[manifest_parsed$hospital == "ALESSANDRIA"]),
    n_distinct(manifest_parsed$patient_id[manifest_parsed$hospital == "NOVARA"]),
    sum(manifest_parsed$alias_letter == "A"),
    sum(manifest_parsed$alias_letter == "B"),
    sum(manifest_parsed$alias_letter == "C"),
    n_distinct(manifest_parsed$patient_id[manifest_parsed$has_T0 &
                                            manifest_parsed$has_T1]),
    n_distinct(manifest_parsed$patient_id[manifest_parsed$has_T0 &
                                            manifest_parsed$has_T2]),
    n_distinct(manifest_parsed$patient_id[manifest_parsed$complete_T0_T1_T2]),
    sum(violations$order_violation, na.rm = TRUE)
  )
)

write_csv_safe(qc_summary, file.path(
  OUT, "PRJNA1125274_metadata_and_run_QC_summary.csv"))

cat("DONE 98C (pre-sequence freeze stage; FASTQ readiness depends on\n")
cat("supplement download completion - rerun after aria2 finishes if runs missing)\n")
