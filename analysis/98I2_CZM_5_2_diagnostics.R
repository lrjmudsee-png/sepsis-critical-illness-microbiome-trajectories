# Diagnostic follow-up for Step98I: identify patient/run loss under the
# historical addendum's total-count >=5 / prevalence >=2 feature rule.
# No FASTQ, DADA2, patient mapping, or frozen results are modified.

options(stringsAsFactors = FALSE)
suppressPackageStartupMessages({
  library(tidyr)
  library(tibble)
  library(readr)
  library(zCompositions)
  library(dplyr)
})

ROOT <- "E:/sepsis_project"
OUT <- file.path(ROOT, "results", "V2_UPGRADE_20260910_REMAINING_VALIDATION")
OBJ <- file.path(ROOT, "data", "PRJNA1125274", "03_dada2",
                 "PRJNA1125274_analysis_object_external_validation.rds")
DATE_ANOMALY <- "LETTER_ORDER_TIME_CODE_DATE_ANOMALY"
CACHE <- file.path(OUT, ".CZM_5_2_imputation_cache.rds")

obj <- readRDS(OBJ)
counts <- obj$counts
meta <- as.data.frame(obj$metadata)
meta$run_id <- as.character(meta$run_accession)
meta$patient_id <- as.character(meta$patient_id)
meta$depth <- rowSums(counts[meta$run_id, , drop = FALSE])

fw <- meta |>
  filter(timepoint %in% c("T0", "T1", "T2"), !is.na(patient_id), depth >= 2000) |>
  select(run_id, patient_id, hospital, timepoint, mapping_status, depth)

cnt <- counts[fw$run_id, , drop = FALSE]
keep <- colSums(cnt) >= 5 & colSums(cnt > 0) >= 2
cnt <- cnt[, keep, drop = FALSE]

warning_messages <- character()
if (file.exists(CACHE)) {
  cached <- readRDS(CACHE)
  imp <- cached$imp
  warning_messages <- cached$warning_messages
} else {
  imp <- withCallingHandlers(
    zCompositions::cmultRepl(cnt, label = 0, method = "CZM",
                            output = "p-counts", z.warning = 1,
                            z.delete = FALSE, suppress.print = TRUE),
    warning = function(w) {
      warning_messages <<- c(warning_messages, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )
  saveRDS(list(imp = imp, warning_messages = warning_messages), CACHE)
}

clr <- log(imp) - rowMeans(log(imp))
aitch <- as.matrix(dist(clr, method = "euclidean"))

base <- fw |>
  filter(timepoint == "T0") |>
  distinct(patient_id, .keep_all = TRUE) |>
  select(patient_id, base_run = run_id)

diag <- fw |>
  left_join(base, by = "patient_id") |>
  mutate(
    retained_reads = rowSums(cnt[run_id, , drop = FALSE]),
    n_present_features = rowSums(cnt[run_id, , drop = FALSE] > 0),
    zero_fraction = 1 - n_present_features / ncol(cnt),
    imputation_row_all_finite = apply(imp[run_id, , drop = FALSE], 1,
                                      function(x) all(is.finite(x))),
    clr_row_all_finite = apply(clr[run_id, , drop = FALSE], 1,
                               function(x) all(is.finite(x))),
    aitchison_from_T0_CZM = ifelse(is.na(base_run), NA_real_,
                                   aitch[cbind(match(run_id, rownames(aitch)),
                                                match(base_run, colnames(aitch)))]),
    distance_finite = is.finite(aitchison_from_T0_CZM)
  )

complete_ids <- diag |>
  group_by(patient_id) |>
  summarise(complete = all(c("T0", "T1", "T2") %in% timepoint),
            date_anomaly = any(mapping_status == DATE_ANOMALY), .groups = "drop") |>
  filter(complete)

patient_diag <- diag |>
  filter(patient_id %in% complete_ids$patient_id) |>
  left_join(complete_ids, by = "patient_id") |>
  group_by(patient_id, hospital, date_anomaly) |>
  summarise(
    T0_run = run_id[match("T0", timepoint)],
    T1_run = run_id[match("T1", timepoint)],
    T2_run = run_id[match("T2", timepoint)],
    T0_distance_finite = distance_finite[match("T0", timepoint)],
    T1_distance_finite = distance_finite[match("T1", timepoint)],
    T2_distance_finite = distance_finite[match("T2", timepoint)],
    T0_zero_fraction = zero_fraction[match("T0", timepoint)],
    T1_zero_fraction = zero_fraction[match("T1", timepoint)],
    T2_zero_fraction = zero_fraction[match("T2", timepoint)],
    T0_retained_reads = retained_reads[match("T0", timepoint)],
    T1_retained_reads = retained_reads[match("T1", timepoint)],
    T2_retained_reads = retained_reads[match("T2", timepoint)],
    complete_aitchison_pair = T1_distance_finite & T2_distance_finite,
    in_date_clean_primary = !first(date_anomaly),
    .groups = "drop"
  ) |>
  arrange(complete_aitchison_pair, hospital, patient_id)

write_excel_csv(patient_diag,
                file.path(OUT, "PRJNA1125274_CZM_5_2_patient_diagnostics.csv"), na = "")
write_excel_csv(diag |> filter(!distance_finite | !imputation_row_all_finite | !clr_row_all_finite),
                file.path(OUT, "PRJNA1125274_CZM_5_2_nonfinite_run_diagnostics.csv"), na = "")

warning_summary <- tibble(message = warning_messages) |>
  count(message, sort = TRUE, name = "n_occurrences")
write_excel_csv(warning_summary,
                file.path(OUT, "PRJNA1125274_CZM_5_2_warning_summary.csv"), na = "")
unlink(CACHE, force = TRUE)

cat("CZM 5/2 diagnostic complete\n")
cat("Features:", ncol(cnt), "\n")
cat("Complete patients with finite T1/T2:", sum(patient_diag$complete_aitchison_pair), "/", nrow(patient_diag), "\n")
cat("Date-clean primary with finite T1/T2:",
    sum(patient_diag$complete_aitchison_pair & patient_diag$in_date_clean_primary), "/",
    sum(patient_diag$in_date_clean_primary), "\n")
