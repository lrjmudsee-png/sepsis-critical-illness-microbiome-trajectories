# ============================================================
# Sepsis V2 - Step93T10B
# FORCED RUN-LEVEL BRIDGE
#
# Core strategy:
#   EII sample/run ID
#        -> Step93Q Run
#        -> Patient_true
#        -> Step93S2 trajectory
#
# This intentionally does NOT trust EII patient_id naming.
# It also compares the corrected Step92B2 EII table against the older
# Step92B table and selects the table/identifier column with the
# highest exact Run overlap with PRJNA691455 Step93Q.
# ============================================================

options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(stringr)
  library(tibble)
})

ROOT <- "E:/sepsis_project"
PROJECT <- "PRJNA691455"

OUT <- file.path(
  ROOT,
  "results",
  "V2_33T10B_FORCED_RUN_BRIDGE"
)

dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

logfile <- file.path(OUT, "_STEP93T10B_runtime.txt")
writeLines(paste("START", Sys.time()), logfile)

logmsg <- function(...) {
  z <- paste0(...)
  cat(z, "\n")
  cat(z, "\n", file = logfile, append = TRUE)
}

norm_id <- function(x) {
  x <- toupper(str_trim(as.character(x)))
  x[x %in% c("", "NA", "N/A", "NULL")] <- NA_character_
  x
}

norm_patient <- function(x) {
  x <- norm_id(x)
  x <- str_replace(x, regex("^PRJNA691455[_:\\-]*", ignore_case = TRUE), "")
  x <- str_replace(x, regex("FFO.*$", ignore_case = TRUE), "")
  x
}

pick_exact <- function(nms, priority) {
  x <- priority[priority %in% nms]
  if (length(x) == 0) return(NA_character_)
  x[1]
}

safe_read <- function(f) {
  tryCatch(
    read_csv(f, show_col_types = FALSE, guess_max = 20000),
    error = function(e) NULL
  )
}

# ------------------------------------------------------------
# 1. Step93Q Run -> Patient_true bridge
# ------------------------------------------------------------

q_file <- file.path(
  ROOT,
  "results",
  "V2_33Q_TRAJECTORY_PATIENT_ID_RECONSTRUCTION_FIX",
  "V2_STEP93Q_fixed_longitudinal_genus_profile.csv"
)

if (!file.exists(q_file)) {
  q_candidates <- list.files(
    file.path(ROOT, "results"),
    recursive = TRUE,
    full.names = TRUE,
    pattern = "V2_STEP93Q_fixed_longitudinal_genus_profile\\.csv$"
  )
  if (length(q_candidates) == 0) stop("Cannot find Step93Q fixed longitudinal genus profile.")
  q_file <- q_candidates[1]
}

q <- safe_read(q_file)
if (is.null(q)) stop("Cannot read Step93Q file.")

q_run_col <- pick_exact(
  names(q),
  c("Run", "Run_ID", "run_id", "run", "RunID", "run_accession")
)

if (is.na(q_run_col)) {
  tmp <- names(q)[str_detect(names(q), regex("^run|run_", ignore_case = TRUE))]
  if (length(tmp) > 0) q_run_col <- tmp[1]
}

q_patient_col <- pick_exact(
  names(q),
  c("Patient_true", "patient_true", "Patient", "patient_id")
)

if (is.na(q_run_col) || is.na(q_patient_col)) {
  stop(
    paste0(
      "Step93Q bridge columns not found. Columns: ",
      paste(names(q), collapse = ", ")
    )
  )
}

q_bridge <- q %>%
  transmute(
    run_key = norm_id(.data[[q_run_col]]),
    Patient_true = norm_patient(.data[[q_patient_col]])
  ) %>%
  filter(!is.na(run_key), !is.na(Patient_true)) %>%
  distinct()

write_csv(q_bridge, file.path(OUT, "01_STEP93Q_run_to_Patient_true_bridge.csv"))

logmsg("Step93Q bridge rows: ", nrow(q_bridge))
logmsg("Step93Q unique runs: ", n_distinct(q_bridge$run_key))
logmsg("Step93Q unique patients: ", n_distinct(q_bridge$Patient_true))

# ------------------------------------------------------------
# 2. Candidate EII tables
# ------------------------------------------------------------

eii_files <- c(
  file.path(
    ROOT,
    "results",
    "V2_32B2_STEP92B_BASELINE_HARMONIZATION_FIXED",
    "V2_STEP92B2_corrected_timepoint_EII.csv"
  ),
  file.path(
    ROOT,
    "results",
    "V2_32B_STEP92B_EII_CONSTRUCTION",
    "V2_STEP92B_patient_timepoint_EII.csv"
  )
)

eii_files <- eii_files[file.exists(eii_files)]

if (length(eii_files) == 0) {
  stop("No Step92B/Step92B2 EII table found.")
}

# ------------------------------------------------------------
# 3. For each EII file, test ALL plausible identifier columns
#    against Step93Q Run IDs
# ------------------------------------------------------------

audit_list <- list()
loaded <- list()

for (f in eii_files) {

  x <- safe_read(f)
  if (is.null(x)) next

  # Project filter when explicit project column exists
  pcol <- pick_exact(names(x), c("project", "Project", "PROJECT"))
  if (!is.na(pcol)) {
    xp <- x %>%
      filter(toupper(as.character(.data[[pcol]])) == PROJECT)

    # avoid accidentally zeroing a table because project labels differ
    if (nrow(xp) > 0) x <- xp
  }

  id_cols <- names(x)[
    str_detect(
      names(x),
      regex("run|sample|accession|patient|subject|(^id$)|_id$", ignore_case = TRUE)
    )
  ]

  if (length(id_cols) == 0) id_cols <- names(x)

  for (cc in id_cols) {
    vals <- unique(norm_id(x[[cc]]))
    vals <- vals[!is.na(vals)]

    overlap <- intersect(vals, q_bridge$run_key)

    audit_list[[length(audit_list) + 1]] <- tibble(
      eii_file = f,
      eii_rows = nrow(x),
      identifier_column = cc,
      unique_identifier_values = length(vals),
      exact_run_overlap = length(overlap),
      overlap_fraction_of_Q_runs = length(overlap) / n_distinct(q_bridge$run_key)
    )
  }

  loaded[[f]] <- x
}

audit <- bind_rows(audit_list) %>%
  arrange(desc(exact_run_overlap), desc(overlap_fraction_of_Q_runs))

write_csv(audit, file.path(OUT, "02_EII_identifier_run_overlap_audit.csv"))

if (nrow(audit) == 0 || max(audit$exact_run_overlap, na.rm = TRUE) == 0) {

  # Dump identifiers for diagnosis and stop safely
  writeLines(
    c(
      "No EII identifier column had any exact overlap with Step93Q Run IDs.",
      "Inspect 01_STEP93Q_run_to_Patient_true_bridge.csv",
      "Inspect 02_EII_identifier_run_overlap_audit.csv"
    ),
    file.path(OUT, "_STEP93T10B_ZERO_RUN_OVERLAP.txt")
  )

  stop("STEP93T10B found 0 Run-level overlap. Diagnostics written.")
}

best <- audit[1, ]

selected_file <- best$eii_file[[1]]
selected_id_col <- best$identifier_column[[1]]
selected_overlap <- best$exact_run_overlap[[1]]

eii <- loaded[[selected_file]]

logmsg("Selected EII file: ", selected_file)
logmsg("Selected EII identifier column: ", selected_id_col)
logmsg("Exact Step93Q Run overlap: ", selected_overlap)

write_csv(best, file.path(OUT, "03_SELECTED_EII_SOURCE_AND_ID_COLUMN.csv"))

# ------------------------------------------------------------
# 4. Force EII -> Step93Q join by Run
# ------------------------------------------------------------

eii2 <- eii %>%
  mutate(run_key = norm_id(.data[[selected_id_col]]))

mapped_rows <- eii2 %>%
  inner_join(q_bridge, by = "run_key")

unmapped_eii <- eii2 %>%
  anti_join(q_bridge, by = "run_key")

write_csv(mapped_rows, file.path(OUT, "04_EII_rows_mapped_to_Patient_true.csv"))

write_csv(
  unmapped_eii,
  file.path(OUT, "05_EII_rows_not_in_PRJNA691455_Step93Q.csv")
)

logmsg("Mapped EII rows: ", nrow(mapped_rows))
logmsg("Mapped Patient_true: ", n_distinct(mapped_rows$Patient_true))

# ------------------------------------------------------------
# 5. Explicit EII column
# ------------------------------------------------------------

eii_priority <- c(
  "EII_0_100",
  "EII",
  "eii",
  "EII_raw",
  "ecological_instability_index",
  "instability_index"
)

eii_col <- pick_exact(names(mapped_rows), eii_priority)

if (is.na(eii_col)) {
  ec <- names(mapped_rows)[
    str_detect(names(mapped_rows), regex("EII|instability", ignore_case = TRUE))
  ]
  ec <- ec[
    vapply(mapped_rows[ec], is.numeric, logical(1))
  ]
  if (length(ec) == 0) {
    stop("Mapped EII table has no recognizable EII numeric column.")
  }
  eii_col <- ec[1]
}

mapped_rows$EII_value <- suppressWarnings(as.numeric(mapped_rows[[eii_col]]))

writeLines(
  paste0("Selected EII column: ", eii_col),
  file.path(OUT, "06_SELECTED_EII_VALUE_COLUMN.txt")
)

# ------------------------------------------------------------
# 6. Time order
# ------------------------------------------------------------

time_col <- pick_exact(
  names(mapped_rows),
  c(
    "analysis_time_order",
    "time_order",
    "time_numeric",
    "time_day",
    "day",
    "Time_numeric"
  )
)

if (!is.na(time_col)) {
  mapped_rows$time_order_final <- suppressWarnings(
    as.numeric(mapped_rows[[time_col]])
  )
} else {

  factor_col <- pick_exact(
    names(mapped_rows),
    c("time_factor", "time_raw", "timepoint", "Timepoint")
  )

  if (!is.na(factor_col)) {
    tf <- tolower(as.character(mapped_rows[[factor_col]]))
    ord <- suppressWarnings(as.numeric(str_extract(tf, "\\d+")))
    ord[str_detect(tf, "base|inclusion|admission")] <- 0
    mapped_rows$time_order_final <- ord
  } else {
    mapped_rows <- mapped_rows %>%
      group_by(Patient_true) %>%
      mutate(time_order_final = row_number()) %>%
      ungroup()
  }
}

# If some time orders remain NA, preserve source row order within patient
mapped_rows <- mapped_rows %>%
  group_by(Patient_true) %>%
  mutate(
    time_order_final = ifelse(
      is.na(time_order_final),
      1000 + row_number(),
      time_order_final
    )
  ) %>%
  ungroup()

# ------------------------------------------------------------
# 7. Patient-level EII
# ------------------------------------------------------------

patient_eii <- mapped_rows %>%
  filter(!is.na(Patient_true), !is.na(EII_value)) %>%
  arrange(Patient_true, time_order_final) %>%
  group_by(Patient_true) %>%
  summarise(
    n_EII_timepoints = n(),
    EII_baseline = first(EII_value),
    EII_last = last(EII_value),
    EII_change = EII_last - EII_baseline,
    EII_mean = mean(EII_value, na.rm = TRUE),
    EII_max = max(EII_value, na.rm = TRUE),
    .groups = "drop"
  )

write_csv(
  patient_eii,
  file.path(OUT, "07_PRJNA691455_patient_EII_by_forced_Run_bridge.csv")
)

# ------------------------------------------------------------
# 8. Load trajectory clusters
# ------------------------------------------------------------

cluster_file <- file.path(
  ROOT,
  "results",
  "V2_33S2_TRAJECTORY_DISTANCE_METHOD_FIX",
  "trajectory_clusters.csv"
)

if (!file.exists(cluster_file)) {
  stop("Missing Step93S2 trajectory_clusters.csv")
}

clusters <- read_csv(cluster_file, show_col_types = FALSE)

cp <- pick_exact(
  names(clusters),
  c("Patient_true", "patient_true", "Patient", "patient_id")
)

if (is.na(cp)) stop("Cannot detect trajectory patient column.")

clusters2 <- clusters %>%
  transmute(
    Patient_true = norm_patient(.data[[cp]]),
    cluster = as.integer(cluster)
  ) %>%
  distinct(Patient_true, .keep_all = TRUE)

# ------------------------------------------------------------
# 9. Optional signed trajectory distance summary
# ------------------------------------------------------------

dist_file <- file.path(
  ROOT,
  "results",
  "V2_33S2_TRAJECTORY_DISTANCE_METHOD_FIX",
  "signed_trajectory_euclidean_distance_matrix.csv"
)

dist_summary <- NULL

if (file.exists(dist_file)) {
  dm <- read.csv(dist_file, row.names = 1, check.names = FALSE)
  dm <- as.matrix(dm)
  storage.mode(dm) <- "numeric"

  dist_summary <- tibble(
    Patient_true = norm_patient(rownames(dm)),
    mean_signed_trajectory_distance = rowMeans(dm, na.rm = TRUE),
    median_signed_trajectory_distance = apply(dm, 1, median, na.rm = TRUE)
  )
}

# ------------------------------------------------------------
# 10. Final integration
# ------------------------------------------------------------

integrated <- clusters2 %>%
  left_join(patient_eii, by = "Patient_true")

if (!is.null(dist_summary)) {
  integrated <- integrated %>%
    left_join(dist_summary, by = "Patient_true")
}

write_csv(
  integrated,
  file.path(OUT, "08_FINAL_patient_EII_trajectory_integrated.csv")
)

matched <- integrated %>%
  filter(!is.na(EII_mean))

unmatched <- integrated %>%
  filter(is.na(EII_mean))

write_csv(matched, file.path(OUT, "09_FINAL_matched_patients.csv"))
write_csv(unmatched, file.path(OUT, "10_FINAL_unmatched_patients.csv"))

# ------------------------------------------------------------
# 11. Exploratory cluster summaries/tests
# ------------------------------------------------------------

cluster_summary <- matched %>%
  group_by(cluster) %>%
  summarise(
    n = n(),
    median_EII_change = median(EII_change, na.rm = TRUE),
    median_EII_mean = median(EII_mean, na.rm = TRUE),
    median_EII_max = median(EII_max, na.rm = TRUE),
    .groups = "drop"
  )

write_csv(
  cluster_summary,
  file.path(OUT, "11_cluster_EII_descriptive.csv")
)

wilcox_out <- tibble()

if (
  n_distinct(matched$cluster) == 2 &&
  all(table(matched$cluster) >= 2)
) {
  metrics <- c("EII_change", "EII_mean", "EII_max")

  wilcox_out <- bind_rows(lapply(metrics, function(m) {
    xx <- matched %>%
      filter(!is.na(.data[[m]]))

    wt <- tryCatch(
      wilcox.test(xx[[m]] ~ xx$cluster, exact = FALSE),
      error = function(e) NULL
    )

    tibble(
      metric = m,
      n = nrow(xx),
      p_value = if (is.null(wt)) NA_real_ else wt$p.value
    )
  }))
}

write_csv(
  wilcox_out,
  file.path(OUT, "12_cluster_EII_wilcoxon_exploratory.csv")
)

# ------------------------------------------------------------
# 12. Trajectory distance vs EII exploratory correlations
# ------------------------------------------------------------

cor_out <- tibble()

if ("mean_signed_trajectory_distance" %in% names(matched)) {
  metrics <- c("EII_change", "EII_mean", "EII_max")

  cor_out <- bind_rows(lapply(metrics, function(m) {
    xx <- matched %>%
      select(mean_signed_trajectory_distance, all_of(m)) %>%
      filter(complete.cases(.))

    if (nrow(xx) < 3) {
      return(tibble(
        metric = m,
        n = nrow(xx),
        spearman_rho = NA_real_,
        p_value = NA_real_
      ))
    }

    ct <- suppressWarnings(
      cor.test(
        xx$mean_signed_trajectory_distance,
        xx[[m]],
        method = "spearman",
        exact = FALSE
      )
    )

    tibble(
      metric = m,
      n = nrow(xx),
      spearman_rho = unname(ct$estimate),
      p_value = ct$p.value
    )
  }))
}

write_csv(
  cor_out,
  file.path(OUT, "13_trajectory_distance_EII_spearman_exploratory.csv")
)

# ------------------------------------------------------------
# 13. QC
# ------------------------------------------------------------

qc <- tibble(
  selected_EII_file = selected_file,
  selected_EII_identifier_column = selected_id_col,
  selected_EII_value_column = eii_col,
  Step93Q_unique_runs = n_distinct(q_bridge$run_key),
  exact_run_overlap = selected_overlap,
  mapped_EII_rows = nrow(mapped_rows),
  mapped_EII_patients = n_distinct(mapped_rows$Patient_true),
  trajectory_patients = nrow(clusters2),
  final_matched_patients = nrow(matched),
  final_match_rate = nrow(matched) / nrow(clusters2)
)

write_csv(qc, file.path(OUT, "14_STEP93T10B_QC_SUMMARY.csv"))

if (nrow(matched) == 0) {
  writeLines(
    "Run bridge worked insufficiently for final patient integration. Inspect files 01-14.",
    file.path(OUT, "_STEP93T10B_FINAL_ZERO_MATCH.txt")
  )
  stop("STEP93T10B final patient match = 0. Diagnostics written.")
}

writeLines(
  c(
    paste0("Completed: ", Sys.time()),
    paste0("Selected EII: ", selected_file),
    paste0("Run overlap: ", selected_overlap),
    paste0("Mapped EII patients: ", n_distinct(mapped_rows$Patient_true)),
    paste0("Final matched trajectory patients: ", nrow(matched), "/", nrow(clusters2)),
    "STEP93T10B COMPLETE"
  ),
  file.path(OUT, "_STEP93T10B_COMPLETE.ok")
)

cat("STEP93T10B COMPLETE\n")
