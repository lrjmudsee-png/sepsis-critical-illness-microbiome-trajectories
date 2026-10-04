# ============================================================
# Sepsis V2 - Step93T10
# FINAL EII + TAXONOMIC TRAJECTORY INTEGRATION
#
# Fixes the root problem:
# Previous scripts used a generic "id" detector and could select run_id
# instead of the true patient_id from the EII table.
#
# This version:
# 1) locks analysis to PRJNA691455
# 2) explicitly prioritizes true patient_id columns
# 3) uses the frozen PRJNA691455 manifest as a Run -> Patient fallback
# 4) integrates trajectory clusters + signed trajectory distance + EII
# 5) performs exploratory cluster/EII and distance/EII statistics
#
# IMPORTANT:
# n is small (core sepsis longitudinal cohort); interpret statistics
# as exploratory/supportive, not definitive validation.
# ============================================================

options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(stringr)
})

ROOT <- "E:/sepsis_project"
PROJECT <- "PRJNA691455"

OUT <- file.path(
  ROOT,
  "results",
  "V2_33T10_EII_TRAJECTORY_FINAL_INTEGRATION"
)

dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

checkpoint <- function(txt) {
  cat(txt, "\n")
  cat(txt, "\n",
      file = file.path(OUT, "_STEP93T10_runtime_checkpoints.txt"),
      append = TRUE)
}

checkpoint("STEP93T10 START")


# ------------------------------------------------------------
# Helper functions
# ------------------------------------------------------------

pick_exact <- function(nms, priority) {
  hit <- priority[priority %in% nms]
  if (length(hit) == 0) return(NA_character_)
  hit[1]
}

normalize_patient <- function(x) {
  x <- as.character(x)
  x <- str_trim(x)
  x <- toupper(x)

  # remove explicit project prefix if present
  x <- str_replace(
    x,
    regex("^PRJNA691455[_:\\-]*", ignore_case = TRUE),
    ""
  )

  # remove sample suffix if a sample-level ID leaked in
  x <- str_replace(
    x,
    regex("FFO.*$", ignore_case = TRUE),
    ""
  )

  x[x %in% c("", "NA", "N/A", "NULL")] <- NA_character_
  x
}

safe_num <- function(x) {
  suppressWarnings(as.numeric(x))
}


# ------------------------------------------------------------
# 1. Load trajectory clusters from the corrected Step93S2
# ------------------------------------------------------------

cluster_file <- file.path(
  ROOT,
  "results",
  "V2_33S2_TRAJECTORY_DISTANCE_METHOD_FIX",
  "trajectory_clusters.csv"
)

if (!file.exists(cluster_file)) {
  stop("Missing corrected Step93S2 trajectory_clusters.csv")
}

clusters <- read_csv(cluster_file, show_col_types = FALSE)

cluster_patient_col <- pick_exact(
  names(clusters),
  c("Patient_true", "patient_true", "Patient", "patient_id")
)

if (is.na(cluster_patient_col)) {
  stop("Cannot identify patient column in trajectory_clusters.csv")
}

clusters2 <- clusters %>%
  transmute(
    trajectory_patient_raw = as.character(.data[[cluster_patient_col]]),
    patient_key = normalize_patient(.data[[cluster_patient_col]]),
    cluster = as.integer(cluster)
  ) %>%
  distinct(patient_key, .keep_all = TRUE)

write_csv(
  clusters2,
  file.path(OUT, "01_trajectory_cluster_patient_ids.csv")
)

checkpoint(
  paste0("Trajectory patients loaded: ", nrow(clusters2))
)


# ------------------------------------------------------------
# 2. Load signed trajectory distance matrix
# ------------------------------------------------------------

distance_file <- file.path(
  ROOT,
  "results",
  "V2_33S2_TRAJECTORY_DISTANCE_METHOD_FIX",
  "signed_trajectory_euclidean_distance_matrix.csv"
)

distance_centrality <- NULL

if (file.exists(distance_file)) {

  dm <- read.csv(
    distance_file,
    row.names = 1,
    check.names = FALSE
  )

  dm <- as.matrix(dm)
  storage.mode(dm) <- "numeric"

  distance_centrality <- tibble(
    patient_key = normalize_patient(rownames(dm)),
    mean_signed_trajectory_distance = rowMeans(dm, na.rm = TRUE),
    median_signed_trajectory_distance = apply(
      dm, 1, median, na.rm = TRUE
    )
  )

  write_csv(
    distance_centrality,
    file.path(OUT, "02_patient_signed_trajectory_distance_summary.csv")
  )
}


# ------------------------------------------------------------
# 3. Locate the REAL Step92B patient-timepoint EII table
# ------------------------------------------------------------

preferred_eii <- file.path(
  ROOT,
  "results",
  "V2_32B_STEP92B_EII_CONSTRUCTION",
  "V2_STEP92B_patient_timepoint_EII.csv"
)

if (file.exists(preferred_eii)) {
  eii_file <- preferred_eii
} else {
  found <- list.files(
    file.path(ROOT, "results"),
    recursive = TRUE,
    full.names = TRUE,
    pattern = "^V2_STEP92B_patient_timepoint_EII\\.csv$"
  )

  if (length(found) == 0) {
    stop("Cannot find V2_STEP92B_patient_timepoint_EII.csv")
  }

  eii_file <- found[1]
}

eii <- read_csv(
  eii_file,
  show_col_types = FALSE,
  guess_max = 10000
)

write_csv(
  tibble(
    selected_eii_file = eii_file,
    n_rows = nrow(eii),
    n_cols = ncol(eii),
    columns = paste(names(eii), collapse = ";")
  ),
  file.path(OUT, "03_EII_input_audit.csv")
)

checkpoint(
  paste0("EII rows loaded: ", nrow(eii))
)


# ------------------------------------------------------------
# 4. Explicit project filtering
# ------------------------------------------------------------

project_col <- pick_exact(
  names(eii),
  c("project", "Project", "PROJECT")
)

if (!is.na(project_col)) {

  eii_project <- eii %>%
    filter(
      toupper(as.character(.data[[project_col]])) == PROJECT
    )

} else {

  # If project column absent, use frozen manifest runs to identify project.
  eii_project <- eii
}

checkpoint(
  paste0("EII rows after PRJNA691455 project filter: ", nrow(eii_project))
)


# ------------------------------------------------------------
# 5. Pick TRUE patient ID column -- NEVER generic "id"
# ------------------------------------------------------------

patient_priority <- c(
  "patient_id",
  "master_patient_id",
  "Patient_ID",
  "patient_uid",
  "Patient_true",
  "subject_id",
  "Subject_ID",
  "patient"
)

patient_col <- pick_exact(
  names(eii_project),
  patient_priority
)

run_col <- pick_exact(
  names(eii_project),
  c("run_id", "Run", "run", "run_uid", "Run_ID")
)

mapping_method <- NA_character_


# ------------------------------------------------------------
# 6. Prefer direct patient_id; otherwise bridge via frozen manifest
# ------------------------------------------------------------

if (!is.na(patient_col)) {

  eii_project$patient_key <- normalize_patient(
    eii_project[[patient_col]]
  )

  mapping_method <- paste0(
    "DIRECT_EII_PATIENT_COLUMN:",
    patient_col
  )

} else {

  if (is.na(run_col)) {
    stop(
      "EII table has neither a recognized patient_id nor run_id column."
    )
  }

  manifest_file <- file.path(
    ROOT,
    "data",
    "_V2_ANALYSIS_READY",
    "01_PROJECTS",
    PROJECT,
    "00_metadata",
    "PRJNA691455_frozen_16S_manifest.csv"
  )

  if (!file.exists(manifest_file)) {
    stop("PRJNA691455 frozen manifest not found.")
  }

  manifest <- read_csv(
    manifest_file,
    show_col_types = FALSE
  )

  manifest_run_col <- pick_exact(
    names(manifest),
    c("run_id", "Run", "run", "run_uid", "Run_ID")
  )

  manifest_patient_col <- pick_exact(
    names(manifest),
    c("patient_id", "master_patient_id", "Patient_ID", "patient_uid")
  )

  if (is.na(manifest_run_col) || is.na(manifest_patient_col)) {
    stop("Frozen manifest lacks run_id or patient_id.")
  }

  bridge <- manifest %>%
    transmute(
      run_key = as.character(.data[[manifest_run_col]]),
      patient_key = normalize_patient(
        .data[[manifest_patient_col]]
      )
    ) %>%
    distinct()

  eii_project <- eii_project %>%
    mutate(
      run_key = as.character(.data[[run_col]])
    ) %>%
    left_join(
      bridge,
      by = "run_key"
    )

  mapping_method <- "FROZEN_MANIFEST_RUN_TO_PATIENT"
}

writeLines(
  c(
    paste0("EII file: ", eii_file),
    paste0("Project column: ", project_col),
    paste0("Patient mapping method: ", mapping_method),
    paste0("Detected patient column: ", patient_col),
    paste0("Detected run column: ", run_col)
  ),
  file.path(OUT, "04_EII_patient_mapping_method.txt")
)


# ------------------------------------------------------------
# 7. Identify the actual EII numeric column
# ------------------------------------------------------------

eii_priority <- c(
  "EII_0_100",
  "EII",
  "eii",
  "EII_raw",
  "ecological_instability_index",
  "instability_index"
)

eii_col <- pick_exact(
  names(eii_project),
  eii_priority
)

if (is.na(eii_col)) {

  eii_like <- names(eii_project)[
    str_detect(
      names(eii_project),
      regex("EII|instability", ignore_case = TRUE)
    )
  ]

  numeric_like <- eii_like[
    vapply(
      eii_project[eii_like],
      is.numeric,
      logical(1)
    )
  ]

  if (length(numeric_like) == 0) {
    stop(
      paste0(
        "No EII numeric column found. Available columns: ",
        paste(names(eii_project), collapse = ", ")
      )
    )
  }

  eii_col <- numeric_like[1]
}

eii_project$EII_value <- safe_num(
  eii_project[[eii_col]]
)

writeLines(
  paste0("Selected EII value column: ", eii_col),
  file.path(OUT, "05_selected_EII_column.txt")
)


# ------------------------------------------------------------
# 8. Determine longitudinal order
# ------------------------------------------------------------

time_col <- pick_exact(
  names(eii_project),
  c(
    "analysis_time_order",
    "time_order",
    "time_day",
    "Time_numeric",
    "time_numeric",
    "day"
  )
)

if (!is.na(time_col)) {

  eii_project$time_order_final <- safe_num(
    eii_project[[time_col]]
  )

} else {

  eii_project <- eii_project %>%
    group_by(patient_key) %>%
    mutate(
      time_order_final = row_number()
    ) %>%
    ungroup()
}


# ------------------------------------------------------------
# 9. Patient-level EII summaries
# ------------------------------------------------------------

eii_clean <- eii_project %>%
  filter(
    !is.na(patient_key),
    !is.na(EII_value)
  )

patient_eii <- eii_clean %>%
  arrange(
    patient_key,
    time_order_final
  ) %>%
  group_by(patient_key) %>%
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
  file.path(OUT, "06_PRJNA691455_patient_level_EII.csv")
)


# ------------------------------------------------------------
# 10. Join EII to trajectory patients
# ------------------------------------------------------------

integrated <- clusters2 %>%
  left_join(
    patient_eii,
    by = "patient_key"
  )

if (!is.null(distance_centrality)) {

  integrated <- integrated %>%
    left_join(
      distance_centrality,
      by = "patient_key"
    )
}

write_csv(
  integrated,
  file.path(OUT, "07_patient_EII_trajectory_integrated.csv")
)


matched <- integrated %>%
  filter(!is.na(EII_mean))

unmatched <- integrated %>%
  filter(is.na(EII_mean))

write_csv(
  matched,
  file.path(OUT, "08_matched_trajectory_EII_patients.csv")
)

write_csv(
  unmatched,
  file.path(OUT, "09_unmatched_trajectory_patients.csv")
)

write_csv(
  anti_join(
    patient_eii,
    clusters2,
    by = "patient_key"
  ),
  file.path(OUT, "10_EII_patients_not_in_trajectory_set.csv")
)

checkpoint(
  paste0(
    "Matched trajectory patients to EII: ",
    nrow(matched),
    " / ",
    nrow(clusters2)
  )
)


# ------------------------------------------------------------
# 11. Hard sanity guard
# ------------------------------------------------------------

if (nrow(matched) == 0) {

  writeLines(
    paste0(
      "No patients matched. Inspect files 03-10 in: ",
      OUT
    ),
    file.path(OUT, "_STEP93T10_NO_MATCH_ERROR.txt")
  )

  stop(
    "STEP93T10 found 0 EII/trajectory patient matches. Diagnostics were written."
  )
}


# ------------------------------------------------------------
# 12. Cluster vs EII exploratory comparisons
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
  file.path(OUT, "11_cluster_EII_descriptive_summary.csv")
)


wilcox_results <- tibble()

if (
  n_distinct(matched$cluster) == 2 &&
  all(table(matched$cluster) >= 2)
) {

  metrics <- c("EII_change", "EII_mean", "EII_max")

  wilcox_results <- bind_rows(
    lapply(
      metrics,
      function(metric) {

        x <- matched %>%
          filter(!is.na(.data[[metric]]))

        wt <- tryCatch(
          wilcox.test(
            x[[metric]] ~ x$cluster,
            exact = FALSE
          ),
          error = function(e) NULL
        )

        tibble(
          metric = metric,
          n = nrow(x),
          p_value = if (is.null(wt)) NA_real_ else wt$p.value
        )
      }
    )
  )
}

write_csv(
  wilcox_results,
  file.path(OUT, "12_cluster_EII_wilcoxon_results.csv")
)


# ------------------------------------------------------------
# 13. Trajectory divergence vs EII correlations
# ------------------------------------------------------------

cor_results <- tibble()

if (
  "mean_signed_trajectory_distance" %in% names(matched)
) {

  metrics <- c("EII_change", "EII_mean", "EII_max")

  cor_results <- bind_rows(
    lapply(
      metrics,
      function(metric) {

        x <- matched %>%
          select(
            mean_signed_trajectory_distance,
            all_of(metric)
          ) %>%
          filter(
            complete.cases(.)
          )

        if (nrow(x) < 3) {
          return(
            tibble(
              metric = metric,
              n = nrow(x),
              spearman_rho = NA_real_,
              p_value = NA_real_
            )
          )
        }

        ct <- suppressWarnings(
          cor.test(
            x$mean_signed_trajectory_distance,
            x[[metric]],
            method = "spearman",
            exact = FALSE
          )
        )

        tibble(
          metric = metric,
          n = nrow(x),
          spearman_rho = unname(ct$estimate),
          p_value = ct$p.value
        )
      }
    )
  )
}

write_csv(
  cor_results,
  file.path(OUT, "13_trajectory_distance_EII_spearman_results.csv")
)


# ------------------------------------------------------------
# 14. Final QC summary
# ------------------------------------------------------------

qc <- tibble(
  project = PROJECT,
  trajectory_patients = nrow(clusters2),
  EII_patients_project = nrow(patient_eii),
  matched_patients = nrow(matched),
  match_rate = nrow(matched) / nrow(clusters2),
  mapping_method = mapping_method,
  EII_column = eii_col,
  time_column = ifelse(is.na(time_col), "ROW_ORDER_FALLBACK", time_col)
)

write_csv(
  qc,
  file.path(OUT, "14_STEP93T10_QC_SUMMARY.csv")
)


writeLines(
  c(
    paste0("Completed: ", Sys.time()),
    paste0("Project: ", PROJECT),
    paste0("Matched patients: ", nrow(matched), "/", nrow(clusters2)),
    paste0("Mapping method: ", mapping_method),
    paste0("EII column: ", eii_col),
    "STEP93T10 COMPLETE"
  ),
  file.path(OUT, "_STEP93T10_COMPLETE.ok")
)

cat("STEP93T10 COMPLETE\n")
