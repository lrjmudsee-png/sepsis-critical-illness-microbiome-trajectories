#!/usr/bin/env Rscript

# A-stage feasibility audit for the mechanism/clinical upgrade direction.
# This script reads frozen analysis objects and registries only. It does not
# rerun sequence processing and does not estimate the proposed new effects.

options(stringsAsFactors = FALSE, warn = 1)

args <- commandArgs(trailingOnly = TRUE)
arg_value <- function(flag) {
  idx <- match(flag, args)
  if (is.na(idx) || idx == length(args)) stop("Missing required argument: ", flag)
  args[[idx + 1L]]
}
out_dir <- normalizePath(arg_value("--output"), winslash = "/", mustWork = TRUE)

write_csv_once <- function(x, filename) {
  path <- file.path(out_dir, filename)
  if (file.exists(path)) stop("Refusing to overwrite existing audit output: ", path)
  if (!is.data.frame(x) || nrow(x) == 0L) stop("Required output has no rows: ", filename)
  write.csv(x, path, row.names = FALSE, na = "", fileEncoding = "UTF-8")
}

write_lines_once <- function(x, filename) {
  path <- file.path(out_dir, filename)
  if (file.exists(path)) stop("Refusing to overwrite existing audit output: ", path)
  writeLines(enc2utf8(x), path, useBytes = TRUE)
}

read_csv <- function(path) {
  read.csv(path, check.names = FALSE, stringsAsFactors = FALSE, fileEncoding = "UTF-8-BOM")
}

as_flag <- function(x) {
  if (is.logical(x)) return(!is.na(x) & x)
  toupper(trimws(as.character(x))) %in% c("TRUE", "T", "1", "YES", "Y")
}

first_column <- function(df, candidates, default = NA_character_) {
  hit <- candidates[candidates %in% names(df)]
  if (!length(hit)) return(rep(default, nrow(df)))
  as.character(df[[hit[[1L]]]])
}

resolve_time <- function(df, target_labels, candidates) {
  resolved <- rep(NA_character_, nrow(df))
  for (column in candidates[candidates %in% names(df)]) {
    values <- as.character(df[[column]])
    take <- is.na(resolved) & !is.na(values) & values %in% target_labels
    resolved[take] <- values[take]
  }
  if (anyNA(resolved)) {
    fallback <- first_column(df, candidates)
    resolved[is.na(resolved)] <- fallback[is.na(resolved)]
  }
  resolved
}

normalize_family <- function(x) {
  value <- trimws(as.character(x))
  value[is.na(value) | value == ""] <- "<UNRESOLVED>"
  value <- sub("^[fF]__", "", value)
  value <- trimws(value)
  value[value == ""] <- "<UNRESOLVED>"
  value
}

family_group <- function(normalized_family) {
  result <- rep("OTHER", length(normalized_family))
  result[normalized_family %in% c("Enterobacteriaceae", "Enterococcaceae")] <- "P"
  result[normalized_family %in% c("Lachnospiraceae", "Ruminococcaceae", "Oscillospiraceae")] <- "C"
  result
}

collapse_values <- function(x) {
  x <- unique(as.character(x[!is.na(x) & as.character(x) != ""]))
  if (!length(x)) "" else paste(sort(x), collapse = "|")
}

safe_file_info <- function(path, category, role, notes = "") {
  info <- file.info(path)
  data.frame(
    category = category,
    input_role = role,
    path = normalizePath(path, winslash = "/", mustWork = FALSE),
    exists = file.exists(path),
    bytes = if (file.exists(path)) unname(info$size) else NA_real_,
    modified_time = if (file.exists(path)) format(info$mtime, "%Y-%m-%d %H:%M:%S %Z") else "",
    read_status = if (file.exists(path)) "READ_ONLY_INPUT" else "MISSING",
    notes = notes,
    stringsAsFactors = FALSE
  )
}

registry_path <- "E:/sepsis_project/results/V2_27B_ANALYSIS_OBJECTS_AND_REPLICATE_FREEZE/V2_STEP87B_analysis_object_registry.csv"
paired_registry_path <- "E:/sepsis_project/results/V2_UPGRADE_20260907/98A_AITCHISON_LONGITUDINAL_ROBUSTNESS/03_PAIRED_CONTRASTS/V2_98A_paired_early_vs_late_contrasts.csv"
meta_eligibility_path <- "E:/sepsis_project/results/V2_UPGRADE_20260907/98B_RANDOM_EFFECTS_META_ANALYSIS/META_ELIGIBILITY_FREEZE.csv"
external_population_path <- "E:/sepsis_project/results/V2_UPGRADE_20260910_FINAL_FREEZE/PRJNA1125274_FINAL_population_membership.csv"
external_primary_path <- "E:/sepsis_project/results/V2_UPGRADE_20260910_FINAL_FREEZE/PRJNA1125274_FINAL_primary_and_sensitivity_results.csv"
external_displacement_path <- "E:/sepsis_project/results/V2_UPGRADE_20260907/98E_PRJNA1125274_EXTERNAL_VALIDATION/PRJNA1125274_external_validation_sample_displacement.csv"
clinical_summary_path <- "E:/sepsis_project/results/V2_36G2_PRJNA851469_DAY3_CLINICAL_OUTCOME_MODELS_FIXED/10_MANUSCRIPT_SAFE_SUMMARY.txt"
eii_code_path <- "E:/sepsis_project/code/03_data_processing/92B2_EII_baseline_harmonization_fixed.R"
supplement_index_path <- file.path(out_dir, "A03_supplement_sheet_index.csv")
supplement_ids_path <- file.path(out_dir, "A03B_supplement_table2_sample_ids.csv")

required_inputs <- c(
  registry_path, paired_registry_path, meta_eligibility_path,
  external_population_path, external_primary_path, external_displacement_path,
  clinical_summary_path, eii_code_path, supplement_index_path, supplement_ids_path
)
missing_inputs <- required_inputs[!file.exists(required_inputs)]
if (length(missing_inputs)) stop("Missing required inputs:\n", paste(missing_inputs, collapse = "\n"))

registry <- read_csv(registry_path)
names(registry)[1L] <- "project"
selected_projects <- c("PRJNA691455", "PRJNA851469", "PRJNA516701")
registry <- registry[registry$project %in% selected_projects, , drop = FALSE]
if (!setequal(registry$project, selected_projects)) stop("Analysis object registry is missing one or more required projects")

configs <- data.frame(
  project = c("PRJNA691455", "PRJNA851469", "PRJNA516701", "PRJNA1125274"),
  anchor_label = c("Day 1", "Day-1", "DAY_1", "T0"),
  early_label = c("Day 3", "Day-3", "DAY_3", "T1"),
  late_label = c("Day 7", "Day-7", "DAY_7", "T2"),
  analysis_role = c("CORE_SEPSIS_LONGITUDINAL", "ICU_BACKGROUND_LONGITUDINAL", "ICU_BACKGROUND_LONGITUDINAL", "EXTERNAL_VALIDATION"),
  stringsAsFactors = FALSE
)
configs$object_path <- vapply(configs$project, function(project) {
  if (project == "PRJNA1125274") {
    "E:/sepsis_project/data/PRJNA1125274/03_dada2/PRJNA1125274_analysis_object_external_validation.rds"
  } else {
    registry$analysis_object_path[match(project, registry$project)]
  }
}, character(1L))

paired_registry <- read_csv(paired_registry_path)
paired_expected <- paired_registry[
  paired_registry$project %in% selected_projects &
    paired_registry$metric == "aitchison_from_patient_baseline_CZM" &
    paired_registry$subset_type == "ALL_PAIRED",
  c("project", "required_anchor", "early_label", "late_label", "n_pairs"),
  drop = FALSE
]
if (nrow(paired_expected) != 3L) stop("Could not resolve one frozen Aitchison pair definition per natural-history cohort")
for (project in selected_projects) {
  cfg <- configs[configs$project == project, , drop = FALSE]
  frozen <- paired_expected[paired_expected$project == project, , drop = FALSE]
  if (!identical(as.character(cfg$anchor_label), as.character(frozen$required_anchor)) ||
      !identical(as.character(cfg$early_label), as.character(frozen$early_label)) ||
      !identical(as.character(cfg$late_label), as.character(frozen$late_label))) {
    stop("Frozen label mismatch for ", project)
  }
}

external_population <- read_csv(external_population_path)
external_population$in_complete_case_sensitivity <- as_flag(external_population$in_complete_case_sensitivity)
external_population$in_prespecified_primary <- as_flag(external_population$in_prespecified_primary)

inventory_rows <- list(
  safe_file_info(registry_path, "registry", "analysis object registry"),
  safe_file_info(paired_registry_path, "frozen_results", "natural-cohort pair definition"),
  safe_file_info(meta_eligibility_path, "frozen_results", "meta-analysis eligibility freeze"),
  safe_file_info(external_population_path, "frozen_results", "external population membership"),
  safe_file_info(external_primary_path, "frozen_results", "external primary results read for overlap only"),
  safe_file_info(external_displacement_path, "frozen_results", "external sample displacement mapping"),
  safe_file_info(clinical_summary_path, "frozen_results", "existing clinical outcome result"),
  safe_file_info(eii_code_path, "code", "existing EII definition audit"),
  safe_file_info(supplement_index_path, "retrieved_source", "official supplement sheet index"),
  safe_file_info(supplement_ids_path, "retrieved_source", "official Table 2 sample identifiers")
)

sample_audits <- list()
pair_audits <- list()
family_audits <- list()
coverage_audits <- list()
population_audits <- list()
prjna851_metadata <- NULL

for (i in seq_len(nrow(configs))) {
  cfg <- configs[i, , drop = FALSE]
  project <- cfg$project
  object_path <- cfg$object_path
  if (!file.exists(object_path)) stop("Missing analysis object: ", object_path)
  message("Auditing ", project, " from ", object_path)
  obj <- readRDS(object_path)
  if (!all(c("counts", "taxonomy", "metadata") %in% names(obj))) stop("Object lacks required components: ", project)
  counts <- obj$counts
  taxonomy <- obj$taxonomy
  metadata <- obj$metadata
  if (!is.matrix(counts)) counts <- as.matrix(counts)
  if (is.null(rownames(counts)) || is.null(colnames(counts))) stop("Counts matrix lacks row/column keys: ", project)

  run_id <- if (project == "PRJNA1125274") {
    first_column(metadata, c("run_accession", "run_id"))
  } else {
    first_column(metadata, c("run_id", "run_accession"))
  }
  patient_id <- first_column(metadata, c("patient_id", "master_patient_id"))
  time_label <- resolve_time(
    metadata,
    c(cfg$anchor_label, cfg$early_label, cfg$late_label),
    if (project == "PRJNA1125274") c("timepoint", "time_label", "time_raw") else c("time_label", "time_raw", "master_time_raw", "timepoint_key")
  )
  hospital <- first_column(metadata, c("hospital", "center", "site"), "NOT_ENCODED")
  hospital[is.na(hospital) | hospital == ""] <- "NOT_ENCODED"
  primary_sample <- if ("primary_patient_time_include" %in% names(metadata)) {
    as_flag(metadata$primary_patient_time_include)
  } else {
    rep(TRUE, nrow(metadata))
  }

  feature_key <- if (project == "PRJNA1125274") "ASV_ID" else "ASV_sequence"
  if (!feature_key %in% names(taxonomy)) stop("Missing taxonomy feature key ", feature_key, " for ", project)
  if (!"Family" %in% names(taxonomy)) stop("Missing Family taxonomy for ", project)
  taxonomy_key <- as.character(taxonomy[[feature_key]])
  count_feature_key <- as.character(colnames(counts))
  tax_index <- match(count_feature_key, taxonomy_key)
  matched <- !is.na(tax_index)
  family_raw <- rep("<UNRESOLVED>", length(count_feature_key))
  family_raw[matched] <- as.character(taxonomy$Family[tax_index[matched]])
  family_raw[is.na(family_raw) | family_raw == ""] <- "<UNRESOLVED>"
  family_normalized <- normalize_family(family_raw)
  group <- family_group(family_normalized)

  p_feature <- group == "P"
  c_feature <- group == "C"
  sample_total <- rowSums(counts)
  sample_p <- if (any(p_feature)) rowSums(counts[, p_feature, drop = FALSE]) else rep(0, nrow(counts))
  sample_c <- if (any(c_feature)) rowSums(counts[, c_feature, drop = FALSE]) else rep(0, nrow(counts))
  count_run_id <- as.character(rownames(counts))
  count_index <- match(run_id, count_run_id)
  in_counts <- !is.na(count_index)
  p_by_meta <- rep(NA_real_, nrow(metadata)); p_by_meta[in_counts] <- sample_p[count_index[in_counts]]
  c_by_meta <- rep(NA_real_, nrow(metadata)); c_by_meta[in_counts] <- sample_c[count_index[in_counts]]
  total_by_meta <- rep(NA_real_, nrow(metadata)); total_by_meta[in_counts] <- sample_total[count_index[in_counts]]

  frozen_pair_patients <- character()
  in_complete_external <- rep(FALSE, nrow(metadata))
  in_primary_external <- rep(FALSE, nrow(metadata))
  if (project != "PRJNA1125274") {
    difference_path <- sprintf(
      "E:/sepsis_project/results/V2_UPGRADE_20260907/98A_AITCHISON_LONGITUDINAL_ROBUSTNESS/03_PAIRED_CONTRASTS/PATIENT_LEVEL_DIFFERENCES/%s_patient_paired_differences_aitchison_CZM.csv",
      project
    )
    if (!file.exists(difference_path)) stop("Missing patient-level frozen pair file: ", difference_path)
    difference_table <- read_csv(difference_path)
    frozen_pair_patients <- unique(as.character(difference_table$patient_id))
    inventory_rows[[length(inventory_rows) + 1L]] <- safe_file_info(
      difference_path, "frozen_results", paste0(project, " patient-level pair membership")
    )
  } else {
    population_index <- match(run_id, external_population$run_id)
    has_population <- !is.na(population_index)
    in_complete_external[has_population] <- external_population$in_complete_case_sensitivity[population_index[has_population]]
    in_primary_external[has_population] <- external_population$in_prespecified_primary[population_index[has_population]]
    frozen_pair_patients <- unique(external_population$patient_id[external_population$in_prespecified_primary])
  }

  patient_ids <- sort(unique(patient_id[!is.na(patient_id) & patient_id != ""]))
  project_pair_rows <- vector("list", length(patient_ids))
  for (j in seq_along(patient_ids)) {
    pid <- patient_ids[[j]]
    idx <- which(patient_id == pid & primary_sample)
    anchor_idx <- idx[!is.na(time_label[idx]) & time_label[idx] == cfg$anchor_label]
    early_idx <- idx[!is.na(time_label[idx]) & time_label[idx] == cfg$early_label]
    late_idx <- idx[!is.na(time_label[idx]) & time_label[idx] == cfg$late_label]
    exact_anchor <- length(anchor_idx) == 1L
    exact_pair <- length(early_idx) == 1L && length(late_idx) == 1L
    anchor_input_available <- exact_anchor && all(in_counts[anchor_idx])
    early_late_inputs_available <- exact_pair && all(in_counts[c(early_idx, late_idx)])
    # The frozen natural-history ALL_PAIRED estimand requires early+late only.
    # PRJNA516701 deliberately contains one ALL_PAIRED patient without strict
    # DAY_1 and handles that difference in the common-anchor sensitivity.
    # The external validation population is explicitly a complete T0/T1/T2 set.
    object_inputs_available <- early_late_inputs_available &&
      (project != "PRJNA1125274" || anchor_input_available)
    if (project == "PRJNA1125274") {
      frozen_primary <- pid %in% unique(external_population$patient_id[external_population$in_prespecified_primary])
      frozen_sensitivity <- pid %in% unique(external_population$patient_id[external_population$in_complete_case_sensitivity])
    } else {
      frozen_primary <- pid %in% frozen_pair_patients
      frozen_sensitivity <- frozen_primary
    }
    eligibility <- object_inputs_available && frozen_primary
    reason <- if (length(early_idx) != 1L) {
      paste0("early_count_", length(early_idx))
    } else if (length(late_idx) != 1L) {
      paste0("late_count_", length(late_idx))
    } else if (project == "PRJNA1125274" && length(anchor_idx) != 1L) {
      paste0("external_required_anchor_count_", length(anchor_idx))
    } else if (!object_inputs_available) {
      "one_or_more_runs_missing_from_counts"
    } else if (!frozen_primary) {
      if (frozen_sensitivity) "complete_case_sensitivity_only_not_prespecified_primary" else "not_in_frozen_pair_population"
    } else if (!anchor_input_available) {
      "eligible_frozen_ALL_PAIRED_without_strict_anchor_common_anchor_sensitivity_excludes"
    } else {
      "eligible_exact_frozen_pair"
    }
    runs_at <- function(indices) collapse_values(run_id[indices])
    p_at <- function(indices) if (length(indices) == 1L) p_by_meta[indices] else NA_real_
    c_at <- function(indices) if (length(indices) == 1L) c_by_meta[indices] else NA_real_
    any_component_zero <- exact_pair && any(c(
      p_at(early_idx) == 0, c_at(early_idx) == 0,
      p_at(late_idx) == 0, c_at(late_idx) == 0
    ), na.rm = TRUE)
    any_double_zero <- exact_pair && any(c(
      p_at(early_idx) + c_at(early_idx) == 0,
      p_at(late_idx) + c_at(late_idx) == 0
    ), na.rm = TRUE)
    project_pair_rows[[j]] <- data.frame(
      project = project,
      patient_id = pid,
      hospital = collapse_values(hospital[idx]),
      anchor_label = cfg$anchor_label,
      early_label = cfg$early_label,
      late_label = cfg$late_label,
      n_anchor_primary_samples = length(anchor_idx),
      n_early_primary_samples = length(early_idx),
      n_late_primary_samples = length(late_idx),
      anchor_run_id = runs_at(anchor_idx),
      early_run_id = runs_at(early_idx),
      late_run_id = runs_at(late_idx),
      anchor_input_available = anchor_input_available,
      early_late_inputs_available = early_late_inputs_available,
      object_inputs_available = object_inputs_available,
      in_frozen_primary_pair_population = frozen_primary,
      in_frozen_sensitivity_pair_population = frozen_sensitivity,
      A_stage_eligible = eligibility,
      early_P_reads = p_at(early_idx),
      early_C_reads = c_at(early_idx),
      late_P_reads = p_at(late_idx),
      late_C_reads = c_at(late_idx),
      any_P_or_C_zero_at_early_or_late = any_component_zero,
      any_P_plus_C_double_zero_at_early_or_late = any_double_zero,
      eligibility_reason = reason,
      stringsAsFactors = FALSE
    )
  }
  project_pairs <- do.call(rbind, project_pair_rows)
  pair_audits[[project]] <- project_pairs

  sample_role <- rep("other_time", length(time_label))
  sample_role[time_label == cfg$anchor_label] <- "anchor"
  sample_role[time_label == cfg$early_label] <- "early"
  sample_role[time_label == cfg$late_label] <- "late"
  meta_run_duplicates <- duplicated(run_id) | duplicated(run_id, fromLast = TRUE)
  count_run_duplicates <- duplicated(count_run_id) | duplicated(count_run_id, fromLast = TRUE)
  sample_audits[[project]] <- data.frame(
    project = project,
    patient_id = patient_id,
    hospital = hospital,
    run_id = run_id,
    resolved_time_label = time_label,
    prespecified_time_role = sample_role,
    primary_patient_time_include = primary_sample,
    run_present_exactly_in_counts = in_counts,
    metadata_run_id_duplicated = meta_run_duplicates,
    counts_run_id_duplicated = ifelse(in_counts, count_run_duplicates[count_index], NA),
    patient_in_frozen_pair_population = patient_id %in% frozen_pair_patients,
    in_external_complete_case_sensitivity = in_complete_external,
    in_external_prespecified_primary = in_primary_external,
    total_reads = total_by_meta,
    P_reads = p_by_meta,
    C_reads = c_by_meta,
    P_zero = !is.na(p_by_meta) & p_by_meta == 0,
    C_zero = !is.na(c_by_meta) & c_by_meta == 0,
    P_plus_C_double_zero = !is.na(p_by_meta) & !is.na(c_by_meta) & (p_by_meta + c_by_meta == 0),
    audit_status = ifelse(in_counts & !meta_run_duplicates, "EXACT_SAMPLE_COUNT_JOIN", "REVIEW_SAMPLE_COUNT_JOIN"),
    stringsAsFactors = FALSE
  )

  feature_abundance <- colSums(counts)
  family_df <- data.frame(
    original_family_label = family_raw,
    normalized_family_label = family_normalized,
    proposed_group = group,
    n_features = 1,
    total_reads = as.numeric(feature_abundance),
    stringsAsFactors = FALSE
  )
  family_df <- aggregate(
    cbind(n_features, total_reads) ~ original_family_label + normalized_family_label + proposed_group,
    data = family_df,
    FUN = sum
  )
  candidate_stem <- grepl(
    "enterobacter|enterococc|lachnospir|ruminococc|oscillospir",
    tolower(family_df$normalized_family_label)
  )
  family_df$mapping_evidence <- ifelse(
    family_df$proposed_group %in% c("P", "C"),
    "PREDECLARED_EXACT_LABEL_AFTER_ONLY_f_PREFIX_TRIM",
    "RETAINED_IN_OTHER_DENOMINATOR"
  )
  family_df$needs_manual_review <- family_df$normalized_family_label == "<UNRESOLVED>" |
    (candidate_stem & family_df$proposed_group == "OTHER")
  family_df$normalization_rule <- "trim whitespace; remove leading f__/F__ only; no synonym or result-driven remapping"
  family_df$project <- project
  family_df <- family_df[, c(
    "project", "original_family_label", "normalized_family_label", "proposed_group",
    "n_features", "total_reads", "mapping_evidence", "needs_manual_review", "normalization_rule"
  )]
  family_audits[[project]] <- family_df

  eligible_pairs <- project_pairs[project_pairs$A_stage_eligible, , drop = FALSE]
  coverage_audits[[project]] <- data.frame(
    project = project,
    n_samples = nrow(counts),
    n_features = ncol(counts),
    count_feature_key = if (project == "PRJNA1125274") "ASV_ID" else "ASV_sequence",
    taxonomy_feature_key = feature_key,
    n_feature_keys_matched = sum(matched),
    n_feature_keys_unmatched = sum(!matched),
    n_duplicated_taxonomy_keys = sum(duplicated(taxonomy_key)),
    n_P_features = sum(p_feature),
    n_C_features = sum(c_feature),
    n_OTHER_features = sum(group == "OTHER"),
    n_unresolved_family_features = sum(family_normalized == "<UNRESOLVED>"),
    unresolved_family_total_reads = sum(feature_abundance[family_normalized == "<UNRESOLVED>"]),
    unresolved_family_read_fraction = sum(feature_abundance[family_normalized == "<UNRESOLVED>"]) / sum(sample_total),
    P_total_reads = sum(sample_p),
    C_total_reads = sum(sample_c),
    OTHER_total_reads = sum(sample_total) - sum(sample_p) - sum(sample_c),
    P_read_fraction = sum(sample_p) / sum(sample_total),
    C_read_fraction = sum(sample_c) / sum(sample_total),
    n_samples_P_zero = sum(sample_p == 0),
    n_samples_C_zero = sum(sample_c == 0),
    n_samples_P_plus_C_double_zero = sum(sample_p + sample_c == 0),
    n_A_stage_eligible_pairs = nrow(eligible_pairs),
    n_eligible_pairs_with_any_P_or_C_zero = sum(eligible_pairs$any_P_or_C_zero_at_early_or_late),
    n_eligible_pairs_with_P_plus_C_double_zero = sum(eligible_pairs$any_P_plus_C_double_zero_at_early_or_late),
    n_eligible_pairs_fully_nonzero = sum(!eligible_pairs$any_P_or_C_zero_at_early_or_late),
    effect_estimation_performed = FALSE,
    stringsAsFactors = FALSE
  )

  sampled_counts <- counts[
    seq_len(min(nrow(counts), 20L)),
    seq_len(min(ncol(counts), 2000L)),
    drop = FALSE
  ]
  sampled_integer_like <- all(is.finite(sampled_counts)) && all(sampled_counts >= 0) &&
    all(abs(sampled_counts - round(sampled_counts)) < 1e-8)

  if (project == "PRJNA1125274") {
    overall <- data.frame(
      project = project,
      center = "ALL_HOSPITALS",
      n_samples_object = nrow(metadata),
      n_patients_object = length(unique(patient_id)),
      n_complete_case_patients = length(unique(external_population$patient_id[external_population$in_complete_case_sensitivity])),
      n_prespecified_primary_patients = length(unique(external_population$patient_id[external_population$in_prespecified_primary])),
      n_A_stage_eligible_pairs = sum(project_pairs$A_stage_eligible),
      stringsAsFactors = FALSE
    )
    by_center <- lapply(sort(unique(hospital)), function(center) {
      idx <- hospital == center
      pop_center <- external_population[external_population$hospital == center, , drop = FALSE]
      data.frame(
        project = project,
        center = center,
        n_samples_object = sum(idx),
        n_patients_object = length(unique(patient_id[idx])),
        n_complete_case_patients = length(unique(pop_center$patient_id[pop_center$in_complete_case_sensitivity])),
        n_prespecified_primary_patients = length(unique(pop_center$patient_id[pop_center$in_prespecified_primary])),
        n_A_stage_eligible_pairs = sum(project_pairs$A_stage_eligible & project_pairs$hospital == center),
        stringsAsFactors = FALSE
      )
    })
    population_audits[[project]] <- rbind(overall, do.call(rbind, by_center))
  } else {
    population_audits[[project]] <- data.frame(
      project = project,
      center = "STUDY_LEVEL_CENTER_NOT_ENCODED",
      n_samples_object = nrow(metadata),
      n_patients_object = length(unique(patient_id)),
      n_complete_case_patients = NA_integer_,
      n_prespecified_primary_patients = length(frozen_pair_patients),
      n_A_stage_eligible_pairs = sum(project_pairs$A_stage_eligible),
      stringsAsFactors = FALSE
    )
  }

  inventory_rows[[length(inventory_rows) + 1L]] <- safe_file_info(
    object_path, "analysis_object", paste0(project, " counts+taxonomy+metadata"),
    sprintf(
      "components=%s; samples=%d; features=%d; metadata_rows=%d; taxonomy_rows=%d; integer_like_counts=%s",
      paste(names(obj), collapse = "|"), nrow(counts), ncol(counts), nrow(metadata), nrow(taxonomy),
      sampled_integer_like
    )
  )
  if (project == "PRJNA851469") prjna851_metadata <- metadata
  rm(obj, counts, taxonomy, metadata, feature_abundance)
  invisible(gc())
}

input_inventory <- do.call(rbind, inventory_rows)
sample_audit <- do.call(rbind, sample_audits)
pair_audit <- do.call(rbind, pair_audits)
population_audit <- do.call(rbind, population_audits)
family_audit <- do.call(rbind, family_audits)
coverage_audit <- do.call(rbind, coverage_audits)

write_csv_once(input_inventory, "A01_input_inventory.csv")
write_csv_once(sample_audit, "A04_patient_time_join_audit.csv")
write_csv_once(pair_audit, "A05_pair_eligibility.csv")
write_csv_once(population_audit, "A06_population_counts.csv")
write_csv_once(family_audit, "A07_taxonomy_mapping_proposal.csv")
write_csv_once(coverage_audit, "A08_feature_and_zero_coverage.csv")

supplement_ids <- read_csv(supplement_ids_path)
healthy_supplement <- as_flag(supplement_ids$healthy_like)
patient_supplement <- as_flag(supplement_ids$patient_like)
candidate_columns <- intersect(
  c("source_sample_id", "sample_id", "patient_id", "phenotype", "control_type", "sample_class", "analysis_role"),
  names(prjna851_metadata)
)
metadata_text <- apply(prjna851_metadata[, candidate_columns, drop = FALSE], 1L, function(row) paste(row, collapse = "|"))
healthy_object <- grepl("healthy|volunteer|non[-_ ]?icu|control", metadata_text, ignore.case = TRUE)
object_source_ids <- first_column(prjna851_metadata, c("source_sample_id", "sample_id"))
matched_ids <- intersect(unique(object_source_ids), unique(as.character(supplement_ids$sample_id)))
supplement_only_ids <- setdiff(unique(as.character(supplement_ids$sample_id)), unique(object_source_ids))
supplement_healthy_ids <- unique(as.character(supplement_ids$sample_id[healthy_supplement]))
supplement_only_healthy_ids <- intersect(supplement_only_ids, supplement_healthy_ids)
healthy_reference_status <- if (sum(healthy_object) > 0L) {
  "FEASIBLE_FROM_13_FROZEN_HEALTHY_SAMPLES_WITH_OFFICIAL_TABLE2_SOURCE_CHECK"
} else if (sum(healthy_supplement) > 0L) {
  "CONDITIONALLY_FEASIBLE_FROM_OFFICIAL_TABLE2_ONLY"
} else {
  "NO_HEALTHY_SAMPLE_IDS_DETECTED"
}
healthy_audit <- data.frame(
  project = "PRJNA851469",
  frozen_object_samples = nrow(prjna851_metadata),
  frozen_object_healthy_like_samples = sum(healthy_object),
  official_table2_unique_sample_ids = nrow(supplement_ids),
  official_table2_healthy_like_sample_ids = sum(healthy_supplement),
  official_table2_patient_like_sample_ids = sum(patient_supplement),
  exact_sample_id_overlap_object_vs_table2 = length(matched_ids),
  exact_overlap_examples = paste(head(sort(matched_ids), 12L), collapse = "|"),
  official_table2_ids_not_in_frozen_object = length(supplement_only_ids),
  official_table2_healthy_ids_not_in_frozen_object = length(supplement_only_healthy_ids),
  official_table2_healthy_ids_not_in_frozen_object_examples = paste(sort(supplement_only_healthy_ids), collapse = "|"),
  feasibility_status = healthy_reference_status,
  independence_status = "SAME_MICRO_ICU_COHORT_NOT_INDEPENDENT_VALIDATION",
  batch_processing_warning = "Use the 13 healthy samples already co-processed in the frozen RDS for the main healthy reference. Use official Table 2 only to audit the 15-source-sample roster; do not naively merge abundance scales.",
  permitted_A_stage_conclusion = "A healthy-reference analysis is technically feasible from 13 co-processed frozen samples. It remains a secondary same-cohort reference, not independent external validation; two official healthy samples are absent from the frozen RDS.",
  stringsAsFactors = FALSE
)
write_csv_once(healthy_audit, "A09_healthy_reference_feasibility.csv")

coverage_lines <- vapply(seq_len(nrow(coverage_audit)), function(i) {
  row <- coverage_audit[i, ]
  sprintf(
    "- %s: taxonomy keys matched %d/%d; P/C features %d/%d; A-stage eligible pairs %d; pairs with any P or C zero %d.",
    row$project, row$n_feature_keys_matched, row$n_features, row$n_P_features, row$n_C_features,
    row$n_A_stage_eligible_pairs, row$n_eligible_pairs_with_any_P_or_C_zero
  )
}, character(1L))

overlap_md <- c(
  "# A11 Existing-analysis overlap audit",
  "",
  "## Conclusion",
  "",
  "The proposed P/C ecological balance can proceed as a new, predeclared aggregate metric. It is nevertheless derived from the same 16S counts and is not an independent mechanistic validation. Relations with Bray, Aitchison, alpha diversity, ASV instability, and EII must be described as consistency or interpretive analyses rather than additional independent evidence.",
  "",
  "## Confirmed overlap",
  "",
  "- The existing EII combines Bray displacement, Shannon/Simpson, and ASV instability. A P/C--EII association therefore has structural overlap in its data source.",
  "- Step89A2 already tested genus- and family-level CLR trajectories. P/C is a meaningful ecological-axis compression only if its family dictionary, pseudocount sensitivity, and single primary contrast are fixed before effect estimation.",
  "- PRJNA1125274 remains external validation of a microbiome trajectory. It has no direct host molecular readout and cannot be upgraded to host-mechanistic validation.",
  "- The existing PRJNA851469 Day-3 clinical outcome primary model is negative: landmark n=42, events=30, HR=0.929, 95% CI 0.624-1.381, p=0.7145; after adjustment for admission antibiotics HR=1.155, p=0.5087. This negative evidence must be retained.",
  "",
  "## Reporting governance",
  "",
  "- Permitted: present P/C as a predeclared candidate ecological axis assessed for cross-cohort directional consistency.",
  "- Not permitted: call a Bray/EII correlation independent validation, or call 16S functional inference direct functional evidence.",
  "- Not permitted: modify the P/C family list because of coverage or future significance.",
  "",
  "## Audited sources",
  "",
  paste0("- EII code: `", eii_code_path, "`"),
  paste0("- Existing clinical-outcome summary: `", clinical_summary_path, "`"),
  "- Official PRJNA851469 supplement: A02/A03."
)
write_lines_once(overlap_md, "A11_existing_analysis_overlap.md")

all_taxonomy_exact <- all(coverage_audit$n_feature_keys_unmatched == 0L) &&
  all(coverage_audit$n_duplicated_taxonomy_keys == 0L)
all_have_pc <- all(coverage_audit$n_P_features > 0L & coverage_audit$n_C_features > 0L)
pair_counts_text <- paste(
  sprintf("%s=%d", coverage_audit$project, coverage_audit$n_A_stage_eligible_pairs),
  collapse = "; "
)
decision <- if (all_taxonomy_exact && all_have_pc) "GO_WITH_GOVERNANCE" else "HOLD_FOR_DATA_REPAIR"
decision_md <- c(
  "# A12 Feasibility decision",
  "",
  paste0("**Overall decision: ", decision, ".**"),
  "",
  "## Components that may enter stage B",
  "",
  paste0("- P/C balance: ", if (all_taxonomy_exact && all_have_pc) "proceed" else "hold", ". All four objects have strict feature-key joins and the candidate families were not selected from results."),
  paste0("- Frozen pair populations: ", pair_counts_text, ". The external-validation number refers only to 24 prespecified date-clean patients; 30 complete cases remain a sensitivity population."),
  "- Transition architecture/stage phenotypes: may proceed as a secondary exploratory module with multiplicity control, separate from the single primary contrast.",
  "- Healthy-reference resilience: feasible from 13 healthy samples already co-processed in the frozen PRJNA851469 RDS. Official Table 2 lists 15 healthy samples, two of which are absent from the frozen object. This is a same-cohort secondary reference, not independent validation.",
  "",
  "## Components that cannot currently enter",
  "",
  "- Patient-level host/immune mechanism: no public host measurements with a strict sample join were found. Stage C remains on hold unless the corresponding author supplies de-identified host data, sampling times, and linkage keys.",
  "- Causal mechanistic language: current evidence supports ecological mechanistic consistency and clinical association at most, not a causal host mechanism.",
  "",
  "## Governance conditions to freeze before stage B",
  "",
  "1. P={Enterobacteriaceae, Enterococcaceae}; C={Lachnospiraceae, Ruminococcaceae, Oscillospiraceae}. Remove a leading f__ from Family labels only; do not add synonym mappings.",
  "2. Use frozen pairs and fixed early/late times. PRJNA1125274 primary analysis is limited to 24 patients; 30 patients form a sensitivity analysis.",
  "3. Report effect size, 95% CI, p value, zero/pseudocount sensitivity, and all negative results.",
  "4. Label Bray/EII analyses as same-source consistency analyses, not independent validation.",
  "5. Do not alter taxonomy, population, or time contrasts in response to results.",
  "6. PRJNA516701 requires explicit zero governance: 8/15 frozen pairs contain at least one P or C zero and 3/15 contain a P+C double-zero sample. Follow the candidate SAP: add 0.5 only after group aggregation; use 0.1 and 1 as sensitivities; treat double-zero samples as uninterpretable missing values rather than log(0.5/0.5)=0; disclose the resulting pair attrition.",
  "7. Unresolved Family labels remain in OTHER and in the total-read denominator. This affects about 23.4% of PRJNA1125274 reads (versus about 1-2% in the three natural-history cohorts) and must be disclosed as an external-cohort taxonomy-resolution limitation.",
  "",
  "## Coverage snapshot",
  "",
  coverage_lines,
  "",
  "Stage A did not calculate any new delta_B effects, confidence intervals, or p values."
)
write_lines_once(decision_md, "A12_feasibility_decision.md")

questions_md <- c(
  "# A13 Questions requiring scientific approval",
  "",
  "These questions do not block stage A, but they should be approved in writing before stage B:",
  "",
  "1. Approve the fixed P/C family dictionary, including both Ruminococcaceae and Oscillospiraceae in C, with no result-driven revision?",
  "2. Freeze the only primary stage-B contrast as early-to-late delta_B, with transition types, Bray/EII relations, and healthy-reference analyses all secondary/exploratory?",
  "3. Keep PRJNA1125274 strictly at 24 date-clean patients for primary analysis and 30 complete cases for sensitivity?",
  "4. Accept the 13-sample healthy reference as a secondary within-MICRO-ICU reference, not independent validation?",
  "5. Authorize a request to the PRJNA851469/MICRO-ICU corresponding author for de-identified patient-level antibody/host measurements, sampling times, and linkage keys? Stage C should remain paused until those data are obtained.",
  "6. Keep the manuscript's primary narrative as cross-cohort ecological destabilization plus external validation, with P/C balance as the mechanistic-deepening secondary line? This is recommended.",
  "7. Approve a hierarchical multiplicity strategy: one primary P/C contrast, with remaining modules controlled by module-level FDR and labelled exploratory?"
)
write_lines_once(questions_md, "A13_questions_for_scientific_review.md")

session_path <- file.path(out_dir, "session_info.txt")
if (file.exists(session_path)) stop("Refusing to overwrite existing audit output: ", session_path)
session_lines <- c(
  paste0("audit_timestamp=", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
  "script=E:/sepsis_project/code/03_data_processing/99_DIRECTION_UPGRADE/99A2_direction_feasibility_audit.R",
  "scope=A-stage feasibility only; no new effect estimates; no DADA2/FASTQ/98D rerun",
  capture.output(sessionInfo())
)
writeLines(session_lines, session_path, useBytes = TRUE)

message("A-stage R audit completed. Decision: ", decision)
