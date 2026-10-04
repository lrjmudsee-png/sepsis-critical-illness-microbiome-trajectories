# ============================================================
# Sepsis V2 - Step 87B
# Freeze same-patient/same-time replicate rule
# + build cohort-specific portable analysis objects
#
# Current frozen input state (Step87A2):
#   786 sequence-QC samples
#   8 projects
#   390 project-specific patients
#   258 patients with >=2 timepoints
#   132 patients with >=3 timepoints
#
# Known duplicate:
#   PRJNA516701 / patient 4200 / DAY_5
#   SRR8488668 = 45,312 target reads
#   SRR8488665 = 10,823 target reads
#
# Primary rule:
#   - Do NOT delete either canonical sample.
#   - Do NOT merge their count vectors.
#   - For models requiring one observation per patient-timepoint,
#     retain the higher-depth sample as the primary representative.
#   - Hold the other sample for sensitivity analysis.
#
# R 4.4.x / Windows
# ============================================================

options(stringsAsFactors = FALSE)

pkgs <- c("readr", "dplyr", "stringr", "tibble", "purrr")
missing_pkgs <- pkgs[
  !vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_pkgs)) {
  install.packages(
    missing_pkgs,
    repos = "https://cloud.r-project.org"
  )
}

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(stringr)
  library(tibble)
  library(purrr)
})

# ------------------------------------------------------------
# Paths
# ------------------------------------------------------------

ROOT <- "E:/sepsis_project"

CANON <- file.path(
  ROOT,
  "data",
  "_V2_ANALYSIS_READY"
)

PROJECT_ROOT <- file.path(
  CANON,
  "01_PROJECTS"
)

STEP86B <- file.path(
  ROOT,
  "results",
  "V2_26B_DEPTH2000_FINAL_ANALYSIS_READY"
)

STEP87A2 <- file.path(
  ROOT,
  "results",
  "V2_27A2_ANALYSIS_METADATA_REBUILD_AND_REPLICATE_AUDIT_FIXED"
)

OUT <- file.path(
  ROOT,
  "results",
  "V2_27B_ANALYSIS_OBJECTS_AND_REPLICATE_FREEZE"
)

dir.create(
  OUT,
  recursive = TRUE,
  showWarnings = FALSE
)

LOG <- file.path(
  OUT,
  "_STEP87B_runtime_checkpoints.txt"
)

ERR <- file.path(
  OUT,
  "_STEP87B_FATAL_ERROR.txt"
)

if (file.exists(ERR)) {
  unlink(ERR)
}

# ------------------------------------------------------------
# Helpers
# ------------------------------------------------------------

ck <- function(x) {
  cat(
    paste0(x, ": ", Sys.time(), "\n"),
    file = LOG,
    append = TRUE
  )
}

safe_csv <- function(p) {
  suppressMessages(
    read_csv(
      p,
      show_col_types = FALSE,
      progress = FALSE,
      name_repair = "unique"
    )
  )
}

atomic_saveRDS <- function(obj, p) {
  tmp <- paste0(p, ".tmp")

  if (file.exists(tmp)) {
    unlink(tmp)
  }

  saveRDS(obj, tmp)

  if (file.exists(p)) {
    unlink(p)
  }

  if (!file.rename(tmp, p)) {
    stop(
      "Could not atomically save RDS: ",
      p,
      "\nPossible causes: file lock, antivirus, OneDrive sync, or permission problem."
    )
  }
}

taxonomy_path <- function(project) {
  p <- file.path(
    PROJECT_ROOT,
    project,
    "04_taxonomy_qc",
    paste0(
      project,
      "_taxonomy_SILVA1382_minBoot80.csv"
    )
  )

  if (!file.exists(p)) {
    stop(
      project,
      ": taxonomy CSV missing: ",
      p
    )
  }

  p
}

detect_asv_sequence_col <- function(tax) {
  candidates <- c(
    "ASV_sequence",
    "asv_sequence",
    "sequence",
    "Sequence"
  )

  hit <- candidates[candidates %in% names(tax)]

  if (!length(hit)) {
    stop(
      "Taxonomy table has no ASV sequence column. Columns are: ",
      paste(names(tax), collapse = ", ")
    )
  }

  hit[1]
}

patient_time_stats <- function(d) {
  pt <- d |>
    distinct(
      patient_id,
      timepoint_key
    ) |>
    count(
      patient_id,
      name = "n_timepoints"
    )

  tibble(
    patients = n_distinct(d$patient_id),
    patients_GE2 = sum(pt$n_timepoints >= 2),
    patients_GE3 = sum(pt$n_timepoints >= 3)
  )
}

role_registry <- tribble(
  ~project, ~analysis_role,

  "PRJNA691455",
  "CORE_SEPSIS_LONGITUDINAL",

  "PRJEB82425",
  "ICU_INFECTION_LONGITUDINAL_EXTERNAL",

  "PRJNA516701",
  "ICU_BACKGROUND_LONGITUDINAL",

  "PRJNA851469",
  "ICU_BACKGROUND_LONGITUDINAL",

  "PRJNA578267",
  "NONSEPSIS_LONGITUDINAL_CONTROL",

  "PRJNA430161",
  "INTERVENTION_LONGITUDINAL_SUPPORT",

  "PRJNA1166732",
  "INTERVENTION_LONGITUDINAL_SUPPORT",

  "PRJNA978257",
  "STATIC_SUPPORT"
)

# ------------------------------------------------------------
# Main
# ------------------------------------------------------------

main <- function() {

  ck("STEP87B STARTED")

  # ==========================================================
  # 1. Guard Step87A2
  # ==========================================================

  complete_flag <- file.path(
    STEP87A2,
    "_STEP87A_COMPLETE.ok"
  )

  full_meta_path <- file.path(
    STEP87A2,
    "V2_STEP87A2_FULL_analysis_metadata_786samples.csv"
  )

  dup_group_path <- file.path(
    STEP87A2,
    "V2_STEP87A2_duplicate_patient_time_groups.csv"
  )

  dup_detail_path <- file.path(
    STEP87A2,
    "V2_STEP87A2_duplicate_patient_time_detail.csv"
  )

  inconsistent_path <- file.path(
    STEP87A2,
    "V2_STEP87A2_step86_vs_step81D_key_field_inconsistencies.csv"
  )

  unmatched_path <- file.path(
    STEP87A2,
    "V2_STEP87A2_unmatched_final_runs_to_frozen_master.csv"
  )

  required87A2 <- c(
    complete_flag,
    full_meta_path,
    dup_group_path,
    dup_detail_path,
    inconsistent_path,
    unmatched_path
  )

  if (!all(file.exists(required87A2))) {
    stop(
      paste0(
        "Required Step87A2 files are missing:\n",
        paste(
          required87A2[!file.exists(required87A2)],
          collapse = "\n"
        )
      )
    )
  }

  meta786 <- safe_csv(full_meta_path)
  dup_groups <- safe_csv(dup_group_path)
  dup_detail <- safe_csv(dup_detail_path)
  inconsist <- safe_csv(inconsistent_path)
  unmatched <- safe_csv(unmatched_path)

  if (
    nrow(meta786) != 786 ||
    n_distinct(meta786$project) != 8
  ) {
    stop(
      "Step87A2 full metadata guard failed. Expected 786 samples / 8 projects; observed ",
      nrow(meta786),
      " samples / ",
      n_distinct(meta786$project),
      " projects."
    )
  }

  if (
    nrow(inconsist) != 0 ||
    nrow(unmatched) != 0
  ) {
    stop(
      "Step87A2 still contains mapping inconsistencies or unmatched Runs. Do NOT proceed to Step87B."
    )
  }

  if (
    nrow(dup_groups) != 1 ||
    dup_groups$project[1] != "PRJNA516701" ||
    as.character(dup_groups$patient_id[1]) != "4200" ||
    dup_groups$timepoint_key[1] != "DAY_5" ||
    dup_groups$n_runs_same_patient_time[1] != 2
  ) {
    stop(
      "Duplicate patient-time state changed. Expected only PRJNA516701 / patient 4200 / DAY_5 / 2 Runs."
    )
  }

  expected_pair <- c(
    "SRR8488668",
    "SRR8488665"
  )

  if (!setequal(dup_detail$run_id, expected_pair)) {
    stop(
      "Expected duplicate Run pair not found. Observed: ",
      paste(dup_detail$run_id, collapse = ", ")
    )
  }

  dd <- dup_detail |>
    filter(
      run_id %in% expected_pair
    ) |>
    arrange(
      desc(reads_after_taxonomy_preview)
    )

  if (
    nrow(dd) != 2 ||
    any(is.na(dd$reads_after_taxonomy_preview)) ||
    dd$reads_after_taxonomy_preview[1] <=
      dd$reads_after_taxonomy_preview[2]
  ) {
    stop(
      "Unable to determine a unique higher-depth representative for the duplicate pair."
    )
  }

  PRIMARY_REP_KEEP <- dd$run_id[1]
  PRIMARY_REP_HOLD <- dd$run_id[2]

  if (
    PRIMARY_REP_KEEP != "SRR8488668" ||
    PRIMARY_REP_HOLD != "SRR8488665" ||
    dd$reads_after_taxonomy_preview[1] != 45312 ||
    dd$reads_after_taxonomy_preview[2] != 10823 ||
    !all(as.numeric(dd$time_day) == 5) ||
    !all(as.character(dd$collection_date) == "2017-07-17")
  ) {
    stop(
      paste0(
        "Duplicate evidence differs from the frozen Step87A2 result.\n",
        "Do not bypass this guard. Review V2_STEP87A2_duplicate_patient_time_detail.csv first."
      )
    )
  }

  ck("STEP87A2 STATE GUARDED")

  # ==========================================================
  # 2. Freeze primary replicate rule
  # ==========================================================

  replicate_decision <- dup_detail |>
    mutate(
      analysis_replicate_rule =
        "ONE_OBSERVATION_PER_PATIENT_TIMEPOINT",

      source_resolution =
        "BIOLOGICAL_EQUIVALENCE_NOT_ESTABLISHED_FROM_AVAILABLE_METADATA",

      primary_decision = case_when(
        run_id == PRIMARY_REP_KEEP ~
          "KEEP_PRIMARY_HIGHER_TARGET_READ_DEPTH",

        run_id == PRIMARY_REP_HOLD ~
          "HOLD_FROM_PRIMARY_KEEP_FOR_SENSITIVITY",

        TRUE ~
          "UNEXPECTED"
      ),

      canonical_data_deleted = FALSE,

      rationale = case_when(
        run_id == PRIMARY_REP_KEEP ~
          paste0(
            "Same patient/day/date pair; retained as deterministic ",
            "higher-depth representative for one-observation-per-timepoint models."
          ),

        run_id == PRIMARY_REP_HOLD ~
          paste0(
            "Same patient/day/date pair; not declared biologically invalid. ",
            "Held from primary analysis and retained for sensitivity because ",
            "finer replicate provenance is unavailable."
          ),

        TRUE ~ NA_character_
      )
    )

  write_excel_csv(
    replicate_decision,
    file.path(
      OUT,
      "V2_STEP87B_replication_resolution_decision.csv"
    ),
    na = ""
  )

  meta_flagged <- meta786 |>
    mutate(
      primary_patient_time_include =
        run_id != PRIMARY_REP_HOLD,

      replicate_primary_status = case_when(
        run_id == PRIMARY_REP_KEEP ~
          "PRIMARY_REPRESENTATIVE",

        run_id == PRIMARY_REP_HOLD ~
          "SENSITIVITY_ONLY_SAME_DAY_REPLICATE",

        TRUE ~
          "STANDARD_SINGLE_OBSERVATION"
      )
    )

  meta785 <- meta_flagged |>
    filter(
      primary_patient_time_include
    )

  if (nrow(meta785) != 785) {
    stop(
      "Primary metadata row count guard failed. Expected 785; observed ",
      nrow(meta785),
      "."
    )
  }

  duplicate_after <- meta785 |>
    count(
      project,
      patient_id,
      timepoint_key,
      name = "n"
    ) |>
    filter(n > 1)

  if (nrow(duplicate_after) != 0) {
    write_excel_csv(
      duplicate_after,
      file.path(
        OUT,
        "V2_STEP87B_UNEXPECTED_duplicate_patient_time_after_resolution.csv"
      ),
      na = ""
    )

    stop(
      "Patient-time duplicates remain after Step87B replicate resolution."
    )
  }

  st <- meta785 |>
    group_by(project) |>
    group_modify(
      ~patient_time_stats(.x)
    ) |>
    ungroup()

  if (
    sum(st$patients) != 390 ||
    sum(st$patients_GE2) != 258 ||
    sum(st$patients_GE3) != 132
  ) {
    stop(
      paste0(
        "Replicate handling unexpectedly changed longitudinal eligibility.\n",
        "Observed patients / GE2 / GE3 = ",
        sum(st$patients), " / ",
        sum(st$patients_GE2), " / ",
        sum(st$patients_GE3)
      )
    )
  }

  write_excel_csv(
    meta_flagged,
    file.path(
      OUT,
      "V2_STEP87B_metadata_786_WITH_primary_replicate_flag.csv"
    ),
    na = ""
  )

  write_excel_csv(
    meta785,
    file.path(
      OUT,
      "V2_STEP87B_PRIMARY_analysis_metadata_785_patient_time_observations.csv"
    ),
    na = ""
  )

  write_excel_csv(
    st,
    file.path(
      OUT,
      "V2_STEP87B_primary_longitudinal_eligibility_by_project.csv"
    ),
    na = ""
  )

  ck("PRIMARY PATIENT-TIME REPLICATE RULE FROZEN")

  # ==========================================================
  # 3. Load Step86B final sequence registry
  # ==========================================================

  seq_registry_path <- file.path(
    STEP86B,
    "V2_STEP86B_FINAL_analysis_ready_registry.csv"
  )

  if (!file.exists(seq_registry_path)) {
    stop(
      "Step86B sequence registry missing: ",
      seq_registry_path
    )
  }

  seq_registry <- safe_csv(seq_registry_path)

  required_registry_cols <- c(
    "project",
    "samples_final",
    "final_seqtab_path"
  )

  if (!all(required_registry_cols %in% names(seq_registry))) {
    stop(
      "Step86B registry is missing required columns: ",
      paste(
        setdiff(
          required_registry_cols,
          names(seq_registry)
        ),
        collapse = ", "
      )
    )
  }

  if (
    nrow(seq_registry) != 8 ||
    sum(seq_registry$samples_final) != 786
  ) {
    stop(
      "Step86B sequence registry guard failed. Expected 8 projects / 786 samples."
    )
  }

  # ==========================================================
  # 4. Build one portable analysis object per cohort
  # ==========================================================

  object_registry_rows <- list()

  for (proj in role_registry$project) {

    ck(
      paste0(
        proj,
        " OBJECT BUILD START"
      )
    )

    rr <- seq_registry |>
      filter(
        .data$project == .env$proj
      )

    if (nrow(rr) != 1) {
      stop(
        proj,
        ": Step86B registry row missing or non-unique."
      )
    }

    seq_path <- as.character(
      rr$final_seqtab_path[1]
    )

    if (!file.exists(seq_path)) {
      stop(
        proj,
        ": final seqtab missing: ",
        seq_path
      )
    }

    counts_all <- readRDS(seq_path)

    if (
      !is.matrix(counts_all) ||
      nrow(counts_all) != rr$samples_final[1] ||
      ncol(counts_all) == 0
    ) {
      stop(
        proj,
        ": final seqtab dimension/type guard failed."
      )
    }

    if (
      is.null(rownames(counts_all)) ||
      anyDuplicated(rownames(counts_all))
    ) {
      stop(
        proj,
        ": final seqtab rownames are missing or duplicated."
      )
    }

    md_all <- meta_flagged |>
      filter(
        .data$project == .env$proj
      )

    md_primary <- meta785 |>
      filter(
        .data$project == .env$proj
      )

    missing_runs <- setdiff(
      md_all$run_id,
      rownames(counts_all)
    )

    if (length(missing_runs)) {
      writeLines(
        missing_runs,
        file.path(
          OUT,
          paste0(
            proj,
            "_MISSING_RUNS_IN_SEQTAB.txt"
          )
        )
      )

      stop(
        proj,
        ": metadata Runs are not all present in final seqtab. See *_MISSING_RUNS_IN_SEQTAB.txt"
      )
    }

    counts_primary <- counts_all[
      md_primary$run_id,
      ,
      drop = FALSE
    ]

    # Drop ASVs that become zero after excluding the one held replicate.
    counts_primary <- counts_primary[
      ,
      colSums(counts_primary) > 0,
      drop = FALSE
    ]

    if (
      nrow(counts_primary) != nrow(md_primary) ||
      !identical(
        rownames(counts_primary),
        md_primary$run_id
      )
    ) {
      stop(
        proj,
        ": count/metadata row alignment failed."
      )
    }

    tax_path <- taxonomy_path(proj)
    tax <- safe_csv(tax_path)

    seq_col <- detect_asv_sequence_col(tax)

    if (seq_col != "ASV_sequence") {
      names(tax)[names(tax) == seq_col] <- "ASV_sequence"
    }

    if (anyDuplicated(tax$ASV_sequence)) {
      stop(
        proj,
        ": duplicated ASV_sequence values in taxonomy table."
      )
    }

    missing_tax <- setdiff(
      colnames(counts_primary),
      tax$ASV_sequence
    )

    if (length(missing_tax)) {
      writeLines(
        missing_tax,
        file.path(
          OUT,
          paste0(
            proj,
            "_MISSING_ASVS_IN_TAXONOMY.txt"
          )
        )
      )

      stop(
        proj,
        ": taxonomy does not cover all final ASVs. See *_MISSING_ASVS_IN_TAXONOMY.txt"
      )
    }

    tax_primary <- tax[
      match(
        colnames(counts_primary),
        tax$ASV_sequence
      ),
      ,
      drop = FALSE
    ]

    if (
      nrow(tax_primary) != ncol(counts_primary) ||
      !identical(
        as.character(tax_primary$ASV_sequence),
        colnames(counts_primary)
      )
    ) {
      stop(
        proj,
        ": taxonomy-to-ASV exact ordering guard failed."
      )
    }

    pt <- md_primary |>
      distinct(
        patient_id,
        timepoint_key
      ) |>
      count(
        patient_id,
        name =
          "n_distinct_timepoints_primary"
      )

    md_primary <- md_primary |>
      left_join(
        pt,
        by = "patient_id"
      ) |>
      mutate(
        longitudinal_GE2_primary =
          n_distinct_timepoints_primary >= 2,

        longitudinal_GE3_primary =
          n_distinct_timepoints_primary >= 3,

        trajectory_module = case_when(
          analysis_role ==
            "CORE_SEPSIS_LONGITUDINAL" ~
            "CORE_SEPSIS_TRAJECTORY",

          analysis_role ==
            "ICU_INFECTION_LONGITUDINAL_EXTERNAL" ~
            "ICU_INFECTION_EXTERNAL_TRAJECTORY",

          analysis_role ==
            "ICU_BACKGROUND_LONGITUDINAL" ~
            "ICU_BACKGROUND_TRAJECTORY",

          analysis_role ==
            "NONSEPSIS_LONGITUDINAL_CONTROL" ~
            "NONSEPSIS_CONTROL_TRAJECTORY",

          analysis_role ==
            "INTERVENTION_LONGITUDINAL_SUPPORT" ~
            "INTERVENTION_SUPPORT_TRAJECTORY",

          analysis_role ==
            "STATIC_SUPPORT" ~
            "STATIC_SUPPORT_ONLY",

          TRUE ~
            "REVIEW"
        )
      )

    same_day_sensitivity <- NULL

    if (proj == "PRJNA516701") {

      pair_runs <- c(
        PRIMARY_REP_KEEP,
        PRIMARY_REP_HOLD
      )

      pair_metadata <- md_all |>
        filter(
          run_id %in% pair_runs
        )

      pair_counts <- counts_all[
        pair_runs,
        ,
        drop = FALSE
      ]

      same_day_sensitivity <- list(
        pair_metadata = pair_metadata,
        pair_counts = pair_counts,
        primary_keep = PRIMARY_REP_KEEP,
        sensitivity_hold = PRIMARY_REP_HOLD
      )
    }

    provenance <- list(
      project = proj,
      created = as.character(Sys.time()),
      step = "87B",
      sequence_source = seq_path,
      taxonomy_source = tax_path,
      metadata_source = full_meta_path,
      minimum_target_reads = 2000L,
      taxonomy_reference =
        "SILVA 138.2 SSU NR99",
      replicate_rule = paste0(
        "One observation per project/patient/timepoint. ",
        "PRJNA516701 patient 4200 DAY_5: ",
        PRIMARY_REP_KEEP,
        " retained as higher-depth primary representative; ",
        PRIMARY_REP_HOLD,
        " retained for sensitivity."
      ),
      cross_cohort_ASV_merge_allowed = FALSE
    )

    obj <- list(
      counts = counts_primary,
      taxonomy = tax_primary,
      metadata = md_primary,
      provenance = provenance,
      same_day_replicate_sensitivity =
        same_day_sensitivity
    )

    object_dir <- file.path(
      PROJECT_ROOT,
      proj,
      "05_work",
      "step87B_analysis_object"
    )

    dir.create(
      object_dir,
      recursive = TRUE,
      showWarnings = FALSE
    )

    object_path <- file.path(
      object_dir,
      paste0(
        proj,
        "_analysis_object_step87B.rds"
      )
    )

    atomic_saveRDS(
      obj,
      object_path
    )

    write_excel_csv(
      md_primary,
      file.path(
        object_dir,
        paste0(
          proj,
          "_analysis_metadata_step87B.csv"
        )
      ),
      na = ""
    )

    role_this <- role_registry |>
      filter(
        .data$project == .env$proj
      ) |>
      pull(analysis_role)

    object_registry_rows[[
      length(object_registry_rows) + 1
    ]] <- tibble(
      project = proj,
      analysis_role = role_this,
      samples_step86B = nrow(counts_all),
      primary_samples = nrow(counts_primary),
      primary_patients =
        n_distinct(md_primary$patient_id),
      primary_patients_GE2 =
        sum(
          pt$n_distinct_timepoints_primary >= 2
        ),
      primary_patients_GE3 =
        sum(
          pt$n_distinct_timepoints_primary >= 3
        ),
      final_ASVs = ncol(counts_primary),
      final_reads = sum(counts_primary),
      analysis_object_path = object_path
    )

    ck(
      paste0(
        proj,
        " OBJECT BUILD COMPLETE"
      )
    )

    rm(
      counts_all,
      counts_primary,
      tax,
      tax_primary,
      md_all,
      md_primary,
      obj,
      same_day_sensitivity
    )

    gc(verbose = FALSE)
  }

  object_registry <- bind_rows(
    object_registry_rows
  )

  if (
    nrow(object_registry) != 8 ||
    sum(object_registry$primary_samples) != 785 ||
    sum(object_registry$primary_patients) != 390 ||
    sum(object_registry$primary_patients_GE2) != 258 ||
    sum(object_registry$primary_patients_GE3) != 132
  ) {
    write_excel_csv(
      object_registry,
      file.path(
        OUT,
        "V2_STEP87B_FAILED_object_registry_snapshot.csv"
      ),
      na = ""
    )

    stop(
      paste0(
        "Global object-count guard failed.\n",
        "Expected projects/samples/patients/GE2/GE3 = 8/785/390/258/132.\n",
        "Observed = ",
        nrow(object_registry), "/",
        sum(object_registry$primary_samples), "/",
        sum(object_registry$primary_patients), "/",
        sum(object_registry$primary_patients_GE2), "/",
        sum(object_registry$primary_patients_GE3)
      )
    )
  }

  write_excel_csv(
    object_registry,
    file.path(
      OUT,
      "V2_STEP87B_analysis_object_registry.csv"
    ),
    na = ""
  )

  # ==========================================================
  # 5. Freeze analysis modules
  # ==========================================================

  modules <- tribble(
    ~analysis_module,
    ~projects,
    ~purpose,

    "CORE_SEPSIS_TRAJECTORY",
    "PRJNA691455",
    "Primary current sepsis longitudinal trajectory.",

    "ICU_INFECTION_EXTERNAL_TRAJECTORY",
    "PRJEB82425",
    "External ICU infection/event-stage trajectory.",

    "ICU_BACKGROUND_TRAJECTORY",
    "PRJNA516701; PRJNA851469",
    "General critical-illness longitudinal background trajectories.",

    "NONSEPSIS_CONTROL_TRAJECTORY",
    "PRJNA578267",
    "Non-sepsis longitudinal comparator.",

    "INTERVENTION_SUPPORT_TRAJECTORY",
    "PRJNA430161; PRJNA1166732",
    "Supportive intervention trajectories; not pooled as natural-history sepsis cohorts.",

    "STATIC_SUPPORT",
    "PRJNA978257",
    "Static supportive comparison only."
  )

  write_excel_csv(
    modules,
    file.path(
      OUT,
      "V2_STEP87B_analysis_module_registry.csv"
    ),
    na = ""
  )

  # ==========================================================
  # 6. Documentation + completion flag
  # ==========================================================

  readme <- c(
    "SEPSIS V2 - STEP87B ANALYSIS OBJECT FREEZE",
    paste0("Created: ", Sys.time()),
    "",
    "REPLICATE RESOLUTION",
    "PRJNA516701 patient 4200 DAY_5 has two different BioSample/Run records on the same date.",
    "Available frozen metadata do not establish biological equivalence.",
    "",
    "PRIMARY RULE",
    paste0(
      PRIMARY_REP_KEEP,
      " (45,312 target reads) retained as deterministic higher-depth representative."
    ),
    paste0(
      PRIMARY_REP_HOLD,
      " (10,823 target reads) held from primary day-level models but preserved for sensitivity."
    ),
    "Canonical source data are not deleted and the held Run is not labelled biologically invalid.",
    "",
    "PRIMARY STATE",
    "785 unique project/patient/timepoint observations.",
    "390 project-specific patients.",
    "258 patients with >=2 timepoints.",
    "132 patients with >=3 timepoints.",
    "",
    "ANALYSIS OBJECTS",
    "One RDS object per cohort containing counts, SILVA taxonomy, full metadata, provenance, and an optional same-day replicate sensitivity component.",
    "ASV identities remain cohort-specific. No cross-region ASV merge.",
    "",
    "NEXT",
    "Step88: cohort-specific alpha diversity + within-patient beta-diversity change.",
    "Do not pool all eight cohorts into one raw-ASV model."
  )

  writeLines(
    readme,
    file.path(
      OUT,
      "README_STEP87B.txt"
    ),
    useBytes = TRUE
  )

  writeLines(
    c(
      paste0(
        "Completed: ",
        Sys.time()
      ),
      "Status: STEP87B COMPLETE",
      "Primary patient-time observations: 785",
      "Project-specific patients: 390",
      "Patients GE2: 258",
      "Patients GE3: 132",
      paste0(
        "Primary replicate kept: ",
        PRIMARY_REP_KEEP
      ),
      paste0(
        "Sensitivity replicate held: ",
        PRIMARY_REP_HOLD
      ),
      "Eight cohort-specific analysis objects built."
    ),
    file.path(
      OUT,
      "_STEP87B_COMPLETE.ok"
    ),
    useBytes = TRUE
  )

  ck("STEP87B COMPLETE")

  cat("\n============================================================\n")
  cat("SEPSIS V2 - STEP87B COMPLETE\n")
  cat("============================================================\n\n")

  cat("Replicate decision:\n")
  print(
    replicate_decision |>
      select(
        project,
        patient_id,
        timepoint_key,
        run_id,
        reads_after_taxonomy_preview,
        primary_decision
      ),
    n = Inf,
    width = Inf
  )

  cat("\nAnalysis object registry:\n")
  print(
    object_registry,
    n = Inf,
    width = Inf
  )

  cat("\nOutput directory:\n")
  cat(OUT, "\n")
}

tryCatch(
  main(),
  error = function(e) {

    msg <- c(
      paste0(
        "STEP87B FATAL ERROR: ",
        Sys.time()
      ),
      paste0(
        "Message: ",
        conditionMessage(e)
      ),
      paste0(
        "Call: ",
        paste(
          deparse(
            conditionCall(e)
          ),
          collapse = " "
        )
      )
    )

    writeLines(
      msg,
      ERR,
      useBytes = TRUE
    )

    ck("STEP87B FAILED")

    message(
      paste(
        msg,
        collapse = "\n"
      )
    )

    quit(
      save = "no",
      status = 1,
      runLast = FALSE
    )
  }
)
