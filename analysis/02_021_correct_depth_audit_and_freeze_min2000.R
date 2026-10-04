# ============================================================
# Sepsis V2 - Step 86B
# Correct Step86A longitudinal depth audit bug
# + freeze 2,000-read minimum for the CURRENT 8 taxonomy-ready cohorts
#
# IMPORTANT BUG FIX
#   Step86A used:
#       for (project in registry$project) {
#         d <- depth_meta |> filter(.data$project == project)
#       }
#   Inside dplyr's data mask, RHS `project` resolved to the
#   project COLUMN, so the condition became project == project
#   and every project incorrectly contained all 796 samples.
#
#   This script uses a distinct loop variable `proj` and
#   `.data$project == .env$proj`.
#
# WHY 2,000 READS
#   Corrected audit from the 796 selected samples:
#
#   threshold   samples kept   patients kept   GE2 kept   GE3 kept
#      100         795/796        392/393       263/263    133/133
#      500         792/796        391/393       262/263    132/133
#     1000         790/796        391/393       260/263    132/133
#     2000         786/796        390/393       258/263    132/133
#     5000         784/796        389/393       257/263    132/133
#    10000         780/796        386/393       256/263    132/133
#
#   2,000 reads is selected as the practical elbow:
#     - retains 98.74% of samples;
#     - removes extremely low-depth libraries;
#     - also removes PRJNA516701 SRR8488562, which retained only
#       1,461 target reads after losing 83.3% to taxonomy filtering;
#     - going from 2,000 to 5,000 removes additional samples and
#       longitudinal information with little QC gain.
#
# SCOPE
#   This cutoff is frozen for the CURRENT 8 taxonomy-ready cohorts.
#   Special-route cohorts (PacBio/full-length/FASTA/Ion/etc.) must
#   undergo their own platform-aware audit before harmonization.
#
# R 4.4.0 / Windows
# ============================================================

options(stringsAsFactors = FALSE)

pkgs <- c("readr","dplyr","tibble","purrr")
missing <- pkgs[
  !vapply(pkgs, requireNamespace, logical(1), quietly=TRUE)
]

if (length(missing)) {
  install.packages(
    missing,
    repos="https://cloud.r-project.org"
  )
}

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tibble)
  library(purrr)
})

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

STEP85C86A <- file.path(
  ROOT,
  "results",
  "V2_25C_86A_TAXONOMY_FILTER_AND_LONGITUDINAL_DEPTH_AUDIT"
)

OUT <- file.path(
  ROOT,
  "results",
  "V2_26B_DEPTH2000_FINAL_ANALYSIS_READY"
)

dir.create(
  OUT,
  recursive=TRUE,
  showWarnings=FALSE
)

LOG <- file.path(
  OUT,
  "_STEP86B_runtime_checkpoints.txt"
)

ERR <- file.path(
  OUT,
  "_STEP86B_FATAL_ERROR.txt"
)

if (file.exists(ERR)) unlink(ERR)

ck <- function(x) {
  cat(
    paste0(x, ": ", Sys.time(), "\n"),
    file=LOG,
    append=TRUE
  )
}

safe_csv <- function(p) {
  suppressMessages(
    read_csv(
      p,
      show_col_types=FALSE,
      progress=FALSE,
      name_repair="unique"
    )
  )
}

atomic_saveRDS <- function(obj,p) {
  tmp <- paste0(p,".tmp")
  if (file.exists(tmp)) unlink(tmp)
  saveRDS(obj,tmp)
  if (file.exists(p)) unlink(p)
  if (!file.rename(tmp,p)) {
    stop("Could not atomically save: ",p)
  }
}

registry <- tribble(
  ~project, ~expected_samples_before,

  "PRJNA430161", 26L,
  "PRJNA691455", 49L,
  "PRJNA978257", 13L,
  "PRJEB82425", 92L,
  "PRJNA516701", 170L,
  "PRJNA578267", 207L,
  "PRJNA851469", 119L,
  "PRJNA1166732", 120L
)

thresholds <- c(
  100L,
  500L,
  1000L,
  2000L,
  5000L,
  10000L
)

FROZEN_DEPTH <- 2000L

patient_time_stats <- function(d) {

  pt <- d |>
    distinct(
      patient_id,
      timepoint_key
    ) |>
    count(
      patient_id,
      name="n_timepoints"
    )

  tibble(
    patients=n_distinct(
      d$patient_id
    ),
    GE2=sum(
      pt$n_timepoints>=2
    ),
    GE3=sum(
      pt$n_timepoints>=3
    )
  )
}

main <- function() {

  ck("STEP86B STARTED")

  map_path <- file.path(
    STEP85C86A,
    "V2_STEP86A_selected_run_patient_time_map.csv"
  )

  freeze_registry_path <- file.path(
    STEP85C86A,
    "V2_STEP85C_taxonomy_filter_freeze_registry.csv"
  )

  if (
    !file.exists(map_path) ||
    !file.exists(freeze_registry_path)
  ) {
    stop(
      "Required Step85C/86A input files are missing."
    )
  }

  depth_meta <- safe_csv(
    map_path
  )

  freeze_registry <- safe_csv(
    freeze_registry_path
  )

  if (
    nrow(depth_meta)!=796 ||
    n_distinct(depth_meta$project)!=8 ||
    any(
      is.na(depth_meta$patient_id)
    )
  ) {
    stop(
      "Selected run/patient/time map hard guard failed."
    )
  }

  if (
    nrow(freeze_registry)!=8 ||
    sum(freeze_registry$samples)!=796
  ) {
    stop(
      "Taxonomy-filter freeze registry hard guard failed."
    )
  }

  # ==========================================================
  # 1. Recompute the longitudinal threshold audit correctly
  # ==========================================================

  impact_rows <- list()
  excluded_rows <- list()

  for (proj in registry$project) {

    d <- depth_meta |>
      filter(
        .data$project == .env$proj
      )

    expected <- registry |>
      filter(
        .data$project == .env$proj
      ) |>
      pull(
        expected_samples_before
      )

    if (
      nrow(d)!=expected
    ) {
      stop(
        proj,
        ": project sample-count guard failed. Expected ",
        expected,
        "; observed ",
        nrow(d)
      )
    }

    before_stats <- patient_time_stats(
      d
    )

    for (thr in thresholds) {

      kept <- d |>
        filter(
          reads_after_taxonomy_preview>=thr
        )

      excluded <- d |>
        filter(
          reads_after_taxonomy_preview<thr
        )

      after_stats <- patient_time_stats(
        kept
      )

      impact_rows[[
        length(impact_rows)+1
      ]] <- tibble(
        project=proj,
        threshold=thr,

        samples_before=nrow(d),
        samples_after=nrow(kept),
        samples_removed=nrow(excluded),

        patients_before=before_stats$patients,
        patients_after=after_stats$patients,
        patients_lost=
          before_stats$patients-
          after_stats$patients,

        GE2_before=before_stats$GE2,
        GE2_after=after_stats$GE2,
        GE2_lost=
          before_stats$GE2-
          after_stats$GE2,

        GE3_before=before_stats$GE3,
        GE3_after=after_stats$GE3,
        GE3_lost=
          before_stats$GE3-
          after_stats$GE3
      )

      if (nrow(excluded)) {

        excluded_rows[[
          length(excluded_rows)+1
        ]] <- excluded |>
          mutate(
            threshold=thr
          )
      }
    }
  }

  impact <- bind_rows(
    impact_rows
  )

  global_impact <- impact |>
    group_by(
      threshold
    ) |>
    summarise(
      samples_before=sum(samples_before),
      samples_after=sum(samples_after),
      samples_removed=sum(samples_removed),

      patients_before=sum(patients_before),
      patients_after=sum(patients_after),
      patients_lost=sum(patients_lost),

      GE2_before=sum(GE2_before),
      GE2_after=sum(GE2_after),
      GE2_lost=sum(GE2_lost),

      GE3_before=sum(GE3_before),
      GE3_after=sum(GE3_after),
      GE3_lost=sum(GE3_lost),

      .groups="drop"
    )

  expected_global <- tribble(
    ~threshold, ~samples_after, ~patients_after, ~GE2_after, ~GE3_after,
    100L,   795L, 392L, 263L, 133L,
    500L,   792L, 391L, 262L, 132L,
    1000L,  790L, 391L, 260L, 132L,
    2000L,  786L, 390L, 258L, 132L,
    5000L,  784L, 389L, 257L, 132L,
    10000L, 780L, 386L, 256L, 132L
  )

  check <- global_impact |>
    select(
      threshold,
      samples_after,
      patients_after,
      GE2_after,
      GE3_after
    ) |>
    left_join(
      expected_global,
      by="threshold",
      suffix=c("_observed","_expected")
    )

  if (
    any(
      check$samples_after_observed !=
      check$samples_after_expected
    ) ||
    any(
      check$patients_after_observed !=
      check$patients_after_expected
    ) ||
    any(
      check$GE2_after_observed !=
      check$GE2_after_expected
    ) ||
    any(
      check$GE3_after_observed !=
      check$GE3_after_expected
    )
  ) {
    stop(
      "Corrected global longitudinal-impact hard guard failed."
    )
  }

  write_excel_csv(
    impact,
    file.path(
      OUT,
      "V2_STEP86B_CORRECTED_depth_threshold_impact_BY_PROJECT.csv"
    ),
    na=""
  )

  write_excel_csv(
    global_impact,
    file.path(
      OUT,
      "V2_STEP86B_CORRECTED_depth_threshold_impact_GLOBAL.csv"
    ),
    na=""
  )

  write_excel_csv(
    bind_rows(
      excluded_rows
    ),
    file.path(
      OUT,
      "V2_STEP86B_candidate_excluded_runs_ALL_THRESHOLDS.csv"
    ),
    na=""
  )

  ck("CORRECTED THRESHOLD AUDIT COMPLETE")

  # ==========================================================
  # 2. Freeze 2,000-read rule
  # ==========================================================

  frozen_exclusions <- depth_meta |>
    filter(
      reads_after_taxonomy_preview <
      FROZEN_DEPTH
    ) |>
    arrange(
      reads_after_taxonomy_preview
    ) |>
    mutate(
      exclusion_stage=
        "POST_TAXONOMY_LIBRARY_DEPTH_QC",

      exclusion_reason=
        paste0(
          "TARGET_PROKARYOTIC_READS_LT_",
          FROZEN_DEPTH
        ),

      final_sample_status=
        "EXCLUDE_FROM_PRIMARY_MICROBIOME_ANALYSIS"
    )

  if (nrow(frozen_exclusions)!=10) {
    stop(
      "Expected exactly 10 samples below 2,000 reads; observed ",
      nrow(frozen_exclusions)
    )
  }

  if (
    !any(
      frozen_exclusions$project=="PRJNA516701" &
      frozen_exclusions$run_id=="SRR8488562" &
      frozen_exclusions$reads_after_taxonomy_preview==1461
    )
  ) {
    stop(
      "Expected SRR8488562 to be excluded by 2,000-read rule."
    )
  }

  write_excel_csv(
    frozen_exclusions,
    file.path(
      OUT,
      "V2_STEP86B_FINAL_excluded_samples_depth2000.csv"
    ),
    na=""
  )

  # ==========================================================
  # 3. Create final analysis-ready tables per cohort
  # ==========================================================

  final_registry_rows <- list()
  final_metadata_rows <- list()

  for (proj in registry$project) {

    fr <- freeze_registry |>
      filter(
        .data$project == .env$proj
      )

    if (nrow(fr)!=1) {
      stop(
        proj,
        ": freeze registry row missing/non-unique."
      )
    }

    seq_path <- fr$frozen_seqtab_path[1]

    if (!file.exists(seq_path)) {
      stop(
        proj,
        ": frozen taxonomy-filtered sequence table missing: ",
        seq_path
      )
    }

    seqtab <- readRDS(
      seq_path
    )

    d <- depth_meta |>
      filter(
        .data$project == .env$proj
      )

    keep_runs <- d |>
      filter(
        reads_after_taxonomy_preview>=
        FROZEN_DEPTH
      ) |>
      pull(
        run_id
      )

    if (
      !all(
        keep_runs %in%
        rownames(seqtab)
      )
    ) {
      stop(
        proj,
        ": kept Run IDs do not all exist in sequence table."
      )
    }

    final_seqtab <- seqtab[
      keep_runs,
      ,
      drop=FALSE
    ]

    # Remove features that become zero after sample filtering.
    feature_keep <- colSums(
      final_seqtab
    )>0

    final_seqtab <- final_seqtab[
      ,
      feature_keep,
      drop=FALSE
    ]

    if (
      any(
        rowSums(final_seqtab) <
        FROZEN_DEPTH
      )
    ) {
      stop(
        proj,
        ": final sequence table contains sample below frozen threshold."
      )
    }

    final_dir <- file.path(
      PROJECT_ROOT,
      proj,
      "05_work",
      "step86B_analysis_ready"
    )

    dir.create(
      final_dir,
      recursive=TRUE,
      showWarnings=FALSE
    )

    final_path <- file.path(
      final_dir,
      paste0(
        proj,
        "_analysis_ready_taxfiltered_min2000.rds"
      )
    )

    atomic_saveRDS(
      final_seqtab,
      final_path
    )

    final_meta <- d |>
      filter(
        run_id %in%
        rownames(final_seqtab)
      ) |>
      mutate(
        final_sample_status=
          "INCLUDE_PRIMARY_ANALYSIS",

        frozen_min_reads=
          FROZEN_DEPTH
      )

    write_excel_csv(
      final_meta,
      file.path(
        final_dir,
        paste0(
          proj,
          "_analysis_ready_metadata_min2000.csv"
        )
      ),
      na=""
    )

    st <- patient_time_stats(
      final_meta
    )

    final_registry_rows[[
      length(final_registry_rows)+1
    ]] <- tibble(
      project=proj,
      samples_before=nrow(seqtab),
      samples_final=nrow(final_seqtab),
      samples_removed=nrow(seqtab)-nrow(final_seqtab),

      patients_final=st$patients,
      patients_GE2_final=st$GE2,
      patients_GE3_final=st$GE3,

      ASVs_final=ncol(final_seqtab),
      reads_final=sum(final_seqtab),

      minimum_reads_final=min(
        rowSums(final_seqtab)
      ),

      median_reads_final=median(
        rowSums(final_seqtab)
      ),

      final_seqtab_path=final_path
    )

    final_metadata_rows[[
      length(final_metadata_rows)+1
    ]] <- final_meta

    ck(
      paste0(
        proj,
        " ANALYSIS-READY TABLE COMPLETE"
      )
    )
  }

  final_registry <- bind_rows(
    final_registry_rows
  )

  final_metadata <- bind_rows(
    final_metadata_rows
  )

  if (
    sum(final_registry$samples_before)!=796 ||
    sum(final_registry$samples_final)!=786 ||
    sum(final_registry$samples_removed)!=10 ||
    nrow(final_metadata)!=786
  ) {
    stop(
      "Final 2,000-read analysis-ready global count guard failed."
    )
  }

  final_patient_stats <- final_metadata |>
    group_by(
      project
    ) |>
    group_modify(
      ~patient_time_stats(.x)
    ) |>
    ungroup()

  if (
    sum(final_patient_stats$patients)!=390 ||
    sum(final_patient_stats$GE2)!=258 ||
    sum(final_patient_stats$GE3)!=132
  ) {
    stop(
      "Final longitudinal eligibility hard guard failed."
    )
  }

  write_excel_csv(
    final_registry,
    file.path(
      OUT,
      "V2_STEP86B_FINAL_analysis_ready_registry.csv"
    ),
    na=""
  )

  write_excel_csv(
    final_metadata,
    file.path(
      OUT,
      "V2_STEP86B_FINAL_analysis_ready_metadata_786samples.csv"
    ),
    na=""
  )

  # ==========================================================
  # 4. Rule freeze documentation
  # ==========================================================

  rules <- tribble(
    ~rule_component, ~frozen_value,

    "Taxonomy_database",
    "SILVA 138.2 SSU NR99",

    "Taxonomy_target",
    "Bacteria + Archaea",

    "Taxonomy_exclude",
    "Chloroplast/Mitochondria + other Kingdoms + Kingdom-unassigned",

    "Minimum_target_prokaryotic_reads",
    "2000",

    "Threshold_scope",
    "Current 8 taxonomy-ready cohorts",

    "Special_route_cohorts",
    "Platform-aware QC required before later harmonization"
  )

  write_excel_csv(
    rules,
    file.path(
      OUT,
      "V2_STEP86B_FINAL_QC_RULES.csv"
    ),
    na=""
  )

  rationale <- c(
    "SEPSIS V2 - STEP86B FINAL DEPTH QC",
    paste0("Created: ",Sys.time()),
    "",
    "BUG CORRECTION",
    "The Step86A per-project audit was invalid because the loop variable `project` was shadowed by the data-frame column inside dplyr::filter(), causing each project to contain all 796 samples.",
    "Step85C taxonomy filtering itself was valid and is retained.",
    "",
    "CORRECTED GLOBAL THRESHOLD IMPACT",
    "100 reads:   795 samples; 392 patients; GE2=263; GE3=133",
    "500 reads:   792 samples; 391 patients; GE2=262; GE3=132",
    "1000 reads:  790 samples; 391 patients; GE2=260; GE3=132",
    "2000 reads:  786 samples; 390 patients; GE2=258; GE3=132",
    "5000 reads:  784 samples; 389 patients; GE2=257; GE3=132",
    "10000 reads: 780 samples; 386 patients; GE2=256; GE3=132",
    "",
    "FROZEN PRIMARY THRESHOLD",
    "Minimum 2,000 target bacterial/archaeal reads after taxonomy filtering.",
    "",
    "RATIONALE",
    "2,000 reads removes clearly under-sequenced libraries while retaining 98.74% of current samples.",
    "It removes 10/796 samples, 3/393 project-specific patients, 5/263 >=2-timepoint patients, and 1/133 >=3-timepoint patients.",
    "It also removes PRJNA516701 SRR8488562, the unique taxonomy-loss outlier with 83.3% non-target reads and only 1,461 retained target reads.",
    "Increasing to 5,000 removes additional samples/longitudinal information for limited incremental QC benefit.",
    "",
    "SCOPE",
    "This rule is frozen for the current eight taxonomy-ready cohorts.",
    "Special-route cohorts (e.g. PacBio/full-length, FASTA-only, Ion Torrent) require platform-aware read-depth review before integration.",
    "",
    "NEXT",
    "Proceed to analysis-object construction and primary longitudinal alpha/beta diversity, within-patient change, and cohort-wise trajectory analyses."
  )

  writeLines(
    rationale,
    file.path(
      OUT,
      "README_STEP86B_FINAL_DEPTH_QC.txt"
    ),
    useBytes=TRUE
  )

  writeLines(
    c(
      paste0("Completed: ",Sys.time()),
      "Status: STEP86B COMPLETE",
      "Corrected Step86A per-project scoping bug.",
      "Frozen minimum target-prokaryotic library size: 2000 reads.",
      "Final current analysis-ready samples: 786.",
      "Current analysis-ready GE2 patients: 258.",
      "Current analysis-ready GE3 patients: 132."
    ),
    file.path(
      OUT,
      "_STEP86B_COMPLETE.ok"
    )
  )

  ck("STEP86B COMPLETE")

  cat("\n============================================================\n")
  cat("SEPSIS V2 - STEP86B COMPLETE\n")
  cat("============================================================\n\n")

  cat("CORRECTED GLOBAL THRESHOLD IMPACT:\n")
  print(
    global_impact,
    n=Inf,
    width=Inf
  )

  cat("\nFINAL 2,000-READ ANALYSIS-READY REGISTRY:\n")
  print(
    final_registry,
    n=Inf,
    width=Inf
  )

  cat("\nOutput:\n")
  cat(
    OUT,
    "\n"
  )
}

tryCatch(
  main(),
  error=function(e) {

    msg <- c(
      paste0(
        "STEP86B FATAL ERROR: ",
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
          collapse=" "
        )
      )
    )

    writeLines(
      msg,
      ERR,
      useBytes=TRUE
    )

    ck(
      "STEP86B FAILED"
    )

    message(
      paste(
        msg,
        collapse="\n"
      )
    )

    quit(
      save="no",
      status=1,
      runLast=FALSE
    )
  }
)
