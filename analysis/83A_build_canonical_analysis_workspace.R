# ============================================================
# Sepsis V2 - Step 83A
# Build ONE canonical analysis workspace
#
# PURPOSE
#   Stop repeated broad scanning and place every file needed for
#   downstream V2 analysis under one canonical root:
#
#   E:/sepsis_project/data/_V2_ANALYSIS_READY
#
# IMPORTANT
#   - ORIGINAL FILES ARE NOT MOVED OR DELETED.
#   - Large raw files are exposed in the canonical workspace
#     using NTFS HARD LINKS, so disk space is not duplicated.
#   - Small freeze/route CSV/TXT files are copied.
#   - Only files under E:/sepsis_project/data are accepted as
#     sequence/ASV sources.
#   - No F drive, Downloads, Desktop, or whole-disk scanning.
#
# INPUTS
#   Step81D final frozen manifest
#   Step82 raw inventory / run coverage
#   Step82B seqtab validation
#   Step82C frozen processing routes
#
# R 4.4.0 / Windows
# ============================================================

options(stringsAsFactors = FALSE)

pkgs <- c("readr","dplyr","tidyr","stringr","purrr","tibble")
missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly=TRUE)]
if (length(missing)) {
  install.packages(missing, repos="https://cloud.r-project.org")
}

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(purrr)
  library(tibble)
})

ROOT <- "E:/sepsis_project"
DATA_ROOT <- file.path(ROOT,"data")
RESULTS_ROOT <- file.path(ROOT,"results")

FREEZE_ROOT <- file.path(
  RESULTS_ROOT,
  "V2_21D_FINAL_FREEZE_CLINICAL_REPAIRED"
)

STEP82_ROOT <- file.path(
  RESULTS_ROOT,
  "V2_22_ASV_RAW_READINESS_AUDIT"
)

STEP82B_ROOT <- file.path(
  RESULTS_ROOT,
  "V2_22B_ASV_VALIDATION"
)

STEP82C_ROOT <- file.path(
  RESULTS_ROOT,
  "V2_22C_PROCESSING_ROUTE_FREEZE"
)

CANON <- file.path(
  DATA_ROOT,
  "_V2_ANALYSIS_READY"
)

OUT <- file.path(
  RESULTS_ROOT,
  "V2_23A_CANONICAL_WORKSPACE_BUILD"
)

dir.create(CANON, recursive=TRUE, showWarnings=FALSE)
dir.create(OUT, recursive=TRUE, showWarnings=FALSE)

LOG <- file.path(OUT,"_STEP83A_runtime_checkpoints.txt")
ERR <- file.path(OUT,"_STEP83A_FATAL_ERROR.txt")
if (file.exists(ERR)) unlink(ERR)

ck <- function(x) {
  cat(paste0(x,": ",Sys.time(),"\n"),file=LOG,append=TRUE)
}

clean_chr <- function(x) {
  x <- trimws(as.character(x))
  x[
    is.na(x) |
    x=="" |
    tolower(x) %in% c("na","nan","n/a","null","none")
  ] <- NA_character_
  x
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

norm_path <- function(p) {
  ifelse(
    is.na(p),
    NA_character_,
    normalizePath(
      p,
      winslash="/",
      mustWork=FALSE
    )
  )
}

inside_root <- function(path, root) {
  p <- tolower(norm_path(path))
  r <- tolower(norm_path(root))
  !is.na(p) & startsWith(p,paste0(r,"/"))
}

file_is_regular <- function(p) {
  file.exists(p) && !dir.exists(p)
}

safe_file_size <- function(p) {
  if (!file_is_regular(p)) return(NA_real_)
  as.numeric(file.info(p)$size)
}

# ------------------------------------------------------------
# Hard-link helper
# ------------------------------------------------------------
link_one <- function(src,dst,category,project,run_id=NA_character_) {

  rec <- tibble(
    project=project,
    run_id=run_id,
    category=category,
    source_path=src,
    target_path=dst,
    source_bytes=NA_real_,
    action=NA_character_,
    status=NA_character_,
    note=NA_character_
  )

  if (is.na(src) || !file_is_regular(src)) {
    rec$status <- "SOURCE_MISSING"
    rec$note <- "Source file does not exist."
    return(rec)
  }

  rec$source_bytes <- safe_file_size(src)

  if (!inside_root(src,DATA_ROOT)) {
    rec$status <- "OUTSIDE_CANONICAL_DATA_ROOT"
    rec$note <- "Source is outside E:/sepsis_project/data; not linked."
    return(rec)
  }

  dir.create(dirname(dst),recursive=TRUE,showWarnings=FALSE)

  if (file.exists(dst)) {

    if (safe_file_size(dst)==safe_file_size(src)) {
      rec$action <- "KEEP_EXISTING_TARGET"
      rec$status <- "OK"
      rec$note <- "Target already exists with same byte size."
      return(rec)
    }

    rec$status <- "TARGET_COLLISION_DIFFERENT_SIZE"
    rec$note <- "Target exists but size differs; no overwrite performed."
    return(rec)
  }

  ok <- tryCatch(
    file.link(src,dst),
    error=function(e) FALSE,
    warning=function(w) FALSE
  )

  if (isTRUE(ok)) {
    rec$action <- "NTFS_HARD_LINK"
    rec$status <- "OK"
    rec$note <- "No additional raw-data disk copy created."
  } else {
    rec$status <- "HARD_LINK_FAILED"
    rec$note <- "No copy fallback was used for large data. Check NTFS/same-volume permissions."
  }

  rec
}

copy_small <- function(src,dst) {

  if (!file_is_regular(src)) {
    return(FALSE)
  }

  dir.create(dirname(dst),recursive=TRUE,showWarnings=FALSE)

  file.copy(
    src,
    dst,
    overwrite=TRUE,
    copy.mode=TRUE,
    copy.date=TRUE
  )
}

# ------------------------------------------------------------
# Pick raw files for one expected Run from Step82 inventory
# ------------------------------------------------------------
pick_run_files <- function(candidates, route) {

  if (!nrow(candidates)) return(candidates)

  candidates <- candidates |>
    mutate(
      path_l=tolower(gsub("\\\\","/",file_path)),

      is_processed_path=str_detect(
        path_l,
        "/filtered/|/03_dada2|/02_qc|/04_taxonomy|/trimmed/|/processed/|/denois"
      ),

      # Prefer original/unprocessed paths.
      raw_priority=case_when(
        is_processed_path ~ 9L,
        str_detect(path_l,"/01_raw|/raw/|/download") ~ 1L,
        TRUE ~ 3L
      )
    ) |>
    filter(
      inside_root(file_path,DATA_ROOT)
    )

  if (!nrow(candidates)) return(candidates)

  # Remove obvious processed files if any unprocessed candidate exists.
  if (any(!candidates$is_processed_path)) {
    candidates <- candidates |>
      filter(!is_processed_path)
  }

  # Select largest file per mate/type because raw files should generally
  # be >= derived/filtered versions.
  if (route=="PAIRED_END_FASTQ_DADA2_ROUTE") {

    keep <- candidates |>
      filter(read_mate %in% c("R1","R2")) |>
      group_by(read_mate) |>
      arrange(raw_priority,desc(bytes),file_path) |>
      slice(1) |>
      ungroup()

    return(keep)
  }

  if (route=="FASTA_NO_QUALITY_SPECIAL_ROUTE_NOT_STANDARD_DADA2") {

    keep <- candidates |>
      filter(sequence_file_type %in% c("FASTA_GZ","FASTA")) |>
      arrange(raw_priority,desc(bytes),file_path) |>
      slice(1)

    return(keep)
  }

  # Single-end / unspecified: pick the largest unprocessed sequence file.
  candidates |>
    arrange(raw_priority,desc(bytes),file_path) |>
    slice(1)
}

main <- function() {

  ck("STEP83A STARTED")

  # ==========================================================
  # 1. Inputs
  # ==========================================================

  manifest_path <- file.path(
    FREEZE_ROOT,
    "V2_FINAL_FROZEN_16S_sequence_manifest.csv"
  )

  master_candidates <- list.files(
    FREEZE_ROOT,
    pattern="^V2_FINAL_FROZEN_master_metadata_CLINICAL_REPAIRED_.*\\.csv$",
    full.names=TRUE
  )

  if (!length(master_candidates)) {
    stop("Final Step81D master not found.")
  }

  master_path <- master_candidates[
    which.max(file.info(master_candidates)$mtime)
  ]

  raw_inventory_path <- file.path(
    STEP82_ROOT,
    "V2_STEP82_raw_sequence_file_inventory.csv"
  )

  coverage_path <- file.path(
    STEP82_ROOT,
    "V2_STEP82_expected_run_local_coverage.csv"
  )

  validation_path <- file.path(
    STEP82B_ROOT,
    "V2_STEP82B_seqtab_validation_all.csv"
  )

  route_path <- file.path(
    STEP82C_ROOT,
    "V2_STEP82C_FROZEN_sequence_processing_routes.csv"
  )

  compare_path <- file.path(
    STEP82C_ROOT,
    "V2_STEP82C_competing_full_ASV_comparison.csv"
  )

  companion_path <- file.path(
    STEP82C_ROOT,
    "V2_STEP82C_ASV_companion_file_audit.csv"
  )

  needed <- c(
    manifest_path,
    master_path,
    raw_inventory_path,
    coverage_path,
    validation_path,
    route_path
  )

  if (!all(file.exists(needed))) {
    stop("One or more Step81D/82/82B/82C inputs are missing.")
  }

  manifest <- safe_csv(manifest_path)
  raw_inv <- safe_csv(raw_inventory_path)
  run_cov <- safe_csv(coverage_path)
  validation <- safe_csv(validation_path)
  routes <- safe_csv(route_path)

  if (nrow(manifest)!=1941 || n_distinct(manifest$run_id)!=1941) {
    stop("Frozen 16S manifest guard failed.")
  }

  if (nrow(routes)!=16) {
    stop("Expected 16 frozen 16S project routes.")
  }

  projects <- sort(unique(manifest$project))

  ck("INPUTS GUARDED")

  # ==========================================================
  # 2. Canonical folder structure
  # ==========================================================

  dirs <- c(
    "00_FREEZE",
    "01_PROJECTS",
    "90_PENDING_REPAIR",
    "91_SPECIAL_ROUTES",
    "99_MANIFESTS_LOGS"
  )

  for (d in dirs) {
    dir.create(file.path(CANON,d),recursive=TRUE,showWarnings=FALSE)
  }

  for (project in projects) {

    pbase <- file.path(CANON,"01_PROJECTS",project)

    subdirs <- c(
      "00_metadata",
      "01_raw",
      "02_asv_selected",
      "03_asv_partial_reference",
      "04_taxonomy_qc",
      "05_work",
      "99_notes"
    )

    for (d in subdirs) {
      dir.create(file.path(pbase,d),recursive=TRUE,showWarnings=FALSE)
    }
  }

  ck("CANONICAL DIRECTORIES CREATED")

  # ==========================================================
  # 3. Copy the small frozen control files
  # ==========================================================

  freeze_files <- c(
    master_path,
    manifest_path,
    route_path,
    compare_path
  )

  freeze_files <- freeze_files[
    !is.na(freeze_files) &
    file.exists(freeze_files)
  ]

  for (p in freeze_files) {
    copy_small(
      p,
      file.path(
        CANON,
        "00_FREEZE",
        basename(p)
      )
    )
  }

  # Per-project metadata slices
  for (project in projects) {

    z <- manifest |>
      filter(project==!!project)

    write_excel_csv(
      z,
      file.path(
        CANON,
        "01_PROJECTS",
        project,
        "00_metadata",
        paste0(project,"_frozen_16S_manifest.csv")
      ),
      na=""
    )
  }

  ck("FREEZE FILES COPIED")

  # ==========================================================
  # 4. Build canonical raw links from exact frozen Runs
  # ==========================================================

  link_records <- list()
  missing_records <- list()
  noncanonical <- list()

  for (project in projects) {

    route_row <- routes |>
      filter(project==!!project)

    if (!nrow(route_row)) next

    raw_route <- route_row$raw_processing_route[1]

    expected <- manifest |>
      filter(project==!!project)

    for (i in seq_len(nrow(expected))) {

      run_id <- expected$run_id[i]

      cand <- raw_inv |>
        filter(
          project==!!project,
          run_id_from_filename==!!run_id
        )

      # Record paths outside the canonical E:/sepsis_project/data root,
      # but never use them.
      outside <- cand |>
        filter(
          !inside_root(file_path,DATA_ROOT)
        )

      if (nrow(outside)) {
        noncanonical[[length(noncanonical)+1]] <- outside |>
          transmute(
            project,
            run_id=run_id_from_filename,
            file_path,
            reason="OUTSIDE_E_SEPSIS_PROJECT_DATA"
          )
      }

      chosen <- pick_run_files(
        cand,
        raw_route
      )

      if (!nrow(chosen)) {

        missing_records[[length(missing_records)+1]] <- tibble(
          project=project,
          patient_id=expected$patient_id[i],
          sample_id=expected$sample_id[i],
          run_id=run_id,
          time_raw=expected$time_raw[i],
          analysis_module=expected$analysis_module[i],
          processing_route=route_row$frozen_processing_route[1],
          reason="NO_CANONICAL_LOCAL_RAW_FILE"
        )

        next
      }

      for (j in seq_len(nrow(chosen))) {

        src <- chosen$file_path[j]
        dst <- file.path(
          CANON,
          "01_PROJECTS",
          project,
          "01_raw",
          basename(src)
        )

        link_records[[length(link_records)+1]] <- link_one(
          src=src,
          dst=dst,
          category="RAW_SEQUENCE",
          project=project,
          run_id=run_id
        )
      }
    }
  }

  links_raw <- if (length(link_records)) {
    bind_rows(link_records)
  } else {
    tibble()
  }

  missing_raw <- if (length(missing_records)) {
    bind_rows(missing_records)
  } else {
    tibble()
  }

  noncanonical_df <- if (length(noncanonical)) {
    bind_rows(noncanonical)
  } else {
    tibble()
  }

  write_excel_csv(
    links_raw,
    file.path(OUT,"V2_STEP83A_raw_hardlink_manifest.csv"),
    na=""
  )

  write_excel_csv(
    missing_raw,
    file.path(OUT,"V2_STEP83A_missing_canonical_raw_runs.csv"),
    na=""
  )

  write_excel_csv(
    noncanonical_df,
    file.path(OUT,"V2_STEP83A_sources_outside_canonical_data_root.csv"),
    na=""
  )

  ck("RAW HARD LINKS BUILT")

  # ==========================================================
  # 5. Select/link ASV candidates
  #
  # PRJNA691455:
  # Step82C comparison strongly supports rerun_v2:
  #   - 49/49 frozen samples
  #   - 2,398,535 total reads
  #   - median 50,835 reads/sample
  #   - median richness 148
  # Old root seqtab had median reads=0 / median richness=0.
  # Therefore rerun_v2 becomes the canonical selected ASV.
  # ==========================================================

  asv_records <- list()
  asv_decisions <- list()

  for (project in projects) {

    route <- routes |>
      filter(project==!!project)

    val <- validation |>
      filter(project==!!project)

    if (!nrow(route) || !nrow(val)) next

    selected <- NA_character_
    decision_note <- NA_character_

    if (project=="PRJNA691455") {

      rerun <- val |>
        filter(
          validation_status=="FULL_FROZEN_COVERAGE_CANDIDATE",
          str_detect(
            file_path,
            regex("03_dada2_rerun_v2",ignore_case=TRUE)
          )
        ) |>
        arrange(desc(nrow))

      if (nrow(rerun)) {
        selected <- rerun$file_path[1]
        decision_note <- paste0(
          "Selected rerun_v2 after Step82C comparison: 49/49 samples; ",
          "2,398,535 total reads; median 50,835 reads/sample; ",
          "median richness 148. Legacy full-coverage file had median reads/richness 0."
        )
      }

    } else if (
      route$frozen_processing_route[1]=="REUSE_EXISTING_ASV"
    ) {

      fullv <- val |>
        filter(
          validation_status=="FULL_FROZEN_COVERAGE_CANDIDATE"
        ) |>
        arrange(desc(nrow),file_path)

      if (nrow(fullv)) {
        selected <- fullv$file_path[1]
        decision_note <- "Single-content full frozen-manifest coverage ASV selected."
      }
    }

    if (!is.na(selected)) {

      dst <- file.path(
        CANON,
        "01_PROJECTS",
        project,
        "02_asv_selected",
        basename(selected)
      )

      asv_records[[length(asv_records)+1]] <- link_one(
        src=selected,
        dst=dst,
        category="SELECTED_ASV",
        project=project
      )

      asv_decisions[[length(asv_decisions)+1]] <- tibble(
        project=project,
        selected_asv=selected,
        canonical_target=dst,
        decision="SELECTED_FOR_ANALYSIS",
        rationale=decision_note
      )
    }

    # Keep the best partial ASV as REFERENCE ONLY.
    if (
      route$asv_decision[1]=="PARTIAL_ASV_NOT_SUFFICIENT_FOR_FROZEN_MANIFEST"
    ) {

      partial <- val |>
        arrange(desc(frozen_coverage_pct),file_path) |>
        slice(1)

      if (nrow(partial)) {

        src <- partial$file_path[1]
        dst <- file.path(
          CANON,
          "01_PROJECTS",
          project,
          "03_asv_partial_reference",
          basename(src)
        )

        asv_records[[length(asv_records)+1]] <- link_one(
          src=src,
          dst=dst,
          category="PARTIAL_ASV_REFERENCE_ONLY",
          project=project
        )

        asv_decisions[[length(asv_decisions)+1]] <- tibble(
          project=project,
          selected_asv=src,
          canonical_target=dst,
          decision="REFERENCE_ONLY_NOT_FINAL",
          rationale=paste0(
            "Partial frozen-manifest coverage: ",
            route$best_coverage_pct[1],
            "%. Do not use as final cohort ASV without route-specific decision."
          )
        )
      }
    }
  }

  asv_links <- if (length(asv_records)) {
    bind_rows(asv_records)
  } else {
    tibble()
  }

  asv_decision_df <- if (length(asv_decisions)) {
    bind_rows(asv_decisions)
  } else {
    tibble()
  }

  write_excel_csv(
    asv_links,
    file.path(OUT,"V2_STEP83A_ASV_hardlink_manifest.csv"),
    na=""
  )

  write_excel_csv(
    asv_decision_df,
    file.path(OUT,"V2_STEP83A_ASV_selection_decisions.csv"),
    na=""
  )

  ck("ASV LINKS BUILT")

  # ==========================================================
  # 6. Link taxonomy/QC companions for selected ASVs only
  # ==========================================================

  companion_records <- list()

  if (file.exists(companion_path)) {

    companions <- safe_csv(companion_path)

    if (nrow(companions) && nrow(asv_decision_df)) {

      for (i in seq_len(nrow(asv_decision_df))) {

        project <- asv_decision_df$project[i]

        if (
          asv_decision_df$decision[i]!="SELECTED_FOR_ANALYSIS"
        ) next

        selected_dir <- dirname(
          asv_decision_df$selected_asv[i]
        )

        cands <- companions |>
          filter(
            directory==selected_dir,
            !is.na(file_path)
          )

        if (nrow(cands)) {

          for (j in seq_len(nrow(cands))) {

            src <- cands$file_path[j]

            if (!file_is_regular(src)) next

            dst <- file.path(
              CANON,
              "01_PROJECTS",
              project,
              "04_taxonomy_qc",
              basename(src)
            )

            companion_records[[length(companion_records)+1]] <- link_one(
              src=src,
              dst=dst,
              category="ASV_TAXONOMY_QC_COMPANION",
              project=project
            )
          }
        }
      }
    }
  }

  companion_links <- if (length(companion_records)) {
    bind_rows(companion_records)
  } else {
    tibble()
  }

  write_excel_csv(
    companion_links,
    file.path(OUT,"V2_STEP83A_taxonomy_QC_hardlink_manifest.csv"),
    na=""
  )

  ck("ASV COMPANIONS BUILT")

  # ==========================================================
  # 7. Pending-repair manifests
  # ==========================================================

  if (nrow(missing_raw)) {

    for (project in unique(missing_raw$project)) {

      z <- missing_raw |>
        filter(project==!!project)

      pdir <- file.path(
        CANON,
        "90_PENDING_REPAIR",
        project
      )

      dir.create(pdir,recursive=TRUE,showWarnings=FALSE)

      write_excel_csv(
        z,
        file.path(
          pdir,
          paste0(project,"_missing_raw_runs.csv")
        ),
        na=""
      )
    }
  }

  # Special-route pointer for CRA002354.
  if ("CRA002354" %in% projects) {

    cra_note <- c(
      "CRA002354 SPECIAL SEQUENCE ROUTE",
      "",
      "The frozen cohort contains 131/131 exact CRR Run files.",
      "Local sequence files are FASTA (.fa.gz), not FASTQ.",
      "FASTA has no per-base quality scores and must NOT be sent directly",
      "through the standard FASTQ DADA2 filterAndTrim/error-learning workflow.",
      "",
      "Step83B/84 will define the FASTA-specific feature harmonization route."
    )

    writeLines(
      cra_note,
      file.path(
        CANON,
        "91_SPECIAL_ROUTES",
        "CRA002354_FASTA_ROUTE.txt"
      ),
      useBytes=TRUE
    )
  }

  # ==========================================================
  # 8. Project workspace summary
  # ==========================================================

  raw_summary <- links_raw |>
    group_by(project) |>
    summarise(
      canonical_raw_files=sum(status=="OK"),
      raw_link_failures=sum(status!="OK"),
      canonical_raw_GiB=round(
        sum(source_bytes[status=="OK"],na.rm=TRUE)/1024^3,
        3
      ),
      .groups="drop"
    )

  missing_summary <- if (nrow(missing_raw)) {
    missing_raw |>
      count(project,name="missing_runs_in_canonical")
  } else {
    tibble(
      project=character(),
      missing_runs_in_canonical=integer()
    )
  }

  asv_summary <- if (nrow(asv_decision_df)) {
    asv_decision_df |>
      group_by(project) |>
      summarise(
        selected_final_asv=sum(decision=="SELECTED_FOR_ANALYSIS"),
        partial_reference_asv=sum(decision=="REFERENCE_ONLY_NOT_FINAL"),
        .groups="drop"
      )
  } else {
    tibble(
      project=character(),
      selected_final_asv=integer(),
      partial_reference_asv=integer()
    )
  }

  workspace_summary <- routes |>
    select(
      project,
      analysis_module,
      expected_runs,
      frozen_processing_route,
      immediate_priority
    ) |>
    left_join(raw_summary,by="project") |>
    left_join(missing_summary,by="project") |>
    left_join(asv_summary,by="project") |>
    mutate(
      across(
        c(
          canonical_raw_files,
          raw_link_failures,
          missing_runs_in_canonical,
          selected_final_asv,
          partial_reference_asv
        ),
        ~coalesce(.x,0L)
      )
    )

  write_excel_csv(
    workspace_summary,
    file.path(
      OUT,
      "V2_STEP83A_project_workspace_summary.csv"
    ),
    na=""
  )

  # ==========================================================
  # 9. Canonical file manifest
  # ==========================================================

  all_link_manifest <- bind_rows(
    links_raw,
    asv_links,
    companion_links
  )

  write_excel_csv(
    all_link_manifest,
    file.path(
      CANON,
      "99_MANIFESTS_LOGS",
      "V2_CANONICAL_FILE_LINK_MANIFEST.csv"
    ),
    na=""
  )

  write_excel_csv(
    workspace_summary,
    file.path(
      CANON,
      "99_MANIFESTS_LOGS",
      "V2_CANONICAL_PROJECT_SUMMARY.csv"
    ),
    na=""
  )

  # ==========================================================
  # 10. README
  # ==========================================================

  readme <- c(
    "SEPSIS V2 - CANONICAL ANALYSIS WORKSPACE",
    paste0("Created: ",Sys.time()),
    "",
    "CANONICAL ROOT",
    CANON,
    "",
    "FROM THIS POINT FORWARD",
    "All Step83B+ sequence-processing scripts should read from this canonical root.",
    "Do not perform broad drive scans again unless a specifically missing Run requires repair.",
    "",
    "DATA SAFETY",
    "- Original files were NOT moved or deleted.",
    "- Large sequence/ASV files were exposed here using NTFS hard links.",
    "- A hard link does not duplicate the raw file's disk payload.",
    "- Small frozen metadata/control files were copied into 00_FREEZE.",
    "",
    "DIRECTORY LAYOUT",
    "00_FREEZE/",
    "  Final frozen master, 16S manifest, processing routes.",
    "01_PROJECTS/<PROJECT>/00_metadata/",
    "  Frozen project-specific manifest.",
    "01_PROJECTS/<PROJECT>/01_raw/",
    "  Canonical raw sequence links for locally available frozen Runs.",
    "01_PROJECTS/<PROJECT>/02_asv_selected/",
    "  ASV selected for downstream use, if validated.",
    "01_PROJECTS/<PROJECT>/03_asv_partial_reference/",
    "  Incomplete ASV retained only as reference; not final analysis input.",
    "01_PROJECTS/<PROJECT>/04_taxonomy_qc/",
    "  Companion taxonomy/tracking/QC files for selected ASVs.",
    "01_PROJECTS/<PROJECT>/05_work/",
    "  New Step83/84 DADA2 or sequence-processing outputs.",
    "90_PENDING_REPAIR/",
    "  Exact missing Run manifests.",
    "91_SPECIAL_ROUTES/",
    "  CRA002354 and other non-standard processing notes.",
    "99_MANIFESTS_LOGS/",
    "  Canonical link and project summary manifests.",
    "",
    "PRJNA691455",
    "The rerun_v2 seqtab is selected as the canonical ASV because Step82C showed",
    "49/49 frozen samples, 2,398,535 total reads, median 50,835 reads/sample,",
    "and median richness 148. The legacy full-coverage file had median reads=0",
    "and median richness=0, so it is not used as the canonical analysis ASV.",
    "",
    "NEXT",
    "Step83B: read-length / quality / primer audit from 01_PROJECTS only.",
    "No additional broad filesystem scanning is needed."
  )

  writeLines(
    readme,
    file.path(
      CANON,
      "README_V2_ANALYSIS_READY.txt"
    ),
    useBytes=TRUE
  )

  copy_small(
    file.path(
      CANON,
      "README_V2_ANALYSIS_READY.txt"
    ),
    file.path(
      OUT,
      "README_STEP83A_CANONICAL_WORKSPACE.txt"
    )
  )

  ck("STEP83A COMPLETE")

  cat("\n============================================================\n")
  cat("SEPSIS V2 - STEP83A CANONICAL WORKSPACE COMPLETE\n")
  cat("============================================================\n\n")

  cat("Canonical root:\n")
  cat(CANON,"\n\n")

  print(
    workspace_summary,
    n=Inf,
    width=Inf
  )

  cat("\nMissing canonical Runs: ",nrow(missing_raw),"\n",sep="")
  cat("Sources found outside E:/sepsis_project/data: ",nrow(noncanonical_df),"\n",sep="")
  cat("============================================================\n")
}

tryCatch(
  main(),
  error=function(e) {

    msg <- c(
      paste0("STEP83A FATAL ERROR: ",Sys.time()),
      paste0("Message: ",conditionMessage(e)),
      paste0(
        "Call: ",
        paste(deparse(conditionCall(e)),collapse=" ")
      )
    )

    writeLines(msg,ERR,useBytes=TRUE)
    ck("STEP83A FAILED")
    message(paste(msg,collapse="\n"))

    quit(
      save="no",
      status=1,
      runLast=FALSE
    )
  }
)
