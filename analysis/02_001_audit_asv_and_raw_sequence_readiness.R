# ============================================================
# Sepsis V2 - Step 82
# ASV / RAW sequence readiness audit
#
# Purpose:
#   Audit the FINAL Step81D primary 16S layer before DADA2.
#
# This script does NOT run DADA2.
#
# It determines, for each 16S project:
#   1) expected Runs from the frozen manifest
#   2) local raw FASTQ/FASTA files
#   3) exact Run accession -> local file coverage
#   4) paired/single-file layout
#   5) candidate reusable ASV/feature tables
#   6) recommended next action
#
# Outputs are audit files only.
#
# R 4.4.0 / Windows
# ============================================================

options(stringsAsFactors = FALSE)

pkgs <- c(
  "readr","dplyr","tidyr",
  "stringr","purrr","tibble"
)

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
  library(tidyr)
  library(stringr)
  library(purrr)
  library(tibble)
})

ROOT <- "E:/sepsis_project"

FREEZE_ROOT <- file.path(
  ROOT,
  "results",
  "V2_21D_FINAL_FREEZE_CLINICAL_REPAIRED"
)

OUT <- file.path(
  ROOT,
  "results",
  "V2_22_ASV_RAW_READINESS_AUDIT"
)

dir.create(
  OUT,
  recursive=TRUE,
  showWarnings=FALSE
)

LOG <- file.path(
  OUT,
  "_STEP82_runtime_checkpoints.txt"
)

ERR <- file.path(
  OUT,
  "_STEP82_FATAL_ERROR.txt"
)

if (file.exists(ERR)) unlink(ERR)

ck <- function(x) {
  cat(
    paste0(x, ": ", Sys.time(), "\n"),
    file=LOG,
    append=TRUE
  )
}

clean_chr <- function(x) {
  x <- trimws(as.character(x))
  x[
    is.na(x) |
    x=="" |
    tolower(x) %in% c(
      "na","nan","n/a","null","none"
    )
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

find_latest <- function(folder,pattern) {
  x <- list.files(
    folder,
    pattern=pattern,
    full.names=TRUE
  )

  if (!length(x)) return(NA_character_)

  x[
    which.max(
      file.info(x)$mtime
    )
  ]
}

extract_run <- function(filename) {
  str_extract(
    basename(filename),
    regex(
      "(SRR|ERR|DRR|CRR)[0-9]+",
      ignore_case=TRUE
    )
  )
}

detect_read_mate <- function(filename) {

  x <- basename(filename)

  case_when(
    str_detect(
      x,
      regex(
        "(_R?1_|_R?1\\.|_1\\.(fastq|fq)|_1\\.(fastq|fq)\\.gz$|\\.1\\.(fastq|fq))",
        ignore_case=TRUE
      )
    ) ~ "R1",

    str_detect(
      x,
      regex(
        "(_R?2_|_R?2\\.|_2\\.(fastq|fq)|_2\\.(fastq|fq)\\.gz$|\\.2\\.(fastq|fq))",
        ignore_case=TRUE
      )
    ) ~ "R2",

    str_detect(
      x,
      regex(
        "\\.(fa|fasta)(\\.gz)?$",
        ignore_case=TRUE
      )
    ) ~ "FASTA",

    TRUE ~ "UNSPECIFIED"
  )
}

file_type <- function(filename) {
  x <- basename(filename)

  case_when(
    str_detect(
      x,
      regex("\\.(fastq|fq)\\.gz$",ignore_case=TRUE)
    ) ~ "FASTQ_GZ",

    str_detect(
      x,
      regex("\\.(fa|fasta)\\.gz$",ignore_case=TRUE)
    ) ~ "FASTA_GZ",

    str_detect(
      x,
      regex("\\.(fastq|fq)$",ignore_case=TRUE)
    ) ~ "FASTQ",

    str_detect(
      x,
      regex("\\.(fa|fasta)$",ignore_case=TRUE)
    ) ~ "FASTA",

    TRUE ~ "OTHER"
  )
}

main <- function() {

  ck("STEP82 STARTED")

  manifest_path <- file.path(
    FREEZE_ROOT,
    "V2_FINAL_FROZEN_16S_sequence_manifest.csv"
  )

  if (!file.exists(manifest_path)) {
    stop(
      "Final Step81D 16S sequence manifest was not found."
    )
  }

  manifest <- safe_csv(manifest_path)

  if (nrow(manifest) != 1941) {
    stop(
      paste0(
        "Step82 guard failed: expected 1941 primary 16S rows; observed ",
        nrow(manifest)
      )
    )
  }

  if (
    sum(is.na(clean_chr(manifest$run_id))) != 0
  ) {
    stop(
      "Primary 16S manifest contains missing Run IDs."
    )
  }

  if (
    n_distinct(manifest$run_id) != 1941
  ) {
    stop(
      "Primary 16S manifest contains duplicated Run IDs."
    )
  }

  projects <- sort(
    unique(
      manifest$project
    )
  )

  ck("FINAL 16S MANIFEST GUARDED")

  # ==========================================================
  # 1. Raw sequence file inventory
  # ==========================================================

  raw_rows <- list()

  for (project in projects) {

    roots <- c(
      file.path(
        ROOT,
        "data",
        project
      )
    )

    # Search other known data roots if they exist.
    extra_roots <- c(
      "F:/sepsis_V2",
      "F:/sepsis",
      "F:/sepsis_V2/02_raw_data"
    )

    for (r in extra_roots) {
      p <- file.path(r,project)
      if (dir.exists(p)) {
        roots <- c(roots,p)
      }
    }

    roots <- unique(
      roots[
        dir.exists(roots)
      ]
    )

    files <- character()

    for (r in roots) {
      files <- c(
        files,
        list.files(
          r,
          recursive=TRUE,
          full.names=TRUE
        )
      )
    }

    files <- unique(files)

    seq_files <- files[
      str_detect(
        basename(files),
        regex(
          "\\.(fastq|fq|fa|fasta)(\\.gz)?$",
          ignore_case=TRUE
        )
      )
    ]

    if (length(seq_files)) {

      info <- file.info(seq_files)

      raw_rows[[length(raw_rows)+1]] <- tibble(
        project=project,
        file_path=seq_files,
        file_name=basename(seq_files),
        bytes=as.numeric(info$size),
        modified=as.character(info$mtime),
        run_id_from_filename=extract_run(seq_files),
        read_mate=map_chr(
          seq_files,
          detect_read_mate
        ),
        sequence_file_type=map_chr(
          seq_files,
          file_type
        )
      )
    }
  }

  raw_inventory <- if (!length(raw_rows)) {
    tibble(
      project=character(),
      file_path=character(),
      file_name=character(),
      bytes=double(),
      modified=character(),
      run_id_from_filename=character(),
      read_mate=character(),
      sequence_file_type=character()
    )
  } else {
    bind_rows(raw_rows)
  }

  write_excel_csv(
    raw_inventory,
    file.path(
      OUT,
      "V2_STEP82_raw_sequence_file_inventory.csv"
    ),
    na=""
  )

  ck("RAW FILE INVENTORY COMPLETE")

  # ==========================================================
  # 2. Exact expected Run -> local file coverage
  # ==========================================================

  expected_runs <- manifest |>
    select(
      project,
      patient_id,
      sample_id,
      run_id,
      time_raw,
      time_class,
      analysis_module
    ) |>
    distinct()

  raw_by_run <- raw_inventory |>
    filter(
      !is.na(run_id_from_filename)
    ) |>
    group_by(
      project,
      run_id_from_filename
    ) |>
    summarise(
      local_file_n=n(),
      local_bytes=sum(bytes,na.rm=TRUE),

      R1_n=sum(
        read_mate=="R1",
        na.rm=TRUE
      ),

      R2_n=sum(
        read_mate=="R2",
        na.rm=TRUE
      ),

      FASTA_n=sum(
        read_mate=="FASTA",
        na.rm=TRUE
      ),

      unspecified_n=sum(
        read_mate=="UNSPECIFIED",
        na.rm=TRUE
      ),

      file_types=paste(
        sort(unique(sequence_file_type)),
        collapse=";"
      ),

      file_paths=paste(
        file_path,
        collapse=" | "
      ),

      .groups="drop"
    )

  run_coverage <- expected_runs |>
    left_join(
      raw_by_run,
      by=c(
        "project",
        "run_id"="run_id_from_filename"
      )
    ) |>
    mutate(
      local_file_n=coalesce(
        local_file_n,
        0L
      ),

      R1_n=coalesce(
        R1_n,
        0L
      ),

      R2_n=coalesce(
        R2_n,
        0L
      ),

      FASTA_n=coalesce(
        FASTA_n,
        0L
      ),

      unspecified_n=coalesce(
        unspecified_n,
        0L
      ),

      raw_run_covered=
        local_file_n>0,

      local_layout=case_when(
        R1_n>=1 & R2_n>=1 ~
          "PAIRED_FASTQ",

        FASTA_n>=1 ~
          "FASTA_RUN_FILE",

        local_file_n>=1 ~
          "SINGLE_OR_UNSPECIFIED",

        TRUE ~
          "NO_LOCAL_RAW"
      )
    )

  write_excel_csv(
    run_coverage,
    file.path(
      OUT,
      "V2_STEP82_expected_run_local_coverage.csv"
    ),
    na=""
  )

  ck("RUN COVERAGE AUDIT COMPLETE")

  # ==========================================================
  # 3. Candidate ASV / feature-table inventory
  # ==========================================================

  search_roots <- c(
    file.path(ROOT,"data"),
    file.path(ROOT,"results")
  )

  all_candidate_files <- character()

  for (r in search_roots) {

    if (!dir.exists(r)) next

    all_candidate_files <- c(
      all_candidate_files,
      list.files(
        r,
        recursive=TRUE,
        full.names=TRUE
      )
    )
  }

  all_candidate_files <- unique(
    all_candidate_files
  )

  # Exclude current audit output and obvious metadata-only files.
  asv_files <- all_candidate_files[
    str_detect(
      basename(all_candidate_files),
      regex(
        "asv|seqtab|feature[_ -]?table|otu[_ -]?table|taxa[_ -]?table|taxonomy|phyloseq",
        ignore_case=TRUE
      )
    ) &
    str_detect(
      basename(all_candidate_files),
      regex(
        "\\.(csv|tsv|txt|rds|rda|rdata|xlsx)$",
        ignore_case=TRUE
      )
    )
  ]

  asv_files <- asv_files[
    !str_detect(
      asv_files,
      regex(
        "V2_22_ASV_RAW_READINESS_AUDIT",
        ignore_case=TRUE
      )
    )
  ]

  asv_inventory <- tibble(
    file_path=asv_files,
    file_name=basename(asv_files),
    bytes=if (length(asv_files)) {
      as.numeric(
        file.info(asv_files)$size
      )
    } else {
      numeric()
    },

    modified=if (length(asv_files)) {
      as.character(
        file.info(asv_files)$mtime
      )
    } else {
      character()
    }
  )

  # Assign project by path/name only when an exact project ID
  # is present. Do not infer otherwise.
  if (nrow(asv_inventory)) {

    asv_inventory <- asv_inventory |>
      rowwise() |>
      mutate(
        project= {
          hits <- projects[
            str_detect(
              file_path,
              fixed(
                projects,
                ignore_case=TRUE
              )
            )
          ]

          if (length(hits)==1) {
            hits
          } else {
            NA_character_
          }
        }
      ) |>
      ungroup()
  } else {
    asv_inventory$project <- character()
  }

  write_excel_csv(
    asv_inventory,
    file.path(
      OUT,
      "V2_STEP82_ASV_candidate_file_inventory.csv"
    ),
    na=""
  )

  ck("ASV CANDIDATE INVENTORY COMPLETE")

  # ==========================================================
  # 4. Per-project readiness
  # ==========================================================

  raw_project <- run_coverage |>
    group_by(project) |>
    summarise(
      expected_runs=n(),

      runs_with_local_raw=sum(
        raw_run_covered
      ),

      runs_missing_local_raw=sum(
        !raw_run_covered
      ),

      raw_coverage_pct=round(
        100*runs_with_local_raw /
        expected_runs,
        2
      ),

      paired_fastq_runs=sum(
        local_layout=="PAIRED_FASTQ"
      ),

      fasta_run_files=sum(
        local_layout=="FASTA_RUN_FILE"
      ),

      unspecified_raw_runs=sum(
        local_layout=="SINGLE_OR_UNSPECIFIED"
      ),

      total_local_sequence_files=sum(
        local_file_n
      ),

      total_local_raw_GiB=round(
        sum(
          local_bytes,
          na.rm=TRUE
        ) /
        1024^3,
        3
      ),

      .groups="drop"
    )

  asv_project <- tibble(
    project=projects
  ) |>
    left_join(
      asv_inventory |>
        filter(
          !is.na(project)
        ) |>
        count(
          project,
          name="candidate_asv_files"
        ),
      by="project"
    ) |>
    mutate(
      candidate_asv_files=coalesce(
        candidate_asv_files,
        0L
      )
    )

  readiness <- raw_project |>
    left_join(
      asv_project,
      by="project"
    ) |>
    left_join(
      manifest |>
        group_by(
          project,
          analysis_module
        ) |>
        summarise(
          patients=n_distinct(
            patient_id,
            na.rm=TRUE
          ),

          sequence_rows=n(),

          .groups="drop"
        ),
      by="project"
    ) |>
    mutate(
      readiness_class=case_when(

        candidate_asv_files>0 &
        raw_coverage_pct==100 ~
          "ASV_CANDIDATE_FOUND_AND_RAW_COMPLETE",

        candidate_asv_files>0 &
        raw_coverage_pct<100 ~
          "ASV_CANDIDATE_FOUND_RAW_INCOMPLETE_OR_NOT_LOCAL",

        candidate_asv_files==0 &
        raw_coverage_pct==100 ~
          "RAW_COMPLETE_DADA2_REQUIRED",

        candidate_asv_files==0 &
        raw_coverage_pct>0 &
        raw_coverage_pct<100 ~
          "RAW_PARTIAL_REPAIR_REQUIRED",

        candidate_asv_files==0 &
        raw_coverage_pct==0 ~
          "NO_LOCAL_RAW_OR_ASV_CANDIDATE",

        TRUE ~
          "MANUAL_REVIEW"
      ),

      recommended_next_action=case_when(

        candidate_asv_files>0 ~
          "Step82B: inspect ASV candidate schema/sample coverage before deciding whether DADA2 can be skipped",

        raw_coverage_pct==100 ~
          "Step83: validate read layout/amplicon region and prepare cohort-specific DADA2",

        raw_coverage_pct>0 &
        raw_coverage_pct<100 ~
          "Repair missing raw Runs before DADA2",

        raw_coverage_pct==0 ~
          "Check known existing ASV locations or download raw sequence files",

        TRUE ~
          "Manual review"
      )
    ) |>
    arrange(
      analysis_module,
      project
    )

  write_excel_csv(
    readiness,
    file.path(
      OUT,
      "V2_STEP82_project_sequence_readiness.csv"
    ),
    na=""
  )

  # ==========================================================
  # 5. Action queues
  # ==========================================================

  asv_queue <- readiness |>
    filter(
      candidate_asv_files>0
    )

  dada2_queue <- readiness |>
    filter(
      candidate_asv_files==0,
      raw_coverage_pct==100
    )

  repair_queue <- readiness |>
    filter(
      raw_coverage_pct<100
    )

  write_excel_csv(
    asv_queue,
    file.path(
      OUT,
      "V2_STEP82_queue_ASV_validation.csv"
    ),
    na=""
  )

  write_excel_csv(
    dada2_queue,
    file.path(
      OUT,
      "V2_STEP82_queue_DADA2_ready.csv"
    ),
    na=""
  )

  write_excel_csv(
    repair_queue,
    file.path(
      OUT,
      "V2_STEP82_queue_raw_repair.csv"
    ),
    na=""
  )

  # ==========================================================
  # 6. Summary
  # ==========================================================

  summary <- tibble(
    metric=c(
      "Primary_16S_projects",
      "Primary_16S_rows",
      "Expected_unique_Runs",
      "Runs_with_local_raw",
      "Runs_missing_local_raw",
      "Projects_with_ASV_candidate",
      "Projects_raw_100pct",
      "Projects_raw_partial",
      "Projects_raw_zero"
    ),

    value=c(
      length(projects),
      nrow(manifest),
      n_distinct(manifest$run_id),
      sum(run_coverage$raw_run_covered),
      sum(!run_coverage$raw_run_covered),
      sum(readiness$candidate_asv_files>0),
      sum(readiness$raw_coverage_pct==100),
      sum(
        readiness$raw_coverage_pct>0 &
        readiness$raw_coverage_pct<100
      ),
      sum(readiness$raw_coverage_pct==0)
    )
  )

  write_excel_csv(
    summary,
    file.path(
      OUT,
      "V2_STEP82_global_readiness_summary.csv"
    ),
    na=""
  )

  # ==========================================================
  # 7. README
  # ==========================================================

  readme <- c(
    "SEPSIS V2 - STEP82 ASV / RAW READINESS AUDIT",
    paste0("Created: ",Sys.time()),
    "",
    "THIS STEP DOES NOT RUN DADA2.",
    "",
    "INPUT",
    manifest_path,
    "",
    "PURPOSE",
    "1. Confirm which expected frozen 16S Runs are already present locally.",
    "2. Identify paired FASTQ / FASTA / unspecified layouts.",
    "3. Identify candidate reusable ASV/feature tables.",
    "4. Build separate queues for ASV validation, DADA2-ready cohorts, and raw-data repair.",
    "",
    "IMPORTANT INTERPRETATION",
    "- A file is only an ASV CANDIDATE because its filename/path contains ASV/seqtab/feature/OTU/taxonomy terms.",
    "- Step82 does NOT certify that a candidate ASV table is correct or complete.",
    "- Exact Run accession in the local raw filename is required for automatic raw coverage.",
    "- Raw files without accession in filename remain visible in raw inventory but are not automatically assigned to an expected Run.",
    "",
    "NEXT",
    "Step82B should inspect candidate ASV tables for dimensions, sample identifiers, taxonomy pairing, and frozen-manifest coverage.",
    "Only cohorts without a reusable ASV should proceed to cohort-specific DADA2."
  )

  writeLines(
    readme,
    file.path(
      OUT,
      "README_STEP82_ASV_RAW_READINESS.txt"
    ),
    useBytes=TRUE
  )

  ck("STEP82 COMPLETE")

  cat("\n============================================================\n")
  cat("SEPSIS V2 - STEP82 READINESS AUDIT COMPLETE\n")
  cat("============================================================\n\n")

  print(
    summary,
    n=Inf,
    width=Inf
  )

  cat("\nPROJECT READINESS:\n")
  print(
    readiness,
    n=Inf,
    width=Inf
  )

  cat("\nOutput:\n")
  cat(OUT,"\n")
  cat("============================================================\n")
}

tryCatch(
  main(),
  error=function(e) {

    msg <- c(
      paste0(
        "STEP82 FATAL ERROR: ",
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

    ck("STEP82 FAILED")

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
