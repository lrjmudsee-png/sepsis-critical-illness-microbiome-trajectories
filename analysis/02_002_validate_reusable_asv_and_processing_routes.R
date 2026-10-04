# ============================================================
# Sepsis V2 - Step 82B
# Validate reusable ASV/seqtab candidates against FINAL frozen
# 16S manifest + classify non-ASV raw processing routes.
#
# DOES NOT RUN DADA2.
#
# Main questions:
# 1. Does each seqtab candidate really contain the frozen samples?
# 2. Are sample IDs on rows or columns?
# 3. How many frozen Runs/samples are represented?
# 4. Are candidate files duplicates?
# 5. Which projects can reuse an ASV and skip raw repair?
# 6. For projects without reusable ASV, is raw layout:
#      paired FASTQ / single FASTQ / FASTA without qualities?
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

STEP82_ROOT <- file.path(
  ROOT,
  "results",
  "V2_22_ASV_RAW_READINESS_AUDIT"
)

OUT <- file.path(
  ROOT,
  "results",
  "V2_22B_ASV_VALIDATION"
)

dir.create(
  OUT,
  recursive=TRUE,
  showWarnings=FALSE
)

LOG <- file.path(
  OUT,
  "_STEP82B_runtime_checkpoints.txt"
)

ERR <- file.path(
  OUT,
  "_STEP82B_FATAL_ERROR.txt"
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

norm_id <- function(x) {
  x <- toupper(clean_chr(x))
  str_replace_all(x,"[^A-Z0-9]","")
}

extract_run <- function(x) {
  str_extract(
    clean_chr(x),
    regex(
      "(SRR|ERR|DRR|CRR)[0-9]+",
      ignore_case=TRUE
    )
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

is_dna_sequence <- function(x) {
  x <- clean_chr(x)

  !is.na(x) &
  nchar(x)>=20 &
  str_detect(
    x,
    regex(
      "^[ACGTN]+$",
      ignore_case=TRUE
    )
  )
}

safe_md5 <- function(path) {
  tryCatch(
    unname(
      tools::md5sum(path)
    ),
    error=function(e) NA_character_
  )
}

read_rds_summary <- function(
  path,
  project,
  expected
) {

  obj <- tryCatch(
    readRDS(path),
    error=function(e) e
  )

  if (inherits(obj,"error")) {
    return(list(
      summary=tibble(
        project=project,
        file_path=path,
        file_name=basename(path),
        md5=safe_md5(path),
        object_class="READ_ERROR",
        nrow=NA_integer_,
        ncol=NA_integer_,
        sample_axis=NA_character_,
        matched_expected_runs=0L,
        matched_expected_samples=0L,
        expected_runs=n_distinct(expected$run_id),
        frozen_coverage_pct=0,
        feature_axis_dna_pct=NA_real_,
        duplicate_sample_keys=NA_integer_,
        validation_status="READ_ERROR",
        error_message=conditionMessage(obj)
      ),
      coverage=tibble()
    ))
  }

  if (
    !is.matrix(obj) &&
    !is.data.frame(obj)
  ) {
    return(list(
      summary=tibble(
        project=project,
        file_path=path,
        file_name=basename(path),
        md5=safe_md5(path),
        object_class=paste(class(obj),collapse=";"),
        nrow=ifelse(
          is.null(dim(obj)),
          NA_integer_,
          dim(obj)[1]
        ),
        ncol=ifelse(
          is.null(dim(obj)),
          NA_integer_,
          dim(obj)[2]
        ),
        sample_axis=NA_character_,
        matched_expected_runs=0L,
        matched_expected_samples=0L,
        expected_runs=n_distinct(expected$run_id),
        frozen_coverage_pct=0,
        feature_axis_dna_pct=NA_real_,
        duplicate_sample_keys=NA_integer_,
        validation_status="NOT_MATRIX_OR_DATAFRAME",
        error_message=NA_character_
      ),
      coverage=tibble()
    ))
  }

  nr <- nrow(obj)
  nc <- ncol(obj)

  rn <- rownames(obj)
  cn <- colnames(obj)

  if (is.null(rn)) rn <- rep(NA_character_,nr)
  if (is.null(cn)) cn <- rep(NA_character_,nc)

  expected_runs <- clean_chr(expected$run_id)
  expected_samples <- clean_chr(expected$sample_id)
  expected_biosamples <- clean_chr(expected$biosample)

  expected_run_norm <- norm_id(expected_runs)
  expected_sample_norm <- norm_id(expected_samples)
  expected_biosample_norm <- norm_id(expected_biosamples)

  score_axis <- function(keys) {

    keys_chr <- clean_chr(keys)
    keys_run <- extract_run(keys_chr)
    keys_norm <- norm_id(keys_chr)

    run_overlap <- sum(
      !is.na(keys_run) &
      keys_run %in% expected_runs
    )

    sample_overlap <- sum(
      keys_norm %in%
      expected_sample_norm,
      na.rm=TRUE
    )

    biosample_overlap <- sum(
      keys_norm %in%
      expected_biosample_norm,
      na.rm=TRUE
    )

    c(
      run=run_overlap,
      sample=sample_overlap,
      biosample=biosample_overlap,
      total=max(
        run_overlap,
        sample_overlap,
        biosample_overlap
      )
    )
  }

  row_score <- score_axis(rn)
  col_score <- score_axis(cn)

  # DADA2 seqtabs should normally have samples on rows.
  # Still detect orientation instead of assuming.
  if (
    row_score["total"] >
    col_score["total"]
  ) {
    sample_axis <- "ROWS"
    sample_keys <- rn
    feature_keys <- cn
  } else if (
    col_score["total"] >
    row_score["total"]
  ) {
    sample_axis <- "COLUMNS"
    sample_keys <- cn
    feature_keys <- rn
  } else if (
    row_score["total"]>0
  ) {
    # Tie: prefer rows because DADA2 seqtab convention is rows=samples.
    sample_axis <- "ROWS_TIE_PREFERRED"
    sample_keys <- rn
    feature_keys <- cn
  } else {
    sample_axis <- "UNRESOLVED"
    sample_keys <- rn
    feature_keys <- cn
  }

  key_chr <- clean_chr(sample_keys)
  key_run <- extract_run(key_chr)
  key_norm <- norm_id(key_chr)

  # Match each frozen row using strict hierarchy:
  # 1. Run accession extracted from seqtab key
  # 2. normalized exact sample_id
  # 3. normalized exact biosample
  coverage <- expected |>
    transmute(
      project,
      patient_id,
      sample_id,
      biosample,
      run_id,
      time_raw,
      frozen_run_norm=norm_id(run_id),
      frozen_sample_norm=norm_id(sample_id),
      frozen_biosample_norm=norm_id(biosample)
    ) |>
    rowwise() |>
    mutate(
      matched_by_run=
        any(
          !is.na(key_run) &
          key_run==run_id
        ),

      matched_by_sample=
        any(
          !is.na(key_norm) &
          key_norm==frozen_sample_norm
        ),

      matched_by_biosample=
        any(
          !is.na(key_norm) &
          !is.na(frozen_biosample_norm) &
          key_norm==frozen_biosample_norm
        ),

      asv_sample_present=
        matched_by_run |
        matched_by_sample |
        matched_by_biosample,

      match_method=case_when(
        matched_by_run ~ "RUN_ACCESSION",
        matched_by_sample ~ "SAMPLE_ID_EXACT_NORMALIZED",
        matched_by_biosample ~ "BIOSAMPLE_EXACT_NORMALIZED",
        TRUE ~ "UNMATCHED"
      )
    ) |>
    ungroup()

  matched_runs <- sum(
    coverage$matched_by_run,
    na.rm=TRUE
  )

  matched_samples <- sum(
    coverage$asv_sample_present,
    na.rm=TRUE
  )

  expected_n <- n_distinct(
    expected$run_id
  )

  cov_pct <- round(
    100*matched_samples /
    expected_n,
    2
  )

  dna_pct <- if (
    length(feature_keys)
  ) {
    round(
      100*mean(
        is_dna_sequence(feature_keys),
        na.rm=TRUE
      ),
      2
    )
  } else {
    NA_real_
  }

  key_identity <- ifelse(
    !is.na(key_run),
    key_run,
    key_norm
  )

  dup_keys <- sum(
    duplicated(
      key_identity[
        !is.na(key_identity)
      ]
    )
  )

  validation <- case_when(
    sample_axis=="UNRESOLVED" ~
      "SAMPLE_AXIS_UNRESOLVED",

    matched_samples==expected_n &
    dup_keys==0 ~
      "FULL_FROZEN_COVERAGE_CANDIDATE",

    matched_samples>0 ~
      "PARTIAL_FROZEN_COVERAGE",

    TRUE ~
      "NO_FROZEN_ID_MATCH"
  )

  summary <- tibble(
    project=project,
    file_path=path,
    file_name=basename(path),
    md5=safe_md5(path),
    object_class=paste(class(obj),collapse=";"),
    nrow=nr,
    ncol=nc,
    sample_axis=sample_axis,
    matched_expected_runs=matched_runs,
    matched_expected_samples=matched_samples,
    expected_runs=expected_n,
    frozen_coverage_pct=cov_pct,
    feature_axis_dna_pct=dna_pct,
    duplicate_sample_keys=dup_keys,
    validation_status=validation,
    error_message=NA_character_
  )

  list(
    summary=summary,
    coverage=coverage
  )
}

main <- function() {

  ck("STEP82B STARTED")

  manifest_path <- file.path(
    FREEZE_ROOT,
    "V2_FINAL_FROZEN_16S_sequence_manifest.csv"
  )

  asv_inventory_path <- file.path(
    STEP82_ROOT,
    "V2_STEP82_ASV_candidate_file_inventory.csv"
  )

  raw_inventory_path <- file.path(
    STEP82_ROOT,
    "V2_STEP82_raw_sequence_file_inventory.csv"
  )

  readiness_path <- file.path(
    STEP82_ROOT,
    "V2_STEP82_project_sequence_readiness.csv"
  )

  needed <- c(
    manifest_path,
    asv_inventory_path,
    raw_inventory_path,
    readiness_path
  )

  if (!all(file.exists(needed))) {
    stop(
      "One or more Step82/Step81D input files are missing."
    )
  }

  manifest <- safe_csv(manifest_path)
  asv_inv <- safe_csv(asv_inventory_path)
  raw_inv <- safe_csv(raw_inventory_path)
  readiness82 <- safe_csv(readiness_path)

  if (
    nrow(manifest)!=1941 ||
    n_distinct(manifest$run_id)!=1941
  ) {
    stop(
      "Final frozen 16S manifest guard failed."
    )
  }

  projects <- sort(
    unique(manifest$project)
  )

  ck("INPUTS GUARDED")

  # ==========================================================
  # 1. Select actual seqtab / ASV abundance candidates.
  #
  # Exclude taxonomy, prechim, merged and generic inventory.
  # Prefer final abundance tables.
  # ==========================================================

  abundance_candidates <- asv_inv |>
    filter(
      !is.na(project),
      project %in% projects,
      str_detect(
        file_name,
        regex(
          "seqtab.*final|asv[_ -]?table|feature[_ -]?table|otu[_ -]?table",
          ignore_case=TRUE
        )
      ),
      !str_detect(
        file_name,
        regex(
          "taxonomy|prechim|pre[_ -]?chimera|merged",
          ignore_case=TRUE
        )
      )
    ) |>
    mutate(
      extension=tolower(
        tools::file_ext(file_path)
      ),
      md5=map_chr(
        file_path,
        safe_md5
      )
    )

  write_excel_csv(
    abundance_candidates,
    file.path(
      OUT,
      "V2_STEP82B_abundance_table_candidates.csv"
    ),
    na=""
  )

  # ==========================================================
  # 2. Validate RDS seqtabs
  # ==========================================================

  validation_results <- list()
  coverage_results <- list()

  rds_candidates <- abundance_candidates |>
    filter(
      extension=="rds"
    )

  for (i in seq_len(nrow(rds_candidates))) {

    project <- rds_candidates$project[i]
    path <- rds_candidates$file_path[i]

    expected <- manifest |>
      filter(
        project==!!project
      )

    z <- read_rds_summary(
      path,
      project,
      expected
    )

    validation_results[[
      length(validation_results)+1
    ]] <- z$summary

    if (
      nrow(z$coverage)
    ) {
      coverage_results[[
        length(coverage_results)+1
      ]] <- z$coverage |>
        mutate(
          candidate_file=path,
          .before=1
        )
    }
  }

  validation <- if (
    length(validation_results)
  ) {
    bind_rows(
      validation_results
    )
  } else {
    tibble()
  }

  coverage_all <- if (
    length(coverage_results)
  ) {
    bind_rows(
      coverage_results
    )
  } else {
    tibble()
  }

  write_excel_csv(
    validation,
    file.path(
      OUT,
      "V2_STEP82B_seqtab_validation_all.csv"
    ),
    na=""
  )

  write_excel_csv(
    coverage_all,
    file.path(
      OUT,
      "V2_STEP82B_seqtab_frozen_sample_coverage_all.csv"
    ),
    na=""
  )

  ck("RDS SEQTAB VALIDATION COMPLETE")

  # ==========================================================
  # 3. De-duplicate identical candidates by MD5
  # ==========================================================

  unique_validation <- validation

  if (
    nrow(unique_validation)
  ) {

    unique_validation <- unique_validation |>
      mutate(
        path_priority=case_when(
          str_detect(
            file_path,
            regex(
              "rerun_v2",
              ignore_case=TRUE
            )
          ) ~ 1L,

          str_detect(
            file_path,
            regex(
              "/03_dada2/|\\\\03_dada2\\\\",
              ignore_case=TRUE
            )
          ) ~ 2L,

          TRUE ~ 3L
        )
      ) |>
      arrange(
        project,
        md5,
        path_priority,
        file_path
      ) |>
      group_by(
        project,
        md5
      ) |>
      mutate(
        identical_copy_n=n(),
        identical_copy_rank=row_number()
      ) |>
      ungroup()

    write_excel_csv(
      unique_validation,
      file.path(
        OUT,
        "V2_STEP82B_seqtab_validation_with_MD5_duplicates.csv"
      ),
      na=""
    )
  }

  # ==========================================================
  # 4. Project-level best evidence
  #
  # This does NOT silently approve a seqtab if several
  # non-identical full-coverage versions exist.
  # ==========================================================

  project_asv <- tibble(
    project=projects
  )

  if (
    nrow(validation)
  ) {

    val2 <- validation |>
      mutate(
        full_coverage=
          validation_status==
          "FULL_FROZEN_COVERAGE_CANDIDATE",

        path_priority=case_when(
          str_detect(
            file_path,
            regex(
              "rerun_v2",
              ignore_case=TRUE
            )
          ) ~ 1L,

          str_detect(
            file_path,
            regex(
              "/03_dada2/|\\\\03_dada2\\\\",
              ignore_case=TRUE
            )
          ) ~ 2L,

          TRUE ~ 3L
        )
      )

    project_asv <- project_asv |>
      left_join(
        val2 |>
          group_by(project) |>
          summarise(
            validated_seqtab_candidates=n(),

            full_coverage_candidates=sum(
              full_coverage,
              na.rm=TRUE
            ),

            unique_full_coverage_MD5=n_distinct(
              md5[
                full_coverage
              ],
              na.rm=TRUE
            ),

            best_coverage_pct=max(
              frozen_coverage_pct,
              na.rm=TRUE
            ),

            best_candidate_file={
              x <- cur_data_all()

              ord <- order(
                -x$frozen_coverage_pct,
                x$path_priority,
                -x$nrow,
                x$file_path
              )

              x$file_path[
                ord[1]
              ]
            },

            best_candidate_status={
              x <- cur_data_all()

              ord <- order(
                -x$frozen_coverage_pct,
                x$path_priority,
                -x$nrow,
                x$file_path
              )

              x$validation_status[
                ord[1]
              ]
            },

            .groups="drop"
          ),
        by="project"
      )
  }

  project_asv <- project_asv |>
    mutate(
      across(
        c(
          validated_seqtab_candidates,
          full_coverage_candidates,
          unique_full_coverage_MD5
        ),
        ~coalesce(.x,0L)
      ),

      best_coverage_pct=ifelse(
        is.na(best_coverage_pct),
        0,
        best_coverage_pct
      ),

      asv_decision=case_when(
        full_coverage_candidates>=1 &
        unique_full_coverage_MD5==1 ~
          "FULL_COVERAGE_SINGLE_CONTENT_ASV_REUSABLE_CANDIDATE",

        full_coverage_candidates>=2 &
        unique_full_coverage_MD5>1 ~
          "MULTIPLE_NONIDENTICAL_FULL_COVERAGE_ASV_REVIEW",

        best_coverage_pct>0 &
        best_coverage_pct<100 ~
          "PARTIAL_ASV_NOT_SUFFICIENT_FOR_FROZEN_MANIFEST",

        TRUE ~
          "NO_VALIDATED_ASV"
      )
    )

  # ==========================================================
  # 5. Raw-layout classification
  # ==========================================================

  raw_project <- raw_inv |>
    group_by(project) |>
    summarise(
      local_seq_files=n(),

      fastq_gz_n=sum(
        sequence_file_type=="FASTQ_GZ"
      ),

      fasta_gz_n=sum(
        sequence_file_type=="FASTA_GZ"
      ),

      R1_files=sum(
        read_mate=="R1"
      ),

      R2_files=sum(
        read_mate=="R2"
      ),

      fasta_files=sum(
        read_mate=="FASTA"
      ),

      unspecified_files=sum(
        read_mate=="UNSPECIFIED"
      ),

      total_GiB=round(
        sum(bytes,na.rm=TRUE)/
        1024^3,
        3
      ),

      .groups="drop"
    )

  final_decision <- readiness82 |>
    select(
      project,
      expected_runs,
      runs_with_local_raw,
      runs_missing_local_raw,
      raw_coverage_pct,
      paired_fastq_runs,
      fasta_run_files,
      unspecified_raw_runs,
      analysis_module
    ) |>
    left_join(
      project_asv,
      by="project"
    ) |>
    left_join(
      raw_project,
      by="project"
    ) |>
    mutate(
      raw_processing_route=case_when(

        fasta_run_files==expected_runs &
        expected_runs>0 ~
          "FASTA_NO_QUALITY_SPECIAL_ROUTE_NOT_STANDARD_DADA2",

        paired_fastq_runs==expected_runs &
        expected_runs>0 ~
          "PAIRED_END_FASTQ_DADA2_ROUTE",

        unspecified_raw_runs==expected_runs &
        expected_runs>0 &
        fastq_gz_n>=expected_runs ~
          "ONE_FASTQ_PER_RUN_LIKELY_SINGLE_END_VALIDATE_THEN_DADA2",

        raw_coverage_pct==100 ~
          "RAW_COMPLETE_LAYOUT_MIXED_REVIEW",

        raw_coverage_pct>0 &
        raw_coverage_pct<100 ~
          "RAW_PARTIAL",

        TRUE ~
          "NO_LOCAL_RAW"
      ),

      final_next_action=case_when(

        asv_decision==
        "FULL_COVERAGE_SINGLE_CONTENT_ASV_REUSABLE_CANDIDATE" ~
          "REUSE_ASV_AFTER_FINAL_FEATURE/TAXONOMY_PROVENANCE_CHECK; DO_NOT REPAIR RAW YET",

        asv_decision==
        "MULTIPLE_NONIDENTICAL_FULL_COVERAGE_ASV_REVIEW" ~
          "COMPARE FULL-COVERAGE ASV VERSIONS AND SELECT DOCUMENTED PIPELINE VERSION",

        asv_decision==
        "PARTIAL_ASV_NOT_SUFFICIENT_FOR_FROZEN_MANIFEST" ~
          "ASV PARTIAL; RAW PROCESSING/REPAIR STILL REQUIRED",

        raw_processing_route==
        "FASTA_NO_QUALITY_SPECIAL_ROUTE_NOT_STANDARD_DADA2" ~
          "DO NOT RUN STANDARD DADA2; DEFINE FASTA-SPECIFIC HARMONIZATION/FEATURE ROUTE",

        raw_processing_route==
        "PAIRED_END_FASTQ_DADA2_ROUTE" ~
          "READY FOR COHORT-SPECIFIC PAIRED-END DADA2",

        raw_processing_route==
        "ONE_FASTQ_PER_RUN_LIKELY_SINGLE_END_VALIDATE_THEN_DADA2" ~
          "VALIDATE FASTQ CONTENT/PRIMER REGION THEN RUN SINGLE-END DADA2",

        raw_processing_route==
        "RAW_PARTIAL" ~
          "REPAIR MISSING RAW RUNS BEFORE DADA2",

        TRUE ~
          "MANUAL REVIEW"
      )
    ) |>
    arrange(
      analysis_module,
      project
    )

  write_excel_csv(
    final_decision,
    file.path(
      OUT,
      "V2_STEP82B_FINAL_project_processing_decision.csv"
    ),
    na=""
  )

  # ==========================================================
  # 6. Missing frozen samples for the BEST candidate per project
  # ==========================================================

  missing_best <- list()

  for (project in projects) {

    best_file <- final_decision$best_candidate_file[
      final_decision$project==project
    ]

    if (
      length(best_file)==1 &&
      !is.na(best_file) &&
      nrow(coverage_all)
    ) {

      z <- coverage_all |>
        filter(
          project==!!project,
          candidate_file==best_file,
          !asv_sample_present
        )

      if (nrow(z)) {
        missing_best[[
          length(missing_best)+1
        ]] <- z
      }
    }
  }

  missing_best_df <- if (
    length(missing_best)
  ) {
    bind_rows(missing_best)
  } else {
    tibble()
  }

  write_excel_csv(
    missing_best_df,
    file.path(
      OUT,
      "V2_STEP82B_missing_frozen_samples_in_best_ASV.csv"
    ),
    na=""
  )

  # ==========================================================
  # 7. Processing queues
  # ==========================================================

  reuse_queue <- final_decision |>
    filter(
      str_detect(
        asv_decision,
        "^FULL_COVERAGE|^MULTIPLE_NONIDENTICAL_FULL_COVERAGE"
      )
    )

  dada2_queue <- final_decision |>
    filter(
      asv_decision=="NO_VALIDATED_ASV",
      str_detect(
        raw_processing_route,
        "DADA2"
      )
    )

  special_queue <- final_decision |>
    filter(
      raw_processing_route==
      "FASTA_NO_QUALITY_SPECIAL_ROUTE_NOT_STANDARD_DADA2"
    )

  repair_queue <- final_decision |>
    filter(
      asv_decision %in% c(
        "NO_VALIDATED_ASV",
        "PARTIAL_ASV_NOT_SUFFICIENT_FOR_FROZEN_MANIFEST"
      ),
      raw_processing_route=="RAW_PARTIAL"
    )

  write_excel_csv(
    reuse_queue,
    file.path(
      OUT,
      "V2_STEP82B_queue_REUSE_ASV.csv"
    ),
    na=""
  )

  write_excel_csv(
    dada2_queue,
    file.path(
      OUT,
      "V2_STEP82B_queue_DADA2.csv"
    ),
    na=""
  )

  write_excel_csv(
    special_queue,
    file.path(
      OUT,
      "V2_STEP82B_queue_SPECIAL_FASTA.csv"
    ),
    na=""
  )

  write_excel_csv(
    repair_queue,
    file.path(
      OUT,
      "V2_STEP82B_queue_RAW_REPAIR.csv"
    ),
    na=""
  )

  # ==========================================================
  # 8. Summary
  # ==========================================================

  summary <- tibble(
    metric=c(
      "Primary_16S_projects",
      "Projects_with_validated_full_ASV_candidate",
      "Projects_with_multiple_nonidentical_full_ASV",
      "Projects_with_partial_ASV_only",
      "Projects_without_validated_ASV",
      "Projects_standard_paired_DADA2_route",
      "Projects_single_FASTQ_validate_DADA2_route",
      "Projects_FASTA_special_route",
      "Projects_raw_partial"
    ),

    value=c(
      nrow(final_decision),

      sum(
        final_decision$asv_decision==
        "FULL_COVERAGE_SINGLE_CONTENT_ASV_REUSABLE_CANDIDATE"
      ),

      sum(
        final_decision$asv_decision==
        "MULTIPLE_NONIDENTICAL_FULL_COVERAGE_ASV_REVIEW"
      ),

      sum(
        final_decision$asv_decision==
        "PARTIAL_ASV_NOT_SUFFICIENT_FOR_FROZEN_MANIFEST"
      ),

      sum(
        final_decision$asv_decision==
        "NO_VALIDATED_ASV"
      ),

      sum(
        final_decision$raw_processing_route==
        "PAIRED_END_FASTQ_DADA2_ROUTE"
      ),

      sum(
        final_decision$raw_processing_route==
        "ONE_FASTQ_PER_RUN_LIKELY_SINGLE_END_VALIDATE_THEN_DADA2"
      ),

      sum(
        final_decision$raw_processing_route==
        "FASTA_NO_QUALITY_SPECIAL_ROUTE_NOT_STANDARD_DADA2"
      ),

      sum(
        final_decision$raw_processing_route==
        "RAW_PARTIAL"
      )
    )
  )

  write_excel_csv(
    summary,
    file.path(
      OUT,
      "V2_STEP82B_summary.csv"
    ),
    na=""
  )

  # ==========================================================
  # README
  # ==========================================================

  readme <- c(
    "SEPSIS V2 - STEP82B ASV VALIDATION",
    paste0("Created: ",Sys.time()),
    "",
    "THIS STEP DOES NOT RUN DADA2 AND DOES NOT ALTER ANY ASV TABLE.",
    "",
    "WHY THIS STEP IS NEEDED",
    "Step82 discovered candidate seqtab/ASV files by filename only. A candidate is not automatically reusable.",
    "Step82B reads RDS seqtabs, detects the sample axis, extracts Run accessions from sample names, and compares them against the Step81D frozen manifest.",
    "",
    "IMPORTANT",
    "- Full frozen-manifest coverage is required before an existing ASV can replace raw reprocessing.",
    "- Multiple non-identical full-coverage seqtabs are NOT automatically resolved.",
    "- Taxonomy can be reassigned later; absence of an old taxonomy file does not invalidate a correct seqtab.",
    "- CRA002354 .fa.gz files contain no FASTQ quality scores and must NOT be sent directly into the standard DADA2 FASTQ pipeline.",
    "- One FASTQ file per Run is treated as likely single-end only; primer/read content must still be checked before DADA2.",
    "",
    "KEY OUTPUT",
    "V2_STEP82B_FINAL_project_processing_decision.csv",
    "",
    "NEXT",
    "After reviewing this output, freeze one processing route per cohort:",
    "1) reuse validated ASV;",
    "2) paired-end DADA2;",
    "3) single-end DADA2 after validation;",
    "4) raw repair;",
    "5) FASTA-specific route."
  )

  writeLines(
    readme,
    file.path(
      OUT,
      "README_STEP82B_ASV_VALIDATION.txt"
    ),
    useBytes=TRUE
  )

  ck("STEP82B COMPLETE")

  cat("\n============================================================\n")
  cat("SEPSIS V2 - STEP82B ASV VALIDATION COMPLETE\n")
  cat("============================================================\n\n")

  print(
    summary,
    n=Inf,
    width=Inf
  )

  cat("\nFINAL PROJECT DECISIONS:\n")
  print(
    final_decision,
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
        "STEP82B FATAL ERROR: ",
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

    ck("STEP82B FAILED")

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
