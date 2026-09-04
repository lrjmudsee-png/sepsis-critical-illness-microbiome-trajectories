# ============================================================
# Sepsis V2 - Step 82C
# Freeze sequence-processing routes + compare competing ASVs
#
# Search scope is STRICTLY limited to:
#   E:/sepsis_project/data
#   E:/sepsis_project/results
#
# No F drive, Downloads, Desktop, or whole-disk scanning.
#
# R 4.4.0 / Windows
# ============================================================

options(stringsAsFactors = FALSE)

pkgs <- c("readr","dplyr","tidyr","stringr","purrr","tibble")
missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly=TRUE)]
if (length(missing)) install.packages(missing, repos="https://cloud.r-project.org")

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(purrr)
  library(tibble)
})

ROOT <- "E:/sepsis_project"

STEP82B_ROOT <- file.path(ROOT,"results","V2_22B_ASV_VALIDATION")
FREEZE_ROOT <- file.path(ROOT,"results","V2_21D_FINAL_FREEZE_CLINICAL_REPAIRED")

OUT <- file.path(ROOT,"results","V2_22C_PROCESSING_ROUTE_FREEZE")
dir.create(OUT,recursive=TRUE,showWarnings=FALSE)

LOG <- file.path(OUT,"_STEP82C_runtime_checkpoints.txt")
ERR <- file.path(OUT,"_STEP82C_FATAL_ERROR.txt")
if (file.exists(ERR)) unlink(ERR)

ck <- function(x) cat(paste0(x,": ",Sys.time(),"\n"),file=LOG,append=TRUE)

clean_chr <- function(x) {
  x <- trimws(as.character(x))
  x[is.na(x)|x==""|tolower(x)%in%c("na","nan","n/a","null","none")] <- NA_character_
  x
}

safe_csv <- function(p) {
  suppressMessages(read_csv(p,show_col_types=FALSE,progress=FALSE,name_repair="unique"))
}

safe_md5 <- function(p) {
  tryCatch(unname(tools::md5sum(p)),error=function(e) NA_character_)
}

# ------------------------------------------------------------
# Compare two seqtab RDS objects from the same project
# ------------------------------------------------------------
compare_seqtabs <- function(project, p1, p2) {

  a <- readRDS(p1)
  b <- readRDS(p2)

  if (!is.matrix(a) || !is.matrix(b)) {
    return(tibble(
      project=project,
      file_A=p1,
      file_B=p2,
      comparison_status="NON_MATRIX_OBJECT"
    ))
  }

  sampA <- rownames(a)
  sampB <- rownames(b)
  featA <- colnames(a)
  featB <- colnames(b)

  common_samples <- intersect(sampA,sampB)
  common_features <- intersect(featA,featB)

  # sample-level library sizes
  libA <- rowSums(a)
  libB <- rowSums(b)

  lib_cor <- NA_real_
  if (length(common_samples)>=3) {
    lib_cor <- suppressWarnings(
      cor(
        libA[common_samples],
        libB[common_samples],
        method="spearman",
        use="complete.obs"
      )
    )
  }

  richnessA <- rowSums(a>0)
  richnessB <- rowSums(b>0)

  richness_cor <- NA_real_
  if (length(common_samples)>=3) {
    richness_cor <- suppressWarnings(
      cor(
        richnessA[common_samples],
        richnessB[common_samples],
        method="spearman",
        use="complete.obs"
      )
    )
  }

  tibble(
    project=project,
    file_A=p1,
    file_B=p2,
    md5_A=safe_md5(p1),
    md5_B=safe_md5(p2),
    samples_A=nrow(a),
    samples_B=nrow(b),
    features_A=ncol(a),
    features_B=ncol(b),
    common_samples=length(common_samples),
    same_sample_set=setequal(sampA,sampB),
    common_features=length(common_features),
    feature_jaccard=round(
      length(common_features)/
      length(union(featA,featB)),
      6
    ),
    total_reads_A=sum(a),
    total_reads_B=sum(b),
    median_reads_A=median(libA),
    median_reads_B=median(libB),
    sample_library_size_spearman=round(lib_cor,6),
    median_richness_A=median(richnessA),
    median_richness_B=median(richnessB),
    sample_richness_spearman=round(richness_cor,6),
    comparison_status="COMPARED"
  )
}

main <- function() {

  ck("STEP82C STARTED")

  decision_path <- file.path(
    STEP82B_ROOT,
    "V2_STEP82B_FINAL_project_processing_decision.csv"
  )

  coverage_path <- file.path(
    STEP82B_ROOT,
    "V2_STEP82B_seqtab_frozen_sample_coverage_all.csv"
  )

  validation_path <- file.path(
    STEP82B_ROOT,
    "V2_STEP82B_seqtab_validation_all.csv"
  )

  manifest_path <- file.path(
    FREEZE_ROOT,
    "V2_FINAL_FROZEN_16S_sequence_manifest.csv"
  )

  needed <- c(decision_path,coverage_path,validation_path,manifest_path)

  if (!all(file.exists(needed))) {
    stop("One or more Step82B / Step81D inputs are missing.")
  }

  decision <- safe_csv(decision_path)
  coverage <- safe_csv(coverage_path)
  validation <- safe_csv(validation_path)
  manifest <- safe_csv(manifest_path)

  if (nrow(decision)!=16) stop("Expected 16 primary 16S projects.")
  if (nrow(manifest)!=1941) stop("Expected 1941 frozen 16S rows.")

  ck("INPUTS GUARDED")

  # ==========================================================
  # 1. Longitudinal impact of incomplete ASV candidates
  # ==========================================================

  impact_rows <- list()

  for (project in unique(validation$project)) {

    d <- decision |> filter(project==!!project)
    if (!nrow(d)) next

    best <- d$best_candidate_file[1]
    if (is.na(best)) next

    c <- coverage |>
      filter(
        project==!!project,
        candidate_file==best
      )

    if (!nrow(c)) next

    all_patients <- unique(clean_chr(c$patient_id))
    all_patients <- all_patients[!is.na(all_patients)]

    frozen_tp <- c |>
      filter(!is.na(patient_id),!is.na(time_raw)) |>
      distinct(patient_id,time_raw) |>
      count(patient_id,name="n_frozen_tp")

    asv_tp <- c |>
      filter(
        asv_sample_present,
        !is.na(patient_id),
        !is.na(time_raw)
      ) |>
      distinct(patient_id,time_raw) |>
      count(patient_id,name="n_asv_tp")

    pt <- tibble(patient_id=all_patients) |>
      left_join(frozen_tp,by="patient_id") |>
      left_join(asv_tp,by="patient_id") |>
      mutate(
        n_frozen_tp=coalesce(n_frozen_tp,0L),
        n_asv_tp=coalesce(n_asv_tp,0L)
      )

    impact_rows[[length(impact_rows)+1]] <- tibble(
      project=project,
      frozen_rows=nrow(c),
      asv_present_rows=sum(c$asv_sample_present,na.rm=TRUE),
      missing_rows=sum(!c$asv_sample_present,na.rm=TRUE),
      frozen_patients_GE2=sum(pt$n_frozen_tp>=2),
      asv_patients_GE2=sum(pt$n_asv_tp>=2),
      GE2_patients_lost=sum(pt$n_frozen_tp>=2)-sum(pt$n_asv_tp>=2),
      frozen_patients_GE3=sum(pt$n_frozen_tp>=3),
      asv_patients_GE3=sum(pt$n_asv_tp>=3),
      GE3_patients_lost=sum(pt$n_frozen_tp>=3)-sum(pt$n_asv_tp>=3)
    )
  }

  impact <- bind_rows(impact_rows)

  write_excel_csv(
    impact,
    file.path(
      OUT,
      "V2_STEP82C_incomplete_ASV_longitudinal_impact.csv"
    ),
    na=""
  )

  # ==========================================================
  # 2. Compare competing non-identical seqtabs
  # ==========================================================

  comparisons <- list()

  multi_projects <- validation |>
    filter(
      validation_status=="FULL_FROZEN_COVERAGE_CANDIDATE"
    ) |>
    group_by(project) |>
    summarise(
      unique_md5=n_distinct(md5),
      .groups="drop"
    ) |>
    filter(unique_md5>1) |>
    pull(project)

  for (project in multi_projects) {

    v <- validation |>
      filter(
        project==!!project,
        validation_status=="FULL_FROZEN_COVERAGE_CANDIDATE"
      ) |>
      distinct(md5,.keep_all=TRUE)

    if (nrow(v)>=2) {

      pairs <- combn(seq_len(nrow(v)),2,simplify=FALSE)

      for (pair in pairs) {

        p1 <- v$file_path[pair[1]]
        p2 <- v$file_path[pair[2]]

        if (file.exists(p1) && file.exists(p2)) {
          comparisons[[length(comparisons)+1]] <-
            compare_seqtabs(project,p1,p2)
        }
      }
    }
  }

  comparison_df <- if (length(comparisons)) {
    bind_rows(comparisons)
  } else {
    tibble()
  }

  write_excel_csv(
    comparison_df,
    file.path(
      OUT,
      "V2_STEP82C_competing_full_ASV_comparison.csv"
    ),
    na=""
  )

  ck("ASV VERSION COMPARISON COMPLETE")

  # ==========================================================
  # 3. Companion-file audit only inside candidate directories
  # ==========================================================

  companion_rows <- list()

  candidate_dirs <- unique(
    dirname(
      validation$file_path[
        !is.na(validation$file_path)
      ]
    )
  )

  for (d in candidate_dirs) {

    if (!dir.exists(d)) next

    # HARD SAFETY: only project data root is allowed.
    d_norm <- normalizePath(d,winslash="/",mustWork=FALSE)
    allowed <- normalizePath(
      file.path(ROOT,"data"),
      winslash="/",
      mustWork=FALSE
    )

    if (!startsWith(tolower(d_norm),tolower(allowed))) next

    files <- list.files(
      d,
      full.names=TRUE,
      recursive=FALSE
    )

    files <- files[
      str_detect(
        basename(files),
        regex(
          "taxonomy|taxa|track|dada|filter|quality|qc|log|script|\\.R$",
          ignore_case=TRUE
        )
      )
    ]

    if (length(files)) {
      info <- file.info(files)

      companion_rows[[length(companion_rows)+1]] <- tibble(
        directory=d,
        file_path=files,
        file_name=basename(files),
        bytes=as.numeric(info$size),
        modified=as.character(info$mtime)
      )
    }
  }

  companions <- if (length(companion_rows)) {
    bind_rows(companion_rows)
  } else {
    tibble()
  }

  write_excel_csv(
    companions,
    file.path(
      OUT,
      "V2_STEP82C_ASV_companion_file_audit.csv"
    ),
    na=""
  )

  # ==========================================================
  # 4. Freeze processing routes
  # ==========================================================

  route <- decision |>
    left_join(
      impact,
      by="project"
    ) |>
    mutate(
      immediate_priority=case_when(

        project=="PRJEB33360" ~
          "HIGH_CORE_RAW_REPAIR",

        project=="PRJNA595346" ~
          "HIGH_VALIDATION_RAW_REPAIR",

        project=="PRJNA691455" ~
          "HIGH_CORE_ASV_VERSION_SELECTION",

        project=="CRA002354" ~
          "HIGH_CORE_SPECIAL_FASTA",

        analysis_module %in% c(
          "CORE_SEPSIS_LONGITUDINAL",
          "ICU_BACKGROUND_LONGITUDINAL",
          "ICU_INFECTION_LONGITUDINAL_EXTERNAL",
          "NONSEPSIS_LONGITUDINAL_CONTROL"
        ) ~
          "HIGH_MAIN_ANALYSIS",

        analysis_module %in% c(
          "INTERVENTION_LONGITUDINAL_SUPPORT",
          "ORGAN_DYSFUNCTION_LONGITUDINAL_SUPPORT"
        ) ~
          "MEDIUM_SUPPORT",

        TRUE ~
          "LOW_STATIC_SUPPORT"
      ),

      frozen_processing_route=case_when(

        project %in% c(
          "PRJNA430161",
          "PRJNA978257"
        ) ~
          "REUSE_EXISTING_ASV",

        project=="PRJNA691455" ~
          "SELECT_BETWEEN_TWO_FULL_ASV_VERSIONS",

        project=="CRA002354" ~
          "SPECIAL_FASTA_FEATURE_HARMONIZATION",

        project=="PRJEB33360" ~
          "REPAIR_11_RUNS_THEN_REPROCESS_OR_APPEND_COHORT",

        project=="PRJNA595346" ~
          "REPAIR_347_RUNS_THEN_PAIRED_DADA2",

        project=="PRJNA912621" ~
          "SUPPORT_REPAIR_2_RUNS_OR_USE_58_SAMPLE_ASV",

        project=="PRJNA1010969" ~
          "STATIC_SUPPORT_DEFER_1_MISSING_RUN",

        project=="PRJNA797231" ~
          "STATIC_SUPPORT_DEFER_2_MISSING_RUNS",

        raw_processing_route==
        "PAIRED_END_FASTQ_DADA2_ROUTE" ~
          "PAIRED_END_DADA2_AFTER_PARAMETER_AUDIT",

        raw_processing_route==
        "ONE_FASTQ_PER_RUN_LIKELY_SINGLE_END_VALIDATE_THEN_DADA2" ~
          "SINGLE_END_DADA2_AFTER_READ_VALIDATION",

        TRUE ~
          "MANUAL_REVIEW"
      ),

      blocks_primary_longitudinal_analysis=case_when(
        project %in% c(
          "PRJEB33360",
          "PRJNA691455",
          "CRA002354"
        ) ~ TRUE,

        TRUE ~ FALSE
      )
    ) |>
    arrange(
      factor(
        immediate_priority,
        levels=c(
          "HIGH_CORE_RAW_REPAIR",
          "HIGH_CORE_ASV_VERSION_SELECTION",
          "HIGH_CORE_SPECIAL_FASTA",
          "HIGH_VALIDATION_RAW_REPAIR",
          "HIGH_MAIN_ANALYSIS",
          "MEDIUM_SUPPORT",
          "LOW_STATIC_SUPPORT"
        )
      ),
      project
    )

  write_excel_csv(
    route,
    file.path(
      OUT,
      "V2_STEP82C_FROZEN_sequence_processing_routes.csv"
    ),
    na=""
  )

  # ==========================================================
  # 5. Main-analysis queue
  # ==========================================================

  main_queue <- route |>
    filter(
      immediate_priority %in% c(
        "HIGH_CORE_RAW_REPAIR",
        "HIGH_CORE_ASV_VERSION_SELECTION",
        "HIGH_CORE_SPECIAL_FASTA",
        "HIGH_VALIDATION_RAW_REPAIR",
        "HIGH_MAIN_ANALYSIS"
      )
    )

  support_queue <- route |>
    filter(
      immediate_priority %in% c(
        "MEDIUM_SUPPORT",
        "LOW_STATIC_SUPPORT"
      )
    )

  write_excel_csv(
    main_queue,
    file.path(
      OUT,
      "V2_STEP82C_MAIN_ANALYSIS_processing_queue.csv"
    ),
    na=""
  )

  write_excel_csv(
    support_queue,
    file.path(
      OUT,
      "V2_STEP82C_SUPPORT_processing_queue.csv"
    ),
    na=""
  )

  # ==========================================================
  # 6. README
  # ==========================================================

  readme <- c(
    "SEPSIS V2 - STEP82C PROCESSING ROUTE FREEZE",
    paste0("Created: ",Sys.time()),
    "",
    "SEARCH SCOPE",
    "Only E:/sepsis_project/data and previously generated result files are used.",
    "No F drive, Downloads, Desktop, or whole-disk scanning.",
    "",
    "PURPOSE",
    "1. Quantify whether incomplete existing ASVs materially reduce longitudinal eligibility.",
    "2. Compare multiple non-identical full-coverage ASV tables (especially PRJNA691455).",
    "3. Freeze one processing route per cohort.",
    "4. Separate main-analysis blockers from support/static work.",
    "",
    "IMPORTANT",
    "- PRJNA430161 and PRJNA978257 already have single-content full-coverage ASV candidates.",
    "- PRJNA691455 requires explicit selection between two non-identical full-coverage seqtabs.",
    "- PRJEB33360 is a core cohort; its 11 missing Runs should be repaired for the final analysis.",
    "- CRA002354 remains a FASTA-without-quality special route.",
    "- Static support cohorts do not block the main longitudinal pipeline.",
    "",
    "NEXT STEP",
    "Use V2_STEP82C_FROZEN_sequence_processing_routes.csv to start Step83.",
    "Step83 should perform read/primer/quality parameter audits for the cohorts that require new DADA2, not run one generic DADA2 configuration across all cohorts."
  )

  writeLines(
    readme,
    file.path(
      OUT,
      "README_STEP82C_PROCESSING_ROUTE_FREEZE.txt"
    ),
    useBytes=TRUE
  )

  ck("STEP82C COMPLETE")

  cat("\n============================================================\n")
  cat("SEPSIS V2 - STEP82C PROCESSING ROUTE FREEZE COMPLETE\n")
  cat("============================================================\n\n")

  cat("LONGITUDINAL IMPACT OF EXISTING ASV GAPS:\n")
  print(impact,n=Inf,width=Inf)

  cat("\nCOMPETING FULL ASV COMPARISON:\n")
  print(comparison_df,n=Inf,width=Inf)

  cat("\nFROZEN ROUTES:\n")
  print(
    route |>
      select(
        project,
        analysis_module,
        immediate_priority,
        frozen_processing_route,
        blocks_primary_longitudinal_analysis
      ),
    n=Inf,
    width=Inf
  )

  cat("\nOutput:\n",OUT,"\n",sep="")
}

tryCatch(
  main(),
  error=function(e) {

    msg <- c(
      paste0("STEP82C FATAL ERROR: ",Sys.time()),
      paste0("Message: ",conditionMessage(e)),
      paste0("Call: ",paste(deparse(conditionCall(e)),collapse=" "))
    )

    writeLines(msg,ERR,useBytes=TRUE)
    ck("STEP82C FAILED")
    message(paste(msg,collapse="\n"))

    quit(
      save="no",
      status=1,
      runLast=FALSE
    )
  }
)
