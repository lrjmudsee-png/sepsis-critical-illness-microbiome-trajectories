# ============================================================
# Sepsis V2 - Step 83C
# Platform-aware sequence route + DADA2 PILOT parameter freeze
#
# STRICT INPUT:
#   E:/sepsis_project/data/_V2_ANALYSIS_READY
#   E:/sepsis_project/results/V2_23B_READ_QUALITY_PRIMER_AUDIT
#
# NO broad filesystem scan.
# NO full DADA2 run in this step.
#
# This step converts Step83B evidence into:
#   1) one platform-aware route per cohort
#   2) a pilot DADA2 parameter sheet for cohorts that are ready
#   3) explicit HOLD routes for PacBio / Ion Torrent anomalies
#   4) repair / ASV-reuse / static-defer queues
#
# R 4.4.0 / Windows
# ============================================================

options(stringsAsFactors = FALSE)

pkgs <- c("readr","dplyr","stringr","tibble")
missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly=TRUE)]
if (length(missing)) {
  install.packages(missing, repos="https://cloud.r-project.org")
}

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(stringr)
  library(tibble)
})

ROOT <- "E:/sepsis_project"

CANON <- file.path(
  ROOT,
  "data",
  "_V2_ANALYSIS_READY"
)

FREEZE <- file.path(
  CANON,
  "00_FREEZE"
)

STEP83B <- file.path(
  ROOT,
  "results",
  "V2_23B_READ_QUALITY_PRIMER_AUDIT"
)

OUT <- file.path(
  ROOT,
  "results",
  "V2_23C_PLATFORM_PARAMETER_FREEZE"
)

dir.create(OUT, recursive=TRUE, showWarnings=FALSE)

LOG <- file.path(OUT, "_STEP83C_runtime_checkpoints.txt")
ERR <- file.path(OUT, "_STEP83C_FATAL_ERROR.txt")
if (file.exists(ERR)) unlink(ERR)

ck <- function(x) {
  cat(paste0(x, ": ", Sys.time(), "\n"), file=LOG, append=TRUE)
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

main <- function() {

  ck("STEP83C STARTED")

  audit_path <- file.path(
    STEP83B,
    "V2_STEP83B_project_parameter_audit.csv"
  )

  mate_path <- file.path(
    STEP83B,
    "V2_STEP83B_project_mate_summary.csv"
  )

  if (!file.exists(audit_path) || !file.exists(mate_path)) {
    stop("Step83B audit inputs are missing.")
  }

  audit <- safe_csv(audit_path)
  mate <- safe_csv(mate_path)

  if (nrow(audit) != 16) {
    stop("Expected 16 primary 16S projects in Step83B audit.")
  }

  ck("STEP83B INPUTS GUARDED")

  # ==========================================================
  # 1. Platform-aware route registry
  #
  # Notes on evidence encoded here:
  # - PRJEB68229: PacBio Sequel IIe, full-length V1-V9,
  #   27F/1492R; paper reports primers removed and orientation
  #   normalized before public FASTQ. DADA2 long-read support
  #   requires PacBioErrfun/BAND_SIZE=32.
  #
  # - PRJEB67798: Ion PGM V5/V6, expected ~350 bp; original
  #   paper used QIIME2 cutadapt trim-single + dada2 denoise-single.
  #   Our local FASTQ median length ~70 bp is inconsistent with
  #   expected amplicon size, so this cohort is HOLD pending
  #   source/processed-feature validation.
  #
  # - PRJEB33360: MiSeq V3-V4; local files are single merged-like
  #   ~478 bp reads with 341F present and 11 frozen Runs missing.
  #
  # - PRJNA851469: Illumina MiSeq V4, 2x250.
  # - PRJNA578267: Illumina MiSeq V3-V4, 341F/806R, 2x300.
  # - PRJNA1166732: V4 515F/806R, 2x150; Step83B shows R1 begins
  #   with 806R and R2 begins with 515F, so primer orientation
  #   must follow observed read orientation rather than labels.
  # ==========================================================

  params <- tribble(
    ~project, ~platform_route, ~target_region,
    ~primer_R1, ~primer_R2, ~primer_handling,
    ~dada2_mode, ~truncLen_R1, ~truncLen_R2,
    ~maxEE_R1, ~maxEE_R2, ~minLen, ~maxLen,
    ~error_function, ~band_size,
    ~parameter_status, ~next_action,

    "PRJEB33360",
    "ILLUMINA_MISEQ_PREMERGED_SINGLE",
    "V3-V4",
    "341F","reverse_V3V4_end",
    "HOLD_UNTIL_11_RUN_REPAIR_THEN_REMOVE_TERMINAL_PRIMERS",
    "SINGLE_END_AFTER_REPAIR",
    NA,NA,2,NA,350,550,
    "DEFAULT_DADA2",NA,
    "HOLD_RAW_REPAIR",
    "Repair 11 missing Runs; verify all 249 files have same merged-read structure, then run one consistent single-end denoising pipeline.",

    "PRJNA691455",
    "ASV_REUSE",
    "V3-V4",
    NA,NA,"NOT_APPLICABLE",
    "NO_NEW_DADA2",
    NA,NA,NA,NA,NA,NA,
    NA,NA,
    "FROZEN_ASV_SELECTED",
    "Use Step83A selected rerun_v2 ASV; do not reprocess raw reads.",

    "CRA002354",
    "FASTA_NO_QUALITY_SPECIAL",
    "amplicon_16S",
    NA,NA,"NO_FASTQ_QUALITY_AVAILABLE",
    "NO_STANDARD_DADA2",
    NA,NA,NA,NA,NA,NA,
    NA,NA,
    "SPECIAL_ROUTE_HOLD",
    "Do not run standard FASTQ DADA2. Define FASTA-specific feature/taxonomy harmonization route.",

    "PRJNA595346",
    "ILLUMINA_PAIRED_RAW_PARTIAL",
    "16S_amplicon",
    NA,NA,"REASSESS_AFTER_RAW_REPAIR",
    "PAIRED_END_AFTER_REPAIR",
    240,220,2,4,NA,NA,
    "DEFAULT_DADA2",NA,
    "HOLD_RAW_REPAIR",
    "Repair 347 missing Runs first; rerun read/primer audit on complete set before freezing final parameters.",

    "PRJEB67798",
    "ION_TORRENT_PGM_SINGLE",
    "V5-V6_expected_~350bp",
    "published_forward_V5V6","published_reverse_V5V6",
    "ORIGINAL_STUDY_USED_CUTADAPT_TRIM_SINGLE",
    "ION_TORRENT_SINGLE_END_HOLD",
    NA,NA,2,NA,200,450,
    "DEFAULT_DADA2",NA,
    "HOLD_LENGTH_ANOMALY",
    "Local reads have median length ~70 bp despite published ~350-bp V5/V6 amplicon. Validate source/processed feature table before denoising.",

    "PRJEB68229",
    "PACBIO_HIFI_FULL_LENGTH",
    "V1-V9_full_length",
    "27F","1492R",
    "PUBLISHED_PRIMERS_ALREADY_REMOVED_AND_READS_ORIENTED",
    "PACBIO_DADA2_HOLD_QUALITY_ENCODING",
    NA,NA,2,NA,1000,1600,
    "PacBioErrfun",32,
    "HOLD_QUALITY_ENCODING",
    "Reads are full-length PacBio HiFi. Step83B Phred+33-derived mean Q is inconsistent with paper-reported Q29-Q40; resolve quality encoding before PacBio DADA2.",

    "PRJEB82425",
    "ILLUMINA_MISEQ_PAIRED",
    "V3-V4",
    "341F","785R",
    "CUTADAPT_REMOVE_R1_341F_R2_785R",
    "PAIRED_END_DADA2_PILOT",
    270,220,2,4,NA,NA,
    "DEFAULT_DADA2",NA,
    "READY_FOR_PILOT",
    "Run 3-5 sample pilot; verify post-filter retention and merged-read success before full cohort.",

    "PRJNA516701",
    "ILLUMINA_PAIRED",
    "16S_short_amplicon",
    NA,NA,
    "NO_COMMON_PRIMER_DETECTED_TREAT_AS_ALREADY_TRIMMED_FOR_PILOT",
    "PAIRED_END_DADA2_PILOT",
    240,220,2,4,NA,NA,
    "DEFAULT_DADA2",NA,
    "READY_FOR_PILOT",
    "Pilot without trimLeft/cutadapt; verify merge success and inferred ASV length distribution.",

    "PRJNA578267",
    "ILLUMINA_MISEQ_PAIRED",
    "V3-V4",
    "341F","806R",
    "CUTADAPT_REMOVE_R1_341F_R2_806R",
    "PAIRED_END_DADA2_PILOT",
    270,220,2,4,NA,NA,
    "DEFAULT_DADA2",NA,
    "READY_FOR_PILOT",
    "R2 quality deteriorates after ~220-240 cycles; pilot truncLen 270/220 and verify overlap/merge.",

    "PRJNA851469",
    "ILLUMINA_MISEQ_PAIRED",
    "V4",
    NA,NA,
    "NO_COMMON_PRIMER_DETECTED_TREAT_AS_ALREADY_TRIMMED_FOR_PILOT",
    "PAIRED_END_DADA2_PILOT",
    240,220,2,4,NA,NA,
    "DEFAULT_DADA2",NA,
    "READY_FOR_PILOT",
    "Published V4 2x250 MiSeq; pilot 240/220, then verify merge success and V4 ASV lengths.",

    "PRJNA1166732",
    "ILLUMINA_NOVASEQ_PAIRED",
    "V4",
    "806R","515F",
    "CUTADAPT_OBSERVED_ORIENTATION_R1_806R_R2_515F",
    "PAIRED_END_DADA2_PILOT",
    0,0,2,2,100,NA,
    "DEFAULT_DADA2",NA,
    "READY_FOR_PILOT",
    "Reads are 2x150 and high quality. Keep full post-primer reads; observed primer orientation is reversed relative to conventional F/R naming.",

    "PRJNA430161",
    "ASV_REUSE",
    "amplicon_16S",
    NA,NA,"NOT_APPLICABLE",
    "NO_NEW_DADA2",
    NA,NA,NA,NA,NA,NA,
    NA,NA,
    "FROZEN_ASV_SELECTED",
    "Use validated full-coverage existing ASV.",

    "PRJNA912621",
    "SINGLE_END_SUPPORT_PARTIAL",
    "V3-V4_like",
    "341F",NA,
    "SUPPORT_COHORT_DECISION_PENDING_2_RUN_REPAIR",
    "SINGLE_END_AFTER_REPAIR_OR_USE_PARTIAL_ASV",
    NA,NA,2,NA,180,NA,
    "DEFAULT_DADA2",NA,
    "SUPPORT_NOT_BLOCKING",
    "Only 2 Runs missing; 58-sample ASV preserves >=2-timepoint eligibility but loses two >=3 patients. Repair later if needed.",

    "PRJNA1010969",
    "STATIC_SUPPORT_DEFER",
    "V4",
    "515F",NA,
    "DEFER",
    "NO_CURRENT_REPROCESSING",
    NA,NA,NA,NA,NA,NA,
    NA,NA,
    "STATIC_DEFER",
    "Static support; 1 missing Run does not block V2 longitudinal analysis.",

    "PRJNA797231",
    "STATIC_SUPPORT_DEFER",
    "V3-V4",
    "341F",NA,
    "DEFER",
    "NO_CURRENT_REPROCESSING",
    NA,NA,NA,NA,NA,NA,
    NA,NA,
    "STATIC_DEFER",
    "Static support; 2 missing Runs do not block V2 longitudinal analysis.",

    "PRJNA978257",
    "ASV_REUSE",
    "V3-V4",
    NA,NA,"NOT_APPLICABLE",
    "NO_NEW_DADA2",
    NA,NA,NA,NA,NA,NA,
    NA,NA,
    "FROZEN_ASV_SELECTED",
    "Use validated full-coverage existing ASV."
  )

  # ==========================================================
  # 2. Attach Step83B evidence
  # ==========================================================

  evidence <- audit |>
    select(
      project,
      analysis_module,
      expected_runs,
      runs_missing_local_raw,
      canonical_raw_files,
      R1_files,
      R2_files,
      single_files,
      fasta_files,
      primer_state,
      parameter_audit_status
    )

  final <- params |>
    left_join(evidence, by="project") |>
    arrange(
      factor(
        parameter_status,
        levels=c(
          "READY_FOR_PILOT",
          "FROZEN_ASV_SELECTED",
          "HOLD_RAW_REPAIR",
          "HOLD_LENGTH_ANOMALY",
          "HOLD_QUALITY_ENCODING",
          "SPECIAL_ROUTE_HOLD",
          "SUPPORT_NOT_BLOCKING",
          "STATIC_DEFER"
        )
      ),
      project
    )

  if (nrow(final)!=16 || any(is.na(final$analysis_module))) {
    stop("Platform parameter registry did not map cleanly to all 16 cohorts.")
  }

  write_excel_csv(
    final,
    file.path(
      OUT,
      "V2_STEP83C_PLATFORM_AWARE_PARAMETER_FREEZE.csv"
    ),
    na=""
  )

  # ==========================================================
  # 3. Pilot queue
  # ==========================================================

  pilot <- final |>
    filter(parameter_status=="READY_FOR_PILOT") |>
    mutate(
      pilot_samples_per_project=5L,
      pilot_success_min_filter_retention_pct=60,
      pilot_success_min_merge_retention_pct=50,
      pilot_require_no_systematic_zero_sample=TRUE,
      pilot_require_asv_length_plausible=TRUE,
      full_run_allowed_only_after_pilot_pass=TRUE
    )

  hold <- final |>
    filter(
      str_detect(parameter_status,"^HOLD") |
      parameter_status=="SPECIAL_ROUTE_HOLD"
    )

  reuse <- final |>
    filter(parameter_status=="FROZEN_ASV_SELECTED")

  support <- final |>
    filter(
      parameter_status %in% c(
        "SUPPORT_NOT_BLOCKING",
        "STATIC_DEFER"
      )
    )

  write_excel_csv(
    pilot,
    file.path(
      OUT,
      "V2_STEP83C_queue_DADA2_PILOT.csv"
    ),
    na=""
  )

  write_excel_csv(
    hold,
    file.path(
      OUT,
      "V2_STEP83C_queue_HOLD_AND_REPAIR.csv"
    ),
    na=""
  )

  write_excel_csv(
    reuse,
    file.path(
      OUT,
      "V2_STEP83C_queue_ASV_REUSE.csv"
    ),
    na=""
  )

  write_excel_csv(
    support,
    file.path(
      OUT,
      "V2_STEP83C_queue_SUPPORT_DEFER.csv"
    ),
    na=""
  )

  # ==========================================================
  # 4. Source-method evidence registry
  # ==========================================================

  source_registry <- tribble(
    ~project, ~method_evidence, ~source_identifier,

    "PRJEB68229",
    "Full-length 16S V1-V9 with 27F/1492R; PacBio Sequel IIe HiFi; primers removed and reverse reads re-oriented before deposited analysis; original study used UNOISE3.",
    "Early reduction in gut microbiota diversity in critically ill patients is associated with mortality",

    "PRJEB67798",
    "16S V5/V6 expected ~350 bp; Ion PGM 400 chemistry; original processing used QIIME2 cutadapt demux/trim-single and dada2 denoise-single.",
    "Intestinal dysbiosis as an intraoperative predictor of septic complications",

    "PRJEB33360",
    "Illumina MiSeq V3-V4; published primers include 341F/785R.",
    "Gut microbiota profiles in critically ill patients, potential biomarkers and risk variables for sepsis",

    "PRJNA578267",
    "Illumina MiSeq 2x300; V3-V4; modified 341F/806R.",
    "Marked Changes in Gut Microbiota in Cardio-Surgical Intensive Care Patients",

    "PRJNA851469",
    "Illumina MiSeq 2x250; 16S V4; original study processed paired FASTQ with DADA2.",
    "Dysbiosis of a microbiota-immune metasystem in critical illness is associated with nosocomial infections",

    "PRJNA1166732",
    "16S V4 515F/806R; paired-end 2x150 NovaSeq 6000.",
    "Effects of fecal microbiota transplantation and probiotics on the gut microbiome in antibiotic-treated septic patients"
  )

  write_excel_csv(
    source_registry,
    file.path(
      OUT,
      "V2_STEP83C_published_method_evidence_registry.csv"
    ),
    na=""
  )

  # ==========================================================
  # 5. Freeze into canonical workspace
  # ==========================================================

  dir.create(FREEZE, recursive=TRUE, showWarnings=FALSE)

  file.copy(
    file.path(
      OUT,
      "V2_STEP83C_PLATFORM_AWARE_PARAMETER_FREEZE.csv"
    ),
    file.path(
      FREEZE,
      "V2_STEP83C_PLATFORM_AWARE_PARAMETER_FREEZE.csv"
    ),
    overwrite=TRUE
  )

  file.copy(
    file.path(
      OUT,
      "V2_STEP83C_queue_DADA2_PILOT.csv"
    ),
    file.path(
      FREEZE,
      "V2_STEP83C_queue_DADA2_PILOT.csv"
    ),
    overwrite=TRUE
  )

  # ==========================================================
  # 6. Summary
  # ==========================================================

  summary <- tibble(
    metric=c(
      "Primary_16S_projects",
      "Ready_for_DADA2_pilot",
      "Existing_ASV_frozen",
      "Hold_raw_repair",
      "Hold_platform_or_quality_issue",
      "FASTA_special",
      "Support_not_blocking",
      "Static_deferred"
    ),
    value=c(
      nrow(final),
      sum(final$parameter_status=="READY_FOR_PILOT"),
      sum(final$parameter_status=="FROZEN_ASV_SELECTED"),
      sum(final$parameter_status=="HOLD_RAW_REPAIR"),
      sum(final$parameter_status %in% c(
        "HOLD_LENGTH_ANOMALY",
        "HOLD_QUALITY_ENCODING"
      )),
      sum(final$parameter_status=="SPECIAL_ROUTE_HOLD"),
      sum(final$parameter_status=="SUPPORT_NOT_BLOCKING"),
      sum(final$parameter_status=="STATIC_DEFER")
    )
  )

  write_excel_csv(
    summary,
    file.path(
      OUT,
      "V2_STEP83C_summary.csv"
    ),
    na=""
  )

  readme <- c(
    "SEPSIS V2 - STEP83C PLATFORM-AWARE PARAMETER FREEZE",
    paste0("Created: ",Sys.time()),
    "",
    "IMPORTANT",
    "Step83B demonstrated that the project is not composed only of standard Illumina paired-end data.",
    "Therefore no single generic DADA2 pipeline will be applied across all cohorts.",
    "",
    "READY FOR PILOT DADA2",
    paste(pilot$project, collapse=", "),
    "",
    "PLATFORM-SPECIFIC HOLDS",
    "PRJEB68229 = PacBio HiFi full-length; resolve quality-score interpretation before PacBioErrfun DADA2.",
    "PRJEB67798 = Ion Torrent V5/V6; local read-length distribution is inconsistent with the published ~350-bp amplicon and needs source validation.",
    "CRA002354 = FASTA without base-quality scores; not standard FASTQ DADA2.",
    "",
    "WHY PILOT BEFORE FULL RUN",
    "truncLen values are evidence-based starting parameters from sampled quality profiles.",
    "A 5-sample pilot must pass filtering, merging, and ASV-length QC before the whole cohort is processed.",
    "",
    "NEXT STEP",
    "Step84A should run pilot DADA2 only on the READY_FOR_PILOT projects and write tracking/merge/ASV-length QC.",
    "Do not process HOLD cohorts until their route-specific issue is resolved."
  )

  writeLines(
    readme,
    file.path(
      OUT,
      "README_STEP83C_PLATFORM_PARAMETER_FREEZE.txt"
    ),
    useBytes=TRUE
  )

  ck("STEP83C COMPLETE")

  cat("\n============================================================\n")
  cat("SEPSIS V2 - STEP83C PLATFORM PARAMETER FREEZE COMPLETE\n")
  cat("============================================================\n\n")

  print(summary,n=Inf,width=Inf)

  cat("\nREADY FOR PILOT:\n")
  print(
    pilot |>
      select(
        project,
        target_region,
        primer_handling,
        truncLen_R1,
        truncLen_R2,
        maxEE_R1,
        maxEE_R2
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
      paste0("STEP83C FATAL ERROR: ",Sys.time()),
      paste0("Message: ",conditionMessage(e)),
      paste0("Call: ",paste(deparse(conditionCall(e)),collapse=" "))
    )

    writeLines(msg,ERR,useBytes=TRUE)
    ck("STEP83C FAILED")
    message(paste(msg,collapse="\n"))

    quit(
      save="no",
      status=1,
      runLast=FALSE
    )
  }
)
