# ============================================================
# Sepsis V2 - Step 84A2
# Targeted PRJNA1166732 merge-overlap diagnostic / retry
#
# Why:
#   Step84A failed only for PRJNA1166732 at chimera removal with:
#   "Input must be a valid sequence table."
#
# This script does NOT rerun all 5 cohorts.
# It reuses the existing Step84A filtered pilot FASTQs and learned
# error models for PRJNA1166732, then tests conservative merge
# overlap thresholds: 12, 10, 8, 6 nt.
#
# STRICT INPUT:
#   E:/sepsis_project/data/_V2_ANALYSIS_READY
#   E:/sepsis_project/results/V2_24A_DADA2_PILOT
#
# NO broad filesystem scan.
# NO full-cohort processing.
#
# R 4.4.0 / Windows
# ============================================================

options(stringsAsFactors = FALSE)

pkgs <- c("readr","dplyr","stringr","purrr","tibble")
missing <- pkgs[
  !vapply(pkgs, requireNamespace, logical(1), quietly=TRUE)
]

if (length(missing)) {
  install.packages(missing, repos="https://cloud.r-project.org")
}

if (!requireNamespace("dada2", quietly=TRUE)) {
  if (!requireNamespace("BiocManager", quietly=TRUE)) {
    install.packages("BiocManager", repos="https://cloud.r-project.org")
  }
  BiocManager::install("dada2", ask=FALSE, update=FALSE)
}

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(stringr)
  library(purrr)
  library(tibble)
  library(dada2)
})

ROOT <- "E:/sepsis_project"

CANON <- file.path(
  ROOT,
  "data",
  "_V2_ANALYSIS_READY"
)

PWORK <- file.path(
  CANON,
  "01_PROJECTS",
  "PRJNA1166732",
  "05_work",
  "step84A_pilot"
)

STEP84A <- file.path(
  ROOT,
  "results",
  "V2_24A_DADA2_PILOT"
)

OUT <- file.path(
  ROOT,
  "results",
  "V2_24A2_PRJNA1166732_MERGE_DIAGNOSTIC"
)

dir.create(OUT, recursive=TRUE, showWarnings=FALSE)

LOG <- file.path(OUT, "_STEP84A2_runtime_checkpoints.txt")
ERR <- file.path(OUT, "_STEP84A2_FATAL_ERROR.txt")
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

getN <- function(x) {
  sum(getUniques(x))
}

main <- function() {

  ck("STEP84A2 STARTED")

  selected_path <- file.path(
    STEP84A,
    "PRJNA1166732_pilot_selected_runs.csv"
  )

  errF_path <- file.path(
    PWORK,
    "PRJNA1166732_pilot_errF.rds"
  )

  errR_path <- file.path(
    PWORK,
    "PRJNA1166732_pilot_errR.rds"
  )

  needed <- c(
    selected_path,
    errF_path,
    errR_path
  )

  if (!all(file.exists(needed))) {
    stop(
      "Required PRJNA1166732 Step84A pilot files are missing."
    )
  }

  selected <- safe_csv(selected_path)

  if (nrow(selected)!=5) {
    stop("Expected exactly 5 selected PRJNA1166732 pilot Runs.")
  }

  run_ids <- selected$run_id

  filtFs <- file.path(
    PWORK,
    "filtered",
    paste0(run_ids,"_R1_filt.fastq.gz")
  )

  filtRs <- file.path(
    PWORK,
    "filtered",
    paste0(run_ids,"_R2_filt.fastq.gz")
  )

  if (!all(file.exists(filtFs)) || !all(file.exists(filtRs))) {
    stop(
      "One or more Step84A filtered pilot FASTQ files are missing."
    )
  }

  names(filtFs) <- run_ids
  names(filtRs) <- run_ids

  errF <- readRDS(errF_path)
  errR <- readRDS(errR_path)

  ck("EXISTING FILTERED FASTQ AND ERROR MODELS LOADED")

  # ==========================================================
  # 1. Re-run only derep + dada
  # ==========================================================

  derepFs <- derepFastq(
    filtFs,
    verbose=FALSE
  )

  derepRs <- derepFastq(
    filtRs,
    verbose=FALSE
  )

  names(derepFs) <- run_ids
  names(derepRs) <- run_ids

  dadaFs <- dada(
    derepFs,
    err=errF,
    multithread=FALSE,
    pool=FALSE
  )

  dadaRs <- dada(
    derepRs,
    err=errR,
    multithread=FALSE,
    pool=FALSE
  )

  denoised <- tibble(
    run_id=run_ids,
    denoised_F=as.numeric(
      sapply(dadaFs,getN)[run_ids]
    ),
    denoised_R=as.numeric(
      sapply(dadaRs,getN)[run_ids]
    )
  )

  write_excel_csv(
    denoised,
    file.path(
      OUT,
      "PRJNA1166732_denoised_read_counts.csv"
    ),
    na=""
  )

  ck("DENOISING COMPLETE")

  # ==========================================================
  # 2. Test merge overlap thresholds
  # ==========================================================

  overlap_grid <- c(12L,10L,8L,6L)

  threshold_results <- list()
  sample_results <- list()
  rejected_diagnostics <- list()
  seqtabs <- list()
  nochims <- list()

  for (ov in overlap_grid) {

    ck(paste0("MERGE TEST minOverlap=",ov," STARTED"))

    mergers_all <- mergePairs(
      dadaFs,
      derepFs,
      dadaRs,
      derepRs,
      minOverlap=ov,
      maxMismatch=0,
      returnRejects=TRUE,
      verbose=FALSE
    )

    names(mergers_all) <- run_ids

    # Per-sample accepted/rejected abundance.
    sm <- map_dfr(
      run_ids,
      function(run) {

        m <- mergers_all[[run]]

        if (is.null(m) || !nrow(m)) {
          return(
            tibble(
              run_id=run,
              minOverlap=ov,
              total_candidate_pair_reads=0,
              accepted_reads=0,
              rejected_reads=0,
              accepted_unique_pairs=0,
              rejected_unique_pairs=0
            )
          )
        }

        if (!"accept" %in% names(m)) {
          m$accept <- TRUE
        }

        tibble(
          run_id=run,
          minOverlap=ov,

          total_candidate_pair_reads=sum(
            m$abundance,
            na.rm=TRUE
          ),

          accepted_reads=sum(
            m$abundance[
              m$accept %in% TRUE
            ],
            na.rm=TRUE
          ),

          rejected_reads=sum(
            m$abundance[
              m$accept %in% FALSE
            ],
            na.rm=TRUE
          ),

          accepted_unique_pairs=sum(
            m$accept %in% TRUE,
            na.rm=TRUE
          ),

          rejected_unique_pairs=sum(
            m$accept %in% FALSE,
            na.rm=TRUE
          )
        )
      }
    ) |>
      left_join(
        denoised,
        by="run_id"
      ) |>
      mutate(
        merge_of_forward_denoised_pct=round(
          ifelse(
            denoised_F>0,
            100*accepted_reads/denoised_F,
            NA_real_
          ),
          2
        )
      )

    sample_results[[as.character(ov)]] <- sm

    # Rejected-pair diagnostics.
    rej <- map_dfr(
      run_ids,
      function(run) {

        m <- mergers_all[[run]]

        if (
          is.null(m) ||
          !nrow(m) ||
          !"accept" %in% names(m)
        ) {
          return(tibble())
        }

        z <- m |>
          filter(
            accept %in% FALSE
          )

        if (!nrow(z)) return(tibble())

        keep_cols <- intersect(
          c(
            "abundance",
            "nmatch",
            "nmismatch",
            "nindel"
          ),
          names(z)
        )

        z |>
          select(all_of(keep_cols)) |>
          mutate(
            run_id=run,
            minOverlap=ov,
            .before=1
          )
      }
    )

    if (nrow(rej)) {
      rejected_diagnostics[[as.character(ov)]] <- rej
    }

    # Keep only accepted mergers for sequence-table construction.
    mergers_accept <- lapply(
      mergers_all,
      function(m) {

        if (is.null(m) || !nrow(m)) {
          return(
            data.frame(
              abundance=numeric(),
              sequence=character(),
              forward=integer(),
              reverse=integer(),
              nmatch=integer(),
              nmismatch=integer(),
              nindel=integer(),
              prefer=integer(),
              accept=logical()
            )
          )
        }

        if ("accept" %in% names(m)) {
          m <- m[
            m$accept %in% TRUE,
            ,
            drop=FALSE
          ]
        }

        m
      }
    )

    names(mergers_accept) <- run_ids

    accepted_total <- sum(
      sm$accepted_reads,
      na.rm=TRUE
    )

    zero_merge_samples <- sum(
      sm$accepted_reads==0
    )

    seqtab_valid <- FALSE
    chimera_valid <- FALSE
    n_asv_prechim <- 0L
    n_asv_nochim <- 0L
    nonchim_reads <- 0
    median_asv_length <- NA_real_
    length_ok_pct <- 0

    if (
      accepted_total>0 &&
      zero_merge_samples<5
    ) {

      seqtab_try <- tryCatch(
        makeSequenceTable(
          mergers_accept
        ),
        error=function(e) NULL
      )

      if (
        !is.null(seqtab_try) &&
        is.matrix(seqtab_try) &&
        nrow(seqtab_try)>0 &&
        ncol(seqtab_try)>0
      ) {

        seqtab_valid <- TRUE
        n_asv_prechim <- ncol(seqtab_try)
        seqtabs[[as.character(ov)]] <- seqtab_try

        nochim_try <- tryCatch(
          removeBimeraDenovo(
            seqtab_try,
            method="consensus",
            multithread=FALSE,
            verbose=FALSE
          ),
          error=function(e) NULL
        )

        if (
          !is.null(nochim_try) &&
          is.matrix(nochim_try) &&
          nrow(nochim_try)>0 &&
          ncol(nochim_try)>0
        ) {

          chimera_valid <- TRUE
          n_asv_nochim <- ncol(nochim_try)
          nonchim_reads <- sum(nochim_try)

          lens <- nchar(
            colnames(nochim_try)
          )

          median_asv_length <- median(lens)

          length_ok_pct <- round(
            100*mean(
              lens>=200 &
              lens<=330
            ),
            2
          )

          nochims[[as.character(ov)]] <- nochim_try
        }
      }
    }

    median_merge <- median(
      sm$merge_of_forward_denoised_pct,
      na.rm=TRUE
    )

    pass <- isTRUE(seqtab_valid) &&
      isTRUE(chimera_valid) &&
      zero_merge_samples==0 &&
      is.finite(median_merge) &&
      median_merge>=50 &&
      n_asv_nochim>0 &&
      length_ok_pct>=80

    threshold_results[[as.character(ov)]] <- tibble(
      minOverlap=ov,
      accepted_reads=accepted_total,
      zero_merge_samples=zero_merge_samples,
      median_merge_of_forward_denoised_pct=round(
        median_merge,
        2
      ),
      seqtab_valid=seqtab_valid,
      chimera_step_valid=chimera_valid,
      ASVs_prechim=n_asv_prechim,
      ASVs_nochim=n_asv_nochim,
      nonchim_reads=nonchim_reads,
      median_ASV_length=median_asv_length,
      ASVs_200_330bp_pct=length_ok_pct,
      diagnostic_pass=pass
    )

    ck(paste0("MERGE TEST minOverlap=",ov," COMPLETE"))
  }

  thresholds <- bind_rows(
    threshold_results
  )

  samples <- bind_rows(
    sample_results
  )

  rejects <- if (length(rejected_diagnostics)) {
    bind_rows(
      rejected_diagnostics
    )
  } else {
    tibble()
  }

  write_excel_csv(
    thresholds,
    file.path(
      OUT,
      "PRJNA1166732_overlap_threshold_comparison.csv"
    ),
    na=""
  )

  write_excel_csv(
    samples,
    file.path(
      OUT,
      "PRJNA1166732_overlap_threshold_sample_tracking.csv"
    ),
    na=""
  )

  write_excel_csv(
    rejects,
    file.path(
      OUT,
      "PRJNA1166732_rejected_pair_diagnostics.csv"
    ),
    na=""
  )

  # ==========================================================
  # 3. Choose most conservative threshold that passes
  # ==========================================================

  passed <- thresholds |>
    filter(
      diagnostic_pass
    ) |>
    arrange(
      desc(minOverlap)
    )

  if (nrow(passed)) {

    chosen <- passed$minOverlap[1]
    seqtab_chosen <- seqtabs[[as.character(chosen)]]
    nochim_chosen <- nochims[[as.character(chosen)]]

    saveRDS(
      seqtab_chosen,
      file.path(
        OUT,
        "PRJNA1166732_pilot_seqtab_prechim_RETRY.rds"
      )
    )

    saveRDS(
      nochim_chosen,
      file.path(
        OUT,
        "PRJNA1166732_pilot_seqtab_nochim_RETRY.rds"
      )
    )

    decision <- tibble(
      project="PRJNA1166732",
      retry_status="PASS",
      selected_minOverlap=chosen,
      maxMismatch=0L,
      recommendation=
        "Use the highest tested minOverlap that passes all pilot QC; freeze it for the subsequent full-cohort parameter sheet."
    )

  } else {

    decision <- tibble(
      project="PRJNA1166732",
      retry_status="NO_THRESHOLD_PASSED",
      selected_minOverlap=NA_integer_,
      maxMismatch=0L,
      recommendation=
        "Do not lower overlap further automatically. Audit read orientation / amplicon length / primer handling before any full-cohort processing."
    )
  }

  write_excel_csv(
    decision,
    file.path(
      OUT,
      "PRJNA1166732_FINAL_retry_decision.csv"
    ),
    na=""
  )

  # ==========================================================
  # 4. README
  # ==========================================================

  readme <- c(
    "SEPSIS V2 - STEP84A2 PRJNA1166732 MERGE DIAGNOSTIC",
    paste0("Created: ",Sys.time()),
    "",
    "WHY THIS EXISTS",
    "The original Step84A failed only for PRJNA1166732 at removeBimeraDenovo with 'Input must be a valid sequence table'.",
    "This targeted diagnostic tests whether paired reads were failing to produce a valid sequence table because the default 12-nt overlap requirement is too strict after removal of the 806R/515F primer sequences from 2x150 reads.",
    "",
    "TESTED",
    "minOverlap = 12, 10, 8, 6",
    "maxMismatch = 0 for all tests",
    "",
    "SELECTION RULE",
    "The HIGHEST minOverlap that passes all QC is selected.",
    "The script never lowers overlap below 6 and never uses justConcatenate automatically.",
    "",
    "PASS CRITERIA",
    "- valid non-empty sequence table",
    "- chimera removal succeeds",
    "- no pilot sample has zero merged reads",
    "- median accepted merge / denoised-forward reads >= 50%",
    "- >=80% of non-chimeric ASVs are 200-330 bp",
    "",
    "IMPORTANT",
    "If no tested threshold passes, the project remains HOLD. Do not proceed to full-cohort DADA2 until read orientation / amplicon span / primer handling has been audited."
  )

  writeLines(
    readme,
    file.path(
      OUT,
      "README_STEP84A2_PRJNA1166732.txt"
    ),
    useBytes=TRUE
  )

  ck("STEP84A2 COMPLETE")

  cat("\n============================================================\n")
  cat("SEPSIS V2 - STEP84A2 COMPLETE\n")
  cat("============================================================\n\n")

  print(
    thresholds,
    n=Inf,
    width=Inf
  )

  cat("\nFINAL DECISION:\n")
  print(
    decision,
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
        "STEP84A2 FATAL ERROR: ",
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

    ck("STEP84A2 FAILED")

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
