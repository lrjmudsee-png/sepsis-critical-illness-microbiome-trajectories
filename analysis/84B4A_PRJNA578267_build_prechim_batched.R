# ============================================================
# Sepsis V2 - Step 84B4A
# PRJNA578267 BATCHED PRE-CHIMERA SEQUENCE TABLE BUILDER
#
# Purpose:
#   Diagnose/replace the very slow final aggregation stage after
#   all 207 sample-wise MERGE checkpoints are complete.
#
# IMPORTANT:
#   - Does NOT rerun filtering, error learning, DADA or merging.
#   - Does NOT run chimera removal.
#   - Builds sequence tables in small batches with checkpoints.
#   - Saves the final PRE-CHIMERA sequence table first.
#   - Safe to stop/restart.
#
# R 4.4.0 / Windows
# ============================================================

options(stringsAsFactors = FALSE)

if (!requireNamespace("dada2", quietly=TRUE)) {
  stop("Package 'dada2' is required.")
}
if (!requireNamespace("readr", quietly=TRUE)) {
  stop("Package 'readr' is required.")
}
if (!requireNamespace("dplyr", quietly=TRUE)) {
  stop("Package 'dplyr' is required.")
}
if (!requireNamespace("tibble", quietly=TRUE)) {
  stop("Package 'tibble' is required.")
}

suppressPackageStartupMessages({
  library(dada2)
  library(readr)
  library(dplyr)
  library(tibble)
})

PROJECT <- "PRJNA578267"
ROOT <- "E:/sepsis_project"

PWORK <- file.path(
  ROOT, "data", "_V2_ANALYSIS_READY",
  "01_PROJECTS", PROJECT,
  "05_work", "step84B_full"
)

MERGE_DIR <- file.path(
  PWORK,
  "samplewise_checkpoints",
  "MERGE"
)

SOURCE_ORDER <- file.path(
  ROOT, "results",
  "V2_24B3_PRJNA578267_SAMPLEWISE",
  "PRJNA578267_samplewise_run_order.csv"
)

OUT <- file.path(
  ROOT, "results",
  "V2_24B4A_PRJNA578267_PRECHIM_BATCHED"
)

BATCH_DIR <- file.path(
  PWORK,
  "prechim_batch_checkpoints"
)

MERGED_DIR <- file.path(
  PWORK,
  "prechim_merged_checkpoints"
)

dir.create(OUT, recursive=TRUE, showWarnings=FALSE)
dir.create(BATCH_DIR, recursive=TRUE, showWarnings=FALSE)
dir.create(MERGED_DIR, recursive=TRUE, showWarnings=FALSE)

LOG <- file.path(OUT, "_STEP84B4A_runtime.log")
ERR <- file.path(OUT, "_STEP84B4A_FATAL_ERROR.txt")
if (file.exists(ERR)) unlink(ERR)

logmsg <- function(x) {
  z <- paste0(Sys.time(), " | ", x)
  cat(z, "\n")
  cat(z, "\n", file=LOG, append=TRUE)
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

valid_rds <- function(p) {
  if (!file.exists(p)) return(FALSE)
  tryCatch({
    readRDS(p)
    TRUE
  }, error=function(e) FALSE)
}

atomic_saveRDS <- function(obj, p) {
  tmp <- paste0(p, ".tmp")
  if (file.exists(tmp)) unlink(tmp)
  saveRDS(obj, tmp)
  if (file.exists(p)) unlink(p)
  if (!file.rename(tmp, p)) {
    stop("Could not atomically save: ", p)
  }
}

main <- function() {

  logmsg("STEP84B4A STARTED")

  if (!file.exists(SOURCE_ORDER)) {
    stop("Run-order file is missing: ", SOURCE_ORDER)
  }

  order_df <- safe_csv(SOURCE_ORDER)

  if (nrow(order_df) != 207 || anyDuplicated(order_df$run_id)) {
    stop("Expected exactly 207 unique Runs in samplewise run order.")
  }

  runs <- order_df$run_id

  merge_paths <- file.path(
    MERGE_DIR,
    paste0(runs, ".rds")
  )

  good <- vapply(
    merge_paths,
    valid_rds,
    logical(1)
  )

  if (!all(good)) {
    stop(
      "MERGE checkpoints are incomplete/corrupt. Valid=",
      sum(good),
      "/207."
    )
  }

  logmsg("ALL 207 MERGE CHECKPOINTS VALID")

  # ----------------------------------------------------------
  # 1. Per-sample merger inventory without building seqtab yet
  # ----------------------------------------------------------
  inv <- vector("list", length(runs))

  for (i in seq_along(runs)) {
    m <- readRDS(merge_paths[i])

    inv[[i]] <- tibble(
      sample_index=i,
      run_id=runs[i],
      merged_unique_sequences=if (is.null(m)) 0L else nrow(m),
      merged_reads=if (is.null(m) || !nrow(m)) 0 else sum(m$abundance, na.rm=TRUE)
    )

    rm(m)
    if (i %% 20 == 0) gc(verbose=FALSE)
  }

  inv <- bind_rows(inv)

  write_excel_csv(
    inv,
    file.path(OUT, "PRJNA578267_merge_checkpoint_inventory.csv"),
    na=""
  )

  logmsg(
    paste0(
      "MERGER INVENTORY: summed per-sample unique rows=",
      sum(inv$merged_unique_sequences),
      "; merged reads=",
      sum(inv$merged_reads)
    )
  )

  # ----------------------------------------------------------
  # 2. Batch makeSequenceTable
  # ----------------------------------------------------------
  BATCH_SIZE <- 10L
  batch_id <- ceiling(seq_along(runs) / BATCH_SIZE)
  batches <- split(seq_along(runs), batch_id)

  batch_summary <- vector("list", length(batches))

  for (b in seq_along(batches)) {

    idx <- batches[[b]]

    batch_path <- file.path(
      BATCH_DIR,
      sprintf("batch_%03d_seqtab.rds", b)
    )

    if (valid_rds(batch_path)) {

      tab <- readRDS(batch_path)

      logmsg(
        paste0(
          "BATCH ", b, "/", length(batches),
          " CHECKPOINT EXISTS: ",
          nrow(tab), " samples x ", ncol(tab), " sequences"
        )
      )

    } else {

      logmsg(
        paste0(
          "BATCH ", b, "/", length(batches),
          " START (samples ", min(idx), "-", max(idx), ")"
        )
      )

      mergers <- lapply(
        merge_paths[idx],
        readRDS
      )

      names(mergers) <- runs[idx]

      tab <- makeSequenceTable(mergers)

      if (
        !is.matrix(tab) ||
        nrow(tab) != length(idx) ||
        ncol(tab) == 0
      ) {
        stop(
          "Invalid batch sequence table at batch ",
          b
        )
      }

      atomic_saveRDS(tab, batch_path)

      logmsg(
        paste0(
          "BATCH ", b, "/", length(batches),
          " COMPLETE: ",
          nrow(tab), " samples x ", ncol(tab), " sequences"
        )
      )

      rm(mergers)
    }

    batch_summary[[b]] <- tibble(
      batch=b,
      first_sample=min(idx),
      last_sample=max(idx),
      samples=nrow(tab),
      sequences=ncol(tab),
      reads=sum(tab)
    )

    rm(tab)
    gc(verbose=FALSE)
  }

  batch_summary <- bind_rows(batch_summary)

  write_excel_csv(
    batch_summary,
    file.path(OUT, "PRJNA578267_batch_seqtab_summary.csv"),
    na=""
  )

  logmsg("ALL BATCH SEQUENCE TABLES COMPLETE")

  # ----------------------------------------------------------
  # 3. Incrementally merge batch sequence tables
  # ----------------------------------------------------------
  accumulator <- NULL
  last_done <- 0L

  # Recover the latest valid merged checkpoint.
  for (b in rev(seq_along(batches))) {
    p <- file.path(
      MERGED_DIR,
      sprintf("merged_through_batch_%03d.rds", b)
    )
    if (valid_rds(p)) {
      accumulator <- readRDS(p)
      last_done <- b
      logmsg(
        paste0(
          "RECOVERED MERGED CHECKPOINT THROUGH BATCH ",
          b, ": ",
          nrow(accumulator), " samples x ",
          ncol(accumulator), " sequences"
        )
      )
      break
    }
  }

  start_b <- last_done + 1L

  if (start_b <= length(batches)) {

    for (b in seq.int(start_b, length(batches))) {

      tab <- readRDS(
        file.path(
          BATCH_DIR,
          sprintf("batch_%03d_seqtab.rds", b)
        )
      )

      logmsg(
        paste0(
          "MERGING BATCH TABLE ", b, "/", length(batches),
          " INTO ACCUMULATOR"
        )
      )

      if (is.null(accumulator)) {
        accumulator <- tab
      } else {
        accumulator <- mergeSequenceTables(
          accumulator,
          tab,
          repeats="sum"
        )
      }

      cp <- file.path(
        MERGED_DIR,
        sprintf("merged_through_batch_%03d.rds", b)
      )

      atomic_saveRDS(accumulator, cp)

      logmsg(
        paste0(
          "MERGED CHECKPOINT ", b,
          ": ", nrow(accumulator),
          " samples x ", ncol(accumulator),
          " sequences; reads=", sum(accumulator)
        )
      )

      rm(tab)
      gc(verbose=FALSE)
    }
  }

  if (
    is.null(accumulator) ||
    !is.matrix(accumulator) ||
    nrow(accumulator) != 207 ||
    ncol(accumulator) == 0
  ) {
    stop("Final pre-chimera accumulator is invalid.")
  }

  final_path <- file.path(
    PWORK,
    "PRJNA578267_seqtab_prechim_step84B4A_BATCHED.rds"
  )

  atomic_saveRDS(accumulator, final_path)

  # ----------------------------------------------------------
  # 4. Pre-chimera diagnostics
  # ----------------------------------------------------------
  lens <- nchar(colnames(accumulator))
  abund <- colSums(accumulator)

  diagnostic <- tibble(
    metric=c(
      "samples",
      "prechim_ASVs",
      "total_reads",
      "median_ASV_length",
      "ASV_length_min",
      "ASV_length_max",
      "ASVs_350_500bp_pct",
      "singletons_by_total_abundance",
      "ASVs_abundance_le_2",
      "ASVs_abundance_le_5"
    ),
    value=c(
      nrow(accumulator),
      ncol(accumulator),
      sum(accumulator),
      median(lens),
      min(lens),
      max(lens),
      round(100*mean(lens>=350 & lens<=500),2),
      sum(abund==1),
      sum(abund<=2),
      sum(abund<=5)
    )
  )

  write_excel_csv(
    diagnostic,
    file.path(OUT, "PRJNA578267_PRECHIM_DIAGNOSTIC_SUMMARY.csv"),
    na=""
  )

  asv_len_dist <- tibble(
    asv_length=lens,
    abundance=abund
  ) |>
    count(
      asv_length,
      wt=abundance,
      name="reads"
    ) |>
    arrange(asv_length)

  write_excel_csv(
    asv_len_dist,
    file.path(OUT, "PRJNA578267_PRECHIM_length_distribution.csv"),
    na=""
  )

  writeLines(
    c(
      paste0("Completed: ", Sys.time()),
      paste0("Final pre-chimera table: ", final_path),
      paste0("Samples: ", nrow(accumulator)),
      paste0("Pre-chimera ASVs: ", ncol(accumulator)),
      paste0("Total reads: ", sum(accumulator)),
      "Chimera removal was NOT run in Step84B4A."
    ),
    file.path(OUT, "_STEP84B4A_COMPLETE.ok")
  )

  logmsg("STEP84B4A COMPLETE")

  cat("\n============================================================\n")
  cat("PRJNA578267 PRE-CHIMERA TABLE BUILT\n")
  cat("============================================================\n")
  print(diagnostic, n=Inf, width=Inf)
  cat("\nNo chimera removal was run.\n")
  cat("Upload V2_24B4A_PRJNA578267_PRECHIM_BATCHED for review.\n")
}

tryCatch(
  main(),
  error=function(e) {

    msg <- c(
      paste0("STEP84B4A FATAL ERROR: ", Sys.time()),
      paste0("Message: ", conditionMessage(e)),
      paste0(
        "Call: ",
        paste(deparse(conditionCall(e)), collapse=" ")
      )
    )

    writeLines(msg, ERR, useBytes=TRUE)
    message(paste(msg, collapse="\n"))

    quit(
      save="no",
      status=1,
      runLast=FALSE
    )
  }
)
