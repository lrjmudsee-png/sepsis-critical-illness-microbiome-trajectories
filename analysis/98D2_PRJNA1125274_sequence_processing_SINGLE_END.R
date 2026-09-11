# ============================================================
# Sepsis V2 Upgrade - Step 98D (single-end R1 strategy)
# PRJNA1125274 (SURVEIL) DADA2 processing
# ------------------------------------------------------------
# RATIONALE (frozen): the SURVEIL 16S libraries were sequenced
# 2x250 (MiSeq V2 500 cycles) on a V3-V6 amplicon kit. Paired-end
# overlap assembly is not achievable for the majority of reads
# (merge fraction <5% in pilots), because the insert exceeds the
# combined read coverage. We therefore freeze a SINGLE-END R1
# analysis for this cohort (reads begin immediately after the
# forward primer; 5' primer-free, quality Q30 median). All
# downstream ecological metrics are within-cohort so region
# choice does not change the estimand or geometry rules.
#
# Usage: Rscript 98D2_....R pilot | full
#   pilot -> 8-sample SE pilot (can run before downloads finish)
#   full  -> all INCLUDE runs of the final 289 manifest
# Output: analysis object with Step87B-compatible structure.
# ============================================================

options(stringsAsFactors = FALSE)
suppressPackageStartupMessages({
  library(dada2); library(ShortRead); library(Biostrings)
  library(readr); library(dplyr); library(tidyr); library(tibble)
  library(stringr); library(tools)
})

args <- commandArgs(trailingOnly = TRUE)
MODE <- if (length(args)) args[1] else "pilot"
stopifnot(MODE %in% c("pilot", "full"))

ROOT <- "E:/sepsis_project"
DATA112 <- file.path(ROOT, "data", "PRJNA1125274")
RAW <- file.path(DATA112, "01_correct_289_raw")
OUT_D <- file.path(DATA112, "03_dada2")
dir.create(OUT_D, recursive = TRUE, showWarnings = FALSE)
SE_DIR <- file.path(OUT_D, "se")
dir.create(SE_DIR, recursive = TRUE, showWarnings = FALSE)
FILT_SE <- file.path(OUT_D, "se_filtered")
dir.create(FILT_SE, recursive = TRUE, showWarnings = FALSE)

silva_train <- "E:/sepsis_project/data/_V2_ANALYSIS_READY/00_REFERENCE/SILVA_138.2_DADA2_SSU/silva_nr99_v138.2_toGenus_trainset.fa.gz"
silva_species <- "E:/sepsis_project/data/_V2_ANALYSIS_READY/00_REFERENCE/SILVA_138.2_DADA2_SSU/silva_v138.2_assignSpecies.fa.gz"
stopifnot(file.exists(silva_train), file.exists(silva_species))

PARAM_FILE <- file.path(OUT_D, "PRJNA1125274_SE_parameter_freeze.csv")
if (!file.exists(PARAM_FILE)) {
  params <- data.frame(
    parameter = c("strategy", "rationale", "platform", "read_length",
                  "chemistry", "insert_region", "primer_state",
                  "trimLeft", "truncLen", "maxEE", "truncQ", "rm.phix",
                  "error_learning", "denoise", "chimera", "taxonomy",
                  "min_reads_keep", "seed"),
    value = c(
      "SINGLE_END_R1",
      "Paired-end overlap unachievable on V3-V6 2x250 (merge <5%); R1 is primer-free 16S (5' starts at amplicon start)",
      "Illumina MiSeq", "250", "V2 500 cycles",
      "16S V3-V6 (SURVEIL/Arrow Microbiota Solution B)",
      "PRIMERS_ALREADY_TRIMMED_BY_PROVIDER_5P_CLEAN",
      0, 0, 6, 2, TRUE,
      "dada2::learnErrors (R1)", "dada2::dada (R1)",
      "removeBimeraDenovo consensus", "SILVA 138.2 trainset + assignSpecies",
      2000, 20260907
    ), stringsAsFactors = FALSE)
  write_csv(params, PARAM_FILE)
  cat("SE parameter freeze written:", PARAM_FILE, "\n")
}

manifest_path <- file.path(DATA112, "PRJNA1125274_FINAL_289_run_manifest.csv")
stopifnot(file.exists(manifest_path))
man <- suppressMessages(read_csv(manifest_path, show_col_types = FALSE))

run_pipeline <- function(runs_df, tag) {
  r1s <- runs_df$fastq_r1_path
  stopifnot(all(file.exists(r1s)))
  outF <- file.path(FILT_SE, gsub("\\.fastq\\.gz$", "_filt.fastq.gz",
                                  basename(r1s)))

  todo <- which(!file.exists(outF))
  if (length(todo)) {
    cat("Filtering", length(todo), "R1 samples\n")
    filt <- dada2::filterAndTrim(r1s[todo], outF[todo], truncLen = 0,
                                 trimLeft = 0, maxEE = 6, truncQ = 2,
                                 rm.phix = TRUE, compress = TRUE,
                                 multithread = 1, verbose = TRUE)
  } else {
    cat("Reusing existing filtered SE files\n")
  }

  errF <- dada2::learnErrors(outF, multithread = 1, randomize = TRUE)
  ddF <- dada2::dada(outF, err = errF, multithread = 1)
  st <- dada2::makeSequenceTable(ddF)
  st_nochim <- dada2::removeBimeraDenovo(st, method = "consensus",
                                         multithread = 1, verbose = TRUE)
  # align row names to run accession in same order as runs_df
  run_ids <- runs_df$run_accession
  if (nrow(st_nochim) == length(run_ids)) {
    rownames(st_nochim) <- run_ids
  } else {
    nm <- sub("\\.fastq\\.gz$", "", basename(outF))
    newids <- runs_df$run_accession[match(rownames(st_nochim), nm)]
    if (any(is.na(newids))) {
      stop("Cannot map seqtab rows to run accessions")
    }
    rownames(st_nochim) <- newids
  }

  saveRDS(st_nochim, file.path(OUT_D, paste0("PRJNA1125274_", tag,
                                             "_seqtab_se.rds")))
  cat("Seqtab dim:", dim(st_nochim), "\n")
  cat("Read length distribution (nt):\n")
  print(summary(nchar(getSequences(st_nochim))))

  if (MODE == "full") {
    taxa <- dada2::assignTaxonomy(st_nochim, refFasta = silva_train,
                                  multithread = 1, tryRC = TRUE, verbose = TRUE)
    taxa_sp <- dada2::addSpecies(taxa, refFasta = silva_species, verbose = TRUE)
    saveRDS(taxa_sp, file.path(OUT_D, "PRJNA1125274_taxa_final_se.rds"))

    meta <- runs_df[match(rownames(st_nochim), runs_df$run_accession), , drop = FALSE]
    counts <- st_nochim
    colnames(counts) <- paste0("ASV_", seq_len(ncol(counts)))
    asv_tax <- as.data.frame(taxa_sp, stringsAsFactors = FALSE) |>
      mutate(ASV_ID = colnames(counts), ASV_sequence = colnames(st_nochim))
    obj <- list(
      counts = counts,
      taxonomy = as_tibble(asv_tax),
      metadata = as.data.frame(meta),
      provenance = list(
        project = "PRJNA1125274",
        created = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
        step = "98D",
        strategy = "SINGLE_END_R1 (see SE parameter freeze)",
        sequence_source = "local correct-289 ENA fastq (01_correct_289_raw)",
        taxonomy_source = "SILVA 138.2 (dada2 assignTaxonomy/addSpecies)",
        metadata_source = "ENA manifest + BioSample XML + SURVEIL paper",
        minimum_target_reads = 2000,
        processing_parameters = PARAM_FILE,
        cross_cohort_ASV_merge_allowed = FALSE
      ),
      same_day_replicate_sensitivity = NULL
    )
    out_obj <- file.path(OUT_D,
                         "PRJNA1125274_analysis_object_external_validation.rds")
    saveRDS(obj, out_obj)
    write_csv(meta, file.path(OUT_D,
                              "PRJNA1125274_external_validation_metadata.csv"))
    cat("Analysis object written:", out_obj, "\n")
  }
  cat(tag, "COMPLETE\n")
  st_nochim
}

if (MODE == "pilot") {
  set.seed(20260907)
  keep <- man |>
    filter(file.exists(fastq_r1_path)) |>
    group_by(hospital, alias_letter) |>
    slice_sample(n = 1) |>
    ungroup() |>
    distinct(run_accession, .keep_all = TRUE) |>
    slice_head(n = 8)
  cat("SE pilot runs:", nrow(keep), "\n")
  invisible(run_pipeline(keep, "pilot"))
} else {
  inc <- man |> filter(inclusion == "INCLUDE", file.exists(fastq_r1_path))
  missing_fq <- inc |> filter(!file.exists(fastq_r1_path))
  if (nrow(missing_fq)) stop("INCLUDE runs missing R1: ",
                             paste(missing_fq$run_accession, collapse = ","))
  cat("Full SE run over", nrow(inc), "runs\n")
  invisible(run_pipeline(inc, "full"))
}
