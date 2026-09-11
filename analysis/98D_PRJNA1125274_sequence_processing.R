# ============================================================
# Sepsis V2 Upgrade - Step 98D
# PRJNA1125274 (SURVEIL) sequence processing (DADA2 + SILVA 138.2)
# ------------------------------------------------------------
# Pipeline: local FASTQ (verified 289-run manifest)
#   -> parameter freeze (pilot mode, one-time)
#   -> filterAndTrim -> learnErrors -> dada -> mergePairs
#   -> removeBimeraDenovo -> SILVA 138.2 assignTaxonomy
#   -> sequencing-depth / duplicate QC
#   -> PRJNA1125274_analysis_object_external_validation.rds
#
# Usage:
#   Rscript 98D_PRJNA1125274_sequence_processing.R pilot   (parameter
#     freeze + small pilot; safe to run before download completes)
#   Rscript 98D_PRJNA1125274_sequence_processing.R full     (after all
#     289 runs verified INCLUDE)
# ============================================================

options(stringsAsFactors = FALSE)
suppressPackageStartupMessages({
  library(dada2); library(ShortRead); library(Biostrings)
  library(readr); library(dplyr); library(tidyr); library(tibble)
  library(stringr); library(tools)
})

args <- commandArgs(trailingOnly = TRUE)
MODE <- if (length(args)) args[1] else "pilot"
MERGE_MISMATCH <- if (length(args) > 1) as.numeric(args[2]) else 0
stopifnot(MODE %in% c("pilot", "full"))

ROOT <- "E:/sepsis_project"
DATA112 <- file.path(ROOT, "data", "PRJNA1125274")
RAW <- file.path(DATA112, "01_correct_289_raw")
FREEZE_DIR <- file.path(DATA112, "03_dada2")
dir.create(FREEZE_DIR, recursive = TRUE, showWarnings = FALSE)
PARAM_FREEZE <- file.path(
  FREEZE_DIR, "PRJNA1125274_DADA2_parameter_freeze.csv")
PILOT_DIR <- file.path(FREEZE_DIR, "pilot")
dir.create(PILOT_DIR, recursive = TRUE, showWarnings = FALSE)
FILT_DIR <- file.path(FREEZE_DIR, "filtered")
dir.create(FILT_DIR, recursive = TRUE, showWarnings = FALSE)

silva_train <- "E:/sepsis_project/data/_V2_ANALYSIS_READY/00_REFERENCE/SILVA_138.2_DADA2_SSU/silva_nr99_v138.2_toGenus_trainset.fa.gz"
silva_species <- "E:/sepsis_project/data/_V2_ANALYSIS_READY/00_REFERENCE/SILVA_138.2_DADA2_SSU/silva_v138.2_assignSpecies.fa.gz"
stopifnot(file.exists(silva_train), file.exists(silva_species))

# manifest (authoritative, data-folder copy from Step98C)
manifest_path <- file.path(DATA112, "PRJNA1125274_FINAL_289_run_manifest.csv")
stopifnot(file.exists(manifest_path))
man <- suppressMessages(read_csv(manifest_path, show_col_types = FALSE))

# ---------- parameter freeze (pilot) ----------
freeze_params <- function(probe_r1_files) {
  set.seed(20260907)
  sr1 <- readFastq(probe_r1_files[1], withIds = TRUE, n = 40000)
  len1 <- as.integer(width(sr1))
  q1 <- as(quality(sr1), "matrix")
  res <- data.frame(
    parameter = c("platform", "read_length", "chemistry", "insert_region",
                  "primer_state", "trimLeft_R1", "trimLeft_R2",
                  "truncLen_R1", "truncLen_R2", "maxEE_R1", "maxEE_R2",
                  "truncQ", "rm.phix", "minOverlap", "maxMismatch",
                  "min_read_length_observed_R1",
                  "max_read_length_observed_R1",
                  "meanQ_R1", "q30_fraction_R1",
                  "filter_nthreads", "seed", "taxonomy_reference",
                  "assigned_region", "min_depth_keep"),
    value = c("Illumina MiSeq", "250 (2x250 V2 chemistry)", "V2 500 cycles",
              "16S V3-V6 (per SURVEIL paper; kit: Arrow Microbiota Solution B)",
              "PRIMERS_ALREADY_TRIMMED_BY_PROVIDER_5P_CLEAN_NO_PRIMER_IN_READ",
              0, 0, 0, 0, 6, 6, 2, TRUE,
              "12", as.character(MERGE_MISMATCH),  # minOverlap/maxMismatch
              as.character(min(len1)), as.character(max(len1)),
              sprintf("%.1f", mean(q1, na.rm = TRUE)),
              sprintf("%.3f", mean(q1 >= 30, na.rm = TRUE)),
              4, 20260907, basename(silva_train), "V3-V6",
              2000))
  write_csv(res, PARAM_FREEZE)
  cat("Parameter freeze written:", PARAM_FREEZE, "\n")
}

pilot_runs <- function() {
  # choose complete runs spanning letters/hospitals that exist locally
  keep <- man |>
    filter(!is.na(fastq_r1_path), !is.na(fastq_r2_path),
           file.exists(fastq_r1_path), file.exists(fastq_r2_path))
  set.seed(20260907)
  # balance over hospital x letter
  chosen <- keep |>
    group_by(hospital, alias_letter) |>
    slice_sample(n = 1) |>
    ungroup() |>
    slice_head(n = 12)
  chosen <- chosen |> distinct(run_accession, .keep_all = TRUE)
  cat("Pilot runs:", nrow(chosen), "\n")

  fqR1 <- chosen$fastq_r1_path
  fqR2 <- chosen$fastq_r2_path
  # quality profile for freeze (first R1)
  freeze_params(fqR1)

  outR1 <- file.path(PILOT_DIR, paste0(basename(fqR1)))
  outR2 <- file.path(PILOT_DIR, paste0(basename(fqR2)))
  outR1 <- gsub("\\.fastq\\.gz$", "_filt.fastq.gz", outR1)
  outR2 <- gsub("\\.fastq\\.gz$", "_filt.fastq.gz", outR2)

  if (!all(file.exists(outR1)) || !all(file.exists(outR2))) {
    filt <- dada2::filterAndTrim(
      fwd = fqR1, rev = fqR2, filt = outR1, filt.rev = outR2,
      truncLen = c(0, 0), trimLeft = c(0, 0),
      maxEE = c(6, 6), truncQ = 2, rm.phix = TRUE,
      compress = TRUE, multithread = 4, verbose = TRUE)
  } else {
    cat("Reusing existing filtered pilot files\n")
    filt <- NULL
  }

  errF <- dada2::learnErrors(outR1, multithread = 4, randomize = TRUE,
                             nbases = 1e8)
  errR <- dada2::learnErrors(outR2, multithread = 4, randomize = TRUE,
                             nbases = 1e8)
  ddF <- dada2::dada(outR1, err = errF, multithread = 4)
  ddR <- dada2::dada(outR2, err = errR, multithread = 4)
  mergers <- dada2::mergePairs(ddF, outR1, ddR, outR2, verbose = TRUE,
                              maxMismatch = MERGE_MISMATCH)
  seqtab <- dada2::makeSequenceTable(mergers)
  seqtab_nochim <- dada2::removeBimeraDenovo(
    seqtab, method = "consensus", multithread = 4, verbose = TRUE)

  cat("Pilot retention table (per sample):\n")
  print(filt)
  cat("Unique merged variants:", ncol(seqtab), "\n")
  cat("After chimera removal:", ncol(seqtab_nochim), "\n")
  cat("Merged read length distribution (nt):\n")
  print(summary(nchar(getSequences(seqtab_nochim))))

  pdf(file.path(PILOT_DIR, "PRJNA1125274_pilot_quality_profiles.pdf"))
  try(dada2::plotQualityProfile(fqR1[1:6]), silent = TRUE)
  dev.off()
  saveRDS(list(filt = filt, seqtab_nochim = seqtab_nochim),
          file.path(PILOT_DIR, "PRJNA1125274_pilot_results.rds"))
  cat("PILOT COMPLETE\n")
}

full_run <- function() {
  stopifnot(file.exists(PARAM_FREEZE))
  inc <- man |> filter(inclusion == "INCLUDE")
  missing_fq <- inc |> filter(!file.exists(fastq_r1_path) |
                                !file.exists(fastq_r2_path))
  if (nrow(missing_fq)) {
    stop("Runs marked INCLUDE but FASTQ missing: ",
         paste(missing_fq$run_accession, collapse = ", "))
  }
  cat("Full run over", nrow(inc), "runs\n")

  fqR1 <- inc$fastq_r1_path
  fqR2 <- inc$fastq_r2_path
  outR1 <- file.path(FILT_DIR, gsub("\\.fastq\\.gz$", "_filt.fastq.gz",
                                    basename(fqR1)))
  outR2 <- file.path(FILT_DIR, gsub("\\.fastq\\.gz$", "_filt.fastq.gz",
                                    basename(fqR2)))

  # filter only samples not already filtered
  todo <- which(!file.exists(outR1) | !file.exists(outR2))
  if (length(todo)) {
    filt <- dada2::filterAndTrim(
      fwd = fqR1[todo], rev = fqR2[todo],
      filt = outR1[todo], filt.rev = outR2[todo],
      truncLen = c(0, 0), trimLeft = c(0, 0),
      maxEE = c(6, 6), truncQ = 2, rm.phix = TRUE,
      compress = TRUE, multithread = 4, verbose = FALSE)
    cat("filterAndTrim done for", length(todo), "samples\n")
    write_csv(as.data.frame(filt), file.path(
      FREEZE_DIR, "PRJNA1125274_filter_stats.csv"))
  }

  errF <- dada2::learnErrors(outR1, multithread = 4, randomize = TRUE)
  errR <- dada2::learnErrors(outR2, multithread = 4, randomize = TRUE)
  ddF <- dada2::dada(outR1, err = errF, multithread = 4)
  ddR <- dada2::dada(outR2, err = errR, multithread = 4)
  mergers <- dada2::mergePairs(ddF, outR1, ddR, outR2, verbose = TRUE,
                              maxMismatch = MERGE_MISMATCH)
  seqtab <- dada2::makeSequenceTable(mergers)
  seqtab_nochim <- dada2::removeBimeraDenovo(
    seqtab, method = "consensus", multithread = 4, verbose = TRUE)
  rownames(seqtab_nochim) <- inc$run_accession

  saveRDS(seqtab_nochim, file.path(
    FREEZE_DIR, "PRJNA1125274_seqtab_final.rds"))
  cat("Seqtab dim:", dim(seqtab_nochim), "\n")

  taxa <- dada2::assignTaxonomy(
    seqtab_nochim, refFasta = silva_train, multithread = 4,
    tryRC = TRUE, verbose = TRUE)
  taxa_sp <- dada2::addSpecies(taxa, refFasta = silva_species, verbose = TRUE)
  saveRDS(taxa_sp, file.path(
    FREEZE_DIR, "PRJNA1125274_taxa_final.rds"))
  cat("Taxonomy assigned\n")

  # ---- build analysis object ----
  st <- seqtab_nochim
  sample_names <- rownames(st)
  meta <- inc |>
    mutate(run_id = run_accession, sample_id = run_accession) |>
    select(-fastq_r1_path, -fastq_r2_path)
  # align meta to samples
  meta <- meta[match(sample_names, meta$run_accession), ]
  stopifnot(all(sample_names == meta$run_accession))

  counts <- st
  colnames(counts) <- paste0("ASV_", seq_len(ncol(st)))
  asv_tax <- as.data.frame(taxa_sp, stringsAsFactors = FALSE) |>
    mutate(ASV_ID = colnames(counts),
           ASV_sequence = colnames(st))

  obj <- list(
    counts = counts,
    taxonomy = as_tibble(asv_tax),
    metadata = as.data.frame(meta),
    provenance = list(
      project = "PRJNA1125274",
      created = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
      step = "98D",
      sequence_source = "local correct-289 ENA fastq (01_correct_289_raw)",
      taxonomy_source = "SILVA 138.2 (dada2 assignTaxonomy/addSpecies)",
      metadata_source = "ENA manifest + BioSample XML + SURVEIL paper",
      minimum_target_reads = 2000,
      processing_parameters = PARAM_FREEZE,
      cross_cohort_ASV_merge_allowed = FALSE
    ),
    same_day_replicate_sensitivity = NULL
  )

  out_obj <- file.path(
    FREEZE_DIR, "PRJNA1125274_analysis_object_external_validation.rds")
  saveRDS(obj, out_obj)
  write_csv(meta, file.path(
    FREEZE_DIR, "PRJNA1125274_external_validation_metadata.csv"))
  cat("Analysis object written:", out_obj, "\n")
  cat("FULL RUN COMPLETE\n")
}

if (MODE == "pilot") {
  pilot_runs()
} else {
  full_run()
}
