# ============================================================
# Sepsis V2 - Step94B1
# CRA002354 FASTA QC + 97% OTU + SILVA 138.2 GENUS PROFILING
#
# Scientific rationale
# --------------------
# CRA002354 is deposited as per-sample FASTA (.fa.gz), not FASTQ.
# Therefore base-quality-aware DADA2 denoising cannot be applied
# legitimately to these deposited files.
#
# The original study used:
#   - V3/V4 16S
#   - merged/quality-filtered reads
#   - 300-500 bp length window
#   - chimera removal
#   - 97% OTU clustering
#
# This rescue pipeline therefore:
# 1) audits all 131 FASTA files;
# 2) requires strong evidence that they are merged/filtered amplicons;
# 3) pools reads while embedding sample labels;
# 4) VSEARCH dereplication + de novo chimera removal + 97% clustering;
# 5) maps reads back to OTUs to obtain a sample x OTU table;
# 6) assigns taxonomy to OTU centroids with SILVA 138.2 via dada2;
# 7) aggregates to genus-level counts/relative abundances;
# 8) joins frozen clinical/source metadata;
# 9) writes >=2000 and >=4000 read-depth analysis-ready sets.
#
# IMPORTANT:
# - These are OTUs, NOT ASVs.
# - Do not merge CRA002354 OTUs with ASVs from other cohorts.
# - Cross-cohort synthesis must remain effect/direction based.
# ============================================================

options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(stringr)
  library(tibble)
})

ROOT <- "E:/sepsis_project"

RAW_DIR <- file.path(
  ROOT,
  "data",
  "_V2_ANALYSIS_READY",
  "01_PROJECTS",
  "CRA002354",
  "01_raw"
)

MASTER <- file.path(
  ROOT,
  "data",
  "_V2_ANALYSIS_READY",
  "00_FREEZE",
  "V2_FINAL_FROZEN_master_metadata_CLINICAL_REPAIRED_20260814_131620.csv"
)

DATA_OUT <- file.path(
  ROOT,
  "data",
  "CRA002354",
  "03_vsearch97_silva1382"
)

OUT <- file.path(
  ROOT,
  "results",
  "V2_34B1_CRA002354_FASTA_OTU97_SILVA1382"
)

dir.create(DATA_OUT, recursive = TRUE, showWarnings = FALSE)
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

logfile <- file.path(OUT, "_STEP94B1_runtime.txt")
writeLines(paste("START", Sys.time()), logfile)

logmsg <- function(...) {
  z <- paste0(...)
  cat(z, "\n")
  cat(z, "\n", file = logfile, append = TRUE)
}

norm_chr <- function(x) {
  z <- str_squish(as.character(x))
  z[z %in% c("", "NA", "N/A", "NULL", "None", "none", ".", "-")] <- NA_character_
  z
}

run_id_from_path <- function(f) {
  str_extract(toupper(basename(f)), "CRR\\d+")
}

# ------------------------------------------------------------
# 1. Locate VSEARCH
# ------------------------------------------------------------

find_vsearch <- function() {

  candidates <- c(
    Sys.which("vsearch"),
    Sys.which("vsearch.exe"),
    file.path(
      ROOT,
      "tools",
      "vsearch-2.31.0-win-x86_64",
      "bin",
      "vsearch.exe"
    ),
    file.path(
      ROOT,
      "tools",
      "vsearch",
      "vsearch.exe"
    )
  )

  candidates <- unique(candidates[nzchar(candidates)])
  candidates <- candidates[file.exists(candidates)]

  if (length(candidates) == 0) return(NA_character_)
  normalizePath(candidates[1], winslash = "/", mustWork = TRUE)
}

VSEARCH <- find_vsearch()

if (is.na(VSEARCH)) {
  stop(
    paste0(
      "VSEARCH not found. Run GET_VSEARCH_2.31.0.ps1 first, ",
      "or place vsearch.exe on PATH."
    )
  )
}

vsearch_version <- tryCatch(
  system2(VSEARCH, "--version", stdout = TRUE, stderr = TRUE),
  error = function(e) character()
)

writeLines(
  vsearch_version,
  file.path(OUT, "00_VSEARCH_VERSION.txt")
)

logmsg("VSEARCH: ", VSEARCH)

# ------------------------------------------------------------
# 2. Locate SILVA 138.2 training set
# ------------------------------------------------------------

find_silva <- function() {

  search_roots <- c(
    file.path(ROOT, "reference"),
    file.path(ROOT, "references"),
    file.path(ROOT, "ref"),
    file.path(ROOT, "data"),
    ROOT
  )

  search_roots <- unique(search_roots[dir.exists(search_roots)])

  hits <- character()

  for (rr in search_roots) {
    h <- list.files(
      rr,
      recursive = TRUE,
      full.names = TRUE,
      pattern = "silva.*train.*set.*\\.fa(\\.gz)?$",
      ignore.case = TRUE
    )
    hits <- c(hits, h)

    # Avoid repeatedly scanning the entire project once hits exist in
    # a preferred reference directory.
    if (
      length(hits) > 0 &&
      basename(rr) %in% c("reference", "references", "ref")
    ) break
  }

  hits <- unique(hits)

  if (length(hits) == 0) return(character())

  score <- rep(0, length(hits))
  z <- tolower(hits)

  score <- score +
    100 * str_detect(z, "138[._-]?2|1382|v138[._-]?2")
  score <- score +
    50 * str_detect(z, "nr99")
  score <- score +
    20 * str_detect(z, "train")
  score <- score -
    50 * str_detect(z, "138[._-]?1|v138[._-]?1")

  hits[order(score, decreasing = TRUE)]
}

silva_hits <- find_silva()

write_csv(
  tibble(
    candidate = silva_hits,
    preferred = seq_along(silva_hits) == 1
  ),
  file.path(OUT, "01_SILVA_REFERENCE_CANDIDATES.csv")
)

if (length(silva_hits) == 0) {
  stop(
    "No SILVA training-set FASTA found under E:/sepsis_project."
  )
}

SILVA_TRAIN <- silva_hits[1]

if (!str_detect(
  tolower(SILVA_TRAIN),
  "138[._-]?2|1382|v138[._-]?2"
)) {
  stop(
    paste0(
      "A taxonomy training set was found, but it is not clearly SILVA 138.2: ",
      SILVA_TRAIN,
      ". Inspect 01_SILVA_REFERENCE_CANDIDATES.csv."
    )
  )
}

logmsg("SILVA train set: ", SILVA_TRAIN)

# ------------------------------------------------------------
# 3. Identify exactly 131 canonical FASTA files
# ------------------------------------------------------------

if (!dir.exists(RAW_DIR)) {
  stop(paste0("Missing canonical CRA002354 raw directory: ", RAW_DIR))
}

fa_files <- list.files(
  RAW_DIR,
  full.names = TRUE,
  pattern = "^CRR\\d+\\.fa\\.gz$",
  ignore.case = TRUE
)

fa_files <- fa_files[
  !is.na(run_id_from_path(fa_files))
]

fa_files <- fa_files[
  !duplicated(run_id_from_path(fa_files))
]

fa_files <- fa_files[
  order(run_id_from_path(fa_files))
]

if (length(fa_files) != 131) {
  stop(
    paste0(
      "Expected 131 unique CRR*.fa.gz files in canonical raw directory, found ",
      length(fa_files),
      "."
    )
  )
}

# ------------------------------------------------------------
# 4. FASTA streaming QC + pooled sample-labelled FASTA
# ------------------------------------------------------------

POOLED <- file.path(DATA_OUT, "CRA002354_pooled_sample_labeled.fa")

if (file.exists(POOLED)) file.remove(POOLED)

out_con <- file(POOLED, open = "wt")
on.exit(try(close(out_con), silent = TRUE), add = TRUE)

qc_list <- vector("list", length(fa_files))

global_seq_counter <- 0L

for (i in seq_along(fa_files)) {

  f <- fa_files[i]
  run <- run_id_from_path(f)

  con <- gzfile(f, open = "rt")

  lengths <- integer()
  ambiguous <- logical()

  current_len <- 0L
  current_ambig <- FALSE
  nseq <- 0L

  repeat {

    lines <- readLines(
      con,
      n = 100000,
      warn = FALSE
    )

    if (length(lines) == 0) break

    for (ln in lines) {

      if (startsWith(ln, ">")) {

        if (nseq > 0L) {
          lengths <- c(lengths, current_len)
          ambiguous <- c(ambiguous, current_ambig)
        }

        nseq <- nseq + 1L
        global_seq_counter <- global_seq_counter + 1L

        current_len <- 0L
        current_ambig <- FALSE

        writeLines(
          paste0(
            ">SEQ",
            global_seq_counter,
            ";sample=",
            run,
            ";"
          ),
          out_con
        )

      } else {

        s <- toupper(
          gsub("\\s+", "", ln)
        )

        if (nzchar(s)) {
          current_len <- current_len + nchar(s)

          if (grepl("[^ACGT]", s)) {
            current_ambig <- TRUE
          }

          writeLines(s, out_con)
        }
      }
    }
  }

  close(con)

  if (nseq > 0L) {
    lengths <- c(lengths, current_len)
    ambiguous <- c(ambiguous, current_ambig)
  }

  if (length(lengths) != nseq) {
    stop(
      paste0(
        "FASTA parser mismatch for ",
        basename(f),
        ": headers=",
        nseq,
        " lengths=",
        length(lengths)
      )
    )
  }

  qc_list[[i]] <- tibble(
    run_id = run,
    file = f,
    compressed_bytes = file.info(f)$size,
    reads = nseq,
    min_length = ifelse(nseq > 0, min(lengths), NA_integer_),
    q1_length = ifelse(nseq > 0, unname(quantile(lengths, 0.25)), NA_real_),
    median_length = ifelse(nseq > 0, median(lengths), NA_real_),
    q3_length = ifelse(nseq > 0, unname(quantile(lengths, 0.75)), NA_real_),
    max_length = ifelse(nseq > 0, max(lengths), NA_integer_),
    fraction_300_500 =
      ifelse(nseq > 0, mean(lengths >= 300 & lengths <= 500), NA_real_),
    ambiguous_reads =
      ifelse(nseq > 0, sum(ambiguous), NA_integer_),
    ambiguous_fraction =
      ifelse(nseq > 0, mean(ambiguous), NA_real_)
  )

  if (i %% 10 == 0 || i == length(fa_files)) {
    logmsg(
      "FASTA QC: ",
      i,
      "/",
      length(fa_files),
      " samples"
    )
  }
}

close(out_con)

qc <- bind_rows(qc_list)

write_csv(
  qc,
  file.path(OUT, "02_FASTA_SAMPLE_QC.csv")
)

global_qc <- qc %>%
  summarise(
    samples = n(),
    total_reads = sum(reads),
    min_reads = min(reads),
    median_reads = median(reads),
    max_reads = max(reads),
    weighted_fraction_300_500 =
      weighted.mean(
        fraction_300_500,
        reads
      ),
    total_ambiguous_reads =
      sum(ambiguous_reads),
    weighted_ambiguous_fraction =
      sum(ambiguous_reads) /
      sum(reads)
  )

write_csv(
  global_qc,
  file.path(OUT, "03_FASTA_GLOBAL_QC.csv")
)

# Gate the rescue pipeline.
length_gate <- (
  global_qc$weighted_fraction_300_500 >= 0.95
)

ambiguity_gate <- (
  global_qc$weighted_ambiguous_fraction <= 0.01
)

sample_depth_gate <- (
  sum(qc$reads >= 2000) >= 120
)

qc_gate <- tibble(
  criterion = c(
    ">=95% reads length 300-500 bp",
    "<=1% reads contain ambiguous bases",
    ">=120/131 samples have >=2000 deposited reads"
  ),
  passed = c(
    length_gate,
    ambiguity_gate,
    sample_depth_gate
  ),
  observed = c(
    global_qc$weighted_fraction_300_500,
    global_qc$weighted_ambiguous_fraction,
    sum(qc$reads >= 2000)
  )
)

write_csv(
  qc_gate,
  file.path(OUT, "04_FASTA_PROCESSING_GATE.csv")
)

if (!all(qc_gate$passed)) {

  writeLines(
    c(
      "CRA002354 FASTA QC gate failed.",
      "Do not continue to OTU clustering until 02-04 outputs are reviewed."
    ),
    file.path(OUT, "_STEP94B1_FASTA_QC_REVIEW_REQUIRED.txt")
  )

  stop("CRA002354 FASTA QC gate failed.")
}

# ------------------------------------------------------------
# 5. VSEARCH pipeline
# ------------------------------------------------------------

UNIQUES <- file.path(DATA_OUT, "CRA002354_uniques.fa")
NONCHIM <- file.path(DATA_OUT, "CRA002354_uniques_nonchimeric.fa")
OTUS <- file.path(DATA_OUT, "CRA002354_OTU97_centroids.fa")
OTU_TABLE <- file.path(DATA_OUT, "CRA002354_OTU97_table.tsv")

THREADS <- max(
  1L,
  min(
    12L,
    parallel::detectCores(logical = TRUE) - 1L
  )
)

run_vsearch <- function(args, label) {

  logmsg("VSEARCH START: ", label)

  stdout_file <- file.path(
    OUT,
    paste0(
      "VSEARCH_",
      gsub("[^A-Za-z0-9]+", "_", label),
      ".stdout.txt"
    )
  )

  stderr_file <- file.path(
    OUT,
    paste0(
      "VSEARCH_",
      gsub("[^A-Za-z0-9]+", "_", label),
      ".stderr.txt"
    )
  )

  status <- system2(
    VSEARCH,
    args = args,
    stdout = stdout_file,
    stderr = stderr_file
  )

  if (!identical(status, 0L)) {
    stop(
      paste0(
        "VSEARCH failed at ",
        label,
        " with exit code ",
        status,
        ". See ",
        stderr_file
      )
    )
  }

  logmsg("VSEARCH COMPLETE: ", label)
}

# A. exact-sequence dereplication with abundance sizes
run_vsearch(
  c(
    "--fastx_uniques", shQuote(POOLED),
    "--fastaout", shQuote(UNIQUES),
    "--sizeout",
    "--relabel", "Uniq",
    "--threads", as.character(THREADS)
  ),
  "01_dereplication"
)

# B. de novo chimera removal
run_vsearch(
  c(
    "--uchime3_denovo", shQuote(UNIQUES),
    "--nonchimeras", shQuote(NONCHIM),
    "--sizein",
    "--sizeout"
  ),
  "02_chimera_removal"
)

# C. abundance-guided 97% OTU clustering
run_vsearch(
  c(
    "--cluster_size", shQuote(NONCHIM),
    "--id", "0.97",
    "--centroids", shQuote(OTUS),
    "--sizein",
    "--sizeout",
    "--relabel", "OTU_",
    "--threads", as.character(THREADS)
  ),
  "03_cluster_97pct"
)

# D. map original labelled reads back to nonchimeric OTU centroids
run_vsearch(
  c(
    "--usearch_global", shQuote(POOLED),
    "--db", shQuote(OTUS),
    "--id", "0.97",
    "--strand", "plus",
    "--otutabout", shQuote(OTU_TABLE),
    "--threads", as.character(THREADS),
    "--maxaccepts", "1",
    "--maxrejects", "64"
  ),
  "04_map_reads_to_OTUs"
)

# ------------------------------------------------------------
# 6. Read OTU table and depth audit
# ------------------------------------------------------------

otu <- read_tsv(
  OTU_TABLE,
  show_col_types = FALSE,
  name_repair = "minimal"
)

if (nrow(otu) == 0 || ncol(otu) < 2) {
  stop("VSEARCH OTU table is empty.")
}

names(otu)[1] <- "otu_id"

otu_ids <- as.character(otu$otu_id)

counts <- as.matrix(
  otu[, -1, drop = FALSE]
)

storage.mode(counts) <- "numeric"

rownames(counts) <- otu_ids

# VSEARCH OTU table is OTU x sample.
depth <- colSums(counts)

depth_df <- tibble(
  run_id = names(depth),
  assigned_OTU97_reads = as.numeric(depth),
  keep_ge2000 = assigned_OTU97_reads >= 2000,
  keep_ge4000 = assigned_OTU97_reads >= 4000
)

write_csv(
  depth_df,
  file.path(OUT, "05_OTU97_ASSIGNED_DEPTH.csv")
)

# ------------------------------------------------------------
# 7. SILVA 138.2 taxonomy on OTU centroids
# ------------------------------------------------------------

if (!requireNamespace("dada2", quietly = TRUE)) {
  stop("R package 'dada2' is required for SILVA taxonomy assignment.")
}

if (!requireNamespace("Biostrings", quietly = TRUE)) {
  stop("Bioconductor package 'Biostrings' is required.")
}

dna <- Biostrings::readDNAStringSet(OTUS)

otu_seq_ids <- sub(
  ";.*$",
  "",
  names(dna)
)

seqs <- as.character(dna)
names(seqs) <- otu_seq_ids

logmsg(
  "Assigning SILVA taxonomy to ",
  length(seqs),
  " OTU centroids"
)

tax <- dada2::assignTaxonomy(
  seqs,
  refFasta = SILVA_TRAIN,
  multithread = TRUE,
  tryRC = TRUE,
  minBoot = 50
)

tax_df <- as.data.frame(
  tax,
  stringsAsFactors = FALSE
) %>%
  rownames_to_column("otu_id")

write_csv(
  tax_df,
  file.path(DATA_OUT, "CRA002354_OTU97_SILVA1382_taxonomy.csv")
)

# ------------------------------------------------------------
# 8. Build stable genus labels
# ------------------------------------------------------------

make_genus_label <- function(row) {

  ranks <- c(
    "Kingdom",
    "Phylum",
    "Class",
    "Order",
    "Family",
    "Genus"
  )

  ranks <- ranks[ranks %in% names(row)]

  vals <- as.character(row[ranks])
  names(vals) <- ranks

  genus <- vals["Genus"]

  if (
    length(genus) == 1 &&
    !is.na(genus) &&
    nzchar(genus)
  ) {
    return(genus)
  }

  # Preserve unclassified abundance rather than silently dropping it.
  nonmissing <- vals[
    !is.na(vals) &
      nzchar(vals)
  ]

  if (length(nonmissing) == 0) {
    return("Unclassified_Bacteria")
  }

  last_rank <- names(nonmissing)[length(nonmissing)]
  last_val <- nonmissing[length(nonmissing)]

  paste0(
    "Unclassified_",
    last_rank,
    "_",
    gsub("[^A-Za-z0-9]+", "_", last_val)
  )
}

tax_df$genus_label <- apply(
  tax_df,
  1,
  make_genus_label
)

# Align taxonomy to OTU count table.
tax_map <- tax_df %>%
  select(
    otu_id,
    genus_label
  )

genus_for_otu <- tax_map$genus_label[
  match(
    rownames(counts),
    tax_map$otu_id
  )
]

genus_for_otu[
  is.na(genus_for_otu) |
    genus_for_otu == ""
] <- "Unclassified_Bacteria"

# Aggregate OTU counts to genus.
genus_counts <- rowsum(
  counts,
  group = genus_for_otu,
  reorder = FALSE
)

# genus x sample -> sample x genus
genus_counts_t <- t(genus_counts)

genus_rel <- genus_counts_t /
  rowSums(genus_counts_t)

genus_count_df <- as.data.frame(
  genus_counts_t,
  check.names = FALSE
) %>%
  rownames_to_column("Run_ID")

genus_rel_df <- as.data.frame(
  genus_rel,
  check.names = FALSE
) %>%
  rownames_to_column("Run_ID")

write_csv(
  genus_count_df,
  file.path(
    DATA_OUT,
    "CRA002354_genus_counts_OTU97_SILVA1382.csv"
  )
)

write_csv(
  genus_rel_df,
  file.path(
    DATA_OUT,
    "CRA002354_genus_relative_abundance_OTU97_SILVA1382.csv"
  )
)

# ------------------------------------------------------------
# 9. Join frozen clinical + infection-source metadata
# ------------------------------------------------------------

master <- read_csv(
  MASTER,
  show_col_types = FALSE,
  guess_max = 50000,
  name_repair = "unique"
)

meta <- master %>%
  filter(
    toupper(project) == "CRA002354"
  ) %>%
  transmute(
    Run_ID = norm_chr(run_id),
    patient_id = norm_chr(patient_id),
    sample_id = norm_chr(sample_id),
    time_day = suppressWarnings(as.numeric(time_day)),
    time_raw = norm_chr(time_raw),
    baseline_sepsis_status =
      norm_chr(baseline_sepsis_status),
    sample_sepsis_status =
      norm_chr(sample_sepsis_status),
    infection_source_raw =
      norm_chr(infection_source_standard),
    outcome_28d = norm_chr(outcome_28d),
    baseline_sofa =
      suppressWarnings(as.numeric(baseline_sofa)),
    baseline_apache_ii =
      suppressWarnings(as.numeric(baseline_apache_ii)),
    baseline_lactate =
      suppressWarnings(as.numeric(baseline_lactate))
  ) %>%
  mutate(
    infection_source_group = case_when(
      tolower(infection_source_raw) == "lung" ~ "RESPIRATORY",
      tolower(infection_source_raw) %in%
        c("abdominal", "intestinal") ~ "ABDOMINAL_GI",
      tolower(infection_source_raw) == "blood" ~ "BLOODSTREAM",
      tolower(infection_source_raw) == "urinary" ~ "URINARY",
      tolower(infection_source_raw) ==
        "surgical_site" ~ "SKIN_SOFT_TISSUE",
      tolower(infection_source_raw) == "other" ~ "OTHER_UNKNOWN",
      TRUE ~ "UNMAPPED"
    ),
    pulmonary_binary = ifelse(
      infection_source_group == "RESPIRATORY",
      "PULMONARY",
      "NONPULMONARY_RECORDED"
    )
  ) %>%
  distinct(Run_ID, .keep_all = TRUE)

analysis_meta <- meta %>%
  inner_join(
    depth_df,
    by = c("Run_ID" = "run_id")
  )

write_csv(
  analysis_meta,
  file.path(
    DATA_OUT,
    "CRA002354_analysis_metadata_with_depth_and_source.csv"
  )
)

# ------------------------------------------------------------
# 10. Genus-level ecological summaries
# ------------------------------------------------------------

calc_alpha <- function(M) {

  observed <- rowSums(M > 0)

  shannon <- apply(M, 1, function(p) {
    p <- p[p > 0]
    -sum(p * log(p))
  })

  simpson <- apply(M, 1, function(p) {
    1 - sum(p^2)
  })

  tibble(
    Run_ID = rownames(M),
    Observed_Genera = observed,
    Shannon = shannon,
    Simpson = simpson
  )
}

alpha <- calc_alpha(genus_rel)

alpha_meta <- alpha %>%
  left_join(
    analysis_meta,
    by = "Run_ID"
  )

write_csv(
  alpha_meta,
  file.path(OUT, "06_GENUS_ALPHA_WITH_METADATA.csv")
)

# ------------------------------------------------------------
# 11. Depth-filtered analysis-ready exports
# ------------------------------------------------------------

keep2000 <- analysis_meta$Run_ID[
  analysis_meta$keep_ge2000
]

keep4000 <- analysis_meta$Run_ID[
  analysis_meta$keep_ge4000
]

rel2000 <- genus_rel_df %>%
  filter(Run_ID %in% keep2000)

rel4000 <- genus_rel_df %>%
  filter(Run_ID %in% keep4000)

meta2000 <- analysis_meta %>%
  filter(keep_ge2000)

meta4000 <- analysis_meta %>%
  filter(keep_ge4000)

write_csv(
  rel2000,
  file.path(
    DATA_OUT,
    "CRA002354_genus_relative_abundance_min2000.csv"
  )
)

write_csv(
  rel4000,
  file.path(
    DATA_OUT,
    "CRA002354_genus_relative_abundance_min4000.csv"
  )
)

write_csv(
  meta2000,
  file.path(
    DATA_OUT,
    "CRA002354_analysis_metadata_min2000.csv"
  )
)

write_csv(
  meta4000,
  file.path(
    DATA_OUT,
    "CRA002354_analysis_metadata_min4000.csv"
  )
)

# ------------------------------------------------------------
# 12. Final QC / readiness
# ------------------------------------------------------------

source_depth_qc <- analysis_meta %>%
  group_by(
    pulmonary_binary
  ) %>%
  summarise(
    samples_total = n(),
    patients_total = n_distinct(patient_id),
    samples_ge2000 = sum(keep_ge2000),
    patients_ge2000 =
      n_distinct(patient_id[keep_ge2000]),
    samples_ge4000 = sum(keep_ge4000),
    patients_ge4000 =
      n_distinct(patient_id[keep_ge4000]),
    median_assigned_reads =
      median(assigned_OTU97_reads),
    .groups = "drop"
  )

write_csv(
  source_depth_qc,
  file.path(OUT, "07_SOURCE_STRATIFIED_DEPTH_QC.csv")
)

patient_long_qc <- meta2000 %>%
  group_by(
    patient_id,
    pulmonary_binary
  ) %>%
  summarise(
    n_samples = n(),
    min_day = min(time_day, na.rm = TRUE),
    max_day = max(time_day, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  group_by(
    pulmonary_binary
  ) %>%
  summarise(
    patients = n(),
    patients_GE2 = sum(n_samples >= 2),
    patients_GE3 = sum(n_samples >= 3),
    .groups = "drop"
  )

write_csv(
  patient_long_qc,
  file.path(OUT, "08_SOURCE_STRATIFIED_LONGITUDINAL_QC_MIN2000.csv")
)

tax_summary <- tibble(
  OTU97_centroids = nrow(tax_df),
  OTUs_with_genus =
    sum(
      !is.na(tax_df$Genus) &
        tax_df$Genus != ""
    ),
  OTU_genus_assignment_fraction =
    mean(
      !is.na(tax_df$Genus) &
        tax_df$Genus != ""
    ),
  genus_labels_after_unclassified_preservation =
    n_distinct(tax_df$genus_label)
)

write_csv(
  tax_summary,
  file.path(OUT, "09_TAXONOMY_QC_SUMMARY.csv")
)

readiness <- tibble(
  unique_raw_samples = length(fa_files),
  total_deposited_reads = global_qc$total_reads,
  samples_ge2000_OTU_assigned =
    sum(depth_df$keep_ge2000),
  samples_ge4000_OTU_assigned =
    sum(depth_df$keep_ge4000),
  patients_ge2000 =
    n_distinct(meta2000$patient_id),
  pulmonary_patients_ge2000 =
    n_distinct(
      meta2000$patient_id[
        meta2000$pulmonary_binary == "PULMONARY"
      ]
    ),
  nonpulmonary_patients_ge2000 =
    n_distinct(
      meta2000$patient_id[
        meta2000$pulmonary_binary ==
          "NONPULMONARY_RECORDED"
      ]
    ),
  pulmonary_GE2_patients_min2000 =
    patient_long_qc$patients_GE2[
      patient_long_qc$pulmonary_binary ==
        "PULMONARY"
    ][1],
  nonpulmonary_GE2_patients_min2000 =
    patient_long_qc$patients_GE2[
      patient_long_qc$pulmonary_binary ==
        "NONPULMONARY_RECORDED"
    ][1],
  ready_for_source_longitudinal_analysis =
    (
      sum(depth_df$keep_ge2000) >= 120 &&
      n_distinct(meta2000$patient_id) >= 55 &&
      n_distinct(
        meta2000$patient_id[
          meta2000$pulmonary_binary == "PULMONARY"
        ]
      ) >= 30 &&
      n_distinct(
        meta2000$patient_id[
          meta2000$pulmonary_binary ==
            "NONPULMONARY_RECORDED"
        ]
      ) >= 15
    )
)

write_csv(
  readiness,
  file.path(OUT, "10_STEP94B1_ANALYSIS_READINESS.csv")
)

interpretation <- c(
  "STEP94B1 CRA002354 FASTA -> OTU97 -> SILVA138.2",
  "",
  paste0(
    "Canonical raw FASTA samples: ",
    length(fa_files)
  ),
  paste0(
    "Deposited reads audited: ",
    global_qc$total_reads
  ),
  paste0(
    "Weighted fraction 300-500 bp: ",
    round(global_qc$weighted_fraction_300_500, 5)
  ),
  paste0(
    "Weighted ambiguous-read fraction: ",
    round(global_qc$weighted_ambiguous_fraction, 6)
  ),
  "",
  "These CRA002354 data are FASTA without per-base quality scores, so the rescue pipeline intentionally does NOT call them DADA2 ASVs.",
  "The processed features are 97% OTUs, used only for within-cohort ecological/source analyses.",
  "Cross-cohort inference must continue to harmonize effect sizes/directions rather than merge feature tables.",
  "",
  paste0(
    "Ready for pulmonary-vs-nonpulmonary longitudinal analysis: ",
    readiness$ready_for_source_longitudinal_analysis
  )
)

writeLines(
  interpretation,
  file.path(OUT, "11_STEP94B1_INTERPRETATION.txt")
)

writeLines(
  c(
    paste0("Completed: ", Sys.time()),
    paste0("VSEARCH: ", VSEARCH),
    paste0("SILVA: ", SILVA_TRAIN),
    paste0(
      "Ready for source longitudinal analysis: ",
      readiness$ready_for_source_longitudinal_analysis
    ),
    "STEP94B1 COMPLETE"
  ),
  file.path(OUT, "_STEP94B1_COMPLETE.ok")
)

cat("STEP94B1 COMPLETE\n")
