# ============================================================
# Sepsis V2 - Step94B1A
# CRA002354 ambiguity resolution + publication-aligned FASTA QC
# + 97% OTU + SILVA 138.2
#
# Step94B1 stopped because 37% of deposited reads contained
# non-ACGT symbols. The original publication required:
#   - 300-500 bp
#   - no ambiguous bases
#   - homopolymers <8 bp
#   - Q20
#
# FASTA has no per-base Q scores, so Q20 cannot be re-applied.
# This script reproduces every sequence-level rule that is
# recoverable from FASTA:
#   1) length 300-500 bp
#   2) A/C/G/T only
#   3) no homopolymer run >=8 bp
#
# It first audits exact ambiguous IUPAC symbols, then filters.
# Only filtered reads enter VSEARCH OTU clustering.
#
# Features are 97% OTUs, NOT ASVs.
# ============================================================

options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(stringr)
  library(tibble)
})

if (!requireNamespace("Biostrings", quietly = TRUE)) {
  stop("Bioconductor package 'Biostrings' is required.")
}

if (!requireNamespace("dada2", quietly = TRUE)) {
  stop("R package 'dada2' is required for SILVA taxonomy assignment.")
}

ROOT <- "E:/sepsis_project"

RAW_DIR <- file.path(
  ROOT, "data", "_V2_ANALYSIS_READY",
  "01_PROJECTS", "CRA002354", "01_raw"
)

MASTER <- file.path(
  ROOT, "data", "_V2_ANALYSIS_READY", "00_FREEZE",
  "V2_FINAL_FROZEN_master_metadata_CLINICAL_REPAIRED_20260814_131620.csv"
)

DATA_OUT <- file.path(
  ROOT, "data", "CRA002354",
  "03_vsearch97_silva1382"
)

OUT <- file.path(
  ROOT, "results",
  "V2_34B1A_CRA002354_AMBIGUITY_RESOLUTION_OTU97"
)

dir.create(DATA_OUT, recursive = TRUE, showWarnings = FALSE)
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

logfile <- file.path(OUT, "_STEP94B1A_runtime.txt")
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
# 1. VSEARCH detection
# ------------------------------------------------------------

find_vsearch <- function() {

  candidates <- c(
    Sys.which("vsearch"),
    Sys.which("vsearch.exe"),
    file.path(
      ROOT, "tools",
      "vsearch-2.31.0-win-x86_64",
      "bin", "vsearch.exe"
    ),
    file.path(
      ROOT, "tools", "vsearch", "vsearch.exe"
    )
  )

  candidates <- unique(candidates[nzchar(candidates)])
  candidates <- candidates[file.exists(candidates)]

  if (length(candidates) == 0) return(NA_character_)

  normalizePath(
    candidates[1],
    winslash = "/",
    mustWork = TRUE
  )
}

VSEARCH <- find_vsearch()

if (is.na(VSEARCH)) {
  stop("VSEARCH not found. Run GET_VSEARCH_2.31.0.ps1 first.")
}

writeLines(
  system2(VSEARCH, "--version", stdout = TRUE, stderr = TRUE),
  file.path(OUT, "00_VSEARCH_VERSION.txt")
)

# ------------------------------------------------------------
# 2. SILVA 138.2 detection
# ------------------------------------------------------------

find_silva <- function() {

  roots <- c(
    file.path(ROOT, "data", "_V2_ANALYSIS_READY", "00_REFERENCE"),
    file.path(ROOT, "reference"),
    file.path(ROOT, "references"),
    file.path(ROOT, "SILVA_138_2")
  )

  roots <- unique(roots[dir.exists(roots)])

  hits <- character()

  for (rr in roots) {
    hits <- c(
      hits,
      list.files(
        rr,
        recursive = TRUE,
        full.names = TRUE,
        pattern = "silva.*train.*set.*\\.fa(\\.gz)?$",
        ignore.case = TRUE
      )
    )
  }

  hits <- unique(hits)

  if (length(hits) == 0) return(character())

  z <- tolower(hits)

  score <-
    100 * str_detect(z, "138[._-]?2|1382|v138[._-]?2") +
    50 * str_detect(z, "nr99") +
    20 * str_detect(z, "train") -
    100 * str_detect(z, "138[._-]?1|v138[._-]?1")

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
  stop("No SILVA training set found.")
}

SILVA_TRAIN <- silva_hits[1]

if (!str_detect(
  tolower(SILVA_TRAIN),
  "138[._-]?2|1382|v138[._-]?2"
)) {
  stop("Preferred taxonomy reference is not clearly SILVA 138.2.")
}

# ------------------------------------------------------------
# 3. Canonical raw files
# ------------------------------------------------------------

fa_files <- list.files(
  RAW_DIR,
  full.names = TRUE,
  pattern = "^CRR\\d+\\.fa\\.gz$",
  ignore.case = TRUE
)

fa_files <- fa_files[
  !duplicated(run_id_from_path(fa_files))
]

fa_files <- fa_files[
  order(run_id_from_path(fa_files))
]

if (length(fa_files) != 131) {
  stop(
    paste0(
      "Expected 131 unique CRR*.fa.gz files; found ",
      length(fa_files)
    )
  )
}

# ------------------------------------------------------------
# 4. Exact ambiguity audit + publication-aligned filtering
# ------------------------------------------------------------

POOLED_CLEAN <- file.path(
  DATA_OUT,
  "CRA002354_publication_aligned_clean_pooled.fa"
)

if (file.exists(POOLED_CLEAN)) {
  file.remove(POOLED_CLEAN)
}

sample_qc <- vector("list", length(fa_files))
symbol_rows <- list()

first_write <- TRUE

for (i in seq_along(fa_files)) {

  f <- fa_files[i]
  run <- run_id_from_path(f)

  dna <- Biostrings::readDNAStringSet(f)

  n0 <- length(dna)

  if (n0 == 0) {
    stop(paste0("No FASTA records in ", basename(f)))
  }

  widths <- Biostrings::width(dna)

  # Exact IUPAC alphabet counts.
  af_all <- Biostrings::alphabetFrequency(
    dna,
    baseOnly = FALSE
  )

  totals <- colSums(af_all)

  symbol_rows[[length(symbol_rows) + 1]] <- tibble(
    run_id = run,
    symbol = names(totals),
    base_count = as.numeric(totals)
  ) %>%
    filter(base_count > 0)

  # "other" in baseOnly=TRUE includes all non-ACGT symbols.
  af_base <- Biostrings::alphabetFrequency(
    dna,
    baseOnly = TRUE
  )

  no_ambiguous <- af_base[, "other"] == 0

  length_ok <- widths >= 300 & widths <= 500

  # Publication criterion: homopolymers shorter than 8 bp.
  hp8 <-
    Biostrings::vcountPattern(
      "AAAAAAAA",
      dna,
      fixed = TRUE
    ) > 0 |
    Biostrings::vcountPattern(
      "CCCCCCCC",
      dna,
      fixed = TRUE
    ) > 0 |
    Biostrings::vcountPattern(
      "GGGGGGGG",
      dna,
      fixed = TRUE
    ) > 0 |
    Biostrings::vcountPattern(
      "TTTTTTTT",
      dna,
      fixed = TRUE
    ) > 0

  homopolymer_ok <- !hp8

  keep <-
    length_ok &
    no_ambiguous &
    homopolymer_ok

  kept <- sum(keep)

  # Create explicit VSEARCH sample tags.
  clean <- dna[keep]

  names(clean) <- paste0(
    "SEQ",
    seq_len(length(clean)),
    "_",
    run,
    ";sample=",
    run,
    ";"
  )

  Biostrings::writeXStringSet(
    clean,
    filepath = POOLED_CLEAN,
    format = "fasta",
    append = !first_write,
    compress = FALSE,
    width = 0
  )

  first_write <- FALSE

  sample_qc[[i]] <- tibble(
    run_id = run,
    original_reads = n0,
    reads_length_300_500 = sum(length_ok),
    reads_no_ambiguous = sum(no_ambiguous),
    reads_no_homopolymer8 = sum(homopolymer_ok),
    reads_publication_aligned_clean = kept,
    retention_fraction = kept / n0,
    discarded_length = sum(!length_ok),
    discarded_ambiguous = sum(!no_ambiguous),
    discarded_homopolymer8 = sum(hp8),
    clean_ge2000 = kept >= 2000,
    clean_ge4000 = kept >= 4000
  )

  rm(dna, clean, af_all, af_base)
  gc(FALSE)

  if (i %% 10 == 0 || i == length(fa_files)) {
    logmsg("Cleaned ", i, "/", length(fa_files), " samples")
  }
}

sample_qc <- bind_rows(sample_qc)
symbol_audit <- bind_rows(symbol_rows)

write_csv(
  sample_qc,
  file.path(OUT, "02_PUBLICATION_ALIGNED_FILTER_QC_BY_SAMPLE.csv")
)

write_csv(
  symbol_audit,
  file.path(OUT, "03_IUPAC_SYMBOL_COUNTS_BY_SAMPLE.csv")
)

symbol_global <- symbol_audit %>%
  group_by(symbol) %>%
  summarise(
    base_count = sum(base_count),
    .groups = "drop"
  ) %>%
  mutate(
    fraction_of_all_bases =
      base_count / sum(base_count)
  ) %>%
  arrange(desc(base_count))

write_csv(
  symbol_global,
  file.path(OUT, "04_GLOBAL_IUPAC_SYMBOL_AUDIT.csv")
)

filter_global <- sample_qc %>%
  summarise(
    samples = n(),
    original_reads = sum(original_reads),
    clean_reads = sum(reads_publication_aligned_clean),
    overall_retention_fraction =
      clean_reads / original_reads,
    min_clean_reads =
      min(reads_publication_aligned_clean),
    median_clean_reads =
      median(reads_publication_aligned_clean),
    max_clean_reads =
      max(reads_publication_aligned_clean),
    samples_clean_ge2000 =
      sum(clean_ge2000),
    samples_clean_ge4000 =
      sum(clean_ge4000)
  )

write_csv(
  filter_global,
  file.path(OUT, "05_GLOBAL_FILTER_QC.csv")
)

# Gate after REAL filtering, not estimated filtering.
filter_gate <- tibble(
  criterion = c(
    ">=120 samples retain >=4000 publication-aligned reads",
    "All samples retain >=2000 publication-aligned reads",
    "Overall clean-read retention >=50%"
  ),
  passed = c(
    filter_global$samples_clean_ge4000 >= 120,
    filter_global$samples_clean_ge2000 == 131,
    filter_global$overall_retention_fraction >= 0.50
  ),
  observed = c(
    filter_global$samples_clean_ge4000,
    filter_global$samples_clean_ge2000,
    filter_global$overall_retention_fraction
  )
)

write_csv(
  filter_gate,
  file.path(OUT, "06_FILTER_GATE.csv")
)

if (!all(filter_gate$passed)) {
  writeLines(
    c(
      "Publication-aligned FASTA filter gate failed.",
      "Review outputs 02-06 before OTU processing."
    ),
    file.path(OUT, "_STEP94B1A_FILTER_REVIEW_REQUIRED.txt")
  )
  stop("Publication-aligned filter gate failed.")
}

# ------------------------------------------------------------
# 5. VSEARCH 97% OTU pipeline
# ------------------------------------------------------------

UNIQUES <- file.path(
  DATA_OUT,
  "CRA002354_clean_uniques.fa"
)

NONCHIM <- file.path(
  DATA_OUT,
  "CRA002354_clean_uniques_nonchimeric.fa"
)

OTUS <- file.path(
  DATA_OUT,
  "CRA002354_OTU97_centroids.fa"
)

OTU_TABLE <- file.path(
  DATA_OUT,
  "CRA002354_OTU97_table.tsv"
)

THREADS <- max(
  1L,
  min(
    12L,
    parallel::detectCores(logical = TRUE) - 1L
  )
)

run_vsearch <- function(args, label) {

  stdout_file <- file.path(
    OUT,
    paste0("VSEARCH_", label, ".stdout.txt")
  )

  stderr_file <- file.path(
    OUT,
    paste0("VSEARCH_", label, ".stderr.txt")
  )

  logmsg("VSEARCH START: ", label)

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
        "; exit code ",
        status
      )
    )
  }

  logmsg("VSEARCH COMPLETE: ", label)
}

# Exact dereplication. Sample labels are not needed in unique centroids;
# sample identity is retained in POOLED_CLEAN for later mapping.
run_vsearch(
  c(
    "--fastx_uniques", shQuote(POOLED_CLEAN),
    "--fastaout", shQuote(UNIQUES),
    "--sizeout",
    "--relabel", "Uniq",
    "--fasta_width", "0"
  ),
  "01_DEREPLICATION"
)

# De novo chimera removal.
run_vsearch(
  c(
    "--uchime3_denovo", shQuote(UNIQUES),
    "--nonchimeras", shQuote(NONCHIM),
    "--sizein",
    "--sizeout",
    "--fasta_width", "0"
  ),
  "02_CHIMERA_REMOVAL"
)

# 97% OTU clustering.
run_vsearch(
  c(
    "--cluster_size", shQuote(NONCHIM),
    "--id", "0.97",
    "--centroids", shQuote(OTUS),
    "--sizein",
    "--sizeout",
    "--relabel", "OTU_",
    "--threads", as.character(THREADS),
    "--qmask", "none",
    "--fasta_width", "0"
  ),
  "03_OTU97_CLUSTERING"
)

# Map filtered sample-labelled reads back to OTUs.
run_vsearch(
  c(
    "--usearch_global", shQuote(POOLED_CLEAN),
    "--db", shQuote(OTUS),
    "--id", "0.97",
    "--strand", "plus",
    "--otutabout", shQuote(OTU_TABLE),
    "--threads", as.character(THREADS),
    "--maxaccepts", "1",
    "--maxrejects", "64",
    "--qmask", "none",
    "--dbmask", "none"
  ),
  "04_MAP_READS_TO_OTUS"
)

# ------------------------------------------------------------
# 6. OTU table
# ------------------------------------------------------------

otu <- read_tsv(
  OTU_TABLE,
  show_col_types = FALSE,
  name_repair = "minimal"
)

if (nrow(otu) == 0 || ncol(otu) < 2) {
  stop("OTU table is empty.")
}

names(otu)[1] <- "otu_id"

# Defensive removal if VSEARCH adds taxonomy column.
sample_cols <- setdiff(
  names(otu),
  c("otu_id", "taxonomy")
)

counts <- as.matrix(
  otu[, sample_cols, drop = FALSE]
)

storage.mode(counts) <- "numeric"
rownames(counts) <- otu$otu_id

assigned_depth <- colSums(counts)

depth_df <- tibble(
  Run_ID = names(assigned_depth),
  assigned_OTU97_reads =
    as.numeric(assigned_depth),
  keep_ge2000 =
    assigned_OTU97_reads >= 2000,
  keep_ge4000 =
    assigned_OTU97_reads >= 4000
)

write_csv(
  depth_df,
  file.path(OUT, "07_OTU97_ASSIGNED_DEPTH.csv")
)

# ------------------------------------------------------------
# 7. SILVA 138.2 taxonomy, bootstrap >=80
# ------------------------------------------------------------

dna_otus <- Biostrings::readDNAStringSet(OTUS)

otu_seq_ids <- sub(
  ";.*$",
  "",
  names(dna_otus)
)

seqs <- as.character(dna_otus)
names(seqs) <- otu_seq_ids

logmsg(
  "Assigning SILVA 138.2 taxonomy to ",
  length(seqs),
  " OTU centroids"
)

tax <- dada2::assignTaxonomy(
  seqs,
  refFasta = SILVA_TRAIN,
  multithread = TRUE,
  tryRC = TRUE,
  minBoot = 80
)

tax_df <- as.data.frame(
  tax,
  stringsAsFactors = FALSE
) %>%
  rownames_to_column("otu_id")

write_csv(
  tax_df,
  file.path(
    DATA_OUT,
    "CRA002354_OTU97_SILVA1382_taxonomy_minBoot80.csv"
  )
)

# ------------------------------------------------------------
# 8. Genus aggregation
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

  if (
    "Genus" %in% names(vals) &&
    !is.na(vals["Genus"]) &&
    nzchar(vals["Genus"])
  ) {
    return(vals["Genus"])
  }

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

genus_for_otu <- tax_df$genus_label[
  match(
    rownames(counts),
    tax_df$otu_id
  )
]

genus_for_otu[
  is.na(genus_for_otu) |
    genus_for_otu == ""
] <- "Unclassified_Bacteria"

genus_counts <- rowsum(
  counts,
  group = genus_for_otu,
  reorder = FALSE
)

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
# 9. Join frozen metadata
# ------------------------------------------------------------

master <- read_csv(
  MASTER,
  show_col_types = FALSE,
  guess_max = 50000,
  name_repair = "unique"
)

meta <- master %>%
  filter(toupper(project) == "CRA002354") %>%
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
    outcome_28d =
      norm_chr(outcome_28d),
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
      tolower(infection_source_raw) == "surgical_site" ~ "SKIN_SOFT_TISSUE",
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
    by = "Run_ID"
  )

write_csv(
  analysis_meta,
  file.path(
    DATA_OUT,
    "CRA002354_analysis_metadata_with_depth_and_source.csv"
  )
)

# ------------------------------------------------------------
# 10. Primary >=4000 + sensitivity >=2000 exports
# ------------------------------------------------------------

keep4000 <- analysis_meta$Run_ID[
  analysis_meta$keep_ge4000
]

keep2000 <- analysis_meta$Run_ID[
  analysis_meta$keep_ge2000
]

write_csv(
  genus_rel_df %>%
    filter(Run_ID %in% keep4000),
  file.path(
    DATA_OUT,
    "CRA002354_genus_relative_abundance_PRIMARY_min4000.csv"
  )
)

write_csv(
  analysis_meta %>%
    filter(keep_ge4000),
  file.path(
    DATA_OUT,
    "CRA002354_analysis_metadata_PRIMARY_min4000.csv"
  )
)

write_csv(
  genus_rel_df %>%
    filter(Run_ID %in% keep2000),
  file.path(
    DATA_OUT,
    "CRA002354_genus_relative_abundance_SENSITIVITY_min2000.csv"
  )
)

write_csv(
  analysis_meta %>%
    filter(keep_ge2000),
  file.path(
    DATA_OUT,
    "CRA002354_analysis_metadata_SENSITIVITY_min2000.csv"
  )
)

# ------------------------------------------------------------
# 11. Taxonomy / source / longitudinal readiness QC
# ------------------------------------------------------------

tax_summary <- tibble(
  OTU97_centroids = nrow(tax_df),
  OTUs_with_genus =
    sum(
      !is.na(tax_df$Genus) &
        tax_df$Genus != ""
    ),
  genus_assignment_fraction =
    mean(
      !is.na(tax_df$Genus) &
        tax_df$Genus != ""
    ),
  genus_labels_preserving_unclassified =
    n_distinct(tax_df$genus_label)
)

write_csv(
  tax_summary,
  file.path(OUT, "08_TAXONOMY_QC_SUMMARY.csv")
)

source_depth <- analysis_meta %>%
  group_by(pulmonary_binary) %>%
  summarise(
    samples_total = n(),
    patients_total =
      n_distinct(patient_id),
    samples_ge4000 =
      sum(keep_ge4000),
    patients_ge4000 =
      n_distinct(patient_id[keep_ge4000]),
    samples_ge2000 =
      sum(keep_ge2000),
    patients_ge2000 =
      n_distinct(patient_id[keep_ge2000]),
    median_assigned_reads =
      median(assigned_OTU97_reads),
    .groups = "drop"
  )

write_csv(
  source_depth,
  file.path(OUT, "09_SOURCE_STRATIFIED_DEPTH_QC.csv")
)

patient_long <- analysis_meta %>%
  filter(keep_ge4000) %>%
  group_by(
    patient_id,
    pulmonary_binary
  ) %>%
  summarise(
    n_samples = n(),
    min_day =
      min(time_day, na.rm = TRUE),
    max_day =
      max(time_day, na.rm = TRUE),
    .groups = "drop"
  )

long_summary <- patient_long %>%
  group_by(pulmonary_binary) %>%
  summarise(
    patients = n(),
    patients_GE2 =
      sum(n_samples >= 2),
    patients_GE3 =
      sum(n_samples >= 3),
    .groups = "drop"
  )

write_csv(
  long_summary,
  file.path(OUT, "10_PRIMARY_MIN4000_LONGITUDINAL_DEPTH_BY_SOURCE.csv")
)

pulm_n <- source_depth$patients_ge4000[
  source_depth$pulmonary_binary == "PULMONARY"
]

nonpulm_n <- source_depth$patients_ge4000[
  source_depth$pulmonary_binary == "NONPULMONARY_RECORDED"
]

pulm_ge2 <- long_summary$patients_GE2[
  long_summary$pulmonary_binary == "PULMONARY"
]

nonpulm_ge2 <- long_summary$patients_GE2[
  long_summary$pulmonary_binary == "NONPULMONARY_RECORDED"
]

readiness <- tibble(
  FASTA_samples = 131,
  clean_reads =
    filter_global$clean_reads,
  clean_read_retention =
    filter_global$overall_retention_fraction,
  samples_clean_ge4000 =
    filter_global$samples_clean_ge4000,
  OTU97_centroids =
    nrow(tax_df),
  OTUs_with_genus =
    tax_summary$OTUs_with_genus,
  samples_assigned_ge4000 =
    sum(depth_df$keep_ge4000),
  patients_assigned_ge4000 =
    n_distinct(
      analysis_meta$patient_id[
        analysis_meta$keep_ge4000
      ]
    ),
  pulmonary_patients_ge4000 =
    ifelse(length(pulm_n) == 0, 0, pulm_n),
  nonpulmonary_patients_ge4000 =
    ifelse(length(nonpulm_n) == 0, 0, nonpulm_n),
  pulmonary_GE2_patients =
    ifelse(length(pulm_ge2) == 0, 0, pulm_ge2),
  nonpulmonary_GE2_patients =
    ifelse(length(nonpulm_ge2) == 0, 0, nonpulm_ge2),
  ready_for_step94B2 =
    (
      sum(depth_df$keep_ge4000) >= 115 &&
      n_distinct(
        analysis_meta$patient_id[
          analysis_meta$keep_ge4000
        ]
      ) >= 55 &&
      ifelse(length(pulm_n) == 0, 0, pulm_n) >= 30 &&
      ifelse(length(nonpulm_n) == 0, 0, nonpulm_n) >= 15 &&
      ifelse(length(pulm_ge2) == 0, 0, pulm_ge2) >= 20 &&
      ifelse(length(nonpulm_ge2) == 0, 0, nonpulm_ge2) >= 10
    )
)

write_csv(
  readiness,
  file.path(OUT, "11_STEP94B1A_ANALYSIS_READINESS.csv")
)

interpretation <- c(
  "STEP94B1A CRA002354 AMBIGUITY RESOLUTION",
  "",
  paste0(
    "Original deposited reads: ",
    filter_global$original_reads
  ),
  paste0(
    "Publication-aligned clean reads retained: ",
    filter_global$clean_reads
  ),
  paste0(
    "Retention fraction: ",
    round(
      filter_global$overall_retention_fraction,
      4
    )
  ),
  paste0(
    "Samples retaining >=4000 clean reads before OTU mapping: ",
    filter_global$samples_clean_ge4000,
    "/131"
  ),
  "",
  "The deposited FASTA cannot reproduce the original Q20 filtering because quality scores are absent.",
  "Accordingly, only FASTA-recoverable publication criteria were re-applied: 300-500 bp, no ambiguous bases, and homopolymers shorter than 8 bp.",
  "",
  paste0(
    "Ready for Step94B2 pulmonary-vs-nonpulmonary longitudinal ecology: ",
    readiness$ready_for_step94B2
  )
)

writeLines(
  interpretation,
  file.path(OUT, "12_STEP94B1A_INTERPRETATION.txt")
)

writeLines(
  c(
    paste0("Completed: ", Sys.time()),
    paste0("SILVA: ", SILVA_TRAIN),
    paste0(
      "Ready for Step94B2: ",
      readiness$ready_for_step94B2
    ),
    "STEP94B1A COMPLETE"
  ),
  file.path(OUT, "_STEP94B1A_COMPLETE.ok")
)

cat("STEP94B1A COMPLETE\n")
