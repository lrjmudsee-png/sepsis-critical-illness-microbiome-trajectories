# ============================================================
# Sepsis V2 - Step94B1B
# CRA002354 FASTA WRITE FIX + OTU97 RESUME
#
# Step94B1A filtering itself was successful:
# - 3,684,641 publication-aligned clean reads retained (62.31%)
# - 131/131 samples retained >=4000 clean reads
# - only ambiguous symbol detected was N
# - no >=8-bp homopolymer was detected
#
# But the pooled FASTA written by Step94B1A was empty to VSEARCH.
# Root cause: Biostrings::writeXStringSet(..., width=0).
#
# Biostrings defines width as the maximum number of letters per FASTA
# sequence line. Here we explicitly use width=80 and add hard file-size
# and VSEARCH-output sanity checks before proceeding.
#
# Pipeline:
# publication-aligned filtering
# -> correctly written pooled FASTA
# -> exact dereplication
# -> de novo chimera removal
# -> 97% OTU clustering
# -> read-to-OTU mapping
# -> SILVA 138.2 taxonomy (bootstrap >=80)
# -> genus counts/relative abundance
# -> source + longitudinal metadata
# -> >=4000 primary and >=2000 sensitivity datasets
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
  stop("R package 'dada2' is required.")
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
  "V2_34B1B_CRA002354_FASTA_WRITE_FIX_AND_OTU97_RESUME"
)

dir.create(DATA_OUT, recursive = TRUE, showWarnings = FALSE)
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

logfile <- file.path(OUT, "_STEP94B1B_runtime.txt")
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
# 1. VSEARCH
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
      ROOT, "tools", "vsearch",
      "vsearch.exe"
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

if (is.na(VSEARCH)) stop("VSEARCH not found.")

writeLines(
  system2(VSEARCH, "--version", stdout = TRUE, stderr = TRUE),
  file.path(OUT, "00_VSEARCH_VERSION.txt")
)

# ------------------------------------------------------------
# 2. SILVA 138.2
# ------------------------------------------------------------

find_silva <- function() {
  roots <- c(
    file.path(ROOT, "data", "_V2_ANALYSIS_READY", "00_REFERENCE"),
    file.path(ROOT, "reference"),
    file.path(ROOT, "references")
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

if (length(silva_hits) == 0) stop("No SILVA training set found.")

SILVA_TRAIN <- silva_hits[1]

if (!str_detect(
  tolower(SILVA_TRAIN),
  "138[._-]?2|1382|v138[._-]?2"
)) {
  stop("Preferred taxonomy reference is not clearly SILVA 138.2.")
}

# ------------------------------------------------------------
# 3. Raw files
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
# 4. Re-filter and CORRECTLY write pooled FASTA
# ------------------------------------------------------------

POOLED <- file.path(
  DATA_OUT,
  "CRA002354_publication_aligned_clean_pooled_WIDTH80.fa"
)

if (file.exists(POOLED)) file.remove(POOLED)

sample_qc <- vector("list", length(fa_files))
symbol_rows <- list()

first_write <- TRUE
global_counter <- 0L

for (i in seq_along(fa_files)) {

  f <- fa_files[i]
  run <- run_id_from_path(f)

  dna <- Biostrings::readDNAStringSet(f)

  n0 <- length(dna)
  widths <- Biostrings::width(dna)

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

  af_base <- Biostrings::alphabetFrequency(
    dna,
    baseOnly = TRUE
  )

  no_ambiguous <- af_base[, "other"] == 0
  length_ok <- widths >= 300 & widths <= 500

  hp8 <-
    Biostrings::vcountPattern("AAAAAAAA", dna, fixed = TRUE) > 0 |
    Biostrings::vcountPattern("CCCCCCCC", dna, fixed = TRUE) > 0 |
    Biostrings::vcountPattern("GGGGGGGG", dna, fixed = TRUE) > 0 |
    Biostrings::vcountPattern("TTTTTTTT", dna, fixed = TRUE) > 0

  keep <- length_ok & no_ambiguous & !hp8

  clean <- dna[keep]

  if (length(clean) > 0) {

    ids <- global_counter + seq_len(length(clean))
    global_counter <- global_counter + length(clean)

    names(clean) <- paste0(
      "SEQ", ids,
      ";sample=", run, ";"
    )

    # CRITICAL FIX: width=80, not width=0.
    Biostrings::writeXStringSet(
      clean,
      filepath = POOLED,
      format = "fasta",
      append = !first_write,
      compress = FALSE,
      width = 80
    )

    first_write <- FALSE
  }

  sample_qc[[i]] <- tibble(
    run_id = run,
    original_reads = n0,
    clean_reads = sum(keep),
    retention_fraction = sum(keep) / n0,
    discarded_length = sum(!length_ok),
    discarded_ambiguous = sum(!no_ambiguous),
    discarded_homopolymer8 = sum(hp8),
    clean_ge2000 = sum(keep) >= 2000,
    clean_ge4000 = sum(keep) >= 4000
  )

  rm(dna, clean, af_all, af_base)
  gc(FALSE)

  if (i %% 10 == 0 || i == length(fa_files)) {
    logmsg(
      "Re-filtered/wrote ",
      i,
      "/",
      length(fa_files),
      " samples"
    )
  }
}

sample_qc <- bind_rows(sample_qc)
symbol_audit <- bind_rows(symbol_rows)

write_csv(
  sample_qc,
  file.path(OUT, "02_FILTER_QC_BY_SAMPLE.csv")
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
  file.path(OUT, "03_GLOBAL_IUPAC_SYMBOL_AUDIT.csv")
)

global_qc <- sample_qc %>%
  summarise(
    samples = n(),
    original_reads = sum(original_reads),
    clean_reads = sum(clean_reads),
    retention_fraction = clean_reads / original_reads,
    min_clean_reads = min(clean_reads),
    median_clean_reads = median(clean_reads),
    max_clean_reads = max(clean_reads),
    samples_ge2000 = sum(clean_ge2000),
    samples_ge4000 = sum(clean_ge4000)
  )

pooled_bytes <- file.info(POOLED)$size

write_csv(
  global_qc %>%
    mutate(
      pooled_fasta_bytes = pooled_bytes,
      clean_records_written = global_counter
    ),
  file.path(OUT, "04_GLOBAL_FILTER_AND_FASTA_WRITE_QC.csv")
)

# Hard FASTA writer gate:
# 3.68M x ~300 nt alone should exceed 1 GB.
# Use a conservative lower bound to catch header-only/empty FASTA.
writer_gate <- tibble(
  criterion = c(
    "131 samples retain >=4000 clean reads",
    "clean FASTA records equal retained reads",
    "pooled FASTA > 500 MB"
  ),
  passed = c(
    global_qc$samples_ge4000 == 131,
    global_counter == global_qc$clean_reads,
    pooled_bytes > 500 * 1024^2
  ),
  observed = c(
    global_qc$samples_ge4000,
    global_counter,
    pooled_bytes
  )
)

write_csv(
  writer_gate,
  file.path(OUT, "05_FASTA_WRITER_GATE.csv")
)

if (!all(writer_gate$passed)) {
  writeLines(
    "Corrected pooled FASTA still failed the writer sanity gate.",
    file.path(OUT, "_STEP94B1B_FASTA_WRITE_FAILED.txt")
  )
  stop("Corrected FASTA writer gate failed.")
}

# ------------------------------------------------------------
# 5. VSEARCH runner
# ------------------------------------------------------------

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
        " with exit code ",
        status
      )
    )
  }

  logmsg("VSEARCH COMPLETE: ", label)

  invisible(
    list(
      stdout = stdout_file,
      stderr = stderr_file
    )
  )
}

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

run_vsearch(
  c(
    "--fastx_uniques", shQuote(POOLED),
    "--fastaout", shQuote(UNIQUES),
    "--sizeout",
    "--relabel", "Uniq",
    "--fasta_width", "80"
  ),
  "01_DEREPLICATION"
)

if (
  !file.exists(UNIQUES) ||
  file.info(UNIQUES)$size < 1024
) {
  stop(
    "Dereplication output is empty. Inspect VSEARCH_01_DEREPLICATION.stderr.txt."
  )
}

run_vsearch(
  c(
    "--uchime3_denovo", shQuote(UNIQUES),
    "--nonchimeras", shQuote(NONCHIM),
    "--sizein",
    "--sizeout",
    "--fasta_width", "80"
  ),
  "02_CHIMERA_REMOVAL"
)

if (
  !file.exists(NONCHIM) ||
  file.info(NONCHIM)$size < 1024
) {
  stop(
    "Non-chimeric FASTA is empty. Inspect chimera-removal log."
  )
}

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
    "--fasta_width", "80"
  ),
  "03_OTU97_CLUSTERING"
)

if (
  !file.exists(OTUS) ||
  file.info(OTUS)$size < 1024
) {
  stop("OTU centroid FASTA is empty.")
}

run_vsearch(
  c(
    "--usearch_global", shQuote(POOLED),
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

if (
  !file.exists(OTU_TABLE) ||
  file.info(OTU_TABLE)$size < 100
) {
  stop("OTU table is empty after read mapping.")
}

# ------------------------------------------------------------
# 6. Read OTU table
# ------------------------------------------------------------

otu <- read_tsv(
  OTU_TABLE,
  show_col_types = FALSE,
  name_repair = "minimal"
)

if (nrow(otu) == 0 || ncol(otu) < 2) {
  stop("OTU table has no usable data.")
}

names(otu)[1] <- "otu_id"

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
  assigned_OTU97_reads = as.numeric(assigned_depth),
  keep_ge2000 = assigned_OTU97_reads >= 2000,
  keep_ge4000 = assigned_OTU97_reads >= 4000
)

write_csv(
  depth_df,
  file.path(OUT, "06_OTU97_ASSIGNED_DEPTH.csv")
)

# ------------------------------------------------------------
# 7. SILVA taxonomy
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
    "Kingdom", "Phylum", "Class",
    "Order", "Family", "Genus"
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
# 9. Frozen metadata join
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
    baseline_sepsis_status = norm_chr(baseline_sepsis_status),
    sample_sepsis_status = norm_chr(sample_sepsis_status),
    infection_source_raw = norm_chr(infection_source_standard),
    outcome_28d = norm_chr(outcome_28d),
    baseline_sofa = suppressWarnings(as.numeric(baseline_sofa)),
    baseline_apache_ii = suppressWarnings(as.numeric(baseline_apache_ii)),
    baseline_lactate = suppressWarnings(as.numeric(baseline_lactate))
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
# 10. Primary and sensitivity exports
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
# 11. QC and readiness
# ------------------------------------------------------------

tax_qc <- tibble(
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
  tax_qc,
  file.path(OUT, "07_TAXONOMY_QC_SUMMARY.csv")
)

source_depth <- analysis_meta %>%
  group_by(pulmonary_binary) %>%
  summarise(
    samples_total = n(),
    patients_total = n_distinct(patient_id),
    samples_ge4000 = sum(keep_ge4000),
    patients_ge4000 = n_distinct(patient_id[keep_ge4000]),
    samples_ge2000 = sum(keep_ge2000),
    patients_ge2000 = n_distinct(patient_id[keep_ge2000]),
    median_assigned_reads = median(assigned_OTU97_reads),
    .groups = "drop"
  )

write_csv(
  source_depth,
  file.path(OUT, "08_SOURCE_STRATIFIED_DEPTH_QC.csv")
)

patient_long <- analysis_meta %>%
  filter(keep_ge4000) %>%
  group_by(
    patient_id,
    pulmonary_binary
  ) %>%
  summarise(
    n_samples = n(),
    min_day = min(time_day, na.rm = TRUE),
    max_day = max(time_day, na.rm = TRUE),
    .groups = "drop"
  )

long_qc <- patient_long %>%
  group_by(pulmonary_binary) %>%
  summarise(
    patients = n(),
    patients_GE2 = sum(n_samples >= 2),
    patients_GE3 = sum(n_samples >= 3),
    .groups = "drop"
  )

write_csv(
  long_qc,
  file.path(OUT, "09_PRIMARY_MIN4000_LONGITUDINAL_DEPTH_BY_SOURCE.csv")
)

get0 <- function(df, group, col) {
  x <- df[df$pulmonary_binary == group, col, drop = TRUE]
  if (length(x) == 0) 0 else x[1]
}

readiness <- tibble(
  clean_reads = global_qc$clean_reads,
  clean_read_retention = global_qc$retention_fraction,
  samples_assigned_ge4000 = sum(depth_df$keep_ge4000),
  samples_assigned_ge2000 = sum(depth_df$keep_ge2000),
  patients_assigned_ge4000 =
    n_distinct(
      analysis_meta$patient_id[
        analysis_meta$keep_ge4000
      ]
    ),
  pulmonary_patients_ge4000 =
    get0(source_depth, "PULMONARY", "patients_ge4000"),
  nonpulmonary_patients_ge4000 =
    get0(source_depth, "NONPULMONARY_RECORDED", "patients_ge4000"),
  pulmonary_GE2_patients =
    get0(long_qc, "PULMONARY", "patients_GE2"),
  nonpulmonary_GE2_patients =
    get0(long_qc, "NONPULMONARY_RECORDED", "patients_GE2"),
  ready_for_step94B2 =
    (
      sum(depth_df$keep_ge4000) >= 115 &&
      n_distinct(
        analysis_meta$patient_id[
          analysis_meta$keep_ge4000
        ]
      ) >= 55 &&
      get0(source_depth, "PULMONARY", "patients_ge4000") >= 30 &&
      get0(source_depth, "NONPULMONARY_RECORDED", "patients_ge4000") >= 15 &&
      get0(long_qc, "PULMONARY", "patients_GE2") >= 20 &&
      get0(long_qc, "NONPULMONARY_RECORDED", "patients_GE2") >= 10
    )
)

write_csv(
  readiness,
  file.path(OUT, "10_STEP94B1B_ANALYSIS_READINESS.csv")
)

writeLines(
  c(
    "STEP94B1B FASTA WRITE FIX + OTU97 RESUME",
    "",
    "Step94B1A biological filtering was valid; the failure was technical.",
    "Only N was detected as an ambiguous IUPAC symbol.",
    "The pooled FASTA is now written with width=80 and validated before VSEARCH.",
    "",
    paste0(
      "Clean reads: ",
      global_qc$clean_reads
    ),
    paste0(
      "Samples >=4000 OTU-assigned reads: ",
      sum(depth_df$keep_ge4000)
    ),
    paste0(
      "Ready for Step94B2: ",
      readiness$ready_for_step94B2
    )
  ),
  file.path(OUT, "11_STEP94B1B_INTERPRETATION.txt")
)

writeLines(
  c(
    paste0("Completed: ", Sys.time()),
    paste0("Ready for Step94B2: ", readiness$ready_for_step94B2),
    "STEP94B1B COMPLETE"
  ),
  file.path(OUT, "_STEP94B1B_COMPLETE.ok")
)

cat("STEP94B1B COMPLETE\n")
