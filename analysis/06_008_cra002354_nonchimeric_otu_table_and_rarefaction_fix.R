# ============================================================
# Sepsis V2 - Step94B1D
# CRA002354 NON-CHIMERIC OTU TABLE + RAREFACTION FIX
#
# Why this step is required:
# Step94B1B correctly identified/removes chimeric UNIQUE sequences before
# OTU clustering, but then mapped ALL publication-aligned clean reads
# back to non-chimeric OTU centroids.
#
# B1B logs show:
#   total clean reads       = 3,684,641
#   abundance-weighted chimera reads = 207,074 (5.6%)
#   expected non-chimera reads       = 3,477,567
#   B1B OTU mapping counted          = 3,684,290 (99.99% of ALL clean reads)
#
# Therefore many reads whose exact sequence was called chimera could map
# at >=97% to a non-chimeric OTU and re-enter the abundance table.
#
# This correction:
# 1) exact-matches sample-labelled clean reads to the NONCHIMERIC UNIQUE
#    sequence database;
# 2) retains only exact non-chimeric reads;
# 3) maps those reads to the EXISTING 97% non-chimeric OTU centroids;
# 4) rebuilds sample x OTU counts;
# 5) audits corrected depth/source balance;
# 6) calculates abundance-weighted SILVA taxonomy coverage;
# 7) rebuilds genus counts/relative abundance for exploratory taxonomy;
# 8) rarefies OTU counts to 4000 reads for PRIMARY alpha diversity;
# 9) rarefies to 2000 for V2 depth sensitivity.
#
# Existing OTU centroids and SILVA taxonomy are NOT recomputed.
# ============================================================

options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(stringr)
  library(tibble)
})

if (!requireNamespace("vegan", quietly = TRUE)) {
  stop("R package 'vegan' is required.")
}

ROOT <- "E:/sepsis_project"

DATA_DIR <- file.path(
  ROOT,
  "data",
  "CRA002354",
  "03_vsearch97_silva1382"
)

B1B_OUT <- file.path(
  ROOT,
  "results",
  "V2_34B1B_CRA002354_FASTA_WRITE_FIX_AND_OTU97_RESUME"
)

OUT <- file.path(
  ROOT,
  "results",
  "V2_34B1D_CRA002354_NONCHIMERIC_OTU_TABLE_AND_RAREFACTION_FIX"
)

dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

POOLED <- file.path(
  DATA_DIR,
  "CRA002354_publication_aligned_clean_pooled_WIDTH80.fa"
)

NONCHIM_UNIQUE <- file.path(
  DATA_DIR,
  "CRA002354_clean_uniques_nonchimeric.fa"
)

OTUS <- file.path(
  DATA_DIR,
  "CRA002354_OTU97_centroids.fa"
)

TAX_FILE <- file.path(
  DATA_DIR,
  "CRA002354_OTU97_SILVA1382_taxonomy_minBoot80.csv"
)

META_FILE <- file.path(
  DATA_DIR,
  "CRA002354_analysis_metadata_with_depth_and_source.csv"
)

NONCHIM_READS <- file.path(
  DATA_DIR,
  "CRA002354_sample_labeled_EXACT_NONCHIMERIC_reads.fa"
)

CORRECTED_OTU_TABLE <- file.path(
  DATA_DIR,
  "CRA002354_OTU97_table_NONCHIMERIC_CORRECTED.tsv"
)

VSEARCH <- file.path(
  ROOT,
  "tools",
  "vsearch-2.31.0-win-x86_64",
  "bin",
  "vsearch.exe"
)

if (!file.exists(VSEARCH)) {
  p <- Sys.which("vsearch")
  if (!nzchar(p)) p <- Sys.which("vsearch.exe")
  if (!nzchar(p)) stop("VSEARCH not found.")
  VSEARCH <- p
}

required_files <- c(
  POOLED,
  NONCHIM_UNIQUE,
  OTUS,
  TAX_FILE,
  META_FILE
)

missing <- required_files[!file.exists(required_files)]

if (length(missing) > 0) {
  stop(
    paste0(
      "Missing required B1B product(s): ",
      paste(missing, collapse = "; ")
    )
  )
}

# Confirm B1B completed.
if (!file.exists(file.path(B1B_OUT, "_STEP94B1B_COMPLETE.ok"))) {
  stop("Step94B1B complete flag is missing.")
}

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

  cat("VSEARCH START:", label, "\n")

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

  cat("VSEARCH COMPLETE:", label, "\n")
}

# ------------------------------------------------------------
# 1. Keep ONLY reads whose exact sequence is present in the
#    non-chimeric unique sequence database.
# ------------------------------------------------------------

if (file.exists(NONCHIM_READS)) {
  file.remove(NONCHIM_READS)
}

run_vsearch(
  c(
    "--search_exact", shQuote(POOLED),
    "--db", shQuote(NONCHIM_UNIQUE),
    "--strand", "plus",
    "--matched", shQuote(NONCHIM_READS),
    "--threads", as.character(THREADS),
    "--qmask", "none",
    "--dbmask", "none",
    "--fasta_width", "80"
  ),
  "01_EXACT_NONCHIMERIC_READ_EXTRACTION"
)

if (
  !file.exists(NONCHIM_READS) ||
  file.info(NONCHIM_READS)$size < 500 * 1024^2
) {
  stop("Exact non-chimeric matched FASTA is unexpectedly small/empty.")
}

# ------------------------------------------------------------
# 2. Map ONLY exact non-chimeric reads to existing OTU97 centroids
# ------------------------------------------------------------

if (file.exists(CORRECTED_OTU_TABLE)) {
  file.remove(CORRECTED_OTU_TABLE)
}

run_vsearch(
  c(
    "--usearch_global", shQuote(NONCHIM_READS),
    "--db", shQuote(OTUS),
    "--id", "0.97",
    "--strand", "plus",
    "--otutabout", shQuote(CORRECTED_OTU_TABLE),
    "--threads", as.character(THREADS),
    "--maxaccepts", "1",
    "--maxrejects", "64",
    "--qmask", "none",
    "--dbmask", "none"
  ),
  "02_MAP_EXACT_NONCHIMERIC_READS_TO_OTU97"
)

if (
  !file.exists(CORRECTED_OTU_TABLE) ||
  file.info(CORRECTED_OTU_TABLE)$size < 100
) {
  stop("Corrected OTU table is empty.")
}

# ------------------------------------------------------------
# 3. Read corrected OTU table
# ------------------------------------------------------------

otu <- read_tsv(
  CORRECTED_OTU_TABLE,
  show_col_types = FALSE,
  name_repair = "minimal"
)

names(otu)[1] <- "otu_id"

sample_cols <- setdiff(
  names(otu),
  c("otu_id", "taxonomy")
)

otu_mat <- as.matrix(
  otu[, sample_cols, drop = FALSE]
)

storage.mode(otu_mat) <- "numeric"
rownames(otu_mat) <- otu$otu_id

# OTU x sample -> sample x OTU
count_mat <- t(otu_mat)

corrected_depth <- rowSums(count_mat)

depth_df <- tibble(
  Run_ID = rownames(count_mat),
  corrected_nonchimera_OTU97_reads =
    as.numeric(corrected_depth),
  keep_ge4000 =
    corrected_nonchimera_OTU97_reads >= 4000,
  keep_ge2000 =
    corrected_nonchimera_OTU97_reads >= 2000
)

write_csv(
  depth_df,
  file.path(OUT, "01_CORRECTED_NONCHIMERIC_SAMPLE_DEPTH.csv")
)

# ------------------------------------------------------------
# 4. Reconcile totals against B1B chimera report
# ------------------------------------------------------------

expected_clean <- 3684641
expected_chimera_abundance <- 207074
expected_nonchimera <- 3477567

total_corrected <- sum(corrected_depth)

reconciliation <- tibble(
  B1B_clean_reads = expected_clean,
  B1B_abundance_weighted_chimera_reads =
    expected_chimera_abundance,
  expected_nonchimera_reads =
    expected_nonchimera,
  corrected_OTU_assigned_reads =
    total_corrected,
  difference_from_expected_nonchimera =
    total_corrected - expected_nonchimera,
  corrected_fraction_of_clean =
    total_corrected / expected_clean,
  expected_nonchimera_fraction =
    expected_nonchimera / expected_clean,
  passes_reconciliation =
    abs(total_corrected - expected_nonchimera) /
      expected_nonchimera < 0.005
)

write_csv(
  reconciliation,
  file.path(OUT, "02_CHIMERA_COUNT_RECONCILIATION.csv")
)

if (!reconciliation$passes_reconciliation) {
  stop(
    paste0(
      "Corrected read total does not reconcile with UCHIME abundance report. ",
      "Expected about ",
      expected_nonchimera,
      ", got ",
      total_corrected,
      "."
    )
  )
}

# ------------------------------------------------------------
# 5. Update metadata/depth/source balance
# ------------------------------------------------------------

meta_old <- read_csv(
  META_FILE,
  show_col_types = FALSE
)

# Drop B1B depth columns and replace them.
meta <- meta_old %>%
  select(
    -any_of(
      c(
        "assigned_OTU97_reads",
        "keep_ge2000",
        "keep_ge4000"
      )
    )
  ) %>%
  inner_join(
    depth_df,
    by = "Run_ID"
  )

write_csv(
  meta,
  file.path(
    DATA_DIR,
    "CRA002354_analysis_metadata_NONCHIMERIC_CORRECTED.csv"
  )
)

source_depth <- meta %>%
  group_by(pulmonary_binary) %>%
  summarise(
    samples_total = n(),
    patients_total = n_distinct(patient_id),
    samples_ge4000 = sum(keep_ge4000),
    patients_ge4000 =
      n_distinct(patient_id[keep_ge4000]),
    samples_ge2000 = sum(keep_ge2000),
    patients_ge2000 =
      n_distinct(patient_id[keep_ge2000]),
    min_corrected_reads =
      min(corrected_nonchimera_OTU97_reads),
    median_corrected_reads =
      median(corrected_nonchimera_OTU97_reads),
    .groups = "drop"
  )

write_csv(
  source_depth,
  file.path(OUT, "03_CORRECTED_SOURCE_STRATIFIED_DEPTH_QC.csv")
)

patient_long_4000 <- meta %>%
  filter(keep_ge4000) %>%
  group_by(patient_id, pulmonary_binary) %>%
  summarise(
    n_samples = n(),
    first_day = min(time_day),
    last_day = max(time_day),
    .groups = "drop"
  ) %>%
  group_by(pulmonary_binary) %>%
  summarise(
    patients = n(),
    patients_GE2 = sum(n_samples >= 2),
    patients_GE3 = sum(n_samples >= 3),
    .groups = "drop"
  )

write_csv(
  patient_long_4000,
  file.path(
    OUT,
    "04_CORRECTED_PRIMARY4000_LONGITUDINAL_DEPTH_BY_SOURCE.csv"
  )
)

# ------------------------------------------------------------
# 6. Taxonomy abundance-weighted coverage audit
# ------------------------------------------------------------

tax <- read_csv(
  TAX_FILE,
  show_col_types = FALSE
)

otu_totals <- tibble(
  otu_id = colnames(count_mat),
  abundance = as.numeric(colSums(count_mat))
)

tax_ab <- otu_totals %>%
  left_join(
    tax,
    by = "otu_id"
  )

ranks <- intersect(
  c("Kingdom","Phylum","Class","Order","Family","Genus"),
  names(tax_ab)
)

tax_coverage <- bind_rows(
  lapply(ranks, function(rk) {
    assigned <- !is.na(tax_ab[[rk]]) &
      tax_ab[[rk]] != ""

    tibble(
      rank = rk,
      OTUs_assigned = sum(assigned),
      OTUs_total = nrow(tax_ab),
      unweighted_OTU_fraction =
        mean(assigned),
      reads_in_assigned_OTUs =
        sum(tax_ab$abundance[assigned]),
      total_reads =
        sum(tax_ab$abundance),
      abundance_weighted_fraction =
        sum(tax_ab$abundance[assigned]) /
        sum(tax_ab$abundance)
    )
  })
)

write_csv(
  tax_coverage,
  file.path(OUT, "05_TAXONOMY_COVERAGE_WEIGHTED_BY_ABUNDANCE.csv")
)

# ------------------------------------------------------------
# 7. Corrected genus table (exploratory taxonomy only)
# ------------------------------------------------------------

make_genus_label <- function(row) {

  ranks0 <- c(
    "Kingdom","Phylum","Class",
    "Order","Family","Genus"
  )

  ranks0 <- ranks0[ranks0 %in% names(row)]

  vals <- as.character(row[ranks0])
  names(vals) <- ranks0

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

  lr <- names(nonmissing)[length(nonmissing)]
  lv <- nonmissing[length(nonmissing)]

  paste0(
    "Unclassified_",
    lr,
    "_",
    gsub("[^A-Za-z0-9]+", "_", lv)
  )
}

tax$genus_label <- apply(
  tax,
  1,
  make_genus_label
)

genus_for_otu <- tax$genus_label[
  match(
    colnames(count_mat),
    tax$otu_id
  )
]

genus_for_otu[
  is.na(genus_for_otu) |
  genus_for_otu == ""
] <- "Unclassified_Bacteria"

# rowsum works on rows; transpose OTU dimension into rows first.
genus_counts_t <- rowsum(
  t(count_mat),
  group = genus_for_otu,
  reorder = FALSE
)

genus_counts <- t(genus_counts_t)

genus_rel <- genus_counts /
  rowSums(genus_counts)

write_csv(
  as.data.frame(
    genus_counts,
    check.names = FALSE
  ) %>%
    rownames_to_column("Run_ID"),
  file.path(
    DATA_DIR,
    "CRA002354_genus_counts_NONCHIMERIC_CORRECTED.csv"
  )
)

write_csv(
  as.data.frame(
    genus_rel,
    check.names = FALSE
  ) %>%
    rownames_to_column("Run_ID"),
  file.path(
    DATA_DIR,
    "CRA002354_genus_relative_abundance_NONCHIMERIC_CORRECTED.csv"
  )
)

# ------------------------------------------------------------
# 8. Corrected OTU counts + relative abundance
# ------------------------------------------------------------

otu_count_df <- as.data.frame(
  count_mat,
  check.names = FALSE
) %>%
  rownames_to_column("Run_ID")

otu_rel <- count_mat / rowSums(count_mat)

otu_rel_df <- as.data.frame(
  otu_rel,
  check.names = FALSE
) %>%
  rownames_to_column("Run_ID")

write_csv(
  otu_count_df,
  file.path(
    DATA_DIR,
    "CRA002354_OTU97_counts_NONCHIMERIC_CORRECTED.csv"
  )
)

write_csv(
  otu_rel_df,
  file.path(
    DATA_DIR,
    "CRA002354_OTU97_relative_abundance_NONCHIMERIC_CORRECTED.csv"
  )
)

# ------------------------------------------------------------
# 9. Primary 4000-read rarefaction at OTU level
# ------------------------------------------------------------

primary_ids <- meta$Run_ID[
  meta$keep_ge4000
]

M4000 <- count_mat[
  primary_ids,
  ,
  drop = FALSE
]

set.seed(20260824)

rare4000 <- vegan::rrarefy(
  M4000,
  sample = 4000
)

alpha_from_counts <- function(M) {

  depth <- rowSums(M)

  rel <- M / depth

  observed <- rowSums(M > 0)

  shannon <- apply(rel, 1, function(p) {
    p <- p[p > 0]
    -sum(p * log(p))
  })

  simpson <- apply(rel, 1, function(p) {
    1 - sum(p^2)
  })

  tibble(
    Run_ID = rownames(M),
    rarefied_depth = depth,
    Observed_OTU97 = observed,
    Shannon_OTU97 = shannon,
    Simpson_OTU97 = simpson
  )
}

alpha4000 <- alpha_from_counts(
  rare4000
) %>%
  left_join(
    meta,
    by = "Run_ID"
  )

write_csv(
  as.data.frame(
    rare4000,
    check.names = FALSE
  ) %>%
    rownames_to_column("Run_ID"),
  file.path(
    DATA_DIR,
    "CRA002354_OTU97_rarefied_PRIMARY_4000.csv"
  )
)

write_csv(
  alpha4000,
  file.path(
    DATA_DIR,
    "CRA002354_OTU97_alpha_PRIMARY_rarefied4000.csv"
  )
)

write_csv(
  meta %>%
    filter(keep_ge4000),
  file.path(
    DATA_DIR,
    "CRA002354_analysis_metadata_PRIMARY_NONCHIMERIC_min4000.csv"
  )
)

# ------------------------------------------------------------
# 10. 2000-read V2 sensitivity rarefaction
# ------------------------------------------------------------

sens_ids <- meta$Run_ID[
  meta$keep_ge2000
]

M2000 <- count_mat[
  sens_ids,
  ,
  drop = FALSE
]

set.seed(20260824)

rare2000 <- vegan::rrarefy(
  M2000,
  sample = 2000
)

alpha2000 <- alpha_from_counts(
  rare2000
) %>%
  left_join(
    meta,
    by = "Run_ID"
  )

write_csv(
  as.data.frame(
    rare2000,
    check.names = FALSE
  ) %>%
    rownames_to_column("Run_ID"),
  file.path(
    DATA_DIR,
    "CRA002354_OTU97_rarefied_SENSITIVITY_2000.csv"
  )
)

write_csv(
  alpha2000,
  file.path(
    DATA_DIR,
    "CRA002354_OTU97_alpha_SENSITIVITY_rarefied2000.csv"
  )
)

write_csv(
  meta %>%
    filter(keep_ge2000),
  file.path(
    DATA_DIR,
    "CRA002354_analysis_metadata_SENSITIVITY_NONCHIMERIC_min2000.csv"
  )
)

# ------------------------------------------------------------
# 11. Final readiness
# ------------------------------------------------------------

get0 <- function(df, group, col) {
  x <- df[
    df$pulmonary_binary == group,
    col,
    drop = TRUE
  ]
  if (length(x) == 0) 0 else x[1]
}

final_ready <- tibble(
  corrected_total_nonchimera_reads =
    total_corrected,
  corrected_samples_ge4000 =
    sum(meta$keep_ge4000),
  corrected_samples_ge2000 =
    sum(meta$keep_ge2000),
  corrected_patients_ge4000 =
    n_distinct(meta$patient_id[meta$keep_ge4000]),
  pulmonary_patients_ge4000 =
    get0(
      source_depth,
      "PULMONARY",
      "patients_ge4000"
    ),
  nonpulmonary_patients_ge4000 =
    get0(
      source_depth,
      "NONPULMONARY_RECORDED",
      "patients_ge4000"
    ),
  pulmonary_GE2_ge4000 =
    get0(
      patient_long_4000,
      "PULMONARY",
      "patients_GE2"
    ),
  nonpulmonary_GE2_ge4000 =
    get0(
      patient_long_4000,
      "NONPULMONARY_RECORDED",
      "patients_GE2"
    ),
  genus_abundance_weighted_assignment =
    tax_coverage$abundance_weighted_fraction[
      tax_coverage$rank == "Genus"
    ][1],
  ready_for_corrected_step94B2 =
    (
      reconciliation$passes_reconciliation &&
      sum(meta$keep_ge4000) >= 110 &&
      n_distinct(meta$patient_id[meta$keep_ge4000]) >= 55 &&
      get0(source_depth, "PULMONARY", "patients_ge4000") >= 30 &&
      get0(source_depth, "NONPULMONARY_RECORDED", "patients_ge4000") >= 14 &&
      get0(patient_long_4000, "PULMONARY", "patients_GE2") >= 20 &&
      get0(patient_long_4000, "NONPULMONARY_RECORDED", "patients_GE2") >= 10
    )
)

write_csv(
  final_ready,
  file.path(
    OUT,
    "06_STEP94B1D_CORRECTED_ANALYSIS_READINESS.csv"
  )
)

writeLines(
  c(
    "STEP94B1D NON-CHIMERIC OTU TABLE + RAREFACTION FIX",
    "",
    paste0(
      "Corrected non-chimeric reads: ",
      total_corrected
    ),
    paste0(
      "Expected from UCHIME abundance report: ",
      expected_nonchimera
    ),
    paste0(
      "Samples >=4000 after true chimera exclusion: ",
      sum(meta$keep_ge4000)
    ),
    paste0(
      "Samples >=2000 after true chimera exclusion: ",
      sum(meta$keep_ge2000)
    ),
    paste0(
      "Abundance-weighted genus assignment fraction: ",
      signif(
        final_ready$genus_abundance_weighted_assignment,
        4
      )
    ),
    paste0(
      "Ready for corrected Step94B2: ",
      final_ready$ready_for_corrected_step94B2
    ),
    "",
    "Primary alpha metrics are now OTU97-level metrics after deterministic rarefaction to 4000 reads.",
    "Genus abundance is retained for exploratory taxonomic analyses, not for the primary richness metric."
  ),
  file.path(
    OUT,
    "07_STEP94B1D_INTERPRETATION.txt"
  )
)

writeLines(
  c(
    paste0("Completed: ", Sys.time()),
    paste0(
      "Ready for corrected Step94B2: ",
      final_ready$ready_for_corrected_step94B2
    ),
    "STEP94B1D COMPLETE"
  ),
  file.path(
    OUT,
    "_STEP94B1D_COMPLETE.ok"
  )
)

cat("STEP94B1D COMPLETE\n")
