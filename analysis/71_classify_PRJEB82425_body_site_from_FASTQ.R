# ============================================================
# Sepsis V2 - Step 71
# Classify PRJEB82425 rectal vs tracheal/ETA using raw FASTQ
# R 4.4.0 / Windows
#
# Study design:
#   ETA/tracheal aspirates: V1-V2 16S region
#   Rectal swabs:          V3-V4 16S region
#
# This script:
# 1) reads the Step70 patient-time map
# 2) finds paired FASTQ.gz recursively under PRJEB82425
# 3) scans a small number of reads only (no full FASTQ processing)
# 4) detects characteristic primer-region sequence motifs
# 5) classifies each Run conservatively:
#       rectal
#       tracheal_ETA
#       unresolved
# 6) creates the gut-only metadata map
#
# IMPORTANT:
# Classification is conservative. If primer motifs were already
# trimmed from raw reads, the run remains unresolved rather than
# being guessed.
# ============================================================

options(stringsAsFactors = FALSE)

pkgs <- c(
  "readr",
  "dplyr",
  "stringr",
  "tibble",
  "purrr"
)

missing <- pkgs[
  !vapply(
    pkgs,
    requireNamespace,
    logical(1),
    quietly = TRUE
  )
]

if (length(missing) > 0) {

  message(
    "Installing missing packages: ",
    paste(missing, collapse = ", ")
  )

  install.packages(
    missing,
    repos = "https://cloud.r-project.org"
  )
}

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(stringr)
  library(tibble)
  library(purrr)
})


# ------------------------------------------------------------
# 1. Paths
# ------------------------------------------------------------
PROJECT <- "PRJEB82425"

PROJECT_ROOT <- file.path(
  "E:/sepsis_project/data",
  PROJECT
)

STEP70_ROOT <- paste0(
  "E:/sepsis_project/results/",
  "V2_10_exact_metadata_mapping"
)

MAP_FILE <- file.path(
  STEP70_ROOT,
  "PRJEB82425_explicit_sample_patient_time_map.csv"
)

OUT_ROOT <- paste0(
  "E:/sepsis_project/results/",
  "V2_11_PRJEB82425_body_site"
)

dir.create(
  OUT_ROOT,
  recursive = TRUE,
  showWarnings = FALSE
)


# ------------------------------------------------------------
# 2. Parameters
# ------------------------------------------------------------
# Number of reads sampled from each FASTQ mate.
N_READS <- 5000

# Only inspect the beginning of each read.
# Primers/adapters should occur near the 5' end.
PREFIX_NT <- 80

# Minimum proportion of inspected reads carrying
# region-specific evidence.
MIN_SIGNAL <- 0.02

# Winning region needs to be this much stronger than the other.
MIN_RATIO <- 3


# ------------------------------------------------------------
# 3. Sequence helpers
# ------------------------------------------------------------
clean_seq <- function(x) {

  toupper(
    gsub(
      "[^ACGTN]",
      "",
      x
    )
  )
}


revcomp <- function(x) {

  x <- clean_seq(x)

  paste0(
    rev(
      chartr(
        "ACGT",
        "TGCA",
        strsplit(
          x,
          "",
          fixed = TRUE
        )[[1]]
      )
    ),
    collapse = ""
  )
}


# Stable sequence fragments from commonly used primer sets
# for these 16S regions.
#
# V1-V2:
#   forward region around 27F: AGAGTTTGAT...
#   reverse region around 338R: GCTGCCTCC...
#
# V3-V4:
#   forward region around 341F: CCTACGGG...
#   reverse region around 785R/805R: GACTAC...GGGTATCTAATCC
#
# We use short, stable fragments rather than requiring a
# full degenerate-primer exact sequence.
V12_BASE_MOTIFS <- c(
  "AGAGTTTGAT",
  "GCTGCCTCC",
  "CCCGTAGGAG"
)

V34_BASE_MOTIFS <- c(
  "CCTACGGG",
  "GGACTAC",
  "GACTAC",
  "GGGTATCTAATCC"
)


expand_with_revcomp <- function(motifs) {

  unique(
    c(
      motifs,
      vapply(
        motifs,
        revcomp,
        character(1)
      )
    )
  )
}


V12_MOTIFS <- expand_with_revcomp(
  V12_BASE_MOTIFS
)

V34_MOTIFS <- expand_with_revcomp(
  V34_BASE_MOTIFS
)


# ------------------------------------------------------------
# 4. Read first N FASTQ sequences
# ------------------------------------------------------------
read_fastq_sequences <- function(
  path,
  n_reads = N_READS
) {

  if (
    is.na(path) ||
    !file.exists(path)
  ) {
    return(character())
  }

  con <- gzfile(
    path,
    open = "rt"
  )

  on.exit(
    close(con),
    add = TRUE
  )


  # FASTQ = 4 lines/read
  lines <- tryCatch(
    readLines(
      con,
      n = n_reads * 4,
      warn = FALSE
    ),
    error = function(e) character()
  )


  if (length(lines) < 4) {
    return(character())
  }


  seq_idx <- seq(
    from = 2,
    to = length(lines),
    by = 4
  )

  seq_idx <- seq_idx[
    seq_idx <= length(lines)
  ]


  seqs <- lines[
    seq_idx
  ]

  clean_seq(
    seqs
  )
}


prefixes <- function(
  seqs,
  n = PREFIX_NT
) {

  if (length(seqs) == 0) {
    return(character())
  }

  substr(
    seqs,
    1,
    n
  )
}


count_motif_reads <- function(
  seqs,
  motifs
) {

  if (length(seqs) == 0) {
    return(0L)
  }


  hit <- rep(
    FALSE,
    length(seqs)
  )


  for (m in motifs) {

    hit <- hit |
      grepl(
        m,
        seqs,
        fixed = TRUE
      )
  }


  sum(hit)
}


scan_fastq <- function(path) {

  seqs <- read_fastq_sequences(
    path
  )

  seqs <- prefixes(
    seqs
  )

  n <- length(seqs)


  if (n == 0) {

    return(
      tibble(
        Reads_Scanned = 0,
        V12_Hits = 0,
        V34_Hits = 0,
        V12_Prop = NA_real_,
        V34_Prop = NA_real_
      )
    )
  }


  v12 <- count_motif_reads(
    seqs,
    V12_MOTIFS
  )

  v34 <- count_motif_reads(
    seqs,
    V34_MOTIFS
  )


  tibble(
    Reads_Scanned = n,
    V12_Hits = v12,
    V34_Hits = v34,
    V12_Prop = v12 / n,
    V34_Prop = v34 / n
  )
}


# ------------------------------------------------------------
# 5. Load metadata
# ------------------------------------------------------------
if (!file.exists(MAP_FILE)) {

  stop(
    paste0(
      "Cannot find Step70 map: ",
      MAP_FILE
    )
  )
}


meta <- suppressMessages(
  read_csv(
    MAP_FILE,
    show_col_types = FALSE
  )
)


if (!"Run_ID" %in% names(meta)) {

  stop(
    "Run_ID column missing from PRJEB82425 map."
  )
}


# ------------------------------------------------------------
# 6. Locate FASTQ files recursively
# ------------------------------------------------------------
message(
  "Searching FASTQ files under: ",
  PROJECT_ROOT
)


fastqs <- list.files(
  PROJECT_ROOT,
  pattern = "\\.fastq\\.gz$",
  recursive = TRUE,
  full.names = TRUE,
  ignore.case = TRUE
)


message(
  "FASTQ.gz files found: ",
  length(fastqs)
)


fastq_index <- tibble(
  FASTQ_Path = fastqs,
  FASTQ_Name = basename(fastqs)
) |>
  mutate(
    Run_ID = str_extract(
      FASTQ_Name,
      "(ERR|SRR|DRR|CRR)[0-9]+"
    ),

    Mate = case_when(

      str_detect(
        FASTQ_Name,
        "_1\\.fastq\\.gz$"
      ) ~ "R1",

      str_detect(
        FASTQ_Name,
        "_2\\.fastq\\.gz$"
      ) ~ "R2",

      TRUE ~ "single_or_unknown"
    )
  ) |>
  filter(
    !is.na(Run_ID)
  )


write_excel_csv(
  fastq_index,
  file.path(
    OUT_ROOT,
    "PRJEB82425_FASTQ_inventory.csv"
  ),
  na = ""
)


# ------------------------------------------------------------
# 7. Scan each Run
# ------------------------------------------------------------
runs <- sort(
  unique(
    meta$Run_ID[
      !is.na(
        meta$Run_ID
      )
    ]
  )
)


results <- vector(
  "list",
  length(runs)
)


for (i in seq_along(runs)) {

  run <- runs[i]


  f <- fastq_index |>
    filter(
      Run_ID == run
    )


  r1 <- f$FASTQ_Path[
    f$Mate == "R1"
  ]

  r2 <- f$FASTQ_Path[
    f$Mate == "R2"
  ]


  if (length(r1) == 0) {
    r1 <- NA_character_
  } else {
    r1 <- r1[1]
  }


  if (length(r2) == 0) {
    r2 <- NA_character_
  } else {
    r2 <- r2[1]
  }


  s1 <- scan_fastq(
    r1
  )

  s2 <- scan_fastq(
    r2
  )


  total_reads <- sum(
    c(
      s1$Reads_Scanned,
      s2$Reads_Scanned
    ),
    na.rm = TRUE
  )


  total_v12 <- sum(
    c(
      s1$V12_Hits,
      s2$V12_Hits
    ),
    na.rm = TRUE
  )


  total_v34 <- sum(
    c(
      s1$V34_Hits,
      s2$V34_Hits
    ),
    na.rm = TRUE
  )


  prop_v12 <- if (
    total_reads > 0
  ) {
    total_v12 / total_reads
  } else {
    NA_real_
  }


  prop_v34 <- if (
    total_reads > 0
  ) {
    total_v34 / total_reads
  } else {
    NA_real_
  }


  ratio_v12 <- (
    prop_v12 + 1e-8
  ) / (
    prop_v34 + 1e-8
  )


  ratio_v34 <- (
    prop_v34 + 1e-8
  ) / (
    prop_v12 + 1e-8
  )


  body_site <- case_when(

    !is.na(prop_v34) &&
      prop_v34 >= MIN_SIGNAL &&
      ratio_v34 >= MIN_RATIO ~
      "rectal",

    !is.na(prop_v12) &&
      prop_v12 >= MIN_SIGNAL &&
      ratio_v12 >= MIN_RATIO ~
      "tracheal_ETA",

    TRUE ~
      "unresolved"
  )


  confidence <- case_when(

    body_site == "rectal" ~
      ratio_v34,

    body_site == "tracheal_ETA" ~
      ratio_v12,

    TRUE ~
      pmax(
        ratio_v12,
        ratio_v34,
        na.rm = TRUE
      )
  )


  results[[i]] <- tibble(
    Run_ID = run,

    R1_Path = r1,
    R2_Path = r2,

    R1_Reads_Scanned =
      s1$Reads_Scanned,

    R2_Reads_Scanned =
      s2$Reads_Scanned,

    V12_Hits = total_v12,
    V34_Hits = total_v34,

    V12_Prop = prop_v12,
    V34_Prop = prop_v34,

    Classification_Ratio =
      confidence,

    Body_Site_FASTQ =
      body_site,

    Classification_Method =
      paste(
        "FASTQ 5-prime motif scan;",
        "V1-V2 => ETA;",
        "V3-V4 => rectal"
      )
  )


  if (
    i %% 10 == 0 ||
    i == length(runs)
  ) {

    message(
      "Scanned ",
      i,
      "/",
      length(runs),
      " runs"
    )
  }
}


classification <- bind_rows(
  results
)


write_excel_csv(
  classification,
  file.path(
    OUT_ROOT,
    "PRJEB82425_FASTQ_body_site_classification.csv"
  ),
  na = ""
)


# ------------------------------------------------------------
# 8. Join back to metadata
# ------------------------------------------------------------
final_map <- meta |>
  left_join(
    classification |>
      select(
        Run_ID,
        Body_Site_FASTQ,
        V12_Prop,
        V34_Prop,
        Classification_Ratio,
        Classification_Method
      ),
    by = "Run_ID"
  ) |>
  mutate(
    Body_Site = case_when(

      Body_Site_FASTQ == "rectal" ~
        "rectal",

      Body_Site_FASTQ == "tracheal_ETA" ~
        "tracheal_ETA",

      TRUE ~
        NA_character_
    )
  )


write_excel_csv(
  final_map,
  file.path(
    OUT_ROOT,
    "PRJEB82425_sample_patient_time_BODY_SITE_map.csv"
  ),
  na = ""
)


gut_map <- final_map |>
  filter(
    Body_Site == "rectal"
  )


write_excel_csv(
  gut_map,
  file.path(
    OUT_ROOT,
    "PRJEB82425_GUT_ONLY_final_map.csv"
  ),
  na = ""
)


eta_map <- final_map |>
  filter(
    Body_Site == "tracheal_ETA"
  )


write_excel_csv(
  eta_map,
  file.path(
    OUT_ROOT,
    "PRJEB82425_ETA_ONLY_map.csv"
  ),
  na = ""
)


unresolved_map <- final_map |>
  filter(
    is.na(Body_Site)
  )


write_excel_csv(
  unresolved_map,
  file.path(
    OUT_ROOT,
    "PRJEB82425_BODY_SITE_unresolved.csv"
  ),
  na = ""
)


# ------------------------------------------------------------
# 9. QC
# ------------------------------------------------------------
qc <- tibble(
  Total_Metadata_Runs =
    n_distinct(
      meta$Run_ID
    ),

  FASTQ_Runs_Found =
    n_distinct(
      fastq_index$Run_ID
    ),

  Rectal_Runs =
    n_distinct(
      gut_map$Run_ID
    ),

  ETA_Runs =
    n_distinct(
      eta_map$Run_ID
    ),

  Unresolved_Runs =
    n_distinct(
      unresolved_map$Run_ID
    ),

  Rectal_Patients =
    n_distinct(
      gut_map$Patient_ID,
      na.rm = TRUE
    ),

  ETA_Patients =
    n_distinct(
      eta_map$Patient_ID,
      na.rm = TRUE
    )
)


write_excel_csv(
  qc,
  file.path(
    OUT_ROOT,
    "PRJEB82425_body_site_QC.csv"
  ),
  na = ""
)


# ------------------------------------------------------------
# 10. Longitudinal gut QC
# ------------------------------------------------------------
gut_long <- gut_map |>
  filter(
    !is.na(Patient_ID),
    !is.na(Time_Raw)
  ) |>
  distinct(
    Patient_ID,
    Time_Raw
  ) |>
  count(
    Patient_ID,
    name = "N_Timepoints"
  )


gut_long_qc <- tibble(
  Gut_Runs =
    n_distinct(
      gut_map$Run_ID
    ),

  Gut_Samples =
    n_distinct(
      gut_map$Sample_ID
    ),

  Gut_Patients =
    n_distinct(
      gut_map$Patient_ID,
      na.rm = TRUE
    ),

  Gut_Patients_GE2_Timepoints =
    sum(
      gut_long$N_Timepoints >= 2
    ),

  Gut_Patients_GE3_Timepoints =
    sum(
      gut_long$N_Timepoints >= 3
    )
)


write_excel_csv(
  gut_long_qc,
  file.path(
    OUT_ROOT,
    "PRJEB82425_GUT_longitudinal_QC.csv"
  ),
  na = ""
)


# ------------------------------------------------------------
# 11. Console
# ------------------------------------------------------------
cat(
  "\n============================================================\n"
)

cat(
  "SEPSIS V2 - STEP 71 COMPLETE\n"
)

cat(
  "R version: ",
  R.version.string,
  "\n",
  sep = ""
)

cat(
  "============================================================\n\n"
)


cat(
  "BODY SITE QC:\n"
)

print(
  qc,
  width = Inf
)


cat(
  "\nGUT LONGITUDINAL QC:\n"
)

print(
  gut_long_qc,
  width = Inf
)


cat(
  "\nExpected from publication:",
  "93 rectal swabs + 73 ETA samples = 166 samples total.\n"
)

cat(
  "The ENA metadata currently contains 165 mapped runs,",
  "so a one-sample difference from the publication total may remain",
  "and should be documented rather than guessed.\n"
)


cat(
  "\nOutput folder:\n",
  OUT_ROOT,
  "\n",
  sep = ""
)

cat(
  "============================================================\n"
)
