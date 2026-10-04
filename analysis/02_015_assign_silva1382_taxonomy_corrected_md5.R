# ============================================================
# Sepsis V2 - Step 85A
# Assign genus-level taxonomy to all CURRENT FINAL selected ASVs
# using SILVA 138.2 SSU NR99 DADA2 training set
#
# CURRENT TAXONOMY-READY COHORTS (8):
#   Existing validated ASV:
#     PRJNA430161
#     PRJNA691455
#     PRJNA978257
#
#   Step84 paired-end DADA2 closed:
#     PRJEB82425
#     PRJNA516701
#     PRJNA578267
#     PRJNA851469
#     PRJNA1166732
#
# IMPORTANT
#   - Cohort ASV tables remain SEPARATE.
#   - We do NOT merge ASV identities across different 16S regions.
#   - Taxonomy is assigned with the SAME SILVA version and settings.
#   - Species assignment is NOT performed in Step85A.
#   - Output includes bootstrap confidence matrices.
#   - Strict input root only:
#       E:/sepsis_project/data/_V2_ANALYSIS_READY
#
# Reference:
#   SILVA 138.2 SSU NR99 DADA2-formatted toGenus training set.
#   minBoot = 80
#   tryRC = TRUE
#   multithread = FALSE (Windows-safe)
#
# R 4.4.0 / Windows
# ============================================================

options(stringsAsFactors = FALSE)

# ------------------------------------------------------------
# Packages
# ------------------------------------------------------------
cran_pkgs <- c("readr","dplyr","stringr","purrr","tibble")
missing_cran <- cran_pkgs[
  !vapply(cran_pkgs, requireNamespace, logical(1), quietly=TRUE)
]

if (length(missing_cran)) {
  install.packages(
    missing_cran,
    repos="https://cloud.r-project.org"
  )
}

if (!requireNamespace("dada2", quietly=TRUE)) {
  if (!requireNamespace("BiocManager", quietly=TRUE)) {
    install.packages(
      "BiocManager",
      repos="https://cloud.r-project.org"
    )
  }

  BiocManager::install(
    "dada2",
    ask=FALSE,
    update=FALSE
  )
}

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(stringr)
  library(purrr)
  library(tibble)
  library(dada2)
})

# ------------------------------------------------------------
# Paths
# ------------------------------------------------------------
ROOT <- "E:/sepsis_project"

CANON <- file.path(
  ROOT,
  "data",
  "_V2_ANALYSIS_READY"
)

PROJECT_ROOT <- file.path(
  CANON,
  "01_PROJECTS"
)

REFERENCE_ROOT <- file.path(
  CANON,
  "00_REFERENCE",
  "SILVA_138.2_DADA2_SSU"
)

OUT <- file.path(
  ROOT,
  "results",
  "V2_25A2_SILVA1382_TAXONOMY"
)

dir.create(
  REFERENCE_ROOT,
  recursive=TRUE,
  showWarnings=FALSE
)

dir.create(
  OUT,
  recursive=TRUE,
  showWarnings=FALSE
)

LOG <- file.path(
  OUT,
  "_STEP85A2_runtime_checkpoints.txt"
)

ERR <- file.path(
  OUT,
  "_STEP85A2_FATAL_ERROR.txt"
)

if (file.exists(ERR)) unlink(ERR)

ck <- function(x) {
  cat(
    paste0(x,": ",Sys.time(),"\n"),
    file=LOG,
    append=TRUE
  )
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

md5_one <- function(p) {
  if (!file.exists(p)) return(NA_character_)
  unname(tools::md5sum(p))
}

valid_rds_matrix <- function(p) {

  if (!file.exists(p)) return(FALSE)

  tryCatch(
    {
      x <- readRDS(p)

      is.matrix(x) &&
        nrow(x)>0 &&
        ncol(x)>0
    },
    error=function(e) FALSE
  )
}

# ------------------------------------------------------------
# Project registry
# ------------------------------------------------------------
registry <- tribble(
  ~project, ~expected_samples, ~source_stage, ~preferred_pattern,

  "PRJNA430161",
  26L,
  "STEP83A_EXISTING_ASV_REUSE",
  "seqtab_final",

  "PRJNA691455",
  49L,
  "STEP83A_EXISTING_ASV_REUSE_RERUN_V2_SELECTED",
  "seqtab_final",

  "PRJNA978257",
  13L,
  "STEP83A_EXISTING_ASV_REUSE",
  "seqtab_final",

  "PRJEB82425",
  92L,
  "STEP84B_FULL_DADA2",
  "seqtab_final_step84B",

  "PRJNA516701",
  170L,
  "STEP84C_QC_REVIEWED",
  "seqtab_final_step84C_QC_REVIEWED",

  "PRJNA578267",
  207L,
  "STEP84B4B_MINPARENT8",
  "seqtab_final_step84B4B_minParent8",

  "PRJNA851469",
  119L,
  "STEP84B_FULL_DADA2",
  "seqtab_final_step84B",

  "PRJNA1166732",
  120L,
  "STEP84B_FULL_DADA2_MINOVERLAP8",
  "seqtab_final_step84B"
)

# ------------------------------------------------------------
# SILVA official reference files
# ------------------------------------------------------------
BASE_URL <- paste0(
  "https://www.arb-silva.de/",
  "current-release/DADA2/1.36.0/SSU/"
)

# Direct download target behind the official SILVA file listing.
# Using this for forced re-download avoids HTML/redirect responses being
# mistaken for small checksum files.
FILEADMIN_BASE_URL <- paste0(
  "https://www.arb-silva.de/fileadmin/silva_databases/",
  "current/DADA2/1.36.0/SSU/"
)

TRAIN_FILE <- "silva_nr99_v138.2_toGenus_trainset.fa.gz"
TRAIN_MD5  <- paste0(TRAIN_FILE,".md5")
TRAIN_DOI  <- paste0(TRAIN_FILE,".doi")

TRAIN_PATH <- file.path(
  REFERENCE_ROOT,
  TRAIN_FILE
)

MD5_PATH <- file.path(
  REFERENCE_ROOT,
  TRAIN_MD5
)

DOI_PATH <- file.path(
  REFERENCE_ROOT,
  TRAIN_DOI
)

download_reference <- function(filename) {

  dest <- file.path(
    REFERENCE_ROOT,
    filename
  )

  if (
    file.exists(dest) &&
    file.info(dest)$size>0
  ) {
    message(
      "Reference file already exists: ",
      filename
    )
    return(dest)
  }

  url <- paste0(
    BASE_URL,
    filename
  )

  message(
    "Downloading official SILVA file:\n",
    url
  )

  download.file(
    url=url,
    destfile=dest,
    mode="wb",
    method="libcurl",
    quiet=FALSE
  )

  if (
    !file.exists(dest) ||
    file.info(dest)$size==0
  ) {
    stop(
      "Download failed or empty file: ",
      filename
    )
  }

  dest
}

extract_md5_from_file <- function(path) {

  if (!file.exists(path) || file.info(path)$size<=0) {
    return(NA_character_)
  }

  # Read as raw bytes rather than readLines().  This is deliberately
  # tolerant of BOMs, unusual line endings and checksum-file formatting.
  n <- as.integer(file.info(path)$size)
  con <- file(path, open="rb")
  on.exit(close(con), add=TRUE)
  raw <- readBin(con, what="raw", n=n)

  txt <- tryCatch(
    rawToChar(raw),
    error=function(e) paste(as.character(raw), collapse=" ")
  )

  # Base-R PCRE extraction, independent of stringr/regex masking.
  hit <- regmatches(
    txt,
    regexpr(
      "[0-9A-Fa-f]{32}",
      txt,
      perl=TRUE
    )
  )

  if (!length(hit) || identical(hit, character(0)) || nchar(hit)!=32) {
    return(NA_character_)
  }

  tolower(hit)
}

force_redownload_official <- function(filename) {

  dest <- file.path(
    REFERENCE_ROOT,
    filename
  )

  if (file.exists(dest)) {
    backup <- paste0(
      dest,
      ".replaced_",
      format(Sys.time(), "%Y%m%d_%H%M%S")
    )
    file.rename(dest, backup)
  }

  url <- paste0(
    FILEADMIN_BASE_URL,
    filename
  )

  message(
    "Force-downloading official SILVA file:\n",
    url
  )

  download.file(
    url=url,
    destfile=dest,
    mode="wb",
    method="libcurl",
    quiet=FALSE
  )

  if (!file.exists(dest) || file.info(dest)$size<=0) {
    stop("Forced SILVA download failed: ", filename)
  }

  dest
}

verify_reference_md5 <- function() {

  if (!file.exists(TRAIN_PATH)) {
    stop("Training FASTA is missing: ", TRAIN_PATH)
  }

  # The first 85A attempt may have downloaded a malformed/redirect checksum
  # response.  If the local checksum file contains no 32-hex digest, replace
  # ONLY that tiny checksum file from the resolved official file URL.
  expected <- extract_md5_from_file(MD5_PATH)

  if (is.na(expected)) {

    message(
      "Local SILVA .md5 file contains no parsable 32-hex checksum; ",
      "re-downloading the official checksum file."
    )

    force_redownload_official(TRAIN_MD5)
    expected <- extract_md5_from_file(MD5_PATH)
  }

  if (is.na(expected)) {
    # Persist a diagnostic copy of the unparseable checksum response.
    if (file.exists(MD5_PATH)) {
      file.copy(
        MD5_PATH,
        file.path(
          OUT,
          "SILVA_MD5_UNPARSEABLE_RESPONSE.bin"
        ),
        overwrite=TRUE
      )
    }
    stop(
      "Official SILVA checksum file still contains no parsable 32-hex MD5."
    )
  }

  observed <- tolower(md5_one(TRAIN_PATH))

  # If the checksum is valid but does not match the training file, replace
  # the training file once from the resolved official download URL and retry.
  if (is.na(observed) || expected!=observed) {

    message(
      "SILVA MD5 mismatch. Re-downloading the training FASTA once from ",
      "the resolved official SILVA file URL."
    )

    force_redownload_official(TRAIN_FILE)
    observed <- tolower(md5_one(TRAIN_PATH))
  }

  if (is.na(observed) || expected!=observed) {
    stop(
      paste0(
        "SILVA training-set MD5 verification failed after one clean retry.\n",
        "Expected: ", expected, "\n",
        "Observed: ", observed
      )
    )
  }

  # Light gzip/readability guard before taxonomy assignment.
  con <- gzfile(TRAIN_PATH, open="rt")
  first_lines <- tryCatch(
    readLines(con, n=4, warn=FALSE),
    error=function(e) character(0),
    finally=close(con)
  )

  if (!length(first_lines) || !startsWith(first_lines[1], ">")) {
    stop(
      "SILVA training FASTA passed MD5 but failed gzip/FASTA readability guard."
    )
  }

  tibble(
    file=TRAIN_FILE,
    expected_md5=expected,
    observed_md5=observed,
    md5_match=TRUE,
    gzip_fasta_readable=TRUE,
    verification_route="ROBUST_RAW_MD5_PARSE_AND_OFFICIAL_FILEADMIN_RETRY"
  )
}

# ------------------------------------------------------------
# Select final ASV in a project's 02_asv_selected directory
# ------------------------------------------------------------
choose_selected_asv <- function(
  project,
  preferred_pattern
) {

  d <- file.path(
    PROJECT_ROOT,
    project,
    "02_asv_selected"
  )

  if (!dir.exists(d)) {
    stop(
      project,
      ": 02_asv_selected directory does not exist."
    )
  }

  files <- list.files(
    d,
    full.names=TRUE,
    recursive=FALSE,
    pattern="\\.rds$",
    ignore.case=TRUE
  )

  files <- files[
    vapply(
      files,
      valid_rds_matrix,
      logical(1)
    )
  ]

  if (!length(files)) {
    stop(
      project,
      ": no valid matrix RDS in 02_asv_selected."
    )
  }

  preferred <- files[
    str_detect(
      basename(files),
      regex(
        preferred_pattern,
        ignore_case=TRUE
      )
    )
  ]

  candidates <- if (
    length(preferred)
  ) {
    preferred
  } else {
    files
  }

  # If several candidates have identical content, choosing any
  # is safe. If content differs, stop instead of guessing.
  md5s <- vapply(
    candidates,
    md5_one,
    character(1)
  )

  unique_md5 <- unique(
    md5s[
      !is.na(md5s)
    ]
  )

  if (
    length(candidates)>1 &&
    length(unique_md5)>1
  ) {
    stop(
      paste0(
        project,
        ": multiple non-identical selected ASV RDS candidates remain:\n",
        paste(
          basename(candidates),
          collapse="\n"
        )
      )
    )
  }

  candidates[1]
}

# ------------------------------------------------------------
# Taxonomic completeness helper
# ------------------------------------------------------------
taxonomy_qc <- function(
  project,
  taxa
) {

  ranks <- colnames(
    taxa
  )

  map_dfr(
    ranks,
    function(rank) {

      x <- taxa[,rank]

      tibble(
        project=project,
        rank=rank,
        total_ASVs=nrow(taxa),
        assigned_ASVs=sum(
          !is.na(x) &
          x!=""
        ),
        assigned_pct=round(
          100*mean(
            !is.na(x) &
            x!=""
          ),
          2
        )
      )
    }
  )
}

# ------------------------------------------------------------
# Main
# ------------------------------------------------------------
main <- function() {

  ck("STEP85A2 STARTED")

  # ==========================================================
  # 1. Download + verify SILVA
  # ==========================================================

  download_reference(
    TRAIN_FILE
  )

  download_reference(
    TRAIN_MD5
  )

  # DOI is useful for provenance but not required for execution.
  try(
    download_reference(
      TRAIN_DOI
    ),
    silent=TRUE
  )

  md5_qc <- verify_reference_md5()

  write_excel_csv(
    md5_qc,
    file.path(
      OUT,
      "V2_STEP85A_SILVA_reference_MD5_QC.csv"
    ),
    na=""
  )

  reference_registry <- tibble(
    database="SILVA",
    release="138.2",
    molecule="SSU",
    reference_subset="NR99",
    dada2_reference_file=TRAIN_FILE,
    local_path=TRAIN_PATH,
    md5=md5_one(TRAIN_PATH),
    assignment_method="DADA2 assignTaxonomy naive Bayesian classifier",
    minBoot=80,
    tryRC=TRUE,
    multithread=FALSE,
    taxonomic_target="Kingdom_to_Genus",
    species_assignment="NOT_RUN_IN_STEP85A"
  )

  write_excel_csv(
    reference_registry,
    file.path(
      OUT,
      "V2_STEP85A_reference_registry.csv"
    ),
    na=""
  )

  ck("SILVA 138.2 REFERENCE READY AND MD5 VERIFIED")

  # ==========================================================
  # 2. Resolve all eight final selected ASV tables
  # ==========================================================

  selected_rows <- list()

  for (
    i in seq_len(
      nrow(registry)
    )
  ) {

    project <- registry$project[i]

    p <- choose_selected_asv(
      project,
      registry$preferred_pattern[i]
    )

    x <- readRDS(p)

    if (
      nrow(x)!=
      registry$expected_samples[i]
    ) {
      stop(
        paste0(
          project,
          ": sample-count guard failed. Expected ",
          registry$expected_samples[i],
          "; observed ",
          nrow(x)
        )
      )
    }

    if (
      is.null(colnames(x)) ||
      anyDuplicated(colnames(x))
    ) {
      stop(
        project,
        ": ASV sequence columns are missing or duplicated."
      )
    }

    if (
      any(rowSums(x)==0)
    ) {
      stop(
        project,
        ": selected ASV table contains zero-read samples."
      )
    }

    selected_rows[[
      length(selected_rows)+1
    ]] <- tibble(
      project=project,
      source_stage=
        registry$source_stage[i],
      expected_samples=
        registry$expected_samples[i],
      observed_samples=
        nrow(x),
      ASVs=
        ncol(x),
      total_reads=
        sum(x),
      selected_asv_path=
        p,
      selected_asv_md5=
        md5_one(p)
    )
  }

  selected_registry <- bind_rows(
    selected_rows
  )

  if (
    nrow(selected_registry)!=8 ||
    sum(
      selected_registry$observed_samples
    )!=796
  ) {
    stop(
      "Eight-cohort selected-ASV registry guard failed."
    )
  }

  write_excel_csv(
    selected_registry,
    file.path(
      OUT,
      "V2_STEP85A_selected_ASV_registry.csv"
    ),
    na=""
  )

  ck("EIGHT SELECTED ASV TABLES GUARDED")

  # ==========================================================
  # 3. Assign taxonomy per cohort independently
  # ==========================================================

  tax_qc_rows <- list()
  cohort_summary_rows <- list()

  for (
    i in seq_len(
      nrow(selected_registry)
    )
  ) {

    project <- selected_registry$project[i]

    ck(
      paste0(
        project,
        " TAXONOMY STARTED"
      )
    )

    seqtab <- readRDS(
      selected_registry$selected_asv_path[i]
    )

    seqs <- colnames(
      seqtab
    )

    # --------------------------------------------------------
    # Taxonomy, with bootstrap confidence output
    # --------------------------------------------------------
    res <- assignTaxonomy(
      seqs,
      TRAIN_PATH,
      minBoot=80,
      tryRC=TRUE,
      outputBootstraps=TRUE,
      multithread=FALSE,
      verbose=TRUE
    )

    taxa <- res$tax
    boot <- res$boot

    if (
      !is.matrix(taxa) ||
      nrow(taxa)!=length(seqs)
    ) {
      stop(
        project,
        ": taxonomy output dimension mismatch."
      )
    }

    if (
      !identical(
        rownames(taxa),
        seqs
      )
    ) {
      stop(
        project,
        ": taxonomy row sequences do not match ASV sequences."
      )
    }

    # Cohort-local stable ASV IDs.
    asv_ids <- sprintf(
      "%s_ASV%06d",
      project,
      seq_along(seqs)
    )

    tax_df <- as.data.frame(
      taxa,
      stringsAsFactors=FALSE
    ) |>
      mutate(
        project=project,
        ASV_ID=asv_ids,
        ASV_sequence=seqs,
        total_abundance=as.numeric(
          colSums(
            seqtab[
              ,
              seqs,
              drop=FALSE
            ]
          )
        ),
        .before=1
      )

    boot_df <- as.data.frame(
      boot,
      stringsAsFactors=FALSE
    )

    names(boot_df) <- paste0(
      names(boot_df),
      "_bootstrap"
    )

    boot_df <- boot_df |>
      mutate(
        project=project,
        ASV_ID=asv_ids,
        ASV_sequence=seqs,
        .before=1
      )

    tax_dir <- file.path(
      PROJECT_ROOT,
      project,
      "04_taxonomy_qc"
    )

    dir.create(
      tax_dir,
      recursive=TRUE,
      showWarnings=FALSE
    )

    tax_csv <- file.path(
      tax_dir,
      paste0(
        project,
        "_taxonomy_SILVA1382_minBoot80.csv"
      )
    )

    boot_csv <- file.path(
      tax_dir,
      paste0(
        project,
        "_taxonomy_bootstrap_SILVA1382.csv"
      )
    )

    write_excel_csv(
      tax_df,
      tax_csv,
      na=""
    )

    write_excel_csv(
      boot_df,
      boot_csv,
      na=""
    )

    saveRDS(
      taxa,
      file.path(
        tax_dir,
        paste0(
          project,
          "_taxonomy_SILVA1382_minBoot80.rds"
        )
      )
    )

    saveRDS(
      boot,
      file.path(
        tax_dir,
        paste0(
          project,
          "_taxonomy_bootstrap_SILVA1382.rds"
        )
      )
    )

    # ASV ID map keeps sequence identity explicit.
    write_excel_csv(
      tibble(
        project=project,
        ASV_ID=asv_ids,
        ASV_sequence=seqs,
        total_abundance=as.numeric(
          colSums(seqtab)
        )
      ),
      file.path(
        tax_dir,
        paste0(
          project,
          "_ASV_ID_sequence_map.csv"
        )
      ),
      na=""
    )

    tq <- taxonomy_qc(
      project,
      taxa
    )

    tax_qc_rows[[
      length(tax_qc_rows)+1
    ]] <- tq

    kingdom <- if (
      "Kingdom" %in%
      colnames(taxa)
    ) {
      taxa[,"Kingdom"]
    } else {
      rep(
        NA_character_,
        nrow(taxa)
      )
    }

    genus <- if (
      "Genus" %in%
      colnames(taxa)
    ) {
      taxa[,"Genus"]
    } else {
      rep(
        NA_character_,
        nrow(taxa)
      )
    }

    cohort_summary_rows[[
      length(cohort_summary_rows)+1
    ]] <- tibble(
      project=project,
      samples=nrow(seqtab),
      ASVs=ncol(seqtab),
      total_reads=sum(seqtab),
      Bacteria_ASVs=sum(
        kingdom=="Bacteria",
        na.rm=TRUE
      ),
      Archaea_ASVs=sum(
        kingdom=="Archaea",
        na.rm=TRUE
      ),
      Kingdom_unassigned_ASVs=sum(
        is.na(kingdom)
      ),
      Genus_assigned_ASVs=sum(
        !is.na(genus)
      ),
      Genus_assigned_pct=round(
        100*mean(
          !is.na(genus)
        ),
        2
      ),
      taxonomy_csv=tax_csv,
      bootstrap_csv=boot_csv
    )

    ck(
      paste0(
        project,
        " TAXONOMY COMPLETE"
      )
    )

    rm(
      seqtab,
      res,
      taxa,
      boot,
      tax_df,
      boot_df
    )

    gc(
      verbose=FALSE
    )
  }

  tax_qc <- bind_rows(
    tax_qc_rows
  )

  cohort_summary <- bind_rows(
    cohort_summary_rows
  )

  write_excel_csv(
    tax_qc,
    file.path(
      OUT,
      "V2_STEP85A_taxonomy_rank_completeness.csv"
    ),
    na=""
  )

  write_excel_csv(
    cohort_summary,
    file.path(
      OUT,
      "V2_STEP85A_cohort_taxonomy_summary.csv"
    ),
    na=""
  )

  # ==========================================================
  # 4. Global closeout
  # ==========================================================

  global <- tibble(
    metric=c(
      "Taxonomy_ready_projects",
      "Selected_samples",
      "Total_cohort_specific_ASVs",
      "Projects_taxonomy_completed",
      "Reference_MD5_verified",
      "minBoot",
      "tryRC"
    ),
    value=c(
      8,
      sum(
        selected_registry$observed_samples
      ),
      sum(
        selected_registry$ASVs
      ),
      nrow(
        cohort_summary
      ),
      1,
      80,
      1
    )
  )

  write_excel_csv(
    global,
    file.path(
      OUT,
      "V2_STEP85A_global_summary.csv"
    ),
    na=""
  )

  # ==========================================================
  # README
  # ==========================================================

  readme <- c(
    "SEPSIS V2 - STEP85A SILVA 138.2 TAXONOMY",
    paste0(
      "Created: ",
      Sys.time()
    ),
    "",
    "REFERENCE",
    "SILVA 138.2 SSU NR99 DADA2-formatted toGenus training set.",
    paste0(
      "Local reference: ",
      TRAIN_PATH
    ),
    "",
    "SETTINGS",
    "assignTaxonomy",
    "minBoot = 80",
    "tryRC = TRUE",
    "outputBootstraps = TRUE",
    "multithread = FALSE",
    "",
    "WHY tryRC=TRUE",
    "The selected ASV tables originate from multiple published cohorts and some reused ASV tables may differ in sequence orientation. tryRC allows DADA2 to compare the reverse complement and use the better-supported classification.",
    "",
    "IMPORTANT",
    "Cohort ASV tables remain separate.",
    "ASV IDs are cohort-specific and must not be merged across different 16S variable regions.",
    "Species-level assignment is not performed in Step85A.",
    "",
    "CURRENT TAXONOMY-READY PROJECTS",
    paste(
      registry$project,
      collapse=", "
    ),
    "",
    "NEXT",
    "Step85B/86 should review taxonomy completeness and remove obvious non-target features (e.g. non-prokaryotic / organelle classifications where applicable), then freeze one global sample-depth QC rule before alpha/beta diversity and longitudinal modeling."
  )

  writeLines(
    readme,
    file.path(
      OUT,
      "README_STEP85A_SILVA1382_TAXONOMY.txt"
    ),
    useBytes=TRUE
  )

  writeLines(
    c(
      paste0(
        "Completed: ",
        Sys.time()
      ),
      "Status: STEP85A2 COMPLETE",
      "Reference: SILVA 138.2 SSU NR99",
      "minBoot: 80",
      "tryRC: TRUE",
      "Eight selected ASV cohorts assigned independently."
    ),
    file.path(
      OUT,
      "_STEP85A_COMPLETE.ok"
    )
  )

  ck("STEP85A2 COMPLETE")

  cat("\n============================================================\n")
  cat("SEPSIS V2 - STEP85A2 COMPLETE\n")
  cat("============================================================\n\n")

  print(
    global,
    n=Inf,
    width=Inf
  )

  cat("\nCOHORT TAXONOMY SUMMARY:\n")

  print(
    cohort_summary,
    n=Inf,
    width=Inf
  )

  cat("\nOutput:\n")
  cat(
    OUT,
    "\n"
  )
}

tryCatch(
  main(),
  error=function(e) {

    msg <- c(
      paste0(
        "STEP85A FATAL ERROR: ",
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

    ck(
      "STEP85A2 FAILED"
    )

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
