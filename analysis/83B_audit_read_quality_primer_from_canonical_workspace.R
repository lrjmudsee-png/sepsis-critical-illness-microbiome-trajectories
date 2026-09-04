# ============================================================
# Sepsis V2 - Step 83B
# Read length / quality / primer audit from CANONICAL workspace
#
# STRICT INPUT ROOT:
#   E:/sepsis_project/data/_V2_ANALYSIS_READY
#
# NO broad filesystem scan.
# NO DADA2 is run in this step.
#
# Purpose:
#   - sample a small number of canonical raw files per cohort
#   - determine read length distributions
#   - inspect FASTQ quality profiles
#   - detect common 16S primer motifs if still present
#   - distinguish paired / single / FASTA special routes
#   - generate evidence for cohort-specific DADA2 parameters
#
# R 4.4.0 / Windows
# ============================================================

options(stringsAsFactors = FALSE)

pkgs <- c(
  "readr","dplyr","tidyr",
  "stringr","purrr","tibble"
)

missing <- pkgs[
  !vapply(pkgs, requireNamespace, logical(1), quietly=TRUE)
]

if (length(missing)) {
  install.packages(
    missing,
    repos="https://cloud.r-project.org"
  )
}

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(purrr)
  library(tibble)
})

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

FREEZE_ROOT <- file.path(
  CANON,
  "00_FREEZE"
)

OUT <- file.path(
  ROOT,
  "results",
  "V2_23B_READ_QUALITY_PRIMER_AUDIT"
)

dir.create(
  OUT,
  recursive=TRUE,
  showWarnings=FALSE
)

LOG <- file.path(
  OUT,
  "_STEP83B_runtime_checkpoints.txt"
)

ERR <- file.path(
  OUT,
  "_STEP83B_FATAL_ERROR.txt"
)

if (file.exists(ERR)) unlink(ERR)

ck <- function(x) {
  cat(
    paste0(x,": ",Sys.time(),"\n"),
    file=LOG,
    append=TRUE
  )
}

clean_chr <- function(x) {
  x <- trimws(as.character(x))
  x[
    is.na(x) |
    x=="" |
    tolower(x) %in% c(
      "na","nan","n/a","null","none"
    )
  ] <- NA_character_
  x
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

# ============================================================
# File helpers
# ============================================================

sequence_file_type <- function(path) {
  x <- basename(path)

  case_when(
    str_detect(
      x,
      regex("\\.(fastq|fq)\\.gz$",ignore_case=TRUE)
    ) ~ "FASTQ_GZ",

    str_detect(
      x,
      regex("\\.(fastq|fq)$",ignore_case=TRUE)
    ) ~ "FASTQ",

    str_detect(
      x,
      regex("\\.(fa|fasta)\\.gz$",ignore_case=TRUE)
    ) ~ "FASTA_GZ",

    str_detect(
      x,
      regex("\\.(fa|fasta)$",ignore_case=TRUE)
    ) ~ "FASTA",

    TRUE ~ "OTHER"
  )
}

detect_mate <- function(path) {

  x <- basename(path)

  case_when(
    str_detect(
      x,
      regex(
        "(_R?1_|_R?1\\.|_1\\.(fastq|fq)|_1\\.(fastq|fq)\\.gz$|\\.1\\.(fastq|fq))",
        ignore_case=TRUE
      )
    ) ~ "R1",

    str_detect(
      x,
      regex(
        "(_R?2_|_R?2\\.|_2\\.(fastq|fq)|_2\\.(fastq|fq)\\.gz$|\\.2\\.(fastq|fq))",
        ignore_case=TRUE
      )
    ) ~ "R2",

    str_detect(
      x,
      regex("\\.(fa|fasta)(\\.gz)?$",ignore_case=TRUE)
    ) ~ "FASTA",

    TRUE ~ "SINGLE_OR_UNSPECIFIED"
  )
}

open_text_connection <- function(path) {
  if (
    str_detect(
      path,
      regex("\\.gz$",ignore_case=TRUE)
    )
  ) {
    gzfile(path,"rt")
  } else {
    file(path,"rt")
  }
}

pick_evenly <- function(x,n=3) {

  x <- sort(unique(x))

  if (length(x)<=n) return(x)

  idx <- unique(
    round(
      seq(
        1,
        length(x),
        length.out=n
      )
    )
  )

  x[idx]
}

# ============================================================
# FASTQ parser
# ============================================================

read_fastq_sample <- function(
  path,
  max_reads=1000
) {

  con <- open_text_connection(path)
  on.exit(close(con),add=TRUE)

  z <- readLines(
    con,
    n=max_reads*4,
    warn=FALSE
  )

  usable <- floor(length(z)/4)*4

  if (usable<4) {
    return(
      list(
        seq=character(),
        qual=character()
      )
    )
  }

  z <- z[seq_len(usable)]

  list(
    seq=z[seq(2,usable,by=4)],
    qual=z[seq(4,usable,by=4)]
  )
}

# ============================================================
# FASTA parser
# ============================================================

read_fasta_sample <- function(
  path,
  max_reads=1000
) {

  con <- open_text_connection(path)
  on.exit(close(con),add=TRUE)

  seqs <- character()
  current <- character()

  while (
    length(seqs)<max_reads
  ) {

    line <- readLines(
      con,
      n=1,
      warn=FALSE
    )

    if (!length(line)) {
      if (length(current)) {
        seqs <- c(
          seqs,
          paste0(current,collapse="")
        )
      }
      break
    }

    if (startsWith(line,">")) {

      if (length(current)) {
        seqs <- c(
          seqs,
          paste0(current,collapse="")
        )
        current <- character()

        if (length(seqs)>=max_reads) break
      }

    } else {
      current <- c(
        current,
        trimws(line)
      )
    }
  }

  seqs
}

# ============================================================
# Quality calculations
# ============================================================

quality_matrix <- function(qstrings) {

  if (!length(qstrings)) {
    return(
      matrix(
        numeric(),
        nrow=0,
        ncol=0
      )
    )
  }

  lens <- nchar(qstrings)

  max_len <- max(lens)

  mat <- matrix(
    NA_real_,
    nrow=length(qstrings),
    ncol=max_len
  )

  for (i in seq_along(qstrings)) {

    q <- as.integer(
      charToRaw(
        qstrings[i]
      )
    ) - 33L

    mat[
      i,
      seq_along(q)
    ] <- q
  }

  mat
}

quality_cycle_summary <- function(
  qmat,
  project,
  file_path,
  mate
) {

  if (!nrow(qmat) || !ncol(qmat)) {
    return(tibble())
  }

  tibble(
    project=project,
    file_path=file_path,
    mate=mate,
    cycle=seq_len(ncol(qmat)),

    n_reads=apply(
      qmat,
      2,
      function(x) sum(!is.na(x))
    ),

    q25=apply(
      qmat,
      2,
      function(x) {
        if (all(is.na(x))) return(NA_real_)
        quantile(
          x,
          0.25,
          na.rm=TRUE,
          names=FALSE
        )
      }
    ),

    q50=apply(
      qmat,
      2,
      function(x) {
        if (all(is.na(x))) return(NA_real_)
        median(
          x,
          na.rm=TRUE
        )
      }
    ),

    q75=apply(
      qmat,
      2,
      function(x) {
        if (all(is.na(x))) return(NA_real_)
        quantile(
          x,
          0.75,
          na.rm=TRUE,
          names=FALSE
        )
      }
    ),

    pct_q30=apply(
      qmat,
      2,
      function(x) {
        if (all(is.na(x))) return(NA_real_)
        100*mean(
          x>=30,
          na.rm=TRUE
        )
      }
    )
  )
}

summarise_fastq_file <- function(
  path,
  project,
  mate,
  max_reads=1000
) {

  x <- read_fastq_sample(
    path,
    max_reads=max_reads
  )

  seqs <- x$seq
  quals <- x$qual

  if (!length(seqs)) {

    return(
      list(
        file_summary=tibble(
          project=project,
          file_path=path,
          file_name=basename(path),
          mate=mate,
          file_type=sequence_file_type(path),
          sampled_reads=0L,
          read_length_min=NA_real_,
          read_length_q10=NA_real_,
          read_length_median=NA_real_,
          read_length_q90=NA_real_,
          read_length_max=NA_real_,
          median_read_mean_Q=NA_real_,
          median_last20_Q=NA_real_,
          overall_pct_bases_Q30=NA_real_,
          min_quality_ascii=NA_integer_,
          max_quality_ascii=NA_integer_,
          parse_status="NO_READS_PARSED"
        ),
        cycles=tibble(),
        seqs=character()
      )
    )
  }

  lens <- nchar(seqs)

  qmat <- quality_matrix(quals)

  read_mean_q <- rowMeans(
    qmat,
    na.rm=TRUE
  )

  last20_q <- rep(
    NA_real_,
    nrow(qmat)
  )

  for (i in seq_len(nrow(qmat))) {

    idx <- which(
      !is.na(qmat[i,])
    )

    if (!length(idx)) next

    tail_idx <- tail(
      idx,
      min(
        20,
        length(idx)
      )
    )

    last20_q[i] <- mean(
      qmat[
        i,
        tail_idx
      ],
      na.rm=TRUE
    )
  }

  all_q_chars <- unlist(
    lapply(
      quals,
      function(q) as.integer(charToRaw(q))
    ),
    use.names=FALSE
  )

  file_summary <- tibble(
    project=project,
    file_path=path,
    file_name=basename(path),
    mate=mate,
    file_type=sequence_file_type(path),
    sampled_reads=length(seqs),

    read_length_min=min(lens),

    read_length_q10=as.numeric(
      quantile(
        lens,
        0.10,
        names=FALSE
      )
    ),

    read_length_median=median(lens),

    read_length_q90=as.numeric(
      quantile(
        lens,
        0.90,
        names=FALSE
      )
    ),

    read_length_max=max(lens),

    median_read_mean_Q=round(
      median(
        read_mean_q,
        na.rm=TRUE
      ),
      2
    ),

    median_last20_Q=round(
      median(
        last20_q,
        na.rm=TRUE
      ),
      2
    ),

    overall_pct_bases_Q30=round(
      100*mean(
        qmat>=30,
        na.rm=TRUE
      ),
      2
    ),

    min_quality_ascii=min(
      all_q_chars
    ),

    max_quality_ascii=max(
      all_q_chars
    ),

    parse_status="OK"
  )

  cycles <- quality_cycle_summary(
    qmat,
    project,
    path,
    mate
  )

  list(
    file_summary=file_summary,
    cycles=cycles,
    seqs=seqs
  )
}

summarise_fasta_file <- function(
  path,
  project,
  max_reads=1000
) {

  seqs <- read_fasta_sample(
    path,
    max_reads=max_reads
  )

  if (!length(seqs)) {

    return(
      list(
        file_summary=tibble(
          project=project,
          file_path=path,
          file_name=basename(path),
          mate="FASTA",
          file_type=sequence_file_type(path),
          sampled_reads=0L,
          read_length_min=NA_real_,
          read_length_q10=NA_real_,
          read_length_median=NA_real_,
          read_length_q90=NA_real_,
          read_length_max=NA_real_,
          median_read_mean_Q=NA_real_,
          median_last20_Q=NA_real_,
          overall_pct_bases_Q30=NA_real_,
          min_quality_ascii=NA_integer_,
          max_quality_ascii=NA_integer_,
          parse_status="NO_SEQUENCES_PARSED"
        ),
        cycles=tibble(),
        seqs=character()
      )
    )
  }

  lens <- nchar(seqs)

  list(
    file_summary=tibble(
      project=project,
      file_path=path,
      file_name=basename(path),
      mate="FASTA",
      file_type=sequence_file_type(path),
      sampled_reads=length(seqs),

      read_length_min=min(lens),

      read_length_q10=as.numeric(
        quantile(
          lens,
          0.10,
          names=FALSE
        )
      ),

      read_length_median=median(lens),

      read_length_q90=as.numeric(
        quantile(
          lens,
          0.90,
          names=FALSE
        )
      ),

      read_length_max=max(lens),

      median_read_mean_Q=NA_real_,
      median_last20_Q=NA_real_,
      overall_pct_bases_Q30=NA_real_,
      min_quality_ascii=NA_integer_,
      max_quality_ascii=NA_integer_,
      parse_status="OK_FASTA_NO_QUALITY"
    ),

    cycles=tibble(),
    seqs=seqs
  )
}

# ============================================================
# IUPAC primer detection
# ============================================================

iupac_to_regex <- function(seq) {

  map <- c(
    A="A",
    C="C",
    G="G",
    T="T",
    R="[AG]",
    Y="[CT]",
    S="[GC]",
    W="[AT]",
    K="[GT]",
    M="[AC]",
    B="[CGT]",
    D="[AGT]",
    H="[ACT]",
    V="[ACG]",
    N="[ACGT]"
  )

  x <- strsplit(
    toupper(seq),
    ""
  )[[1]]

  paste0(
    unname(
      map[x]
    ),
    collapse=""
  )
}

primer_catalog <- tribble(
  ~primer_name, ~sequence, ~region_hint,
  "27F",    "AGAGTTTGATCMTGGCTCAG",     "V1-start/full-length",
  "338F",   "ACTCCTACGGGAGGCAGCAG",     "V3",
  "341F",   "CCTACGGGNGGCWGCAG",        "V3-V4",
  "515F",   "GTGYCAGCMGCCGCGGTAA",      "V4",
  "515F_v2","GTGCCAGCMGCCGCGGTAA",      "V4",
  "799F",   "AACMGGATTAGATACCCKG",       "V5",
  "785R",   "GACTACHVGGGTATCTAATCC",     "V3-V4",
  "806R",   "GGACTACHVGGGTWTCTAAT",      "V3-V4/V4",
  "806R_v2","GGACTACNVGGGTWTCTAAT",      "V4",
  "1193R",  "ACGTCATCCCCACCTTCC",        "V5-V7",
  "1492R",  "TACGGYTACCTTGTTACGACTT",    "full-length"
)

detect_primers <- function(
  seqs,
  project,
  file_path,
  mate,
  search_prefix=45
) {

  if (!length(seqs)) return(tibble())

  prefixes <- substr(
    toupper(seqs),
    1,
    search_prefix
  )

  map_dfr(
    seq_len(nrow(primer_catalog)),
    function(i) {

      rx <- iupac_to_regex(
        primer_catalog$sequence[i]
      )

      hit <- str_detect(
        prefixes,
        regex(
          rx,
          ignore_case=TRUE
        )
      )

      tibble(
        project=project,
        file_path=file_path,
        mate=mate,
        primer_name=primer_catalog$primer_name[i],
        primer_sequence=primer_catalog$sequence[i],
        region_hint=primer_catalog$region_hint[i],
        sampled_reads=length(prefixes),
        reads_with_primer=sum(hit,na.rm=TRUE),
        primer_detection_pct=round(
          100*mean(hit,na.rm=TRUE),
          2
        )
      )
    }
  )
}

# ============================================================
# Main
# ============================================================

main <- function() {

  ck("STEP83B STARTED")

  route_path <- file.path(
    FREEZE_ROOT,
    "V2_STEP82C_FROZEN_sequence_processing_routes.csv"
  )

  if (!file.exists(route_path)) {
    stop(
      "Canonical processing-route file was not found in 00_FREEZE."
    )
  }

  routes <- safe_csv(route_path)

  if (nrow(routes)!=16) {
    stop(
      "Expected 16 primary 16S project routes."
    )
  }

  if (!dir.exists(PROJECT_ROOT)) {
    stop(
      "Canonical 01_PROJECTS directory does not exist."
    )
  }

  ck("CANONICAL INPUTS GUARDED")

  # ==========================================================
  # 1. Canonical raw inventory ONLY
  # ==========================================================

  inventory_rows <- list()

  for (
    project in routes$project
  ) {

    raw_dir <- file.path(
      PROJECT_ROOT,
      project,
      "01_raw"
    )

    files <- if (dir.exists(raw_dir)) {
      list.files(
        raw_dir,
        full.names=TRUE,
        recursive=FALSE
      )
    } else {
      character()
    }

    files <- files[
      str_detect(
        basename(files),
        regex(
          "\\.(fastq|fq|fa|fasta)(\\.gz)?$",
          ignore_case=TRUE
        )
      )
    ]

    if (length(files)) {

      inventory_rows[[
        length(inventory_rows)+1
      ]] <- tibble(
        project=project,
        file_path=files,
        file_name=basename(files),
        file_type=map_chr(
          files,
          sequence_file_type
        ),
        mate=map_chr(
          files,
          detect_mate
        ),
        bytes=as.numeric(
          file.info(files)$size
        )
      )
    }
  }

  inventory <- if (length(inventory_rows)) {
    bind_rows(inventory_rows)
  } else {
    tibble()
  }

  write_excel_csv(
    inventory,
    file.path(
      OUT,
      "V2_STEP83B_canonical_raw_inventory.csv"
    ),
    na=""
  )

  ck("CANONICAL RAW INVENTORY COMPLETE")

  # ==========================================================
  # 2. Select a small audit sample per cohort
  # ==========================================================

  sample_plan_rows <- list()

  for (
    project in routes$project
  ) {

    z <- inventory |>
      filter(
        project==!!project
      )

    if (!nrow(z)) next

    for (
      mate in unique(z$mate)
    ) {

      files <- z$file_path[
        z$mate==mate
      ]

      n_pick <- if (
        mate %in% c(
          "R1","R2"
        )
      ) {
        3
      } else {
        4
      }

      chosen <- pick_evenly(
        files,
        n=n_pick
      )

      if (length(chosen)) {

        sample_plan_rows[[
          length(sample_plan_rows)+1
        ]] <- tibble(
          project=project,
          mate=mate,
          file_path=chosen
        )
      }
    }
  }

  sample_plan <- if (
    length(sample_plan_rows)
  ) {
    bind_rows(sample_plan_rows)
  } else {
    tibble()
  }

  write_excel_csv(
    sample_plan,
    file.path(
      OUT,
      "V2_STEP83B_sampled_file_plan.csv"
    ),
    na=""
  )

  # ==========================================================
  # 3. Parse sampled sequence files
  # ==========================================================

  file_summaries <- list()
  cycle_summaries <- list()
  primer_summaries <- list()

  for (
    i in seq_len(
      nrow(sample_plan)
    )
  ) {

    project <- sample_plan$project[i]
    mate <- sample_plan$mate[i]
    path <- sample_plan$file_path[i]

    ftype <- sequence_file_type(path)

    if (
      ftype %in% c(
        "FASTQ_GZ",
        "FASTQ"
      )
    ) {

      z <- summarise_fastq_file(
        path,
        project,
        mate,
        max_reads=1000
      )

    } else if (
      ftype %in% c(
        "FASTA_GZ",
        "FASTA"
      )
    ) {

      z <- summarise_fasta_file(
        path,
        project,
        max_reads=1000
      )

    } else {
      next
    }

    file_summaries[[
      length(file_summaries)+1
    ]] <- z$file_summary

    if (
      nrow(z$cycles)
    ) {
      cycle_summaries[[
        length(cycle_summaries)+1
      ]] <- z$cycles
    }

    if (
      length(z$seqs)
    ) {

      primer_summaries[[
        length(primer_summaries)+1
      ]] <- detect_primers(
        z$seqs,
        project,
        path,
        mate
      )
    }
  }

  file_summary <- if (
    length(file_summaries)
  ) {
    bind_rows(
      file_summaries
    )
  } else {
    tibble()
  }

  cycle_summary <- if (
    length(cycle_summaries)
  ) {
    bind_rows(
      cycle_summaries
    )
  } else {
    tibble()
  }

  primer_summary <- if (
    length(primer_summaries)
  ) {
    bind_rows(
      primer_summaries
    )
  } else {
    tibble()
  }

  write_excel_csv(
    file_summary,
    file.path(
      OUT,
      "V2_STEP83B_sampled_file_read_quality_summary.csv"
    ),
    na=""
  )

  write_excel_csv(
    cycle_summary,
    file.path(
      OUT,
      "V2_STEP83B_quality_by_cycle.csv"
    ),
    na=""
  )

  write_excel_csv(
    primer_summary,
    file.path(
      OUT,
      "V2_STEP83B_primer_detection_all.csv"
    ),
    na=""
  )

  ck("READ QUALITY AND PRIMER AUDIT COMPLETE")

  # ==========================================================
  # 4. Best primer per project/mate
  # ==========================================================

  best_primer <- primer_summary |>
    group_by(
      project,
      mate
    ) |>
    arrange(
      desc(primer_detection_pct),
      primer_name
    ) |>
    slice(1) |>
    ungroup() |>
    select(
      project,
      mate,
      best_primer=primer_name,
      primer_region_hint=region_hint,
      best_primer_detection_pct=
        primer_detection_pct
    )

  write_excel_csv(
    best_primer,
    file.path(
      OUT,
      "V2_STEP83B_best_primer_by_project_mate.csv"
    ),
    na=""
  )

  # ==========================================================
  # 5. Project / mate summary
  # ==========================================================

  mate_summary <- file_summary |>
    group_by(
      project,
      mate,
      file_type
    ) |>
    summarise(
      sampled_files=n(),

      sampled_reads=sum(
        sampled_reads
      ),

      median_read_length=median(
        read_length_median,
        na.rm=TRUE
      ),

      min_sampled_q10_length=min(
        read_length_q10,
        na.rm=TRUE
      ),

      max_sampled_q90_length=max(
        read_length_q90,
        na.rm=TRUE
      ),

      median_read_mean_Q=if (
        all(
          is.na(
            median_read_mean_Q
          )
        )
      ) {
        NA_real_
      } else {
        median(
          median_read_mean_Q,
          na.rm=TRUE
        )
      },

      median_tail20_Q=if (
        all(
          is.na(
            median_last20_Q
          )
        )
      ) {
        NA_real_
      } else {
        median(
          median_last20_Q,
          na.rm=TRUE
        )
      },

      median_pct_bases_Q30=if (
        all(
          is.na(
            overall_pct_bases_Q30
          )
        )
      ) {
        NA_real_
      } else {
        median(
          overall_pct_bases_Q30,
          na.rm=TRUE
        )
      },

      .groups="drop"
    ) |>
    left_join(
      best_primer,
      by=c(
        "project",
        "mate"
      )
    )

  write_excel_csv(
    mate_summary,
    file.path(
      OUT,
      "V2_STEP83B_project_mate_summary.csv"
    ),
    na=""
  )

  # ==========================================================
  # 6. Create one project-level parameter-readiness table
  # ==========================================================

  project_raw <- inventory |>
    group_by(project) |>
    summarise(
      canonical_raw_files=n(),
      R1_files=sum(mate=="R1"),
      R2_files=sum(mate=="R2"),
      single_files=sum(
        mate=="SINGLE_OR_UNSPECIFIED"
      ),
      fasta_files=sum(
        mate=="FASTA"
      ),
      .groups="drop"
    )

  # Pivot mate summaries to wide.
  mate_wide <- mate_summary |>
    mutate(
      mate2=case_when(
        mate=="R1" ~ "R1",
        mate=="R2" ~ "R2",
        mate=="FASTA" ~ "FASTA",
        TRUE ~ "SINGLE"
      )
    ) |>
    select(
      project,
      mate2,
      median_read_length,
      median_read_mean_Q,
      median_tail20_Q,
      median_pct_bases_Q30,
      best_primer,
      primer_region_hint,
      best_primer_detection_pct
    ) |>
    pivot_wider(
      names_from=mate2,
      values_from=c(
        median_read_length,
        median_read_mean_Q,
        median_tail20_Q,
        median_pct_bases_Q30,
        best_primer,
        primer_region_hint,
        best_primer_detection_pct
      ),
      names_sep="_"
    )

  project_audit <- routes |>
    select(
      project,
      analysis_module,
      expected_runs,
      runs_missing_local_raw,
      frozen_processing_route,
      immediate_priority
    ) |>
    left_join(
      project_raw,
      by="project"
    ) |>
    left_join(
      mate_wide,
      by="project"
    ) |>
    mutate(
      canonical_raw_files=coalesce(
        canonical_raw_files,
        0L
      ),

      R1_files=coalesce(
        R1_files,
        0L
      ),

      R2_files=coalesce(
        R2_files,
        0L
      ),

      single_files=coalesce(
        single_files,
        0L
      ),

      fasta_files=coalesce(
        fasta_files,
        0L
      ),

      primer_state=case_when(

        coalesce(
          best_primer_detection_pct_R1,
          0
        )>=50 |
        coalesce(
          best_primer_detection_pct_R2,
          0
        )>=50 |
        coalesce(
          best_primer_detection_pct_SINGLE,
          0
        )>=50 |
        coalesce(
          best_primer_detection_pct_FASTA,
          0
        )>=50 ~
          "COMMON_16S_PRIMER_STRONGLY_DETECTED",

        coalesce(
          best_primer_detection_pct_R1,
          0
        )>=5 |
        coalesce(
          best_primer_detection_pct_R2,
          0
        )>=5 |
        coalesce(
          best_primer_detection_pct_SINGLE,
          0
        )>=5 |
        coalesce(
          best_primer_detection_pct_FASTA,
          0
        )>=5 ~
          "COMMON_16S_PRIMER_PARTIALLY_DETECTED",

        TRUE ~
          "NO_COMMON_PRIMER_AT_READ_START_OR_ALREADY_TRIMMED"
      ),

      parameter_audit_status=case_when(

        frozen_processing_route==
        "REUSE_EXISTING_ASV" ~
          "ASV_REUSE_NO_NEW_DADA2_REQUIRED",

        frozen_processing_route==
        "SELECT_BETWEEN_TWO_FULL_ASV_VERSIONS" ~
          "ASV_SELECTED_IN_STEP83A_NO_NEW_DADA2_REQUIRED",

        str_detect(
          frozen_processing_route,
          "STATIC_SUPPORT_DEFER"
        ) ~
          "DEFER_STATIC_SUPPORT",

        frozen_processing_route==
        "SPECIAL_FASTA_FEATURE_HARMONIZATION" ~
          "FASTA_SPECIAL_ROUTE_NO_STANDARD_FASTQ_DADA2",

        runs_missing_local_raw>0 ~
          "RAW_REPAIR_REQUIRED_BEFORE_FINAL_DADA2",

        canonical_raw_files>0 ~
          "READ_AUDIT_COMPLETE_READY_FOR_COHORT_PARAMETER_SELECTION",

        TRUE ~
          "NO_CANONICAL_RAW"
      )
    )

  write_excel_csv(
    project_audit,
    file.path(
      OUT,
      "V2_STEP83B_project_parameter_audit.csv"
    ),
    na=""
  )

  # ==========================================================
  # 7. Step83C queues
  # ==========================================================

  ready_dada2 <- project_audit |>
    filter(
      parameter_audit_status==
      "READ_AUDIT_COMPLETE_READY_FOR_COHORT_PARAMETER_SELECTION"
    )

  repair_first <- project_audit |>
    filter(
      parameter_audit_status==
      "RAW_REPAIR_REQUIRED_BEFORE_FINAL_DADA2"
    )

  special <- project_audit |>
    filter(
      parameter_audit_status==
      "FASTA_SPECIAL_ROUTE_NO_STANDARD_FASTQ_DADA2"
    )

  reuse <- project_audit |>
    filter(
      str_detect(
        parameter_audit_status,
        "ASV_"
      )
    )

  write_excel_csv(
    ready_dada2,
    file.path(
      OUT,
      "V2_STEP83B_queue_DADA2_PARAMETER_SELECTION.csv"
    ),
    na=""
  )

  write_excel_csv(
    repair_first,
    file.path(
      OUT,
      "V2_STEP83B_queue_RAW_REPAIR_FIRST.csv"
    ),
    na=""
  )

  write_excel_csv(
    special,
    file.path(
      OUT,
      "V2_STEP83B_queue_FASTA_SPECIAL.csv"
    ),
    na=""
  )

  write_excel_csv(
    reuse,
    file.path(
      OUT,
      "V2_STEP83B_queue_ASV_REUSE.csv"
    ),
    na=""
  )

  # ==========================================================
  # 8. Summary
  # ==========================================================

  summary <- tibble(
    metric=c(
      "Projects",
      "Canonical_sequence_files",
      "Sampled_files",
      "Projects_ready_for_parameter_selection",
      "Projects_raw_repair_required",
      "Projects_ASV_reuse_no_DADA2",
      "Projects_FASTA_special",
      "Projects_static_deferred",
      "Sampled_files_parse_failed"
    ),

    value=c(
      nrow(project_audit),
      nrow(inventory),
      nrow(sample_plan),

      sum(
        project_audit$parameter_audit_status==
        "READ_AUDIT_COMPLETE_READY_FOR_COHORT_PARAMETER_SELECTION"
      ),

      sum(
        project_audit$parameter_audit_status==
        "RAW_REPAIR_REQUIRED_BEFORE_FINAL_DADA2"
      ),

      sum(
        str_detect(
          project_audit$parameter_audit_status,
          "ASV_"
        )
      ),

      sum(
        project_audit$parameter_audit_status==
        "FASTA_SPECIAL_ROUTE_NO_STANDARD_FASTQ_DADA2"
      ),

      sum(
        project_audit$parameter_audit_status==
        "DEFER_STATIC_SUPPORT"
      ),

      sum(
        file_summary$parse_status %in%
        c(
          "NO_READS_PARSED",
          "NO_SEQUENCES_PARSED"
        )
      )
    )
  )

  write_excel_csv(
    summary,
    file.path(
      OUT,
      "V2_STEP83B_summary.csv"
    ),
    na=""
  )

  # ==========================================================
  # README
  # ==========================================================

  readme <- c(
    "SEPSIS V2 - STEP83B READ QUALITY / PRIMER AUDIT",
    paste0("Created: ",Sys.time()),
    "",
    "STRICT DATA ROOT",
    CANON,
    "",
    "NO BROAD FILESYSTEM SCAN WAS PERFORMED.",
    "",
    "WHAT THIS STEP DOES",
    "- Samples a small number of canonical sequence files per project.",
    "- Reads up to 1000 sequences per sampled file.",
    "- Measures read length distributions.",
    "- Measures Phred+33 FASTQ quality profiles.",
    "- Calculates quality-by-cycle summaries.",
    "- Searches the first 45 bases for common 16S primer motifs.",
    "- Generates project-level readiness for cohort-specific DADA2 parameter selection.",
    "",
    "WHAT THIS STEP DOES NOT DO",
    "- It does not run DADA2.",
    "- It does not alter raw reads.",
    "- It does not assume one truncLen/trimLeft configuration for all cohorts.",
    "- It does not force primer identity if primers were already removed.",
    "",
    "IMPORTANT",
    "- CRA002354 FASTA files have no per-base quality scores and remain a separate route.",
    "- Primer not detected can mean primers were already trimmed; it is not automatically an error.",
    "- DADA2 parameters should be chosen cohort by cohort from this audit.",
    "",
    "KEY OUTPUT",
    "V2_STEP83B_project_parameter_audit.csv",
    "",
    "NEXT",
    "Step83C will convert the read/quality/primer evidence into a frozen processing parameter sheet for each cohort that actually needs new DADA2."
  )

  writeLines(
    readme,
    file.path(
      OUT,
      "README_STEP83B_READ_QUALITY_PRIMER_AUDIT.txt"
    ),
    useBytes=TRUE
  )

  ck("STEP83B COMPLETE")

  cat("\n============================================================\n")
  cat("SEPSIS V2 - STEP83B READ/QUALITY/PRIMER AUDIT COMPLETE\n")
  cat("============================================================\n\n")

  print(
    summary,
    n=Inf,
    width=Inf
  )

  cat("\nPROJECT PARAMETER AUDIT:\n")

  print(
    project_audit,
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
        "STEP83B FATAL ERROR: ",
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

    ck("STEP83B FAILED")

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
