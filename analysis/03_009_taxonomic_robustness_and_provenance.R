# ============================================================
# Sepsis V2 - Step 89B
# Taxonomic robustness, patient-level coherence, and
# candidate ASV provenance audit
#
# PURPOSE
# 1) Quantify cross-cohort concordance of taxonomic effect sizes.
# 2) Test whether patients move in a common taxonomic direction
#    or reach ecological displacement through heterogeneous routes.
# 3) Audit the ASV-level provenance of prioritized sepsis taxa
#    before biological interpretation.
#
# INPUTS
# - Step87B cohort-specific analysis objects
# - Step88B frozen contrast registry
# - Step89A2 common-anchor taxonomic results
#
# IMPORTANT
# - No raw ASV merging across cohorts.
# - Common-anchor paired contrasts remain primary.
# - Genus/Family effect directions can be compared across cohorts,
#   but raw relative-abundance magnitudes are not pooled.
# - Candidate ASV sequences are exported for later external
#   sequence-identity verification if needed.
# ============================================================

options(stringsAsFactors = FALSE)

pkgs <- c(
  "readr",
  "dplyr",
  "tidyr",
  "tibble",
  "purrr",
  "stringr"
)

missing_pkgs <- pkgs[
  !vapply(
    pkgs,
    requireNamespace,
    logical(1),
    quietly = TRUE
  )
]

if (length(missing_pkgs)) {
  install.packages(
    missing_pkgs,
    repos = "https://cloud.r-project.org"
  )
}

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(purrr)
  library(stringr)
})

set.seed(20260820)

ROOT <- "E:/sepsis_project"

STEP87B <- file.path(
  ROOT,
  "results",
  "V2_27B_ANALYSIS_OBJECTS_AND_REPLICATE_FREEZE"
)

STEP88B <- file.path(
  ROOT,
  "results",
  "V2_28B_STEP88B_ANCHOR_ROBUSTNESS_AND_PAIRED_CONTRASTS"
)

STEP89A2 <- file.path(
  ROOT,
  "results",
  "V2_29A2_STEP89A_TAXONOMIC_PAIRED_TRAJECTORIES_FIXED"
)

OUT <- file.path(
  ROOT,
  "results",
  "V2_29B_STEP89B_TAXONOMIC_ROBUSTNESS_AND_PROVENANCE"
)

DIR_CONC <- file.path(
  OUT,
  "01_CROSS_COHORT_CONCORDANCE"
)

DIR_COH <- file.path(
  OUT,
  "02_PATIENT_LEVEL_COHERENCE"
)

DIR_PROV <- file.path(
  OUT,
  "03_CANDIDATE_ASV_PROVENANCE"
)

DIR_SYN <- file.path(
  OUT,
  "04_SYNTHESIS"
)

for (d in c(
  OUT,
  DIR_CONC,
  DIR_COH,
  DIR_PROV,
  DIR_SYN
)) {
  dir.create(
    d,
    recursive = TRUE,
    showWarnings = FALSE
  )
}

LOG <- file.path(
  OUT,
  "_STEP89B_runtime_checkpoints.txt"
)

ERR <- file.path(
  OUT,
  "_STEP89B_FATAL_ERROR.txt"
)

if (file.exists(ERR)) {
  unlink(ERR)
}

ck <- function(x) {
  cat(
    paste0(
      x,
      ": ",
      Sys.time(),
      "\n"
    ),
    file = LOG,
    append = TRUE
  )
}

safe_csv <- function(p) {
  suppressMessages(
    read_csv(
      p,
      show_col_types = FALSE,
      progress = FALSE,
      name_repair = "unique"
    )
  )
}

write_csv_safe <- function(
  x,
  p
) {
  write_excel_csv(
    x,
    p,
    na = ""
  )
}

clean_chr <- function(x) {

  y <- as.character(x)

  y[
    is.na(y) |
    trimws(y) == "" |
    tolower(
      trimws(y)
    ) %in%
      c(
        "na",
        "nan",
        "null"
      )
  ] <- NA_character_

  y
}

resolve_time_labels <- function(d) {

  tr <- clean_chr(
    d$time_raw
  )

  tk <- clean_chr(
    d$timepoint_key
  )

  label <- ifelse(
    !is.na(tr),
    tr,
    tk
  )

  map <- tibble(
    label = label,
    ord = as.numeric(
      d$analysis_time_order
    )
  ) |>
    filter(
      !is.na(label),
      !is.na(ord)
    ) |>
    distinct()

  collision <- map |>
    count(
      label,
      name = "n_orders"
    ) |>
    filter(
      n_orders > 1
    )

  if (
    nrow(collision) > 0
  ) {
    label <- ifelse(
      !is.na(tk),
      tk,
      label
    )
  }

  label
}

normalize_tax_value <- function(x) {

  y <- clean_chr(x)

  y <- str_replace(
    y,
    "^[a-zA-Z]__",
    ""
  )

  y <- str_trim(y)

  bad <- is.na(y) |
    tolower(y) %in%
      c(
        "",
        "unknown",
        "unclassified",
        "uncultured",
        "uncultured bacterium",
        "uncultured_bacterium",
        "unidentified",
        "norank",
        "no rank",
        "metagenome"
      ) |
    str_detect(
      tolower(
        ifelse(
          is.na(y),
          "",
          y
        )
      ),
      "^uncultured($|[ _])"
    )

  y[bad] <- NA_character_

  y
}

detect_tax_col <- function(
  tax,
  target
) {

  nm <- names(tax)

  hit <- nm[
    tolower(nm) ==
      tolower(target)
  ]

  if (!length(hit)) {
    stop(
      "Taxonomy table has no ",
      target,
      " column. Columns: ",
      paste(
        nm,
        collapse = ", "
      )
    )
  }

  hit[1]
}

make_tax_labels <- function(
  tax,
  rank
) {

  if (
    rank ==
      "Family"
  ) {

    fam_col <- detect_tax_col(
      tax,
      "Family"
    )

    fam <- normalize_tax_value(
      tax[[fam_col]]
    )

    return(fam)
  }

  if (
    rank ==
      "Genus"
  ) {

    gen_col <- detect_tax_col(
      tax,
      "Genus"
    )

    fam_col <- detect_tax_col(
      tax,
      "Family"
    )

    gen <- normalize_tax_value(
      tax[[gen_col]]
    )

    fam <- normalize_tax_value(
      tax[[fam_col]]
    )

    out <- gen

    fallback <- is.na(out) &
      !is.na(fam)

    out[fallback] <- paste0(
      "Unclassified_",
      fam[fallback]
    )

    return(out)
  }

  stop(
    "Unsupported rank: ",
    rank
  )
}

collapse_counts <- function(
  counts,
  tax,
  rank
) {

  labels <- make_tax_labels(
    tax,
    rank
  )

  keep <- !is.na(labels)

  if (!any(keep)) {
    stop(
      rank,
      ": no classified taxa."
    )
  }

  x <- t(
    rowsum(
      t(
        counts[
          ,
          keep,
          drop = FALSE
        ]
      ),
      group =
        labels[keep],
      reorder = FALSE
    )
  )

  storage.mode(x) <- "numeric"

  x
}

relative_abundance <- function(
  counts
) {

  rs <- rowSums(counts)

  if (
    any(
      rs <= 0
    )
  ) {
    stop(
      "Zero row sum after taxonomic collapse."
    )
  }

  counts /
    rs
}

clr_matrix <- function(
  counts,
  pseudocount = 0.5
) {

  logx <- log(
    counts +
      pseudocount
  )

  logx -
    rowMeans(logx)
}

build_common_anchor_pair <- function(
  md,
  required_anchor,
  early_label,
  late_label
) {

  d <- md |>
    filter(
      !is.na(
        analysis_time_order
      )
    )

  d$time_order_numeric <- as.numeric(
    d$analysis_time_order
  )

  d$time_label <- resolve_time_labels(d)

  patient_n <- d |>
    distinct(
      patient_id,
      timepoint_key
    ) |>
    count(
      patient_id,
      name = "n_time"
    )

  eligible <- patient_n |>
    filter(
      n_time >= 2
    ) |>
    pull(
      patient_id
    )

  d <- d |>
    filter(
      patient_id %in%
        eligible
    )

  ref <- d |>
    group_by(
      patient_id
    ) |>
    arrange(
      time_order_numeric,
      run_id,
      .by_group = TRUE
    ) |>
    slice(1) |>
    ungroup() |>
    select(
      patient_id,
      reference_time_label =
        time_label
    )

  d <- d |>
    left_join(
      ref,
      by = "patient_id"
    ) |>
    filter(
      reference_time_label ==
        required_anchor
    )

  paired_ids <- d |>
    filter(
      time_label %in%
        c(
          early_label,
          late_label
        )
    ) |>
    distinct(
      patient_id,
      time_label
    ) |>
    count(
      patient_id,
      name =
        "n_contrast_times"
    ) |>
    filter(
      n_contrast_times ==
        2
    ) |>
    pull(
      patient_id
    )

  d |>
    filter(
      patient_id %in%
        paired_ids,
      time_label %in%
        c(
          early_label,
          late_label
        )
    ) |>
    arrange(
      patient_id,
      time_order_numeric
    )
}

filter_tested_taxa <- function(
  counts,
  md_pair,
  early_label,
  late_label,
  min_prevalence = 0.20,
  min_mean_ra = 0.001
) {

  x <- counts[
    md_pair$run_id,
    ,
    drop = FALSE
  ]

  ra <- relative_abundance(x)

  er <- which(
    md_pair$time_label ==
      early_label
  )

  lr <- which(
    md_pair$time_label ==
      late_label
  )

  prev_e <- colMeans(
    x[
      er,
      ,
      drop = FALSE
    ] > 0
  )

  prev_l <- colMeans(
    x[
      lr,
      ,
      drop = FALSE
    ] > 0
  )

  mean_e <- colMeans(
    ra[
      er,
      ,
      drop = FALSE
    ]
  )

  mean_l <- colMeans(
    ra[
      lr,
      ,
      drop = FALSE
    ]
  )

  names(
    which(
      pmax(
        prev_e,
        prev_l
      ) >=
        min_prevalence &
      pmax(
        mean_e,
        mean_l
      ) >=
        min_mean_ra
    )
  )
}

patient_delta_matrix <- function(
  collapsed,
  md_pair,
  early_label,
  late_label,
  tested_taxa
) {

  x <- collapsed[
    md_pair$run_id,
    ,
    drop = FALSE
  ]

  clr <- clr_matrix(
    x,
    pseudocount = 0.5
  )

  if (
    length(
      tested_taxa
    )
  ) {
    clr <- clr[
      ,
      tested_taxa,
      drop = FALSE
    ]
  }

  pats <- unique(
    md_pair$patient_id
  )

  rows <- vector(
    "list",
    length(pats)
  )

  for (
    i in seq_along(pats)
  ) {

    p <- pats[i]

    mp <- md_pair |>
      filter(
        patient_id ==
          p
      )

    er <- mp |>
      filter(
        time_label ==
          early_label
      ) |>
      pull(run_id)

    lr <- mp |>
      filter(
        time_label ==
          late_label
      ) |>
      pull(run_id)

    if (
      length(er) != 1 ||
      length(lr) != 1
    ) {
      stop(
        "Patient ",
        p,
        " does not have exactly one early and one late sample."
      )
    }

    rows[[i]] <- clr[
      lr,
      ,
      drop = FALSE
    ] -
      clr[
        er,
        ,
        drop = FALSE
      ]

    rownames(
      rows[[i]]
    ) <- p
  }

  do.call(
    rbind,
    rows
  )
}

cosine_summary <- function(
  delta
) {

  if (
    nrow(delta) < 2 ||
    ncol(delta) < 1
  ) {
    return(
      tibble(
        n_patients =
          nrow(delta),
        n_taxa =
          ncol(delta),
        coherence_ratio =
          NA_real_,
        median_pairwise_cosine =
          NA_real_,
        q1_pairwise_cosine =
          NA_real_,
        q3_pairwise_cosine =
          NA_real_,
        fraction_pairwise_cosine_positive =
          NA_real_
      )
    )
  }

  norms <- sqrt(
    rowSums(
      delta^2
    )
  )

  valid <- norms > 0

  delta2 <- delta[
    valid,
    ,
    drop = FALSE
  ]

  norms2 <- norms[valid]

  if (
    nrow(delta2) < 2
  ) {
    return(
      tibble(
        n_patients =
          nrow(delta),
        n_taxa =
          ncol(delta),
        coherence_ratio =
          NA_real_,
        median_pairwise_cosine =
          NA_real_,
        q1_pairwise_cosine =
          NA_real_,
        q3_pairwise_cosine =
          NA_real_,
        fraction_pairwise_cosine_positive =
          NA_real_
      )
    )
  }

  unit <- delta2 /
    norms2

  sim <- unit %*%
    t(unit)

  vals <- sim[
    upper.tri(sim)
  ]

  mean_vector <- colMeans(
    delta2
  )

  coherence_ratio <- sqrt(
    sum(
      mean_vector^2
    )
  ) /
    mean(norms2)

  tibble(
    n_patients =
      nrow(delta),
    n_taxa =
      ncol(delta),
    coherence_ratio =
      coherence_ratio,
    median_pairwise_cosine =
      median(
        vals,
        na.rm = TRUE
      ),
    q1_pairwise_cosine =
      quantile(
        vals,
        0.25,
        na.rm = TRUE
      ),
    q3_pairwise_cosine =
      quantile(
        vals,
        0.75,
        na.rm = TRUE
      ),
    fraction_pairwise_cosine_positive =
      mean(
        vals > 0,
        na.rm = TRUE
      )
  )
}

safe_spearman <- function(
  x,
  y
) {

  ok <- complete.cases(
    x,
    y
  )

  x <- x[ok]
  y <- y[ok]

  if (
    length(x) < 5 ||
    length(
      unique(x)
    ) < 2 ||
    length(
      unique(y)
    ) < 2
  ) {
    return(
      c(
        rho =
          NA_real_,
        p =
          NA_real_
      )
    )
  }

  tt <- suppressWarnings(
    cor.test(
      x,
      y,
      method =
        "spearman",
      exact =
        FALSE
    )
  )

  c(
    rho =
      unname(
        tt$estimate
      ),
    p =
      tt$p.value
  )
}

topn_set <- function(
  d,
  n = 20
) {

  d |>
    filter(
      !is.na(
        paired_effect_dz
      )
    ) |>
    arrange(
      desc(
        abs(
          paired_effect_dz
        )
      )
    ) |>
    slice_head(
      n = n
    ) |>
    pull(
      taxon
    )
}

jaccard <- function(
  a,
  b
) {

  u <- union(
    a,
    b
  )

  if (!length(u)) {
    return(
      NA_real_
    )
  }

  length(
    intersect(
      a,
      b
    )
  ) /
    length(u)
}

cross_cohort_concordance <- function(
  primary,
  rank,
  projects
) {

  d <- primary |>
    filter(
      tax_rank ==
        rank,
      project %in%
        projects
    )

  pairs <- combn(
    projects,
    2,
    simplify = FALSE
  )

  rows <- list()

  for (
    i in seq_along(pairs)
  ) {

    p1 <- pairs[[i]][1]
    p2 <- pairs[[i]][2]

    d1 <- d |>
      filter(
        project ==
          p1
      ) |>
      select(
        taxon,
        dz1 =
          paired_effect_dz,
        clr1 =
          mean_clr_difference,
        dir1 =
          direction
      )

    d2 <- d |>
      filter(
        project ==
          p2
      ) |>
      select(
        taxon,
        dz2 =
          paired_effect_dz,
        clr2 =
          mean_clr_difference,
        dir2 =
          direction
      )

    m <- inner_join(
      d1,
      d2,
      by = "taxon"
    )

    s_dz <- safe_spearman(
      m$dz1,
      m$dz2
    )

    s_clr <- safe_spearman(
      m$clr1,
      m$clr2
    )

    direction_agreement <- if (
      nrow(m)
    ) {
      mean(
        m$dir1 ==
          m$dir2,
        na.rm = TRUE
      )
    } else {
      NA_real_
    }

    top1 <- topn_set(
      d |>
        filter(
          project ==
            p1
        ),
      n = 20
    )

    top2 <- topn_set(
      d |>
        filter(
          project ==
            p2
        ),
      n = 20
    )

    rows[[i]] <- tibble(
      tax_rank =
        rank,
      project_1 =
        p1,
      project_2 =
        p2,
      n_shared_tested_taxa =
        nrow(m),
      spearman_dz_rho =
        s_dz["rho"],
      spearman_dz_p =
        s_dz["p"],
      spearman_mean_clr_rho =
        s_clr["rho"],
      spearman_mean_clr_p =
        s_clr["p"],
      direction_agreement_fraction =
        direction_agreement,
      top20_jaccard =
        jaccard(
          top1,
          top2
        )
    )
  }

  bind_rows(rows)
}

candidate_provenance <- function(
  obj,
  candidates
) {

  counts <- obj$counts
  tax <- as.data.frame(
    obj$taxonomy,
    stringsAsFactors = FALSE
  )

  if (
    ncol(counts) !=
      nrow(tax)
  ) {
    stop(
      "Candidate provenance: count/taxonomy dimension mismatch."
    )
  }

  seqs <- colnames(counts)

  lib <- rowSums(counts)
  asv_ra <- counts /
    lib

  lineage_cols <- names(tax)[
    tolower(
      names(tax)
    ) %in%
      c(
        "kingdom",
        "domain",
        "phylum",
        "class",
        "order",
        "family",
        "genus",
        "species"
      )
  ]

  summary_rows <- list()
  detail_rows <- list()
  fasta_lines <- character(0)

  k <- 0L

  for (
    i in seq_len(
      nrow(candidates)
    )
  ) {

    rank <- candidates$tax_rank[i]
    taxon <- candidates$taxon[i]

    labels <- make_tax_labels(
      tax,
      rank
    )

    idx <- which(
      labels ==
        taxon
    )

    if (!length(idx)) {
      next
    }

    totals <- colSums(
      counts[
        ,
        idx,
        drop = FALSE
      ]
    )

    total_taxon_reads <- sum(totals)

    dominant_local <- which.max(
      totals
    )

    dominant_idx <- idx[
      dominant_local
    ]

    dominant_fraction <- if (
      total_taxon_reads > 0
    ) {
      totals[
        dominant_local
      ] /
        total_taxon_reads
    } else {
      NA_real_
    }

    taxon_counts <- rowSums(
      counts[
        ,
        idx,
        drop = FALSE
      ]
    )

    taxon_ra <- rowSums(
      asv_ra[
        ,
        idx,
        drop = FALSE
      ]
    )

    summary_rows[[
      length(
        summary_rows
      ) + 1
    ]] <- tibble(
      tax_rank =
        rank,
      taxon =
        taxon,
      n_contributing_ASVs =
        length(idx),
      total_reads_all_primary_samples =
        total_taxon_reads,
      prevalence_all_primary_samples =
        mean(
          taxon_counts > 0
        ),
      mean_relative_abundance_all_primary_samples =
        mean(
          taxon_ra
        ),
      median_relative_abundance_all_primary_samples =
        median(
          taxon_ra
        ),
      dominant_ASV_sequence =
        seqs[
          dominant_idx
        ],
      dominant_ASV_read_fraction_within_taxon =
        dominant_fraction,
      single_ASV_dominant_flag =
        !is.na(
          dominant_fraction
        ) &&
        dominant_fraction >=
          0.80,
      manual_sequence_identity_check_recommended =
        TRUE
    )

    for (
      j in seq_along(idx)
    ) {

      ii <- idx[j]

      row <- tibble(
        tax_rank =
          rank,
        taxon =
          taxon,
        ASV_sequence =
          seqs[ii],
        total_reads =
          totals[j],
        read_fraction_within_taxon =
          if (
            total_taxon_reads > 0
          ) {
            totals[j] /
              total_taxon_reads
          } else {
            NA_real_
          },
        prevalence =
          mean(
            counts[
              ,
              ii
            ] > 0
          ),
        mean_relative_abundance =
          mean(
            asv_ra[
              ,
              ii
            ]
          )
      )

      for (
        lc in lineage_cols
      ) {
        row[[lc]] <- as.character(
          tax[
            ii,
            lc
          ]
        )
      }

      detail_rows[[
        length(
          detail_rows
        ) + 1
      ]] <- row

      k <- k + 1L

      fasta_lines <- c(
        fasta_lines,
        paste0(
          ">candidate_",
          k,
          "|",
          rank,
          "|",
          str_replace_all(
            taxon,
            "[^A-Za-z0-9._-]",
            "_"
          ),
          "|reads=",
          totals[j]
        ),
        seqs[ii]
      )
    }
  }

  list(
    summary =
      bind_rows(
        summary_rows
      ),
    detail =
      bind_rows(
        detail_rows
      ),
    fasta =
      fasta_lines
  )
}

main <- function() {

  ck(
    "STEP89B STARTED"
  )

  registry_path <- file.path(
    STEP87B,
    "V2_STEP87B_analysis_object_registry.csv"
  )

  step87_complete <- file.path(
    STEP87B,
    "_STEP87B_COMPLETE.ok"
  )

  step88_complete <- file.path(
    STEP88B,
    "_STEP88B_COMPLETE.ok"
  )

  contrast_path <- file.path(
    STEP88B,
    "02_PAIRED_CONTRASTS",
    "V2_STEP88B_paired_early_vs_late_followup_contrasts.csv"
  )

  step89_complete <- file.path(
    STEP89A2,
    "_STEP89A2_COMPLETE.ok"
  )

  primary_path <- file.path(
    STEP89A2,
    "04_CROSS_COHORT_SYNTHESIS",
    "V2_STEP89A2_PRIMARY_common_anchor_all_taxa.csv"
  )

  genus_syn_path <- file.path(
    STEP89A2,
    "04_CROSS_COHORT_SYNTHESIS",
    "V2_STEP89A2_GENUS_cross_cohort_trajectory_synthesis.csv"
  )

  family_syn_path <- file.path(
    STEP89A2,
    "04_CROSS_COHORT_SYNTHESIS",
    "V2_STEP89A2_FAMILY_cross_cohort_trajectory_synthesis.csv"
  )

  core_path <- file.path(
    STEP89A2,
    "04_CROSS_COHORT_SYNTHESIS",
    "V2_STEP89A2_CORE_SEPSIS_taxonomic_shortlist.csv"
  )

  required <- c(
    registry_path,
    step87_complete,
    step88_complete,
    contrast_path,
    step89_complete,
    primary_path,
    genus_syn_path,
    family_syn_path,
    core_path
  )

  if (
    !all(
      file.exists(
        required
      )
    )
  ) {
    stop(
      paste0(
        "Required upstream files are missing:\n",
        paste(
          required[
            !file.exists(
              required
            )
          ],
          collapse = "\n"
        )
      )
    )
  }

  registry <- safe_csv(
    registry_path
  )

  contrasts_raw <- safe_csv(
    contrast_path
  )

  primary <- safe_csv(
    primary_path
  )

  genus_syn <- safe_csv(
    genus_syn_path
  )

  family_syn <- safe_csv(
    family_syn_path
  )

  core <- safe_csv(
    core_path
  )

  if (
    nrow(registry) != 8 ||
    sum(
      registry$primary_samples
    ) != 785
  ) {
    stop(
      "Step87B registry guard failed."
    )
  }

  if (
    nrow(
      primary |>
        filter(
          project ==
            "PRJNA691455",
          tax_rank ==
            "Genus",
          fdr <
            0.10
        )
    ) != 0
  ) {
    stop(
      "Step89A2 core sepsis Genus q<0.10 state changed unexpectedly."
    )
  }

  if (
    nrow(
      primary |>
        filter(
          project ==
            "PRJNA691455",
          tax_rank ==
            "Family",
          fdr <
            0.10
        )
    ) != 4
  ) {
    stop(
      "Step89A2 core sepsis Family q<0.10 state changed unexpectedly."
    )
  }

  ck(
    "UPSTREAM STATE GUARDED"
  )

  contrast_registry <- contrasts_raw |>
    filter(
      metric ==
        "bray_from_patient_baseline",
      subset_type ==
        "ALL_PAIRED"
    ) |>
    distinct(
      project,
      analysis_role_group,
      required_anchor,
      early_label,
      late_label
    )

  natural_history_projects <- c(
    "PRJNA691455",
    "PRJNA851469",
    "PRJNA516701",
    "PRJNA578267"
  )

  expected_pairs <- tribble(
    ~project,
    ~n_pairs_expected,

    "PRJNA691455",
    9L,

    "PRJNA851469",
    14L,

    "PRJNA516701",
    14L,

    "PRJNA578267",
    32L
  )

  # ----------------------------------------------------------
  # 1. Cross-cohort effect concordance
  # ----------------------------------------------------------

  concordance <- bind_rows(
    cross_cohort_concordance(
      primary,
      "Genus",
      natural_history_projects
    ),
    cross_cohort_concordance(
      primary,
      "Family",
      natural_history_projects
    )
  )

  write_csv_safe(
    concordance,
    file.path(
      DIR_CONC,
      "V2_STEP89B_cross_cohort_effect_concordance.csv"
    )
  )

  # ----------------------------------------------------------
  # 2. Patient-level delta-vector coherence
  # ----------------------------------------------------------

  coherence_rows <- list()

  for (
    proj in natural_history_projects
  ) {

    ck(
      paste0(
        proj,
        " PATIENT COHERENCE START"
      )
    )

    cr <- contrast_registry |>
      filter(
        project ==
          proj
      )

    if (
      nrow(cr) != 1
    ) {
      stop(
        proj,
        ": contrast registry row missing/non-unique."
      )
    }

    rr <- registry |>
      filter(
        project ==
          proj
      )

    object_path <- as.character(
      rr$analysis_object_path[1]
    )

    if (
      !file.exists(
        object_path
      )
    ) {
      stop(
        proj,
        ": analysis object missing: ",
        object_path
      )
    }

    obj <- readRDS(
      object_path
    )

    counts <- obj$counts

    tax <- as.data.frame(
      obj$taxonomy,
      stringsAsFactors = FALSE
    )

    md <- as_tibble(
      obj$metadata
    )

    md_pair <- build_common_anchor_pair(
      md =
        md,
      required_anchor =
        cr$required_anchor,
      early_label =
        cr$early_label,
      late_label =
        cr$late_label
    )

    n_pairs <- n_distinct(
      md_pair$patient_id
    )

    expected_n <- expected_pairs |>
      filter(
        project ==
          proj
      ) |>
      pull(
        n_pairs_expected
      )

    if (
      n_pairs !=
        expected_n
    ) {
      stop(
        proj,
        ": common-anchor pair count mismatch. Expected ",
        expected_n,
        "; observed ",
        n_pairs,
        "."
      )
    }

    for (
      rank in c(
        "Genus",
        "Family"
      )
    ) {

      collapsed <- collapse_counts(
        counts,
        tax,
        rank
      )

      tested <- filter_tested_taxa(
        counts =
          collapsed,
        md_pair =
          md_pair,
        early_label =
          cr$early_label,
        late_label =
          cr$late_label,
        min_prevalence =
          0.20,
        min_mean_ra =
          0.001
      )

      delta <- patient_delta_matrix(
        collapsed =
          collapsed,
        md_pair =
          md_pair,
        early_label =
          cr$early_label,
        late_label =
          cr$late_label,
        tested_taxa =
          tested
      )

      coh <- cosine_summary(
        delta
      ) |>
        mutate(
          project =
            proj,
          tax_rank =
            rank,
          required_anchor =
            cr$required_anchor,
          early_label =
            cr$early_label,
          late_label =
            cr$late_label,
          .before = 1
        )

      coherence_rows[[
        length(
          coherence_rows
        ) + 1
      ]] <- coh

      delta_df <- as.data.frame(
        delta,
        check.names = FALSE
      ) |>
        rownames_to_column(
          "patient_id"
        )

      write_csv_safe(
        delta_df,
        file.path(
          DIR_COH,
          paste0(
            proj,
            "_",
            toupper(rank),
            "_patient_CLR_delta_matrix.csv"
          )
        )
      )
    }

    ck(
      paste0(
        proj,
        " PATIENT COHERENCE COMPLETE"
      )
    )

    rm(
      obj,
      counts,
      tax,
      md,
      md_pair
    )

    gc(
      verbose = FALSE
    )
  }

  coherence <- bind_rows(
    coherence_rows
  )

  write_csv_safe(
    coherence,
    file.path(
      DIR_COH,
      "V2_STEP89B_patient_taxonomic_direction_coherence.csv"
    )
  )

  # ----------------------------------------------------------
  # 3. Candidate definition for provenance audit
  # ----------------------------------------------------------

  core_candidates <- core |>
    filter(
      subset_type ==
        "COMMON_ANCHOR_SENSITIVITY",
      fdr <
        0.10 |
        wilcoxon_p <
          0.05
    ) |>
    select(
      tax_rank,
      taxon
    ) |>
    distinct()

  syn_candidates <- bind_rows(
    genus_syn |>
      select(
        tax_rank,
        taxon,
        synthesis_pattern
      ),
    family_syn |>
      select(
        tax_rank,
        taxon,
        synthesis_pattern
      )
  ) |>
    filter(
      synthesis_pattern %in%
        c(
          "CRITICAL_ILLNESS_CONSISTENT_NONSEPSIS_OPPOSITE_DIRECTIONAL",
          "CRITICAL_ILLNESS_CONSISTENT_NONSEPSIS_OPPOSITE_STRONG",
          "CRITICAL_ILLNESS_CONSISTENT_DIRECTIONAL",
          "CRITICAL_ILLNESS_CONSISTENT_STRONG"
        )
    ) |>
    select(
      tax_rank,
      taxon
    ) |>
    distinct()

  candidates <- bind_rows(
    core_candidates,
    syn_candidates
  ) |>
    distinct() |>
    arrange(
      tax_rank,
      taxon
    )

  write_csv_safe(
    candidates,
    file.path(
      DIR_PROV,
      "V2_STEP89B_candidate_taxa_for_ASV_audit.csv"
    )
  )

  # ----------------------------------------------------------
  # 4. ASV provenance in the core sepsis cohort
  # ----------------------------------------------------------

  rr_core <- registry |>
    filter(
      project ==
        "PRJNA691455"
    )

  if (
    nrow(rr_core) != 1
  ) {
    stop(
      "Core sepsis registry row missing/non-unique."
    )
  }

  core_object_path <- as.character(
    rr_core$analysis_object_path[1]
  )

  if (
    !file.exists(
      core_object_path
    )
  ) {
    stop(
      "Core sepsis analysis object missing: ",
      core_object_path
    )
  }

  core_obj <- readRDS(
    core_object_path
  )

  prov <- candidate_provenance(
    core_obj,
    candidates
  )

  write_csv_safe(
    prov$summary,
    file.path(
      DIR_PROV,
      "V2_STEP89B_CORE_SEPSIS_candidate_ASV_provenance_summary.csv"
    )
  )

  write_csv_safe(
    prov$detail,
    file.path(
      DIR_PROV,
      "V2_STEP89B_CORE_SEPSIS_candidate_ASV_provenance_detail.csv"
    )
  )

  writeLines(
    prov$fasta,
    file.path(
      DIR_PROV,
      "V2_STEP89B_CORE_SEPSIS_candidate_ASVs_for_sequence_identity_check.fasta"
    ),
    useBytes = TRUE
  )

  # ----------------------------------------------------------
  # 5. Integrated synthesis
  # ----------------------------------------------------------

  tax_signal_summary <- primary |>
    filter(
      project %in%
        natural_history_projects
    ) |>
    group_by(
      project,
      tax_rank
    ) |>
    summarise(
      n_tested =
        n(),
      n_q_lt_005 =
        sum(
          fdr <
            0.05,
          na.rm = TRUE
        ),
      n_q_lt_010 =
        sum(
          fdr <
            0.10,
          na.rm = TRUE
        ),
      n_nominal_p_lt_005 =
        sum(
          wilcoxon_p <
            0.05,
          na.rm = TRUE
        ),
      .groups = "drop"
    ) |>
    left_join(
      coherence,
      by = c(
        "project",
        "tax_rank"
      )
    )

  write_csv_safe(
    tax_signal_summary,
    file.path(
      DIR_SYN,
      "V2_STEP89B_taxonomic_signal_and_coherence_summary.csv"
    )
  )

  # Flag whether the core finding is ecological-without-universal-taxon.
  core_genus <- tax_signal_summary |>
    filter(
      project ==
        "PRJNA691455",
      tax_rank ==
        "Genus"
    )

  core_family <- tax_signal_summary |>
    filter(
      project ==
        "PRJNA691455",
      tax_rank ==
        "Family"
    )

  interpretation_state <- case_when(
    core_genus$n_q_lt_010 == 0 &
      core_family$n_q_lt_005 == 0 ~
      "ROBUST_ECOLOGICAL_DISPLACEMENT_WITHOUT_STRONG_CORE_TAXONOMIC_SIGNATURE",

    TRUE ~
      "CORE_TAXONOMIC_SIGNATURE_REQUIRES_REVIEW"
  )

  readme <- c(
    "SEPSIS V2 - STEP89B TAXONOMIC ROBUSTNESS AND PROVENANCE",
    paste0(
      "Created: ",
      Sys.time()
    ),
    "",
    "WHY THIS STEP",
    "Step89A2 found robust ecological displacement but no core sepsis Genus q<0.10 and no Family q<0.05.",
    "Step89B quantifies whether taxonomic effect directions generalize and whether individual patients share a common taxonomic route.",
    "",
    "CROSS-COHORT CONCORDANCE",
    "Spearman correlations, direction agreement, and top-20 effect overlap are calculated for common tested Genus/Family taxa.",
    "No raw ASVs are pooled across cohorts.",
    "",
    "PATIENT COHERENCE",
    "Each patient is represented by a late-minus-early CLR taxonomic delta vector.",
    "coherence_ratio approaches 1 when patients move in a common taxonomic direction and approaches 0 when their taxonomic routes cancel across patients.",
    "Pairwise cosine similarity is exported as a complementary measure.",
    "",
    "CANDIDATE PROVENANCE",
    "Core sepsis nominal/FDR candidates and cross-cohort directional candidates are traced back to their contributing ASVs.",
    "Dominant-ASV fraction and ASV sequences are exported before biological interpretation.",
    "",
    "INTERPRETATION STATE",
    interpretation_state,
    "",
    "NEXT",
    "Use the provenance and coherence results to decide whether the manuscript should emphasize specific taxa or instead emphasize conserved ecological instability with heterogeneous taxonomic realization."
  )

  writeLines(
    readme,
    file.path(
      OUT,
      "README_STEP89B.txt"
    ),
    useBytes = TRUE
  )

  writeLines(
    c(
      paste0(
        "Completed: ",
        Sys.time()
      ),
      "Status: STEP89B COMPLETE",
      "Four natural-history cohorts audited for taxonomic concordance.",
      "Patient-level CLR delta-vector coherence calculated at Genus and Family levels.",
      paste0(
        "Candidate taxa sent to core-sepsis ASV provenance audit: ",
        nrow(candidates)
      ),
      paste0(
        "Core-sepsis candidate ASV detail rows: ",
        nrow(
          prov$detail
        )
      ),
      paste0(
        "Interpretation state: ",
        interpretation_state
      )
    ),
    file.path(
      OUT,
      "_STEP89B_COMPLETE.ok"
    ),
    useBytes = TRUE
  )

  ck(
    "STEP89B COMPLETE"
  )

  cat(
    "\n============================================================\n"
  )
  cat(
    "SEPSIS V2 - STEP89B COMPLETE\n"
  )
  cat(
    "============================================================\n\n"
  )

  cat(
    "Taxonomic signal/coherence summary:\n"
  )

  print(
    tax_signal_summary,
    n = Inf,
    width = Inf
  )

  cat(
    "\nCross-cohort concordance:\n"
  )

  print(
    concordance,
    n = Inf,
    width = Inf
  )

  cat(
    "\nCandidate provenance summary:\n"
  )

  print(
    prov$summary,
    n = Inf,
    width = Inf
  )

  cat(
    "\nOutput directory:\n"
  )

  cat(
    OUT,
    "\n"
  )
}

tryCatch(
  main(),
  error = function(e) {

    msg <- c(
      paste0(
        "STEP89B FATAL ERROR: ",
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
          collapse = " "
        )
      )
    )

    writeLines(
      msg,
      ERR,
      useBytes = TRUE
    )

    ck(
      "STEP89B FAILED"
    )

    message(
      paste(
        msg,
        collapse = "\n"
      )
    )

    quit(
      save = "no",
      status = 1,
      runLast = FALSE
    )
  }
)
