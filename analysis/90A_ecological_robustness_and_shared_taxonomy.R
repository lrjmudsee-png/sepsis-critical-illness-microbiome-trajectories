# ============================================================
# Sepsis V2 - Step 90A
# Ecological robustness + shared-taxonomy coherence audit
#
# PURPOSE
# 1) Stress-test the common-anchor Bray trajectory result using:
#    - paired sign consistency
#    - bootstrap CIs
#    - leave-one-patient-out (LOO) analysis
# 2) Recalculate patient taxonomic coherence in a SHARED
#    Genus/Family feature space across the four natural-history
#    cohorts, so cohort differences are not driven simply by
#    different numbers of tested taxa.
# 3) Audit whether core sepsis nominal taxonomic candidates are
#    broad patient-level signals or leverage/outlier-driven.
#
# This is a robustness step; it does not rerun upstream
# microbiome processing and does not merge raw ASVs.
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

STEP88A2 <- file.path(
  ROOT,
  "results",
  "V2_28A2_STEP88A_LONGITUDINAL_DIVERSITY_AND_DISPLACEMENT_FIXED"
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

STEP89B <- file.path(
  ROOT,
  "results",
  "V2_29B_STEP89B_TAXONOMIC_ROBUSTNESS_AND_PROVENANCE"
)

OUT <- file.path(
  ROOT,
  "results",
  "V2_30A_STEP90A_ECOLOGICAL_ROBUSTNESS_AND_SHARED_TAXONOMY"
)

DIR_BRAY <- file.path(
  OUT,
  "01_BRAY_ROBUSTNESS"
)

DIR_SHARED <- file.path(
  OUT,
  "02_SHARED_TAXONOMY_COHERENCE"
)

DIR_TAXON <- file.path(
  OUT,
  "03_CORE_SEPSIS_TAXON_BREADTH"
)

DIR_SYN <- file.path(
  OUT,
  "04_SYNTHESIS"
)

for (d in c(
  OUT,
  DIR_BRAY,
  DIR_SHARED,
  DIR_TAXON,
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
  "_STEP90A_runtime_checkpoints.txt"
)

ERR <- file.path(
  OUT,
  "_STEP90A_FATAL_ERROR.txt"
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
      name_repair = "minimal"
    )
  )
}

write_csv_safe <- function(x, p) {
  write_excel_csv(
    x,
    p,
    na = ""
  )
}

safe_wilcox <- function(x) {

  x <- as.numeric(
    x[
      !is.na(x)
    ]
  )

  if (!length(x)) {
    return(NA_real_)
  }

  if (
    all(
      abs(x) <
        1e-15
    )
  ) {
    return(1)
  }

  tryCatch(
    suppressWarnings(
      wilcox.test(
        x,
        mu = 0,
        exact = FALSE,
        correct = FALSE
      )$p.value
    ),
    error = function(e)
      NA_real_
  )
}

exact_sign_p <- function(x) {

  x <- as.numeric(
    x[
      !is.na(x)
    ]
  )

  pos <- sum(
    x > 0
  )

  neg <- sum(
    x < 0
  )

  n <- pos + neg

  if (!n) {
    return(NA_real_)
  }

  binom.test(
    pos,
    n,
    p = 0.5,
    alternative = "two.sided"
  )$p.value
}

boot_stat_ci <- function(
  x,
  stat = c(
    "mean",
    "median"
  ),
  B = 3000,
  seed = 20260820
) {

  stat <- match.arg(stat)

  x <- as.numeric(
    x[
      !is.na(x)
    ]
  )

  if (
    length(x) < 3
  ) {
    return(
      c(
        low = NA_real_,
        high = NA_real_
      )
    )
  }

  set.seed(seed)

  vals <- replicate(
    B,
    {
      z <- sample(
        x,
        length(x),
        replace = TRUE
      )

      if (
        stat == "mean"
      ) {
        mean(z)
      } else {
        median(z)
      }
    }
  )

  q <- quantile(
    vals,
    c(
      0.025,
      0.975
    ),
    na.rm = TRUE,
    names = FALSE
  )

  c(
    low = q[1],
    high = q[2]
  )
}

paired_summary <- function(
  pair,
  metric,
  early_label,
  late_label,
  project
) {

  wide <- pair |>
    select(
      patient_id,
      time_factor,
      all_of(metric)
    ) |>
    distinct(
      patient_id,
      time_factor,
      .keep_all = TRUE
    ) |>
    pivot_wider(
      names_from =
        time_factor,
      values_from =
        all_of(metric)
    )

  if (
    !(early_label %in%
        names(wide)) ||
    !(late_label %in%
        names(wide))
  ) {
    stop(
      project,
      ": missing contrast labels for ",
      metric,
      "."
    )
  }

  wide <- wide |>
    filter(
      !is.na(
        .data[[early_label]]
      ),
      !is.na(
        .data[[late_label]]
      )
    )

  delta <- wide[[late_label]] -
    wide[[early_label]]

  n <- length(delta)

  mean_ci <- boot_stat_ci(
    delta,
    "mean",
    seed =
      20260820 + n
  )

  med_ci <- boot_stat_ci(
    delta,
    "median",
    seed =
      20261820 + n
  )

  sd_delta <- if (
    n >= 2
  ) {
    sd(delta)
  } else {
    NA_real_
  }

  dz <- if (
    !is.na(sd_delta) &&
    sd_delta > 0
  ) {
    mean(delta) /
      sd_delta
  } else {
    NA_real_
  }

  tibble(
    project =
      project,
    metric =
      metric,
    n_pairs =
      n,
    early_median =
      median(
        wide[[early_label]],
        na.rm = TRUE
      ),
    late_median =
      median(
        wide[[late_label]],
        na.rm = TRUE
      ),
    mean_delta =
      mean(
        delta,
        na.rm = TRUE
      ),
    median_delta =
      median(
        delta,
        na.rm = TRUE
      ),
    mean_delta_ci_low =
      mean_ci["low"],
    mean_delta_ci_high =
      mean_ci["high"],
    median_delta_ci_low =
      med_ci["low"],
    median_delta_ci_high =
      med_ci["high"],
    paired_effect_dz =
      dz,
    n_positive =
      sum(
        delta > 0,
        na.rm = TRUE
      ),
    n_negative =
      sum(
        delta < 0,
        na.rm = TRUE
      ),
    n_tie =
      sum(
        abs(delta) <
          1e-15,
        na.rm = TRUE
      ),
    positive_fraction_non_ties =
      ifelse(
        sum(delta != 0) > 0,
        sum(delta > 0) /
          sum(delta != 0),
        NA_real_
      ),
    sign_test_p =
      exact_sign_p(delta),
    wilcoxon_p =
      safe_wilcox(delta)
  )
}

loo_bray <- function(
  pair,
  early_label,
  late_label,
  project
) {

  metric <-
    "bray_from_patient_baseline"

  wide <- pair |>
    select(
      patient_id,
      time_factor,
      all_of(metric)
    ) |>
    distinct(
      patient_id,
      time_factor,
      .keep_all = TRUE
    ) |>
    pivot_wider(
      names_from =
        time_factor,
      values_from =
        all_of(metric)
    ) |>
    filter(
      !is.na(
        .data[[early_label]]
      ),
      !is.na(
        .data[[late_label]]
      )
    ) |>
    mutate(
      delta =
        .data[[late_label]] -
        .data[[early_label]]
    )

  full_mean <- mean(
    wide$delta
  )

  rows <- vector(
    "list",
    nrow(wide)
  )

  for (
    i in seq_len(
      nrow(wide)
    )
  ) {

    d <- wide$delta[-i]

    rows[[i]] <- tibble(
      project =
        project,
      omitted_patient_id =
        wide$patient_id[i],
      n_remaining =
        length(d),
      mean_delta =
        mean(d),
      median_delta =
        median(d),
      n_positive =
        sum(d > 0),
      n_negative =
        sum(d < 0),
      positive_fraction_non_ties =
        ifelse(
          sum(d != 0) > 0,
          sum(d > 0) /
            sum(d != 0),
          NA_real_
        ),
      wilcoxon_p =
        safe_wilcox(d),
      sign_test_p =
        exact_sign_p(d),
      direction_preserved =
        sign(
          mean(d)
        ) ==
          sign(
            full_mean
          )
    )
  }

  bind_rows(rows)
}

cosine_metrics <- function(mat) {

  x <- as.matrix(mat)

  if (
    nrow(x) < 2 ||
    ncol(x) < 1
  ) {
    return(
      c(
        coherence_ratio =
          NA_real_,
        median_pairwise_cosine =
          NA_real_,
        fraction_pairwise_cosine_positive =
          NA_real_
      )
    )
  }

  norms <- sqrt(
    rowSums(
      x^2
    )
  )

  keep <- norms > 0

  x <- x[
    keep,
    ,
    drop = FALSE
  ]

  norms <- norms[keep]

  if (
    nrow(x) < 2
  ) {
    return(
      c(
        coherence_ratio =
          NA_real_,
        median_pairwise_cosine =
          NA_real_,
        fraction_pairwise_cosine_positive =
          NA_real_
      )
    )
  }

  unit <- x /
    norms

  sim <- unit %*%
    t(unit)

  vals <- sim[
    upper.tri(sim)
  ]

  mean_vec <- colMeans(x)

  coherence <- sqrt(
    sum(
      mean_vec^2
    )
  ) /
    mean(norms)

  c(
    coherence_ratio =
      coherence,
    median_pairwise_cosine =
      median(
        vals,
        na.rm = TRUE
      ),
    fraction_pairwise_cosine_positive =
      mean(
        vals > 0,
        na.rm = TRUE
      )
  )
}

bootstrap_cosine <- function(
  mat,
  B = 3000,
  seed = 20260820
) {

  x <- as.matrix(mat)

  n <- nrow(x)

  if (
    n < 3
  ) {
    return(
      tibble(
        coherence_ci_low =
          NA_real_,
        coherence_ci_high =
          NA_real_,
        median_cosine_ci_low =
          NA_real_,
        median_cosine_ci_high =
          NA_real_
      )
    )
  }

  set.seed(seed)

  vals <- replicate(
    B,
    {
      idx <- sample(
        seq_len(n),
        n,
        replace = TRUE
      )

      cosine_metrics(
        x[
          idx,
          ,
          drop = FALSE
        ]
      )[
        c(
          "coherence_ratio",
          "median_pairwise_cosine"
        )
      ]
    }
  )

  tibble(
    coherence_ci_low =
      quantile(
        vals[
          1,
        ],
        0.025,
        na.rm = TRUE
      ),
    coherence_ci_high =
      quantile(
        vals[
          1,
        ],
        0.975,
        na.rm = TRUE
      ),
    median_cosine_ci_low =
      quantile(
        vals[
          2,
        ],
        0.025,
        na.rm = TRUE
      ),
    median_cosine_ci_high =
      quantile(
        vals[
          2,
        ],
        0.975,
        na.rm = TRUE
      )
  )
}

candidate_breadth <- function(
  delta_df,
  candidates,
  rank
) {

  taxa <- intersect(
    candidates$taxon[
      candidates$tax_rank ==
        rank
    ],
    setdiff(
      names(delta_df),
      "patient_id"
    )
  )

  rows <- list()

  for (
    taxon in taxa
  ) {

    d <- delta_df[[taxon]]

    n <- sum(
      !is.na(d)
    )

    d2 <- d[
      !is.na(d)
    ]

    total_abs <- sum(
      abs(d2)
    )

    max_leverage <- if (
      total_abs > 0
    ) {
      max(
        abs(d2)
      ) /
        total_abs
    } else {
      NA_real_
    }

    full_mean <- mean(
      d2
    )

    loo_mean <- if (
      length(d2) >= 3
    ) {
      sapply(
        seq_along(d2),
        function(i)
          mean(
            d2[-i]
          )
      )
    } else {
      NA_real_
    }

    loo_direction_preserved <- if (
      length(loo_mean) > 1
    ) {
      mean(
        sign(
          loo_mean
        ) ==
          sign(
            full_mean
          ),
        na.rm = TRUE
      )
    } else {
      NA_real_
    }

    rows[[
      length(rows) + 1
    ]] <- tibble(
      tax_rank =
        rank,
      taxon =
        taxon,
      n_patients =
        n,
      mean_delta =
        mean(d2),
      median_delta =
        median(d2),
      n_positive =
        sum(
          d2 > 0
        ),
      n_negative =
        sum(
          d2 < 0
        ),
      n_tie =
        sum(
          abs(d2) <
            1e-15
        ),
      positive_fraction_non_ties =
        ifelse(
          sum(d2 != 0) > 0,
          sum(d2 > 0) /
            sum(d2 != 0),
          NA_real_
        ),
      sign_test_p =
        exact_sign_p(d2),
      wilcoxon_p =
        safe_wilcox(d2),
      max_single_patient_absolute_contribution =
        max_leverage,
      loo_mean_direction_preserved_fraction =
        loo_direction_preserved
    )
  }

  bind_rows(rows)
}

main <- function() {

  ck(
    "STEP90A STARTED"
  )

  disp_path <- file.path(
    STEP88A2,
    "03_BETA_DISPLACEMENT",
    "V2_STEP88A2_ALL_within_patient_bray_displacement.csv"
  )

  step88a_complete <- file.path(
    STEP88A2,
    "_STEP88A2_COMPLETE.ok"
  )

  contrast_path <- file.path(
    STEP88B,
    "02_PAIRED_CONTRASTS",
    "V2_STEP88B_paired_early_vs_late_followup_contrasts.csv"
  )

  step88b_complete <- file.path(
    STEP88B,
    "_STEP88B_COMPLETE.ok"
  )

  core_path <- file.path(
    STEP89A2,
    "04_CROSS_COHORT_SYNTHESIS",
    "V2_STEP89A2_CORE_SEPSIS_taxonomic_shortlist.csv"
  )

  step89a_complete <- file.path(
    STEP89A2,
    "_STEP89A2_COMPLETE.ok"
  )

  step89b_complete <- file.path(
    STEP89B,
    "_STEP89B_COMPLETE.ok"
  )

  required <- c(
    disp_path,
    step88a_complete,
    contrast_path,
    step88b_complete,
    core_path,
    step89a_complete,
    step89b_complete
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
        "Required upstream files missing:\n",
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

  disp <- safe_csv(
    disp_path
  )

  contrast_raw <- safe_csv(
    contrast_path
  )

  core <- safe_csv(
    core_path
  )

  natural_projects <- c(
    "PRJNA691455",
    "PRJNA851469",
    "PRJNA516701",
    "PRJNA578267"
  )

  expected_pairs <- tribble(
    ~project,
    ~n_expected,

    "PRJNA691455",
    9L,

    "PRJNA851469",
    14L,

    "PRJNA516701",
    14L,

    "PRJNA578267",
    32L
  )

  contrasts <- contrast_raw |>
    filter(
      metric ==
        "bray_from_patient_baseline",
      subset_type ==
        "ALL_PAIRED",
      project %in%
        natural_projects
    ) |>
    distinct(
      project,
      required_anchor,
      early_label,
      late_label
    )

  if (
    nrow(contrasts) != 4
  ) {
    stop(
      "Natural-history contrast registry does not contain four unique cohorts."
    )
  }

  ck(
    "UPSTREAM STATE GUARDED"
  )

  # ----------------------------------------------------------
  # 1. Bray/alpha common-anchor robustness
  # ----------------------------------------------------------

  paired_rows <- list()
  loo_rows <- list()
  loo_summary_rows <- list()

  for (
    i in seq_len(
      nrow(
        contrasts
      )
    )
  ) {

    cr <- contrasts[
      i,
      ,
      drop = FALSE
    ]

    proj <- cr$project

    anchored_ids <- disp |>
      filter(
        project ==
          proj,
        is_patient_reference,
        time_factor ==
          cr$required_anchor
      ) |>
      distinct(
        patient_id
      ) |>
      pull(
        patient_id
      )

    pair <- disp |>
      filter(
        project ==
          proj,
        patient_id %in%
          anchored_ids,
        time_factor %in%
          c(
            cr$early_label,
            cr$late_label
          )
      )

    valid_ids <- pair |>
      distinct(
        patient_id,
        time_factor
      ) |>
      count(
        patient_id,
        name =
          "n_times"
      ) |>
      filter(
        n_times == 2
      ) |>
      pull(
        patient_id
      )

    pair <- pair |>
      filter(
        patient_id %in%
          valid_ids
      )

    n_pair <- n_distinct(
      pair$patient_id
    )

    expected_n <- expected_pairs |>
      filter(
        project ==
          proj
      ) |>
      pull(
        n_expected
      )

    if (
      n_pair !=
        expected_n
    ) {
      stop(
        proj,
        ": pair-count guard failed. Expected ",
        expected_n,
        "; observed ",
        n_pair,
        "."
      )
    }

    for (
      metric in c(
        "bray_from_patient_baseline",
        "Shannon",
        "Simpson",
        "Observed_ASV"
      )
    ) {

      paired_rows[[
        length(
          paired_rows
        ) + 1
      ]] <- paired_summary(
        pair =
          pair,
        metric =
          metric,
        early_label =
          cr$early_label,
        late_label =
          cr$late_label,
        project =
          proj
      ) |>
        mutate(
          required_anchor =
            cr$required_anchor,
          early_label =
            cr$early_label,
          late_label =
            cr$late_label,
          .after =
            project
        )
    }

    loo <- loo_bray(
      pair =
        pair,
      early_label =
        cr$early_label,
      late_label =
        cr$late_label,
      project =
        proj
    )

    loo_rows[[
      length(
        loo_rows
      ) + 1
    ]] <- loo

    loo_summary_rows[[
      length(
        loo_summary_rows
      ) + 1
    ]] <- tibble(
      project =
        proj,
      n_leave_one_out_runs =
        nrow(loo),
      direction_preserved_fraction =
        mean(
          loo$direction_preserved
        ),
      min_LOO_mean_delta =
        min(
          loo$mean_delta
        ),
      max_LOO_mean_delta =
        max(
          loo$mean_delta
        ),
      min_LOO_median_delta =
        min(
          loo$median_delta
        ),
      max_LOO_median_delta =
        max(
          loo$median_delta
        ),
      min_LOO_wilcoxon_p =
        min(
          loo$wilcoxon_p,
          na.rm = TRUE
        ),
      max_LOO_wilcoxon_p =
        max(
          loo$wilcoxon_p,
          na.rm = TRUE
        ),
      fraction_LOO_wilcoxon_p_lt_005 =
        mean(
          loo$wilcoxon_p <
            0.05,
          na.rm = TRUE
        )
    )
  }

  paired_summary_all <- bind_rows(
    paired_rows
  )

  loo_all <- bind_rows(
    loo_rows
  )

  loo_summary <- bind_rows(
    loo_summary_rows
  )

  write_csv_safe(
    paired_summary_all,
    file.path(
      DIR_BRAY,
      "V2_STEP90A_common_anchor_paired_metric_robustness.csv"
    )
  )

  write_csv_safe(
    loo_all,
    file.path(
      DIR_BRAY,
      "V2_STEP90A_Bray_leave_one_patient_out_detail.csv"
    )
  )

  write_csv_safe(
    loo_summary,
    file.path(
      DIR_BRAY,
      "V2_STEP90A_Bray_leave_one_patient_out_summary.csv"
    )
  )

  # ----------------------------------------------------------
  # 2. Shared-taxonomy feature-space coherence
  # ----------------------------------------------------------

  shared_rows <- list()
  shared_taxa_rows <- list()

  for (
    rank in c(
      "GENUS",
      "FAMILY"
    )
  ) {

    matrices <- list()

    for (
      proj in natural_projects
    ) {

      p <- file.path(
        STEP89B,
        "02_PATIENT_LEVEL_COHERENCE",
        paste0(
          proj,
          "_",
          rank,
          "_patient_CLR_delta_matrix.csv"
        )
      )

      if (
        !file.exists(p)
      ) {
        stop(
          "Missing Step89B patient delta matrix: ",
          p
        )
      }

      matrices[[proj]] <-
        safe_csv(p)
    }

    taxa_sets <- lapply(
      matrices,
      function(x)
        setdiff(
          names(x),
          "patient_id"
        )
    )

    shared_taxa <- Reduce(
      intersect,
      taxa_sets
    )

    if (
      length(
        shared_taxa
      ) < 5
    ) {
      stop(
        rank,
        ": fewer than 5 shared taxa across the four natural-history cohorts."
      )
    }

    shared_taxa_rows[[
      length(
        shared_taxa_rows
      ) + 1
    ]] <- tibble(
      tax_rank =
        rank,
      taxon =
        shared_taxa
    )

    for (
      j in seq_along(
        natural_projects
      )
    ) {

      proj <- natural_projects[j]

      d <- matrices[[proj]]

      mat <- as.matrix(
        d[
          ,
          shared_taxa,
          drop = FALSE
        ]
      )

      rownames(mat) <- d$patient_id

      met <- cosine_metrics(mat)

      ci <- bootstrap_cosine(
        mat,
        B = 3000,
        seed =
          20260820 + j +
          ifelse(
            rank == "GENUS",
            0,
            100
          )
      )

      shared_rows[[
        length(
          shared_rows
        ) + 1
      ]] <- tibble(
        project =
          proj,
        tax_rank =
          rank,
        n_patients =
          nrow(mat),
        n_shared_taxa =
          ncol(mat),
        coherence_ratio =
          met[
            "coherence_ratio"
          ],
        coherence_ci_low =
          ci$coherence_ci_low,
        coherence_ci_high =
          ci$coherence_ci_high,
        median_pairwise_cosine =
          met[
            "median_pairwise_cosine"
          ],
        median_cosine_ci_low =
          ci$median_cosine_ci_low,
        median_cosine_ci_high =
          ci$median_cosine_ci_high,
        fraction_pairwise_cosine_positive =
          met[
            "fraction_pairwise_cosine_positive"
          ]
      )
    }
  }

  shared_coherence <- bind_rows(
    shared_rows
  )

  shared_taxa_table <- bind_rows(
    shared_taxa_rows
  )

  write_csv_safe(
    shared_coherence,
    file.path(
      DIR_SHARED,
      "V2_STEP90A_shared_taxonomy_patient_coherence.csv"
    )
  )

  write_csv_safe(
    shared_taxa_table,
    file.path(
      DIR_SHARED,
      "V2_STEP90A_taxa_shared_across_four_natural_history_cohorts.csv"
    )
  )

  # ----------------------------------------------------------
  # 3. Core-sepsis taxon breadth / leverage audit
  # ----------------------------------------------------------

  candidates <- core |>
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
      taxon,
      fdr,
      wilcoxon_p,
      paired_effect_dz,
      direction
    ) |>
    distinct()

  genus_delta_path <- file.path(
    STEP89B,
    "02_PATIENT_LEVEL_COHERENCE",
    "PRJNA691455_GENUS_patient_CLR_delta_matrix.csv"
  )

  family_delta_path <- file.path(
    STEP89B,
    "02_PATIENT_LEVEL_COHERENCE",
    "PRJNA691455_FAMILY_patient_CLR_delta_matrix.csv"
  )

  genus_delta <- safe_csv(
    genus_delta_path
  )

  family_delta <- safe_csv(
    family_delta_path
  )

  genus_breadth <- candidate_breadth(
    genus_delta,
    candidates,
    "Genus"
  )

  family_breadth <- candidate_breadth(
    family_delta,
    candidates,
    "Family"
  )

  breadth <- bind_rows(
    genus_breadth,
    family_breadth
  ) |>
    left_join(
      candidates,
      by = c(
        "tax_rank",
        "taxon"
      ),
      suffix = c(
        "_robustness",
        "_step89A2"
      )
    ) |>
    mutate(
      broad_patient_support_flag =
        case_when(
          loo_mean_direction_preserved_fraction ==
            1 &
            max_single_patient_absolute_contribution <
              0.35 &
            (
              positive_fraction_non_ties >=
                0.70 |
              positive_fraction_non_ties <=
                0.30
            ) ~
            "BROAD_DIRECTIONAL_SUPPORT",

          loo_mean_direction_preserved_fraction ==
            1 ~
            "DIRECTION_ROBUST_BUT_HETEROGENEOUS",

          TRUE ~
            "POTENTIALLY_LEVERAGE_SENSITIVE"
        )
    )

  write_csv_safe(
    breadth,
    file.path(
      DIR_TAXON,
      "V2_STEP90A_CORE_SEPSIS_candidate_patient_breadth_and_leverage.csv"
    )
  )

  # ----------------------------------------------------------
  # 4. Integrated interpretation state
  # ----------------------------------------------------------

  sepsis_bray <- paired_summary_all |>
    filter(
      project ==
        "PRJNA691455",
      metric ==
        "bray_from_patient_baseline"
    )

  sepsis_loo <- loo_summary |>
    filter(
      project ==
        "PRJNA691455"
    )

  sepsis_shared_genus <- shared_coherence |>
    filter(
      project ==
        "PRJNA691455",
      tax_rank ==
        "GENUS"
    )

  interpretation_state <- case_when(
    sepsis_bray$mean_delta >
      0 &
      sepsis_loo$direction_preserved_fraction ==
        1 &
      sepsis_shared_genus$median_pairwise_cosine <
        0.20 ~
      "ROBUST_PATIENT_LEVEL_ECOLOGICAL_DISPLACEMENT_WITH_TAXONOMIC_ROUTE_HETEROGENEITY",

    TRUE ~
      "MIXED_ROBUSTNESS_STATE_REQUIRES_REVIEW"
  )

  synthesis <- tibble(
    interpretation_state =
      interpretation_state,
    sepsis_bray_n_pairs =
      sepsis_bray$n_pairs,
    sepsis_bray_mean_delta =
      sepsis_bray$mean_delta,
    sepsis_bray_mean_delta_ci_low =
      sepsis_bray$mean_delta_ci_low,
    sepsis_bray_mean_delta_ci_high =
      sepsis_bray$mean_delta_ci_high,
    sepsis_bray_positive_fraction =
      sepsis_bray$positive_fraction_non_ties,
    sepsis_bray_sign_test_p =
      sepsis_bray$sign_test_p,
    sepsis_bray_wilcoxon_p =
      sepsis_bray$wilcoxon_p,
    sepsis_bray_LOO_direction_preserved_fraction =
      sepsis_loo$direction_preserved_fraction,
    sepsis_shared_genus_n_taxa =
      sepsis_shared_genus$n_shared_taxa,
    sepsis_shared_genus_median_pairwise_cosine =
      sepsis_shared_genus$median_pairwise_cosine,
    sepsis_shared_genus_median_cosine_ci_low =
      sepsis_shared_genus$median_cosine_ci_low,
    sepsis_shared_genus_median_cosine_ci_high =
      sepsis_shared_genus$median_cosine_ci_high
  )

  write_csv_safe(
    synthesis,
    file.path(
      DIR_SYN,
      "V2_STEP90A_integrated_interpretation_state.csv"
    )
  )

  readme <- c(
    "SEPSIS V2 - STEP90A ECOLOGICAL ROBUSTNESS AND SHARED TAXONOMY",
    paste0(
      "Created: ",
      Sys.time()
    ),
    "",
    "BRAY ROBUSTNESS",
    "Common-anchor paired early-vs-late changes were stress-tested with bootstrap confidence intervals, exact sign tests, and leave-one-patient-out analysis.",
    "",
    "SHARED TAXONOMY COHERENCE",
    "Patient CLR delta-vector coherence was recalculated using only Genus/Family taxa present in all four natural-history cohorts.",
    "This reduces the chance that apparent coherence differences are caused merely by different taxonomic feature counts.",
    "",
    "CORE SEPSIS CANDIDATE BREADTH",
    "Nominal/FDR candidate taxa were audited for patient-level sign consistency, single-patient leverage, and leave-one-out direction stability.",
    "",
    "INTERPRETATION STATE",
    interpretation_state,
    "",
    "NEXT",
    "If the ecological displacement remains robust while shared-feature taxonomic coherence remains low, freeze the principal V2 biological narrative around ecological instability rather than a universal differential-taxon signature.",
    "Sequence-identity verification of unusual environmental/thermophilic candidates should be completed before any taxon-specific mechanistic claim."
  )

  writeLines(
    readme,
    file.path(
      OUT,
      "README_STEP90A.txt"
    ),
    useBytes = TRUE
  )

  writeLines(
    c(
      paste0(
        "Completed: ",
        Sys.time()
      ),
      "Status: STEP90A COMPLETE",
      "Four natural-history cohorts underwent common-anchor paired robustness analysis.",
      "Bray leave-one-patient-out analysis completed.",
      "Shared Genus/Family coherence analysis completed.",
      "Core sepsis candidate breadth/leverage audit completed.",
      paste0(
        "Interpretation state: ",
        interpretation_state
      )
    ),
    file.path(
      OUT,
      "_STEP90A_COMPLETE.ok"
    ),
    useBytes = TRUE
  )

  ck(
    "STEP90A COMPLETE"
  )

  cat(
    "\n============================================================\n"
  )
  cat(
    "SEPSIS V2 - STEP90A COMPLETE\n"
  )
  cat(
    "============================================================\n\n"
  )

  cat(
    "Paired robustness:\n"
  )

  print(
    paired_summary_all,
    n = Inf,
    width = Inf
  )

  cat(
    "\nLOO summary:\n"
  )

  print(
    loo_summary,
    n = Inf,
    width = Inf
  )

  cat(
    "\nShared-taxonomy coherence:\n"
  )

  print(
    shared_coherence,
    n = Inf,
    width = Inf
  )

  cat(
    "\nCore candidate breadth:\n"
  )

  print(
    breadth,
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
        "STEP90A FATAL ERROR: ",
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
      "STEP90A FAILED"
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
