# ============================================================
# Sepsis V2 - Step 90B
# Corrected shared-taxonomy coherence inference
#
# WHY:
# Step90A point estimates are valid, but bootstrap CIs for
# pairwise cosine are not suitable because patient resampling
# with replacement can duplicate the same patient and create
# artificial cosine=1 "pairs".
#
# THIS STEP REPLACES THOSE CIs WITH:
# 1) patient-level leave-one-out (LOO) ranges
# 2) whole-vector sign-flip permutation tests
#
# The sign-flip null preserves each patient's taxonomic delta
# magnitude and internal covariance, but destroys a common
# cohort-level direction.
#
# INPUT:
# Step89B patient CLR delta matrices.
# No upstream microbiome processing is rerun.
# ============================================================

options(stringsAsFactors = FALSE)

pkgs <- c(
  "readr",
  "dplyr",
  "tidyr",
  "tibble",
  "purrr"
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
})

set.seed(20260820)

ROOT <- "E:/sepsis_project"

STEP89B <- file.path(
  ROOT,
  "results",
  "V2_29B_STEP89B_TAXONOMIC_ROBUSTNESS_AND_PROVENANCE"
)

STEP90A <- file.path(
  ROOT,
  "results",
  "V2_30A_STEP90A_ECOLOGICAL_ROBUSTNESS_AND_SHARED_TAXONOMY"
)

OUT <- file.path(
  ROOT,
  "results",
  "V2_30B_STEP90B_CORRECTED_COHERENCE_INFERENCE"
)

dir.create(
  OUT,
  recursive = TRUE,
  showWarnings = FALSE
)

LOG <- file.path(
  OUT,
  "_STEP90B_runtime_checkpoints.txt"
)

ERR <- file.path(
  OUT,
  "_STEP90B_FATAL_ERROR.txt"
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
        mean_pairwise_cosine =
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
        mean_pairwise_cosine =
          NA_real_,
        median_pairwise_cosine =
          NA_real_,
        fraction_pairwise_cosine_positive =
          NA_real_
      )
    )
  }

  unit <- x / norms

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
    mean_pairwise_cosine =
      mean(
        vals,
        na.rm = TRUE
      ),
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

pairwise_cosine_detail <- function(
  mat,
  project,
  rank
) {

  x <- as.matrix(mat)

  ids <- rownames(x)

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

  ids <- ids[keep]

  norms <- norms[keep]

  unit <- x / norms

  sim <- unit %*%
    t(unit)

  idx <- which(
    upper.tri(sim),
    arr.ind = TRUE
  )

  tibble(
    project =
      project,
    tax_rank =
      rank,
    patient_1 =
      ids[
        idx[, 1]
      ],
    patient_2 =
      ids[
        idx[, 2]
      ],
    cosine =
      sim[idx]
  )
}

loo_metrics <- function(
  mat,
  project,
  rank
) {

  x <- as.matrix(mat)

  ids <- rownames(x)

  full <- cosine_metrics(x)

  rows <- vector(
    "list",
    nrow(x)
  )

  for (
    i in seq_len(
      nrow(x)
    )
  ) {

    m <- cosine_metrics(
      x[
        -i,
        ,
        drop = FALSE
      ]
    )

    rows[[i]] <- tibble(
      project =
        project,
      tax_rank =
        rank,
      omitted_patient_id =
        ids[i],
      n_remaining =
        nrow(x) - 1,
      coherence_ratio =
        m[
          "coherence_ratio"
        ],
      mean_pairwise_cosine =
        m[
          "mean_pairwise_cosine"
        ],
      median_pairwise_cosine =
        m[
          "median_pairwise_cosine"
        ],
      fraction_pairwise_cosine_positive =
        m[
          "fraction_pairwise_cosine_positive"
        ],
      coherence_direction_relative_to_full =
        sign(
          m[
            "mean_pairwise_cosine"
          ]
        ) ==
          sign(
            full[
              "mean_pairwise_cosine"
            ]
          )
    )
  }

  bind_rows(rows)
}

all_sign_patterns <- function(n) {

  # One global sign reversal produces the same pairwise cosines
  # and coherence norm as its opposite pattern, so fixing the
  # first patient to +1 removes redundant mirror patterns.
  if (
    n <= 16
  ) {

    k <- n - 1

    grid <- expand.grid(
      rep(
        list(
          c(
            -1,
            1
          )
        ),
        k
      )
    )

    cbind(
      1,
      as.matrix(grid)
    )

  } else {

    NULL
  }
}

sign_flip_test <- function(
  mat,
  B = 20000,
  seed = 20260820
) {

  x <- as.matrix(mat)

  n <- nrow(x)

  obs <- cosine_metrics(x)

  exact_patterns <- all_sign_patterns(n)

  if (
    !is.null(
      exact_patterns
    )
  ) {

    signs <- exact_patterns
    method <- "EXACT_SIGN_FLIP"

  } else {

    set.seed(seed)

    signs <- matrix(
      sample(
        c(
          -1,
          1
        ),
        size =
          B * n,
        replace =
          TRUE
      ),
      nrow = B,
      ncol = n
    )

    signs[, 1] <- 1

    method <-
      "MONTE_CARLO_SIGN_FLIP"
  }

  nperm <- nrow(signs)

  perm_coh <- numeric(nperm)

  perm_mean_cos <- numeric(nperm)

  perm_median_cos <- numeric(nperm)

  for (
    i in seq_len(nperm)
  ) {

    xp <- x *
      signs[
        i,
      ]

    m <- cosine_metrics(xp)

    perm_coh[i] <-
      m[
        "coherence_ratio"
      ]

    perm_mean_cos[i] <-
      m[
        "mean_pairwise_cosine"
      ]

    perm_median_cos[i] <-
      m[
        "median_pairwise_cosine"
      ]
  }

  # One-sided: evidence for greater common-direction coherence
  # than expected under random patient-level vector reversal.
  p_coh <- (
    sum(
      perm_coh >=
        obs[
          "coherence_ratio"
        ]
    ) +
      1
  ) /
    (
      nperm +
        1
    )

  p_mean <- (
    sum(
      perm_mean_cos >=
        obs[
          "mean_pairwise_cosine"
        ]
    ) +
      1
  ) /
    (
      nperm +
        1
    )

  p_median <- (
    sum(
      perm_median_cos >=
        obs[
          "median_pairwise_cosine"
        ]
    ) +
      1
  ) /
    (
      nperm +
        1
    )

  tibble(
    n_patients =
      n,
    n_permutations =
      nperm,
    permutation_method =
      method,
    coherence_ratio =
      obs[
        "coherence_ratio"
      ],
    coherence_signflip_p =
      p_coh,
    mean_pairwise_cosine =
      obs[
        "mean_pairwise_cosine"
      ],
    mean_cosine_signflip_p =
      p_mean,
    median_pairwise_cosine =
      obs[
        "median_pairwise_cosine"
      ],
    median_cosine_signflip_p =
      p_median,
    fraction_pairwise_cosine_positive =
      obs[
        "fraction_pairwise_cosine_positive"
      ]
  )
}

main <- function() {

  ck(
    "STEP90B STARTED"
  )

  step89b_complete <- file.path(
    STEP89B,
    "_STEP89B_COMPLETE.ok"
  )

  step90a_complete <- file.path(
    STEP90A,
    "_STEP90A_COMPLETE.ok"
  )

  shared_taxa_path <- file.path(
    STEP90A,
    "02_SHARED_TAXONOMY_COHERENCE",
    "V2_STEP90A_taxa_shared_across_four_natural_history_cohorts.csv"
  )

  if (
    !all(
      file.exists(
        c(
          step89b_complete,
          step90a_complete,
          shared_taxa_path
        )
      )
    )
  ) {
    stop(
      "Required Step89B/Step90A files are missing."
    )
  }

  shared_taxa_tbl <- safe_csv(
    shared_taxa_path
  )

  projects <- c(
    "PRJNA691455",
    "PRJNA851469",
    "PRJNA516701",
    "PRJNA578267"
  )

  expected_pairs <- c(
    PRJNA691455 = 9,
    PRJNA851469 = 14,
    PRJNA516701 = 14,
    PRJNA578267 = 32
  )

  result_rows <- list()
  loo_rows <- list()
  pair_rows <- list()

  for (
    rank in c(
      "GENUS",
      "FAMILY"
    )
  ) {

    taxa <- shared_taxa_tbl |>
      filter(
        tax_rank ==
          rank
      ) |>
      pull(
        taxon
      )

    expected_taxa <- if (
      rank ==
        "GENUS"
    ) {
      22L
    } else {
      20L
    }

    if (
      length(taxa) !=
        expected_taxa
    ) {
      stop(
        rank,
        ": expected ",
        expected_taxa,
        " shared taxa; observed ",
        length(taxa),
        "."
      )
    }

    for (
      j in seq_along(
        projects
      )
    ) {

      proj <- projects[j]

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
          "Missing patient delta matrix: ",
          p
        )
      }

      d <- safe_csv(p)

      if (
        nrow(d) !=
          expected_pairs[
            proj
          ]
      ) {
        stop(
          proj,
          " ",
          rank,
          ": expected ",
          expected_pairs[
            proj
          ],
          " patients; observed ",
          nrow(d),
          "."
        )
      }

      missing_taxa <- setdiff(
        taxa,
        names(d)
      )

      if (
        length(
          missing_taxa
        )
      ) {
        stop(
          proj,
          " ",
          rank,
          ": shared taxa missing from matrix: ",
          paste(
            missing_taxa,
            collapse = ", "
          )
        )
      }

      mat <- as.matrix(
        d[
          ,
          taxa,
          drop = FALSE
        ]
      )

      rownames(mat) <-
        d$patient_id

      perm <- sign_flip_test(
        mat,
        B = 20000,
        seed =
          20260820 +
          j +
          ifelse(
            rank ==
              "GENUS",
            0,
            100
          )
      ) |>
        mutate(
          project =
            proj,
          tax_rank =
            rank,
          n_shared_taxa =
            ncol(mat),
          .before = 1
        )

      result_rows[[
        length(
          result_rows
        ) + 1
      ]] <- perm

      loo_rows[[
        length(
          loo_rows
        ) + 1
      ]] <- loo_metrics(
        mat,
        project =
          proj,
        rank =
          rank
      )

      pair_rows[[
        length(
          pair_rows
        ) + 1
      ]] <- pairwise_cosine_detail(
        mat,
        project =
          proj,
        rank =
          rank
      )
    }
  }

  results <- bind_rows(
    result_rows
  )

  loo <- bind_rows(
    loo_rows
  )

  pair_detail <- bind_rows(
    pair_rows
  )

  loo_summary <- loo |>
    group_by(
      project,
      tax_rank
    ) |>
    summarise(
      n_LOO =
        n(),
      coherence_ratio_LOO_min =
        min(
          coherence_ratio,
          na.rm = TRUE
        ),
      coherence_ratio_LOO_max =
        max(
          coherence_ratio,
          na.rm = TRUE
        ),
      mean_cosine_LOO_min =
        min(
          mean_pairwise_cosine,
          na.rm = TRUE
        ),
      mean_cosine_LOO_max =
        max(
          mean_pairwise_cosine,
          na.rm = TRUE
        ),
      median_cosine_LOO_min =
        min(
          median_pairwise_cosine,
          na.rm = TRUE
        ),
      median_cosine_LOO_max =
        max(
          median_pairwise_cosine,
          na.rm = TRUE
        ),
      fraction_positive_LOO_min =
        min(
          fraction_pairwise_cosine_positive,
          na.rm = TRUE
        ),
      fraction_positive_LOO_max =
        max(
          fraction_pairwise_cosine_positive,
          na.rm = TRUE
        ),
      .groups =
        "drop"
    )

  results_final <- results |>
    left_join(
      loo_summary,
      by = c(
        "project",
        "tax_rank"
      )
    ) |>
    mutate(
      coherence_interpretation =
        case_when(
          median_cosine_signflip_p <
            0.05 &
            median_pairwise_cosine >
              0.15 ~
            "EVIDENCE_OF_SHARED_TAXONOMIC_DIRECTION",

          median_cosine_signflip_p <
            0.10 &
            median_pairwise_cosine >
              0.10 ~
            "WEAK_TO_MODERATE_SHARED_DIRECTION",

          TRUE ~
            "NO_STRONG_SHARED_TAXONOMIC_DIRECTION"
        )
    )

  write_csv_safe(
    results_final,
    file.path(
      OUT,
      "V2_STEP90B_corrected_shared_taxonomy_coherence_inference.csv"
    )
  )

  write_csv_safe(
    loo,
    file.path(
      OUT,
      "V2_STEP90B_shared_taxonomy_coherence_LOO_detail.csv"
    )
  )

  write_csv_safe(
    pair_detail,
    file.path(
      OUT,
      "V2_STEP90B_shared_taxonomy_pairwise_cosine_detail.csv"
    )
  )

  core_genus <- results_final |>
    filter(
      project ==
        "PRJNA691455",
      tax_rank ==
        "GENUS"
    )

  nonsepsis_genus <- results_final |>
    filter(
      project ==
        "PRJNA578267",
      tax_rank ==
        "GENUS"
    )

  final_state <- case_when(
    core_genus$median_pairwise_cosine <
      0.20 &
      core_genus$median_cosine_signflip_p >=
        0.05 &
      nonsepsis_genus$median_pairwise_cosine >
        core_genus$median_pairwise_cosine ~
      "ECOLOGICAL_DISPLACEMENT_ROBUST_TAXONOMIC_DIRECTION_NOT_UNIVERSALLY_COHERENT",

    TRUE ~
      "COHERENCE_INFERENCE_REQUIRES_REVIEW"
  )

  readme <- c(
    "SEPSIS V2 - STEP90B CORRECTED COHERENCE INFERENCE",
    paste0(
      "Created: ",
      Sys.time()
    ),
    "",
    "CORRECTION",
    "Step90A shared-taxonomy point estimates remain valid.",
    "Step90A bootstrap CIs for pairwise cosine should not be used because resampling patients with replacement can duplicate one patient and create artificial same-patient cosine=1 pairs.",
    "",
    "REPLACEMENT INFERENCE",
    "Step90B uses whole-patient taxonomic delta-vector sign flips to test whether observed cohort coherence exceeds random direction.",
    "For n<=16 patients, the sign-flip test is exact after removing the redundant global mirror pattern.",
    "For larger cohorts, Monte Carlo sign flips are used.",
    "",
    "ROBUSTNESS",
    "Leave-one-patient-out ranges are exported for coherence ratio, mean cosine, median cosine, and fraction of positive patient-pair cosines.",
    "",
    "FINAL STATE",
    final_state,
    "",
    "USAGE",
    "Use Step90B inference rather than Step90A bootstrap coherence CIs in any manuscript or figure.",
    "Step90A Bray paired/LOO results are unaffected and remain valid."
  )

  writeLines(
    readme,
    file.path(
      OUT,
      "README_STEP90B.txt"
    ),
    useBytes = TRUE
  )

  writeLines(
    c(
      paste0(
        "Completed: ",
        Sys.time()
      ),
      "Status: STEP90B COMPLETE",
      "Shared-taxonomy point estimates independently recalculated.",
      "Sign-flip coherence inference completed.",
      "Leave-one-patient-out coherence ranges completed.",
      "Step90A bootstrap pairwise-cosine CIs superseded.",
      paste0(
        "Final state: ",
        final_state
      )
    ),
    file.path(
      OUT,
      "_STEP90B_COMPLETE.ok"
    ),
    useBytes = TRUE
  )

  ck(
    "STEP90B COMPLETE"
  )

  cat(
    "\n============================================================\n"
  )

  cat(
    "SEPSIS V2 - STEP90B COMPLETE\n"
  )

  cat(
    "============================================================\n\n"
  )

  print(
    results_final,
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
        "STEP90B FATAL ERROR: ",
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
      "STEP90B FAILED"
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
