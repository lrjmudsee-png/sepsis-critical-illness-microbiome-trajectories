# ============================================================
# Sepsis V2 - Step 89A2
# FIXED paired taxonomic trajectory screen (Genus + Family)
#
# Scientific design:
# - Uses the Step88B pre-specified early-vs-late FOLLOW-UP contrast.
# - Primary analysis = common-anchor paired patients.
# - Secondary analysis = all paired patients.
# - Collapses ASVs to Genus and Family WITHIN each cohort.
# - Does NOT merge raw ASVs across cohorts.
# - Uses CLR paired differences for compositional inference.
# - Relative abundance is exported for interpretation only.
# - BH-FDR is applied within project x rank x subset.
#
# Main natural-history contrast:
#   PRJNA691455  Day 7 vs Day 3, Day 1 common anchor
#   PRJNA851469  Day-7 vs Day-3, Day-1 common anchor
#   PRJNA516701  DAY_7 vs DAY_3, DAY_1 common anchor
#   PRJNA578267  third vs second, first-visit common anchor
#
# Supportive contrasts:
#   PRJNA1166732 Week 2 vs Week 1
#   PRJNA430161  Day 28 vs Day 7
#   PRJEB82425    Infection_D5 vs Infection_D1
#
# IMPORTANT:
# - PRJEB82425 common-anchor n is small and remains supportive.
# - Species-level inference is intentionally NOT performed here.
# ============================================================

options(stringsAsFactors = FALSE)

pkgs <- c(
  "readr",
  "dplyr",
  "tidyr",
  "tibble",
  "purrr",
  "stringr",
  "ggplot2"
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
  library(ggplot2)
})

set.seed(20260820)

# ------------------------------------------------------------
# Paths
# ------------------------------------------------------------

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

OUT <- file.path(
  ROOT,
  "results",
  "V2_29A2_STEP89A_TAXONOMIC_PAIRED_TRAJECTORIES_FIXED"
)

DIR_QC <- file.path(OUT, "01_QC")
DIR_GENUS <- file.path(OUT, "02_GENUS")
DIR_FAMILY <- file.path(OUT, "03_FAMILY")
DIR_SYN <- file.path(OUT, "04_CROSS_COHORT_SYNTHESIS")
DIR_PLOT <- file.path(OUT, "05_EXPLORATORY_PLOTS")

for (d in c(
  OUT,
  DIR_QC,
  DIR_GENUS,
  DIR_FAMILY,
  DIR_SYN,
  DIR_PLOT
)) {
  dir.create(
    d,
    recursive = TRUE,
    showWarnings = FALSE
  )
}

LOG <- file.path(
  OUT,
  "_STEP89A2_runtime_checkpoints.txt"
)

ERR <- file.path(
  OUT,
  "_STEP89A2_FATAL_ERROR.txt"
)

if (file.exists(ERR)) {
  unlink(ERR)
}

ck <- function(x) {
  cat(
    paste0(x, ": ", Sys.time(), "\n"),
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

write_csv_safe <- function(x, p) {
  write_excel_csv(
    x,
    p,
    na = ""
  )
}

bh_adjust <- function(p) {

  out <- rep(
    NA_real_,
    length(p)
  )

  ok <- !is.na(p)

  if (any(ok)) {
    out[ok] <- p.adjust(
      p[ok],
      method = "BH"
    )
  }

  out
}

clean_chr <- function(x) {

  y <- as.character(x)

  y[
    is.na(y) |
    trimws(y) == "" |
    tolower(trimws(y)) %in%
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
        "na",
        "nan",
        "null",
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

  if (
    ncol(counts) !=
      nrow(tax)
  ) {
    stop(
      rank,
      ": counts/taxonomy dimensions do not match."
    )
  }

  labels <- make_tax_labels(
    tax,
    rank
  )

  keep <- !is.na(labels)

  if (
    sum(keep) == 0
  ) {
    stop(
      rank,
      ": no classified taxa remain."
    )
  }

  collapsed <- t(
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
      reorder =
        FALSE
    )
  )

  storage.mode(
    collapsed
  ) <- "numeric"

  collapsed
}

clr_matrix <- function(
  counts,
  pseudocount = 0.5
) {

  if (
    any(
      counts < 0
    )
  ) {
    stop(
      "Negative counts encountered before CLR."
    )
  }

  logx <- log(
    counts +
      pseudocount
  )

  clr <- logx -
    rowMeans(
      logx
    )

  clr
}

relative_abundance <- function(
  counts
) {

  denom <- rowSums(
    counts
  )

  if (
    any(
      denom <= 0
    )
  ) {
    stop(
      "Zero library after taxonomic collapse."
    )
  }

  counts /
    denom
}

paired_boot_ci <- function(
  diff,
  stat = c(
    "mean",
    "median"
  ),
  B = 2000,
  seed = 20260820
) {

  stat <- match.arg(
    stat
  )

  x <- as.numeric(
    diff[
      !is.na(diff)
    ]
  )

  n <- length(x)

  if (
    n < 3
  ) {
    return(
      c(
        low =
          NA_real_,
        high =
          NA_real_
      )
    )
  }

  set.seed(seed)

  vals <- replicate(
    B,
    {

      z <- sample(
        x,
        size = n,
        replace = TRUE
      )

      if (
        stat ==
          "mean"
      ) {
        mean(z)
      } else {
        median(z)
      }
    }
  )

  q <- quantile(
    vals,
    probs = c(
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

safe_wilcox <- function(
  x
) {

  x <- as.numeric(
    x[
      !is.na(x)
    ]
  )

  if (!length(x)) {
    return(
      NA_real_
    )
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
        paired = FALSE,
        exact = FALSE,
        correct = FALSE
      )$p.value
    ),
    error = function(e)
      NA_real_
  )
}

build_pair_ids <- function(
  md,
  required_anchor,
  early_label,
  late_label,
  subset_type
) {

  md2 <- md |>
    filter(
      !is.na(
        analysis_time_order
      )
    )

  patient_n_time <- md2 |>
    distinct(
      patient_id,
      timepoint_key
    ) |>
    count(
      patient_id,
      name = "n_time"
    )

  eligible <- patient_n_time |>
    filter(
      n_time >= 2
    ) |>
    pull(
      patient_id
    )

  md2 <- md2 |>
    filter(
      patient_id %in%
        eligible
    )

  md2$time_order_numeric <- as.numeric(
    md2$analysis_time_order
  )

  md2$time_label <- resolve_time_labels(
    md2
  )

  reference <- md2 |>
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

  md2 <- md2 |>
    left_join(
      reference,
      by = "patient_id"
    )

  if (
    subset_type ==
      "COMMON_ANCHOR_SENSITIVITY"
  ) {

    md2 <- md2 |>
      filter(
        reference_time_label ==
          required_anchor
      )
  }

  pairs <- md2 |>
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
      name = "n_contrast_times"
    ) |>
    filter(
      n_contrast_times ==
        2
    ) |>
    pull(
      patient_id
    )

  list(
    metadata = md2 |>
      filter(
        patient_id %in%
          pairs,
        time_label %in%
          c(
            early_label,
            late_label
          )
      ),
    pair_ids = pairs
  )
}

taxon_stats <- function(
  project,
  role_group,
  rank,
  subset_type,
  required_anchor,
  early_label,
  late_label,
  counts_collapsed,
  md_pair,
  min_prevalence = 0.20,
  min_mean_ra = 0.001
) {

  if (
    nrow(md_pair) == 0
  ) {
    return(
      tibble()
    )
  }

  if (
    anyDuplicated(
      md_pair$run_id
    )
  ) {
    stop(
      project,
      ": duplicated run_id in paired metadata."
    )
  }

  x <- counts_collapsed[
    md_pair$run_id,
    ,
    drop = FALSE
  ]

  if (
    !identical(
      rownames(x),
      as.character(
        md_pair$run_id
      )
    )
  ) {
    stop(
      project,
      ": paired count/metadata alignment failed."
    )
  }

  ra <- relative_abundance(
    x
  )

  clr <- clr_matrix(
    x,
    pseudocount = 0.5
  )

  early_rows <- which(
    md_pair$time_label ==
      early_label
  )

  late_rows <- which(
    md_pair$time_label ==
      late_label
  )

  prevalence_early <- colMeans(
    x[
      early_rows,
      ,
      drop = FALSE
    ] > 0
  )

  prevalence_late <- colMeans(
    x[
      late_rows,
      ,
      drop = FALSE
    ] > 0
  )

  mean_ra_early <- colMeans(
    ra[
      early_rows,
      ,
      drop = FALSE
    ]
  )

  mean_ra_late <- colMeans(
    ra[
      late_rows,
      ,
      drop = FALSE
    ]
  )

  keep_taxa <- names(
    which(
      pmax(
        prevalence_early,
        prevalence_late
      ) >=
        min_prevalence &
      pmax(
        mean_ra_early,
        mean_ra_late
      ) >=
        min_mean_ra
    )
  )

  if (!length(keep_taxa)) {
    return(
      tibble()
    )
  }

  rows <- vector(
    "list",
    length(
      keep_taxa
    )
  )

  for (
    j in seq_along(
      keep_taxa
    )
  ) {

    taxon <- keep_taxa[j]

    dd <- tibble(
      patient_id =
        md_pair$patient_id,
      time_label =
        md_pair$time_label,
      clr =
        clr[
          ,
          taxon
        ],
      relative_abundance =
        ra[
          ,
          taxon
        ]
    )

    wide_clr <- dd |>
      select(
        patient_id,
        time_label,
        clr
      ) |>
      pivot_wider(
        names_from =
          time_label,
        values_from =
          clr
      )

    wide_ra <- dd |>
      select(
        patient_id,
        time_label,
        relative_abundance
      ) |>
      pivot_wider(
        names_from =
          time_label,
        values_from =
          relative_abundance
      )

    if (
      !(early_label %in%
          names(wide_clr)) ||
      !(late_label %in%
          names(wide_clr))
    ) {
      next
    }

    wide_clr <- wide_clr |>
      filter(
        !is.na(
          .data[[early_label]]
        ),
        !is.na(
          .data[[late_label]]
        )
      )

    wide_ra <- wide_ra |>
      filter(
        patient_id %in%
          wide_clr$patient_id
      )

    diff <- wide_clr[[late_label]] -
      wide_clr[[early_label]]

    n <- length(diff)

    if (!n) {
      next
    }

    sd_diff <- if (
      n >= 2
    ) {
      sd(
        diff,
        na.rm = TRUE
      )
    } else {
      NA_real_
    }

    dz <- if (
      !is.na(
        sd_diff
      ) &&
      sd_diff > 0
    ) {
      mean(
        diff,
        na.rm = TRUE
      ) /
        sd_diff
    } else {
      NA_real_
    }

    ci_mean <- paired_boot_ci(
      diff,
      stat = "mean",
      seed =
        20260820 + j
    )

    ci_median <- paired_boot_ci(
      diff,
      stat = "median",
      seed =
        20261820 + j
    )

    early_ra <- wide_ra[[early_label]]

    late_ra <- wide_ra[[late_label]]

    rows[[j]] <- tibble(
      project = project,
      analysis_role_group =
        role_group,
      tax_rank =
        rank,
      subset_type =
        subset_type,
      required_anchor =
        required_anchor,
      early_label =
        early_label,
      late_label =
        late_label,
      taxon =
        taxon,
      n_pairs =
        n,
      prevalence_early =
        prevalence_early[
          taxon
        ],
      prevalence_late =
        prevalence_late[
          taxon
        ],
      mean_ra_early =
        mean_ra_early[
          taxon
        ],
      mean_ra_late =
        mean_ra_late[
          taxon
        ],
      median_ra_early =
        median(
          early_ra,
          na.rm = TRUE
        ),
      median_ra_late =
        median(
          late_ra,
          na.rm = TRUE
        ),
      mean_clr_difference =
        mean(
          diff,
          na.rm = TRUE
        ),
      median_clr_difference =
        median(
          diff,
          na.rm = TRUE
        ),
      mean_clr_diff_ci_low =
        ci_mean[
          "low"
        ],
      mean_clr_diff_ci_high =
        ci_mean[
          "high"
        ],
      median_clr_diff_ci_low =
        ci_median[
          "low"
        ],
      median_clr_diff_ci_high =
        ci_median[
          "high"
        ],
      paired_effect_dz =
        dz,
      wilcoxon_p =
        safe_wilcox(
          diff
        ),
      direction = case_when(
        mean(
          diff,
          na.rm = TRUE
        ) > 0 ~
          "INCREASE",
        mean(
          diff,
          na.rm = TRUE
        ) < 0 ~
          "DECREASE",
        TRUE ~
          "NO_CHANGE"
      )
    )
  }

  bind_rows(
    rows
  ) |>
    mutate(
      fdr =
        bh_adjust(
          wilcoxon_p
        ),
      evidence_flag =
        case_when(
          fdr < 0.05 ~
            "FDR_LT_0.05",
          fdr < 0.10 ~
            "FDR_LT_0.10",
          wilcoxon_p < 0.05 ~
            "NOMINAL_ONLY",
          TRUE ~
            "NO_STATISTICAL_SIGNAL"
        )
    ) |>
    arrange(
      fdr,
      wilcoxon_p,
      desc(
        abs(
          paired_effect_dz
        )
      )
    )
}

make_cross_cohort_synthesis <- function(
  primary_results,
  rank
) {

  x <- primary_results |>
    filter(
      tax_rank ==
        rank,
      subset_type ==
        "COMMON_ANCHOR_SENSITIVITY",
      project %in%
        c(
          "PRJNA691455",
          "PRJNA851469",
          "PRJNA516701",
          "PRJNA578267"
        )
    ) |>
    select(
      project,
      taxon,
      paired_effect_dz,
      mean_clr_difference,
      fdr,
      direction
    )

  if (!nrow(x)) {
    return(
      tibble()
    )
  }

  w <- x |>
    pivot_wider(
      names_from =
        project,
      values_from =
        c(
          paired_effect_dz,
          mean_clr_difference,
          fdr,
          direction
        ),
      names_sep =
        "__"
    )

  get_col <- function(
    prefix,
    project
  ) {
    paste0(
      prefix,
      "__",
      project
    )
  }

  direction_value <- function(
    data,
    project
  ) {

    nm <- get_col(
      "direction",
      project
    )

    if (
      nm %in%
        names(data)
    ) {
      data[[nm]]
    } else {
      rep(
        NA_character_,
        nrow(data)
      )
    }
  }

  fdr_value <- function(
    data,
    project
  ) {

    nm <- get_col(
      "fdr",
      project
    )

    if (
      nm %in%
        names(data)
    ) {
      data[[nm]]
    } else {
      rep(
        NA_real_,
        nrow(data)
      )
    }
  }

  d_sepsis <- direction_value(
    w,
    "PRJNA691455"
  )

  d_icu1 <- direction_value(
    w,
    "PRJNA851469"
  )

  d_icu2 <- direction_value(
    w,
    "PRJNA516701"
  )

  d_non <- direction_value(
    w,
    "PRJNA578267"
  )

  q_sepsis <- fdr_value(
    w,
    "PRJNA691455"
  )

  q_icu1 <- fdr_value(
    w,
    "PRJNA851469"
  )

  q_icu2 <- fdr_value(
    w,
    "PRJNA516701"
  )

  q_non <- fdr_value(
    w,
    "PRJNA578267"
  )

  same_sepsis_icu1 <-
    !is.na(d_sepsis) &
    !is.na(d_icu1) &
    d_sepsis ==
      d_icu1

  same_sepsis_icu2 <-
    !is.na(d_sepsis) &
    !is.na(d_icu2) &
    d_sepsis ==
      d_icu2

  all_critical_same <-
    !is.na(d_sepsis) &
    !is.na(d_icu1) &
    !is.na(d_icu2) &
    d_sepsis ==
      d_icu1 &
    d_sepsis ==
      d_icu2

  recovery_opposite <-
    all_critical_same &
    !is.na(d_non) &
    d_non !=
      d_sepsis

  n_q10_critical <-
    rowSums(
      cbind(
        !is.na(q_sepsis) &
          q_sepsis < 0.10,
        !is.na(q_icu1) &
          q_icu1 < 0.10,
        !is.na(q_icu2) &
          q_icu2 < 0.10
      )
    )

  w |>
    mutate(
      tax_rank =
        rank,
      n_critical_cohorts_q_lt_0.10 =
        n_q10_critical,
      all_critical_same_direction =
        all_critical_same,
      nonsepsis_opposite_direction =
        recovery_opposite,
      synthesis_pattern =
        case_when(
          recovery_opposite &
            n_q10_critical >= 2 ~
            "CRITICAL_ILLNESS_CONSISTENT_NONSEPSIS_OPPOSITE_STRONG",

          recovery_opposite ~
            "CRITICAL_ILLNESS_CONSISTENT_NONSEPSIS_OPPOSITE_DIRECTIONAL",

          all_critical_same &
            n_q10_critical >= 2 ~
            "CRITICAL_ILLNESS_CONSISTENT_STRONG",

          all_critical_same ~
            "CRITICAL_ILLNESS_CONSISTENT_DIRECTIONAL",

          !is.na(d_sepsis) &
            (
              same_sepsis_icu1 |
              same_sepsis_icu2
            ) ~
            "SEPSIS_PLUS_ONE_ICU_DIRECTIONAL_SUPPORT",

          !is.na(q_sepsis) &
            q_sepsis < 0.10 ~
            "SEPSIS_ASSOCIATED_CANDIDATE",

          TRUE ~
            "MIXED_OR_WEAK"
        )
    ) |>
    arrange(
      desc(
        n_critical_cohorts_q_lt_0.10
      ),
      synthesis_pattern,
      taxon
    )
}

make_effect_heatmap <- function(
  synthesis,
  rank,
  out_path
) {

  if (!nrow(synthesis)) {
    return(
      invisible(NULL)
    )
  }

  effect_cols <- grep(
    "^paired_effect_dz__",
    names(synthesis),
    value = TRUE
  )

  if (!length(effect_cols)) {
    return(
      invisible(NULL)
    )
  }

  candidates <- synthesis |>
    filter(
      synthesis_pattern !=
        "MIXED_OR_WEAK"
    ) |>
    slice_head(
      n = 30
    )

  if (!nrow(candidates)) {
    return(
      invisible(NULL)
    )
  }

  long <- candidates |>
    select(
      taxon,
      all_of(
        effect_cols
      )
    ) |>
    pivot_longer(
      cols =
        all_of(
          effect_cols
        ),
      names_to =
        "project",
      values_to =
        "paired_effect_dz"
    ) |>
    mutate(
      project =
        str_remove(
          project,
          "^paired_effect_dz__"
        )
    )

  g <- ggplot(
    long,
    aes(
      x = project,
      y = taxon,
      fill =
        paired_effect_dz
    )
  ) +
    geom_tile() +
    geom_text(
      aes(
        label =
          ifelse(
            is.na(
              paired_effect_dz
            ),
            "",
            sprintf(
              "%.2f",
              paired_effect_dz
            )
          )
      ),
      size = 2.5
    ) +
    labs(
      title = paste0(
        rank,
        " paired CLR trajectory effect sizes"
      ),
      x = NULL,
      y = NULL,
      fill = "paired dz"
    ) +
    theme_bw(
      base_size = 10
    ) +
    theme(
      axis.text.x =
        element_text(
          angle = 35,
          hjust = 1
        )
    )

  ggsave(
    out_path,
    g,
    width = 9,
    height = max(
      5,
      0.27 *
        n_distinct(
          long$taxon
        ) +
        2
    )
  )
}

main <- function() {

  ck(
    "STEP89A2 STARTED"
  )

  step87_complete <- file.path(
    STEP87B,
    "_STEP87B_COMPLETE.ok"
  )

  object_registry_path <- file.path(
    STEP87B,
    "V2_STEP87B_analysis_object_registry.csv"
  )

  step88b_complete <- file.path(
    STEP88B,
    "_STEP88B_COMPLETE.ok"
  )

  paired88b_path <- file.path(
    STEP88B,
    "02_PAIRED_CONTRASTS",
    "V2_STEP88B_paired_early_vs_late_followup_contrasts.csv"
  )

  anchor88b_path <- file.path(
    STEP88B,
    "01_ANCHOR_AUDIT",
    "V2_STEP88B_anchor_audit.csv"
  )

  required <- c(
    step87_complete,
    object_registry_path,
    step88b_complete,
    paired88b_path,
    anchor88b_path
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
        "Required Step87B/Step88B files are missing:\n",
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
    object_registry_path
  )

  paired88b <- safe_csv(
    paired88b_path
  )

  anchor88b <- safe_csv(
    anchor88b_path
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

  contrast_registry <- paired88b |>
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

  expected_projects <- c(
    "PRJNA691455",
    "PRJNA851469",
    "PRJNA516701",
    "PRJNA578267",
    "PRJNA1166732",
    "PRJNA430161",
    "PRJEB82425"
  )

  if (
    !setequal(
      contrast_registry$project,
      expected_projects
    )
  ) {
    stop(
      "Step88B contrast registry differs from expected seven cohorts."
    )
  }

  expected_common_pairs <- tribble(
    ~project,
    ~n_expected,

    "PRJNA691455",
    9L,

    "PRJNA851469",
    14L,

    "PRJNA516701",
    14L,

    "PRJNA578267",
    32L,

    "PRJNA1166732",
    40L,

    "PRJNA430161",
    8L,

    "PRJEB82425",
    6L
  )

  ck(
    "INPUT STATE GUARDED"
  )

  result_rows <- list()
  qc_rows <- list()

  for (
    i in seq_len(
      nrow(
        contrast_registry
      )
    )
  ) {

    cr <- contrast_registry[
      i,
      ,
      drop = FALSE
    ]

    proj <- cr$project

    ck(
      paste0(
        proj,
        " START"
      )
    )

    rr <- registry |>
      filter(
        project ==
          proj
      )

    if (
      nrow(rr) != 1
    ) {
      stop(
        proj,
        ": Step87B registry row missing/non-unique."
      )
    }

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

    if (
      !is.matrix(
        counts
      ) ||
      nrow(counts) !=
        nrow(md) ||
      ncol(counts) !=
        nrow(tax) ||
      !identical(
        rownames(
          counts
        ),
        as.character(
          md$run_id
        )
      )
    ) {
      stop(
        proj,
        ": Step87B object alignment guard failed."
      )
    }

    for (
      subset_type in c(
        "COMMON_ANCHOR_SENSITIVITY",
        "ALL_PAIRED"
      )
    ) {

      pair_info <- build_pair_ids(
        md =
          md,
        required_anchor =
          cr$required_anchor,
        early_label =
          cr$early_label,
        late_label =
          cr$late_label,
        subset_type =
          subset_type
      )

      md_pair <- pair_info$metadata

      n_pairs <- n_distinct(
        md_pair$patient_id
      )

      if (
        subset_type ==
          "COMMON_ANCHOR_SENSITIVITY"
      ) {

        expected_n <- expected_common_pairs |>
          filter(
            project ==
              proj
          ) |>
          pull(
            n_expected
          )

        if (
          length(expected_n) != 1 ||
          n_pairs !=
            expected_n
        ) {
          stop(
            paste0(
              proj,
              ": common-anchor pair-count guard failed. Expected ",
              expected_n,
              "; observed ",
              n_pairs,
              "."
            )
          )
        }
      }

      for (
        rank in c(
          "Genus",
          "Family"
        )
      ) {

        collapsed <- collapse_counts(
          counts =
            counts,
          tax =
            tax,
          rank =
            rank
        )

        stats <- taxon_stats(
          project =
            proj,
          role_group =
            cr$analysis_role_group,
          rank =
            rank,
          subset_type =
            subset_type,
          required_anchor =
            cr$required_anchor,
          early_label =
            cr$early_label,
          late_label =
            cr$late_label,
          counts_collapsed =
            collapsed,
          md_pair =
            md_pair,
          min_prevalence =
            0.20,
          min_mean_ra =
            0.001
        )

        result_rows[[
          length(
            result_rows
          ) + 1
        ]] <- stats

        qc_rows[[
          length(
            qc_rows
          ) + 1
        ]] <- tibble(
          project =
            proj,
          analysis_role_group =
            cr$analysis_role_group,
          tax_rank =
            rank,
          subset_type =
            subset_type,
          required_anchor =
            cr$required_anchor,
          early_label =
            cr$early_label,
          late_label =
            cr$late_label,
          n_pairs =
            n_pairs,
          n_collapsed_taxa =
            ncol(
              collapsed
            ),
          n_tested_taxa =
            nrow(
              stats
            ),
          n_fdr_lt_005 =
            sum(
              stats$fdr <
                0.05,
              na.rm = TRUE
            ),
          n_fdr_lt_010 =
            sum(
              stats$fdr <
                0.10,
              na.rm = TRUE
            )
        )
      }
    }

    ck(
      paste0(
        proj,
        " COMPLETE"
      )
    )

    rm(
      obj,
      counts,
      tax,
      md
    )

    gc(
      verbose = FALSE
    )
  }

  results <- bind_rows(
    result_rows
  )

  qc <- bind_rows(
    qc_rows
  )

  if (
    nrow(qc) !=
      7 * 2 * 2
  ) {
    stop(
      "Step89A2 QC row-count guard failed."
    )
  }

  if (
    !all(
      expected_common_pairs$project %in%
        qc$project
    )
  ) {
    stop(
      "Not all expected projects were processed."
    )
  }

  # ----------------------------------------------------------
  # Export project/rank results
  # ----------------------------------------------------------

  genus <- results |>
    filter(
      tax_rank ==
        "Genus"
    )

  family <- results |>
    filter(
      tax_rank ==
        "Family"
    )

  write_csv_safe(
    qc,
    file.path(
      DIR_QC,
      "V2_STEP89A2_taxonomic_analysis_QC.csv"
    )
  )

  write_csv_safe(
    genus,
    file.path(
      DIR_GENUS,
      "V2_STEP89A2_ALL_genus_paired_CLR_results.csv"
    )
  )

  write_csv_safe(
    family,
    file.path(
      DIR_FAMILY,
      "V2_STEP89A2_ALL_family_paired_CLR_results.csv"
    )
  )

  for (
    proj in expected_projects
  ) {

    write_csv_safe(
      genus |>
        filter(
          project ==
            proj
        ),
      file.path(
        DIR_GENUS,
        paste0(
          proj,
          "_genus_paired_CLR_results.csv"
        )
      )
    )

    write_csv_safe(
      family |>
        filter(
          project ==
            proj
        ),
      file.path(
        DIR_FAMILY,
        paste0(
          proj,
          "_family_paired_CLR_results.csv"
        )
      )
    )
  }

  # ----------------------------------------------------------
  # Primary common-anchor candidate tables
  # ----------------------------------------------------------

  primary_candidates <- results |>
    filter(
      subset_type ==
        "COMMON_ANCHOR_SENSITIVITY"
    ) |>
    mutate(
      priority =
        case_when(
          fdr < 0.05 ~
            "HIGH",
          fdr < 0.10 ~
            "MODERATE",
          wilcoxon_p < 0.05 ~
            "EXPLORATORY",
          TRUE ~
            "LOW"
        )
    )

  write_csv_safe(
    primary_candidates,
    file.path(
      DIR_SYN,
      "V2_STEP89A2_PRIMARY_common_anchor_all_taxa.csv"
    )
  )

  # ----------------------------------------------------------
  # Cross-cohort natural-history synthesis
  # ----------------------------------------------------------

  genus_syn <- make_cross_cohort_synthesis(
    primary_candidates,
    "Genus"
  )

  family_syn <- make_cross_cohort_synthesis(
    primary_candidates,
    "Family"
  )

  write_csv_safe(
    genus_syn,
    file.path(
      DIR_SYN,
      "V2_STEP89A2_GENUS_cross_cohort_trajectory_synthesis.csv"
    )
  )

  write_csv_safe(
    family_syn,
    file.path(
      DIR_SYN,
      "V2_STEP89A2_FAMILY_cross_cohort_trajectory_synthesis.csv"
    )
  )

  make_effect_heatmap(
    genus_syn,
    "Genus",
    file.path(
      DIR_PLOT,
      "V2_STEP89A2_GENUS_effect_size_heatmap.pdf"
    )
  )

  make_effect_heatmap(
    family_syn,
    "Family",
    file.path(
      DIR_PLOT,
      "V2_STEP89A2_FAMILY_effect_size_heatmap.pdf"
    )
  )

  # ----------------------------------------------------------
  # Core sepsis shortlist
  # ----------------------------------------------------------

  core_sepsis <- primary_candidates |>
    filter(
      project ==
        "PRJNA691455"
    ) |>
    arrange(
      tax_rank,
      fdr,
      wilcoxon_p,
      desc(
        abs(
          paired_effect_dz
        )
      )
    )

  write_csv_safe(
    core_sepsis,
    file.path(
      DIR_SYN,
      "V2_STEP89A2_CORE_SEPSIS_taxonomic_shortlist.csv"
    )
  )

  # ----------------------------------------------------------
  # README + complete
  # ----------------------------------------------------------

  core_genus_q05 <- core_sepsis |>
    filter(
      tax_rank ==
        "Genus",
      fdr <
        0.05
    ) |>
    nrow()

  core_genus_q10 <- core_sepsis |>
    filter(
      tax_rank ==
        "Genus",
      fdr <
        0.10
    ) |>
    nrow()

  core_family_q05 <- core_sepsis |>
    filter(
      tax_rank ==
        "Family",
      fdr <
        0.05
    ) |>
    nrow()

  readme <- c(
    "SEPSIS V2 - STEP89A2 TAXONOMIC PAIRED TRAJECTORIES",
    paste0(
      "Created: ",
      Sys.time()
    ),
    "",
    "PRIMARY ANALYSIS",
    "Common-anchor paired early-vs-late follow-up contrast within each cohort.",
    "Genus and Family are analyzed separately.",
    "CLR paired differences are used for compositional inference.",
    "Relative abundance summaries are descriptive only.",
    "",
    "FILTER",
    "Taxa are tested if prevalence >=20% in either contrast timepoint AND mean relative abundance >=0.1% in either timepoint.",
    "CLR is calculated using all classified taxa at the relevant rank before the testing filter.",
    "",
    "MULTIPLE TESTING",
    "BH-FDR is controlled within project x taxonomic rank x subset.",
    "q<0.05 = high-priority signal; q<0.10 = moderate-priority candidate.",
    "",
    "CROSS-COHORT SYNTHESIS",
    "Only Genus/Family labels are harmonized across cohorts.",
    "Raw ASVs are never merged.",
    "Cross-cohort synthesis emphasizes direction and paired standardized effect size, not raw abundance magnitude.",
    "",
    "CORE SEPSIS",
    paste0(
      "Genus q<0.05: ",
      core_genus_q05
    ),
    paste0(
      "Genus q<0.10: ",
      core_genus_q10
    ),
    paste0(
      "Family q<0.05: ",
      core_family_q05
    ),
    "",
    "INTERPRETATION",
    "Taxa shared in direction across sepsis and ICU-background cohorts represent critical-illness-associated trajectory candidates.",
    "Taxa showing the opposite direction in the non-sepsis surgical control are prioritized as destabilization-versus-recovery candidates.",
    "PRJEB82425 remains supportive because its common Inclusion anchor subset is small.",
    "",
    "NEXT",
    "Audit Step89A2 candidates. If coherent, proceed to focused longitudinal models/visualization for prioritized taxa and manuscript-level cross-cohort evidence synthesis."
  )

  writeLines(
    readme,
    file.path(
      OUT,
      "README_STEP89A2.txt"
    ),
    useBytes = TRUE
  )

  writeLines(
    c(
      paste0(
        "Completed: ",
        Sys.time()
      ),
      "Status: STEP89A2 COMPLETE",
      "Seven longitudinal cohorts analyzed.",
      "Genus and Family paired CLR analyses completed.",
      "Common-anchor analysis treated as primary.",
      "All-paired analysis retained as sensitivity.",
      paste0(
        "Core sepsis Genus q<0.05: ",
        core_genus_q05
      ),
      paste0(
        "Core sepsis Genus q<0.10: ",
        core_genus_q10
      ),
      paste0(
        "Core sepsis Family q<0.05: ",
        core_family_q05
      ),
      "Cross-cohort natural-history synthesis generated."
    ),
    file.path(
      OUT,
      "_STEP89A2_COMPLETE.ok"
    ),
    useBytes = TRUE
  )

  ck(
    "STEP89A2 COMPLETE"
  )

  cat(
    "\n============================================================\n"
  )
  cat(
    "SEPSIS V2 - STEP89A2 COMPLETE\n"
  )
  cat(
    "============================================================\n\n"
  )

  cat(
    "QC summary:\n"
  )

  print(
    qc,
    n = Inf,
    width = Inf
  )

  cat(
    "\nCore sepsis top signals:\n"
  )

  print(
    core_sepsis |>
      select(
        tax_rank,
        taxon,
        n_pairs,
        mean_clr_difference,
        paired_effect_dz,
        wilcoxon_p,
        fdr,
        direction,
        evidence_flag
      ) |>
      slice_head(
        n = 30
      ),
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
        "STEP89A2 FATAL ERROR: ",
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
      "STEP89A2 FAILED"
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
