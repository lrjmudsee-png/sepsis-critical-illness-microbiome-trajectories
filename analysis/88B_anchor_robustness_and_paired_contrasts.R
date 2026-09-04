# ============================================================
# Sepsis V2 - Step 88B
# Anchor audit + paired trajectory contrasts + robustness layer
#
# INPUT:
#   Step88A2 completed outputs
#
# PURPOSE:
#   1) Audit whether each cohort uses a common patient baseline anchor.
#   2) Define one pre-specified early-vs-late FOLLOW-UP contrast per cohort.
#   3) Compute paired within-patient changes for:
#        - Bray displacement from patient baseline
#        - Shannon
#        - Simpson
#        - Observed ASV
#   4) Repeat contrasts in a common-anchor sensitivity subset.
#   5) Add BH-FDR to the primary categorical mixed-model global tests.
#   6) Produce a direction/robustness registry for later cross-cohort synthesis.
#
# IMPORTANT:
# - No cross-cohort ASV merge.
# - Raw Bray magnitudes are NOT pooled across cohorts.
# - Cross-cohort interpretation is based on within-cohort direction,
#   paired standardized effect size, and consistency.
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

OUT <- file.path(
  ROOT,
  "results",
  "V2_28B_STEP88B_ANCHOR_ROBUSTNESS_AND_PAIRED_CONTRASTS"
)

DIR_AUDIT <- file.path(OUT, "01_ANCHOR_AUDIT")
DIR_PAIR <- file.path(OUT, "02_PAIRED_CONTRASTS")
DIR_MODEL <- file.path(OUT, "03_MODEL_FDR")
DIR_SYN <- file.path(OUT, "04_SYNTHESIS")

for (d in c(
  OUT,
  DIR_AUDIT,
  DIR_PAIR,
  DIR_MODEL,
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
  "_STEP88B_runtime_checkpoints.txt"
)

ERR <- file.path(
  OUT,
  "_STEP88B_FATAL_ERROR.txt"
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

paired_boot_ci <- function(
  diff,
  stat = c("mean", "median"),
  B = 2000,
  seed = 20260820
) {

  stat <- match.arg(stat)

  x <- as.numeric(
    diff[
      !is.na(diff)
    ]
  )

  n <- length(x)

  if (n < 3) {
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
        size = n,
        replace = TRUE
      )

      if (stat == "mean") {
        mean(z)
      } else {
        median(z)
      }
    }
  )

  as.numeric(
    quantile(
      vals,
      probs = c(
        0.025,
        0.975
      ),
      na.rm = TRUE,
      names = FALSE
    )
  ) |>
    setNames(
      c(
        "low",
        "high"
      )
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

  if (all(abs(x) < 1e-15)) {
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
    error = function(e) NA_real_
  )
}

paired_effect <- function(
  d,
  project,
  early_label,
  late_label,
  metric,
  subset_type = "ALL_PAIRED",
  required_anchor = NA_character_
) {

  dd <- d |>
    filter(
      .data$project == .env$project
    )

  if (
    subset_type ==
      "COMMON_ANCHOR_SENSITIVITY"
  ) {

    if (
      is.na(required_anchor) ||
      required_anchor == ""
    ) {
      stop(
        project,
        ": common-anchor sensitivity requested without an anchor."
      )
    }

    anchored_ids <- dd |>
      filter(
        is_patient_reference,
        time_factor ==
          required_anchor
      ) |>
      distinct(
        patient_id
      ) |>
      pull(
        patient_id
      )

    dd <- dd |>
      filter(
        patient_id %in%
          anchored_ids
      )
  }

  pair <- dd |>
    filter(
      time_factor %in%
        c(
          early_label,
          late_label
        )
    ) |>
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
    !(early_label %in% names(pair)) ||
    !(late_label %in% names(pair))
  ) {
    return(
      tibble(
        project = project,
        metric = metric,
        subset_type =
          subset_type,
        required_anchor =
          required_anchor,
        early_label =
          early_label,
        late_label =
          late_label,
        n_pairs = 0L,
        early_median =
          NA_real_,
        late_median =
          NA_real_,
        median_difference =
          NA_real_,
        mean_difference =
          NA_real_,
        mean_diff_ci_low =
          NA_real_,
        mean_diff_ci_high =
          NA_real_,
        median_diff_ci_low =
          NA_real_,
        median_diff_ci_high =
          NA_real_,
        paired_sd_difference =
          NA_real_,
        paired_effect_dz =
          NA_real_,
        wilcoxon_p =
          NA_real_,
        direction =
          "NO_DATA"
      )
    )
  }

  pair <- pair |>
    filter(
      !is.na(
        .data[[early_label]]
      ),
      !is.na(
        .data[[late_label]]
      )
    )

  diff <- pair[[late_label]] -
    pair[[early_label]]

  n <- length(diff)

  if (!n) {
    return(
      tibble(
        project = project,
        metric = metric,
        subset_type =
          subset_type,
        required_anchor =
          required_anchor,
        early_label =
          early_label,
        late_label =
          late_label,
        n_pairs = 0L,
        early_median =
          NA_real_,
        late_median =
          NA_real_,
        median_difference =
          NA_real_,
        mean_difference =
          NA_real_,
        mean_diff_ci_low =
          NA_real_,
        mean_diff_ci_high =
          NA_real_,
        median_diff_ci_low =
          NA_real_,
        median_diff_ci_high =
          NA_real_,
        paired_sd_difference =
          NA_real_,
        paired_effect_dz =
          NA_real_,
        wilcoxon_p =
          NA_real_,
        direction =
          "NO_DATA"
      )
    )
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
    !is.na(sd_diff) &&
    sd_diff > 0
  ) {
    mean(
      diff,
      na.rm = TRUE
    ) / sd_diff
  } else {
    NA_real_
  }

  ci_mean <- paired_boot_ci(
    diff,
    stat = "mean"
  )

  ci_med <- paired_boot_ci(
    diff,
    stat = "median"
  )

  mean_diff <- mean(
    diff,
    na.rm = TRUE
  )

  direction <- case_when(
    mean_diff > 0 ~
      "INCREASE",
    mean_diff < 0 ~
      "DECREASE",
    TRUE ~
      "NO_CHANGE"
  )

  tibble(
    project = project,
    metric = metric,
    subset_type =
      subset_type,
    required_anchor =
      required_anchor,
    early_label =
      early_label,
    late_label =
      late_label,
    n_pairs =
      n,
    early_median =
      median(
        pair[[early_label]],
        na.rm = TRUE
      ),
    late_median =
      median(
        pair[[late_label]],
        na.rm = TRUE
      ),
    median_difference =
      median(
        diff,
        na.rm = TRUE
      ),
    mean_difference =
      mean_diff,
    mean_diff_ci_low =
      ci_mean["low"],
    mean_diff_ci_high =
      ci_mean["high"],
    median_diff_ci_low =
      ci_med["low"],
    median_diff_ci_high =
      ci_med["high"],
    paired_sd_difference =
      sd_diff,
    paired_effect_dz =
      dz,
    wilcoxon_p =
      safe_wilcox(diff),
    direction =
      direction
  )
}

main <- function() {

  ck("STEP88B STARTED")

  complete_path <- file.path(
    STEP88A2,
    "_STEP88A2_COMPLETE.ok"
  )

  disp_path <- file.path(
    STEP88A2,
    "03_BETA_DISPLACEMENT",
    "V2_STEP88A2_ALL_within_patient_bray_displacement.csv"
  )

  model_global_path <- file.path(
    STEP88A2,
    "05_MODELS",
    "V2_STEP88A2_model_global_time_tests.csv"
  )

  model_status_path <- file.path(
    STEP88A2,
    "05_MODELS",
    "V2_STEP88A2_model_status.csv"
  )

  qc_path <- file.path(
    STEP88A2,
    "01_QC",
    "V2_STEP88A2_cohort_QC_summary.csv"
  )

  required <- c(
    complete_path,
    disp_path,
    model_global_path,
    model_status_path,
    qc_path
  )

  if (!all(file.exists(required))) {
    stop(
      paste0(
        "Required Step88A2 files are missing:\n",
        paste(
          required[
            !file.exists(required)
          ],
          collapse = "\n"
        )
      )
    )
  }

  disp <- safe_csv(
    disp_path
  )

  model_global <- safe_csv(
    model_global_path
  )

  model_status <- safe_csv(
    model_status_path
  )

  qc <- safe_csv(
    qc_path
  )

  if (
    nrow(qc) != 8 ||
    sum(qc$primary_samples) != 785 ||
    sum(qc$eligible_GE2_samples) != 653 ||
    sum(qc$eligible_GE2_patients) != 258
  ) {
    stop(
      "Step88A2 QC state differs from expected frozen state."
    )
  }

  if (
    nrow(model_status) != 48 ||
    any(
      model_status$status !=
        "OK"
    )
  ) {
    stop(
      "Step88A2 model status is not fully OK."
    )
  }

  ck("STEP88A2 STATE GUARDED")

  # ----------------------------------------------------------
  # Pre-specified contrast registry
  # ----------------------------------------------------------

  contrasts <- tribble(
    ~project,
    ~analysis_role_group,
    ~required_anchor,
    ~early_followup,
    ~late_followup,

    "PRJNA691455",
    "CORE_SEPSIS",
    "Day 1",
    "Day 3",
    "Day 7",

    "PRJNA851469",
    "ICU_BACKGROUND",
    "Day-1",
    "Day-3",
    "Day-7",

    "PRJNA516701",
    "ICU_BACKGROUND",
    "DAY_1",
    "DAY_3",
    "DAY_7",

    "PRJNA578267",
    "NONSEPSIS_SURGICAL_CONTROL",
    "Time first",
    "Time second",
    "Time third",

    "PRJNA1166732",
    "INTERVENTION_SUPPORT",
    "Baseline",
    "Week 1",
    "Week 2",

    "PRJNA430161",
    "INTERVENTION_SUPPORT",
    "Baseline",
    "Day 7",
    "Day 28",

    "PRJEB82425",
    "ICU_INFECTION_EXTERNAL",
    "Inclusion",
    "Infection_D1",
    "Infection_D5"
  )

  # ----------------------------------------------------------
  # Anchor audit
  # ----------------------------------------------------------

  refs <- disp |>
    filter(
      is_patient_reference
    ) |>
    count(
      project,
      time_factor,
      name =
        "n_reference_patients"
    )

  project_ref_total <- refs |>
    group_by(
      project
    ) |>
    summarise(
      n_reference_patients_total =
        sum(
          n_reference_patients
        ),
      dominant_reference_n =
        max(
          n_reference_patients
        ),
      dominant_reference_stage =
        time_factor[
          which.max(
            n_reference_patients
          )
        ],
      .groups = "drop"
    ) |>
    mutate(
      dominant_reference_fraction =
        dominant_reference_n /
        n_reference_patients_total,
      anchor_structure =
        case_when(
          dominant_reference_fraction ==
            1 ~
            "STRICT_COMMON_ANCHOR",
          dominant_reference_fraction >=
            0.90 ~
            "MOSTLY_COMMON_ANCHOR",
          TRUE ~
            "HETEROGENEOUS_ANCHOR"
        )
    )

  anchor_audit <- contrasts |>
    left_join(
      project_ref_total,
      by = "project"
    ) |>
    left_join(
      refs |>
        rename(
          required_anchor =
            time_factor,
          n_required_anchor_patients =
            n_reference_patients
        ),
      by = c(
        "project",
        "required_anchor"
      )
    ) |>
    mutate(
      n_required_anchor_patients =
        replace_na(
          n_required_anchor_patients,
          0L
        ),
      required_anchor_fraction =
        n_required_anchor_patients /
        n_reference_patients_total,
      required_anchor_flag =
        case_when(
          required_anchor_fraction ==
            1 ~
            "STRICT",
          required_anchor_fraction >=
            0.90 ~
            "MOSTLY",
          required_anchor_fraction >=
            0.50 ~
            "PARTIAL",
          TRUE ~
            "LIMITED"
        )
    )

  write_csv_safe(
    refs,
    file.path(
      DIR_AUDIT,
      "V2_STEP88B_reference_stage_counts.csv"
    )
  )

  write_csv_safe(
    anchor_audit,
    file.path(
      DIR_AUDIT,
      "V2_STEP88B_anchor_audit.csv"
    )
  )

  # ----------------------------------------------------------
  # Paired contrasts
  # ----------------------------------------------------------

  metrics <- c(
    "bray_from_patient_baseline",
    "Shannon",
    "Simpson",
    "Observed_ASV"
  )

  paired_rows <- list()

  for (
    i in seq_len(
      nrow(
        contrasts
      )
    )
  ) {

    rr <- contrasts[
      i,
      ,
      drop = FALSE
    ]

    for (
      metric in metrics
    ) {

      paired_rows[[
        length(paired_rows) + 1
      ]] <- paired_effect(
        d = disp,
        project =
          rr$project,
        early_label =
          rr$early_followup,
        late_label =
          rr$late_followup,
        metric =
          metric,
        subset_type =
          "ALL_PAIRED",
        required_anchor =
          rr$required_anchor
      ) |>
        mutate(
          analysis_role_group =
            rr$analysis_role_group,
          .after = project
        )

      paired_rows[[
        length(paired_rows) + 1
      ]] <- paired_effect(
        d = disp,
        project =
          rr$project,
        early_label =
          rr$early_followup,
        late_label =
          rr$late_followup,
        metric =
          metric,
        subset_type =
          "COMMON_ANCHOR_SENSITIVITY",
        required_anchor =
          rr$required_anchor
      ) |>
        mutate(
          analysis_role_group =
            rr$analysis_role_group,
          .after = project
        )
    }
  }

  paired <- bind_rows(
    paired_rows
  ) |>
    group_by(
      subset_type,
      metric
    ) |>
    mutate(
      wilcoxon_fdr_within_metric =
        bh_adjust(
          wilcoxon_p
        )
    ) |>
    ungroup()

  write_csv_safe(
    paired,
    file.path(
      DIR_PAIR,
      "V2_STEP88B_paired_early_vs_late_followup_contrasts.csv"
    )
  )

  # ----------------------------------------------------------
  # Primary global model FDR
  # Categorical-time model is the primary flexible model.
  # Continuous-day model remains a secondary trend analysis.
  # ----------------------------------------------------------

  primary_global <- model_global |>
    filter(
      model_type %in%
        c(
          "CATEGORICAL_TIME_MIXED_MODEL",
          "FOLLOWUP_CATEGORICAL_TIME_MIXED_MODEL"
        )
    ) |>
    group_by(
      outcome
    ) |>
    mutate(
      global_time_fdr_within_outcome =
        bh_adjust(
          global_time_p
        )
    ) |>
    ungroup()

  secondary_trend <- model_global |>
    filter(
      model_type %in%
        c(
          "CONTINUOUS_DAY_MIXED_MODEL",
          "FOLLOWUP_CONTINUOUS_DAY_MIXED_MODEL"
        )
    ) |>
    group_by(
      outcome
    ) |>
    mutate(
      global_time_fdr_within_outcome =
        bh_adjust(
          global_time_p
        )
    ) |>
    ungroup()

  write_csv_safe(
    primary_global,
    file.path(
      DIR_MODEL,
      "V2_STEP88B_PRIMARY_categorical_global_models_with_FDR.csv"
    )
  )

  write_csv_safe(
    secondary_trend,
    file.path(
      DIR_MODEL,
      "V2_STEP88B_SECONDARY_continuous_day_models_with_FDR.csv"
    )
  )

  # ----------------------------------------------------------
  # Bray synthesis registry
  # ----------------------------------------------------------

  bray_all <- paired |>
    filter(
      metric ==
        "bray_from_patient_baseline",
      subset_type ==
        "ALL_PAIRED"
    ) |>
    select(
      project,
      analysis_role_group,
      n_pairs_all =
        n_pairs,
      bray_mean_diff_all =
        mean_difference,
      bray_dz_all =
        paired_effect_dz,
      bray_p_all =
        wilcoxon_p,
      bray_fdr_all =
        wilcoxon_fdr_within_metric,
      direction_all =
        direction
    )

  bray_anchor <- paired |>
    filter(
      metric ==
        "bray_from_patient_baseline",
      subset_type ==
        "COMMON_ANCHOR_SENSITIVITY"
    ) |>
    select(
      project,
      n_pairs_common_anchor =
        n_pairs,
      bray_mean_diff_common_anchor =
        mean_difference,
      bray_dz_common_anchor =
        paired_effect_dz,
      bray_p_common_anchor =
        wilcoxon_p,
      bray_fdr_common_anchor =
        wilcoxon_fdr_within_metric,
      direction_common_anchor =
        direction
    )

  synthesis <- contrasts |>
    select(
      project,
      analysis_role_group
    ) |>
    left_join(
      anchor_audit |>
        select(
          project,
          required_anchor,
          required_anchor_fraction,
          required_anchor_flag
        ),
      by = "project"
    ) |>
    left_join(
      bray_all,
      by = c(
        "project",
        "analysis_role_group"
      )
    ) |>
    left_join(
      bray_anchor,
      by = "project"
    ) |>
    mutate(
      direction_robust =
        direction_all ==
          direction_common_anchor,
      formal_interpretation_tier =
        case_when(
          required_anchor_flag %in%
            c(
              "STRICT",
              "MOSTLY"
            ) &
            n_pairs_common_anchor >=
              8 ~
              "FORMAL_TRAJECTORY_EVIDENCE",

          n_pairs_common_anchor >=
            5 ~
              "SUPPORTIVE_SENSITIVITY_ONLY",

          TRUE ~
              "DESCRIPTIVE_ONLY"
        )
    )

  write_csv_safe(
    synthesis,
    file.path(
      DIR_SYN,
      "V2_STEP88B_Bray_direction_and_robustness_registry.csv"
    )
  )

  # Hard guards for the scientifically important planned contrasts.
  expected_pairs <- tribble(
    ~project,
    ~n_all,
    ~n_anchor,

    "PRJNA691455",
    9L,
    9L,

    "PRJNA851469",
    14L,
    14L,

    "PRJNA516701",
    15L,
    14L,

    "PRJNA578267",
    35L,
    32L,

    "PRJNA1166732",
    40L,
    40L,

    "PRJNA430161",
    8L,
    8L,

    "PRJEB82425",
    22L,
    6L
  )

  observed_pairs <- paired |>
    filter(
      metric ==
        "bray_from_patient_baseline"
    ) |>
    select(
      project,
      subset_type,
      n_pairs
    ) |>
    pivot_wider(
      names_from =
        subset_type,
      values_from =
        n_pairs
    ) |>
    rename(
      n_all =
        ALL_PAIRED,
      n_anchor =
        COMMON_ANCHOR_SENSITIVITY
    )

  guard <- expected_pairs |>
    left_join(
      observed_pairs,
      by = "project",
      suffix = c(
        "_expected",
        "_observed"
      )
    )

  if (
    any(
      guard$n_all_expected !=
        guard$n_all_observed
    ) ||
    any(
      guard$n_anchor_expected !=
        guard$n_anchor_observed
    )
  ) {
    write_csv_safe(
      guard,
      file.path(
        OUT,
        "V2_STEP88B_FAILED_pair_count_guard.csv"
      )
    )

    stop(
      "Planned paired contrast counts differ from the audited Step88A2 state."
    )
  }

  readme <- c(
    "SEPSIS V2 - STEP88B ANCHOR ROBUSTNESS AND PAIRED CONTRASTS",
    paste0(
      "Created: ",
      Sys.time()
    ),
    "",
    "WHY THIS STEP EXISTS",
    "Step88A2 established cohort-specific time effects but some cohorts use heterogeneous earliest-sample anchors.",
    "This step formalizes early-vs-late follow-up paired contrasts and repeats them in patients sharing the intended common anchor.",
    "",
    "PRIMARY INTERPRETATION RULE",
    "Raw Bray-Curtis values are not compared or pooled directly across cohorts because ASV feature spaces and 16S regions differ.",
    "Cross-cohort evidence is interpreted using within-cohort direction, paired standardized effect size (dz), FDR, anchor robustness, and role of each cohort.",
    "",
    "CORE SEPSIS CONTRAST",
    "PRJNA691455: Day 7 versus Day 3, with Day 1 as common anchor.",
    "",
    "ICU BACKGROUND CONTRASTS",
    "PRJNA851469: Day 7 versus Day 3, Day 1 anchor.",
    "PRJNA516701: Day 7 versus Day 3, Day 1 anchor sensitivity.",
    "",
    "NON-SEPSIS CONTROL CONTRAST",
    "PRJNA578267: third versus second visit, first-visit anchor sensitivity.",
    "",
    "EXTERNAL INFECTION COHORT",
    "PRJEB82425 has heterogeneous patient reference stages. The Inclusion-anchored subset is therefore supportive rather than definitive.",
    "",
    "NEXT",
    "If Step88B confirms robust direction patterns, proceed to cohort-specific taxonomic longitudinal trajectories and cross-cohort taxonomic effect synthesis."
  )

  writeLines(
    readme,
    file.path(
      OUT,
      "README_STEP88B.txt"
    ),
    useBytes = TRUE
  )

  writeLines(
    c(
      paste0(
        "Completed: ",
        Sys.time()
      ),
      "Status: STEP88B COMPLETE",
      "Step88A2 state audited: 785 primary observations / 653 longitudinal samples.",
      "Seven planned cohort contrasts completed.",
      "Common-anchor sensitivity completed.",
      "Primary categorical mixed-model global tests BH-FDR adjusted.",
      "No cross-cohort raw Bray pooling performed."
    ),
    file.path(
      OUT,
      "_STEP88B_COMPLETE.ok"
    ),
    useBytes = TRUE
  )

  ck("STEP88B COMPLETE")

  cat(
    "\n============================================================\n"
  )
  cat(
    "SEPSIS V2 - STEP88B COMPLETE\n"
  )
  cat(
    "============================================================\n\n"
  )

  cat(
    "Anchor audit:\n"
  )

  print(
    anchor_audit,
    n = Inf,
    width = Inf
  )

  cat(
    "\nBray robustness registry:\n"
  )

  print(
    synthesis,
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
        "STEP88B FATAL ERROR: ",
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

    ck("STEP88B FAILED")

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
