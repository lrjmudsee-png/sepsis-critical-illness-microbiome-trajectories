# ============================================================
# Sepsis V2 - Step 88A
# Cohort-specific longitudinal diversity + within-patient
# Bray-Curtis displacement from each patient's earliest
# longitudinal sample.
#
# Input frozen state from Step87B:
#   785 primary patient-time observations
#   390 project-specific patients
#   258 patients with >=2 timepoints
#   132 patients with >=3 timepoints
#   8 cohort-specific RDS analysis objects
#
# Longitudinal repeated-measures analysis state expected here:
#   653 samples from 258 patients
#
# IMPORTANT:
# - Cohorts are analyzed separately.
# - No cross-cohort ASV merging.
# - Static controls/singletons are NOT forced into time models.
# - PRJNA691455 longitudinal model includes only its 29
#   time-resolved sepsis samples (10 patients).
# - Bray-Curtis is computed on relative abundance.
# - Shannon/Simpson are computed from the count table.
# - Observed richness is retained descriptively and modeled
#   with log10 library size adjustment.
# ============================================================

options(stringsAsFactors = FALSE)

required_pkgs <- c(
  "readr",
  "dplyr",
  "tidyr",
  "purrr",
  "tibble",
  "stringr",
  "vegan",
  "lme4",
  "lmerTest",
  "ggplot2"
)

missing_pkgs <- required_pkgs[
  !vapply(
    required_pkgs,
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
  library(purrr)
  library(tibble)
  library(stringr)
  library(vegan)
  library(lme4)
  library(lmerTest)
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

OUT <- file.path(
  ROOT,
  "results",
  "V2_28A_STEP88A_LONGITUDINAL_DIVERSITY_AND_DISPLACEMENT"
)

DIR_QC <- file.path(OUT, "01_QC")
DIR_ALPHA <- file.path(OUT, "02_ALPHA")
DIR_BETA <- file.path(OUT, "03_BETA_DISPLACEMENT")
DIR_PCOA <- file.path(OUT, "04_PCOA")
DIR_MODEL <- file.path(OUT, "05_MODELS")
DIR_PLOT <- file.path(OUT, "06_EXPLORATORY_PLOTS")

for (d in c(
  OUT,
  DIR_QC,
  DIR_ALPHA,
  DIR_BETA,
  DIR_PCOA,
  DIR_MODEL,
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
  "_STEP88A_runtime_checkpoints.txt"
)

ERR <- file.path(
  OUT,
  "_STEP88A_FATAL_ERROR.txt"
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

# ------------------------------------------------------------
# Expected frozen state
# ------------------------------------------------------------

expected_longitudinal <- tribble(
  ~project,       ~eligible_samples_GE2, ~eligible_patients_GE2,

  "PRJEB82425",                  89L, 35L,
  "PRJNA1166732",               120L, 40L,
  "PRJNA430161",                 26L,  9L,
  "PRJNA516701",                 97L, 41L,
  "PRJNA578267",                192L, 80L,
  "PRJNA691455",                 29L, 10L,
  "PRJNA851469",                100L, 43L,
  "PRJNA978257",                  0L,  0L
)

project_order <- expected_longitudinal$project

# ------------------------------------------------------------
# Helpers
# ------------------------------------------------------------

clean_chr <- function(x) {
  y <- as.character(x)
  y[
    is.na(y) |
    trimws(y) == "" |
    tolower(trimws(y)) %in% c("na", "nan", "null")
  ] <- NA_character_
  y
}

make_time_label <- function(md) {
  tr <- clean_chr(md$time_raw)
  tk <- clean_chr(md$timepoint_key)
  ifelse(
    !is.na(tr),
    tr,
    tk
  )
}

resolve_adjustment <- function(project, d) {

  candidate <- switch(
    project,

    "PRJEB82425" = "phenotype",

    "PRJNA1166732" = "intervention_arm",

    "PRJNA430161" = "intervention_arm",

    NA_character_
  )

  if (is.na(candidate)) {
    return(NA_character_)
  }

  if (!(candidate %in% names(d))) {
    return(NA_character_)
  }

  vals <- clean_chr(d[[candidate]])

  if (length(unique(na.omit(vals))) < 2) {
    return(NA_character_)
  }

  candidate
}

safe_lmer <- function(
  d,
  outcome,
  time_mode = c("factor", "numeric"),
  include_depth = TRUE,
  adjustment = NA_character_
) {

  time_mode <- match.arg(time_mode)

  result_empty <- list(
    status = "NOT_RUN",
    global = tibble(),
    coefficients = tibble()
  )

  needed <- c(
    outcome,
    "patient_id"
  )

  if (time_mode == "factor") {
    needed <- c(needed, "time_factor")
  } else {
    needed <- c(needed, "time_day")
  }

  if (include_depth) {
    needed <- c(
      needed,
      "log10_library_size"
    )
  }

  if (!is.na(adjustment)) {
    needed <- c(
      needed,
      adjustment
    )
  }

  needed <- unique(needed)

  if (!all(needed %in% names(d))) {
    result_empty$status <- paste0(
      "MISSING_COLUMNS: ",
      paste(
        setdiff(
          needed,
          names(d)
        ),
        collapse = ", "
      )
    )
    return(result_empty)
  }

  dd <- d |>
    select(
      all_of(needed)
    ) |>
    drop_na()

  if (
    nrow(dd) < 8 ||
    n_distinct(dd$patient_id) < 4
  ) {
    result_empty$status <- "TOO_FEW_COMPLETE_OBSERVATIONS"
    return(result_empty)
  }

  if (time_mode == "factor") {
    dd$time_factor <- droplevels(
      factor(dd$time_factor)
    )

    if (nlevels(dd$time_factor) < 2) {
      result_empty$status <- "TIME_FACTOR_HAS_LT2_LEVELS"
      return(result_empty)
    }

    time_term <- "time_factor"

  } else {

    if (
      n_distinct(dd$time_day) < 2
    ) {
      result_empty$status <- "TIME_NUMERIC_HAS_LT2_VALUES"
      return(result_empty)
    }

    time_term <- "time_day"
  }

  fixed_base <- character(0)

  if (include_depth) {
    fixed_base <- c(
      fixed_base,
      "log10_library_size"
    )
  }

  if (!is.na(adjustment)) {
    dd[[adjustment]] <- factor(
      dd[[adjustment]]
    )

    if (nlevels(dd[[adjustment]]) >= 2) {
      fixed_base <- c(
        fixed_base,
        adjustment
      )
    }
  }

  full_terms <- c(
    time_term,
    fixed_base
  )

  reduced_terms <- fixed_base

  rhs_full <- paste(
    c(
      full_terms,
      "(1 | patient_id)"
    ),
    collapse = " + "
  )

  rhs_reduced <- paste(
    c(
      if (length(reduced_terms)) {
        reduced_terms
      } else {
        "1"
      },
      "(1 | patient_id)"
    ),
    collapse = " + "
  )

  f_full <- as.formula(
    paste(
      outcome,
      "~",
      rhs_full
    )
  )

  f_reduced <- as.formula(
    paste(
      outcome,
      "~",
      rhs_reduced
    )
  )

  fit_full <- tryCatch(
    lmerTest::lmer(
      f_full,
      data = dd,
      REML = FALSE,
      control = lme4::lmerControl(
        optimizer = "bobyqa",
        optCtrl = list(
          maxfun = 2e5
        )
      )
    ),
    error = function(e) e
  )

  if (inherits(fit_full, "error")) {
    result_empty$status <- paste0(
      "FULL_MODEL_ERROR: ",
      conditionMessage(fit_full)
    )
    return(result_empty)
  }

  fit_reduced <- tryCatch(
    lmerTest::lmer(
      f_reduced,
      data = dd,
      REML = FALSE,
      control = lme4::lmerControl(
        optimizer = "bobyqa",
        optCtrl = list(
          maxfun = 2e5
        )
      )
    ),
    error = function(e) e
  )

  if (inherits(fit_reduced, "error")) {
    result_empty$status <- paste0(
      "REDUCED_MODEL_ERROR: ",
      conditionMessage(fit_reduced)
    )
    return(result_empty)
  }

  lrt <- tryCatch(
    anova(
      fit_reduced,
      fit_full
    ),
    error = function(e) e
  )

  global_p <- NA_real_
  chisq <- NA_real_
  df_diff <- NA_real_

  if (!inherits(lrt, "error")) {

    lrt_df <- as.data.frame(lrt)

    if (nrow(lrt_df) >= 2) {

      if ("Pr(>Chisq)" %in% names(lrt_df)) {
        global_p <- lrt_df[2, "Pr(>Chisq)"]
      }

      if ("Chisq" %in% names(lrt_df)) {
        chisq <- lrt_df[2, "Chisq"]
      }

      if ("Chi Df" %in% names(lrt_df)) {
        df_diff <- lrt_df[2, "Chi Df"]
      }
    }
  }

  cc <- as.data.frame(
    coef(
      summary(fit_full)
    )
  )

  cc <- rownames_to_column(
    cc,
    "term"
  )

  names(cc) <- make.names(
    names(cc)
  )

  p_col <- grep(
    "^Pr",
    names(cc),
    value = TRUE
  )

  coef_out <- tibble(
    term = cc$term,
    estimate = cc$Estimate,
    std_error = cc$Std..Error,
    df = if ("df" %in% names(cc)) cc$df else NA_real_,
    statistic = if ("t.value" %in% names(cc)) cc$t.value else NA_real_,
    p_value = if (length(p_col)) cc[[p_col[1]]] else NA_real_
  )

  singular <- lme4::isSingular(
    fit_full,
    tol = 1e-4
  )

  conv_messages <- fit_full@optinfo$conv$lme4$messages

  if (is.null(conv_messages)) {
    conv_messages <- ""
  } else {
    conv_messages <- paste(
      conv_messages,
      collapse = "; "
    )
  }

  global <- tibble(
    n_observations = nrow(dd),
    n_patients = n_distinct(dd$patient_id),
    time_mode = time_mode,
    outcome = outcome,
    adjustment = ifelse(
      is.na(adjustment),
      "",
      adjustment
    ),
    include_depth = include_depth,
    likelihood_ratio_chisq = chisq,
    df_difference = df_diff,
    global_time_p = global_p,
    singular_fit = singular,
    convergence_message = conv_messages
  )

  list(
    status = "OK",
    global = global,
    coefficients = coef_out
  )
}

make_alpha_plot <- function(
  d,
  project,
  p
) {

  dd <- d |>
    select(
      run_id,
      patient_id,
      time_factor,
      time_order_numeric,
      Observed_ASV,
      Shannon,
      Simpson
    ) |>
    pivot_longer(
      cols = c(
        Observed_ASV,
        Shannon,
        Simpson
      ),
      names_to = "metric",
      values_to = "value"
    )

  g <- ggplot(
    dd,
    aes(
      x = time_factor,
      y = value
    )
  ) +
    geom_boxplot(
      outlier.shape = NA
    ) +
    geom_jitter(
      width = 0.08,
      alpha = 0.45,
      size = 1.2
    ) +
    facet_wrap(
      ~metric,
      scales = "free_y",
      ncol = 1
    ) +
    labs(
      title = paste0(
        project,
        " longitudinal alpha diversity"
      ),
      x = "Time",
      y = NULL
    ) +
    theme_bw(
      base_size = 11
    ) +
    theme(
      axis.text.x = element_text(
        angle = 35,
        hjust = 1
      )
    )

  ggsave(
    p,
    g,
    width = 7.5,
    height = 9
  )
}

make_displacement_plot <- function(
  d,
  project,
  p
) {

  g <- ggplot(
    d,
    aes(
      x = time_factor,
      y = bray_from_patient_baseline,
      group = patient_id
    )
  ) +
    geom_line(
      alpha = 0.22
    ) +
    geom_point(
      alpha = 0.55,
      size = 1.4
    ) +
    stat_summary(
      aes(group = 1),
      fun = median,
      geom = "line",
      linewidth = 1.1
    ) +
    stat_summary(
      aes(group = 1),
      fun = median,
      geom = "point",
      size = 2.4
    ) +
    labs(
      title = paste0(
        project,
        " within-patient Bray-Curtis displacement"
      ),
      x = "Time",
      y = "Bray-Curtis distance from earliest longitudinal sample"
    ) +
    theme_bw(
      base_size = 11
    ) +
    theme(
      axis.text.x = element_text(
        angle = 35,
        hjust = 1
      )
    )

  ggsave(
    p,
    g,
    width = 7.5,
    height = 5.5
  )
}

make_pcoa_plot <- function(
  d,
  project,
  p
) {

  g <- ggplot(
    d,
    aes(
      x = Axis1,
      y = Axis2,
      shape = time_factor
    )
  ) +
    geom_point(
      alpha = 0.72,
      size = 2
    ) +
    labs(
      title = paste0(
        project,
        " Bray-Curtis PCoA"
      ),
      x = "PCoA Axis 1",
      y = "PCoA Axis 2",
      shape = "Time"
    ) +
    theme_bw(
      base_size = 11
    )

  ggsave(
    p,
    g,
    width = 7.2,
    height = 5.5
  )
}

# ------------------------------------------------------------
# Main
# ------------------------------------------------------------

main <- function() {

  ck("STEP88A STARTED")

  complete87b <- file.path(
    STEP87B,
    "_STEP87B_COMPLETE.ok"
  )

  object_registry_path <- file.path(
    STEP87B,
    "V2_STEP87B_analysis_object_registry.csv"
  )

  primary_meta_path <- file.path(
    STEP87B,
    "V2_STEP87B_PRIMARY_analysis_metadata_785_patient_time_observations.csv"
  )

  required_input <- c(
    complete87b,
    object_registry_path,
    primary_meta_path
  )

  if (!all(file.exists(required_input))) {
    stop(
      paste0(
        "Required Step87B files are missing:\n",
        paste(
          required_input[
            !file.exists(required_input)
          ],
          collapse = "\n"
        )
      )
    )
  }

  registry <- safe_csv(
    object_registry_path
  )

  frozen_meta <- safe_csv(
    primary_meta_path
  )

  if (
    nrow(registry) != 8 ||
    sum(registry$primary_samples) != 785 ||
    nrow(frozen_meta) != 785
  ) {
    stop(
      "Step87B state guard failed. Expected 8 objects and 785 primary observations."
    )
  }

  if (
    sum(registry$primary_patients) != 390 ||
    sum(registry$primary_patients_GE2) != 258 ||
    sum(registry$primary_patients_GE3) != 132
  ) {
    stop(
      "Step87B patient eligibility guard failed. Expected 390 / 258 / 132."
    )
  }

  ck("STEP87B STATE GUARDED")

  qc_rows <- list()
  alpha_rows <- list()
  alpha_summary_rows <- list()
  displacement_rows <- list()
  displacement_summary_rows <- list()
  pcoa_rows <- list()
  model_global_rows <- list()
  model_coef_rows <- list()
  model_status_rows <- list()

  for (proj in project_order) {

    ck(
      paste0(
        proj,
        " START"
      )
    )

    rr <- registry |>
      filter(
        project == proj
      )

    if (nrow(rr) != 1) {
      stop(
        proj,
        ": object registry row missing or non-unique."
      )
    }

    object_path <- as.character(
      rr$analysis_object_path[1]
    )

    if (!file.exists(object_path)) {
      stop(
        paste0(
          proj,
          ": Step87B analysis object does not exist:\n",
          object_path,
          "\nThe Step87B result ZIP contains the registry, but the actual RDS objects remain under data/_V2_ANALYSIS_READY/01_PROJECTS/..."
        )
      )
    }

    obj <- readRDS(
      object_path
    )

    if (
      !all(
        c(
          "counts",
          "taxonomy",
          "metadata",
          "provenance"
        ) %in% names(obj)
      )
    ) {
      stop(
        proj,
        ": malformed Step87B analysis object."
      )
    }

    counts <- obj$counts
    md <- as_tibble(
      obj$metadata
    )

    if (
      !is.matrix(counts) ||
      nrow(counts) != nrow(md) ||
      !identical(
        rownames(counts),
        as.character(md$run_id)
      )
    ) {
      stop(
        proj,
        ": count/metadata alignment failed."
      )
    }

    lib <- rowSums(counts)

    if (any(lib < 2000)) {
      stop(
        proj,
        ": found library size below frozen 2,000-read threshold."
      )
    }

    # ========================================================
    # Alpha metrics for every primary sample
    # ========================================================

    alpha <- md |>
      mutate(
        library_size = as.numeric(lib),
        log10_library_size = log10(
          library_size
        ),
        Observed_ASV = as.numeric(
          rowSums(
            counts > 0
          )
        ),
        Shannon = as.numeric(
          vegan::diversity(
            counts,
            index = "shannon"
          )
        ),
        Simpson = as.numeric(
          vegan::diversity(
            counts,
            index = "simpson"
          )
        )
      )

    alpha_rows[[
      length(alpha_rows) + 1
    ]] <- alpha

    write_csv_safe(
      alpha,
      file.path(
        DIR_ALPHA,
        paste0(
          proj,
          "_sample_alpha_diversity.csv"
        )
      )
    )

    # ========================================================
    # Build true longitudinal GE2 subset
    # ========================================================

    long0 <- alpha |>
      filter(
        !is.na(
          analysis_time_order
        )
      )

    patient_n_time <- long0 |>
      distinct(
        patient_id,
        timepoint_key
      ) |>
      count(
        patient_id,
        name = "n_longitudinal_timepoints"
      )

    eligible_patients <- patient_n_time |>
      filter(
        n_longitudinal_timepoints >= 2
      ) |>
      pull(
        patient_id
      )

    long <- long0 |>
      filter(
        patient_id %in% eligible_patients
      )

    ee <- expected_longitudinal |>
      filter(
        project == proj
      )

    if (
      nrow(long) != ee$eligible_samples_GE2[1] ||
      n_distinct(long$patient_id) !=
        ee$eligible_patients_GE2[1]
    ) {
      stop(
        paste0(
          proj,
          ": longitudinal eligibility guard failed. ",
          "Expected ",
          ee$eligible_samples_GE2[1],
          " samples / ",
          ee$eligible_patients_GE2[1],
          " patients; observed ",
          nrow(long),
          " / ",
          n_distinct(long$patient_id),
          "."
        )
      )
    }

    if (nrow(long) > 0) {

      long <- long |>
        mutate(
          time_label = make_time_label(
            cur_data()
          ),
          time_order_numeric = as.numeric(
            analysis_time_order
          )
        )

      time_levels <- long |>
        distinct(
          time_order_numeric,
          time_label
        ) |>
        arrange(
          time_order_numeric
        ) |>
        pull(
          time_label
        )

      if (
        anyDuplicated(time_levels)
      ) {
        time_levels <- unique(
          time_levels
        )
      }

      long <- long |>
        mutate(
          time_factor = factor(
            time_label,
            levels = time_levels,
            ordered = TRUE
          )
        ) |>
        arrange(
          patient_id,
          time_order_numeric
        )

      # ======================================================
      # Alpha summaries
      # ======================================================

      alpha_summary <- long |>
        group_by(
          project,
          analysis_role,
          time_axis_type,
          time_order_numeric,
          time_factor
        ) |>
        summarise(
          n_samples = n(),
          n_patients = n_distinct(
            patient_id
          ),

          observed_median = median(
            Observed_ASV,
            na.rm = TRUE
          ),
          observed_q1 = quantile(
            Observed_ASV,
            0.25,
            na.rm = TRUE
          ),
          observed_q3 = quantile(
            Observed_ASV,
            0.75,
            na.rm = TRUE
          ),

          shannon_median = median(
            Shannon,
            na.rm = TRUE
          ),
          shannon_q1 = quantile(
            Shannon,
            0.25,
            na.rm = TRUE
          ),
          shannon_q3 = quantile(
            Shannon,
            0.75,
            na.rm = TRUE
          ),

          simpson_median = median(
            Simpson,
            na.rm = TRUE
          ),
          simpson_q1 = quantile(
            Simpson,
            0.25,
            na.rm = TRUE
          ),
          simpson_q3 = quantile(
            Simpson,
            0.75,
            na.rm = TRUE
          ),

          median_library_size = median(
            library_size,
            na.rm = TRUE
          ),

          .groups = "drop"
        )

      alpha_summary_rows[[
        length(alpha_summary_rows) + 1
      ]] <- alpha_summary

      write_csv_safe(
        alpha_summary,
        file.path(
          DIR_ALPHA,
          paste0(
            proj,
            "_longitudinal_alpha_summary.csv"
          )
        )
      )

      make_alpha_plot(
        long,
        proj,
        file.path(
          DIR_PLOT,
          paste0(
            proj,
            "_alpha_diversity_exploratory.pdf"
          )
        )
      )

      # ======================================================
      # Relative abundance + Bray-Curtis
      # ======================================================

      long_counts <- counts[
        long$run_id,
        ,
        drop = FALSE
      ]

      rel <- long_counts /
        rowSums(
          long_counts
        )

      bray <- vegan::vegdist(
        rel,
        method = "bray"
      )

      bray_mat <- as.matrix(
        bray
      )

      baseline_map <- long |>
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
          baseline_run_id = run_id,
          baseline_time_order =
            time_order_numeric,
          baseline_time_label =
            time_label
        )

      displacement <- long |>
        left_join(
          baseline_map,
          by = "patient_id"
        ) |>
        rowwise() |>
        mutate(
          bray_from_patient_baseline =
            as.numeric(
              bray_mat[
                run_id,
                baseline_run_id
              ]
            )
        ) |>
        ungroup()

      if (
        any(
          is.na(
            displacement$
              bray_from_patient_baseline
          )
        )
      ) {
        stop(
          proj,
          ": NA Bray-Curtis displacement generated."
        )
      }

      baseline_rows <- displacement |>
        filter(
          run_id == baseline_run_id
        )

      if (
        any(
          abs(
            baseline_rows$
              bray_from_patient_baseline
          ) > 1e-12
        )
      ) {
        stop(
          proj,
          ": baseline Bray-Curtis distance is not zero."
        )
      }

      displacement_rows[[
        length(displacement_rows) + 1
      ]] <- displacement

      write_csv_safe(
        displacement,
        file.path(
          DIR_BETA,
          paste0(
            proj,
            "_within_patient_bray_displacement.csv"
          )
        )
      )

      displacement_summary <- displacement |>
        group_by(
          project,
          analysis_role,
          time_axis_type,
          time_order_numeric,
          time_factor
        ) |>
        summarise(
          n_samples = n(),
          n_patients = n_distinct(
            patient_id
          ),
          bray_median = median(
            bray_from_patient_baseline,
            na.rm = TRUE
          ),
          bray_q1 = quantile(
            bray_from_patient_baseline,
            0.25,
            na.rm = TRUE
          ),
          bray_q3 = quantile(
            bray_from_patient_baseline,
            0.75,
            na.rm = TRUE
          ),
          bray_mean = mean(
            bray_from_patient_baseline,
            na.rm = TRUE
          ),
          bray_sd = sd(
            bray_from_patient_baseline,
            na.rm = TRUE
          ),
          .groups = "drop"
        )

      displacement_summary_rows[[
        length(displacement_summary_rows) + 1
      ]] <- displacement_summary

      write_csv_safe(
        displacement_summary,
        file.path(
          DIR_BETA,
          paste0(
            proj,
            "_within_patient_bray_displacement_summary.csv"
          )
        )
      )

      make_displacement_plot(
        displacement,
        proj,
        file.path(
          DIR_PLOT,
          paste0(
            proj,
            "_bray_displacement_exploratory.pdf"
          )
        )
      )

      # ======================================================
      # PCoA
      # ======================================================

      pcoa_fit <- stats::cmdscale(
        bray,
        k = 2,
        eig = TRUE,
        add = TRUE
      )

      pcoa_points <- as.data.frame(
        pcoa_fit$points
      )

      pcoa_points <- rownames_to_column(
        pcoa_points,
        "run_id"
      )

      names(pcoa_points)[
        names(pcoa_points) == "V1"
      ] <- "Axis1"

      names(pcoa_points)[
        names(pcoa_points) == "V2"
      ] <- "Axis2"

      pcoa <- long |>
        select(
          project,
          run_id,
          patient_id,
          analysis_role,
          time_axis_type,
          time_day,
          time_order_numeric,
          time_factor
        ) |>
        left_join(
          pcoa_points,
          by = "run_id"
        )

      eig <- pcoa_fit$eig

      positive_eig <- eig[
        eig > 0
      ]

      var1 <- if (
        length(positive_eig) >= 1
      ) {
        100 * positive_eig[1] /
          sum(positive_eig)
      } else {
        NA_real_
      }

      var2 <- if (
        length(positive_eig) >= 2
      ) {
        100 * positive_eig[2] /
          sum(positive_eig)
      } else {
        NA_real_
      }

      pcoa <- pcoa |>
        mutate(
          Axis1_variance_percent =
            var1,
          Axis2_variance_percent =
            var2
        )

      pcoa_rows[[
        length(pcoa_rows) + 1
      ]] <- pcoa

      write_csv_safe(
        pcoa,
        file.path(
          DIR_PCOA,
          paste0(
            proj,
            "_Bray_PCoA_coordinates.csv"
          )
        )
      )

      make_pcoa_plot(
        pcoa,
        proj,
        file.path(
          DIR_PLOT,
          paste0(
            proj,
            "_Bray_PCoA_exploratory.pdf"
          )
        )
      )

      # ======================================================
      # Mixed models
      # ======================================================

      adj <- resolve_adjustment(
        proj,
        long
      )

      for (
        metric in c(
          "Observed_ASV",
          "Shannon",
          "Simpson"
        )
      ) {

        m_factor <- safe_lmer(
          d = long,
          outcome = metric,
          time_mode = "factor",
          include_depth = TRUE,
          adjustment = adj
        )

        status_row <- tibble(
          project = proj,
          analysis_family =
            "ALPHA_DIVERSITY",
          outcome = metric,
          model_type =
            "CATEGORICAL_TIME_MIXED_MODEL",
          status = m_factor$status
        )

        model_status_rows[[
          length(model_status_rows) + 1
        ]] <- status_row

        if (
          m_factor$status == "OK"
        ) {

          model_global_rows[[
            length(model_global_rows) + 1
          ]] <- m_factor$global |>
            mutate(
              project = proj,
              analysis_family =
                "ALPHA_DIVERSITY",
              model_type =
                "CATEGORICAL_TIME_MIXED_MODEL",
              .before = 1
            )

          model_coef_rows[[
            length(model_coef_rows) + 1
          ]] <- m_factor$coefficients |>
            mutate(
              project = proj,
              analysis_family =
                "ALPHA_DIVERSITY",
              outcome = metric,
              model_type =
                "CATEGORICAL_TIME_MIXED_MODEL",
              .before = 1
            )
        }

        if (
          all(
            !is.na(
              long$time_day
            )
          ) &&
          n_distinct(
            long$time_day
          ) >= 2
        ) {

          m_numeric <- safe_lmer(
            d = long,
            outcome = metric,
            time_mode = "numeric",
            include_depth = TRUE,
            adjustment = adj
          )

          model_status_rows[[
            length(model_status_rows) + 1
          ]] <- tibble(
            project = proj,
            analysis_family =
              "ALPHA_DIVERSITY",
            outcome = metric,
            model_type =
              "CONTINUOUS_DAY_MIXED_MODEL",
            status = m_numeric$status
          )

          if (
            m_numeric$status == "OK"
          ) {

            model_global_rows[[
              length(model_global_rows) + 1
            ]] <- m_numeric$global |>
              mutate(
                project = proj,
                analysis_family =
                  "ALPHA_DIVERSITY",
                model_type =
                  "CONTINUOUS_DAY_MIXED_MODEL",
                .before = 1
              )

            model_coef_rows[[
              length(model_coef_rows) + 1
            ]] <- m_numeric$coefficients |>
              mutate(
                project = proj,
                analysis_family =
                  "ALPHA_DIVERSITY",
                outcome = metric,
                model_type =
                  "CONTINUOUS_DAY_MIXED_MODEL",
                .before = 1
              )
          }
        }
      }

      # Bray displacement models
      b_factor <- safe_lmer(
        d = displacement,
        outcome =
          "bray_from_patient_baseline",
        time_mode = "factor",
        include_depth = FALSE,
        adjustment = adj
      )

      model_status_rows[[
        length(model_status_rows) + 1
      ]] <- tibble(
        project = proj,
        analysis_family =
          "WITHIN_PATIENT_BRAY_DISPLACEMENT",
        outcome =
          "bray_from_patient_baseline",
        model_type =
          "CATEGORICAL_TIME_MIXED_MODEL",
        status = b_factor$status
      )

      if (
        b_factor$status == "OK"
      ) {

        model_global_rows[[
          length(model_global_rows) + 1
        ]] <- b_factor$global |>
          mutate(
            project = proj,
            analysis_family =
              "WITHIN_PATIENT_BRAY_DISPLACEMENT",
            model_type =
              "CATEGORICAL_TIME_MIXED_MODEL",
            .before = 1
          )

        model_coef_rows[[
          length(model_coef_rows) + 1
        ]] <- b_factor$coefficients |>
          mutate(
            project = proj,
            analysis_family =
              "WITHIN_PATIENT_BRAY_DISPLACEMENT",
            outcome =
              "bray_from_patient_baseline",
            model_type =
              "CATEGORICAL_TIME_MIXED_MODEL",
            .before = 1
          )
      }

      if (
        all(
          !is.na(
            displacement$time_day
          )
        ) &&
        n_distinct(
          displacement$time_day
        ) >= 2
      ) {

        b_numeric <- safe_lmer(
          d = displacement,
          outcome =
            "bray_from_patient_baseline",
          time_mode = "numeric",
          include_depth = FALSE,
          adjustment = adj
        )

        model_status_rows[[
          length(model_status_rows) + 1
        ]] <- tibble(
          project = proj,
          analysis_family =
            "WITHIN_PATIENT_BRAY_DISPLACEMENT",
          outcome =
            "bray_from_patient_baseline",
          model_type =
            "CONTINUOUS_DAY_MIXED_MODEL",
          status = b_numeric$status
        )

        if (
          b_numeric$status == "OK"
        ) {

          model_global_rows[[
            length(model_global_rows) + 1
          ]] <- b_numeric$global |>
            mutate(
              project = proj,
              analysis_family =
                "WITHIN_PATIENT_BRAY_DISPLACEMENT",
              model_type =
                "CONTINUOUS_DAY_MIXED_MODEL",
              .before = 1
            )

          model_coef_rows[[
            length(model_coef_rows) + 1
          ]] <- b_numeric$coefficients |>
            mutate(
              project = proj,
              analysis_family =
                "WITHIN_PATIENT_BRAY_DISPLACEMENT",
              outcome =
                "bray_from_patient_baseline",
              model_type =
                "CONTINUOUS_DAY_MIXED_MODEL",
              .before = 1
            )
        }
      }
    }

    qc_rows[[
      length(qc_rows) + 1
    ]] <- tibble(
      project = proj,
      primary_samples =
        nrow(md),
      primary_patients =
        n_distinct(
          md$patient_id
        ),
      final_ASVs =
        ncol(counts),
      min_library_size =
        min(lib),
      median_library_size =
        median(lib),
      max_library_size =
        max(lib),
      longitudinal_axis_samples =
        nrow(long0),
      longitudinal_axis_patients =
        n_distinct(
          long0$patient_id
        ),
      eligible_GE2_samples =
        nrow(long),
      eligible_GE2_patients =
        n_distinct(
          long$patient_id
        ),
      time_axis_types = paste(
        sort(
          unique(
            na.omit(
              md$time_axis_type
            )
          )
        ),
        collapse = ";"
      )
    )

    ck(
      paste0(
        proj,
        " COMPLETE"
      )
    )

    rm(
      obj,
      counts,
      md,
      alpha,
      long0,
      long
    )

    gc(
      verbose = FALSE
    )
  }

  qc <- bind_rows(
    qc_rows
  )

  alpha_all <- bind_rows(
    alpha_rows
  )

  alpha_summary_all <- bind_rows(
    alpha_summary_rows
  )

  displacement_all <- bind_rows(
    displacement_rows
  )

  displacement_summary_all <- bind_rows(
    displacement_summary_rows
  )

  pcoa_all <- bind_rows(
    pcoa_rows
  )

  model_global <- bind_rows(
    model_global_rows
  )

  model_coef <- bind_rows(
    model_coef_rows
  )

  model_status <- bind_rows(
    model_status_rows
  )

  # ----------------------------------------------------------
  # Final guards
  # ----------------------------------------------------------

  if (
    nrow(qc) != 8 ||
    sum(qc$primary_samples) != 785 ||
    sum(qc$primary_patients) != 390
  ) {
    stop(
      "Final Step88A primary-state guard failed."
    )
  }

  if (
    sum(
      qc$eligible_GE2_samples
    ) != 653 ||
    sum(
      qc$eligible_GE2_patients
    ) != 258
  ) {
    stop(
      paste0(
        "Final Step88A longitudinal-state guard failed. ",
        "Expected 653 samples / 258 patients; observed ",
        sum(
          qc$eligible_GE2_samples
        ),
        " / ",
        sum(
          qc$eligible_GE2_patients
        ),
        "."
      )
    )
  }

  core_qc <- qc |>
    filter(
      project ==
        "PRJNA691455"
    )

  if (
    core_qc$eligible_GE2_samples != 29 ||
    core_qc$eligible_GE2_patients != 10
  ) {
    stop(
      "Core PRJNA691455 longitudinal guard failed. Expected 29 samples / 10 patients."
    )
  }

  # ----------------------------------------------------------
  # Write global outputs
  # ----------------------------------------------------------

  write_csv_safe(
    qc,
    file.path(
      DIR_QC,
      "V2_STEP88A_cohort_QC_summary.csv"
    )
  )

  write_csv_safe(
    alpha_all,
    file.path(
      DIR_ALPHA,
      "V2_STEP88A_ALL_sample_alpha_diversity.csv"
    )
  )

  write_csv_safe(
    alpha_summary_all,
    file.path(
      DIR_ALPHA,
      "V2_STEP88A_ALL_longitudinal_alpha_summary.csv"
    )
  )

  write_csv_safe(
    displacement_all,
    file.path(
      DIR_BETA,
      "V2_STEP88A_ALL_within_patient_bray_displacement.csv"
    )
  )

  write_csv_safe(
    displacement_summary_all,
    file.path(
      DIR_BETA,
      "V2_STEP88A_ALL_within_patient_bray_displacement_summary.csv"
    )
  )

  write_csv_safe(
    pcoa_all,
    file.path(
      DIR_PCOA,
      "V2_STEP88A_ALL_Bray_PCoA_coordinates.csv"
    )
  )

  write_csv_safe(
    model_status,
    file.path(
      DIR_MODEL,
      "V2_STEP88A_model_status.csv"
    )
  )

  write_csv_safe(
    model_global,
    file.path(
      DIR_MODEL,
      "V2_STEP88A_model_global_time_tests.csv"
    )
  )

  write_csv_safe(
    model_coef,
    file.path(
      DIR_MODEL,
      "V2_STEP88A_model_fixed_effect_coefficients.csv"
    )
  )

  # ----------------------------------------------------------
  # README / completion
  # ----------------------------------------------------------

  n_model_fail <- sum(
    model_status$status != "OK"
  )

  readme <- c(
    "SEPSIS V2 - STEP88A LONGITUDINAL DIVERSITY AND DISPLACEMENT",
    paste0(
      "Created: ",
      Sys.time()
    ),
    "",
    "FROZEN INPUT",
    "785 primary patient-time observations across 8 cohorts.",
    "390 project-specific patients.",
    "",
    "TRUE REPEATED-MEASURES SUBSET",
    "653 samples from 258 patients with >=2 longitudinal timepoints.",
    "Static controls and one-time samples are not forced into trajectory models.",
    "PRJNA691455 core sepsis trajectory = 29 samples from 10 patients.",
    "",
    "ALPHA DIVERSITY",
    "Observed ASVs, Shannon, and Simpson are calculated per sample.",
    "Categorical-time mixed models are fit for all longitudinal cohorts.",
    "Continuous-day mixed models are additionally fit where numeric day is valid.",
    "Alpha mixed models include log10 library size; PRJEB82425 is additionally adjusted for phenotype; intervention cohorts are adjusted for intervention arm.",
    "",
    "BETA DIVERSITY",
    "Bray-Curtis is calculated on within-cohort relative abundance tables.",
    "Each patient's earliest longitudinal observation is used as that patient's displacement reference.",
    "No raw ASV table is merged across cohorts.",
    "",
    "PCoA",
    "Cohort-specific Bray-Curtis PCoA coordinates are exported for exploratory visualization.",
    "",
    "MODEL INTERPRETATION",
    "Step88A models are cohort-specific trajectory tests. They are not yet cross-cohort meta-analysis.",
    "Do not interpret intervention-support cohorts as natural-history sepsis replication cohorts.",
    "",
    paste0(
      "Model fits with non-OK status: ",
      n_model_fail,
      ". Check 05_MODELS/V2_STEP88A_model_status.csv."
    ),
    "",
    "NEXT",
    "After auditing Step88A results, proceed to taxonomic longitudinal trajectories and then cross-cohort direction/effect synthesis."
  )

  writeLines(
    readme,
    file.path(
      OUT,
      "README_STEP88A.txt"
    ),
    useBytes = TRUE
  )

  completion_status <- if (
    n_model_fail == 0
  ) {
    "STEP88A COMPLETE"
  } else {
    "STEP88A COMPLETE WITH MODEL WARNINGS"
  }

  writeLines(
    c(
      paste0(
        "Completed: ",
        Sys.time()
      ),
      paste0(
        "Status: ",
        completion_status
      ),
      "Primary observations audited: 785",
      "Project-specific patients audited: 390",
      "Longitudinal GE2 samples analyzed: 653",
      "Longitudinal GE2 patients analyzed: 258",
      "Core PRJNA691455 longitudinal samples: 29",
      "Core PRJNA691455 longitudinal patients: 10",
      paste0(
        "Non-OK mixed-model fits: ",
        n_model_fail
      )
    ),
    file.path(
      OUT,
      "_STEP88A_COMPLETE.ok"
    ),
    useBytes = TRUE
  )

  ck(
    "STEP88A COMPLETE"
  )

  cat(
    "\n============================================================\n"
  )
  cat(
    "SEPSIS V2 - STEP88A COMPLETE\n"
  )
  cat(
    "============================================================\n\n"
  )

  cat(
    "Cohort QC summary:\n"
  )

  print(
    qc,
    n = Inf,
    width = Inf
  )

  cat(
    "\nModel status:\n"
  )

  print(
    model_status,
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
        "STEP88A FATAL ERROR: ",
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
      "STEP88A FAILED"
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
