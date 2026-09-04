# ============================================================
# Sepsis V2 - Step 91B4
# Literature-validated infection-source freeze and source-stratified
# ecological trajectory analysis
#
# CORE COHORT:
# PRJNA691455
# Source: Yang XJ et al. World J Gastroenterol. 2021;27(19):2376-2393.
# DOI: 10.3748/wjg.v27.i19.2376
# Published Table 2 gives patient-level PRIMARY INFECTION SITE.
#
# S1  Lungs
# S2  Abdominal cavity
# S3  Lungs
# S4  Lungs
# S5  Abdomen
# S6  Lungs
# S7  Lungs
# S8  Abdominal cavity
# S9  Abdominal cavity
# S10 Abdominal cavity
#
# Standardized:
# PULMONARY = S1,S3,S4,S6,S7
# ABDOMINAL = S2,S5,S8,S9,S10
#
# EXTERNAL SUPPORT:
# PRJEB82425
# Published study groups:
# Pneumonia/VAP (n=17), Other infection (n=12), Control (n=9)
# Supplement Table S1 reports aggregate "Other infection" sites:
#   Fever of unknown origin = 9
#   UTI = 1
#   Primary bacteremia/catheter infection = 2
# Individual subtype mapping is NOT publicly recoverable from the
# available frozen metadata / supplement, so those 12 are NOT split.
#
# ANALYSIS PRINCIPLE:
# - Core infection-source analysis is exploratory due n=10.
# - Primary source-modification endpoint:
#     Day7 Bray distance from patient-specific Day1 baseline
#     PULMONARY vs ABDOMINAL, n=5 vs 5.
# - Secondary:
#     Day7-Day3 incremental Bray displacement, complete n=4 vs 5.
# - Exact label-permutation tests are used because sample size is tiny.
# - Alpha-diversity deltas are exploratory only.
# - PRJEB82425 is descriptive external support only.
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

set.seed(20260821)

ROOT <- "E:/sepsis_project"

STEP87B <- file.path(
  ROOT,
  "results",
  "V2_27B_ANALYSIS_OBJECTS_AND_REPLICATE_FREEZE"
)

STEP88A2 <- file.path(
  ROOT,
  "results",
  "V2_28A2_STEP88A_LONGITUDINAL_DIVERSITY_AND_DISPLACEMENT_FIXED"
)

STEP91A <- file.path(
  ROOT,
  "results",
  "V2_31A_STEP91A_INFECTION_SOURCE_FEASIBILITY_AUDIT"
)

OUT <- file.path(
  ROOT,
  "results",
  "V2_31B4_STEP91B_INFECTION_SOURCE_FREEZE_AND_ANALYSIS_ATTRITION_FIXED"
)

DIR_FREEZE <- file.path(
  OUT,
  "01_SOURCE_FREEZE"
)

DIR_CORE <- file.path(
  OUT,
  "02_CORE_SEPSIS_SOURCE_ANALYSIS"
)

DIR_EXT <- file.path(
  OUT,
  "03_EXTERNAL_SUPPORT_PRJEB82425"
)

DIR_SYN <- file.path(
  OUT,
  "04_SYNTHESIS"
)

for (d in c(
  OUT,
  DIR_FREEZE,
  DIR_CORE,
  DIR_EXT,
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
  "_STEP91B4_runtime_checkpoints.txt"
)

ERR <- file.path(
  OUT,
  "_STEP91B4_FATAL_ERROR.txt"
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

exact_two_group_permutation <- function(
  y,
  group
) {

  y <- as.numeric(y)
  group <- as.character(group)

  ok <- complete.cases(
    y,
    group
  )

  y <- y[ok]
  group <- group[ok]

  lev <- sort(
    unique(group)
  )

  if (
    length(lev) != 2
  ) {
    stop(
      "exact_two_group_permutation requires exactly two groups."
    )
  }

  n1 <- sum(
    group ==
      lev[1]
  )

  n <- length(y)

  obs <- mean(
    y[
      group ==
        lev[1]
    ]
  ) -
    mean(
      y[
        group ==
          lev[2]
    ]
  )

  combos <- combn(
    n,
    n1
  )

  perm <- apply(
    combos,
    2,
    function(idx) {
      mean(
        y[idx]
      ) -
        mean(
          y[
            -idx
          ]
        )
    }
  )

  p <- mean(
    abs(perm) >=
      abs(obs) -
      1e-12
  )

  tibble(
    group_1 =
      lev[1],
    group_2 =
      lev[2],
    n_group_1 =
      n1,
    n_group_2 =
      n - n1,
    mean_group_1 =
      mean(
        y[
          group ==
            lev[1]
        ]
      ),
    mean_group_2 =
      mean(
        y[
          group ==
            lev[2]
        ]
      ),
    median_group_1 =
      median(
        y[
          group ==
            lev[1]
        ]
      ),
    median_group_2 =
      median(
        y[
          group ==
            lev[2]
        ]
      ),
    mean_difference_group1_minus_group2 =
      obs,
    exact_permutation_p =
      p,
    n_exact_label_allocations =
      ncol(combos)
  )
}

bootstrap_mean_difference <- function(
  y,
  group,
  B = 5000,
  seed = 20260821
) {

  y <- as.numeric(y)
  group <- as.character(group)

  ok <- complete.cases(
    y,
    group
  )

  y <- y[ok]
  group <- group[ok]

  lev <- sort(
    unique(group)
  )

  if (
    length(lev) != 2
  ) {
    return(
      c(
        low = NA_real_,
        high = NA_real_
      )
    )
  }

  y1 <- y[
    group ==
      lev[1]
  ]

  y2 <- y[
    group ==
      lev[2]
  ]

  if (
    length(y1) < 2 ||
    length(y2) < 2
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
      mean(
        sample(
          y1,
          length(y1),
          replace = TRUE
        )
      ) -
        mean(
          sample(
            y2,
            length(y2),
            replace = TRUE
          )
        )
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

hedges_g <- function(
  y,
  group
) {

  y <- as.numeric(y)
  group <- as.character(group)

  ok <- complete.cases(
    y,
    group
  )

  y <- y[ok]
  group <- group[ok]

  lev <- sort(
    unique(group)
  )

  if (
    length(lev) != 2
  ) {
    return(NA_real_)
  }

  y1 <- y[
    group ==
      lev[1]
  ]

  y2 <- y[
    group ==
      lev[2]
  ]

  n1 <- length(y1)
  n2 <- length(y2)

  if (
    n1 < 2 ||
    n2 < 2
  ) {
    return(NA_real_)
  }

  sp <- sqrt(
    (
      (n1 - 1) *
        var(y1) +
      (n2 - 1) *
        var(y2)
    ) /
      (
        n1 +
          n2 -
          2
      )
  )

  if (
    is.na(sp) ||
    sp == 0
  ) {
    return(NA_real_)
  }

  d <- (
    mean(y1) -
      mean(y2)
  ) /
    sp

  J <- 1 -
    3 /
      (
        4 *
          (
            n1 +
              n2
          ) -
          9
      )

  J * d
}

cliffs_delta <- function(
  y,
  group
) {

  y <- as.numeric(y)
  group <- as.character(group)

  ok <- complete.cases(
    y,
    group
  )

  y <- y[ok]
  group <- group[ok]

  lev <- sort(
    unique(group)
  )

  if (
    length(lev) != 2
  ) {
    return(NA_real_)
  }

  y1 <- y[
    group ==
      lev[1]
  ]

  y2 <- y[
    group ==
      lev[2]
  ]

  comp <- outer(
    y1,
    y2,
    "-"
  )

  (
    sum(
      comp > 0
    ) -
      sum(
        comp < 0
      )
  ) /
    length(comp)
}

source_compare <- function(
  data,
  outcome,
  analysis_label,
  seed_offset = 0L
) {

  d <- data |>
    filter(
      !is.na(
        .data[[outcome]]
      ),
      !is.na(
        infection_source_standard
      )
    )

  if (
    n_distinct(
      d$infection_source_standard
    ) != 2
  ) {
    stop(
      analysis_label,
      ": expected two source groups."
    )
  }

  ex <- exact_two_group_permutation(
    d[[outcome]],
    d$infection_source_standard
  )

  ci <- bootstrap_mean_difference(
    d[[outcome]],
    d$infection_source_standard,
    B = 5000,
    seed =
      20260821 +
      seed_offset
  )

  ex |>
    mutate(
      analysis =
        analysis_label,
      outcome =
        outcome,
      mean_difference_ci_low =
        ci["low"],
      mean_difference_ci_high =
        ci["high"],
      hedges_g =
        hedges_g(
          d[[outcome]],
          d$infection_source_standard
        ),
      cliffs_delta =
        cliffs_delta(
          d[[outcome]],
          d$infection_source_standard
        ),
      .before = 1
    )
}

main <- function() {

  ck(
    "STEP91B4 STARTED"
  )

  registry_path <- file.path(
    STEP87B,
    "V2_STEP87B_analysis_object_registry.csv"
  )

  step87_complete <- file.path(
    STEP87B,
    "_STEP87B_COMPLETE.ok"
  )

  disp_path <- file.path(
    STEP88A2,
    "03_BETA_DISPLACEMENT",
    "V2_STEP88A2_ALL_within_patient_bray_displacement.csv"
  )

  step88_complete <- file.path(
    STEP88A2,
    "_STEP88A2_COMPLETE.ok"
  )

  step91a_complete <- file.path(
    STEP91A,
    "_STEP91A_COMPLETE.ok"
  )

  required <- c(
    registry_path,
    step87_complete,
    disp_path,
    step88_complete,
    step91a_complete
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

  registry <- safe_csv(
    registry_path
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

  disp <- safe_csv(
    disp_path
  )

  if (
    nrow(disp) !=
      653
  ) {
    stop(
      "Step88A2 longitudinal displacement table expected 653 rows; observed ",
      nrow(disp),
      "."
    )
  }

  ck(
    "UPSTREAM STATE GUARDED"
  )

  # ----------------------------------------------------------
  # 1. Freeze literature-validated PRJNA691455 source mapping
  # ----------------------------------------------------------

  source_map <- tribble(
    ~patient_id,
    ~published_patient_no,
    ~primary_infection_site_raw,
    ~infection_source_standard,

    "PRJNA691455_S1",
    "S1",
    "Lungs",
    "PULMONARY",

    "PRJNA691455_S2",
    "S2",
    "Abdominal cavity",
    "ABDOMINAL",

    "PRJNA691455_S3",
    "S3",
    "Lungs",
    "PULMONARY",

    "PRJNA691455_S4",
    "S4",
    "Lungs",
    "PULMONARY",

    "PRJNA691455_S5",
    "S5",
    "Abdomen",
    "ABDOMINAL",

    "PRJNA691455_S6",
    "S6",
    "Lungs",
    "PULMONARY",

    "PRJNA691455_S7",
    "S7",
    "Lungs",
    "PULMONARY",

    "PRJNA691455_S8",
    "S8",
    "Abdominal cavity",
    "ABDOMINAL",

    "PRJNA691455_S9",
    "S9",
    "Abdominal cavity",
    "ABDOMINAL",

    "PRJNA691455_S10",
    "S10",
    "Abdominal cavity",
    "ABDOMINAL"
  ) |>
    mutate(
      project =
        "PRJNA691455",
      mapping_level =
        "PATIENT_LEVEL_PRIMARY_INFECTION_SITE",
      mapping_confidence =
        "HIGH",
      provenance_type =
        "PUBLISHED_TABLE_2",
      provenance_citation =
        "Yang XJ et al. World J Gastroenterol. 2021;27(19):2376-2393.",
      provenance_doi =
        "10.3748/wjg.v27.i19.2376"
    ) |>
    select(
      project,
      patient_id,
      published_patient_no,
      primary_infection_site_raw,
      infection_source_standard,
      mapping_level,
      mapping_confidence,
      provenance_type,
      provenance_citation,
      provenance_doi
    )

  core <- disp |>
    filter(
      project ==
        "PRJNA691455"
    )

  core_ids <- sort(
    unique(
      core$patient_id
    )
  )

  if (
    !setequal(
      core_ids,
      source_map$patient_id
    )
  ) {
    stop(
      paste0(
        "PRJNA691455 patient ID mapping guard failed.\n",
        "Observed: ",
        paste(
          core_ids,
          collapse = ", "
        ),
        "\nExpected: ",
        paste(
          sort(
            source_map$patient_id
          ),
          collapse = ", "
        )
      )
    )
  }

  if (
    sum(
      source_map$infection_source_standard ==
        "PULMONARY"
    ) != 5 ||
    sum(
      source_map$infection_source_standard ==
        "ABDOMINAL"
    ) != 5
  ) {
    stop(
      "Core source mapping is not 5 pulmonary / 5 abdominal."
    )
  }

  write_csv_safe(
    source_map,
    file.path(
      DIR_FREEZE,
      "V2_STEP91B4_PRJNA691455_literature_validated_source_mapping.csv"
    )
  )

  # Overall source feasibility freeze registry
  source_registry <- tribble(
    ~project,
    ~source_resolution,
    ~eligible_role,
    ~patient_level_source_available,
    ~notes,

    "PRJNA691455",
    "PULMONARY_VS_ABDOMINAL",
    "CORE_SOURCE_ANALYSIS",
    TRUE,
    "Published Table 2 maps all 10 sepsis patients: 5 pulmonary, 5 abdominal.",

    "PRJEB82425",
    "PNEUMONIA_VAP_VS_OTHER_INFECTION_VS_CONTROL",
    "COARSE_EXTERNAL_SUPPORT",
    TRUE,
    "Patient-level study group available, but the 12 Other_infection patients cannot be individually split into FUO/UTI/bacteremia from public metadata.",

    "PRJNA516701",
    "NONE",
    "NOT_ELIGIBLE_FOR_SOURCE_ANALYSIS",
    FALSE,
    "No defensible patient-level infection-source variable identified.",

    "PRJNA851469",
    "NONE",
    "NOT_ELIGIBLE_FOR_SOURCE_ANALYSIS",
    FALSE,
    "No defensible patient-level infection-source variable identified.",

    "PRJNA578267",
    "NONE",
    "NONSEPSIS_CONTROL_NOT_SOURCE_ANALYSIS",
    FALSE,
    "Non-sepsis longitudinal control.",

    "PRJNA430161",
    "NONE",
    "INTERVENTION_SUPPORT_NOT_SOURCE_ANALYSIS",
    FALSE,
    "Intervention support cohort.",

    "PRJNA1166732",
    "NONE",
    "INTERVENTION_SUPPORT_NOT_SOURCE_ANALYSIS",
    FALSE,
    "Intervention support cohort.",

    "PRJNA978257",
    "NONE",
    "STATIC_SUPPORT_NOT_SOURCE_ANALYSIS",
    FALSE,
    "Static support cohort."
  )

  write_csv_safe(
    source_registry,
    file.path(
      DIR_FREEZE,
      "V2_STEP91B4_source_analysis_eligibility_registry.csv"
    )
  )

  ck(
    "SOURCE MAPPING FROZEN"
  )

  # ----------------------------------------------------------
  # 2. Core sepsis source-stratified trajectory analysis
  # ----------------------------------------------------------

  # Step88A2 already carries a frozen metadata field named
  # infection_source_standard. Keep it for audit rather than
  # allowing left_join() to create .x/.y suffixes.
  if (
    "infection_source_standard" %in%
      names(core)
  ) {

    core <- core |>
      rename(
        infection_source_standard_frozen_pre91B =
          infection_source_standard
      )

  } else {

    core$infection_source_standard_frozen_pre91B <-
      NA_character_
  }

  core2 <- core |>
    left_join(
      source_map |>
        select(
          patient_id,
          primary_infection_site_raw,
          infection_source_standard
        ),
      by =
        "patient_id"
    )

  if (
    !(
      "infection_source_standard" %in%
        names(core2)
    )
  ) {
    stop(
      "Validated infection_source_standard column missing after source join."
    )
  }

  if (
    any(
      is.na(
        core2$infection_source_standard
      )
    )
  ) {
    stop(
      "Unmapped core sepsis samples remain after source join."
    )
  }

  # Audit any pre-existing nonmissing frozen values instead of
  # silently overwriting them. If present and conflicting, stop.
  source_conflict_audit <- core2 |>
    distinct(
      patient_id,
      infection_source_standard_frozen_pre91B,
      infection_source_standard
    ) |>
    mutate(
      frozen_pre91B_nonmissing =
        !is.na(
          infection_source_standard_frozen_pre91B
        ) &
        trimws(
          as.character(
            infection_source_standard_frozen_pre91B
          )
        ) != "",
      conflict =
        frozen_pre91B_nonmissing &
        toupper(
          trimws(
            as.character(
              infection_source_standard_frozen_pre91B
            )
          )
        ) !=
          toupper(
            trimws(
              as.character(
                infection_source_standard
              )
            )
          )
    )

  write_csv_safe(
    source_conflict_audit,
    file.path(
      DIR_FREEZE,
      "V2_STEP91B4_PRJNA691455_preexisting_vs_literature_source_audit.csv"
    )
  )

  if (
    any(
      source_conflict_audit$conflict,
      na.rm = TRUE
    )
  ) {
    stop(
      paste0(
        "Conflict detected between frozen pre-Step91B infection_source_standard ",
        "and literature-validated source mapping. Review audit CSV."
      )
    )
  }

  # One row per patient with source and timepoint-specific metrics.
  core_wide <- core2 |>
    select(
      patient_id,
      infection_source_standard,
      time_factor,
      bray_from_patient_baseline,
      Shannon,
      Simpson,
      Observed_ASV
    ) |>
    pivot_wider(
      names_from =
        time_factor,
      values_from =
        c(
          bray_from_patient_baseline,
          Shannon,
          Simpson,
          Observed_ASV
        ),
      names_sep =
        "__"
    ) |>
    mutate(
      bray_day7_from_day1 =
        .data[["bray_from_patient_baseline__Day 7"]],

      bray_increment_day7_minus_day3 =
        .data[["bray_from_patient_baseline__Day 7"]] -
          .data[["bray_from_patient_baseline__Day 3"]],

      shannon_delta_day7_minus_day1 =
        .data[["Shannon__Day 7"]] -
          .data[["Shannon__Day 1"]],

      simpson_delta_day7_minus_day1 =
        .data[["Simpson__Day 7"]] -
          .data[["Simpson__Day 1"]],

      observed_asv_delta_day7_minus_day1 =
        .data[["Observed_ASV__Day 7"]] -
          .data[["Observed_ASV__Day 1"]]
    )

  if (
    nrow(
      core_wide |>
        filter(
          !is.na(
            bray_day7_from_day1
          )
        )
    ) != 10
  ) {
    stop(
      "Expected all 10 core sepsis patients to have Day7 Bray displacement."
    )
  }

  if (
    nrow(
      core_wide |>
        filter(
          !is.na(
            bray_increment_day7_minus_day3
          )
        )
    ) != 9
  ) {
    stop(
      "Expected 9 core sepsis patients with complete Day3-Day7 Bray increment."
    )
  }

  complete_increment_counts <- core_wide |>
    filter(
      !is.na(
        bray_increment_day7_minus_day3
      )
    ) |>
    count(
      infection_source_standard
    )

  if (
    !(
      complete_increment_counts$n[
        complete_increment_counts$infection_source_standard ==
          "PULMONARY"
      ] == 4 &&
      complete_increment_counts$n[
        complete_increment_counts$infection_source_standard ==
          "ABDOMINAL"
      ] == 5
    )
  ) {
    stop(
      "Expected Day3-Day7 complete patients: Pulmonary=4, Abdominal=5."
    )
  }

  write_csv_safe(
    core_wide,
    file.path(
      DIR_CORE,
      "V2_STEP91B4_PRJNA691455_patient_level_source_trajectory_metrics.csv"
    )
  )

  analyses <- list(
    source_compare(
      core_wide,
      "bray_day7_from_day1",
      "PRIMARY_Day7_Bray_from_Day1_baseline",
      1
    ),

    source_compare(
      core_wide,
      "bray_increment_day7_minus_day3",
      "SECONDARY_Day7_minus_Day3_Bray_increment",
      2
    ),

    source_compare(
      core_wide,
      "shannon_delta_day7_minus_day1",
      "EXPLORATORY_Shannon_Day7_minus_Day1",
      3
    ),

    source_compare(
      core_wide,
      "simpson_delta_day7_minus_day1",
      "EXPLORATORY_Simpson_Day7_minus_Day1",
      4
    ),

    source_compare(
      core_wide,
      "observed_asv_delta_day7_minus_day1",
      "EXPLORATORY_ObservedASV_Day7_minus_Day1",
      5
    )
  )

  core_tests <- bind_rows(
    analyses
  ) |>
    mutate(
      interpretation_tier =
        case_when(
          str_starts(
            analysis,
            "PRIMARY_"
          ) ~
            "PRIMARY_SOURCE_MODIFICATION_ENDPOINT",

          str_starts(
            analysis,
            "SECONDARY_"
          ) ~
            "SECONDARY_SOURCE_MODIFICATION_ENDPOINT",

          TRUE ~
            "EXPLORATORY_ALPHA_ENDPOINT"
        ),
      alpha_exploratory_fdr =
        ifelse(
          interpretation_tier ==
            "EXPLORATORY_ALPHA_ENDPOINT",
          bh_adjust(
            ifelse(
              interpretation_tier ==
                "EXPLORATORY_ALPHA_ENDPOINT",
              exact_permutation_p,
              NA_real_
            )
          ),
          NA_real_
        )
    )

  # Recalculate exploratory FDR only across the 3 alpha outcomes.
  alpha_idx <- which(
    core_tests$interpretation_tier ==
      "EXPLORATORY_ALPHA_ENDPOINT"
  )

  if (
    length(
      alpha_idx
    )
  ) {
    core_tests$alpha_exploratory_fdr[
      alpha_idx
    ] <- p.adjust(
      core_tests$exact_permutation_p[
        alpha_idx
      ],
      method = "BH"
    )
  }

  write_csv_safe(
    core_tests,
    file.path(
      DIR_CORE,
      "V2_STEP91B4_PRJNA691455_source_stratified_exact_tests.csv"
    )
  )

  core_group_summary <- core_wide |>
    group_by(
      infection_source_standard
    ) |>
    summarise(
      n_patients =
        n(),
      median_Day7_Bray =
        median(
          bray_day7_from_day1,
          na.rm = TRUE
        ),
      mean_Day7_Bray =
        mean(
          bray_day7_from_day1,
          na.rm = TRUE
        ),
      n_complete_Day3_Day7 =
        sum(
          !is.na(
            bray_increment_day7_minus_day3
          )
        ),
      median_Day7_minus_Day3_Bray_increment =
        median(
          bray_increment_day7_minus_day3,
          na.rm = TRUE
        ),
      mean_Day7_minus_Day3_Bray_increment =
        mean(
          bray_increment_day7_minus_day3,
          na.rm = TRUE
        ),
      n_positive_Day7_minus_Day3 =
        sum(
          bray_increment_day7_minus_day3 >
            0,
          na.rm = TRUE
        ),
      .groups =
        "drop"
    )

  write_csv_safe(
    core_group_summary,
    file.path(
      DIR_CORE,
      "V2_STEP91B4_PRJNA691455_source_group_descriptive_summary.csv"
    )
  )

  ck(
    "CORE SOURCE ANALYSIS COMPLETE"
  )

  # ----------------------------------------------------------
  # 3. PRJEB82425 external coarse infection-phenotype support
  # ----------------------------------------------------------

  # First validate the FULL primary cohort using the Step87B object.
  # Published grouping is 17 VAP/Pneumonia, 12 Other infection,
  # and 9 Control = 38 patients total.
  rr_ext <- registry |>
    filter(
      project ==
        "PRJEB82425"
    )

  if (
    nrow(rr_ext) != 1
  ) {
    stop(
      "PRJEB82425 Step87B registry row missing/non-unique."
    )
  }

  ext_object_path <- as.character(
    rr_ext$analysis_object_path[1]
  )

  if (
    !file.exists(
      ext_object_path
    )
  ) {
    stop(
      "PRJEB82425 Step87B analysis object missing: ",
      ext_object_path
    )
  }

  ext_obj <- readRDS(
    ext_object_path
  )

  ext_full_md <- as_tibble(
    ext_obj$metadata
  ) |>
    mutate(
      infection_source_proxy =
        case_when(
          phenotype ==
            "Pneumonia" ~
            "PULMONARY_VAP",

          phenotype ==
            "Other_infection" ~
            "OTHER_INFECTION_MIXED",

          phenotype ==
            "Control" ~
            "UNINFECTED_CONTROL",

          TRUE ~
            NA_character_
        )
    )

  ext_full_patient_groups <- ext_full_md |>
    distinct(
      patient_id,
      phenotype,
      infection_source_proxy
    )

  ext_full_counts <- ext_full_patient_groups |>
    count(
      infection_source_proxy,
      name =
        "n_full_primary_patients"
    )

  expected_full <- tibble(
    infection_source_proxy =
      c(
        "PULMONARY_VAP",
        "OTHER_INFECTION_MIXED",
        "UNINFECTED_CONTROL"
      ),
    n_expected_full =
      c(
        17L,
        12L,
        9L
      )
  )

  full_guard <- expected_full |>
    left_join(
      ext_full_counts,
      by =
        "infection_source_proxy"
    )

  if (
    any(
      full_guard$n_expected_full !=
        full_guard$n_full_primary_patients
    )
  ) {
    write_csv_safe(
      full_guard,
      file.path(
        DIR_EXT,
        "V2_STEP91B4_PRJEB82425_FAILED_FULL_COHORT_group_count_guard.csv"
      )
    )

    stop(
      "PRJEB82425 full Step87B cohort differs from validated 17/12/9 grouping."
    )
  }

  write_csv_safe(
    ext_full_patient_groups,
    file.path(
      DIR_EXT,
      "V2_STEP91B4_PRJEB82425_FULL_primary_patient_coarse_infection_group.csv"
    )
  )

  # Step88A2 displacement is intentionally restricted to true
  # longitudinal GE2 patients. Therefore it contains 35/38 patients:
  # VAP 16, Other infection 11, Control 8. One single-timepoint
  # patient from each published group is excluded from the
  # longitudinal branch.
  ext <- disp |>
    filter(
      project ==
        "PRJEB82425"
    ) |>
    mutate(
      infection_source_proxy =
        case_when(
          phenotype ==
            "Pneumonia" ~
            "PULMONARY_VAP",

          phenotype ==
            "Other_infection" ~
            "OTHER_INFECTION_MIXED",

          phenotype ==
            "Control" ~
            "UNINFECTED_CONTROL",

          TRUE ~
            NA_character_
        )
    )

  ext_patient_groups <- ext |>
    distinct(
      patient_id,
      phenotype,
      infection_source_proxy
    )

  ext_counts <- ext_patient_groups |>
    count(
      infection_source_proxy,
      name =
        "n_longitudinal_GE2_patients"
    )

  expected_longitudinal <- tibble(
    infection_source_proxy =
      c(
        "PULMONARY_VAP",
        "OTHER_INFECTION_MIXED",
        "UNINFECTED_CONTROL"
      ),
    n_expected_longitudinal_GE2 =
      c(
        16L,
        11L,
        8L
      )
  )

  longitudinal_guard <- expected_longitudinal |>
    left_join(
      ext_counts,
      by =
        "infection_source_proxy"
    )

  if (
    any(
      longitudinal_guard$n_expected_longitudinal_GE2 !=
        longitudinal_guard$n_longitudinal_GE2_patients
    ) ||
    sum(
      longitudinal_guard$n_longitudinal_GE2_patients
    ) != 35
  ) {
    write_csv_safe(
      longitudinal_guard,
      file.path(
        DIR_EXT,
        "V2_STEP91B4_PRJEB82425_FAILED_LONGITUDINAL_group_count_guard.csv"
      )
    )

    stop(
      paste0(
        "PRJEB82425 Step88A2 longitudinal GE2 counts differ from ",
        "expected 16/11/8 (35 total)."
      )
    )
  }

  attrition <- expected_full |>
    left_join(
      ext_full_counts,
      by =
        "infection_source_proxy"
    ) |>
    left_join(
      expected_longitudinal,
      by =
        "infection_source_proxy"
    ) |>
    left_join(
      ext_counts,
      by =
        "infection_source_proxy"
    ) |>
    mutate(
      n_excluded_from_longitudinal_GE2 =
        n_full_primary_patients -
          n_longitudinal_GE2_patients
    )

  if (
    any(
      attrition$n_excluded_from_longitudinal_GE2 != 1
    )
  ) {
    write_csv_safe(
      attrition,
      file.path(
        DIR_EXT,
        "V2_STEP91B4_PRJEB82425_FAILED_attrition_guard.csv"
      )
    )

    stop(
      "PRJEB82425 expected exactly one single-timepoint exclusion per group."
    )
  }

  write_csv_safe(
    ext_patient_groups,
    file.path(
      DIR_EXT,
      "V2_STEP91B4_PRJEB82425_LONGITUDINAL_GE2_patient_coarse_infection_group.csv"
    )
  )

  write_csv_safe(
    attrition,
    file.path(
      DIR_EXT,
      "V2_STEP91B4_PRJEB82425_full_to_longitudinal_attrition_by_group.csv"
    )
  )

  # Published aggregate composition of Other Infection group.
  ext_other_infection_aggregate <- tribble(
    ~published_subtype,
    ~n_patients,
    ~patient_level_mapping_publicly_available,

    "Fever of unknown origin",
    9L,
    FALSE,

    "Urinary tract infection",
    1L,
    FALSE,

    "Primary bacteremia/catheter-related infection",
    2L,
    FALSE
  ) |>
    mutate(
      source =
        "Kritikos et al. BMC Infect Dis. 2025;25:468, Supplementary Table S1",
      doi =
        "10.1186/s12879-025-10825-6"
    )

  write_csv_safe(
    ext_other_infection_aggregate,
    file.path(
      DIR_EXT,
      "V2_STEP91B4_PRJEB82425_other_infection_subtypes_AGGREGATE_ONLY.csv"
    )
  )

  # Inclusion-anchored patients are the cleanest subset for
  # patient-specific displacement. Among infected patients,
  # only 6 have both Infection_D1 and Infection_D5:
  # 3 VAP, 3 Other infection. This remains descriptive only.
  inclusion_anchor_ids <- ext |>
    filter(
      is_patient_reference,
      time_factor ==
        "Inclusion"
    ) |>
    pull(
      patient_id
    ) |>
    unique()

  ext_pair <- ext |>
    filter(
      patient_id %in%
        inclusion_anchor_ids,
      infection_source_proxy %in%
        c(
          "PULMONARY_VAP",
          "OTHER_INFECTION_MIXED"
        ),
      time_factor %in%
        c(
          "Infection_D1",
          "Infection_D5"
        )
    ) |>
    select(
      patient_id,
      infection_source_proxy,
      time_factor,
      bray_from_patient_baseline
    ) |>
    distinct() |>
    pivot_wider(
      names_from =
        time_factor,
      values_from =
        bray_from_patient_baseline
    ) |>
    filter(
      !is.na(
        Infection_D1
      ),
      !is.na(
        Infection_D5
      )
    ) |>
    mutate(
      bray_increment_D5_minus_D1 =
        Infection_D5 -
          Infection_D1
    )

  ext_pair_counts <- ext_pair |>
    count(
      infection_source_proxy
    )

  if (
    nrow(ext_pair) != 6 ||
    !(
      ext_pair_counts$n[
        ext_pair_counts$infection_source_proxy ==
          "PULMONARY_VAP"
      ] == 3 &&
      ext_pair_counts$n[
        ext_pair_counts$infection_source_proxy ==
          "OTHER_INFECTION_MIXED"
      ] == 3
    )
  ) {
    stop(
      "PRJEB82425 clean Inclusion-anchor D1-D5 subset expected 3 VAP + 3 Other infection."
    )
  }

  write_csv_safe(
    ext_pair,
    file.path(
      DIR_EXT,
      "V2_STEP91B4_PRJEB82425_Inclusion_anchor_D1_D5_DESCRIPTIVE_ONLY.csv"
    )
  )

  ext_pair_summary <- ext_pair |>
    group_by(
      infection_source_proxy
    ) |>
    summarise(
      n =
        n(),
      median_Infection_D1_Bray =
        median(
          Infection_D1
        ),
      median_Infection_D5_Bray =
        median(
          Infection_D5
        ),
      mean_increment =
        mean(
          bray_increment_D5_minus_D1
        ),
      median_increment =
        median(
          bray_increment_D5_minus_D1
        ),
      n_positive_increment =
        sum(
          bray_increment_D5_minus_D1 >
            0
        ),
      .groups =
        "drop"
    )

  write_csv_safe(
    ext_pair_summary,
    file.path(
      DIR_EXT,
      "V2_STEP91B4_PRJEB82425_Inclusion_anchor_D1_D5_DESCRIPTIVE_summary.csv"
    )
  )

  ck(
    "EXTERNAL SUPPORT COMPLETE"
  )

  # ----------------------------------------------------------
  # 4. Integrated interpretation
  # ----------------------------------------------------------

  primary_row <- core_tests |>
    filter(
      interpretation_tier ==
        "PRIMARY_SOURCE_MODIFICATION_ENDPOINT"
    )

  secondary_row <- core_tests |>
    filter(
      interpretation_tier ==
        "SECONDARY_SOURCE_MODIFICATION_ENDPOINT"
    )

  primary_signal <- if (
    primary_row$exact_permutation_p <
      0.05
  ) {
    "SOURCE_ASSOCIATED"
  } else {
    "NO_CLEAR_SOURCE_ASSOCIATION"
  }

  secondary_signal <- if (
    secondary_row$exact_permutation_p <
      0.05
  ) {
    "SOURCE_ASSOCIATED"
  } else {
    "NO_CLEAR_SOURCE_ASSOCIATION"
  }

  interpretation <- tibble(
    core_source_mapping =
      "HIGH_CONFIDENCE_PUBLISHED_PATIENT_LEVEL",
    core_groups =
      "PULMONARY_5_VS_ABDOMINAL_5",
    primary_endpoint =
      "Day7 Bray distance from patient-specific Day1 baseline",
    primary_exact_permutation_p =
      primary_row$exact_permutation_p,
    primary_source_signal =
      primary_signal,
    secondary_endpoint =
      "Day7 minus Day3 Bray increment",
    secondary_exact_permutation_p =
      secondary_row$exact_permutation_p,
    secondary_source_signal =
      secondary_signal,
    external_PRJEB82425_role =
      "DESCRIPTIVE_COARSE_SUPPORT_ONLY",
    source_branch_conclusion =
      case_when(
        primary_signal ==
          "NO_CLEAR_SOURCE_ASSOCIATION" &
          secondary_signal ==
            "NO_CLEAR_SOURCE_ASSOCIATION" ~
          "NO_EVIDENCE_THAT_CORE_PROGRESSIVE_DISPLACEMENT_IS_STRONGLY_DEPENDENT_ON_PULMONARY_VS_ABDOMINAL_SOURCE; POWER_IS VERY LIMITED",

        TRUE ~
          "POSSIBLE_SOURCE_MODIFICATION_REQUIRES_CAUTIOUS_REVIEW"
      )
  )

  write_csv_safe(
    interpretation,
    file.path(
      DIR_SYN,
      "V2_STEP91B4_integrated_source_interpretation.csv"
    )
  )

  readme <- c(
    "SEPSIS V2 - STEP91B4 INFECTION SOURCE FREEZE AND ANALYSIS",
    paste0(
      "Created: ",
      Sys.time()
    ),
    "",
    "CORE SOURCE RESCUE",
    "Step91A found no defensible infection-source field for core analysis in the frozen metadata.",
    "A pre-existing generic infection_source_standard column is retained as infection_source_standard_frozen_pre91B for audit; the literature-validated mapping is the analysis field.",
    "Manual source-paper validation recovered patient-level primary infection site for all 10 PRJNA691455 sepsis patients from published Table 2.",
    "",
    "CORE MAPPING",
    "Pulmonary: S1, S3, S4, S6, S7 (n=5).",
    "Abdominal: S2, S5, S8, S9, S10 (n=5).",
    "Citation: Yang XJ et al. World J Gastroenterol. 2021;27(19):2376-2393. DOI 10.3748/wjg.v27.i19.2376.",
    "",
    "PRIMARY SOURCE-MODIFICATION ENDPOINT",
    "Day7 Bray distance from each patient's own Day1 baseline, Pulmonary vs Abdominal.",
    "Exact label-permutation testing is used because n=5 vs 5.",
    "",
    "SECONDARY",
    "Day7-Day3 incremental Bray displacement among 9 complete patients (Pulmonary n=4, Abdominal n=5).",
    "",
    "PRJEB82425",
    "Published/full Step87B patient grouping is validated: VAP n=17, Other infection n=12, Control n=9.",
    "The Step88A2 longitudinal GE2 subset contains VAP n=16, Other infection n=11, Control n=8; one single-timepoint patient per group is excluded.",
    "Supplement Table S1 shows aggregate Other infection composition: FUO 9, UTI 1, bacteremia/catheter infection 2.",
    "These subtypes cannot be assigned to individual public patient IDs, so no fabricated detailed source mapping is created.",
    "The clean Inclusion-anchor Infection_D1-to-D5 subset has only 3 VAP and 3 Other-infection patients and is descriptive only.",
    "",
    "INTERPRETATION",
    interpretation$source_branch_conclusion,
    "",
    "IMPORTANT LIMITATION",
    "A nonsignificant pulmonary-vs-abdominal test in 10 patients is not evidence of equivalence.",
    "The source analysis is supportive/exploratory and should not replace the main patient-anchored ecological-displacement result."
  )

  writeLines(
    readme,
    file.path(
      OUT,
      "README_STEP91B4.txt"
    ),
    useBytes = TRUE
  )

  writeLines(
    c(
      paste0(
        "Completed: ",
        Sys.time()
      ),
      "Status: STEP91B4 COMPLETE",
      "PRJNA691455 source mapping frozen from published Table 2.",
      "Core source groups: Pulmonary n=5 / Abdominal n=5.",
      "Primary exact small-sample source comparison completed.",
      "Secondary Day7-Day3 source comparison completed.",
      "PRJEB82425 full grouping (17/12/9) and longitudinal GE2 attrition (16/11/8) validated; detailed Other-infection subtypes remain aggregate-only.",
      paste0(
        "Source branch conclusion: ",
        interpretation$source_branch_conclusion
      )
    ),
    file.path(
      OUT,
      "_STEP91B4_COMPLETE.ok"
    ),
    useBytes = TRUE
  )

  ck(
    "STEP91B4 COMPLETE"
  )

  cat(
    "\n============================================================\n"
  )

  cat(
    "SEPSIS V2 - STEP91B4 COMPLETE\n"
  )

  cat(
    "============================================================\n\n"
  )

  cat(
    "Core source mapping:\n"
  )

  print(
    source_map,
    n = Inf,
    width = Inf
  )

  cat(
    "\nCore source tests:\n"
  )

  print(
    core_tests,
    n = Inf,
    width = Inf
  )

  cat(
    "\nPRJEB82425 descriptive paired subset:\n"
  )

  print(
    ext_pair_summary,
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
        "STEP91B4 FATAL ERROR: ",
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
      "STEP91B4 FAILED"
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
