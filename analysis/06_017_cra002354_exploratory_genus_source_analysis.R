# ============================================================
# Sepsis V2 - Step94B3A
# CRA002354 EXPLORATORY GENUS-LEVEL INFECTION-SOURCE ANALYSIS
#
# Prerequisites:
# - Step94B1D corrected true-nonchimeric OTU ecology
# - Step94B1E taxonomy ID repair
# - Step94B2D infection-source ecological evidence freeze
#
# Role:
# EXPLORATORY / SUPPORTIVE ONLY.
# This step must not overwrite or upgrade the frozen ecological conclusion.
#
# Primary taxonomic strategy:
# 1) start from OTU97 counts rarefied to exactly 4000 reads;
# 2) map OTUs to genus using ID-fixed SILVA 138.2 taxonomy;
# 3) retain explicitly genus-assigned taxa only;
# 4) aggregate rarefied OTUs to genus;
# 5) CLR-transform genus counts with pseudocount 0.5;
# 6) perform baseline source contrasts;
# 7) perform longitudinal source x personal-time mixed models;
# 8) BH-FDR within each prespecified family;
# 9) sensitivity excluding OTHER_UNKNOWN;
# 10) pulmonary vs abdominal/GI is descriptive/exploratory only.
# ============================================================

options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(tibble)
  library(lme4)
})

ROOT <- "E:/sepsis_project"

DATA_DIR <- file.path(
  ROOT,
  "data",
  "CRA002354",
  "03_vsearch97_silva1382"
)

B1E_OUT <- file.path(
  ROOT,
  "results",
  "V2_34B1E_CRA002354_TAXONOMY_ID_REPAIR"
)

B2D_OUT <- file.path(
  ROOT,
  "results",
  "V2_34B2D_INFECTION_SOURCE_EVIDENCE_FREEZE"
)

OUT <- file.path(
  ROOT,
  "results",
  "V2_34B3A_CRA002354_EXPLORATORY_GENUS_SOURCE_ANALYSIS"
)

dir.create(
  OUT,
  recursive = TRUE,
  showWarnings = FALSE
)

# ------------------------------------------------------------
# 0. Hard gates
# ------------------------------------------------------------

TAX_READY <- file.path(
  B1E_OUT,
  "03_STEP94B1E_TAXONOMY_READINESS.csv"
)

if (!file.exists(TAX_READY)) {
  stop("Step94B1E taxonomy readiness file is missing.")
}

tr <- read_csv(
  TAX_READY,
  show_col_types = FALSE
)

if (
  !"taxonomy_ready_for_exploratory_genus_analysis" %in% names(tr) ||
  !isTRUE(
    tr$taxonomy_ready_for_exploratory_genus_analysis[1]
  )
) {
  stop("Taxonomy is not ready for exploratory genus analysis.")
}

if (
  !file.exists(
    file.path(
      B2D_OUT,
      "_STEP94B2D_COMPLETE.ok"
    )
  )
) {
  stop("Step94B2D infection-source ecological branch is not frozen.")
}

RARE4000_FILE <- file.path(
  DATA_DIR,
  "CRA002354_OTU97_rarefied_PRIMARY_4000.csv"
)

META_FILE <- file.path(
  DATA_DIR,
  "CRA002354_analysis_metadata_PRIMARY_NONCHIMERIC_min4000.csv"
)

TAX_FILE <- file.path(
  DATA_DIR,
  "CRA002354_OTU97_SILVA1382_taxonomy_minBoot80_ID_FIXED.csv"
)

needed <- c(
  RARE4000_FILE,
  META_FILE,
  TAX_FILE
)

missing <- needed[
  !file.exists(needed)
]

if (length(missing) > 0) {
  stop(
    paste0(
      "Missing required file(s): ",
      paste(
        missing,
        collapse = "; "
      )
    )
  )
}

# ------------------------------------------------------------
# 1. Load rarefied OTU table + fixed taxonomy + metadata
# ------------------------------------------------------------

otu <- read_csv(
  RARE4000_FILE,
  show_col_types = FALSE,
  name_repair = "minimal"
)

meta <- read_csv(
  META_FILE,
  show_col_types = FALSE
)

tax <- read_csv(
  TAX_FILE,
  show_col_types = FALSE
)

if (
  !"Run_ID" %in% names(otu) ||
  !"Run_ID" %in% names(meta) ||
  !"otu_id" %in% names(tax) ||
  !"Genus" %in% names(tax)
) {
  stop(
    "Required Run_ID / otu_id / Genus columns are missing."
  )
}

otu_cols <- setdiff(
  names(otu),
  "Run_ID"
)

Motu <- as.matrix(
  otu[
    ,
    otu_cols,
    drop = FALSE
  ]
)

storage.mode(Motu) <- "numeric"
rownames(Motu) <- otu$Run_ID

# Keep only samples with metadata.
common_samples <- intersect(
  rownames(Motu),
  meta$Run_ID
)

Motu <- Motu[
  common_samples,
  ,
  drop = FALSE
]

meta <- meta %>%
  filter(
    Run_ID %in%
      common_samples
  ) %>%
  arrange(
    match(
      Run_ID,
      common_samples
    )
  )

Motu <- Motu[
  meta$Run_ID,
  ,
  drop = FALSE
]

# Personal baseline time axis.
meta <- meta %>%
  group_by(patient_id) %>%
  mutate(
    personal_baseline_day =
      min(
        time_day,
        na.rm = TRUE
      ),
    days_since_personal_baseline =
      time_day -
      personal_baseline_day,
    is_personal_baseline =
      time_day ==
      personal_baseline_day
  ) %>%
  ungroup()

# ------------------------------------------------------------
# 2. Map OTUs to explicitly assigned genera
# ------------------------------------------------------------

tax_map <- tax %>%
  transmute(
    otu_id =
      as.character(otu_id),
    Genus =
      as.character(Genus)
  ) %>%
  mutate(
    Genus =
      str_trim(Genus),
    valid_genus =
      !is.na(Genus) &
      Genus != "" &
      !str_detect(
        tolower(Genus),
        "^unclassified$|^uncultured$|^unknown$"
      )
  ) %>%
  filter(valid_genus) %>%
  select(
    otu_id,
    Genus
  ) %>%
  distinct(
    otu_id,
    .keep_all = TRUE
  )

otu_to_genus <- tax_map$Genus[
  match(
    colnames(Motu),
    tax_map$otu_id
  )
]

assigned_otu <- !is.na(
  otu_to_genus
)

if (
  sum(assigned_otu) < 100
) {
  stop(
    "Too few explicitly genus-assigned OTUs."
  )
}

# Rarefied sample x assigned OTU.
Massigned <- Motu[
  ,
  assigned_otu,
  drop = FALSE
]

genus_for_otu <-
  otu_to_genus[
    assigned_otu
  ]

# Aggregate OTU x sample by genus.
genus_counts_t <- rowsum(
  t(Massigned),
  group =
    genus_for_otu,
  reorder = FALSE
)

Mgenus <- t(
  genus_counts_t
)

rownames(Mgenus) <-
  rownames(Massigned)

# ------------------------------------------------------------
# 3. Genus coverage QC
# ------------------------------------------------------------

assigned_depth <-
  rowSums(Mgenus)

coverage_qc <- tibble(
  Run_ID =
    rownames(Mgenus),
  rarefied_total_reads =
    rowSums(Motu),
  genus_assigned_reads =
    assigned_depth,
  genus_assigned_fraction =
    genus_assigned_reads /
    rarefied_total_reads
) %>%
  left_join(
    meta %>%
      select(
        Run_ID,
        patient_id,
        pulmonary_binary,
        infection_source_group,
        time_day
      ),
    by = "Run_ID"
  )

write_csv(
  coverage_qc,
  file.path(
    OUT,
    "01_RAREFIED_GENUS_COVERAGE_QC.csv"
  )
)

coverage_summary <- coverage_qc %>%
  summarise(
    samples = n(),
    min_fraction =
      min(
        genus_assigned_fraction
      ),
    median_fraction =
      median(
        genus_assigned_fraction
      ),
    mean_fraction =
      mean(
        genus_assigned_fraction
      ),
    max_fraction =
      max(
        genus_assigned_fraction
      ),
    samples_ge80pct =
      sum(
        genus_assigned_fraction >=
          0.80
      ),
    samples_ge70pct =
      sum(
        genus_assigned_fraction >=
          0.70
      )
  )

write_csv(
  coverage_summary,
  file.path(
    OUT,
    "02_RAREFIED_GENUS_COVERAGE_SUMMARY.csv"
  )
)

# ------------------------------------------------------------
# 4. CLR transformation
# ------------------------------------------------------------

PSEUDO <- 0.5

logM <- log(
  Mgenus +
    PSEUDO
)

Mclr <- logM -
  rowMeans(
    logM
  )

# Relative abundance only for descriptive summaries.
Mrel <- Mgenus /
  rowSums(Mgenus)

# ------------------------------------------------------------
# Helpers
# ------------------------------------------------------------

rank_biserial <- function(
  x,
  g
) {

  keep <- !is.na(x) &
    !is.na(g)

  x <- x[keep]
  g <- g[keep]

  a <- x[
    g ==
      "PULMONARY"
  ]

  b <- x[
    g ==
      "NONPULMONARY_RECORDED"
  ]

  if (
    length(a) == 0 ||
    length(b) == 0
  ) {
    return(
      tibble(
        n_pulmonary =
          length(a),
        n_nonpulmonary =
          length(b),
        median_pulmonary =
          NA_real_,
        median_nonpulmonary =
          NA_real_,
        median_difference =
          NA_real_,
        rank_biserial =
          NA_real_,
        p_value =
          NA_real_
      )
    )
  }

  wt <- suppressWarnings(
    wilcox.test(
      a,
      b,
      exact = FALSE
    )
  )

  U <- unname(
    wt$statistic
  )

  tibble(
    n_pulmonary =
      length(a),
    n_nonpulmonary =
      length(b),
    median_pulmonary =
      median(a),
    median_nonpulmonary =
      median(b),
    median_difference =
      median(a) -
      median(b),
    rank_biserial =
      2 * U /
      (
        length(a) *
        length(b)
      ) -
      1,
    p_value =
      wt$p.value
  )
}

fit_genus_interaction <- function(
  y,
  d,
  genus_name,
  label
) {

  dd <- d
  dd$y <- y

  dd <- dd %>%
    filter(
      !is.na(y),
      !is.na(
        days_since_personal_baseline
      ),
      !is.na(
        pulmonary_binary
      ),
      !is.na(
        patient_id
      )
    )

  if (
    nrow(dd) < 20 ||
    n_distinct(
      dd$patient_id
    ) < 10
  ) {
    return(
      tibble(
        analysis =
          label,
        genus =
          genus_name,
        nobs =
          nrow(dd),
        npatients =
          n_distinct(
            dd$patient_id
          ),
        interaction_beta =
          NA_real_,
        se =
          NA_real_,
        ci95_low =
          NA_real_,
        ci95_high =
          NA_real_,
        p_LRT =
          NA_real_,
        singular =
          NA,
        status =
          "INSUFFICIENT_DATA"
      )
    )
  }

  dd$pulmonary_binary <- factor(
    dd$pulmonary_binary,
    levels = c(
      "NONPULMONARY_RECORDED",
      "PULMONARY"
    )
  )

  full <- tryCatch(
    lmer(
      y ~
        days_since_personal_baseline *
        pulmonary_binary +
        (1|patient_id),
      data = dd,
      REML = FALSE,
      control =
        lmerControl(
          optimizer =
            "bobyqa"
        )
    ),
    error =
      function(e) NULL
  )

  red <- tryCatch(
    lmer(
      y ~
        days_since_personal_baseline +
        pulmonary_binary +
        (1|patient_id),
      data = dd,
      REML = FALSE,
      control =
        lmerControl(
          optimizer =
            "bobyqa"
        )
    ),
    error =
      function(e) NULL
  )

  if (
    is.null(full) ||
    is.null(red)
  ) {
    return(
      tibble(
        analysis =
          label,
        genus =
          genus_name,
        nobs =
          nrow(dd),
        npatients =
          n_distinct(
            dd$patient_id
          ),
        interaction_beta =
          NA_real_,
        se =
          NA_real_,
        ci95_low =
          NA_real_,
        ci95_high =
          NA_real_,
        p_LRT =
          NA_real_,
        singular =
          NA,
        status =
          "MODEL_ERROR"
      )
    )
  }

  int_name <-
    "days_since_personal_baseline:pulmonary_binaryPULMONARY"

  cn <- names(
    fixef(full)
  )

  if (
    !int_name %in% cn
  ) {
    alt <-
      "pulmonary_binaryPULMONARY:days_since_personal_baseline"

    if (
      alt %in% cn
    ) {
      int_name <- alt
    }
  }

  beta <- if (
    int_name %in% cn
  ) {
    unname(
      fixef(full)[
        int_name
      ]
    )
  } else {
    NA_real_
  }

  se <- if (
    int_name %in% cn
  ) {
    unname(
      sqrt(
        diag(
          vcov(full)
        )
      )[
        int_name
      ]
    )
  } else {
    NA_real_
  }

  p <- tryCatch(
    anova(
      red,
      full
    )$`Pr(>Chisq)`[2],
    error =
      function(e) NA_real_
  )

  tibble(
    analysis =
      label,
    genus =
      genus_name,
    nobs =
      nrow(dd),
    npatients =
      n_distinct(
        dd$patient_id
      ),
    interaction_beta =
      beta,
    se =
      se,
    ci95_low =
      beta -
      1.96 * se,
    ci95_high =
      beta +
      1.96 * se,
    p_LRT =
      p,
    singular =
      isSingular(
        full,
        tol = 1e-4
      ),
    status =
      "OK"
  )
}

# ------------------------------------------------------------
# 5. Baseline genus family
# ------------------------------------------------------------

base_idx <- which(
  meta$is_personal_baseline
)

base_meta <- meta[
  base_idx,
  ,
  drop = FALSE
]

Mbase_clr <- Mclr[
  base_idx,
  ,
  drop = FALSE
]

Mbase_rel <- Mrel[
  base_idx,
  ,
  drop = FALSE
]

# Prevalence = raw rarefied genus count >0.
Mbase_count <- Mgenus[
  base_idx,
  ,
  drop = FALSE
]

base_prev <- colMeans(
  Mbase_count > 0
)

base_patient_prev <- lapply(
  colnames(Mbase_count),
  function(g) {

    present <-
      Mbase_count[
        ,
        g
      ] > 0

    tibble(
      genus = g,
      prevalence =
        mean(present),
      pulmonary_patients_present =
        n_distinct(
          base_meta$patient_id[
            present &
            base_meta$pulmonary_binary ==
              "PULMONARY"
          ]
        ),
      nonpulmonary_patients_present =
        n_distinct(
          base_meta$patient_id[
            present &
            base_meta$pulmonary_binary ==
              "NONPULMONARY_RECORDED"
          ]
        )
    )
  }
) %>%
  bind_rows()

eligible_base <-
  base_patient_prev %>%
  filter(
    prevalence >=
      0.10,
    pulmonary_patients_present >=
      3,
    nonpulmonary_patients_present >=
      3
  )

write_csv(
  base_patient_prev,
  file.path(
    OUT,
    "03_BASELINE_GENUS_PREVALENCE_FILTER.csv"
  )
)

baseline_results <- bind_rows(
  lapply(
    eligible_base$genus,
    function(g) {

      rr <- rank_biserial(
        Mbase_clr[
          ,
          g
        ],
        base_meta$pulmonary_binary
      )

      rr$genus <- g

      rr$mean_rel_abundance_pulmonary <-
        mean(
          Mbase_rel[
            base_meta$pulmonary_binary ==
              "PULMONARY",
            g
          ]
        )

      rr$mean_rel_abundance_nonpulmonary <-
        mean(
          Mbase_rel[
            base_meta$pulmonary_binary ==
              "NONPULMONARY_RECORDED",
            g
          ]
        )

      rr
    }
  )
) %>%
  mutate(
    FDR =
      p.adjust(
        p_value,
        method = "BH"
      ),
    direction =
      case_when(
        median_difference > 0 ~
          "HIGHER_CLR_IN_PULMONARY",
        median_difference < 0 ~
          "LOWER_CLR_IN_PULMONARY",
        TRUE ~
          "NO_DIRECTION"
      )
  ) %>%
  arrange(
    FDR,
    p_value
  )

write_csv(
  baseline_results,
  file.path(
    OUT,
    "04_EXPLORATORY_BASELINE_GENUS_SOURCE_CONTRAST.csv"
  )
)

# ------------------------------------------------------------
# 6. Longitudinal genus source x personal-time family
#    Use repeated patients only.
# ------------------------------------------------------------

patient_n <- meta %>%
  count(
    patient_id,
    name = "n_samples"
  )

repeat_ids <- patient_n %>%
  filter(
    n_samples >= 2
  ) %>%
  pull(
    patient_id
  )

long_idx <- which(
  meta$patient_id %in%
    repeat_ids
)

long_meta <- meta[
  long_idx,
  ,
  drop = FALSE
]

Mlong_count <- Mgenus[
  long_idx,
  ,
  drop = FALSE
]

Mlong_clr <- Mclr[
  long_idx,
  ,
  drop = FALSE
]

long_prev <- lapply(
  colnames(Mlong_count),
  function(g) {

    present <-
      Mlong_count[
        ,
        g
      ] > 0

    tibble(
      genus = g,
      sample_prevalence =
        mean(present),
      pulmonary_patients_present =
        n_distinct(
          long_meta$patient_id[
            present &
            long_meta$pulmonary_binary ==
              "PULMONARY"
          ]
        ),
      nonpulmonary_patients_present =
        n_distinct(
          long_meta$patient_id[
            present &
            long_meta$pulmonary_binary ==
              "NONPULMONARY_RECORDED"
          ]
        )
    )
  }
) %>%
  bind_rows()

eligible_long <- long_prev %>%
  filter(
    sample_prevalence >=
      0.10,
    pulmonary_patients_present >=
      3,
    nonpulmonary_patients_present >=
      3
  )

write_csv(
  long_prev,
  file.path(
    OUT,
    "05_LONGITUDINAL_GENUS_PREVALENCE_FILTER.csv"
  )
)

long_results <- bind_rows(
  lapply(
    eligible_long$genus,
    function(g) {

      fit_genus_interaction(
        Mlong_clr[
          ,
          g
        ],
        long_meta,
        genus_name = g,
        label =
          "PRIMARY_REPEATED_PATIENT_CLR"
      )
    }
  )
) %>%
  mutate(
    FDR =
      p.adjust(
        p_LRT,
        method = "BH"
      ),
    direction =
      case_when(
        interaction_beta > 0 ~
          "MORE_POSITIVE_TIME_SLOPE_IN_PULMONARY",
        interaction_beta < 0 ~
          "MORE_NEGATIVE_TIME_SLOPE_IN_PULMONARY",
        TRUE ~
          "NO_DIRECTION"
      )
  ) %>%
  arrange(
    FDR,
    p_LRT
  )

write_csv(
  long_results,
  file.path(
    OUT,
    "06_EXPLORATORY_LONGITUDINAL_GENUS_SOURCE_TIME_INTERACTIONS.csv"
  )
)

# ------------------------------------------------------------
# 7. Sensitivity: exclude OTHER_UNKNOWN
# ------------------------------------------------------------

known_idx <- which(
  long_meta$infection_source_group !=
    "OTHER_UNKNOWN"
)

known_meta <- long_meta[
  known_idx,
  ,
  drop = FALSE
]

Mknown_clr <- Mlong_clr[
  known_idx,
  ,
  drop = FALSE
]

# Refit same primary eligible genera to preserve comparability.
known_results <- bind_rows(
  lapply(
    eligible_long$genus,
    function(g) {

      fit_genus_interaction(
        Mknown_clr[
          ,
          g
        ],
        known_meta,
        genus_name = g,
        label =
          "SENSITIVITY_EXCLUDE_OTHER_UNKNOWN"
      )
    }
  )
) %>%
  mutate(
    FDR =
      p.adjust(
        p_LRT,
        method = "BH"
      )
  ) %>%
  arrange(
    FDR,
    p_LRT
  )

write_csv(
  known_results,
  file.path(
    OUT,
    "07_SENSITIVITY_EXCLUDE_OTHER_UNKNOWN_GENUS_INTERACTIONS.csv"
  )
)

# ------------------------------------------------------------
# 8. Robustness table for primary longitudinal hits
# ------------------------------------------------------------

primary_hits <- long_results %>%
  filter(
    !is.na(FDR),
    FDR < 0.05
  ) %>%
  select(
    genus,
    primary_beta =
      interaction_beta,
    primary_p =
      p_LRT,
    primary_FDR =
      FDR
  )

hit_robustness <- primary_hits %>%
  left_join(
    known_results %>%
      select(
        genus,
        sensitivity_beta =
          interaction_beta,
        sensitivity_p =
          p_LRT,
        sensitivity_FDR =
          FDR
      ),
    by = "genus"
  ) %>%
  mutate(
    direction_consistent =
      sign(
        primary_beta
      ) ==
      sign(
        sensitivity_beta
      ),
    sensitivity_FDR_significant =
      sensitivity_FDR <
      0.05
  )

write_csv(
  hit_robustness,
  file.path(
    OUT,
    "08_PRIMARY_LONGITUDINAL_HIT_ROBUSTNESS.csv"
  )
)

# ------------------------------------------------------------
# 9. Descriptive pulmonary vs abdominal/GI baseline contrast
# ------------------------------------------------------------

abd_idx <- which(
  base_meta$infection_source_group %in%
    c(
      "RESPIRATORY",
      "ABDOMINAL_GI"
    )
)

abd_meta <- base_meta[
  abd_idx,
  ,
  drop = FALSE
]

Mabd_clr <- Mbase_clr[
  abd_idx,
  ,
  drop = FALSE
]

Mabd_rel <- Mbase_rel[
  abd_idx,
  ,
  drop = FALSE
]

# Re-label for a descriptive two-group comparison.
abd_group <- ifelse(
  abd_meta$infection_source_group ==
    "RESPIRATORY",
  "PULMONARY",
  "ABDOMINAL_GI"
)

abd_results <- bind_rows(
  lapply(
    eligible_base$genus,
    function(g) {

      a <-
        Mabd_clr[
          abd_group ==
            "PULMONARY",
          g
        ]

      b <-
        Mabd_clr[
          abd_group ==
            "ABDOMINAL_GI",
          g
        ]

      if (
        length(a) < 5 ||
        length(b) < 5
      ) {
        return(
          tibble(
            genus = g,
            n_pulmonary =
              length(a),
            n_abdominal_GI =
              length(b),
            median_CLR_pulmonary =
              NA_real_,
            median_CLR_abdominal_GI =
              NA_real_,
            median_difference =
              NA_real_,
            p_value =
              NA_real_
          )
        )
      }

      wt <- suppressWarnings(
        wilcox.test(
          a,
          b,
          exact = FALSE
        )
      )

      tibble(
        genus = g,
        n_pulmonary =
          length(a),
        n_abdominal_GI =
          length(b),
        median_CLR_pulmonary =
          median(a),
        median_CLR_abdominal_GI =
          median(b),
        median_difference =
          median(a) -
          median(b),
        p_value =
          wt$p.value
      )
    }
  )
) %>%
  mutate(
    FDR =
      p.adjust(
        p_value,
        method = "BH"
      )
  ) %>%
  arrange(
    FDR,
    p_value
  )

write_csv(
  abd_results,
  file.path(
    OUT,
    "09_DESCRIPTIVE_PULMONARY_VS_ABDOMINAL_GI_BASELINE_GENUS.csv"
  )
)

# ------------------------------------------------------------
# 10. Top-taxonomy summary
# ------------------------------------------------------------

top_baseline <- baseline_results %>%
  slice_head(
    n = 20
  ) %>%
  mutate(
    family =
      "BASELINE"
  )

top_long <- long_results %>%
  slice_head(
    n = 20
  ) %>%
  transmute(
    genus,
    interaction_beta,
    p_LRT,
    FDR,
    direction,
    family =
      "LONGITUDINAL"
  )

write_csv(
  top_baseline,
  file.path(
    OUT,
    "10_TOP20_BASELINE_GENUS_SIGNALS.csv"
  )
)

write_csv(
  top_long,
  file.path(
    OUT,
    "11_TOP20_LONGITUDINAL_GENUS_SIGNALS.csv"
  )
)

# ------------------------------------------------------------
# 11. Evidence classification
# ------------------------------------------------------------

n_base_hits <- sum(
  baseline_results$FDR <
    0.05,
  na.rm = TRUE
)

n_long_hits <- sum(
  long_results$FDR <
    0.05,
  na.rm = TRUE
)

n_robust_long <- sum(
  hit_robustness$
    sensitivity_FDR_significant &
    hit_robustness$
      direction_consistent,
  na.rm = TRUE
)

tier <- case_when(

  n_long_hits > 0 &
    n_robust_long > 0 ~
    "EXPLORATORY_LONGITUDINAL_GENUS_SIGNALS_WITH_SENSITIVITY_SUPPORT",

  n_long_hits > 0 ~
    "EXPLORATORY_LONGITUDINAL_GENUS_SIGNALS_WITHOUT_FULL_SENSITIVITY_SUPPORT",

  n_long_hits == 0 &
    n_base_hits > 0 ~
    "EXPLORATORY_BASELINE_GENUS_DIFFERENCES_ONLY",

  TRUE ~
    "NO_FDR_SIGNIFICANT_GENUS_LEVEL_SOURCE_SIGNAL"
)

evidence <- tibble(
  taxonomy_mapping_fraction =
    tr$taxonomy_rows_mapped_fraction[1],
  abundance_weighted_genus_assignment =
    tr$abundance_weighted_genus_assignment[1],
  baseline_genera_tested =
    nrow(
      baseline_results
    ),
  baseline_FDR_hits =
    n_base_hits,
  longitudinal_genera_tested =
    nrow(
      long_results
    ),
  longitudinal_FDR_hits =
    n_long_hits,
  longitudinal_hits_robust_to_excluding_OTHER_UNKNOWN =
    n_robust_long,
  evidence_tier =
    tier,
  branch_role =
    "EXPLORATORY_SUPPORTIVE_ONLY"
)

write_csv(
  evidence,
  file.path(
    OUT,
    "12_GENUS_SOURCE_EVIDENCE_SUMMARY.csv"
  )
)

# ------------------------------------------------------------
# 12. Simple diagnostic figures
# ------------------------------------------------------------

pdf(
  file.path(
    OUT,
    "Figure_STEP94B3A_genus_source_signals.pdf"
  ),
  width = 10,
  height = 5
)

par(
  mfrow = c(1,2),
  mar = c(4.5,4.5,3,1)
)

if (
  nrow(
    baseline_results
  ) > 0
) {

  plot(
    baseline_results$median_difference,
    -log10(
      pmax(
        baseline_results$p_value,
        1e-300
      )
    ),
    xlab =
      "Pulmonary - non-pulmonary median CLR",
    ylab =
      "-log10(p)",
    main =
      "A. Baseline genus contrasts",
    pch = 16
  )

  abline(
    h =
      -log10(0.05),
    lty = 2
  )

} else {
  plot.new()
  title(
    "A. No eligible baseline genera"
  )
}

if (
  nrow(
    long_results
  ) > 0
) {

  plot(
    long_results$interaction_beta,
    -log10(
      pmax(
        long_results$p_LRT,
        1e-300
      )
    ),
    xlab =
      "Source x time interaction beta (CLR)",
    ylab =
      "-log10(p)",
    main =
      "B. Longitudinal genus interactions",
    pch = 16
  )

  abline(
    h =
      -log10(0.05),
    lty = 2
  )

} else {
  plot.new()
  title(
    "B. No eligible longitudinal genera"
  )
}

dev.off()

# ------------------------------------------------------------
# 13. Frozen interpretation
# ------------------------------------------------------------

interpretation <- c(
  "STEP94B3A EXPLORATORY GENUS SOURCE ANALYSIS",
  "",
  paste0(
    "Taxonomy mapping fraction: ",
    signif(
      tr$taxonomy_rows_mapped_fraction[1],
      5
    )
  ),
  paste0(
    "Abundance-weighted genus assignment: ",
    signif(
      tr$abundance_weighted_genus_assignment[1],
      5
    )
  ),
  paste0(
    "Baseline genera tested after prevalence filtering: ",
    nrow(
      baseline_results
    )
  ),
  paste0(
    "Baseline FDR-significant genera: ",
    n_base_hits
  ),
  paste0(
    "Longitudinal genera tested after prevalence filtering: ",
    nrow(
      long_results
    )
  ),
  paste0(
    "Longitudinal FDR-significant source x time genera: ",
    n_long_hits
  ),
  paste0(
    "Longitudinal hits robust after excluding OTHER_UNKNOWN: ",
    n_robust_long
  ),
  paste0(
    "Evidence tier: ",
    tier
  ),
  "",
  "Guardrails:",
  "This branch is exploratory/supportive only.",
  "Do not reinterpret taxonomic hits as a reproducible sepsis-specific signature.",
  "Do not let genus-level sensitivity-only significance upgrade the frozen ecological source conclusion.",
  "Pulmonary vs abdominal/GI is descriptive because of the much smaller abdominal/GI sample size."
)

writeLines(
  interpretation,
  file.path(
    OUT,
    "13_STEP94B3A_INTERPRETATION.txt"
  )
)

writeLines(
  c(
    paste0(
      "Completed: ",
      Sys.time()
    ),
    paste0(
      "Evidence tier: ",
      tier
    ),
    "STEP94B3A COMPLETE"
  ),
  file.path(
    OUT,
    "_STEP94B3A_COMPLETE.ok"
  )
)

cat(
  "STEP94B3A COMPLETE\n"
)
