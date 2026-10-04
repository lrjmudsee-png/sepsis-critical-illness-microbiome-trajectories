# ============================================================
# Sepsis V2 - Step94B3B
# CRA002354 GENUS COVERAGE ROBUSTNESS SENSITIVITY
#
# Purpose:
# B3A is the authoritative exploratory genus analysis.
# This step does NOT replace its primary family.
#
# It checks whether B3A signals are robust after excluding samples
# with low sample-specific genus taxonomic coverage.
#
# Sensitivity thresholds:
# - genus-assigned fraction >= 0.70
# - genus-assigned fraction >= 0.80
#
# Same B3A eligible genus families are retained where analyzable.
# ============================================================

options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(stringr)
  library(tibble)
  library(lme4)
})

ROOT <- "E:/sepsis_project"

DATA_DIR <- file.path(
  ROOT, "data", "CRA002354", "03_vsearch97_silva1382"
)

B3A <- file.path(
  ROOT, "results",
  "V2_34B3A_CRA002354_EXPLORATORY_GENUS_SOURCE_ANALYSIS"
)

OUT <- file.path(
  ROOT, "results",
  "V2_34B3B_CRA002354_GENUS_COVERAGE_ROBUSTNESS"
)

dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(file.path(B3A, "_STEP94B3A_COMPLETE.ok"))) {
  stop("Step94B3A is not complete.")
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

COV_FILE <- file.path(
  B3A,
  "01_RAREFIED_GENUS_COVERAGE_QC.csv"
)

BASE_PRIMARY_FILE <- file.path(
  B3A,
  "04_EXPLORATORY_BASELINE_GENUS_SOURCE_CONTRAST.csv"
)

LONG_PRIMARY_FILE <- file.path(
  B3A,
  "06_EXPLORATORY_LONGITUDINAL_GENUS_SOURCE_TIME_INTERACTIONS.csv"
)

needed <- c(
  RARE4000_FILE,
  META_FILE,
  TAX_FILE,
  COV_FILE,
  BASE_PRIMARY_FILE,
  LONG_PRIMARY_FILE
)

missing <- needed[!file.exists(needed)]
if (length(missing) > 0) {
  stop(paste("Missing:", paste(missing, collapse = "; ")))
}

# ------------------------------------------------------------
# 1. Reconstruct B3A genus CLR matrix
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

cov <- read_csv(
  COV_FILE,
  show_col_types = FALSE
)

base_primary <- read_csv(
  BASE_PRIMARY_FILE,
  show_col_types = FALSE
)

long_primary <- read_csv(
  LONG_PRIMARY_FILE,
  show_col_types = FALSE
)

otu_cols <- setdiff(names(otu), "Run_ID")

Motu <- as.matrix(
  otu[, otu_cols, drop = FALSE]
)

storage.mode(Motu) <- "numeric"
rownames(Motu) <- otu$Run_ID

common <- Reduce(
  intersect,
  list(
    rownames(Motu),
    meta$Run_ID,
    cov$Run_ID
  )
)

meta <- meta %>%
  filter(Run_ID %in% common) %>%
  arrange(match(Run_ID, common)) %>%
  group_by(patient_id) %>%
  mutate(
    personal_baseline_day = min(time_day, na.rm = TRUE),
    days_since_personal_baseline = time_day - personal_baseline_day,
    is_personal_baseline = time_day == personal_baseline_day
  ) %>%
  ungroup()

Motu <- Motu[meta$Run_ID, , drop = FALSE]

cov <- cov %>%
  select(Run_ID, genus_assigned_fraction)

meta <- meta %>%
  left_join(cov, by = "Run_ID")

tax_map <- tax %>%
  transmute(
    otu_id = as.character(otu_id),
    Genus = str_trim(as.character(Genus))
  ) %>%
  filter(
    !is.na(Genus),
    Genus != "",
    !str_detect(
      tolower(Genus),
      "^unclassified$|^uncultured$|^unknown$"
    )
  ) %>%
  distinct(otu_id, .keep_all = TRUE)

otu_to_genus <- tax_map$Genus[
  match(colnames(Motu), tax_map$otu_id)
]

assigned <- !is.na(otu_to_genus)

Massigned <- Motu[, assigned, drop = FALSE]
genus_for_otu <- otu_to_genus[assigned]

genus_counts_t <- rowsum(
  t(Massigned),
  group = genus_for_otu,
  reorder = FALSE
)

Mgenus <- t(genus_counts_t)
rownames(Mgenus) <- rownames(Motu)

PSEUDO <- 0.5
logM <- log(Mgenus + PSEUDO)
Mclr <- logM - rowMeans(logM)

Mrel <- Mgenus / rowSums(Mgenus)

# ------------------------------------------------------------
# Helpers
# ------------------------------------------------------------

rank_biserial <- function(x, g) {

  keep <- !is.na(x) & !is.na(g)
  x <- x[keep]
  g <- g[keep]

  a <- x[g == "PULMONARY"]
  b <- x[g == "NONPULMONARY_RECORDED"]

  if (length(a) < 3 || length(b) < 3) {
    return(tibble(
      n_pulmonary = length(a),
      n_nonpulmonary = length(b),
      median_difference = NA_real_,
      rank_biserial = NA_real_,
      p_value = NA_real_,
      status = "INSUFFICIENT_GROUP_N"
    ))
  }

  wt <- suppressWarnings(
    wilcox.test(a, b, exact = FALSE)
  )

  U <- unname(wt$statistic)

  tibble(
    n_pulmonary = length(a),
    n_nonpulmonary = length(b),
    median_difference = median(a) - median(b),
    rank_biserial =
      2 * U / (length(a) * length(b)) - 1,
    p_value = wt$p.value,
    status = "OK"
  )
}

fit_interaction <- function(y, d, genus, threshold) {

  dd <- d
  dd$y <- y

  # Repeated patients only after coverage filtering.
  keep_patients <- dd %>%
    count(patient_id) %>%
    filter(n >= 2) %>%
    pull(patient_id)

  dd <- dd %>%
    filter(patient_id %in% keep_patients)

  if (
    nrow(dd) < 20 ||
    n_distinct(dd$patient_id) < 10 ||
    n_distinct(dd$pulmonary_binary) < 2
  ) {
    return(tibble(
      threshold = threshold,
      genus = genus,
      nobs = nrow(dd),
      npatients = n_distinct(dd$patient_id),
      beta = NA_real_,
      se = NA_real_,
      p_LRT = NA_real_,
      singular = NA,
      status = "INSUFFICIENT_DATA"
    ))
  }

  # Require >=3 patients in each source with genus presence after filtering.
  # Presence is assessed outside this function by caller.

  dd$pulmonary_binary <- factor(
    dd$pulmonary_binary,
    levels = c("NONPULMONARY_RECORDED", "PULMONARY")
  )

  full <- tryCatch(
    lmer(
      y ~
        days_since_personal_baseline *
        pulmonary_binary +
        (1|patient_id),
      data = dd,
      REML = FALSE,
      control = lmerControl(optimizer = "bobyqa")
    ),
    error = function(e) NULL
  )

  red <- tryCatch(
    lmer(
      y ~
        days_since_personal_baseline +
        pulmonary_binary +
        (1|patient_id),
      data = dd,
      REML = FALSE,
      control = lmerControl(optimizer = "bobyqa")
    ),
    error = function(e) NULL
  )

  if (is.null(full) || is.null(red)) {
    return(tibble(
      threshold = threshold,
      genus = genus,
      nobs = nrow(dd),
      npatients = n_distinct(dd$patient_id),
      beta = NA_real_,
      se = NA_real_,
      p_LRT = NA_real_,
      singular = NA,
      status = "MODEL_ERROR"
    ))
  }

  cn <- names(fixef(full))

  intn <- "days_since_personal_baseline:pulmonary_binaryPULMONARY"
  if (!intn %in% cn) {
    intn <- "pulmonary_binaryPULMONARY:days_since_personal_baseline"
  }

  beta <- if (intn %in% cn) unname(fixef(full)[intn]) else NA_real_
  se <- if (intn %in% cn) unname(sqrt(diag(vcov(full)))[intn]) else NA_real_

  p <- tryCatch(
    anova(red, full)$`Pr(>Chisq)`[2],
    error = function(e) NA_real_
  )

  tibble(
    threshold = threshold,
    genus = genus,
    nobs = nrow(dd),
    npatients = n_distinct(dd$patient_id),
    beta = beta,
    se = se,
    p_LRT = p,
    singular = isSingular(full, tol = 1e-4),
    status = "OK"
  )
}

# ------------------------------------------------------------
# 2. Coverage balance audit
# ------------------------------------------------------------

coverage_balance <- bind_rows(
  lapply(c(0.70, 0.80), function(th) {

    base_d <- meta %>%
      filter(is_personal_baseline)

    tibble(
      threshold = th,
      samples_total = nrow(meta),
      samples_retained = sum(meta$genus_assigned_fraction >= th),
      baseline_samples_total = nrow(base_d),
      baseline_samples_retained = sum(base_d$genus_assigned_fraction >= th),
      baseline_pulmonary_retained =
        sum(
          base_d$genus_assigned_fraction >= th &
          base_d$pulmonary_binary == "PULMONARY"
        ),
      baseline_nonpulmonary_retained =
        sum(
          base_d$genus_assigned_fraction >= th &
          base_d$pulmonary_binary == "NONPULMONARY_RECORDED"
        )
    )
  })
)

write_csv(
  coverage_balance,
  file.path(OUT, "01_COVERAGE_THRESHOLD_RETENTION.csv")
)

# ------------------------------------------------------------
# 3. Baseline sensitivity, same B3A genus family
# ------------------------------------------------------------

base_family <- base_primary$genus

baseline_sens_list <- list()

for (th in c(0.70, 0.80)) {

  idx <- which(
    meta$is_personal_baseline &
    meta$genus_assigned_fraction >= th
  )

  d <- meta[idx, , drop = FALSE]
  Mc <- Mclr[idx, , drop = FALSE]
  Mr <- Mrel[idx, , drop = FALSE]
  Mg <- Mgenus[idx, , drop = FALSE]

  res <- bind_rows(
    lapply(base_family, function(g) {

      if (!g %in% colnames(Mc)) {
        return(tibble(
          n_pulmonary = NA_integer_,
          n_nonpulmonary = NA_integer_,
          median_difference = NA_real_,
          rank_biserial = NA_real_,
          p_value = NA_real_,
          status = "GENUS_NOT_FOUND",
          genus = g,
          mean_rel_pulmonary = NA_real_,
          mean_rel_nonpulmonary = NA_real_
        ))
      }

      present <- Mg[, g] > 0

      npp <- n_distinct(
        d$patient_id[
          present &
          d$pulmonary_binary == "PULMONARY"
        ]
      )

      nnp <- n_distinct(
        d$patient_id[
          present &
          d$pulmonary_binary == "NONPULMONARY_RECORDED"
        ]
      )

      if (npp < 3 || nnp < 3) {
        return(tibble(
          n_pulmonary =
            sum(d$pulmonary_binary == "PULMONARY"),
          n_nonpulmonary =
            sum(d$pulmonary_binary == "NONPULMONARY_RECORDED"),
          median_difference = NA_real_,
          rank_biserial = NA_real_,
          p_value = NA_real_,
          status = "PREVALENCE_TOO_LOW_AFTER_FILTER",
          genus = g,
          mean_rel_pulmonary =
            mean(Mr[d$pulmonary_binary == "PULMONARY", g]),
          mean_rel_nonpulmonary =
            mean(Mr[d$pulmonary_binary == "NONPULMONARY_RECORDED", g])
        ))
      }

      rr <- rank_biserial(
        Mc[, g],
        d$pulmonary_binary
      )

      rr$genus <- g
      rr$mean_rel_pulmonary <-
        mean(Mr[d$pulmonary_binary == "PULMONARY", g])
      rr$mean_rel_nonpulmonary <-
        mean(Mr[d$pulmonary_binary == "NONPULMONARY_RECORDED", g])

      rr
    })
  ) %>%
    mutate(
      threshold = th,
      FDR = p.adjust(p_value, method = "BH")
    )

  baseline_sens_list[[as.character(th)]] <- res
}

baseline_sens <- bind_rows(baseline_sens_list)

write_csv(
  baseline_sens,
  file.path(OUT, "02_BASELINE_GENUS_COVERAGE_SENSITIVITY.csv")
)

# ------------------------------------------------------------
# 4. Longitudinal sensitivity, same B3A genus family
# ------------------------------------------------------------

long_family <- long_primary$genus
long_sens_list <- list()

for (th in c(0.70, 0.80)) {

  idx <- which(
    meta$genus_assigned_fraction >= th
  )

  d <- meta[idx, , drop = FALSE]
  Mc <- Mclr[idx, , drop = FALSE]
  Mg <- Mgenus[idx, , drop = FALSE]

  # Repeated patients after sample filtering.
  rep_ids <- d %>%
    count(patient_id) %>%
    filter(n >= 2) %>%
    pull(patient_id)

  keep_idx <- which(d$patient_id %in% rep_ids)
  d2 <- d[keep_idx, , drop = FALSE]
  Mc2 <- Mc[keep_idx, , drop = FALSE]
  Mg2 <- Mg[keep_idx, , drop = FALSE]

  res <- bind_rows(
    lapply(long_family, function(g) {

      if (!g %in% colnames(Mc2)) {
        return(tibble(
          threshold = th,
          genus = g,
          nobs = NA_integer_,
          npatients = NA_integer_,
          beta = NA_real_,
          se = NA_real_,
          p_LRT = NA_real_,
          singular = NA,
          status = "GENUS_NOT_FOUND"
        ))
      }

      present <- Mg2[, g] > 0

      npp <- n_distinct(
        d2$patient_id[
          present &
          d2$pulmonary_binary == "PULMONARY"
        ]
      )

      nnp <- n_distinct(
        d2$patient_id[
          present &
          d2$pulmonary_binary == "NONPULMONARY_RECORDED"
        ]
      )

      if (npp < 3 || nnp < 3) {
        return(tibble(
          threshold = th,
          genus = g,
          nobs = nrow(d2),
          npatients = n_distinct(d2$patient_id),
          beta = NA_real_,
          se = NA_real_,
          p_LRT = NA_real_,
          singular = NA,
          status = "PREVALENCE_TOO_LOW_AFTER_FILTER"
        ))
      }

      fit_interaction(
        Mc2[, g],
        d2,
        genus = g,
        threshold = th
      )
    })
  ) %>%
    mutate(
      FDR = p.adjust(p_LRT, method = "BH")
    )

  long_sens_list[[as.character(th)]] <- res
}

long_sens <- bind_rows(long_sens_list)

write_csv(
  long_sens,
  file.path(OUT, "03_LONGITUDINAL_GENUS_COVERAGE_SENSITIVITY.csv")
)

# ------------------------------------------------------------
# 5. Primary B3A hit robustness
# ------------------------------------------------------------

primary_hits <- base_primary %>%
  filter(FDR < 0.05) %>%
  transmute(
    genus,
    primary_rank_biserial = rank_biserial,
    primary_p = p_value,
    primary_FDR = FDR,
    primary_mean_rel_pulmonary = mean_rel_abundance_pulmonary,
    primary_mean_rel_nonpulmonary = mean_rel_abundance_nonpulmonary,
    max_primary_mean_relative_abundance =
      pmax(
        mean_rel_abundance_pulmonary,
        mean_rel_abundance_nonpulmonary
      ),
    low_abundance_both_groups_lt_0_1pct =
      max_primary_mean_relative_abundance < 0.001
  )

hit_rob <- primary_hits %>%
  left_join(
    baseline_sens %>%
      filter(threshold == 0.70) %>%
      select(
        genus,
        rb_70 = rank_biserial,
        p_70 = p_value,
        FDR_70 = FDR,
        status_70 = status
      ),
    by = "genus"
  ) %>%
  left_join(
    baseline_sens %>%
      filter(threshold == 0.80) %>%
      select(
        genus,
        rb_80 = rank_biserial,
        p_80 = p_value,
        FDR_80 = FDR,
        status_80 = status
      ),
    by = "genus"
  ) %>%
  mutate(
    direction_consistent_70 =
      sign(primary_rank_biserial) == sign(rb_70),
    direction_consistent_80 =
      sign(primary_rank_biserial) == sign(rb_80),
    robust_FDR_70 =
      !is.na(FDR_70) &
      FDR_70 < 0.05 &
      direction_consistent_70,
    robust_FDR_80 =
      !is.na(FDR_80) &
      FDR_80 < 0.05 &
      direction_consistent_80,
    robustness_class = case_when(
      robust_FDR_70 & robust_FDR_80 ~
        "ROBUST_AT_70_AND_80_PERCENT_COVERAGE",
      robust_FDR_70 | robust_FDR_80 ~
        "PARTIAL_COVERAGE_SENSITIVITY_SUPPORT",
      TRUE ~
        "NOT_FDR_ROBUST_TO_COVERAGE_FILTERING"
    )
  )

write_csv(
  hit_rob,
  file.path(OUT, "04_PRIMARY_BASELINE_HIT_COVERAGE_ROBUSTNESS.csv")
)

# ------------------------------------------------------------
# 6. Summary
# ------------------------------------------------------------

n_long70 <- sum(
  long_sens$threshold == 0.70 &
  long_sens$FDR < 0.05,
  na.rm = TRUE
)

n_long80 <- sum(
  long_sens$threshold == 0.80 &
  long_sens$FDR < 0.05,
  na.rm = TRUE
)

summary <- tibble(
  primary_B3A_baseline_FDR_hits = nrow(primary_hits),
  primary_B3A_longitudinal_FDR_hits =
    sum(long_primary$FDR < 0.05, na.rm = TRUE),
  baseline_hits_robust_at_70pct =
    sum(hit_rob$robust_FDR_70, na.rm = TRUE),
  baseline_hits_robust_at_80pct =
    sum(hit_rob$robust_FDR_80, na.rm = TRUE),
  baseline_hits_robust_at_both =
    sum(
      hit_rob$robust_FDR_70 &
      hit_rob$robust_FDR_80,
      na.rm = TRUE
    ),
  longitudinal_FDR_hits_at_70pct = n_long70,
  longitudinal_FDR_hits_at_80pct = n_long80,
  branch_role = "EXPLORATORY_SUPPORTIVE_ONLY"
)

write_csv(
  summary,
  file.path(OUT, "05_GENUS_COVERAGE_ROBUSTNESS_SUMMARY.csv")
)

writeLines(
  c(
    "STEP94B3B GENUS COVERAGE ROBUSTNESS",
    "",
    paste0(
      "Primary B3A baseline FDR hits: ",
      summary$primary_B3A_baseline_FDR_hits
    ),
    paste0(
      "Baseline hits robust at >=70% coverage: ",
      summary$baseline_hits_robust_at_70pct
    ),
    paste0(
      "Baseline hits robust at >=80% coverage: ",
      summary$baseline_hits_robust_at_80pct
    ),
    paste0(
      "Baseline hits robust at both thresholds: ",
      summary$baseline_hits_robust_at_both
    ),
    paste0(
      "Longitudinal FDR hits at >=70% coverage: ",
      n_long70
    ),
    paste0(
      "Longitudinal FDR hits at >=80% coverage: ",
      n_long80
    ),
    "",
    "Interpretation rule:",
    "B3A remains the primary exploratory genus family.",
    "Coverage-filtered results are QC sensitivities only.",
    "Low-abundance genera should be described cautiously even when statistically significant."
  ),
  file.path(OUT, "06_STEP94B3B_INTERPRETATION.txt")
)

writeLines(
  c(
    paste0("Completed: ", Sys.time()),
    "STEP94B3B COMPLETE"
  ),
  file.path(OUT, "_STEP94B3B_COMPLETE.ok")
)

cat("STEP94B3B COMPLETE\n")
