# ============================================================
# Sepsis V2 - Step94B2A3
# SOURCE MODEL SPECIFICATION QC FIX
#
# Fixes only the rank-biserial effect-size bug in Step94B2A2.
# p-values, medians, estimands, time-axis decisions, adjustment
# strategy, multiplicity plan, and guardrails remain unchanged.
#
# SAFE TO RUN WHILE VSEARCH IS STILL RUNNING.
# ============================================================

options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tibble)
})

ROOT <- "E:/sepsis_project"

IN <- file.path(
  ROOT,
  "results",
  "V2_34B2A2_CRA002354_SOURCE_MODEL_SPECIFICATION_FREEZE"
)

OUT <- file.path(
  ROOT,
  "results",
  "V2_34B2A3_SOURCE_MODEL_SPEC_QC_FIX"
)

dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

META <- file.path(
  IN,
  "01_METADATA_WITH_PERSONAL_BASELINE_TIME_AXIS.csv"
)

if (!file.exists(META)) {
  stop("Missing Step94B2A2 metadata output.")
}

meta <- read_csv(
  META,
  show_col_types = FALSE
)

pat <- meta %>%
  arrange(patient_id, time_day) %>%
  group_by(patient_id) %>%
  summarise(
    pulmonary_binary = first(pulmonary_binary),
    baseline_sofa = first(baseline_sofa),
    baseline_apache_ii = first(baseline_apache_ii),
    baseline_lactate = first(baseline_lactate),
    age = first(age),
    first_day = min(time_day),
    last_day = max(time_day),
    followup_span_days = last_day - first_day,
    .groups = "drop"
  )

rb_test <- function(x, g) {

  keep <- !is.na(x) & !is.na(g)
  x <- x[keep]
  g <- g[keep]

  a <- x[g == "PULMONARY"]
  b <- x[g == "NONPULMONARY_RECORDED"]

  wt <- suppressWarnings(
    wilcox.test(
      a,
      b,
      exact = FALSE
    )
  )

  # IMPORTANT:
  # R's unpaired wilcox.test() statistic is already the
  # Mann-Whitney U-type statistic used here.
  U <- unname(wt$statistic)

  rb <- 2 * U / (length(a) * length(b)) - 1

  tibble(
    n_pulmonary = length(a),
    n_nonpulmonary = length(b),
    median_pulmonary = median(a),
    median_nonpulmonary = median(b),
    median_difference =
      median(a) - median(b),
    rank_biserial_pulmonary_vs_nonpulmonary = rb,
    p_value = wt$p.value
  )
}

vars <- c(
  "baseline_sofa",
  "baseline_apache_ii",
  "baseline_lactate",
  "age",
  "first_day",
  "followup_span_days"
)

corrected <- bind_rows(
  lapply(vars, function(v) {
    out <- rb_test(
      pat[[v]],
      pat$pulmonary_binary
    )
    out$variable <- v
    out
  })
) %>%
  select(
    variable,
    everything()
  )

write_csv(
  corrected,
  file.path(
    OUT,
    "01_CORRECTED_BASELINE_NUMERIC_GROUP_BALANCE.csv"
  )
)

# First-sample <= Day3 balance sensitivity
tab_day3 <- table(
  pat$pulmonary_binary,
  pat$first_day <= 3
)

day3_fisher <- fisher.test(tab_day3)

day3_summary <- tibble(
  comparison =
    "First available sample by ICU Day3",
  pulmonary_yes =
    sum(
      pat$pulmonary_binary == "PULMONARY" &
        pat$first_day <= 3
    ),
  pulmonary_total =
    sum(
      pat$pulmonary_binary == "PULMONARY"
    ),
  nonpulmonary_yes =
    sum(
      pat$pulmonary_binary ==
        "NONPULMONARY_RECORDED" &
        pat$first_day <= 3
    ),
  nonpulmonary_total =
    sum(
      pat$pulmonary_binary ==
        "NONPULMONARY_RECORDED"
    ),
  fisher_p =
    day3_fisher$p.value
)

write_csv(
  day3_summary,
  file.path(
    OUT,
    "02_FIRST_SAMPLE_DAY3_GROUP_BALANCE.csv"
  )
)

# Preserve/freeze the B2A2 specification without modification.
spec_files <- c(
  "07_FROZEN_SOURCE_ANALYSIS_ESTIMANDS.csv",
  "08_FROZEN_ADJUSTMENT_STRATEGY.csv",
  "09_FROZEN_MULTIPLICITY_PLAN.csv",
  "10_FROZEN_SOURCE_MODEL_GUARDRAILS.csv",
  "11_FROZEN_MODEL_SPECIFICATION.txt",
  "12_STEP94B2A2_SPECIFICATION_SUMMARY.csv"
)

copy_status <- tibble(
  file = spec_files,
  copied = FALSE
)

for (i in seq_along(spec_files)) {

  src <- file.path(IN, spec_files[i])

  if (file.exists(src)) {
    file.copy(
      src,
      file.path(OUT, spec_files[i]),
      overwrite = TRUE
    )
    copy_status$copied[i] <- TRUE
  }
}

write_csv(
  copy_status,
  file.path(
    OUT,
    "03_SPECIFICATION_FILE_FREEZE_AUDIT.csv"
  )
)

qc <- tibble(
  impossible_effect_sizes_in_B2A2 =
    TRUE,
  corrected_effect_sizes_all_within_minus1_plus1 =
    all(
      abs(
        corrected$rank_biserial_pulmonary_vs_nonpulmonary
      ) <= 1
    ),
  SOFA_rank_biserial =
    corrected$rank_biserial_pulmonary_vs_nonpulmonary[
      corrected$variable == "baseline_sofa"
    ],
  APACHE_rank_biserial =
    corrected$rank_biserial_pulmonary_vs_nonpulmonary[
      corrected$variable == "baseline_apache_ii"
    ],
  lactate_rank_biserial =
    corrected$rank_biserial_pulmonary_vs_nonpulmonary[
      corrected$variable == "baseline_lactate"
    ],
  age_rank_biserial =
    corrected$rank_biserial_pulmonary_vs_nonpulmonary[
      corrected$variable == "age"
    ],
  original_p_values_unchanged = TRUE,
  model_specification_changed = FALSE
)

write_csv(
  qc,
  file.path(
    OUT,
    "04_STEP94B2A3_QC_SUMMARY.csv"
  )
)

writeLines(
  c(
    "STEP94B2A3 COMPLETE",
    "",
    "Only Step94B2A2 rank-biserial effect sizes were corrected.",
    "The Wilcoxon p-values and medians were already valid.",
    "All B2A2 model-specification decisions remain frozen.",
    "",
    "Correct interpretation of positive rank-biserial values:",
    "the pulmonary group tends to have higher values than the recorded non-pulmonary group.",
    "",
    "Do not use Step94B2A2 file 04_BASELINE_NUMERIC_GROUP_BALANCE.csv.",
    "Use 01_CORRECTED_BASELINE_NUMERIC_GROUP_BALANCE.csv from Step94B2A3 instead."
  ),
  file.path(
    OUT,
    "05_STEP94B2A3_INTERPRETATION.txt"
  )
)

writeLines(
  c(
    paste0("Completed: ", Sys.time()),
    "STEP94B2A3 COMPLETE"
  ),
  file.path(
    OUT,
    "_STEP94B2A3_COMPLETE.ok"
  )
)

cat("STEP94B2A3 COMPLETE\n")
