# ============================================================
# Sepsis V2 - Step92B
# Ecological Instability Index (EII) construction
#
# Goal:
# Build a microbiome-only ecological instability metric.
#
# Components:
# 1. Bray-Curtis displacement from patient baseline
# 2. Alpha diversity instability
# 3. Community turnover proxy
#
# This step does NOT use clinical outcomes.
# ============================================================

options(stringsAsFactors = FALSE)

pkgs <- c(
  "readr",
  "dplyr",
  "tidyr",
  "stringr"
)

missing <- pkgs[
  !vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)
]

if(length(missing)){
  install.packages(
    missing,
    repos="https://cloud.r-project.org"
  )
}

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(stringr)
})

ROOT <- "E:/sepsis_project"

STEP88 <- file.path(
  ROOT,
  "results",
  "V2_28A2_STEP88A_LONGITUDINAL_DIVERSITY_AND_DISPLACEMENT_FIXED"
)

OUT <- file.path(
  ROOT,
  "results",
  "V2_32B_STEP92B_EII_CONSTRUCTION"
)

dir.create(
  OUT,
  recursive=TRUE,
  showWarnings=FALSE
)

input <- file.path(
  STEP88,
  "03_BETA_DISPLACEMENT",
  "V2_STEP88A2_ALL_within_patient_bray_displacement.csv"
)

if(!file.exists(input)){
  stop("Step88A2 displacement file missing.")
}

dat <- read_csv(
  input,
  show_col_types=FALSE
)

required <- c(
  "patient_id",
  "project",
  "time_factor",
  "bray_from_patient_baseline",
  "Shannon",
  "Simpson",
  "Observed_ASV"
)

missing_cols <- setdiff(
  required,
  names(dat)
)

if(length(missing_cols)>0){
  stop(
    paste(
      "Missing columns:",
      paste(missing_cols, collapse=", ")
    )
  )
}

# ------------------------------------------------------------
# Normalize time labels
# ------------------------------------------------------------

dat <- dat %>%
  mutate(
    time_numeric =
      case_when(
        str_detect(
          time_factor,
          "Day 1|D1|Inclusion"
        ) ~ 1,

        str_detect(
          time_factor,
          "Day 3|D3"
        ) ~ 3,

        str_detect(
          time_factor,
          "Day 5|D5"
        ) ~ 5,

        str_detect(
          time_factor,
          "Day 7|D7"
        ) ~ 7,

        TRUE ~ NA_real_
      )
  )

if(any(is.na(dat$time_numeric))){
  warning("Some time labels could not be normalized.")
}

# ------------------------------------------------------------
# Calculate instability components
# ------------------------------------------------------------

eii <- dat %>%
  arrange(
    project,
    patient_id,
    time_numeric
  ) %>%
  group_by(
    project,
    patient_id
  ) %>%
  mutate(

    # existing ecological displacement
    bray_instability =
      bray_from_patient_baseline,

    # diversity departure from baseline
    shannon_instability =
      abs(
        Shannon -
          first(
            Shannon[
              time_numeric ==
                min(time_numeric, na.rm=TRUE)
            ]
          )
      ),

    simpson_instability =
      abs(
        Simpson -
          first(
            Simpson[
              time_numeric ==
                min(time_numeric, na.rm=TRUE)
            ]
          )
      ),

    observed_asv_instability =
      abs(
        Observed_ASV -
          first(
            Observed_ASV[
              time_numeric ==
                min(time_numeric, na.rm=TRUE)
            ]
          )
      )
  ) %>%
  ungroup()

# ------------------------------------------------------------
# Standardize within whole analysis object
# ------------------------------------------------------------

z <- function(x){
  if(sd(x, na.rm=TRUE)==0){
    return(rep(0,length(x)))
  }
  as.numeric(
    scale(x)
  )
}

eii <- eii %>%
  mutate(

    z_bray =
      z(
        bray_instability
      ),

    z_shannon =
      z(
        shannon_instability
      ),

    z_simpson =
      z(
        simpson_instability
      ),

    z_asv =
      z(
        observed_asv_instability
      ),

    EII_raw =
      rowMeans(
        cbind(
          z_bray,
          z_shannon,
          z_simpson,
          z_asv
        ),
        na.rm=TRUE
      )
  )

# Rescale 0-100 for interpretation

range_raw <- range(
  eii$EII_raw,
  na.rm=TRUE
)

eii <- eii %>%
  mutate(
    EII_0_100 =
      100 *
      (
        EII_raw -
          range_raw[1]
      ) /
      diff(range_raw)
  )

write_csv(
  eii,
  file.path(
    OUT,
    "V2_STEP92B_patient_timepoint_EII.csv"
  )
)

patient_summary <- eii %>%
  group_by(
    project,
    patient_id
  ) %>%
  summarise(
    n_timepoints =
      n(),
    baseline_EII =
      EII_0_100[
        which.min(time_numeric)
      ][1],
    max_EII =
      max(
        EII_0_100,
        na.rm=TRUE
      ),
    final_EII =
      EII_0_100[
        which.max(time_numeric)
      ][1],
    EII_change =
      final_EII -
      baseline_EII,
    .groups="drop"
  )

write_csv(
  patient_summary,
  file.path(
    OUT,
    "V2_STEP92B_patient_level_EII_summary.csv"
  )
)

cohort_summary <- patient_summary %>%
  group_by(
    project
  ) %>%
  summarise(
    n_patients=n(),
    mean_final_EII=
      mean(
        final_EII,
        na.rm=TRUE
      ),
    median_final_EII=
      median(
        final_EII,
        na.rm=TRUE
      ),
    mean_EII_change=
      mean(
        EII_change,
        na.rm=TRUE
      ),
    .groups="drop"
  )

write_csv(
  cohort_summary,
  file.path(
    OUT,
    "V2_STEP92B_cohort_EII_summary.csv"
  )
)

writeLines(
  c(
    paste0("Completed: ",Sys.time()),
    "STEP92B COMPLETE",
    "Ecological Instability Index constructed."
  ),
  file.path(
    OUT,
    "_STEP92B_COMPLETE.ok"
  )
)

cat("STEP92B COMPLETE\n")
