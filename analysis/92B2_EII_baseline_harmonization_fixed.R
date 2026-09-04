# ============================================================
# Sepsis V2 - Step92B2
# EII baseline harmonization fixed version
#
# Fix:
# Step92B assumed Day1/Day3/Day7 labels universally.
# This version first audits time labels and uses cohort-specific
# baseline anchoring.
#
# Output:
# 1. timepoint audit
# 2. baseline mapping
# 3. corrected EII table
# 4. attrition records
# ============================================================

options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(stringr)
  library(tidyr)
})

ROOT <- "E:/sepsis_project"

INPUT <- file.path(
  ROOT,
  "results",
  "V2_28A2_STEP88A_LONGITUDINAL_DIVERSITY_AND_DISPLACEMENT_FIXED",
  "03_BETA_DISPLACEMENT",
  "V2_STEP88A2_ALL_within_patient_bray_displacement.csv"
)

OUT <- file.path(
  ROOT,
  "results",
  "V2_32B2_STEP92B_BASELINE_HARMONIZATION_FIXED"
)

dir.create(
  OUT,
  recursive = TRUE,
  showWarnings = FALSE
)

if(!file.exists(INPUT)){
  stop("Missing Step88A2 input.")
}

dat <- read_csv(
  INPUT,
  show_col_types = FALSE
)

required <- c(
  "project",
  "patient_id",
  "time_factor",
  "bray_from_patient_baseline",
  "Shannon",
  "Simpson",
  "Observed_ASV"
)

miss <- setdiff(required,names(dat))

if(length(miss)>0){
  stop(
    paste(
      "Missing:",
      paste(miss,collapse=",")
    )
  )
}

# ------------------------------------------------------------
# Audit raw time labels
# ------------------------------------------------------------

time_audit <- dat %>%
  count(
    project,
    time_factor,
    name="n_samples"
  ) %>%
  arrange(
    project,
    time_factor
  )

write_csv(
  time_audit,
  file.path(
    OUT,
    "V2_STEP92B2_raw_timepoint_audit.csv"
  )
)

# ------------------------------------------------------------
# Cohort-specific baseline rules
# ------------------------------------------------------------

baseline_rules <- tribble(
  ~project, ~baseline_regex,

  "PRJNA691455",
  "Day 1|Day1|D1",

  "PRJEB82425",
  "Inclusion|Baseline",

  "PRJNA516701",
  "Baseline|baseline|T0|Day 0|Day0",

  "PRJNA851469",
  "Baseline|baseline|T0|Day 0|Day0",

  "PRJNA578267",
  "Baseline|baseline|T0|Day 0|Day0",

  "PRJNA430161",
  "Baseline|baseline|T0|Day 0|Day0",

  "PRJNA1166732",
  "Baseline|baseline|T0|Day 0|Day0",

  "PRJNA978257",
  "Baseline|baseline|T0|Day 0|Day0"
)

write_csv(
  baseline_rules,
  file.path(
    OUT,
    "V2_STEP92B2_baseline_mapping_dictionary.csv"
  )
)

# ------------------------------------------------------------
# Identify baseline
# ------------------------------------------------------------

dat2 <- dat %>%
  left_join(
    baseline_rules,
    by="project"
  ) %>%
  mutate(
    is_baseline =
      str_detect(
        time_factor,
        baseline_regex
      )
  )

baseline_check <- dat2 %>%
  group_by(
    project,
    patient_id
  ) %>%
  summarise(
    n_baseline =
      sum(is_baseline,na.rm=TRUE),
    .groups="drop"
  )

write_csv(
  baseline_check,
  file.path(
    OUT,
    "V2_STEP92B2_baseline_patient_check.csv"
  )
)

# keep patients with exactly one baseline

eligible <- baseline_check %>%
  filter(
    n_baseline==1
  ) %>%
  select(
    project,
    patient_id
  )

attrition <- baseline_check %>%
  mutate(
    retained =
      n_baseline==1
  )

write_csv(
  attrition,
  file.path(
    OUT,
    "V2_STEP92B2_baseline_attrition.csv"
  )
)

dat3 <- dat2 %>%
  inner_join(
    eligible,
    by=c(
      "project",
      "patient_id"
    )
  )

# ------------------------------------------------------------
# Baseline anchored instability
# ------------------------------------------------------------

eii <- dat3 %>%
  group_by(
    project,
    patient_id
  ) %>%
  mutate(

    baseline_shannon =
      Shannon[
        is_baseline
      ][1],

    baseline_simpson =
      Simpson[
        is_baseline
      ][1],

    baseline_asv =
      Observed_ASV[
        is_baseline
      ][1],

    shannon_instability =
      abs(
        Shannon -
          baseline_shannon
      ),

    simpson_instability =
      abs(
        Simpson -
          baseline_simpson
      ),

    asv_instability =
      abs(
        Observed_ASV -
          baseline_asv
      )

  ) %>%
  ungroup()


zscore <- function(x){
  if(sd(x,na.rm=TRUE)==0){
    return(rep(0,length(x)))
  }
  as.numeric(scale(x))
}


eii <- eii %>%
  mutate(

    z_bray =
      zscore(
        bray_from_patient_baseline
      ),

    z_shannon =
      zscore(
        shannon_instability
      ),

    z_simpson =
      zscore(
        simpson_instability
      ),

    z_asv =
      zscore(
        asv_instability
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


rng <- range(
  eii$EII_raw,
  na.rm=TRUE
)


eii <- eii %>%
  mutate(
    EII_0_100 =
      100*
      (
        EII_raw-rng[1]
      )/
      diff(rng)
  )


write_csv(
  eii,
  file.path(
    OUT,
    "V2_STEP92B2_corrected_timepoint_EII.csv"
  )
)


patient_summary <- eii %>%
  group_by(
    project,
    patient_id
  ) %>%
  summarise(
    n_timepoints=n(),
    max_EII=max(EII_0_100,na.rm=TRUE),
    final_EII=EII_0_100[
      which.max(EII_0_100)
    ][1],
    .groups="drop"
  )


write_csv(
  patient_summary,
  file.path(
    OUT,
    "V2_STEP92B2_patient_EII_summary.csv"
  )
)


writeLines(
  c(
    paste0(
      "Completed: ",
      Sys.time()
    ),
    "STEP92B2 COMPLETE",
    "Baseline harmonization completed."
  ),
  file.path(
    OUT,
    "_STEP92B2_COMPLETE.ok"
  )
)

cat("STEP92B2 COMPLETE\n")
