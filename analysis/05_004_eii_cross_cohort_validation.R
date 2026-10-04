# ============================================================
# Sepsis V2 - Step92C
# EII cross-cohort validation
#
# Purpose:
# Validate whether Ecological Instability Index is reproducible
# across longitudinal cohorts.
#
# Input:
# Step92B2 corrected EII
#
# Output:
# - cohort level summaries
# - trajectory summaries
# - statistical comparison tables
# ============================================================

options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
})

ROOT <- "E:/sepsis_project"

INPUT <- file.path(
  ROOT,
  "results",
  "V2_32B2_STEP92B_BASELINE_HARMONIZATION_FIXED",
  "V2_STEP92B2_corrected_timepoint_EII.csv"
)

OUT <- file.path(
  ROOT,
  "results",
  "V2_32C_STEP92C_EII_CROSS_COHORT_VALIDATION"
)

dir.create(
  OUT,
  recursive = TRUE,
  showWarnings = FALSE
)

if(!file.exists(INPUT)){
  stop("Step92B2 EII file missing.")
}

dat <- read_csv(
  INPUT,
  show_col_types = FALSE
)

# patient trajectory summary

patient <- dat %>%
  group_by(
    project,
    patient_id
  ) %>%
  summarise(
    n_timepoints=n(),
    baseline_EII =
      EII_0_100[is_baseline][1],
    max_EII =
      max(EII_0_100,na.rm=TRUE),
    final_EII =
      EII_0_100[
        which.max(
          as.numeric(
            factor(time_factor)
          )
        )
      ][1],
    EII_change =
      max_EII - baseline_EII,
    .groups="drop"
  )

write_csv(
  patient,
  file.path(
    OUT,
    "V2_STEP92C_patient_EII_trajectory_summary.csv"
  )
)

# cohort summary

cohort <- patient %>%
  group_by(project) %>%
  summarise(
    n_patients=n(),
    mean_EII_change =
      mean(
        EII_change,
        na.rm=TRUE
      ),
    median_EII_change =
      median(
        EII_change,
        na.rm=TRUE
      ),
    proportion_increase =
      mean(
        EII_change>0,
        na.rm=TRUE
      ),
    .groups="drop"
  )

write_csv(
  cohort,
  file.path(
    OUT,
    "V2_STEP92C_cohort_EII_validation_summary.csv"
  )
)

# time trajectory

trajectory <- dat %>%
  group_by(
    project,
    time_factor
  ) %>%
  summarise(
    n=n(),
    mean_EII =
      mean(
        EII_0_100,
        na.rm=TRUE
      ),
    median_EII =
      median(
        EII_0_100,
        na.rm=TRUE
      ),
    .groups="drop"
  )

write_csv(
  trajectory,
  file.path(
    OUT,
    "V2_STEP92C_EII_time_trajectory_summary.csv"
  )
)

# figure-ready plot

p <- ggplot(
  trajectory,
  aes(
    x=time_factor,
    y=mean_EII,
    group=project
  )
)+
  geom_line()+
  geom_point()+
  facet_wrap(~project, scales="free_x")+
  theme_bw()+
  labs(
    title="Ecological Instability Index trajectories",
    x="Time",
    y="Mean EII (0-100)"
  )

ggsave(
  file.path(
    OUT,
    "Figure_STEP92C_EII_trajectory.png"
  ),
  p,
  width=10,
  height=6,
  dpi=300
)

writeLines(
  c(
    paste0(
      "Completed: ",
      Sys.time()
    ),
    "STEP92C COMPLETE"
  ),
  file.path(
    OUT,
    "_STEP92C_COMPLETE.ok"
  )
)

cat("STEP92C COMPLETE\n")
