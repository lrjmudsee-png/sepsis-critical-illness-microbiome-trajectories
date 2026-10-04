
# ============================================================
# Step93R FIX
# Handle NA in genus delta matrix before Bray-Curtis
# ============================================================

options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(vegan)
})

ROOT <- "E:/sepsis_project"

IN <- file.path(
  ROOT,
  "results",
  "V2_33Q_TRAJECTORY_PATIENT_ID_RECONSTRUCTION_FIX",
  "V2_STEP93Q_fixed_longitudinal_genus_profile.csv"
)

OUT <- file.path(
  ROOT,
  "results",
  "V2_33R_TAXONOMIC_TRAJECTORY_DIVERGENCE_ANALYSIS"
)

dir.create(
  OUT,
  recursive=TRUE,
  showWarnings=FALSE
)

dat <- read_csv(
  IN,
  show_col_types=FALSE
)

meta_cols <- c(
  "Project",
  "Run",
  "BioSample",
  "Sample_Name",
  "Source",
  "Group",
  "Patient_ID",
  "Patient_ID_reconstructed",
  "Timepoint_reconstructed",
  "Patient_true"
)

genus_cols <- setdiff(
  names(dat),
  meta_cols
)


# repeated patients
patients <- dat %>%
  group_by(Patient_true) %>%
  summarise(
    n_timepoints=n_distinct(Time_numeric),
    .groups="drop"
  ) %>%
  filter(n_timepoints>=2)


dat2 <- dat %>%
  filter(Patient_true %in% patients$Patient_true)


baseline <- dat2 %>%
  group_by(Patient_true) %>%
  filter(Time_numeric == min(Time_numeric, na.rm=TRUE)) %>%
  slice(1) %>%
  ungroup()


followup <- dat2 %>%
  group_by(Patient_true) %>%
  filter(Time_numeric == max(Time_numeric, na.rm=TRUE)) %>%
  slice(1) %>%
  ungroup()


delta <- as.data.frame(
  t(
    mapply(
      function(i){
        as.numeric(followup[i,genus_cols]) -
          as.numeric(baseline[i,genus_cols])
      },
      seq_len(nrow(baseline))
    )
  )
)

# fix column structure
if(ncol(delta)!=length(genus_cols)){
  delta <- as.data.frame(
    do.call(
      rbind,
      lapply(
        seq_len(nrow(baseline)),
        function(i){
          as.numeric(followup[i,genus_cols]) -
            as.numeric(baseline[i,genus_cols])
        }
      )
    )
  )
}

colnames(delta) <- genus_cols

# critical fix:
# convert missing abundance to zero
delta <- delta %>%
  mutate(
    across(
      everything(),
      ~replace_na(as.numeric(.x),0)
    )
  )

delta$Patient_true <- baseline$Patient_true


write_csv(
  delta,
  file.path(
    OUT,
    "V2_STEP93R_patient_genus_delta_matrix_FIXED.csv"
  )
)


mat <- delta %>%
  select(all_of(genus_cols)) %>%
  as.matrix()

rownames(mat) <- delta$Patient_true


bc <- vegdist(
  abs(mat),
  method="bray"
)

write.csv(
  as.matrix(bc),
  file.path(
    OUT,
    "V2_STEP93R_taxonomic_trajectory_bray_distance_matrix_FIXED.csv"
  )
)


summary <- data.frame(
  mean_pairwise_distance=mean(as.vector(bc)),
  n_patients=nrow(mat)
)

write_csv(
  summary,
  file.path(
    OUT,
    "V2_STEP93R_taxonomic_trajectory_heterogeneity_index_FIXED.csv"
  )
)


writeLines(
  c(
    paste0("Completed: ",Sys.time()),
    "STEP93R FIX COMPLETE"
  ),
  file.path(
    OUT,
    "_STEP93R_FIX_COMPLETE.ok"
  )
)

cat("STEP93R FIX COMPLETE\n")
