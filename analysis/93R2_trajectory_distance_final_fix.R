
# ============================================================
# Step93R2
# Final fix for trajectory distance calculation
#
# Previous issue:
# delta matrix still contained NA and metadata columns.
#
# This version:
# 1. selects only numeric genus delta columns
# 2. replaces NA/NaN/Inf with 0
# 3. calculates Bray-Curtis distance
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
  "V2_33R_TAXONOMIC_TRAJECTORY_DIVERGENCE_ANALYSIS",
  "V2_STEP93R_patient_genus_delta_matrix.csv"
)

OUT <- file.path(
  ROOT,
  "results",
  "V2_33R2_TRAJECTORY_DISTANCE_FINAL_FIX"
)

dir.create(OUT, recursive=TRUE, showWarnings=FALSE)

dat <- read_csv(
  IN,
  show_col_types=FALSE
)

# patient identifier
patient <- dat$Patient_true

# keep numeric genus delta features only
mat <- dat %>%
  select(where(is.numeric))

# remove non-genus numeric metadata if present
drop_cols <- c(
  "Time_numeric"
)

mat <- mat %>%
  select(
    -any_of(drop_cols)
  )

mat <- as.data.frame(
  lapply(
    mat,
    function(x){
      x[!is.finite(x)] <- 0
      x[is.na(x)] <- 0
      as.numeric(x)
    }
  )
)

rownames(mat) <- patient


write_csv(
  data.frame(
    Patient_true=rownames(mat),
    mat
  ),
  file.path(
    OUT,
    "patient_genus_delta_matrix_clean.csv"
  )
)

bc <- vegdist(
  abs(as.matrix(mat)),
  method="bray"
)

write.csv(
  as.matrix(bc),
  file.path(
    OUT,
    "trajectory_bray_distance_matrix.csv"
  )
)

summary <- data.frame(
  n_patients=nrow(mat),
  mean_pairwise_distance=mean(as.vector(bc))
)

write_csv(
  summary,
  file.path(
    OUT,
    "trajectory_heterogeneity_index.csv"
  )
)

writeLines(
  c(
    paste0("Completed: ", Sys.time()),
    "STEP93R2 COMPLETE"
  ),
  file.path(
    OUT,
    "_STEP93R2_COMPLETE.ok"
  )
)

cat("STEP93R2 COMPLETE\n")
