
# ============================================================
# Step93S2
# Method correction:
# Bray-Curtis cannot use signed negative delta vectors.
#
# Use:
# 1. Euclidean distance on signed delta for direction
# 2. CLR/PCA after removing zero variance columns
# 3. clustering on Euclidean trajectory distance
# ============================================================

options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(stats)
})

ROOT <- "E:/sepsis_project"

IN <- file.path(
  ROOT,
  "results",
  "V2_33R2_TRAJECTORY_DISTANCE_FINAL_FIX",
  "patient_genus_delta_matrix_clean.csv"
)

OUT <- file.path(
  ROOT,
  "results",
  "V2_33S2_TRAJECTORY_DISTANCE_METHOD_FIX"
)

dir.create(OUT, recursive=TRUE, showWarnings=FALSE)

dat <- read_csv(
  IN,
  show_col_types=FALSE
)

patient <- dat$Patient_true

mat <- dat %>%
  select(-Patient_true) %>%
  mutate(across(everything(), as.numeric)) %>%
  as.matrix()

rownames(mat) <- patient

mat[is.na(mat)] <- 0
mat[!is.finite(mat)] <- 0


# remove zero variance genus features
keep <- apply(
  mat,
  2,
  function(x) sd(x)!=0
)

mat <- mat[,keep,drop=FALSE]


# save cleaned signed trajectory matrix
write_csv(
  data.frame(
    Patient_true=rownames(mat),
    mat
  ),
  file.path(
    OUT,
    "signed_trajectory_matrix_clean.csv"
  )
)


# signed trajectory distance
dist_mat <- dist(
  mat,
  method="euclidean"
)

write.csv(
  as.matrix(dist_mat),
  file.path(
    OUT,
    "signed_trajectory_euclidean_distance_matrix.csv"
  )
)


# PCA
pca <- prcomp(
  mat,
  scale.=TRUE
)

write_csv(
  data.frame(
    Patient_true=rownames(pca$x),
    PC1=pca$x[,1],
    PC2=pca$x[,2]
  ),
  file.path(
    OUT,
    "trajectory_PCA_coordinates.csv"
  )
)


# clustering
hc <- hclust(
  dist_mat,
  method="ward.D2"
)

cluster_df <- data.frame(
  Patient_true=names(cutree(hc,k=2)),
  cluster=cutree(hc,k=2)
)

write_csv(
  cluster_df,
  file.path(
    OUT,
    "trajectory_clusters.csv"
  )
)


writeLines(
  c(
    paste0("Completed: ",Sys.time()),
    "STEP93S2 COMPLETE"
  ),
  file.path(
    OUT,
    "_STEP93S2_COMPLETE.ok"
  )
)

cat("STEP93S2 COMPLETE\n")
