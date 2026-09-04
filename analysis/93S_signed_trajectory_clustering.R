
# ============================================================
# Sepsis V2 - Step93S
# Signed taxonomic trajectory and patient clustering
#
# Input:
# patient_genus_delta_matrix_clean.csv
#
# Outputs:
# signed distance matrix
# PCA coordinates
# patient clusters
# ============================================================

options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(vegan)
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
  "V2_33S_SIGNED_TRAJECTORY_CLUSTERING"
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

patient <- dat$Patient_true

mat <- dat %>%
  select(-Patient_true) %>%
  as.matrix()

rownames(mat) <- patient

mat[is.na(mat)] <- 0
mat[!is.finite(mat)] <- 0


# signed Bray-Curtis
signed_distance <- vegdist(
  mat,
  method="bray"
)

write.csv(
  as.matrix(signed_distance),
  file.path(
    OUT,
    "signed_trajectory_bray_distance_matrix.csv"
  )
)


# PCA
pca <- prcomp(
  mat,
  scale.=TRUE
)

pca_df <- data.frame(
  Patient_true=rownames(pca$x),
  PC1=pca$x[,1],
  PC2=pca$x[,2]
)

write_csv(
  pca_df,
  file.path(
    OUT,
    "trajectory_PCA_coordinates.csv"
  )
)


# hierarchical clustering
hc <- hclust(
  signed_distance,
  method="ward.D2"
)

clusters <- data.frame(
  Patient_true=names(hc$order),
  cluster=cutree(
    hc,
    k=2
  )
)

write_csv(
  clusters,
  file.path(
    OUT,
    "trajectory_clusters.csv"
  )
)


writeLines(
  c(
    paste0("Completed: ",Sys.time()),
    "STEP93S COMPLETE"
  ),
  file.path(
    OUT,
    "_STEP93S_COMPLETE.ok"
  )
)

cat("STEP93S COMPLETE\n")
