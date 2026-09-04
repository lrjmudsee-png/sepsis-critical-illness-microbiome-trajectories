# ============================================================
# Sepsis V2 - Step93U
# TAXONOMIC TRAJECTORY HETEROGENEITY x ECOLOGICAL INSTABILITY
# ROBUSTNESS + FIGURE-READY OUTPUT
#
# Rationale
# ---------
# Step93T10B achieved 10/10 patient matching for PRJNA691455.
#
# Important methodological decision:
# The previous k=2 cluster solution is 9 vs 1, so it is NOT suitable
# for inferential "trajectory subtype" claims.
#
# Main analysis here therefore treats taxonomic trajectory as a
# continuous/directional phenomenon:
#
#   heterogeneous signed genus-change routes
#          +
#   within-patient ecological displacement
#
# Outputs:
# - patient-level ecological endpoints
# - pairwise trajectory direction heterogeneity
# - cluster validity audit
# - Spearman + permutation + leave-one-out robustness
# - Day3 -> Day7 paired Bray confirmation
# - provisional figure-ready PDF/PNG
#
# Interpretation remains exploratory because n=10.
# ============================================================

options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(stringr)
  library(tibble)
})

ROOT <- "E:/sepsis_project"

IN_T10B <- file.path(
  ROOT,
  "results",
  "V2_33T10B_FORCED_RUN_BRIDGE"
)

IN_S2 <- file.path(
  ROOT,
  "results",
  "V2_33S2_TRAJECTORY_DISTANCE_METHOD_FIX"
)

OUT <- file.path(
  ROOT,
  "results",
  "V2_33U_TRAJECTORY_ECOLOGY_ROBUSTNESS"
)

dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

set.seed(20260822)

# ------------------------------------------------------------
# helpers
# ------------------------------------------------------------

norm_patient <- function(x) {
  x <- toupper(str_trim(as.character(x)))
  x <- str_replace(x, regex("^PRJNA691455[_:\\-]*", ignore_case=TRUE), "")
  x <- str_replace(x, regex("FFO.*$", ignore_case=TRUE), "")
  x[x %in% c("", "NA", "N/A", "NULL")] <- NA_character_
  x
}

safe_spearman <- function(x, y) {
  ok <- complete.cases(x, y)
  x <- x[ok]
  y <- y[ok]

  if (length(x) < 3) {
    return(c(n=length(x), rho=NA, p=NA))
  }

  ct <- suppressWarnings(
    cor.test(x, y, method="spearman", exact=FALSE)
  )

  c(
    n=length(x),
    rho=unname(ct$estimate),
    p=ct$p.value
  )
}

perm_spearman <- function(x, y, B=10000) {
  ok <- complete.cases(x, y)
  x <- x[ok]
  y <- y[ok]

  if (length(x) < 4) return(NA_real_)

  obs <- suppressWarnings(
    cor(x, y, method="spearman")
  )

  perm <- replicate(
    B,
    suppressWarnings(
      cor(x, sample(y), method="spearman")
    )
  )

  (sum(abs(perm) >= abs(obs), na.rm=TRUE) + 1) / (sum(!is.na(perm)) + 1)
}

loo_spearman <- function(x, y) {
  ok <- complete.cases(x, y)
  x <- x[ok]
  y <- y[ok]

  if (length(x) < 5) {
    return(c(min=NA, median=NA, max=NA))
  }

  rhos <- sapply(seq_along(x), function(i) {
    suppressWarnings(
      cor(
        x[-i],
        y[-i],
        method="spearman"
      )
    )
  })

  c(
    min=min(rhos, na.rm=TRUE),
    median=median(rhos, na.rm=TRUE),
    max=max(rhos, na.rm=TRUE)
  )
}

# ------------------------------------------------------------
# 1. Load Step93T10B final integration
# ------------------------------------------------------------

integrated_file <- file.path(
  IN_T10B,
  "08_FINAL_patient_EII_trajectory_integrated.csv"
)

mapped_file <- file.path(
  IN_T10B,
  "04_EII_rows_mapped_to_Patient_true.csv"
)

if (!file.exists(integrated_file) || !file.exists(mapped_file)) {
  stop("Step93T10B outputs missing.")
}

integrated <- read_csv(
  integrated_file,
  show_col_types=FALSE
)

mapped <- read_csv(
  mapped_file,
  show_col_types=FALSE,
  guess_max=10000
)

integrated <- integrated %>%
  mutate(Patient_true = norm_patient(Patient_true))

mapped <- mapped %>%
  mutate(Patient_true = norm_patient(Patient_true))

# ------------------------------------------------------------
# 2. Latest ecological endpoint per patient
# ------------------------------------------------------------

if (!"analysis_time_order" %in% names(mapped)) {
  stop("analysis_time_order missing from Step93T10B mapped rows.")
}

latest <- mapped %>%
  arrange(Patient_true, analysis_time_order) %>%
  group_by(Patient_true) %>%
  slice_tail(n=1) %>%
  ungroup() %>%
  transmute(
    Patient_true,
    latest_day = analysis_time_order,
    latest_Bray = bray_from_patient_baseline,
    latest_Shannon_instability = shannon_instability,
    latest_Simpson_instability = simpson_instability,
    latest_ASV_instability = asv_instability,
    latest_EII = EII_0_100
  )

patient_final <- integrated %>%
  left_join(latest, by="Patient_true")

write_csv(
  patient_final,
  file.path(OUT, "01_patient_trajectory_ecology_endpoints.csv")
)

# ------------------------------------------------------------
# 3. Ecological endpoint summary
# ------------------------------------------------------------

eco_summary <- tibble(
  metric = c(
    "latest_Bray",
    "latest_EII",
    "latest_Shannon_instability",
    "latest_Simpson_instability",
    "latest_ASV_instability"
  )
) %>%
  rowwise() %>%
  mutate(
    n = sum(!is.na(patient_final[[metric]])),
    median = median(patient_final[[metric]], na.rm=TRUE),
    q1 = quantile(patient_final[[metric]], 0.25, na.rm=TRUE),
    q3 = quantile(patient_final[[metric]], 0.75, na.rm=TRUE),
    min = min(patient_final[[metric]], na.rm=TRUE),
    max = max(patient_final[[metric]], na.rm=TRUE)
  ) %>%
  ungroup()

write_csv(
  eco_summary,
  file.path(OUT, "02_ecological_endpoint_summary.csv")
)

# ------------------------------------------------------------
# 4. Day3 -> Day7 paired Bray confirmation
# ------------------------------------------------------------

bray_37 <- mapped %>%
  filter(
    analysis_time_order %in% c(3,7)
  ) %>%
  select(
    Patient_true,
    analysis_time_order,
    bray_from_patient_baseline
  ) %>%
  distinct() %>%
  tidyr::pivot_wider(
    names_from=analysis_time_order,
    values_from=bray_from_patient_baseline,
    names_prefix="Day"
  ) %>%
  filter(
    !is.na(Day3),
    !is.na(Day7)
  ) %>%
  mutate(
    change_Day7_minus_Day3 = Day7 - Day3
  )

write_csv(
  bray_37,
  file.path(OUT, "03_paired_Day3_Day7_Bray.csv")
)

if (nrow(bray_37) >= 3) {
  wt <- wilcox.test(
    bray_37$Day7,
    bray_37$Day3,
    paired=TRUE,
    exact=FALSE
  )

  bray_test <- tibble(
    n_pairs=nrow(bray_37),
    median_Day3=median(bray_37$Day3),
    median_Day7=median(bray_37$Day7),
    median_change=median(bray_37$change_Day7_minus_Day3),
    mean_change=mean(bray_37$change_Day7_minus_Day3),
    wilcoxon_p=wt$p.value
  )
} else {
  bray_test <- tibble(
    n_pairs=nrow(bray_37),
    median_Day3=NA,
    median_Day7=NA,
    median_change=NA,
    mean_change=NA,
    wilcoxon_p=NA
  )
}

write_csv(
  bray_test,
  file.path(OUT, "04_paired_Day3_Day7_Bray_test.csv")
)

# ------------------------------------------------------------
# 5. Read signed trajectory matrix
# ------------------------------------------------------------

signed_file <- file.path(
  IN_S2,
  "signed_trajectory_matrix_clean.csv"
)

if (!file.exists(signed_file)) {
  stop("signed_trajectory_matrix_clean.csv missing.")
}

signed <- read_csv(
  signed_file,
  show_col_types=FALSE
)

patient_candidates <- c(
  "Patient_true",
  "patient_true",
  "Patient",
  "patient_id"
)

pc <- patient_candidates[patient_candidates %in% names(signed)][1]

if (is.na(pc)) {
  pc <- names(signed)[1]
}

patients_signed <- norm_patient(signed[[pc]])

x <- signed %>%
  select(where(is.numeric))

X <- as.matrix(x)
rownames(X) <- patients_signed

# remove non-finite columns
keep <- apply(X, 2, function(z) all(is.finite(z)))
X <- X[, keep, drop=FALSE]

# ------------------------------------------------------------
# 6. Pairwise trajectory distances + cosine similarity
# ------------------------------------------------------------

pairs <- combn(seq_len(nrow(X)), 2)

pairwise <- lapply(seq_len(ncol(pairs)), function(j) {

  i <- pairs[1,j]
  k <- pairs[2,j]

  a <- X[i,]
  b <- X[k,]

  euclid <- sqrt(sum((a-b)^2))

  denom <- sqrt(sum(a^2)) * sqrt(sum(b^2))

  cosine <- ifelse(
    denom == 0,
    NA_real_,
    sum(a*b) / denom
  )

  tibble(
    Patient_A=rownames(X)[i],
    Patient_B=rownames(X)[k],
    signed_Euclidean=euclid,
    cosine_similarity=cosine,
    cosine_dissimilarity=1-cosine
  )
})

pairwise <- bind_rows(pairwise)

write_csv(
  pairwise,
  file.path(OUT, "05_pairwise_signed_trajectory_heterogeneity.csv")
)

heterogeneity <- tibble(
  n_patients=nrow(X),
  n_pairs=nrow(pairwise),
  n_features=ncol(X),
  mean_pairwise_Euclidean=mean(pairwise$signed_Euclidean, na.rm=TRUE),
  median_pairwise_Euclidean=median(pairwise$signed_Euclidean, na.rm=TRUE),
  mean_cosine_similarity=mean(pairwise$cosine_similarity, na.rm=TRUE),
  median_cosine_similarity=median(pairwise$cosine_similarity, na.rm=TRUE),
  proportion_cosine_le_0=mean(pairwise$cosine_similarity <= 0, na.rm=TRUE)
)

write_csv(
  heterogeneity,
  file.path(OUT, "06_trajectory_heterogeneity_summary.csv")
)

# ------------------------------------------------------------
# 7. Cluster validity audit
# ------------------------------------------------------------

cluster_sizes <- patient_final %>%
  count(cluster, name="n")

min_cluster_size <- min(cluster_sizes$n)

cluster_audit <- tibble(
  requested_k=2,
  n_patients=nrow(patient_final),
  cluster_sizes=paste(
    paste0(cluster_sizes$cluster, ":", cluster_sizes$n),
    collapse=";"
  ),
  min_cluster_size=min_cluster_size,
  inferential_cluster_comparison_allowed=min_cluster_size >= 2,
  recommendation=ifelse(
    min_cluster_size >= 2,
    "Exploratory cluster comparison possible",
    "Do not use subtype/cluster inferential claims; retain clusters only as exploratory visualization"
  )
)

write_csv(
  cluster_audit,
  file.path(OUT, "07_cluster_validity_audit.csv")
)

# silhouette audit for k=2..4
silhouette_out <- list()

if (requireNamespace("cluster", quietly=TRUE) && nrow(X) >= 5) {

  D <- dist(X, method="euclidean")
  hc <- hclust(D, method="ward.D2")

  for(k in 2:min(4, nrow(X)-1)) {
    cl <- cutree(hc, k=k)

    sil <- tryCatch(
      cluster::silhouette(cl, D),
      error=function(e) NULL
    )

    silhouette_out[[length(silhouette_out)+1]] <- tibble(
      k=k,
      min_cluster_size=min(table(cl)),
      max_cluster_size=max(table(cl)),
      mean_silhouette=ifelse(
        is.null(sil),
        NA_real_,
        mean(sil[, "sil_width"])
      )
    )
  }
}

silhouette_out <- bind_rows(silhouette_out)

write_csv(
  silhouette_out,
  file.path(OUT, "08_cluster_silhouette_k2_to_k4.csv")
)

# ------------------------------------------------------------
# 8. Continuous trajectory x ecology robustness
# ------------------------------------------------------------

predictors <- c(
  "mean_signed_trajectory_distance",
  "median_signed_trajectory_distance"
)

outcomes <- c(
  "latest_Bray",
  "latest_EII",
  "latest_Shannon_instability",
  "latest_Simpson_instability",
  "latest_ASV_instability"
)

robust_list <- list()

for(pred in predictors) {
  for(out in outcomes) {

    if (!(pred %in% names(patient_final)) ||
        !(out %in% names(patient_final))) next

    xx <- patient_final[[pred]]
    yy <- patient_final[[out]]

    s <- safe_spearman(xx, yy)
    loo <- loo_spearman(xx, yy)
    pp <- perm_spearman(xx, yy, B=10000)

    robust_list[[length(robust_list)+1]] <- tibble(
      predictor=pred,
      outcome=out,
      n=as.integer(s["n"]),
      spearman_rho=as.numeric(s["rho"]),
      asymptotic_p=as.numeric(s["p"]),
      permutation_p_10000=pp,
      LOO_rho_min=as.numeric(loo["min"]),
      LOO_rho_median=as.numeric(loo["median"]),
      LOO_rho_max=as.numeric(loo["max"])
    )
  }
}

robust <- bind_rows(robust_list) %>%
  mutate(
    interpretation=case_when(
      n < 8 ~ "insufficient_n",
      permutation_p_10000 < 0.05 &
        LOO_rho_min > 0 ~ "directionally_robust_positive",
      permutation_p_10000 < 0.05 &
        LOO_rho_max < 0 ~ "directionally_robust_negative",
      permutation_p_10000 < 0.05 ~ "significant_but_LOO_unstable",
      TRUE ~ "exploratory_not_significant"
    )
  )

write_csv(
  robust,
  file.path(OUT, "09_trajectory_ecology_robustness.csv")
)

# ------------------------------------------------------------
# 9. Provisional figure-ready graphics
# ------------------------------------------------------------

# PDF
pdf(
  file.path(OUT, "Figure_STEP93U_trajectory_ecology_robustness.pdf"),
  width=10,
  height=8
)

par(mfrow=c(2,2), mar=c(4.5,4.5,2.5,1))

# A: latest Bray
bp <- barplot(
  patient_final$latest_Bray,
  names.arg=patient_final$Patient_true,
  las=2,
  ylim=c(0,1),
  ylab="Latest Bray-Curtis from baseline",
  main="A. Within-patient ecological displacement"
)
abline(h=median(patient_final$latest_Bray, na.rm=TRUE), lty=2)

# B: trajectory distance vs latest EII
plot(
  patient_final$mean_signed_trajectory_distance,
  patient_final$latest_EII,
  pch=19,
  xlab="Mean signed-trajectory distance to other patients",
  ylab="Latest EII",
  main="B. Trajectory heterogeneity vs EII"
)
text(
  patient_final$mean_signed_trajectory_distance,
  patient_final$latest_EII,
  labels=patient_final$Patient_true,
  pos=3,
  cex=.65
)

# C: trajectory distance vs Bray
plot(
  patient_final$mean_signed_trajectory_distance,
  patient_final$latest_Bray,
  pch=19,
  xlab="Mean signed-trajectory distance to other patients",
  ylab="Latest Bray-Curtis from baseline",
  main="C. Trajectory heterogeneity vs Bray"
)
text(
  patient_final$mean_signed_trajectory_distance,
  patient_final$latest_Bray,
  labels=patient_final$Patient_true,
  pos=3,
  cex=.65
)

# D: Day3 vs Day7 paired
if (nrow(bray_37) > 0) {
  matplot(
    t(as.matrix(bray_37[,c("Day3","Day7")])),
    type="l",
    lty=1,
    xaxt="n",
    xlab="",
    ylab="Bray-Curtis from baseline",
    main="D. Paired Day3 to Day7 displacement"
  )
  axis(1, at=1:2, labels=c("Day3","Day7"))
} else {
  plot.new()
  title("D. Day3-Day7 data unavailable")
}

dev.off()

# PNG
png(
  file.path(OUT, "Figure_STEP93U_trajectory_ecology_robustness.png"),
  width=1800,
  height=1400,
  res=180
)

par(mfrow=c(2,2), mar=c(4.5,4.5,2.5,1))

barplot(
  patient_final$latest_Bray,
  names.arg=patient_final$Patient_true,
  las=2,
  ylim=c(0,1),
  ylab="Latest Bray-Curtis from baseline",
  main="A. Within-patient ecological displacement"
)
abline(h=median(patient_final$latest_Bray, na.rm=TRUE), lty=2)

plot(
  patient_final$mean_signed_trajectory_distance,
  patient_final$latest_EII,
  pch=19,
  xlab="Mean signed-trajectory distance to other patients",
  ylab="Latest EII",
  main="B. Trajectory heterogeneity vs EII"
)
text(
  patient_final$mean_signed_trajectory_distance,
  patient_final$latest_EII,
  labels=patient_final$Patient_true,
  pos=3,
  cex=.65
)

plot(
  patient_final$mean_signed_trajectory_distance,
  patient_final$latest_Bray,
  pch=19,
  xlab="Mean signed-trajectory distance to other patients",
  ylab="Latest Bray-Curtis from baseline",
  main="C. Trajectory heterogeneity vs Bray"
)
text(
  patient_final$mean_signed_trajectory_distance,
  patient_final$latest_Bray,
  labels=patient_final$Patient_true,
  pos=3,
  cex=.65
)

if (nrow(bray_37) > 0) {
  matplot(
    t(as.matrix(bray_37[,c("Day3","Day7")])),
    type="l",
    lty=1,
    xaxt="n",
    xlab="",
    ylab="Bray-Curtis from baseline",
    main="D. Paired Day3 to Day7 displacement"
  )
  axis(1, at=1:2, labels=c("Day3","Day7"))
} else {
  plot.new()
  title("D. Day3-Day7 data unavailable")
}

dev.off()

# ------------------------------------------------------------
# 10. Compact interpretation sheet
# ------------------------------------------------------------

key_rho <- robust %>%
  filter(
    predictor=="mean_signed_trajectory_distance",
    outcome %in% c("latest_Bray","latest_EII")
  )

interpretation <- c(
  "STEP93U INTERPRETATION",
  "",
  paste0(
    "Patients: ",
    nrow(patient_final)
  ),
  paste0(
    "Latest Bray median [IQR]: ",
    round(median(patient_final$latest_Bray, na.rm=TRUE),3),
    " [",
    round(quantile(patient_final$latest_Bray,.25,na.rm=TRUE),3),
    ", ",
    round(quantile(patient_final$latest_Bray,.75,na.rm=TRUE),3),
    "]"
  ),
  paste0(
    "k=2 cluster sizes: ",
    paste(cluster_sizes$n, collapse=" vs "),
    "."
  ),
  "",
  "Methodological conclusion:",
  "The 9-vs-1 cluster solution must NOT be treated as validated trajectory subtypes.",
  "Primary interpretation should use continuous signed-trajectory heterogeneity and patient-level ecological displacement.",
  "",
  "EII remains an exploratory ecological-state index, not a clinical prognostic score.",
  "Because EII contains Bray/alpha-instability components, EII associations are internal ecological consistency evidence rather than independent external validation."
)

writeLines(
  interpretation,
  file.path(OUT, "10_STEP93U_INTERPRETATION.txt")
)

writeLines(
  c(
    paste0("Completed: ", Sys.time()),
    "STEP93U COMPLETE"
  ),
  file.path(OUT, "_STEP93U_COMPLETE.ok")
)

cat("STEP93U COMPLETE\n")
