# ============================================================
# Sepsis V2 Upgrade - Step 98E
# PRJNA1125274 (SURVEIL) external validation
# ------------------------------------------------------------
# Runs the pre-specified endpoints from PRJNA1125274_VALIDATION_SAP.md
# using the Step98D analysis object (counts/taxonomy/metadata).
#
# Primary: per-patient trajectory contrast
#   Delta = D(T0,T2) - D(T0,T1),  complete T0/T1/T2 patients only,
#   separately for Bray-Curtis and Aitchison (CZM primary;
#   pseudocount 0.5 / 1 sensitivity), never pooled together.
#   Discovery-consistent direction (from Step88B/98A) = INCREASE.
#
# Secondary: T0->T1 within-patient beta-diversity change via
#   restricted permutation (patient blocking), Bray & Aitchison.
#
# Sensitivity: excluding the 12 date-anomaly patients; per-hospital.
# Interpretation: SAP Situation A/B/C per metric.
# ============================================================

options(stringsAsFactors = FALSE)
suppressPackageStartupMessages({
  library(dada2); library(vegan)
  library(readr); library(tibble); library(purrr)
  library(stringr); library(zCompositions)
  # dplyr/tidyr last so dplyr::select masks MASS::select (vegan attaches MASS)
  library(dplyr); library(tidyr)
})

set.seed(20260907)

ROOT <- "E:/sepsis_project"
OUT <- file.path(ROOT, "results", "V2_UPGRADE_20260907",
                 "98E_PRJNA1125274_EXTERNAL_VALIDATION")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
DATA112 <- file.path(ROOT, "data", "PRJNA1125274")
OBJ <- file.path(DATA112, "03_dada2",
                 "PRJNA1125274_analysis_object_external_validation.rds")

write_csv_safe <- function(x, p) write_excel_csv(x, p, na = "")
J_paired <- function(n) 1 - 3 / (4 * (n - 1) - 1)

obj <- readRDS(OBJ)
counts <- obj$counts
meta <- obj$metadata
meta$run_id <- as.character(meta$run_accession)
meta$patient_id <- as.character(meta$patient_id)

cat("Runs in object:", nrow(counts), " ASVs:", ncol(counts), "\n")
stopifnot(all(rownames(counts) %in% meta$run_id))

# keep complete framework: run, patient, timepoint, hospital
df <- meta |>
  filter(timepoint %in% c("T0", "T1", "T2"), !is.na(patient_id)) |>
  mutate(depth = counts[run_id, , drop = FALSE] |> rowSums()) |>
  filter(depth >= 2000) |>
  select(run_id, patient_id, hospital, timepoint,
         alias_letter, collection_date, mapping_status,
         complete_T0_T1_T2, depth)

cat("Samples passing depth >= 2000:", nrow(df), "\n")

# ------------------------------------------------------------
# helper: Aitchison distance from personal baseline
# ------------------------------------------------------------
clr_matrix <- function(x) log(x) - rowMeans(log(x))
zero_impute <- function(counts, method = c("CZM", "PC0.5", "PC1")) {
  method <- match.arg(method)
  if (method == "CZM") {
    zCompositions::cmultRepl(counts, label = 0, method = "CZM",
                             output = "p-counts", z.warning = 1,
                             z.delete = FALSE, suppress.print = TRUE)
  } else if (method == "PC0.5") counts + 0.5 else counts + 1
}

displacement_table <- function(fw, runs) {
  # fw: framework data.frame (rows = runs), runs: count matrix subset rows
  out <- fw
  cnt <- counts[fw$run_id, , drop = FALSE]
  feat <- colSums(cnt)
  if (any(feat == 0)) cnt <- cnt[, feat > 0, drop = FALSE]
  # Technical amendment (frozen, see SAP addendum): prevalence prefilter for
  # tractable compositional distances. Geometry and zero-handling unchanged.
  keep_prevalent <- colSums(cnt) >= 20 & colSums(cnt > 0) >= 3
  cnt <- cnt[, keep_prevalent, drop = FALSE]
  cat("    prevalence-filtered ASVs:", ncol(cnt), "\n")
  rel <- cnt / rowSums(cnt)
  bray_mat <- as.matrix(vegan::vegdist(rel, method = "bray"))

  imp <- zero_impute(cnt, "CZM")
  dst <- as.matrix(dist(clr_matrix(imp), method = "euclidean"))
  imp05 <- zero_impute(cnt, "PC0.5")
  dst05 <- as.matrix(dist(clr_matrix(imp05), method = "euclidean"))
  imp1 <- zero_impute(cnt, "PC1")
  dst1 <- as.matrix(dist(clr_matrix(imp1), method = "euclidean"))

  # baseline per patient = T0 run
  base <- fw |> filter(timepoint == "T0") |>
    dplyr::select(patient_id, base_run = run_id) |> distinct(patient_id, .keep_all = TRUE)
  fw2 <- fw |> dplyr::left_join(base, by = "patient_id")
  # patients without a depth-passing T0 cannot anchor a T0-based trajectory
  n_no_base <- sum(is.na(fw2$base_run))
  if (n_no_base) cat("  dropping", n_no_base, "samples without depth-passing T0\n")
  fw2 <- fw2 |> dplyr::filter(!is.na(base_run))
  fw2 <- fw2 |>
    rowwise() |>
    mutate(
      bray_from_T0 = bray_mat[run_id, base_run],
      aitchison_from_T0_CZM = dst[run_id, base_run],
      aitchison_from_T0_PC0.5 = dst05[run_id, base_run],
      aitchison_from_T0_PC1 = dst1[run_id, base_run],
      is_baseline = run_id == base_run
    ) |>
    ungroup()
  # QC baseline distance == 0
  if (any(abs(fw2$bray_from_T0[fw2$is_baseline]) > 1e-9) ||
      any(abs(fw2$aitchison_from_T0_CZM[fw2$is_baseline]) > 1e-9)) {
    stop("baseline distance QC failed")
  }
  fw2
}

# ------------------------------------------------------------
# primary trajectory contrast
# ------------------------------------------------------------
primary_contrast <- function(disp, metric, exclude_date_anomaly = FALSE) {
  d <- disp |> filter(timepoint %in% c("T1", "T2"))
  if (exclude_date_anomaly) {
    d <- d |> filter(mapping_status != "LETTER_ORDER_TIME_CODE_DATE_ANOMALY")
  }
  pair <- d |>
    select(patient_id, timepoint, val = all_of(metric)) |>
    distinct(patient_id, timepoint, .keep_all = TRUE) |>
    pivot_wider(names_from = timepoint, values_from = val) |>
    filter(!is.na(T1), !is.na(T2))
  diffv <- pair$T2 - pair$T1
  n <- length(diffv)
  m <- mean(diffv); s <- sd(diffv)
  dz <- if (s > 0) m / s else NA_real_
  gz <- J_paired(n) * dz
  se_g <- J_paired(n) * sqrt(1 / n + dz^2 / (2 * n))
  tcrit <- qt(0.975, n - 1)
  tstat <- if (s > 0) m / (s / sqrt(n)) else NA_real_
  p_t <- 2 * pt(-abs(tstat), n - 1)
  p_w <- suppressWarnings(wilcox.test(diffv, mu = 0, exact = FALSE,
                                      correct = FALSE)$p.value)
  tibble(
    metric = metric, n_patients = n,
    n_excluded_anomaly = sum(disp$mapping_status ==
                               "LETTER_ORDER_TIME_CODE_DATE_ANOMALY"),
    mean_paired_diff = m, sd_paired_diff = s,
    median_paired_diff = median(diffv),
    cohen_dz = dz, hedges_gz = gz,
    gz_ci_low = gz - tcrit * se_g, gz_ci_high = gz + tcrit * se_g,
    paired_t_p = p_t, wilcoxon_p = p_w,
    direction = case_when(m > 0 ~ "INCREASE", m < 0 ~ "DECREASE",
                          TRUE ~ "NO_CHANGE"),
    discovery_concordant_direction = "INCREASE"
  )
}

# ------------------------------------------------------------
# secondary: T0->T1 restricted-permutation beta diversity
# ------------------------------------------------------------
t0t1_beta <- function(disp) {
  d <- disp |> filter(timepoint %in% c("T0", "T1"))
  # only patients with both T0 and T1
  keep <- d |> group_by(patient_id) |> filter(n_distinct(timepoint) == 2) |>
    ungroup()
  cnt <- counts[keep$run_id, , drop = FALSE]
  feat <- colSums(cnt)
  if (any(feat == 0)) cnt <- cnt[, feat > 0, drop = FALSE]
  keep_prevalent <- colSums(cnt) >= 20 & colSums(cnt > 0) >= 3
  cnt <- cnt[, keep_prevalent, drop = FALSE]
  rel <- cnt / rowSums(cnt)
  bc <- as.matrix(vegan::vegdist(rel, method = "bray"))
  cl <- clr_matrix(zero_impute(cnt, "CZM"))
  ac <- as.matrix(dist(cl, method = "euclidean"))
  time_factor <- factor(keep$timepoint)
  res <- list()
  for (geom in c("bray", "aitchison")) {
    dm <- if (geom == "bray") bc else ac
    set.seed(20260907)
    a <- vegan::adonis2(dm ~ time_factor, strata = keep$patient_id,
                        permutations = 999, by = "terms")
    res[[geom]] <- tibble(
      geometry = geom,
      n_samples = nrow(keep),
      n_patients = n_distinct(keep$patient_id),
      F = a["time_factor", "F"],
      R2 = a["time_factor", "R2"],
      p = a["time_factor", "Pr(>F)"]
    )
  }
  bind_rows(res)
}

cat("Building displacement table...\n")
LOG <- file.path(OUT, "_STEP98E_runtime.txt")
cat(paste0(Sys.time(), " start\n"), file = LOG)
stamp <- function(x) cat(paste0(Sys.time(), " ", x, "\n"), file = LOG, append = TRUE)
disp <- displacement_table(df)
write_csv_safe(disp, file.path(
  OUT, "PRJNA1125274_external_validation_sample_displacement.csv"))
stamp("displacement written")

cat("Primary contrasts...\n")
metrics <- c("bray_from_T0", "aitchison_from_T0_CZM",
             "aitchison_from_T0_PC0.5", "aitchison_from_T0_PC1")
prim_all <- map_dfr(metrics, ~ primary_contrast(disp, .x, FALSE))
prim_noanom <- map_dfr(metrics, ~ primary_contrast(disp, .x, TRUE))
prim_all$sensitivity_set <- "ALL_COMPLETE_T0T1T2"
prim_noanom$sensitivity_set <- "EXCLUDING_DATE_ANOMALY"
prim <- bind_rows(prim_all, prim_noanom)
write_csv_safe(prim, file.path(
  OUT, "PRJNA1125274_external_validation_primary_contrast.csv")); stamp("primary written")

cat("Hospital-stratified primary (bray + aitchison CZM)...\n")
prim_hosp <- map_dfr(unique(disp$hospital), function(h) {
  dd <- disp |> filter(hospital == h)
  map_dfr(c("bray_from_T0", "aitchison_from_T0_CZM"),
          ~ primary_contrast(dd, .x, FALSE)) |>
    mutate(hospital = h)
})
write_csv_safe(prim_hosp, file.path(
  OUT, "PRJNA1125274_external_validation_primary_by_hospital.csv")); stamp("hospital written")

cat("Secondary T0->T1 community change...\n")
beta <- t0t1_beta(disp)
write_csv_safe(beta, file.path(
  OUT, "PRJNA1125274_T0_T1_community_change_permutation.csv")); stamp("beta written")

cat("Patient/run accounting...\n")
accounting <- tibble(
  runs_total = nrow(df),
  runs_T0 = sum(df$timepoint == "T0"),
  runs_T1 = sum(df$timepoint == "T1"),
  runs_T2 = sum(df$timepoint == "T2"),
  patients_total = n_distinct(df$patient_id),
  patients_complete_T0T1T2 = n_distinct(
    df$patient_id[df$complete_T0_T1_T2]),
  patients_after_depth_T0T1T2 = prim_all$n_patients[1]
)
write_csv_safe(accounting, file.path(
  OUT, "PRJNA1125274_external_validation_accounting.csv")); stamp("accounting written")

# ------------------------------------------------------------
# conclusion per SAP rules
# ------------------------------------------------------------
conclude <- function(row) {
  if (is.na(row$paired_t_p)) return("NO_DATA")
  dir_match <- row$direction == row$discovery_concordant_direction
  sig <- row$paired_t_p < 0.05
  if (!dir_match) return("Situation C: discordant replication / transportability limitation")
  if (dir_match & sig) return("Situation A: independent external replication supported")
  return("Situation B: directionally concordant but imprecise replication")
}
conclusion <- prim |>
  filter(metric %in% c("bray_from_T0", "aitchison_from_T0_CZM"),
         sensitivity_set == "ALL_COMPLETE_T0T1T2") |>
  rowwise() |>
  mutate(sap_conclusion = conclude(cur_data())) |>
  ungroup() |>
  select(metric, sap_conclusion, hedges_gz, gz_ci_low, gz_ci_high,
         paired_t_p, direction)
write_csv_safe(conclusion, file.path(
  OUT, "PRJNA1125274_external_validation_conclusion.csv"))
cat("DONE 98E\n")
