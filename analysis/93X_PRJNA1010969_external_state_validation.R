# ============================================================
# Sepsis V2 - Step93X
# PRJNA1010969 CROSS-SECTIONAL EXTERNAL ECOLOGICAL STATE VALIDATION
#
# Scientific question
# -------------------
# After freezing the longitudinal evidence, test whether an independent
# cross-sectional cohort containing Control / Trauma / Sepsis samples
# shows:
#
# 1) group separation in genus-level community structure;
# 2) greater ecological distance from the healthy-control centroid in
#    Sepsis and/or Trauma;
# 3) alpha-diversity differences;
# 4) exploratory genus-level differences;
# 5) directional alignment between PRJNA691455 longitudinal genus change
#    and PRJNA1010969 Sepsis-vs-Control / Sepsis-vs-Trauma differences.
#
# IMPORTANT:
# - This is EXTERNAL CROSS-SECTIONAL STATE VALIDATION.
# - It is NOT a longitudinal replication.
# - "distance from healthy centroid" is not equivalent to within-patient
#   displacement from personal baseline.
# - Genus screens are exploratory and compositional.
# ============================================================

options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(stringr)
  library(tibble)
})

set.seed(20260822)

ROOT <- "E:/sepsis_project"
PROJECT <- "PRJNA1010969"

ABUND_FILE <- file.path(
  ROOT,
  "data",
  PROJECT,
  "03_dada2_rerun_v2",
  "PRJNA1010969_genus_relative_abundance.csv"
)

OUT <- file.path(
  ROOT,
  "results",
  "V2_33X_PRJNA1010969_EXTERNAL_STATE_VALIDATION"
)

dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

logfile <- file.path(OUT, "_STEP93X_runtime.txt")
writeLines(paste("START", Sys.time()), logfile)

logmsg <- function(...) {
  z <- paste0(...)
  cat(z, "\n")
  cat(z, "\n", file = logfile, append = TRUE)
}

norm_id <- function(x) {
  x <- toupper(str_trim(as.character(x)))
  x[x %in% c("", "NA", "N/A", "NULL")] <- NA_character_
  x
}

pick_exact <- function(nms, candidates) {
  x <- candidates[candidates %in% nms]
  if (length(x) == 0) return(NA_character_)
  x[1]
}

bray_vec <- function(a, b) {
  den <- sum(a + b, na.rm = TRUE)
  if (!is.finite(den) || den <= 0) return(NA_real_)
  sum(abs(a - b), na.rm = TRUE) / den
}

bray_matrix <- function(X) {
  n <- nrow(X)
  D <- matrix(0, n, n)
  rownames(D) <- rownames(X)
  colnames(D) <- rownames(X)

  for (i in seq_len(n - 1)) {
    for (j in (i + 1):n) {
      d <- bray_vec(X[i, ], X[j, ])
      D[i, j] <- d
      D[j, i] <- d
    }
  }
  D
}

pseudo_f_distance <- function(D, group) {
  group <- as.factor(group)
  n <- length(group)
  k <- nlevels(group)

  if (n <= k || k < 2) return(NA_real_)

  total_ss <- sum(D[upper.tri(D)]^2, na.rm = TRUE) / n

  within_ss <- 0
  for (g in levels(group)) {
    idx <- which(group == g)
    ng <- length(idx)
    if (ng >= 2) {
      dg <- D[idx, idx, drop = FALSE]
      within_ss <- within_ss +
        sum(dg[upper.tri(dg)]^2, na.rm = TRUE) / ng
    }
  }

  between_ss <- total_ss - within_ss

  (between_ss / (k - 1)) /
    (within_ss / (n - k))
}

permanova_custom <- function(D, group, B = 9999) {
  obs <- pseudo_f_distance(D, group)
  if (!is.finite(obs)) {
    return(tibble(
      n = length(group),
      groups = length(unique(group)),
      pseudo_F = NA_real_,
      permutation_p = NA_real_,
      permutations = B
    ))
  }

  perm <- replicate(
    B,
    pseudo_f_distance(D, sample(group))
  )

  p <- (sum(perm >= obs, na.rm = TRUE) + 1) /
    (sum(is.finite(perm)) + 1)

  tibble(
    n = length(group),
    groups = length(unique(group)),
    pseudo_F = obs,
    permutation_p = p,
    permutations = B
  )
}

cliff_delta <- function(a, b) {
  a <- a[is.finite(a)]
  b <- b[is.finite(b)]
  if (length(a) == 0 || length(b) == 0) return(NA_real_)
  cmp <- outer(a, b, FUN = "-")
  (sum(cmp > 0) - sum(cmp < 0)) / length(cmp)
}

safe_wilcox_p <- function(a, b) {
  a <- a[is.finite(a)]
  b <- b[is.finite(b)]
  if (length(a) < 2 || length(b) < 2) return(NA_real_)
  suppressWarnings(
    tryCatch(
      wilcox.test(a, b, exact = FALSE)$p.value,
      error = function(e) NA_real_
    )
  )
}

# ------------------------------------------------------------
# 1. Read abundance matrix
# ------------------------------------------------------------

if (!file.exists(ABUND_FILE)) {
  stop(paste0("Missing abundance file: ", ABUND_FILE))
}

abund <- read_csv(
  ABUND_FILE,
  show_col_types = FALSE,
  guess_max = 10000,
  name_repair = "unique"
)

run_col_abund <- pick_exact(
  names(abund),
  c("Run_ID", "Run", "run_id", "run", "RunID")
)

if (is.na(run_col_abund)) {
  stop("Cannot identify Run column in genus relative abundance table.")
}

run_ids <- norm_id(abund[[run_col_abund]])

# Numeric genus columns only; abundance table is expected to contain
# Run_ID + genus relative abundance features.
num_cols <- names(abund)[vapply(abund, is.numeric, logical(1))]

if (length(num_cols) < 10) {
  stop("Too few numeric genus columns detected.")
}

X <- as.matrix(abund[, num_cols, drop = FALSE])
storage.mode(X) <- "numeric"
X[!is.finite(X)] <- 0
X[X < 0] <- 0

rs <- rowSums(X)
keep_nonzero <- rs > 0

abund <- abund[keep_nonzero, , drop = FALSE]
X <- X[keep_nonzero, , drop = FALSE]
run_ids <- run_ids[keep_nonzero]

X <- X / rowSums(X)
rownames(X) <- run_ids

logmsg("Abundance samples: ", nrow(X))
logmsg("Genus features: ", ncol(X))

# ------------------------------------------------------------
# 2. Find the manually downloaded SraRunTable by exact Run overlap
# ------------------------------------------------------------

meta_candidates <- list.files(
  file.path(ROOT, "data", PROJECT),
  recursive = TRUE,
  full.names = TRUE,
  pattern = "SraRunTable.*\\.csv$",
  ignore.case = TRUE
)

if (length(meta_candidates) == 0) {
  stop("No PRJNA1010969 SraRunTable CSV found.")
}

meta_audit <- list()
meta_loaded <- list()

for (f in meta_candidates) {
  m <- tryCatch(
    read_csv(
      f,
      show_col_types = FALSE,
      guess_max = 10000,
      name_repair = "unique"
    ),
    error = function(e) NULL
  )

  if (is.null(m)) next

  rc <- pick_exact(
    names(m),
    c("Run", "Run_ID", "run_id", "run", "Run accession")
  )

  if (is.na(rc)) {
    hits <- names(m)[
      str_detect(names(m), regex("^run$|run.*access", ignore_case = TRUE))
    ]
    if (length(hits) > 0) rc <- hits[1]
  }

  if (is.na(rc)) next

  overlap <- length(
    intersect(
      unique(norm_id(m[[rc]])),
      unique(run_ids)
    )
  )

  meta_audit[[length(meta_audit) + 1]] <- tibble(
    file = f,
    rows = nrow(m),
    run_column = rc,
    exact_run_overlap = overlap,
    columns = paste(names(m), collapse = ";")
  )

  meta_loaded[[f]] <- m
}

meta_audit <- bind_rows(meta_audit) %>%
  arrange(desc(exact_run_overlap))

write_csv(
  meta_audit,
  file.path(OUT, "01_metadata_candidate_audit.csv")
)

if (nrow(meta_audit) == 0 || max(meta_audit$exact_run_overlap) == 0) {
  stop("No SraRunTable has exact Run overlap with PRJNA1010969 abundance.")
}

best_meta_file <- meta_audit$file[1]
best_run_col <- meta_audit$run_column[1]
meta <- meta_loaded[[best_meta_file]]

sample_name_col <- pick_exact(
  names(meta),
  c("Sample Name", "Sample_Name", "sample_name", "SampleName")
)

if (is.na(sample_name_col)) {
  hits <- names(meta)[
    str_detect(names(meta), regex("sample.*name", ignore_case = TRUE))
  ]
  if (length(hits) > 0) sample_name_col <- hits[1]
}

if (is.na(sample_name_col)) {
  stop("Cannot identify Sample_Name in PRJNA1010969 metadata.")
}

meta2 <- meta %>%
  transmute(
    Run_ID = norm_id(.data[[best_run_col]]),
    Sample_Name = as.character(.data[[sample_name_col]])
  ) %>%
  mutate(
    prefix = toupper(str_extract(Sample_Name, "^[A-Za-z]+")),
    Group = case_when(
      str_starts(toupper(Sample_Name), "HMB") ~ "Control",
      str_starts(toupper(Sample_Name), "TMB") ~ "Trauma",
      str_starts(toupper(Sample_Name), "SMB") ~ "Sepsis",
      TRUE ~ NA_character_
    )
  ) %>%
  distinct(Run_ID, .keep_all = TRUE)

mapping <- tibble(Run_ID = rownames(X)) %>%
  left_join(meta2, by = "Run_ID")

write_csv(
  mapping,
  file.path(OUT, "02_PRJNA1010969_Run_Group_mapping.csv")
)

group_counts <- mapping %>%
  count(Group, name = "n")

write_csv(
  group_counts,
  file.path(OUT, "03_group_sample_counts.csv")
)

if (sum(!is.na(mapping$Group)) < 45) {
  warning("Fewer than 45 abundance samples received reconstructed group labels.")
}

# Keep only labeled samples
idx <- which(!is.na(mapping$Group))
mapping2 <- mapping[idx, , drop = FALSE]
X2 <- X[mapping2$Run_ID, , drop = FALSE]

group <- factor(
  mapping2$Group,
  levels = c("Control", "Trauma", "Sepsis")
)

# ------------------------------------------------------------
# 3. Alpha diversity at genus level
# ------------------------------------------------------------

Observed_Genera <- rowSums(X2 > 0)

Shannon <- apply(X2, 1, function(p) {
  p <- p[p > 0]
  -sum(p * log(p))
})

Simpson <- apply(X2, 1, function(p) {
  1 - sum(p^2)
})

alpha <- mapping2 %>%
  mutate(
    Observed_Genera = Observed_Genera,
    Shannon = Shannon,
    Simpson = Simpson
  )

write_csv(
  alpha,
  file.path(OUT, "04_genus_alpha_diversity_by_sample.csv")
)

alpha_kw <- bind_rows(
  lapply(
    c("Observed_Genera", "Shannon", "Simpson"),
    function(metric) {
      kt <- kruskal.test(alpha[[metric]] ~ alpha$Group)
      tibble(
        metric = metric,
        n = sum(!is.na(alpha[[metric]])),
        kruskal_chisq = unname(kt$statistic),
        df = unname(kt$parameter),
        p_value = kt$p.value
      )
    }
  )
) %>%
  mutate(FDR = p.adjust(p_value, method = "BH"))

write_csv(
  alpha_kw,
  file.path(OUT, "05_alpha_global_tests.csv")
)

contrasts <- list(
  c("Sepsis", "Control"),
  c("Sepsis", "Trauma"),
  c("Trauma", "Control")
)

alpha_pair <- bind_rows(
  lapply(
    c("Observed_Genera", "Shannon", "Simpson"),
    function(metric) {
      bind_rows(
        lapply(contrasts, function(cc) {
          a <- alpha[[metric]][alpha$Group == cc[1]]
          b <- alpha[[metric]][alpha$Group == cc[2]]

          tibble(
            metric = metric,
            group_A = cc[1],
            group_B = cc[2],
            n_A = sum(is.finite(a)),
            n_B = sum(is.finite(b)),
            median_A = median(a, na.rm = TRUE),
            median_B = median(b, na.rm = TRUE),
            median_difference_A_minus_B =
              median(a, na.rm = TRUE) - median(b, na.rm = TRUE),
            cliff_delta_A_vs_B = cliff_delta(a, b),
            wilcoxon_p = safe_wilcox_p(a, b)
          )
        })
      )
    }
  )
) %>%
  mutate(FDR = p.adjust(wilcoxon_p, method = "BH"))

write_csv(
  alpha_pair,
  file.path(OUT, "06_alpha_pairwise_tests.csv")
)

# ------------------------------------------------------------
# 4. Bray-Curtis + PCoA + PERMANOVA
# ------------------------------------------------------------

D <- bray_matrix(X2)

write_csv(
  as.data.frame(D) %>%
    rownames_to_column("Run_ID"),
  file.path(OUT, "07_Bray_distance_matrix.csv")
)

perm_global <- permanova_custom(
  D,
  group,
  B = 9999
)

write_csv(
  perm_global,
  file.path(OUT, "08_Bray_PERMANOVA_global.csv")
)

perm_pair <- bind_rows(
  lapply(contrasts, function(cc) {
    ii <- which(group %in% cc)
    res <- permanova_custom(
      D[ii, ii, drop = FALSE],
      droplevels(group[ii]),
      B = 9999
    )
    res %>%
      mutate(
        group_A = cc[1],
        group_B = cc[2]
      )
  })
) %>%
  mutate(FDR = p.adjust(permutation_p, method = "BH"))

write_csv(
  perm_pair,
  file.path(OUT, "09_Bray_PERMANOVA_pairwise.csv")
)

# PCoA
pc <- cmdscale(
  as.dist(D),
  k = 2,
  eig = TRUE,
  add = TRUE
)

coords <- as.data.frame(pc$points)
names(coords)[1:2] <- c("PCoA1", "PCoA2")

pcoa <- mapping2 %>%
  bind_cols(coords)

eig <- pc$eig
pos_eig <- eig[eig > 0]

variance1 <- if (length(pos_eig) >= 1) pos_eig[1] / sum(pos_eig) else NA_real_
variance2 <- if (length(pos_eig) >= 2) pos_eig[2] / sum(pos_eig) else NA_real_

write_csv(
  pcoa,
  file.path(OUT, "10_Bray_PCoA_coordinates.csv")
)

write_csv(
  tibble(
    axis = c("PCoA1", "PCoA2"),
    variance_fraction_positive_eigenvalues = c(variance1, variance2)
  ),
  file.path(OUT, "11_Bray_PCoA_variance.csv")
)

# ------------------------------------------------------------
# 5. Distance to healthy-control centroid
# ------------------------------------------------------------

control_centroid <- colMeans(
  X2[group == "Control", , drop = FALSE]
)
control_centroid <- control_centroid / sum(control_centroid)

distance_control_centroid <- apply(
  X2,
  1,
  function(x) bray_vec(x, control_centroid)
)

centroid_df <- mapping2 %>%
  mutate(
    Bray_to_Control_Centroid = distance_control_centroid
  )

write_csv(
  centroid_df,
  file.path(OUT, "12_distance_to_Control_centroid.csv")
)

centroid_kw <- kruskal.test(
  centroid_df$Bray_to_Control_Centroid ~ centroid_df$Group
)

write_csv(
  tibble(
    n = nrow(centroid_df),
    kruskal_chisq = unname(centroid_kw$statistic),
    df = unname(centroid_kw$parameter),
    p_value = centroid_kw$p.value
  ),
  file.path(OUT, "13_Control_centroid_distance_global_test.csv")
)

centroid_pair <- bind_rows(
  lapply(contrasts, function(cc) {
    a <- centroid_df$Bray_to_Control_Centroid[centroid_df$Group == cc[1]]
    b <- centroid_df$Bray_to_Control_Centroid[centroid_df$Group == cc[2]]

    tibble(
      group_A = cc[1],
      group_B = cc[2],
      n_A = sum(is.finite(a)),
      n_B = sum(is.finite(b)),
      median_A = median(a, na.rm = TRUE),
      median_B = median(b, na.rm = TRUE),
      median_difference_A_minus_B =
        median(a, na.rm = TRUE) - median(b, na.rm = TRUE),
      cliff_delta_A_vs_B = cliff_delta(a, b),
      wilcoxon_p = safe_wilcox_p(a, b)
    )
  })
) %>%
  mutate(FDR = p.adjust(wilcoxon_p, method = "BH"))

write_csv(
  centroid_pair,
  file.path(OUT, "14_Control_centroid_distance_pairwise_tests.csv")
)

# ------------------------------------------------------------
# 6. Optional dispersion test using vegan if available
# ------------------------------------------------------------

dispersion_result <- tibble(
  method = "vegan_betadisper_if_available",
  available = FALSE,
  F = NA_real_,
  permutation_p = NA_real_,
  permutations = 9999
)

if (requireNamespace("vegan", quietly = TRUE)) {
  bd <- vegan::betadisper(as.dist(D), group)
  pt <- vegan::permutest(bd, permutations = 9999)

  dispersion_result <- tibble(
    method = "vegan::betadisper + permutest",
    available = TRUE,
    F = unname(pt$tab[1, "F"]),
    permutation_p = unname(pt$tab[1, "Pr(>F)"]),
    permutations = 9999
  )
}

write_csv(
  dispersion_result,
  file.path(OUT, "15_beta_dispersion_test.csv")
)

# ------------------------------------------------------------
# 7. Exploratory genus-level group screen
# ------------------------------------------------------------

prevalence <- colMeans(X2 > 0)

keep_genus <- names(prevalence)[
  prevalence >= 0.10
]

genus_list <- lapply(
  keep_genus,
  function(g) {
    y <- X2[, g]

    kw <- tryCatch(
      kruskal.test(y ~ group),
      error = function(e) NULL
    )

    sc <- y[group == "Sepsis"]
    co <- y[group == "Control"]
    tr <- y[group == "Trauma"]

    tibble(
      genus = g,
      prevalence = prevalence[g],
      median_Control = median(co, na.rm = TRUE),
      median_Trauma = median(tr, na.rm = TRUE),
      median_Sepsis = median(sc, na.rm = TRUE),
      KW_p = if (is.null(kw)) NA_real_ else kw$p.value,

      Sepsis_minus_Control =
        median(sc, na.rm = TRUE) - median(co, na.rm = TRUE),
      Sepsis_vs_Control_cliff = cliff_delta(sc, co),
      Sepsis_vs_Control_p = safe_wilcox_p(sc, co),

      Sepsis_minus_Trauma =
        median(sc, na.rm = TRUE) - median(tr, na.rm = TRUE),
      Sepsis_vs_Trauma_cliff = cliff_delta(sc, tr),
      Sepsis_vs_Trauma_p = safe_wilcox_p(sc, tr),

      Trauma_minus_Control =
        median(tr, na.rm = TRUE) - median(co, na.rm = TRUE),
      Trauma_vs_Control_cliff = cliff_delta(tr, co),
      Trauma_vs_Control_p = safe_wilcox_p(tr, co)
    )
  }
)

genus_screen <- bind_rows(genus_list) %>%
  mutate(
    KW_FDR = p.adjust(KW_p, method = "BH"),
    Sepsis_vs_Control_FDR =
      p.adjust(Sepsis_vs_Control_p, method = "BH"),
    Sepsis_vs_Trauma_FDR =
      p.adjust(Sepsis_vs_Trauma_p, method = "BH"),
    Trauma_vs_Control_FDR =
      p.adjust(Trauma_vs_Control_p, method = "BH")
  ) %>%
  arrange(KW_FDR, KW_p)

write_csv(
  genus_screen,
  file.path(OUT, "16_exploratory_genus_group_screen.csv")
)

# ------------------------------------------------------------
# 8. External alignment with PRJNA691455 longitudinal genus change
# ------------------------------------------------------------

signed_file <- file.path(
  ROOT,
  "results",
  "V2_33S2_TRAJECTORY_DISTANCE_METHOD_FIX",
  "signed_trajectory_matrix_clean.csv"
)

alignment <- NULL
alignment_summary <- NULL

if (file.exists(signed_file)) {

  signed <- read_csv(
    signed_file,
    show_col_types = FALSE,
    name_repair = "unique"
  )

  patient_cols <- c(
    "Patient_true",
    "patient_true",
    "Patient",
    "patient_id"
  )
  pcoll <- patient_cols[patient_cols %in% names(signed)][1]

  if (is.na(pcoll)) pcoll <- names(signed)[1]

  numeric_signed <- signed %>%
    select(where(is.numeric))

  longitudinal_mean_delta <- colMeans(
    as.matrix(numeric_signed),
    na.rm = TRUE
  )

  longitudinal <- tibble(
    genus = names(longitudinal_mean_delta),
    PRJNA691455_mean_signed_Day1_to_Day7_delta =
      as.numeric(longitudinal_mean_delta)
  )

  alignment <- genus_screen %>%
    inner_join(longitudinal, by = "genus")

  write_csv(
    alignment,
    file.path(OUT, "17_longitudinal_external_genus_alignment.csv")
  )

  if (nrow(alignment) >= 5) {

    c_sc <- suppressWarnings(
      cor.test(
        alignment$PRJNA691455_mean_signed_Day1_to_Day7_delta,
        alignment$Sepsis_minus_Control,
        method = "spearman",
        exact = FALSE
      )
    )

    c_st <- suppressWarnings(
      cor.test(
        alignment$PRJNA691455_mean_signed_Day1_to_Day7_delta,
        alignment$Sepsis_minus_Trauma,
        method = "spearman",
        exact = FALSE
      )
    )

    sign_sc <- mean(
      sign(alignment$PRJNA691455_mean_signed_Day1_to_Day7_delta) ==
        sign(alignment$Sepsis_minus_Control),
      na.rm = TRUE
    )

    sign_st <- mean(
      sign(alignment$PRJNA691455_mean_signed_Day1_to_Day7_delta) ==
        sign(alignment$Sepsis_minus_Trauma),
      na.rm = TRUE
    )

    alignment_summary <- tibble(
      comparison = c(
        "Longitudinal delta vs external Sepsis-Control",
        "Longitudinal delta vs external Sepsis-Trauma"
      ),
      overlapping_genera = nrow(alignment),
      spearman_rho = c(
        unname(c_sc$estimate),
        unname(c_st$estimate)
      ),
      p_value = c(
        c_sc$p.value,
        c_st$p.value
      ),
      sign_concordance = c(
        sign_sc,
        sign_st
      )
    )

    write_csv(
      alignment_summary,
      file.path(OUT, "18_longitudinal_external_alignment_summary.csv")
    )
  }
}

# ------------------------------------------------------------
# 9. Figure
# ------------------------------------------------------------

make_figure <- function() {

  par(mfrow = c(2, 2), mar = c(4.5, 4.5, 2.7, 1))

  # A PCoA
  pch_map <- c(Control = 1, Trauma = 17, Sepsis = 19)
  plot(
    pcoa$PCoA1,
    pcoa$PCoA2,
    pch = pch_map[pcoa$Group],
    xlab = paste0(
      "PCoA1 (",
      ifelse(is.na(variance1), "NA", paste0(round(100*variance1,1), "%")),
      ")"
    ),
    ylab = paste0(
      "PCoA2 (",
      ifelse(is.na(variance2), "NA", paste0(round(100*variance2,1), "%")),
      ")"
    ),
    main = "A. PRJNA1010969 genus-level Bray PCoA"
  )
  legend(
    "topright",
    legend = names(pch_map),
    pch = pch_map,
    bty = "n"
  )

  # B distance to healthy centroid
  boxplot(
    Bray_to_Control_Centroid ~ Group,
    data = centroid_df,
    ylab = "Bray-Curtis to healthy-control centroid",
    xlab = "",
    main = "B. Distance from healthy ecological state"
  )

  # C Shannon
  boxplot(
    Shannon ~ Group,
    data = alpha,
    ylab = "Genus-level Shannon",
    xlab = "",
    main = "C. Alpha diversity"
  )

  # D longitudinal/external alignment
  if (!is.null(alignment) && nrow(alignment) >= 5) {
    plot(
      alignment$PRJNA691455_mean_signed_Day1_to_Day7_delta,
      alignment$Sepsis_minus_Control,
      pch = 19,
      xlab = "PRJNA691455 mean signed longitudinal genus delta",
      ylab = "PRJNA1010969 Sepsis - Control median abundance",
      main = "D. Longitudinal / external genus alignment"
    )
    abline(h = 0, v = 0, lty = 2)
  } else {
    plot.new()
    title("D. Insufficient overlapping genera")
  }
}

pdf(
  file.path(OUT, "Figure_STEP93X_external_state_validation.pdf"),
  width = 10,
  height = 8
)
make_figure()
dev.off()

png(
  file.path(OUT, "Figure_STEP93X_external_state_validation.png"),
  width = 1800,
  height = 1400,
  res = 180
)
make_figure()
dev.off()

# ------------------------------------------------------------
# 10. Claim guardrails + interpretation
# ------------------------------------------------------------

centroid_sepsis_control <- centroid_pair %>%
  filter(
    group_A == "Sepsis",
    group_B == "Control"
  )

centroid_sepsis_trauma <- centroid_pair %>%
  filter(
    group_A == "Sepsis",
    group_B == "Trauma"
  )

guardrails <- tibble(
  claim = c(
    "PRJNA1010969 is an independent longitudinal replication",
    "PRJNA1010969 can externally validate a cross-sectional ecological state contrast",
    "Distance to control centroid is equivalent to within-patient displacement from baseline",
    "Sepsis-specific composition is supported if Sepsis differs from both Control and Trauma",
    "Genus-level Wilcoxon screen establishes causal biomarkers",
    "Longitudinal-external genus alignment is mechanistic proof"
  ),
  status = c(
    "NO",
    "YES",
    "NO",
    "CONDITIONAL_ON_RESULTS",
    "NO",
    "NO"
  ),
  reason = c(
    "PRJNA1010969 lacks longitudinal time structure in the current metadata.",
    "The cohort contains independent Control, Trauma, and Sepsis groups.",
    "Healthy-centroid distance is a cross-sectional group-state analogue only.",
    "Trauma provides a critical-illness comparator that helps assess specificity.",
    "Relative-abundance tests are exploratory and compositional.",
    "Alignment is cross-cohort concordance evidence only."
  )
)

write_csv(
  guardrails,
  file.path(OUT, "19_STEP93X_claim_guardrails.csv")
)

interpretation <- c(
  "STEP93X INTERPRETATION FRAMEWORK",
  "",
  paste0(
    "Abundance samples with group labels: ",
    nrow(mapping2),
    "."
  ),
  paste0(
    "Control / Trauma / Sepsis counts: ",
    paste(
      paste0(group_counts$Group, "=", group_counts$n),
      collapse = "; "
    ),
    "."
  ),
  "",
  "Primary external-validation questions:",
  "1. Does genus-level Bray community structure differ among Control, Trauma, and Sepsis?",
  "2. Are Sepsis samples farther from the healthy-control centroid than Control samples?",
  "3. Are Sepsis samples also different from Trauma, which would argue for specificity beyond generic severe illness?",
  "4. Do longitudinal genus changes in PRJNA691455 align directionally with cross-sectional Sepsis-state differences in PRJNA1010969?",
  "",
  "Interpretation guardrail:",
  "This cohort is cross-sectional and therefore cannot validate the temporal progression itself. It can validate whether the late/ill ecological state has an external cross-sectional analogue.",
  "",
  "Do not call genus-level nonparametric screens causal biomarkers; treat them as exploratory taxonomic concordance."
)

writeLines(
  interpretation,
  file.path(OUT, "20_STEP93X_INTERPRETATION.txt")
)

# ------------------------------------------------------------
# 11. QC
# ------------------------------------------------------------

qc <- tibble(
  project = PROJECT,
  abundance_rows = nrow(X),
  labeled_rows = nrow(mapping2),
  metadata_file = best_meta_file,
  metadata_run_overlap = meta_audit$exact_run_overlap[1],
  genus_features = ncol(X2),
  group_levels = paste(levels(group), collapse = ";"),
  Bray_PERMANOVA_p = perm_global$permutation_p[1],
  beta_dispersion_test_available = dispersion_result$available[1]
)

write_csv(
  qc,
  file.path(OUT, "21_STEP93X_QC_SUMMARY.csv")
)

writeLines(
  c(
    paste0("Completed: ", Sys.time()),
    paste0("Labeled samples: ", nrow(mapping2)),
    paste0("Genus features: ", ncol(X2)),
    "STEP93X COMPLETE"
  ),
  file.path(OUT, "_STEP93X_COMPLETE.ok")
)

cat("STEP93X COMPLETE\n")
