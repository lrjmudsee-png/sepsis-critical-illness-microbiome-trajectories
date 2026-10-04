# ============================================================
# Sepsis V2 - Step93X2
# PRJNA1010969 COMPLETE GROUP RECONSTRUCTION + EXTERNAL VALIDATION
#
# FIXES Step93X:
# 1) Step93X only mapped HMB/TMB/SMB and left BMT/P/PC unlabeled.
# 2) This version prioritizes source_material_id labels and uses
#    prefix mapping only as a fallback:
#       Control: HMB + PC
#       Trauma : TMB + BMT
#       Sepsis : SMB + P
# 3) Requires all 54 metadata samples to resolve to 18/18/18.
# 4) Re-runs all external-state analyses on all available abundance rows.
# 5) Canonicalizes genus names before longitudinal/external alignment.
#
# IMPORTANT:
# - This is cross-sectional external state validation.
# - PERMANOVA is interpreted together with beta-dispersion.
# - Genus-level screens and cross-cohort alignment remain exploratory.
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
  ROOT, "data", PROJECT, "03_dada2_rerun_v2",
  "PRJNA1010969_genus_relative_abundance.csv"
)

OUT <- file.path(
  ROOT, "results",
  "V2_33X2_PRJNA1010969_COMPLETE_GROUP_VALIDATION"
)

dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

logfile <- file.path(OUT, "_STEP93X2_runtime.txt")
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
  D <- matrix(0, n, n, dimnames = list(rownames(X), rownames(X)))
  if (n >= 2) {
    for (i in seq_len(n - 1)) {
      for (j in (i + 1):n) {
        d <- bray_vec(X[i, ], X[j, ])
        D[i, j] <- d
        D[j, i] <- d
      }
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

  if (!is.finite(within_ss) || within_ss <= 0) return(NA_real_)

  (between_ss / (k - 1)) / (within_ss / (n - k))
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

  perm <- replicate(B, pseudo_f_distance(D, sample(group)))
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

canon_genus <- function(x) {
  original <- as.character(x)
  z <- tolower(original)

  # Turn common make.names separators into spaces.
  z <- str_replace_all(z, "\\.", " ")
  z <- str_replace_all(z, "_", " ")

  # Keep only the genus component if taxonomy hierarchy is present.
  z <- str_replace(z, "^.*?genus\\s+", "")
  z <- str_replace(z, "\\s*\\|\\s*family.*$", "")
  z <- str_replace(z, "\\s+family\\s+.*$", "")

  # Strip generic taxonomic prefixes.
  z <- str_replace(z, "^genus\\s+", "")
  z <- str_squish(z)

  # Stable join key.
  key <- str_replace_all(z, "[^a-z0-9]+", "")

  key[key %in% c("", "na", "nan", "unknown", "unclassified")] <- NA_character_
  key
}

# ------------------------------------------------------------
# 1. Abundance
# ------------------------------------------------------------

if (!file.exists(ABUND_FILE)) stop("Missing PRJNA1010969 abundance file.")

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

if (is.na(run_col_abund)) stop("Cannot identify abundance Run column.")

run_ids <- norm_id(abund[[run_col_abund]])

num_cols <- names(abund)[vapply(abund, is.numeric, logical(1))]
if (length(num_cols) < 10) stop("Too few numeric genus columns.")

X <- as.matrix(abund[, num_cols, drop = FALSE])
storage.mode(X) <- "numeric"
X[!is.finite(X)] <- 0
X[X < 0] <- 0

keep <- rowSums(X) > 0
X <- X[keep, , drop = FALSE]
run_ids <- run_ids[keep]
X <- X / rowSums(X)
rownames(X) <- run_ids

# ------------------------------------------------------------
# 2. Metadata
# ------------------------------------------------------------

meta_candidates <- list.files(
  file.path(ROOT, "data", PROJECT),
  recursive = TRUE,
  full.names = TRUE,
  pattern = "SraRunTable.*\\.csv$",
  ignore.case = TRUE
)

if (length(meta_candidates) == 0) stop("No SraRunTable found.")

audit <- list()
loaded <- list()

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

  rc <- pick_exact(names(m), c("Run", "Run_ID", "run_id", "run"))
  if (is.na(rc)) next

  overlap <- length(intersect(unique(norm_id(m[[rc]])), unique(run_ids)))

  audit[[length(audit) + 1]] <- tibble(
    file = f,
    rows = nrow(m),
    run_column = rc,
    exact_run_overlap = overlap,
    columns = paste(names(m), collapse = ";")
  )
  loaded[[f]] <- m
}

audit <- bind_rows(audit) %>% arrange(desc(exact_run_overlap))
write_csv(audit, file.path(OUT, "01_metadata_candidate_audit.csv"))

if (nrow(audit) == 0 || audit$exact_run_overlap[1] == 0) {
  stop("No metadata table overlaps abundance runs.")
}

meta_file <- audit$file[1]
run_col <- audit$run_column[1]
meta <- loaded[[meta_file]]

sample_col <- pick_exact(
  names(meta),
  c("Sample Name", "Sample_Name", "sample_name", "SampleName")
)
if (is.na(sample_col)) {
  hit <- names(meta)[str_detect(names(meta), regex("sample.*name", ignore_case = TRUE))]
  if (length(hit) > 0) sample_col <- hit[1]
}
if (is.na(sample_col)) stop("No Sample Name column.")

source_col <- pick_exact(
  names(meta),
  c(
    "source_material_id",
    "source material identifiers",
    "source_material_identifiers",
    "Source material identifiers"
  )
)

if (is.na(source_col)) {
  hit <- names(meta)[
    str_detect(names(meta), regex("source.*material", ignore_case = TRUE))
  ]
  if (length(hit) > 0) source_col <- hit[1]
}

meta2 <- meta %>%
  transmute(
    Run_ID = norm_id(.data[[run_col]]),
    Sample_Name = as.character(.data[[sample_col]]),
    source_material_id = if (!is.na(source_col))
      as.character(.data[[source_col]]) else NA_character_
  ) %>%
  mutate(
    prefix = toupper(str_extract(Sample_Name, "^[A-Za-z]+")),

    Group_from_source = case_when(
      str_detect(tolower(source_material_id), "control|healthy") ~ "Control",
      str_detect(tolower(source_material_id), "trauma|injury") ~ "Trauma",
      str_detect(tolower(source_material_id), "sepsis") ~ "Sepsis",
      TRUE ~ NA_character_
    ),

    Group_from_prefix = case_when(
      prefix %in% c("HMB", "PC") ~ "Control",
      prefix %in% c("TMB", "BMT") ~ "Trauma",
      prefix %in% c("SMB", "P") ~ "Sepsis",
      TRUE ~ NA_character_
    ),

    Group = coalesce(Group_from_source, Group_from_prefix),

    source_prefix_agreement = case_when(
      is.na(Group_from_source) ~ NA,
      is.na(Group_from_prefix) ~ FALSE,
      TRUE ~ Group_from_source == Group_from_prefix
    )
  ) %>%
  distinct(Run_ID, .keep_all = TRUE)

write_csv(
  meta2,
  file.path(OUT, "02_COMPLETE_metadata_group_reconstruction.csv")
)

prefix_audit <- meta2 %>%
  count(prefix, Group, name = "n") %>%
  arrange(Group, prefix)

write_csv(
  prefix_audit,
  file.path(OUT, "03_prefix_group_audit.csv")
)

metadata_counts <- meta2 %>%
  count(Group, name = "n")

write_csv(
  metadata_counts,
  file.path(OUT, "04_FULL_METADATA_group_counts.csv")
)

# Hard validation against documented BioProject design: 18/18/18.
expected <- tibble(
  Group = c("Control", "Trauma", "Sepsis"),
  expected_n = c(18, 18, 18)
)

count_check <- expected %>%
  left_join(metadata_counts, by = "Group") %>%
  mutate(
    n = coalesce(n, 0L),
    matches_expected = n == expected_n
  )

write_csv(
  count_check,
  file.path(OUT, "05_METADATA_18_18_18_validation.csv")
)

if (
  any(is.na(meta2$Group)) ||
  !all(count_check$matches_expected)
) {
  writeLines(
    "Metadata grouping failed the expected 18/18/18 design. Inspect files 02-05.",
    file.path(OUT, "_STEP93X2_GROUP_VALIDATION_FAILED.txt")
  )
  stop("STEP93X2 group reconstruction did not reproduce 18/18/18.")
}

# ------------------------------------------------------------
# 3. Join all abundance samples to complete groups
# ------------------------------------------------------------

mapping <- tibble(Run_ID = rownames(X)) %>%
  left_join(
    meta2 %>% select(
      Run_ID, Sample_Name, source_material_id, prefix, Group
    ),
    by = "Run_ID"
  )

write_csv(
  mapping,
  file.path(OUT, "06_COMPLETE_Run_Group_mapping.csv")
)

if (any(is.na(mapping$Group))) {
  stop("Some abundance samples remain unlabeled after complete mapping.")
}

abundance_counts <- mapping %>%
  count(Group, name = "n")

write_csv(
  abundance_counts,
  file.path(OUT, "07_ABUNDANCE_group_counts.csv")
)

group <- factor(mapping$Group, levels = c("Control", "Trauma", "Sepsis"))
X2 <- X[mapping$Run_ID, , drop = FALSE]

logmsg(
  "Complete abundance groups: ",
  paste(paste0(abundance_counts$Group, "=", abundance_counts$n), collapse = "; ")
)

# ------------------------------------------------------------
# 4. Alpha
# ------------------------------------------------------------

Observed_Genera <- rowSums(X2 > 0)

Shannon <- apply(X2, 1, function(p) {
  p <- p[p > 0]
  -sum(p * log(p))
})

Simpson <- apply(X2, 1, function(p) {
  1 - sum(p^2)
})

alpha <- mapping %>%
  mutate(
    Observed_Genera = Observed_Genera,
    Shannon = Shannon,
    Simpson = Simpson
  )

write_csv(alpha, file.path(OUT, "08_alpha_by_sample.csv"))

alpha_global <- bind_rows(
  lapply(c("Observed_Genera", "Shannon", "Simpson"), function(metric) {
    kt <- kruskal.test(alpha[[metric]] ~ group)
    tibble(
      metric = metric,
      n = nrow(alpha),
      kruskal_chisq = unname(kt$statistic),
      df = unname(kt$parameter),
      p_value = kt$p.value
    )
  })
) %>%
  mutate(FDR = p.adjust(p_value, "BH"))

write_csv(
  alpha_global,
  file.path(OUT, "09_alpha_global_tests.csv")
)

contrasts <- list(
  c("Sepsis", "Control"),
  c("Sepsis", "Trauma"),
  c("Trauma", "Control")
)

alpha_pair <- bind_rows(
  lapply(c("Observed_Genera", "Shannon", "Simpson"), function(metric) {
    bind_rows(lapply(contrasts, function(cc) {
      a <- alpha[[metric]][group == cc[1]]
      b <- alpha[[metric]][group == cc[2]]

      tibble(
        metric = metric,
        group_A = cc[1],
        group_B = cc[2],
        n_A = length(a),
        n_B = length(b),
        median_A = median(a, na.rm = TRUE),
        median_B = median(b, na.rm = TRUE),
        median_difference_A_minus_B =
          median(a, na.rm = TRUE) - median(b, na.rm = TRUE),
        cliff_delta_A_vs_B = cliff_delta(a, b),
        wilcoxon_p = safe_wilcox_p(a, b)
      )
    }))
  })
) %>%
  mutate(FDR = p.adjust(wilcoxon_p, "BH"))

write_csv(
  alpha_pair,
  file.path(OUT, "10_alpha_pairwise_tests.csv")
)

# ------------------------------------------------------------
# 5. Bray + PERMANOVA
# ------------------------------------------------------------

D <- bray_matrix(X2)

perm_global <- permanova_custom(D, group, B = 9999)

write_csv(
  perm_global,
  file.path(OUT, "11_Bray_PERMANOVA_global.csv")
)

perm_pair <- bind_rows(
  lapply(contrasts, function(cc) {
    ii <- which(group %in% cc)
    permanova_custom(
      D[ii, ii, drop = FALSE],
      droplevels(group[ii]),
      B = 9999
    ) %>%
      mutate(group_A = cc[1], group_B = cc[2])
  })
) %>%
  mutate(FDR = p.adjust(permutation_p, "BH"))

write_csv(
  perm_pair,
  file.path(OUT, "12_Bray_PERMANOVA_pairwise.csv")
)

# ------------------------------------------------------------
# 6. PCoA
# ------------------------------------------------------------

pc <- cmdscale(as.dist(D), k = 2, eig = TRUE, add = TRUE)

coords <- as.data.frame(pc$points)
names(coords)[1:2] <- c("PCoA1", "PCoA2")

pcoa <- mapping %>% bind_cols(coords)

eig <- pc$eig
pos_eig <- eig[eig > 0]
v1 <- if (length(pos_eig) >= 1) pos_eig[1] / sum(pos_eig) else NA_real_
v2 <- if (length(pos_eig) >= 2) pos_eig[2] / sum(pos_eig) else NA_real_

write_csv(pcoa, file.path(OUT, "13_Bray_PCoA_coordinates.csv"))

write_csv(
  tibble(
    axis = c("PCoA1", "PCoA2"),
    variance_fraction = c(v1, v2)
  ),
  file.path(OUT, "14_Bray_PCoA_variance.csv")
)

# ------------------------------------------------------------
# 7. Distance to control centroid
# ------------------------------------------------------------

control_centroid <- colMeans(X2[group == "Control", , drop = FALSE])
control_centroid <- control_centroid / sum(control_centroid)

centroid_dist <- apply(
  X2, 1, function(z) bray_vec(z, control_centroid)
)

centroid_df <- mapping %>%
  mutate(Bray_to_Control_Centroid = centroid_dist)

write_csv(
  centroid_df,
  file.path(OUT, "15_distance_to_Control_centroid.csv")
)

kw <- kruskal.test(
  centroid_df$Bray_to_Control_Centroid ~ group
)

write_csv(
  tibble(
    n = nrow(centroid_df),
    kruskal_chisq = unname(kw$statistic),
    df = unname(kw$parameter),
    p_value = kw$p.value
  ),
  file.path(OUT, "16_Control_centroid_global_test.csv")
)

centroid_pair <- bind_rows(
  lapply(contrasts, function(cc) {
    a <- centroid_df$Bray_to_Control_Centroid[group == cc[1]]
    b <- centroid_df$Bray_to_Control_Centroid[group == cc[2]]

    tibble(
      group_A = cc[1],
      group_B = cc[2],
      n_A = length(a),
      n_B = length(b),
      median_A = median(a, na.rm = TRUE),
      median_B = median(b, na.rm = TRUE),
      median_difference_A_minus_B =
        median(a, na.rm = TRUE) - median(b, na.rm = TRUE),
      cliff_delta_A_vs_B = cliff_delta(a, b),
      wilcoxon_p = safe_wilcox_p(a, b)
    )
  })
) %>%
  mutate(FDR = p.adjust(wilcoxon_p, "BH"))

write_csv(
  centroid_pair,
  file.path(OUT, "17_Control_centroid_pairwise_tests.csv")
)

# ------------------------------------------------------------
# 8. Dispersion
# ------------------------------------------------------------

dispersion <- tibble(
  method = "vegan_betadisper_if_available",
  available = FALSE,
  F = NA_real_,
  permutation_p = NA_real_,
  permutations = 9999
)

if (requireNamespace("vegan", quietly = TRUE)) {
  bd <- vegan::betadisper(as.dist(D), group)
  pt <- vegan::permutest(bd, permutations = 9999)

  dispersion <- tibble(
    method = "vegan::betadisper + permutest",
    available = TRUE,
    F = unname(pt$tab[1, "F"]),
    permutation_p = unname(pt$tab[1, "Pr(>F)"]),
    permutations = 9999
  )
}

write_csv(
  dispersion,
  file.path(OUT, "18_beta_dispersion_test.csv")
)

# ------------------------------------------------------------
# 9. Exploratory genus screen
# ------------------------------------------------------------

prevalence <- colMeans(X2 > 0)
keep_genus <- names(prevalence)[prevalence >= 0.10]

genus_screen <- bind_rows(
  lapply(keep_genus, function(g) {
    y <- X2[, g]

    co <- y[group == "Control"]
    tr <- y[group == "Trauma"]
    se <- y[group == "Sepsis"]

    kt <- tryCatch(
      kruskal.test(y ~ group),
      error = function(e) NULL
    )

    tibble(
      genus = g,
      genus_key = canon_genus(g),
      prevalence = prevalence[g],
      median_Control = median(co),
      median_Trauma = median(tr),
      median_Sepsis = median(se),
      mean_Control = mean(co),
      mean_Trauma = mean(tr),
      mean_Sepsis = mean(se),
      KW_p = if (is.null(kt)) NA_real_ else kt$p.value,
      Sepsis_minus_Control = median(se) - median(co),
      Sepsis_minus_Control_mean = mean(se) - mean(co),
      Sepsis_vs_Control_cliff = cliff_delta(se, co),
      Sepsis_vs_Control_p = safe_wilcox_p(se, co),
      Sepsis_minus_Trauma = median(se) - median(tr),
      Sepsis_minus_Trauma_mean = mean(se) - mean(tr),
      Sepsis_vs_Trauma_cliff = cliff_delta(se, tr),
      Sepsis_vs_Trauma_p = safe_wilcox_p(se, tr),
      Trauma_minus_Control = median(tr) - median(co),
      Trauma_minus_Control_mean = mean(tr) - mean(co),
      Trauma_vs_Control_cliff = cliff_delta(tr, co),
      Trauma_vs_Control_p = safe_wilcox_p(tr, co)
    )
  })
) %>%
  mutate(
    KW_FDR = p.adjust(KW_p, "BH"),
    Sepsis_vs_Control_FDR = p.adjust(Sepsis_vs_Control_p, "BH"),
    Sepsis_vs_Trauma_FDR = p.adjust(Sepsis_vs_Trauma_p, "BH"),
    Trauma_vs_Control_FDR = p.adjust(Trauma_vs_Control_p, "BH")
  ) %>%
  arrange(KW_FDR, KW_p)

write_csv(
  genus_screen,
  file.path(OUT, "19_exploratory_genus_group_screen.csv")
)

# ------------------------------------------------------------
# 10. Corrected longitudinal/external genus-name alignment
# ------------------------------------------------------------

signed_file <- file.path(
  ROOT, "results",
  "V2_33S2_TRAJECTORY_DISTANCE_METHOD_FIX",
  "signed_trajectory_matrix_clean.csv"
)

alignment_summary <- tibble()

if (file.exists(signed_file)) {

  signed <- read_csv(
    signed_file,
    show_col_types = FALSE,
    name_repair = "unique"
  )

  pcols <- c("Patient_true", "patient_true", "Patient", "patient_id")
  pcoll <- pcols[pcols %in% names(signed)][1]
  if (is.na(pcoll)) pcoll <- names(signed)[1]

  num_signed <- signed %>% select(where(is.numeric))

  long_mean <- colMeans(
    as.matrix(num_signed),
    na.rm = TRUE
  )

  longitudinal <- tibble(
    longitudinal_feature = names(long_mean),
    genus_key = canon_genus(names(long_mean)),
    PRJNA691455_mean_signed_delta = as.numeric(long_mean)
  ) %>%
    filter(!is.na(genus_key)) %>%
    group_by(genus_key) %>%
    summarise(
      PRJNA691455_mean_signed_delta =
        mean(PRJNA691455_mean_signed_delta, na.rm = TRUE),
      longitudinal_feature =
        paste(unique(longitudinal_feature), collapse = " | "),
      .groups = "drop"
    )

  external_for_join <- genus_screen %>%
    filter(!is.na(genus_key)) %>%
    group_by(genus_key) %>%
    summarise(
      external_genus = paste(unique(genus), collapse = " | "),
      Sepsis_minus_Control =
        mean(Sepsis_minus_Control_mean, na.rm = TRUE),
      Sepsis_minus_Trauma =
        mean(Sepsis_minus_Trauma_mean, na.rm = TRUE),
      .groups = "drop"
    )

  alignment <- longitudinal %>%
    inner_join(external_for_join, by = "genus_key")

  write_csv(
    alignment,
    file.path(OUT, "20_longitudinal_external_genus_alignment.csv")
  )

  if (nrow(alignment) >= 5) {

    sc <- suppressWarnings(
      cor.test(
        alignment$PRJNA691455_mean_signed_delta,
        alignment$Sepsis_minus_Control,
        method = "spearman",
        exact = FALSE
      )
    )

    st <- suppressWarnings(
      cor.test(
        alignment$PRJNA691455_mean_signed_delta,
        alignment$Sepsis_minus_Trauma,
        method = "spearman",
        exact = FALSE
      )
    )

    alignment_summary <- tibble(
      comparison = c(
        "Longitudinal delta vs external Sepsis-Control",
        "Longitudinal delta vs external Sepsis-Trauma"
      ),
      overlapping_genera = nrow(alignment),
      spearman_rho = c(unname(sc$estimate), unname(st$estimate)),
      p_value = c(sc$p.value, st$p.value),
      sign_concordance = c(
        mean(
          sign(alignment$PRJNA691455_mean_signed_delta) ==
            sign(alignment$Sepsis_minus_Control),
          na.rm = TRUE
        ),
        mean(
          sign(alignment$PRJNA691455_mean_signed_delta) ==
            sign(alignment$Sepsis_minus_Trauma),
          na.rm = TRUE
        )
      )
    )
  }

  write_csv(
    alignment_summary,
    file.path(OUT, "21_longitudinal_external_alignment_summary.csv")
  )
}

# ------------------------------------------------------------
# 11. Figure
# ------------------------------------------------------------

make_fig <- function() {
  par(mfrow = c(2, 2), mar = c(4.5, 4.5, 2.7, 1))

  pch_map <- c(Control = 1, Trauma = 17, Sepsis = 19)

  plot(
    pcoa$PCoA1,
    pcoa$PCoA2,
    pch = pch_map[pcoa$Group],
    xlab = paste0("PCoA1 (", round(100*v1, 1), "%)"),
    ylab = paste0("PCoA2 (", round(100*v2, 1), "%)"),
    main = "A. Complete-cohort Bray PCoA"
  )
  legend("topright", legend = names(pch_map), pch = pch_map, bty = "n")

  boxplot(
    Bray_to_Control_Centroid ~ Group,
    data = centroid_df,
    ylab = "Bray-Curtis to healthy centroid",
    xlab = "",
    main = "B. Distance from healthy ecological state"
  )

  boxplot(
    Shannon ~ Group,
    data = alpha,
    ylab = "Genus-level Shannon",
    xlab = "",
    main = "C. Alpha diversity"
  )

  if (exists("alignment") && nrow(alignment) >= 5) {
    plot(
      alignment$PRJNA691455_mean_signed_delta,
      alignment$Sepsis_minus_Control,
      pch = 19,
      xlab = "PRJNA691455 mean longitudinal genus delta",
      ylab = "PRJNA1010969 Sepsis-Control mean abundance difference",
      main = "D. Longitudinal / external alignment"
    )
    abline(h = 0, v = 0, lty = 2)
  } else {
    plot.new()
    title("D. Genus alignment unavailable")
  }
}

pdf(
  file.path(OUT, "Figure_STEP93X2_complete_external_validation.pdf"),
  width = 10, height = 8
)
make_fig()
dev.off()

png(
  file.path(OUT, "Figure_STEP93X2_complete_external_validation.png"),
  width = 1800, height = 1400, res = 180
)
make_fig()
dev.off()

# ------------------------------------------------------------
# 12. Interpretation + QC
# ------------------------------------------------------------

dispersion_warning <- (
  nrow(dispersion) == 1 &&
  isTRUE(dispersion$available[1]) &&
  is.finite(dispersion$permutation_p[1]) &&
  dispersion$permutation_p[1] < 0.05
)

interpretation <- c(
  "STEP93X2 COMPLETE EXTERNAL VALIDATION INTERPRETATION",
  "",
  paste0(
    "Full metadata design: ",
    paste(
      paste0(metadata_counts$Group, "=", metadata_counts$n),
      collapse = "; "
    )
  ),
  paste0(
    "Available abundance samples: ",
    paste(
      paste0(abundance_counts$Group, "=", abundance_counts$n),
      collapse = "; "
    )
  ),
  "",
  "Primary rules:",
  "1. Use Step93X2 only; discard Step93X statistical results because Step93X analyzed only 24/53 labeled abundance samples.",
  "2. PRJNA1010969 is cross-sectional external state validation, not longitudinal replication.",
  "3. Sepsis-vs-Trauma is the key specificity comparison beyond generic critical illness.",
  "4. PERMANOVA must be interpreted together with beta-dispersion.",
  if (dispersion_warning)
    "5. Beta-dispersion is significant; therefore PERMANOVA separation may reflect both centroid/location and dispersion differences."
  else
    "5. No significant beta-dispersion warning was detected.",
  "6. Genus-level screens and cross-cohort taxonomic alignment remain exploratory."
)

writeLines(
  interpretation,
  file.path(OUT, "22_STEP93X2_INTERPRETATION.txt")
)

qc <- tibble(
  project = PROJECT,
  metadata_rows = nrow(meta2),
  metadata_Control = sum(meta2$Group == "Control"),
  metadata_Trauma = sum(meta2$Group == "Trauma"),
  metadata_Sepsis = sum(meta2$Group == "Sepsis"),
  abundance_rows = nrow(mapping),
  abundance_Control = sum(mapping$Group == "Control"),
  abundance_Trauma = sum(mapping$Group == "Trauma"),
  abundance_Sepsis = sum(mapping$Group == "Sepsis"),
  metadata_source_column = ifelse(is.na(source_col), "NONE", source_col),
  source_prefix_disagreements =
    sum(meta2$source_prefix_agreement == FALSE, na.rm = TRUE),
  Bray_PERMANOVA_p = perm_global$permutation_p[1],
  beta_dispersion_p =
    ifelse(isTRUE(dispersion$available[1]),
           dispersion$permutation_p[1], NA_real_),
  aligned_genera =
    ifelse(exists("alignment"), nrow(alignment), NA_integer_)
)

write_csv(
  qc,
  file.path(OUT, "23_STEP93X2_QC_SUMMARY.csv")
)

writeLines(
  c(
    paste0("Completed: ", Sys.time()),
    "Full metadata group reconstruction validated at 18/18/18.",
    paste0("Abundance samples analyzed: ", nrow(mapping)),
    "STEP93X2 COMPLETE"
  ),
  file.path(OUT, "_STEP93X2_COMPLETE.ok")
)

cat("STEP93X2 COMPLETE\n")
