# ============================================================
# Sepsis V2 - Step93X3
# PRJNA1010969 EXTERNAL SPECIFICITY ROBUSTNESS
#
# Why:
# Step93X2 established the correct complete groups:
#   Control 17 abundance samples
#   Trauma  18
#   Sepsis  18
#
# But global Bray beta-dispersion was significant.
# Therefore Bray PERMANOVA alone cannot be interpreted as pure centroid
# separation. This step asks whether Sepsis-vs-Trauma separation is robust
# across:
#
#   1) pairwise Bray PERMANOVA + pairwise beta-dispersion
#   2) Hellinger-Euclidean distance
#   3) CLR/Aitchison distance under 3 pseudocount sensitivities
#   4) alpha-diversity specificity
#   5) distance-from-healthy-centroid specificity
#   6) genus-level specificity counts
#   7) longitudinal/external taxonomic alignment
#
# Output is a manuscript-facing specificity evidence matrix.
# ============================================================

options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(stringr)
  library(tibble)
})

if (!requireNamespace("vegan", quietly = TRUE)) {
  stop("Package 'vegan' is required for Step93X3.")
}

set.seed(20260822)

ROOT <- "E:/sepsis_project"

IN_X2 <- file.path(
  ROOT,
  "results",
  "V2_33X2_PRJNA1010969_COMPLETE_GROUP_VALIDATION"
)

ABUND_FILE <- file.path(
  ROOT,
  "data",
  "PRJNA1010969",
  "03_dada2_rerun_v2",
  "PRJNA1010969_genus_relative_abundance.csv"
)

OUT <- file.path(
  ROOT,
  "results",
  "V2_33X3_EXTERNAL_SPECIFICITY_ROBUSTNESS"
)

dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

logfile <- file.path(OUT, "_STEP93X3_runtime.txt")
writeLines(paste("START", Sys.time()), logfile)

logmsg <- function(...) {
  z <- paste0(...)
  cat(z, "\n")
  cat(z, "\n", file = logfile, append = TRUE)
}

pick_exact <- function(nms, candidates) {
  x <- candidates[candidates %in% nms]
  if (length(x) == 0) return(NA_character_)
  x[1]
}

norm_id <- function(x) {
  x <- toupper(str_trim(as.character(x)))
  x[x %in% c("", "NA", "N/A", "NULL")] <- NA_character_
  x
}

pairwise_adonis <- function(D, group, permutations = 9999) {

  contrasts <- list(
    c("Sepsis", "Control"),
    c("Sepsis", "Trauma"),
    c("Trauma", "Control")
  )

  bind_rows(lapply(contrasts, function(cc) {

    ii <- which(group %in% cc)

    dg <- as.dist(
      as.matrix(D)[ii, ii, drop = FALSE]
    )

    g <- droplevels(group[ii])

    md <- data.frame(Group = g)

    fit <- vegan::adonis2(
      dg ~ Group,
      data = md,
      permutations = permutations
    )

    tibble(
      group_A = cc[1],
      group_B = cc[2],
      n = length(ii),
      pseudo_F = unname(fit$F[1]),
      R2 = unname(fit$R2[1]),
      p_value = unname(fit$`Pr(>F)`[1])
    )
  })) %>%
    mutate(FDR = p.adjust(p_value, "BH"))
}

pairwise_betadisper <- function(D, group, permutations = 9999) {

  contrasts <- list(
    c("Sepsis", "Control"),
    c("Sepsis", "Trauma"),
    c("Trauma", "Control")
  )

  bind_rows(lapply(contrasts, function(cc) {

    ii <- which(group %in% cc)
    dg <- as.dist(
      as.matrix(D)[ii, ii, drop = FALSE]
    )
    g <- droplevels(group[ii])

    bd <- vegan::betadisper(dg, g)
    pt <- vegan::permutest(
      bd,
      permutations = permutations
    )

    med <- tapply(
      bd$distances,
      g,
      median,
      na.rm = TRUE
    )

    tibble(
      group_A = cc[1],
      group_B = cc[2],
      n = length(ii),
      median_dispersion_A = unname(med[cc[1]]),
      median_dispersion_B = unname(med[cc[2]]),
      F = unname(pt$tab[1, "F"]),
      p_value = unname(pt$tab[1, "Pr(>F)"])
    )
  })) %>%
    mutate(FDR = p.adjust(p_value, "BH"))
}

global_adonis <- function(D, group, permutations = 9999) {
  md <- data.frame(Group = group)

  fit <- vegan::adonis2(
    D ~ Group,
    data = md,
    permutations = permutations
  )

  tibble(
    n = length(group),
    groups = nlevels(group),
    pseudo_F = unname(fit$F[1]),
    R2 = unname(fit$R2[1]),
    p_value = unname(fit$`Pr(>F)`[1]),
    permutations = permutations
  )
}

global_betadisper <- function(D, group, permutations = 9999) {
  bd <- vegan::betadisper(D, group)
  pt <- vegan::permutest(
    bd,
    permutations = permutations
  )

  tibble(
    n = length(group),
    F = unname(pt$tab[1, "F"]),
    p_value = unname(pt$tab[1, "Pr(>F)"]),
    permutations = permutations
  )
}

# ------------------------------------------------------------
# 1. Load corrected Step93X2 mapping + abundance
# ------------------------------------------------------------

mapping_file <- file.path(
  IN_X2,
  "06_COMPLETE_Run_Group_mapping.csv"
)

if (!file.exists(mapping_file)) {
  stop("Missing Step93X2 complete Run/Group mapping.")
}

mapping <- read_csv(
  mapping_file,
  show_col_types = FALSE
)

if (!file.exists(ABUND_FILE)) {
  stop("Missing PRJNA1010969 genus abundance.")
}

abund <- read_csv(
  ABUND_FILE,
  show_col_types = FALSE,
  guess_max = 10000,
  name_repair = "unique"
)

run_col <- pick_exact(
  names(abund),
  c("Run_ID", "Run", "run_id", "run", "RunID")
)

if (is.na(run_col)) stop("Cannot detect abundance Run column.")

abund$Run_key <- norm_id(abund[[run_col]])
mapping$Run_key <- norm_id(mapping$Run_ID)

num_cols <- names(abund)[
  vapply(abund, is.numeric, logical(1))
]

X <- as.matrix(
  abund[
    match(mapping$Run_key, abund$Run_key),
    num_cols,
    drop = FALSE
  ]
)

storage.mode(X) <- "numeric"
X[!is.finite(X)] <- 0
X[X < 0] <- 0

if (any(rowSums(X) <= 0)) {
  stop("One or more mapped abundance rows have zero total abundance.")
}

X <- X / rowSums(X)
rownames(X) <- mapping$Run_ID

group <- factor(
  mapping$Group,
  levels = c("Control", "Trauma", "Sepsis")
)

write_csv(
  mapping %>% count(Group, name = "n"),
  file.path(OUT, "01_group_counts_confirmed.csv")
)

# ------------------------------------------------------------
# 2. Bray-Curtis: global + pairwise PERMANOVA and dispersion
# ------------------------------------------------------------

D_bray <- vegan::vegdist(
  X,
  method = "bray"
)

bray_global <- global_adonis(D_bray, group)
bray_disp_global <- global_betadisper(D_bray, group)

bray_pair <- pairwise_adonis(D_bray, group)
bray_disp_pair <- pairwise_betadisper(D_bray, group)

write_csv(
  bray_global,
  file.path(OUT, "02_Bray_global_PERMANOVA.csv")
)

write_csv(
  bray_disp_global,
  file.path(OUT, "03_Bray_global_dispersion.csv")
)

write_csv(
  bray_pair,
  file.path(OUT, "04_Bray_pairwise_PERMANOVA.csv")
)

write_csv(
  bray_disp_pair,
  file.path(OUT, "05_Bray_pairwise_dispersion.csv")
)

# ------------------------------------------------------------
# 3. Hellinger-Euclidean sensitivity
# ------------------------------------------------------------

X_hel <- sqrt(X)
D_hel <- dist(X_hel, method = "euclidean")

hel_global <- global_adonis(D_hel, group)
hel_disp_global <- global_betadisper(D_hel, group)
hel_pair <- pairwise_adonis(D_hel, group)
hel_disp_pair <- pairwise_betadisper(D_hel, group)

write_csv(
  hel_global,
  file.path(OUT, "06_Hellinger_global_PERMANOVA.csv")
)

write_csv(
  hel_disp_global,
  file.path(OUT, "07_Hellinger_global_dispersion.csv")
)

write_csv(
  hel_pair,
  file.path(OUT, "08_Hellinger_pairwise_PERMANOVA.csv")
)

write_csv(
  hel_disp_pair,
  file.path(OUT, "09_Hellinger_pairwise_dispersion.csv")
)

# ------------------------------------------------------------
# 4. CLR / Aitchison sensitivity under multiple pseudocounts
# ------------------------------------------------------------

pseudocounts <- c(
  1e-6,
  1e-5,
  1e-4
)

clr_results <- list()
clr_pair_results <- list()
clr_disp_results <- list()
clr_pair_disp_results <- list()

for (pc in pseudocounts) {

  X_pc <- X + pc
  X_pc <- X_pc / rowSums(X_pc)

  logX <- log(X_pc)
  clr <- logX - rowMeans(logX)

  D_clr <- dist(clr, method = "euclidean")

  g <- global_adonis(D_clr, group) %>%
    mutate(pseudocount = pc)

  gd <- global_betadisper(D_clr, group) %>%
    mutate(pseudocount = pc)

  p <- pairwise_adonis(D_clr, group) %>%
    mutate(pseudocount = pc)

  pd <- pairwise_betadisper(D_clr, group) %>%
    mutate(pseudocount = pc)

  clr_results[[length(clr_results)+1]] <- g
  clr_disp_results[[length(clr_disp_results)+1]] <- gd
  clr_pair_results[[length(clr_pair_results)+1]] <- p
  clr_pair_disp_results[[length(clr_pair_disp_results)+1]] <- pd
}

clr_global <- bind_rows(clr_results)
clr_disp_global <- bind_rows(clr_disp_results)
clr_pair <- bind_rows(clr_pair_results)
clr_pair_disp <- bind_rows(clr_pair_disp_results)

write_csv(
  clr_global,
  file.path(OUT, "10_CLR_global_PERMANOVA_pseudocount_sensitivity.csv")
)

write_csv(
  clr_disp_global,
  file.path(OUT, "11_CLR_global_dispersion_pseudocount_sensitivity.csv")
)

write_csv(
  clr_pair,
  file.path(OUT, "12_CLR_pairwise_PERMANOVA_pseudocount_sensitivity.csv")
)

write_csv(
  clr_pair_disp,
  file.path(OUT, "13_CLR_pairwise_dispersion_pseudocount_sensitivity.csv")
)

# ------------------------------------------------------------
# 5. Import scalar specificity evidence from Step93X2
# ------------------------------------------------------------

alpha_pair <- read_csv(
  file.path(IN_X2, "10_alpha_pairwise_tests.csv"),
  show_col_types = FALSE
)

centroid_pair <- read_csv(
  file.path(IN_X2, "17_Control_centroid_pairwise_tests.csv"),
  show_col_types = FALSE
)

genus_screen <- read_csv(
  file.path(IN_X2, "19_exploratory_genus_group_screen.csv"),
  show_col_types = FALSE
)

alignment_summary <- read_csv(
  file.path(IN_X2, "21_longitudinal_external_alignment_summary.csv"),
  show_col_types = FALSE
)

write_csv(
  alpha_pair,
  file.path(OUT, "14_alpha_pairwise_from_X2.csv")
)

write_csv(
  centroid_pair,
  file.path(OUT, "15_Control_centroid_pairwise_from_X2.csv")
)

write_csv(
  alignment_summary,
  file.path(OUT, "16_longitudinal_external_alignment_from_X2.csv")
)

# ------------------------------------------------------------
# 6. Genus-level specificity counts
# ------------------------------------------------------------

genus_counts <- tibble(
  contrast = c(
    "Sepsis_vs_Control",
    "Sepsis_vs_Trauma",
    "Trauma_vs_Control"
  ),
  FDR_significant_genera = c(
    sum(genus_screen$Sepsis_vs_Control_FDR < 0.05, na.rm = TRUE),
    sum(genus_screen$Sepsis_vs_Trauma_FDR < 0.05, na.rm = TRUE),
    sum(genus_screen$Trauma_vs_Control_FDR < 0.05, na.rm = TRUE)
  ),
  screened_genera = nrow(genus_screen)
)

write_csv(
  genus_counts,
  file.path(OUT, "17_genus_specificity_counts.csv")
)

sepsis_trauma_genus <- genus_screen %>%
  filter(
    !is.na(Sepsis_vs_Trauma_FDR),
    Sepsis_vs_Trauma_FDR < 0.05
  ) %>%
  arrange(Sepsis_vs_Trauma_FDR)

write_csv(
  sepsis_trauma_genus,
  file.path(OUT, "18_Sepsis_vs_Trauma_FDR_significant_genera.csv")
)

# ------------------------------------------------------------
# 7. Build Sepsis-vs-Trauma robustness matrix
# ------------------------------------------------------------

get_pair <- function(df, A = "Sepsis", B = "Trauma") {
  df %>%
    filter(group_A == A, group_B == B) %>%
    slice(1)
}

br_st <- get_pair(bray_pair)
brd_st <- get_pair(bray_disp_pair)

he_st <- get_pair(hel_pair)
hed_st <- get_pair(hel_disp_pair)

clr_st <- clr_pair %>%
  filter(group_A == "Sepsis", group_B == "Trauma") %>%
  arrange(pseudocount)

clrd_st <- clr_pair_disp %>%
  filter(group_A == "Sepsis", group_B == "Trauma") %>%
  arrange(pseudocount)

alpha_st <- alpha_pair %>%
  filter(group_A == "Sepsis", group_B == "Trauma")

centroid_st <- centroid_pair %>%
  filter(group_A == "Sepsis", group_B == "Trauma")

robustness <- bind_rows(
  tibble(
    evidence_type = "Bray PERMANOVA",
    variant = "Bray-Curtis",
    effect_or_R2 = br_st$R2,
    p_value = br_st$p_value,
    FDR = br_st$FDR,
    dispersion_p_or_FDR = brd_st$FDR,
    robust_location_evidence =
      br_st$FDR < 0.05 && brd_st$FDR >= 0.05
  ),

  tibble(
    evidence_type = "Hellinger PERMANOVA",
    variant = "Hellinger-Euclidean",
    effect_or_R2 = he_st$R2,
    p_value = he_st$p_value,
    FDR = he_st$FDR,
    dispersion_p_or_FDR = hed_st$FDR,
    robust_location_evidence =
      he_st$FDR < 0.05 && hed_st$FDR >= 0.05
  ),

  bind_rows(lapply(seq_len(nrow(clr_st)), function(i) {
    tibble(
      evidence_type = "CLR PERMANOVA",
      variant = paste0("Aitchison_pc_", clr_st$pseudocount[i]),
      effect_or_R2 = clr_st$R2[i],
      p_value = clr_st$p_value[i],
      FDR = clr_st$FDR[i],
      dispersion_p_or_FDR =
        clrd_st$FDR[
          match(
            clr_st$pseudocount[i],
            clrd_st$pseudocount
          )
        ],
      robust_location_evidence =
        clr_st$FDR[i] < 0.05 &&
        clrd_st$FDR[
          match(
            clr_st$pseudocount[i],
            clrd_st$pseudocount
          )
        ] >= 0.05
    )
  })),

  alpha_st %>%
    transmute(
      evidence_type = "Alpha diversity",
      variant = metric,
      effect_or_R2 = cliff_delta_A_vs_B,
      p_value = wilcoxon_p,
      FDR = FDR,
      dispersion_p_or_FDR = NA_real_,
      robust_location_evidence = FDR < 0.05
    ),

  centroid_st %>%
    transmute(
      evidence_type = "Distance to healthy centroid",
      variant = "Bray_to_Control_Centroid",
      effect_or_R2 = cliff_delta_A_vs_B,
      p_value = wilcoxon_p,
      FDR = FDR,
      dispersion_p_or_FDR = NA_real_,
      robust_location_evidence = FDR < 0.05
    )
)

write_csv(
  robustness,
  file.path(OUT, "19_SEPSIS_vs_TRAUMA_specificity_robustness_matrix.csv")
)

# ------------------------------------------------------------
# 8. Formal evidence tier
# ------------------------------------------------------------

multi_distance_significant <- sum(
  robustness$evidence_type %in%
    c("Bray PERMANOVA", "Hellinger PERMANOVA", "CLR PERMANOVA") &
  robustness$FDR < 0.05,
  na.rm = TRUE
)

multi_distance_unconfounded <- sum(
  robustness$evidence_type %in%
    c("Bray PERMANOVA", "Hellinger PERMANOVA", "CLR PERMANOVA") &
  robustness$robust_location_evidence,
  na.rm = TRUE
)

all_alpha_specific <- all(
  alpha_st$FDR < 0.05
)

centroid_specific <- (
  nrow(centroid_st) == 1 &&
  centroid_st$FDR[1] < 0.05
)

specificity_tier <- case_when(
  multi_distance_unconfounded >= 2 &&
    all_alpha_specific &&
    centroid_specific ~
    "MODERATE_EXTERNAL_SEPSIS_ASSOCIATED_SPECIFICITY",

  multi_distance_significant >= 2 &&
    all_alpha_specific &&
    centroid_specific ~
    "SUPPORTIVE_SPECIFICITY_WITH_DISPERSION_CAVEAT",

  all_alpha_specific &&
    centroid_specific ~
    "SCALAR_ECOLOGICAL_SPECIFICITY_ONLY",

  TRUE ~
    "NO_ROBUST_SEPSIS_SPECIFIC_EXTERNAL_STATE"
)

tier <- tibble(
  specificity_tier = specificity_tier,
  multivariate_distance_tests_significant = multi_distance_significant,
  multivariate_distance_tests_unconfounded_by_pairwise_dispersion =
    multi_distance_unconfounded,
  all_three_alpha_metrics_Sepsis_vs_Trauma_FDR_lt_005 =
    all_alpha_specific,
  Sepsis_farther_than_Trauma_from_Control_centroid =
    centroid_specific,
  Sepsis_vs_Trauma_FDR_significant_genera =
    sum(genus_screen$Sepsis_vs_Trauma_FDR < 0.05, na.rm = TRUE),
  longitudinal_external_Sepsis_Control_rho =
    alignment_summary$spearman_rho[
      str_detect(alignment_summary$comparison, "Sepsis-Control")
    ][1],
  longitudinal_external_Sepsis_Trauma_rho =
    alignment_summary$spearman_rho[
      str_detect(alignment_summary$comparison, "Sepsis-Trauma")
    ][1]
)

write_csv(
  tier,
  file.path(OUT, "20_EXTERNAL_SPECIFICITY_EVIDENCE_TIER.csv")
)

# ------------------------------------------------------------
# 9. Manuscript interpretation guardrails
# ------------------------------------------------------------

interpretation <- c(
  "STEP93X3 EXTERNAL SPECIFICITY ROBUSTNESS",
  "",
  paste0("Formal external specificity tier: ", specificity_tier),
  "",
  "Facts already established from corrected Step93X2:",
  "- Metadata design is 18 Control / 18 Trauma / 18 Sepsis.",
  "- Available genus abundance is 17 Control / 18 Trauma / 18 Sepsis.",
  "- Sepsis has significantly lower Observed genera, Shannon, and Simpson than Trauma.",
  "- Sepsis is significantly farther from the healthy-control centroid than Trauma.",
  "- Longitudinal genus change aligns with external Sepsis-vs-Control differences but not with Sepsis-vs-Trauma differences.",
  "",
  "Interpretation rules:",
  "1. A significant PERMANOVA accompanied by significant pairwise beta-dispersion is not treated as clean centroid/location evidence.",
  "2. Sepsis-vs-Trauma alpha-diversity differences provide independent scalar ecological specificity evidence.",
  "3. Distance to the healthy-control centroid supports a graded Healthy -> Trauma -> Sepsis ecological displacement axis, but it remains cross-sectional.",
  "4. Genus-level Sepsis-vs-Trauma findings are exploratory because they arise from relative-abundance screening.",
  "5. The weak longitudinal-vs-Sepsis-Trauma genus correlation means the longitudinal taxonomic trajectory should not be presented as a replicated sepsis-specific taxonomic signature.",
  "6. The stronger longitudinal-vs-Sepsis-Control correlation is more consistent with a general healthy-to-critical-illness ecological axis."
)

writeLines(
  interpretation,
  file.path(OUT, "21_STEP93X3_INTERPRETATION.txt")
)

# ------------------------------------------------------------
# 10. Figure
# ------------------------------------------------------------

make_fig <- function() {

  par(mfrow = c(2, 2), mar = c(4.5, 4.5, 2.8, 1))

  # A: multivariate FDR
  mv <- robustness %>%
    filter(
      evidence_type %in%
        c("Bray PERMANOVA", "Hellinger PERMANOVA", "CLR PERMANOVA")
    )

  y <- -log10(pmax(mv$FDR, 1e-6))

  barplot(
    y,
    names.arg = mv$variant,
    las = 2,
    ylab = "-log10(FDR)",
    main = "A. Sepsis vs Trauma multivariate sensitivity"
  )
  abline(h = -log10(0.05), lty = 2)

  # B: pairwise dispersion
  dd <- c(
    Bray = brd_st$FDR,
    Hellinger = hed_st$FDR,
    CLR_1e6 =
      clrd_st$FDR[match(1e-6, clrd_st$pseudocount)],
    CLR_1e5 =
      clrd_st$FDR[match(1e-5, clrd_st$pseudocount)],
    CLR_1e4 =
      clrd_st$FDR[match(1e-4, clrd_st$pseudocount)]
  )

  barplot(
    -log10(pmax(dd, 1e-6)),
    names.arg = names(dd),
    las = 2,
    ylab = "-log10(dispersion FDR)",
    main = "B. Sepsis vs Trauma dispersion"
  )
  abline(h = -log10(0.05), lty = 2)

  # C: alpha Cliff delta
  aa <- alpha_st
  barplot(
    aa$cliff_delta_A_vs_B,
    names.arg = aa$metric,
    las = 2,
    ylab = "Cliff delta (Sepsis vs Trauma)",
    main = "C. Alpha-diversity specificity"
  )
  abline(h = 0, lty = 2)

  # D: cross-cohort taxonomic alignment
  if (nrow(alignment_summary) >= 2) {
    barplot(
      alignment_summary$spearman_rho,
      names.arg = c("Sepsis-Control", "Sepsis-Trauma"),
      ylim = range(c(0, alignment_summary$spearman_rho), na.rm = TRUE),
      ylab = "Spearman rho",
      main = "D. Longitudinal/external genus alignment"
    )
    abline(h = 0, lty = 2)
  } else {
    plot.new()
    title("D. Alignment unavailable")
  }
}

pdf(
  file.path(OUT, "Figure_STEP93X3_external_specificity_robustness.pdf"),
  width = 10,
  height = 8
)
make_fig()
dev.off()

png(
  file.path(OUT, "Figure_STEP93X3_external_specificity_robustness.png"),
  width = 1800,
  height = 1400,
  res = 180
)
make_fig()
dev.off()

# ------------------------------------------------------------
# 11. QC
# ------------------------------------------------------------

qc <- tibble(
  abundance_samples = nrow(X),
  Control = sum(group == "Control"),
  Trauma = sum(group == "Trauma"),
  Sepsis = sum(group == "Sepsis"),
  Bray_Sepsis_Trauma_FDR = br_st$FDR,
  Bray_Sepsis_Trauma_dispersion_FDR = brd_st$FDR,
  Hellinger_Sepsis_Trauma_FDR = he_st$FDR,
  Hellinger_Sepsis_Trauma_dispersion_FDR = hed_st$FDR,
  CLR_variants_tested = length(pseudocounts),
  specificity_tier = specificity_tier
)

write_csv(
  qc,
  file.path(OUT, "22_STEP93X3_QC_SUMMARY.csv")
)

writeLines(
  c(
    paste0("Completed: ", Sys.time()),
    paste0("Specificity tier: ", specificity_tier),
    "STEP93X3 COMPLETE"
  ),
  file.path(OUT, "_STEP93X3_COMPLETE.ok")
)

cat("STEP93X3 COMPLETE\n")
