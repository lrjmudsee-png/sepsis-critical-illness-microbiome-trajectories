# ============================================================
# Sepsis V2 - Step93X4
# EXTERNAL VALIDATION EVIDENCE FREEZE
#
# Purpose:
# Correct Step93X3's evidence-tier counting bug and freeze the
# manuscript-facing interpretation of PRJNA1010969.
#
# Key correction:
# Aitchison pseudocount variants are sensitivity variants of ONE
# distance family, not three independent distance tests.
#
# Distance families:
#   1) Bray-Curtis
#   2) Hellinger-Euclidean
#   3) CLR/Aitchison
#
# Frozen interpretation:
# - Bray and Hellinger Sepsis-vs-Trauma PERMANOVA are significant
#   but pairwise dispersion is also significant.
# - CLR/Aitchison Sepsis-vs-Trauma is significant across all three
#   pseudocounts and pairwise dispersion is non-significant.
# - All three alpha-diversity metrics differ Sepsis vs Trauma.
# - Sepsis is farther from the healthy-control centroid than Trauma.
# - 14 exploratory genera differ by FDR.
# - Longitudinal genus change aligns with Sepsis-Control, but NOT
#   Sepsis-Trauma, so no replicated sepsis-specific taxonomic signature.
# ============================================================

options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(stringr)
  library(tibble)
})

ROOT <- "E:/sepsis_project"

IN_X3 <- file.path(
  ROOT,
  "results",
  "V2_33X3_EXTERNAL_SPECIFICITY_ROBUSTNESS"
)

IN_X2 <- file.path(
  ROOT,
  "results",
  "V2_33X2_PRJNA1010969_COMPLETE_GROUP_VALIDATION"
)

OUT <- file.path(
  ROOT,
  "results",
  "V2_33X4_EXTERNAL_VALIDATION_EVIDENCE_FREEZE"
)

dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

required <- c(
  "04_Bray_pairwise_PERMANOVA.csv",
  "05_Bray_pairwise_dispersion.csv",
  "08_Hellinger_pairwise_PERMANOVA.csv",
  "09_Hellinger_pairwise_dispersion.csv",
  "12_CLR_pairwise_PERMANOVA_pseudocount_sensitivity.csv",
  "13_CLR_pairwise_dispersion_pseudocount_sensitivity.csv",
  "14_alpha_pairwise_from_X2.csv",
  "15_Control_centroid_pairwise_from_X2.csv",
  "16_longitudinal_external_alignment_from_X2.csv",
  "17_genus_specificity_counts.csv",
  "18_Sepsis_vs_Trauma_FDR_significant_genera.csv"
)

for (f in required) {
  if (!file.exists(file.path(IN_X3, f))) {
    stop(paste0("Missing Step93X3 source: ", f))
  }
}

bray <- read_csv(
  file.path(IN_X3, "04_Bray_pairwise_PERMANOVA.csv"),
  show_col_types = FALSE
)

bray_disp <- read_csv(
  file.path(IN_X3, "05_Bray_pairwise_dispersion.csv"),
  show_col_types = FALSE
)

hell <- read_csv(
  file.path(IN_X3, "08_Hellinger_pairwise_PERMANOVA.csv"),
  show_col_types = FALSE
)

hell_disp <- read_csv(
  file.path(IN_X3, "09_Hellinger_pairwise_dispersion.csv"),
  show_col_types = FALSE
)

clr <- read_csv(
  file.path(IN_X3, "12_CLR_pairwise_PERMANOVA_pseudocount_sensitivity.csv"),
  show_col_types = FALSE
)

clr_disp <- read_csv(
  file.path(IN_X3, "13_CLR_pairwise_dispersion_pseudocount_sensitivity.csv"),
  show_col_types = FALSE
)

alpha <- read_csv(
  file.path(IN_X3, "14_alpha_pairwise_from_X2.csv"),
  show_col_types = FALSE
)

centroid <- read_csv(
  file.path(IN_X3, "15_Control_centroid_pairwise_from_X2.csv"),
  show_col_types = FALSE
)

alignment <- read_csv(
  file.path(IN_X3, "16_longitudinal_external_alignment_from_X2.csv"),
  show_col_types = FALSE
)

genus_counts <- read_csv(
  file.path(IN_X3, "17_genus_specificity_counts.csv"),
  show_col_types = FALSE
)

genus_st <- read_csv(
  file.path(IN_X3, "18_Sepsis_vs_Trauma_FDR_significant_genera.csv"),
  show_col_types = FALSE
)

get_pair <- function(df, A="Sepsis", B="Trauma") {
  df %>%
    filter(group_A == A, group_B == B)
}

# ------------------------------------------------------------
# 1. Collapse multivariate evidence by distance FAMILY
# ------------------------------------------------------------

br_st <- get_pair(bray) %>% slice(1)
brd_st <- get_pair(bray_disp) %>% slice(1)

he_st <- get_pair(hell) %>% slice(1)
hed_st <- get_pair(hell_disp) %>% slice(1)

clr_st <- get_pair(clr)
clrd_st <- get_pair(clr_disp)

if (nrow(clr_st) == 0 || nrow(clrd_st) == 0) {
  stop("Missing Sepsis-vs-Trauma CLR sensitivity rows.")
}

aitchison_all_sig <- all(clr_st$FDR < 0.05, na.rm = TRUE)

aitchison_all_disp_nonsig <- all(
  clrd_st$FDR >= 0.05,
  na.rm = TRUE
)

aitchison_r2_range <- range(
  clr_st$R2,
  na.rm = TRUE
)

aitchison_fdr_range <- range(
  clr_st$FDR,
  na.rm = TRUE
)

aitchison_disp_fdr_range <- range(
  clrd_st$FDR,
  na.rm = TRUE
)

family_matrix <- tibble(
  distance_family = c(
    "Bray-Curtis",
    "Hellinger-Euclidean",
    "CLR/Aitchison"
  ),
  variants = c(
    1,
    1,
    nrow(clr_st)
  ),
  Sepsis_vs_Trauma_significant = c(
    br_st$FDR[1] < 0.05,
    he_st$FDR[1] < 0.05,
    aitchison_all_sig
  ),
  pairwise_dispersion_significant = c(
    brd_st$FDR[1] < 0.05,
    hed_st$FDR[1] < 0.05,
    !aitchison_all_disp_nonsig
  ),
  clean_location_evidence = c(
    br_st$FDR[1] < 0.05 && brd_st$FDR[1] >= 0.05,
    he_st$FDR[1] < 0.05 && hed_st$FDR[1] >= 0.05,
    aitchison_all_sig && aitchison_all_disp_nonsig
  ),
  R2_or_range = c(
    formatC(br_st$R2[1], digits=4, format="f"),
    formatC(he_st$R2[1], digits=4, format="f"),
    paste0(
      formatC(aitchison_r2_range[1], digits=4, format="f"),
      "–",
      formatC(aitchison_r2_range[2], digits=4, format="f")
    )
  ),
  PERMANOVA_FDR_or_range = c(
    formatC(br_st$FDR[1], digits=4, format="f"),
    formatC(he_st$FDR[1], digits=4, format="f"),
    paste0(
      formatC(aitchison_fdr_range[1], digits=4, format="f"),
      "–",
      formatC(aitchison_fdr_range[2], digits=4, format="f")
    )
  ),
  dispersion_FDR_or_range = c(
    formatC(brd_st$FDR[1], digits=4, format="f"),
    formatC(hed_st$FDR[1], digits=4, format="f"),
    paste0(
      formatC(aitchison_disp_fdr_range[1], digits=4, format="f"),
      "–",
      formatC(aitchison_disp_fdr_range[2], digits=4, format="f")
    )
  )
)

write_csv(
  family_matrix,
  file.path(OUT, "01_FROZEN_multivariate_distance_family_matrix.csv")
)

# ------------------------------------------------------------
# 2. Scalar ecological evidence
# ------------------------------------------------------------

alpha_st <- alpha %>%
  filter(group_A == "Sepsis", group_B == "Trauma")

centroid_st <- centroid %>%
  filter(group_A == "Sepsis", group_B == "Trauma")

scalar_matrix <- bind_rows(
  alpha_st %>%
    transmute(
      evidence = paste0("Alpha_", metric),
      effect = cliff_delta_A_vs_B,
      FDR = FDR,
      direction = ifelse(
        median_difference_A_minus_B < 0,
        "Sepsis_lower_than_Trauma",
        "Sepsis_higher_than_Trauma"
      )
    ),

  centroid_st %>%
    transmute(
      evidence = "Distance_to_healthy_control_centroid",
      effect = cliff_delta_A_vs_B,
      FDR = FDR,
      direction = ifelse(
        median_difference_A_minus_B > 0,
        "Sepsis_farther_than_Trauma",
        "Sepsis_not_farther_than_Trauma"
      )
    )
)

write_csv(
  scalar_matrix,
  file.path(OUT, "02_FROZEN_scalar_ecological_specificity.csv")
)

# ------------------------------------------------------------
# 3. Taxonomic evidence
# ------------------------------------------------------------

taxonomic_summary <- tibble(
  item = c(
    "Screened genera",
    "Sepsis-vs-Trauma FDR-significant genera",
    "Longitudinal-vs-Sepsis-Control overlapping genera",
    "Longitudinal-vs-Sepsis-Control Spearman rho",
    "Longitudinal-vs-Sepsis-Control p",
    "Longitudinal-vs-Sepsis-Trauma overlapping genera",
    "Longitudinal-vs-Sepsis-Trauma Spearman rho",
    "Longitudinal-vs-Sepsis-Trauma p"
  ),
  value = c(
    genus_counts$screened_genera[
      genus_counts$contrast == "Sepsis_vs_Trauma"
    ][1],
    genus_counts$FDR_significant_genera[
      genus_counts$contrast == "Sepsis_vs_Trauma"
    ][1],
    alignment$overlapping_genera[
      str_detect(alignment$comparison, "Sepsis-Control")
    ][1],
    alignment$spearman_rho[
      str_detect(alignment$comparison, "Sepsis-Control")
    ][1],
    alignment$p_value[
      str_detect(alignment$comparison, "Sepsis-Control")
    ][1],
    alignment$overlapping_genera[
      str_detect(alignment$comparison, "Sepsis-Trauma")
    ][1],
    alignment$spearman_rho[
      str_detect(alignment$comparison, "Sepsis-Trauma")
    ][1],
    alignment$p_value[
      str_detect(alignment$comparison, "Sepsis-Trauma")
    ][1]
  )
)

write_csv(
  taxonomic_summary,
  file.path(OUT, "03_FROZEN_taxonomic_external_summary.csv")
)

# Keep only named genera for a cleaner manuscript candidate list.
named_genus_candidates <- genus_st %>%
  filter(!is.na(genus_key), genus_key != "") %>%
  select(
    genus,
    prevalence,
    median_Control,
    median_Trauma,
    median_Sepsis,
    Sepsis_minus_Trauma_mean,
    Sepsis_vs_Trauma_cliff,
    Sepsis_vs_Trauma_FDR,
    Sepsis_vs_Control_FDR,
    Trauma_vs_Control_FDR
  ) %>%
  arrange(Sepsis_vs_Trauma_FDR)

write_csv(
  named_genus_candidates,
  file.path(OUT, "04_EXPLORATORY_named_Sepsis_vs_Trauma_genus_candidates.csv")
)

# ------------------------------------------------------------
# 4. Corrected evidence tier
# ------------------------------------------------------------

n_distance_families <- nrow(family_matrix)
n_clean_distance_families <- sum(
  family_matrix$clean_location_evidence
)

all_alpha_sig <- (
  nrow(alpha_st) == 3 &&
  all(alpha_st$FDR < 0.05)
)

centroid_sig <- (
  nrow(centroid_st) == 1 &&
  centroid_st$FDR[1] < 0.05
)

sepsis_trauma_taxonomic_replication <- (
  alignment$spearman_rho[
    str_detect(alignment$comparison, "Sepsis-Trauma")
  ][1] > 0.2 &&
  alignment$p_value[
    str_detect(alignment$comparison, "Sepsis-Trauma")
  ][1] < 0.05
)

# Conservative tier:
# one clean distance family + strong scalar evidence =
# supportive state differentiation, not a validated specific signature.
corrected_tier <- case_when(

  n_clean_distance_families >= 2 &&
    all_alpha_sig &&
    centroid_sig &&
    sepsis_trauma_taxonomic_replication ~
    "MODERATE_REPLICATED_SEPSIS_ASSOCIATED_SPECIFICITY",

  n_clean_distance_families >= 1 &&
    all_alpha_sig &&
    centroid_sig ~
    "SUPPORTIVE_EXTERNAL_SEPSIS_ASSOCIATED_ECOLOGICAL_STATE_DIFFERENTIATION",

  all_alpha_sig &&
    centroid_sig ~
    "SCALAR_EXTERNAL_SEPSIS_ASSOCIATED_DIFFERENTIATION",

  TRUE ~
    "NO_ROBUST_EXTERNAL_SEPSIS_ASSOCIATED_DIFFERENTIATION"
)

tier <- tibble(
  corrected_external_evidence_tier = corrected_tier,
  independent_distance_families_tested = n_distance_families,
  clean_multivariate_distance_families = n_clean_distance_families,
  clean_family_names = paste(
    family_matrix$distance_family[
      family_matrix$clean_location_evidence
    ],
    collapse = ";"
  ),
  all_three_alpha_metrics_significant_Sepsis_vs_Trauma =
    all_alpha_sig,
  Sepsis_farther_than_Trauma_from_healthy_centroid =
    centroid_sig,
  Sepsis_vs_Trauma_named_FDR_genera =
    nrow(named_genus_candidates),
  longitudinal_Sepsis_Control_taxonomic_alignment =
    alignment$spearman_rho[
      str_detect(alignment$comparison, "Sepsis-Control")
    ][1],
  longitudinal_Sepsis_Trauma_taxonomic_alignment =
    alignment$spearman_rho[
      str_detect(alignment$comparison, "Sepsis-Trauma")
    ][1],
  replicated_sepsis_specific_taxonomic_signature =
    sepsis_trauma_taxonomic_replication
)

write_csv(
  tier,
  file.path(OUT, "05_CORRECTED_EXTERNAL_EVIDENCE_TIER.csv")
)

# ------------------------------------------------------------
# 5. Frozen claim guardrails
# ------------------------------------------------------------

guardrails <- tibble(
  claim = c(
    "PRJNA1010969 externally replicates the longitudinal temporal trajectory",
    "Sepsis differs ecologically from Trauma in the external cohort",
    "Bray-Curtis alone proves centroid separation between Sepsis and Trauma",
    "Aitchison composition differs between Sepsis and Trauma robustly across pseudocounts",
    "Sepsis shows lower alpha diversity than Trauma",
    "Sepsis is farther from the healthy-control centroid than Trauma",
    "A replicated sepsis-specific taxonomic signature is established",
    "The external cohort supports a graded ecological severity axis",
    "The 14 genus hits are validated sepsis biomarkers"
  ),
  status = c(
    "NO",
    "YES_SUPPORTIVE",
    "NO",
    "YES",
    "YES",
    "YES",
    "NO",
    "YES",
    "NO"
  ),
  reason = c(
    "PRJNA1010969 is cross-sectional.",
    "Aitchison separation, alpha diversity, and healthy-centroid distance all distinguish Sepsis from Trauma.",
    "Bray PERMANOVA is accompanied by significant pairwise dispersion.",
    "All three pseudocount variants are significant and pairwise dispersion is non-significant.",
    "Observed genera, Shannon, and Simpson are all FDR-significant.",
    "The Sepsis-Trauma healthy-centroid distance contrast is FDR-significant.",
    "Longitudinal genus change does not align with external Sepsis-Trauma differences.",
    "Control, Trauma, and Sepsis show ordered ecological deterioration across multiple scalar summaries.",
    "They arise from exploratory compositional relative-abundance screening."
  )
)

write_csv(
  guardrails,
  file.path(OUT, "06_FROZEN_external_claim_guardrails.csv")
)

# ------------------------------------------------------------
# 6. Manuscript-ready results paragraph
# ------------------------------------------------------------

a_obs <- alpha_st %>% filter(metric == "Observed_Genera")
a_sha <- alpha_st %>% filter(metric == "Shannon")
a_sim <- alpha_st %>% filter(metric == "Simpson")

sc_align <- alignment %>%
  filter(str_detect(comparison, "Sepsis-Control")) %>%
  slice(1)

st_align <- alignment %>%
  filter(str_detect(comparison, "Sepsis-Trauma")) %>%
  slice(1)

results_text <- c(
  "FROZEN EXTERNAL VALIDATION RESULTS",
  "",
  paste0(
    "In the independent PRJNA1010969 cohort (17 controls, 18 trauma patients, and 18 sepsis patients with genus-level profiles), sepsis showed lower ecological diversity than trauma across observed genera (FDR=",
    formatC(a_obs$FDR, digits=4, format="f"),
    "), Shannon diversity (FDR=",
    formatC(a_sha$FDR, digits=4, format="f"),
    "), and Simpson diversity (FDR=",
    formatC(a_sim$FDR, digits=4, format="f"),
    ")."
  ),
  paste0(
    "Sepsis samples were also farther from the healthy-control centroid than trauma samples (Cliff's delta=",
    formatC(centroid_st$cliff_delta_A_vs_B[1], digits=3, format="f"),
    ", FDR=",
    formatC(centroid_st$FDR[1], digits=4, format="f"),
    ")."
  ),
  paste0(
    "Bray-Curtis and Hellinger PERMANOVA detected Sepsis-Trauma differences (FDR=",
    formatC(br_st$FDR[1], digits=4, format="f"),
    " and ",
    formatC(he_st$FDR[1], digits=4, format="f"),
    ", respectively), but both contrasts were accompanied by significant pairwise dispersion."
  ),
  paste0(
    "In contrast, CLR/Aitchison analysis remained significant across pseudocounts 10^-6 to 10^-4 (R2=",
    formatC(aitchison_r2_range[1], digits=3, format="f"),
    "–",
    formatC(aitchison_r2_range[2], digits=3, format="f"),
    "; FDR=",
    formatC(aitchison_fdr_range[1], digits=4, format="f"),
    "–",
    formatC(aitchison_fdr_range[2], digits=4, format="f"),
    ") without significant pairwise dispersion."
  ),
  paste0(
    "Across 127 overlapping genera, longitudinal change in PRJNA691455 correlated with the external Sepsis-Control contrast (rho=",
    formatC(sc_align$spearman_rho, digits=3, format="f"),
    ", p=",
    formatC(sc_align$p_value, digits=4, format="f"),
    ") but not with the Sepsis-Trauma contrast (rho=",
    formatC(st_align$spearman_rho, digits=3, format="f"),
    ", p=",
    formatC(st_align$p_value, digits=3, format="f"),
    ")."
  ),
  "",
  paste0(
    "Frozen interpretation tier: ",
    corrected_tier,
    "."
  ),
  "",
  "Interpretation:",
  "The external cohort supports sepsis-associated ecological state differentiation beyond severe trauma, but does not establish a reproducible sepsis-specific taxonomic signature."
)

writeLines(
  results_text,
  file.path(OUT, "07_FROZEN_external_results_paragraph.txt")
)

# ------------------------------------------------------------
# 7. Manuscript-ready methods wording
# ------------------------------------------------------------

methods_text <- c(
  "FROZEN EXTERNAL VALIDATION METHODS",
  "",
  "PRJNA1010969 was treated as an independent cross-sectional external-state cohort rather than as a longitudinal replication cohort.",
  "Group structure was reconstructed as Control, Trauma, and Sepsis using source-material metadata with sample-prefix fallback and was validated against the documented 18/18/18 metadata design.",
  "Genus-level community differences were evaluated using Bray-Curtis, Hellinger-Euclidean, and CLR/Aitchison distances.",
  "For CLR/Aitchison analyses, pseudocounts of 10^-6, 10^-5, and 10^-4 were used as a zero-handling sensitivity analysis.",
  "PERMANOVA results were interpreted jointly with pairwise beta-dispersion tests; distance-family variants were treated as sensitivity analyses within one family rather than as independent replications.",
  "Alpha-diversity and healthy-centroid distance contrasts were tested nonparametrically with false-discovery-rate correction.",
  "Genus-level relative-abundance contrasts and longitudinal/external taxonomic alignment were considered exploratory."
)

writeLines(
  methods_text,
  file.path(OUT, "08_FROZEN_external_methods_wording.txt")
)

# ------------------------------------------------------------
# 8. Figure
# ------------------------------------------------------------

make_fig <- function() {

  par(mfrow=c(2,2), mar=c(4.7,4.5,2.7,1))

  # A: distance families
  vals <- c(
    Bray = br_st$FDR[1],
    Hellinger = he_st$FDR[1],
    Aitchison_min = min(clr_st$FDR)
  )

  barplot(
    -log10(vals),
    names.arg = names(vals),
    las=2,
    ylab="-log10(PERMANOVA FDR)",
    main="A. Sepsis vs Trauma multivariate evidence"
  )
  abline(h=-log10(0.05), lty=2)

  # B: dispersion
  dvals <- c(
    Bray = brd_st$FDR[1],
    Hellinger = hed_st$FDR[1],
    Aitchison_max = max(clrd_st$FDR)
  )

  barplot(
    -log10(dvals),
    names.arg=names(dvals),
    las=2,
    ylab="-log10(pairwise dispersion FDR)",
    main="B. Dispersion sensitivity"
  )
  abline(h=-log10(0.05), lty=2)

  # C: alpha effects
  barplot(
    alpha_st$cliff_delta_A_vs_B,
    names.arg=alpha_st$metric,
    las=2,
    ylab="Cliff's delta (Sepsis vs Trauma)",
    main="C. Alpha-diversity depletion"
  )
  abline(h=0, lty=2)

  # D: alignment
  barplot(
    c(sc_align$spearman_rho, st_align$spearman_rho),
    names.arg=c("Sepsis-Control","Sepsis-Trauma"),
    ylab="Spearman rho",
    main="D. Longitudinal/external taxonomic alignment"
  )
  abline(h=0, lty=2)
}

pdf(
  file.path(OUT, "Figure_STEP93X4_external_validation_freeze.pdf"),
  width=10,
  height=8
)
make_fig()
dev.off()

png(
  file.path(OUT, "Figure_STEP93X4_external_validation_freeze.png"),
  width=1800,
  height=1400,
  res=180
)
make_fig()
dev.off()

# ------------------------------------------------------------
# 9. QC
# ------------------------------------------------------------

qc <- tibble(
  distance_families_tested = n_distance_families,
  clean_distance_families = n_clean_distance_families,
  clean_distance_family = paste(
    family_matrix$distance_family[
      family_matrix$clean_location_evidence
    ],
    collapse=";"
  ),
  all_alpha_Sepsis_vs_Trauma_significant = all_alpha_sig,
  healthy_centroid_Sepsis_vs_Trauma_significant = centroid_sig,
  named_Sepsis_vs_Trauma_genus_candidates = nrow(named_genus_candidates),
  longitudinal_Sepsis_Control_rho = sc_align$spearman_rho,
  longitudinal_Sepsis_Trauma_rho = st_align$spearman_rho,
  corrected_tier = corrected_tier
)

write_csv(
  qc,
  file.path(OUT, "09_STEP93X4_QC_SUMMARY.csv")
)

writeLines(
  c(
    paste0("Completed: ", Sys.time()),
    paste0("Corrected tier: ", corrected_tier),
    "Step93X3 raw statistics retained; X3 automatic tier superseded.",
    "STEP93X4 COMPLETE"
  ),
  file.path(OUT, "_STEP93X4_COMPLETE.ok")
)

cat("STEP93X4 COMPLETE\n")
