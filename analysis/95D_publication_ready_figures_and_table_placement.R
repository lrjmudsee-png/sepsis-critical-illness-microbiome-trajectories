# ============================================================
# Sepsis V2 - Step95D
# PUBLICATION-READY FIGURE 1-5 REBUILD + TABLE PLACEMENT
#
# IMPORTANT:
# - NO new inferential analysis.
# - All values are read from already-frozen Step93W/93U/93X4/94B2D
#   and Step95C sources.
# - Step94B3C genus findings remain supplementary/exploratory.
# ============================================================

options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(stringr)
  library(tibble)
})

ROOT <- "E:/sepsis_project"
RESULTS <- file.path(ROOT, "results")

D93W <- file.path(RESULTS, "V2_33W_LONGITUDINAL_EVIDENCE_FREEZE")
D93U <- file.path(RESULTS, "V2_33U_TRAJECTORY_ECOLOGY_ROBUSTNESS")
D93X4 <- file.path(RESULTS, "V2_33X4_EXTERNAL_VALIDATION_EVIDENCE_FREEZE")
D94B2D <- file.path(RESULTS, "V2_34B2D_INFECTION_SOURCE_EVIDENCE_FREEZE")
D94B3C <- file.path(RESULTS, "V2_34B3C_CRA002354_GENUS_EVIDENCE_FREEZE")
D95C <- file.path(RESULTS, "V2_35C_MANUSCRIPT_ARCHITECTURE_FIGURE1_AND_MAIN_TABLES")

OUT <- file.path(
  RESULTS,
  "V2_35D_PUBLICATION_READY_FIGURES_AND_TABLE_PLACEMENT"
)

dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

FIGDIR <- file.path(OUT, "01_FINAL_FIGURES")
TABDIR <- file.path(OUT, "02_TABLE_PLACEMENT")
CAPDIR <- file.path(OUT, "03_CAPTIONS_AND_NOTES")

dir.create(FIGDIR, recursive = TRUE, showWarnings = FALSE)
dir.create(TABDIR, recursive = TRUE, showWarnings = FALSE)
dir.create(CAPDIR, recursive = TRUE, showWarnings = FALSE)

required_flags <- c(
  file.path(D93W, "_STEP93W_COMPLETE.ok"),
  file.path(D93U, "_STEP93U_COMPLETE.ok"),
  file.path(D93X4, "_STEP93X4_COMPLETE.ok"),
  file.path(D94B2D, "_STEP94B2D_COMPLETE.ok"),
  file.path(D94B3C, "_STEP94B3C_COMPLETE.ok"),
  file.path(D95C, "_STEP95C_COMPLETE.ok")
)

if (!all(file.exists(required_flags))) {
  stop("One or more required frozen branches / Step95C are incomplete.")
}

# ------------------------------------------------------------
# Plot helpers
# ------------------------------------------------------------

open_pdf <- function(path, width = 11, height = 7.5) {
  pdf(path, width = width, height = height, family = "sans", useDingbats = FALSE)
}

open_png <- function(path, width = 2200, height = 1500, res = 200) {
  png(path, width = width, height = height, res = res)
}

panel_label <- function(x, y, lab) {
  text(x, y, lab, font = 2, cex = 1.0, adj = c(0,1), xpd = NA)
}

safe_num <- function(x, digits = 3) {
  ifelse(is.na(x), "NA", formatC(x, digits = digits, format = "f"))
}

# ------------------------------------------------------------
# Load frozen data
# ------------------------------------------------------------

W <- read_csv(
  file.path(D93W, "01_FROZEN_longitudinal_evidence_matrix.csv"),
  show_col_types = FALSE
)

U1 <- read_csv(
  file.path(D93U, "01_patient_trajectory_ecology_endpoints.csv"),
  show_col_types = FALSE
)

U3 <- read_csv(
  file.path(D93U, "03_paired_Day3_Day7_Bray.csv"),
  show_col_types = FALSE
)

U5 <- read_csv(
  file.path(D93U, "05_pairwise_signed_trajectory_heterogeneity.csv"),
  show_col_types = FALSE
)

U9 <- read_csv(
  file.path(D93U, "09_trajectory_ecology_robustness.csv"),
  show_col_types = FALSE
)

X1 <- read_csv(
  file.path(D93X4, "01_FROZEN_multivariate_distance_family_matrix.csv"),
  show_col_types = FALSE
)

X2 <- read_csv(
  file.path(D93X4, "02_FROZEN_scalar_ecological_specificity.csv"),
  show_col_types = FALSE
)

X3 <- read_csv(
  file.path(D93X4, "03_FROZEN_taxonomic_external_summary.csv"),
  show_col_types = FALSE
)

B07 <- read_csv(
  file.path(D94B2D, "SOURCE_07_PRIMARY_BRAY_DISPLACEMENT_INTERACTION.csv"),
  show_col_types = FALSE
)

B08 <- read_csv(
  file.path(D94B2D, "SOURCE_08_PRIMARY_SOURCE_SPECIFIC_BRAY_SLOPES.csv"),
  show_col_types = FALSE
)

B20 <- read_csv(
  file.path(D94B2D, "SOURCE_20_KEY_BRAY_SENSITIVITY_SUMMARY.csv"),
  show_col_types = FALSE
)

B04 <- read_csv(
  file.path(D94B2D, "04_BASELINE_BETA_ROBUSTNESS_FREEZE.csv"),
  show_col_types = FALSE
)

B01 <- read_csv(
  file.path(D94B2D, "01_FINAL_INFECTION_SOURCE_EVIDENCE_REGISTRY.csv"),
  show_col_types = FALSE
)

GFREEZE <- read_csv(
  file.path(D94B3C, "03_FINAL_GENUS_BRANCH_EVIDENCE_FREEZE.csv"),
  show_col_types = FALSE
)

# ------------------------------------------------------------
# FIGURE 1
# Study architecture - remove internal Step IDs from manuscript visual
# and show that core sepsis is nested within longitudinal evidence.
# ------------------------------------------------------------

make_fig1 <- function(device = c("pdf","png")) {

  device <- match.arg(device)

  path <- file.path(
    FIGDIR,
    paste0("Figure_1_Study_Architecture_Publication.", device)
  )

  if (device == "pdf") open_pdf(path, 11.5, 7.2) else open_png(path, 2300, 1450, 200)

  par(mar = c(0.5,0.5,1.0,0.5))

  plot(
    NA, xlim = c(0,12), ylim = c(0,8),
    axes = FALSE, xlab = "", ylab = "",
    xaxs = "i", yaxs = "i"
  )

  rect(4.0, 7.0, 8.0, 7.7, lwd = 1.3)
  text(6.0, 7.35, "Public gut-microbiome cohorts in critical illness", font = 2, cex = 1.0)

  # Longitudinal framework
  rect(0.6, 5.1, 6.0, 6.3, lwd = 1.3)
  text(3.3, 5.92, "Longitudinal critical-illness framework", font = 2, cex = 0.95)
  text(3.3, 5.48, "7 cohorts; within-patient ecological change", cex = 0.86)

  # Nested core sepsis
  rect(1.2, 3.25, 5.4, 4.45, lwd = 1.2)
  text(3.3, 4.05, "Core repeated-sepsis cohort", font = 2, cex = 0.92)
  text(3.3, 3.60, "Taxonomic-route heterogeneity + ecological displacement", cex = 0.80)

  # External
  rect(6.5, 5.1, 11.4, 6.3, lwd = 1.3)
  text(8.95, 5.92, "External ecological-state validation", font = 2, cex = 0.94)
  text(8.95, 5.48, "Control vs severe trauma vs sepsis", cex = 0.86)

  # Source
  rect(6.5, 3.25, 11.4, 4.45, lwd = 1.3)
  text(8.95, 4.05, "Infection-source longitudinal extension", font = 2, cex = 0.92)
  text(8.95, 3.60, "Pulmonary vs recorded non-pulmonary sepsis", cex = 0.82)

  arrows(6, 7.0, 3.3, 6.3, length = 0.08)
  arrows(6, 7.0, 8.95, 6.3, length = 0.08)
  arrows(3.3, 5.1, 3.3, 4.45, length = 0.08)
  arrows(8.95, 5.1, 8.95, 4.45, length = 0.08)

  rect(1.3, 0.75, 10.7, 2.0, lwd = 1.4)
  text(
    6.0, 1.52,
    "Central synthesis: progressive ecological displacement can be shared across cohorts",
    font = 2, cex = 0.92
  )
  text(
    6.0, 1.10,
    "while taxonomic routes remain heterogeneous and infection-source specificity is limited.",
    cex = 0.88
  )

  arrows(3.3, 3.25, 4.8, 2.0, length = 0.08)
  arrows(8.95, 3.25, 7.2, 2.0, length = 0.08)
  arrows(8.95, 5.1, 8.0, 2.0, length = 0.08)

  dev.off()
}

make_fig1("pdf")
make_fig1("png")

# ------------------------------------------------------------
# FIGURE 2
# Cross-cohort primary vs common-anchor standardized Bray change
# ------------------------------------------------------------

make_fig2 <- function(device = c("pdf","png")) {

  device <- match.arg(device)
  path <- file.path(
    FIGDIR,
    paste0("Figure_2_Cross_Cohort_Longitudinal_Displacement.", device)
  )

  if (device == "pdf") open_pdf(path, 10.5, 7.0) else open_png(path, 2100, 1400, 200)

  ord <- c(
    "PRJNA691455",
    "PRJEB82425",
    "PRJNA851469",
    "PRJNA516701",
    "PRJNA578267",
    "PRJNA1166732",
    "PRJNA430161"
  )

  d <- W %>%
    mutate(
      ord = match(project, ord)
    ) %>%
    arrange(ord)

  labels <- c(
    "Core sepsis",
    "External ICU infection",
    "ICU background 1",
    "ICU background 2",
    "Non-sepsis surgical control",
    "Intervention support 1",
    "Intervention support 2"
  )

  y <- rev(seq_len(nrow(d)))

  par(mar = c(4.5, 8.0, 1.0, 2.0))

  xlim <- range(
    c(d$primary_dz, d$common_anchor_dz),
    na.rm = TRUE
  )
  xlim <- c(min(-0.75, xlim[1]-0.1), max(2.25, xlim[2]+0.1))

  plot(
    NA,
    xlim = xlim,
    ylim = c(0.5, nrow(d)+0.5),
    yaxt = "n",
    xlab = "Paired standardized change in Bray-Curtis displacement (dz)",
    ylab = "",
    bty = "n"
  )

  abline(v = 0, lty = 2)

  axis(2, at = y, labels = labels, las = 1, tick = FALSE)

  for (i in seq_len(nrow(d))) {
    segments(
      d$primary_dz[i], y[i],
      d$common_anchor_dz[i], y[i],
      lty = 3
    )
  }

  points(d$primary_dz, y, pch = 16, cex = 1.0)
  points(d$common_anchor_dz, y, pch = 1, cex = 1.15)

  # significance marker using frozen FDR
  sig <- ifelse(d$primary_fdr < 0.05, "*", "")
  text(
    rep(xlim[2]-0.02, nrow(d)),
    y,
    labels = sig,
    adj = c(1,0.5),
    cex = 1.2
  )

  legend(
    "bottomright",
    legend = c("Primary early-to-late", "Common-anchor sensitivity", "Primary FDR < 0.05"),
    pch = c(16,1,NA),
    lty = c(NA,NA,NA),
    text.width = strwidth("Common-anchor sensitivity"),
    bty = "n"
  )
  text(
    xlim[2]-0.02,
    0.55,
    "*",
    adj = c(1,0),
    cex = 1.2
  )

  dev.off()
}

make_fig2("pdf")
make_fig2("png")

# ------------------------------------------------------------
# FIGURE 3
# Core sepsis: shared ecological displacement + heterogeneous taxonomic routes
# ------------------------------------------------------------

get_u9 <- function(outcome) {
  U9 %>%
    filter(
      predictor == "mean_signed_trajectory_distance",
      outcome == !!outcome
    ) %>%
    slice(1)
}

r_eii <- get_u9("latest_EII")
r_sim <- get_u9("latest_Simpson_instability")

make_fig3 <- function(device = c("pdf","png")) {

  device <- match.arg(device)
  path <- file.path(
    FIGDIR,
    paste0("Figure_3_Core_Sepsis_Heterogeneous_Taxonomic_Routes.", device)
  )

  if (device == "pdf") open_pdf(path, 11.0, 8.0) else open_png(path, 2200, 1600, 200)

  par(mfrow = c(2,2), mar = c(4.3,4.4,2.2,1.0), oma = c(0,0,0.5,0))

  # A paired Day3-Day7
  rng <- range(c(U3$Day3, U3$Day7), na.rm = TRUE)
  plot(
    c(1,2), rng,
    type = "n",
    xaxt = "n",
    xlab = "",
    ylab = "Bray-Curtis displacement from baseline",
    main = "A. Paired Day 3 to Day 7 displacement",
    bty = "n"
  )
  axis(1, at = c(1,2), labels = c("Day 3","Day 7"))
  for (i in seq_len(nrow(U3))) {
    segments(1, U3$Day3[i], 2, U3$Day7[i])
    points(c(1,2), c(U3$Day3[i], U3$Day7[i]), pch = 16, cex = 0.75)
  }
  text(
    1.5, rng[1] + 0.03 * diff(rng),
    labels = paste0("n = ", nrow(U3), "; mean change = ",
                    safe_num(mean(U3$change_Day7_minus_Day3), 3)),
    cex = 0.78
  )

  # B heterogeneity distribution
  boxplot(
    U5$signed_Euclidean,
    horizontal = TRUE,
    outline = TRUE,
    xlab = "Pairwise signed-trajectory Euclidean distance",
    main = "B. Between-patient taxonomic-route heterogeneity",
    axes = TRUE
  )
  rug(U5$signed_Euclidean, side = 1)
  text(
    mean(U5$signed_Euclidean),
    1.22,
    paste0("mean = ", safe_num(mean(U5$signed_Euclidean), 3),
           "; median = ", safe_num(median(U5$signed_Euclidean), 3)),
    cex = 0.76
  )

  # C EII
  plot(
    U1$mean_signed_trajectory_distance,
    U1$latest_EII,
    pch = 16,
    xlab = "Mean taxonomic-trajectory distance",
    ylab = "Latest EII",
    main = "C. Trajectory heterogeneity vs ecological injury",
    bty = "n"
  )
  abline(lm(latest_EII ~ mean_signed_trajectory_distance, data = U1), lty = 2)
  legend(
    "topleft",
    legend = paste0(
      "Spearman rho = ", safe_num(r_eii$spearman_rho[1], 3),
      "; permutation p = ", safe_num(r_eii$permutation_p_10000[1], 3)
    ),
    bty = "n",
    cex = 0.78
  )

  # D Simpson instability
  plot(
    U1$mean_signed_trajectory_distance,
    U1$latest_Simpson_instability,
    pch = 16,
    xlab = "Mean taxonomic-trajectory distance",
    ylab = "Latest Simpson instability",
    main = "D. Trajectory heterogeneity vs diversity instability",
    bty = "n"
  )
  abline(
    lm(latest_Simpson_instability ~ mean_signed_trajectory_distance, data = U1),
    lty = 2
  )
  legend(
    "topleft",
    legend = paste0(
      "Spearman rho = ", safe_num(r_sim$spearman_rho[1], 3),
      "; permutation p = ", safe_num(r_sim$permutation_p_10000[1], 3)
    ),
    bty = "n",
    cex = 0.78
  )

  dev.off()
}

make_fig3("pdf")
make_fig3("png")

# ------------------------------------------------------------
# FIGURE 4
# External severe-trauma comparison
# ------------------------------------------------------------

parse_mid <- function(x) {
  if (is.na(x)) return(NA_real_)
  nums <- as.numeric(str_extract_all(x, "[0-9]+\\.?[0-9]*")[[1]])
  if (length(nums) == 0) return(NA_real_)
  mean(nums)
}

make_fig4 <- function(device = c("pdf","png")) {

  device <- match.arg(device)
  path <- file.path(
    FIGDIR,
    paste0("Figure_4_External_Validation_Trauma_Comparison.", device)
  )

  if (device == "pdf") open_pdf(path, 11.0, 8.0) else open_png(path, 2200, 1600, 200)

  par(mfrow = c(2,2), mar = c(4.6,4.4,2.2,1.0), oma = c(0,0,0.5,0))

  # A R2
  r2 <- sapply(X1$R2_or_range, parse_mid)
  bp <- barplot(
    r2,
    names.arg = c("Bray-Curtis","Hellinger","CLR/Aitchison"),
    las = 2,
    ylab = "PERMANOVA R2",
    main = "A. Sepsis vs trauma multivariate separation",
    ylim = c(0, max(r2, na.rm = TRUE)*1.35)
  )
  fdr_lab <- X1$PERMANOVA_FDR_or_range
  text(bp, r2, labels = paste0("FDR ", fdr_lab), pos = 3, cex = 0.70)

  # B dispersion FDR
  disp <- sapply(X1$dispersion_FDR_or_range, parse_mid)
  ydisp <- -log10(disp)
  bp2 <- barplot(
    ydisp,
    names.arg = c("Bray-Curtis","Hellinger","CLR/Aitchison"),
    las = 2,
    ylab = "-log10(dispersion FDR)",
    main = "B. Dispersion sensitivity",
    ylim = c(0, max(ydisp, na.rm = TRUE)*1.25)
  )
  abline(h = -log10(0.05), lty = 2)
  text(
    bp2, ydisp,
    labels = ifelse(X1$clean_location_evidence, "clean", "confounded"),
    pos = 3, cex = 0.72
  )

  # C alpha effect
  alpha <- X2 %>%
    filter(str_detect(evidence, "^Alpha_"))
  yy <- rev(seq_len(nrow(alpha)))
  plot(
    alpha$effect,
    yy,
    pch = 16,
    yaxt = "n",
    xlab = "Effect (Sepsis vs Trauma)",
    ylab = "",
    main = "C. Alpha-diversity depletion",
    xlim = c(min(alpha$effect)-0.1, 0.05),
    bty = "n"
  )
  abline(v = 0, lty = 2)
  axis(
    2, at = yy,
    labels = c("Observed genera","Shannon","Simpson"),
    las = 1, tick = FALSE
  )
  text(
    alpha$effect, yy,
    labels = paste0(" FDR=", formatC(alpha$FDR, format = "g", digits = 2)),
    pos = 4, cex = 0.72
  )

  # D taxonomic alignment (exploratory)
  rho_sc <- X3$value[X3$item == "Longitudinal-vs-Sepsis-Control Spearman rho"][1]
  rho_st <- X3$value[X3$item == "Longitudinal-vs-Sepsis-Trauma Spearman rho"][1]
  p_sc <- X3$value[X3$item == "Longitudinal-vs-Sepsis-Control p"][1]
  p_st <- X3$value[X3$item == "Longitudinal-vs-Sepsis-Trauma p"][1]

  vals <- c(rho_sc, rho_st)
  bp3 <- barplot(
    vals,
    names.arg = c("Sepsis-Control","Sepsis-Trauma"),
    ylab = "Spearman rho",
    main = "D. Exploratory taxonomic alignment",
    ylim = c(0, max(0.4, vals + 0.08))
  )
  abline(h = 0, lty = 2)
  text(
    bp3, vals,
    labels = c(
      paste0("p=", formatC(p_sc, format="g", digits=2)),
      paste0("p=", formatC(p_st, format="g", digits=2))
    ),
    pos = 3,
    cex = 0.74
  )

  dev.off()
}

make_fig4("pdf")
make_fig4("png")

# ------------------------------------------------------------
# FIGURE 5
# CRA002354 infection-source ecology: coefficient-centric final visual
# ------------------------------------------------------------

make_fig5 <- function(device = c("pdf","png")) {

  device <- match.arg(device)
  path <- file.path(
    FIGDIR,
    paste0("Figure_5_Infection_Source_Longitudinal_Ecology.", device)
  )

  if (device == "pdf") open_pdf(path, 11.0, 8.0) else open_png(path, 2200, 1600, 200)

  par(mfrow = c(2,2), mar = c(4.6,5.2,2.2,1.0), oma = c(0,0,0.5,0))

  # A primary interaction
  est <- B07$beta_interaction[1]
  lo <- B07$ci95_low[1]
  hi <- B07$ci95_high[1]

  plot(
    est, 1,
    xlim = range(c(lo,hi,0))*1.15,
    ylim = c(0.5,1.5),
    yaxt = "n",
    pch = 16,
    xlab = "Source-by-time interaction beta",
    ylab = "",
    main = "A. Primary Bray displacement interaction",
    bty = "n"
  )
  abline(v = 0, lty = 2)
  segments(lo,1,hi,1,lwd=2)
  points(est,1,pch=16)
  axis(2, at = 1, labels = "Pulmonary vs non-pulmonary", las = 1, tick = FALSE)
  legend(
    "bottomright",
    legend = paste0("p = ", safe_num(B07$p_LRT[1],3)),
    bty = "n",
    cex = 0.82
  )

  # B source-specific slopes
  yy <- c(2,1)
  plot(
    B08$slope_per_day,
    yy,
    xlim = range(c(B08$ci95_low,B08$ci95_high,0))*1.15,
    ylim = c(0.5,2.5),
    yaxt = "n",
    pch = 16,
    xlab = "Estimated Bray displacement slope per day",
    ylab = "",
    main = "B. Source-specific estimated slopes",
    bty = "n"
  )
  abline(v = 0, lty = 2)
  segments(B08$ci95_low, yy, B08$ci95_high, yy, lwd = 2)
  points(B08$slope_per_day, yy, pch = 16)
  axis(
    2,
    at = yy,
    labels = c("Pulmonary","Recorded non-pulmonary"),
    las = 1, tick = FALSE
  )
  text(
    B08$ci95_high, yy,
    labels = paste0(" p=", safe_num(B08$p_LRT,3)),
    pos = 4, cex = 0.72
  )

  # C sensitivity interaction betas
  d <- B20
  y <- rev(seq_len(nrow(d)))
  plot(
    d$beta, y,
    xlim = c(min(0,d$beta)-0.01, max(d$beta)+0.01),
    ylim = c(0.5,nrow(d)+0.5),
    yaxt = "n",
    pch = 16,
    xlab = "Source-by-time interaction beta",
    ylab = "",
    main = "C. Key sensitivity analyses",
    bty = "n"
  )
  abline(v=0,lty=2)
  axis(2, at = y, labels = gsub("_"," ",d$sensitivity), las = 1, tick = FALSE, cex.axis = 0.75)
  text(
    d$beta, y,
    labels = paste0(" p=", safe_num(d$p,3)),
    pos = 4, cex = 0.68
  )

  # D baseline beta robustness
  vals <- c(
    B04$primary_PERMANOVA_R2[1],
    B04$day3_restricted_PERMANOVA_R2[1]
  )
  bp <- barplot(
    vals,
    names.arg = c("Unrestricted baseline","Baseline <= ICU Day 3"),
    las = 2,
    ylab = "PERMANOVA R2",
    main = "D. Baseline source association robustness",
    ylim = c(0, max(vals)*1.8)
  )
  text(
    bp, vals,
    labels = c(
      paste0(
        "p=", safe_num(B04$primary_PERMANOVA_p[1],3),
        "\ndispersion p=", safe_num(B04$primary_dispersion_p[1],3)
      ),
      paste0(
        "p=", safe_num(B04$day3_restricted_PERMANOVA_p[1],3),
        "\ndispersion p=", safe_num(B04$day3_restricted_dispersion_p[1],3)
      )
    ),
    pos = 3,
    cex = 0.68
  )

  dev.off()
}

make_fig5("pdf")
make_fig5("png")

# ------------------------------------------------------------
# Figure captions
# ------------------------------------------------------------

captions <- c(
  "Figure 1. Study architecture and evidence hierarchy.",
  "Public gut-microbiome cohorts were organized into a longitudinal critical-illness framework, a nested core repeated-sepsis trajectory analysis, an external Control-Trauma-Sepsis ecological-state validation, and an infection-source longitudinal extension. The manuscript emphasizes shared ecological displacement with heterogeneous taxonomic routes rather than a universal taxonomic signature.",
  "",
  "Figure 2. Cross-cohort longitudinal ecological displacement.",
  "Primary early-to-late paired standardized Bray-Curtis displacement changes are shown together with common-anchor sensitivity estimates. Positive values indicate increasing displacement from baseline and negative values indicate re-convergence. Primary FDR-significant cohorts are marked with an asterisk. PRJNA578267 showed a non-significant recovery-direction primary contrast but significant re-convergence in common-anchor sensitivity analysis.",
  "",
  "Figure 3. Heterogeneous taxonomic routes accompany ecological displacement in core sepsis.",
  paste0(
    "Paired Day 3-to-Day 7 Bray-Curtis displacement increased on average in the repeated core-sepsis patients. Pairwise signed-trajectory distances demonstrated substantial between-patient taxonomic heterogeneity. Mean trajectory distance was positively associated with latest EII (Spearman rho ",
    safe_num(r_eii$spearman_rho[1],3),
    ", permutation p ",
    safe_num(r_eii$permutation_p_10000[1],3),
    ") and latest Simpson instability (rho ",
    safe_num(r_sim$spearman_rho[1],3),
    ", permutation p ",
    safe_num(r_sim$permutation_p_10000[1],3),
    "). These analyses are ecological/trajectory interpretations and do not establish stable taxonomic subtypes."
  ),
  "",
  "Figure 4. External ecological-state differentiation beyond severe trauma.",
  "Sepsis-vs-Trauma separation was significant across Bray-Curtis, Hellinger-Euclidean, and CLR/Aitchison distance families; only the CLR/Aitchison family showed clean location evidence without significant dispersion heterogeneity. Sepsis also showed lower alpha diversity than trauma and greater distance from the healthy-control centroid. Taxonomic alignment analyses are exploratory and do not establish a reproducible sepsis-specific genus signature.",
  "",
  "Figure 5. Infection source does not clearly modify longitudinal ecological displacement in CRA002354.",
  paste0(
    "The primary source-by-personal-time Bray interaction was ",
    safe_num(B07$beta_interaction[1],3),
    " (95% CI ",
    safe_num(B07$ci95_low[1],3),
    " to ",
    safe_num(B07$ci95_high[1],3),
    "; p=",
    safe_num(B07$p_LRT[1],3),
    "). Estimated source-specific slopes were positive in both pulmonary and recorded non-pulmonary groups, but neither reached conventional significance. Key sensitivity analyses preserved the positive interaction direction without statistical significance. The unrestricted baseline PERMANOVA signal was dispersion-confounded and disappeared when baseline was restricted to ICU Day 3 or earlier."
  )
)

writeLines(
  captions,
  file.path(CAPDIR, "DRAFT_MAIN_FIGURE_CAPTIONS.txt")
)

# ------------------------------------------------------------
# Table placement plan
# ------------------------------------------------------------

table_plan <- tibble(
  current_item = c(
    "Step95C Main_Table_1_Evidence_Architecture",
    "Step95C Main_Table_2_Cross_Cohort_Longitudinal_Results",
    "Step95C Main_Table_3_Compact_Supporting_Result_Summary"
  ),
  recommended_manuscript_location = c(
    "INTERNAL_OR_SUPPLEMENTARY",
    "MAIN_OR_SUPPLEMENTARY_DEPENDING_JOURNAL_LIMIT",
    "INTERNAL_OR_SUPPLEMENTARY"
  ),
  reason = c(
    "Evidence hierarchy is useful for provenance but is not a conventional manuscript Table 1. A manuscript Table 1 should instead summarize cohort/study characteristics.",
    "This directly reports frozen longitudinal effect estimates and is manuscript-relevant.",
    "This is a writing aid that mixes main and exploratory findings; it should not appear as a main manuscript table in its current form."
  ),
  next_action = c(
    "Build a conventional cohort-characteristics Table 1 from frozen metadata.",
    "Format into journal-ready result table after figure QC.",
    "Use as Results drafting checklist; move detailed exploratory items to Supplementary."
  )
)

write_csv(
  table_plan,
  file.path(TABDIR, "TABLE_PLACEMENT_PLAN.csv")
)

# ------------------------------------------------------------
# Visual QC report
# ------------------------------------------------------------

qc_notes <- c(
  "STEP95D VISUAL QC DECISIONS",
  "",
  "Figure 1:",
  "- Rebuilt because the Step95C version was an internal evidence map and treated core sepsis as a parallel rather than nested longitudinal branch.",
  "- Internal Step IDs were removed from the manuscript-facing visual.",
  "",
  "Figure 2:",
  "- Rebuilt from the frozen Step93W evidence matrix.",
  "- Retains primary and common-anchor estimates while removing internal STRICT/MOSTLY/LIMITED labels from the visual.",
  "",
  "Figure 3:",
  "- Rebuilt because the frozen diagnostic figure had overlapping labels and crowded panels.",
  "- New panels directly support the manuscript interpretation: paired displacement, pairwise taxonomic heterogeneity, EII association, and Simpson-instability association.",
  "",
  "Figure 4:",
  "- Rebuilt because the frozen diagnostic figure had overlapping titles/axis labels.",
  "- CLR/Aitchison is explicitly distinguished as the clean multivariate location evidence; taxonomic alignment remains exploratory.",
  "",
  "Figure 5:",
  "- Rebuilt from coefficient/sensitivity registries rather than raw QC scatterplots.",
  "- The final figure now directly displays the null primary interaction, source-specific slopes, sensitivity consistency, and baseline confounding/robustness.",
  "",
  "Table architecture:",
  "- Step95C evidence-architecture table should not be the manuscript Table 1.",
  "- A conventional cohort-characteristics Table 1 should be constructed next.",
  "- The cross-cohort longitudinal results table is manuscript-relevant.",
  "",
  "No new inferential analysis was performed."
)

writeLines(
  qc_notes,
  file.path(CAPDIR, "STEP95D_VISUAL_QC_REPORT.txt")
)

# ------------------------------------------------------------
# Final readiness
# ------------------------------------------------------------

figure_files <- file.path(
  FIGDIR,
  c(
    "Figure_1_Study_Architecture_Publication.pdf",
    "Figure_2_Cross_Cohort_Longitudinal_Displacement.pdf",
    "Figure_3_Core_Sepsis_Heterogeneous_Taxonomic_Routes.pdf",
    "Figure_4_External_Validation_Trauma_Comparison.pdf",
    "Figure_5_Infection_Source_Longitudinal_Ecology.pdf",
    "Figure_1_Study_Architecture_Publication.png",
    "Figure_2_Cross_Cohort_Longitudinal_Displacement.png",
    "Figure_3_Core_Sepsis_Heterogeneous_Taxonomic_Routes.png",
    "Figure_4_External_Validation_Trauma_Comparison.png",
    "Figure_5_Infection_Source_Longitudinal_Ecology.png"
  )
)

ready <- all(file.exists(figure_files))

write_csv(
  tibble(
    final_figure_files_created = ready,
    no_new_inference = TRUE,
    ready_for_panel_visual_review = ready,
    ready_for_conventional_table1_build = ready
  ),
  file.path(OUT, "STEP95D_READINESS.csv")
)

if (!ready) stop("One or more publication-ready figure files were not created.")

writeLines(
  c(
    paste0("Completed: ", Sys.time()),
    "Final Figure 1-5 PDF+PNG created: TRUE",
    "No new inferential analysis: TRUE",
    "Ready for panel visual review: TRUE",
    "STEP95D COMPLETE"
  ),
  file.path(OUT, "_STEP95D_COMPLETE.ok")
)

cat("STEP95D COMPLETE\n")
