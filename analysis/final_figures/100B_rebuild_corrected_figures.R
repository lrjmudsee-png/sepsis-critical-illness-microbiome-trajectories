# ============================================================
# Step100B: rebuild final manuscript figures from Step100A source pack
#
# No inferential analysis is performed. Scientific values, sample sizes,
# effects, p values, confidence intervals and feature counts are read from
# the Step100A machine-readable CSV files.
#
# Usage
#   Rscript analysis/final_figures/100B_rebuild_corrected_figures.R <source_csv_dir> <new_output_dir>
# ============================================================

options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(ragg)
})

args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 2) stop("Usage: Rscript 100B_rebuild_corrected_figures.R <source_csv_dir> <new_output_dir>")
source_dir <- normalizePath(args[1], winslash = "/", mustWork = TRUE)
out_dir <- normalizePath(args[2], winslash = "/", mustWork = FALSE)
if (file.exists(out_dir) || dir.exists(out_dir)) stop("Output already exists; choose a new output directory.")
source(file.path(dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))), "figure_display_overrides.R"), encoding = "UTF-8")
main_dir <- file.path(out_dir, "02_MAIN_FIGURES")
supp_dir <- file.path(out_dir, "03_SUPPLEMENTARY_FIGURES")
required_source <- c(
  "FIGURE1_architecture.csv",
  "FIGURE2_longitudinal_displacement.csv",
  "FIGURE3a_day3_day7_paired.csv",
  "FIGURE3b_pairwise_trajectory.csv",
  "FIGURE3c_correlations.csv",
  "FIGURE3d_cross_cohort_taxonomic_concordance.csv",
  "FIGURE3_summary.csv",
  "FIGURE4_group_counts.csv",
  "FIGURE4_centroid_comparison.csv",
  "FIGURE4_aitchison_pseudocount_sensitivity.csv",
  "FIGURE4_taxonomic_alignment.csv",
  "FIGURE4a_sample_distances.csv",
  "FIGURE5_primary_interaction.csv",
  "FIGURE5_source_specific_slopes.csv",
  "FIGURE5_sensitivity.csv",
  "FIGURE5_baseline_robustness.csv",
  "SUPPLEMENTARY_FIGURE_S1_taxonomic_concordance_sensitivity.csv",
  "SUPPLEMENTARY_FIGURE_S2_pair_counts.csv",
  "SUPPLEMENTARY_FIGURE_S3_bray_meta.csv",
  "SUPPLEMENTARY_FIGURE_S4_aitchison_meta.csv",
  "SUPPLEMENTARY_FIGURES_S3_S4_pooled.csv",
  "SUPPLEMENTARY_FIGURE_S5_family_balance.csv",
  "SUPPLEMENTARY_FIGURE_S6_healthy_reference_paired.csv",
  "SUPPLEMENTARY_FIGURE_S6_summary.csv"
)
missing_source <- required_source[!file.exists(file.path(source_dir, required_source))]
if (length(missing_source)) stop("Step100A source pack is incomplete: ", paste(missing_source, collapse = ", "))
dir.create(main_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(supp_dir, recursive = TRUE, showWarnings = FALSE)

read_pack <- function(name) read_csv(file.path(source_dir, name), show_col_types = FALSE, progress = FALSE)

blue <- "#0072B2"
orange <- "#D55E00"
green <- "#009E73"
purple <- "#CC79A7"
grey <- "#666666"
light_grey <- "#D9D9D9"

fmt <- function(x, digits = 3) formatC(x, format = "f", digits = digits)
fmt_p <- function(x) {
  ifelse(x < 0.001, format(x, scientific = TRUE, digits = 2), formatC(x, format = "f", digits = 3))
}

panel_label <- function(label) {
  mtext(label, side = 3, adj = 0, line = 0.35, font = 2, cex = 1.15)
}

open_device <- function(path, type, width, height) {
  if (type == "pdf") {
    cairo_pdf(path, width = width, height = height, family = "Arial", onefile = FALSE)
  } else if (type == "png") {
    agg_png(path, width = width, height = height, units = "in", res = 300, background = "white")
  } else if (type == "tiff") {
    agg_tiff(path, width = width, height = height, units = "in", res = 600,
             compression = "lzw", background = "white")
  } else {
    stop("Unknown device type: ", type)
  }
}

render_set <- function(stem, directory, width, height, plot_fun) {
  width <- 7.25
  height <- switch(stem,
    Figure_1_Study_Architecture = 5.5,
    Figure_2_Cross_Cohort_Longitudinal_Displacement = 4.6,
    Figure_3_Core_Sepsis_and_Taxonomic_Concordance = 6.8,
    Figure_4_External_Validation_Trauma_Comparison = 6.0,
    Figure_5_Infection_Source_Longitudinal_Ecology = 7.2,
    Supplementary_Figure_S1_Taxonomic_Concordance_Sensitivity = 6.6,
    Supplementary_Figure_S2_Pair_Counts = 4.4,
    Supplementary_Figure_S3_Bray_Meta_Analysis = 3.5,
    Supplementary_Figure_S4_CZM_Aitchison_Meta_Analysis = 3.5,
    Supplementary_Figure_S5_Family_Balance = 4.5,
    Supplementary_Figure_S6_Healthy_Reference = 4.4)
  revised <- corrected_plot(stem)
  if (is.null(revised)) stop("No reviewed final plot implementation for: ", stem)
  rows <- list()
  for (type in c("pdf", "png", "tiff")) {
    path <- file.path(directory, paste0(stem, ".", type))
    open_device(path, type, width, height)
    tryCatch(if (is.null(revised)) plot_fun() else print(revised), finally = dev.off())
    rows[[type]] <- tibble(
      figure = stem,
      format = toupper(type),
      path = normalizePath(path, winslash = "/", mustWork = TRUE),
      bytes = file.info(path)$size,
      width_inches = width,
      height_inches = height,
      dpi = ifelse(type == "pdf", NA_integer_, ifelse(type == "png", 300L, 600L)),
      scientific_values_source = "Step100A CSV source pack"
    )
  }
  bind_rows(rows)
}

export_manifest <- tibble()

# ------------------------------------------------------------------
# Figure 1
# ------------------------------------------------------------------

f1 <- read_pack("FIGURE1_architecture.csv")
f1_val <- function(branch, field) f1[[field]][f1$branch == branch][1]

plot_figure1 <- function() {
  par(mar = c(0.4, 0.4, 0.4, 0.4), family = "sans")
  plot(NA, xlim = c(0, 12), ylim = c(0, 9), axes = FALSE, xlab = "", ylab = "", xaxs = "i", yaxs = "i")

  rect(3.5, 8.0, 8.5, 8.7, lwd = 1.4, border = grey)
  text(6, 8.35, "Public gut-microbiome cohorts", font = 2, cex = 1.05)
  text(6, 8.08, paste0("n = ", sum(f1$cohort_count[f1$branch %in% c("DISCOVERY", "INDEPENDENT_LONGITUDINAL_VALIDATION", "EXTERNAL_CROSS_SECTIONAL_COMPARISON", "INFECTION_SOURCE_EXTENSION")]), " cohorts"), cex = 0.82)

  x0 <- c(0.25, 3.2, 6.15, 9.1)
  x1 <- c(2.9, 5.85, 8.8, 11.75)
  fill <- c("#E8F2F8", "#E8F5EF", "#FCEFE8", "#F2ECF6")
  borders <- c(blue, green, orange, purple)
  branches <- c("DISCOVERY", "INDEPENDENT_LONGITUDINAL_VALIDATION", "EXTERNAL_CROSS_SECTIONAL_COMPARISON", "INFECTION_SOURCE_EXTENSION")
  titles <- c("Discovery evidence", "Independent longitudinal\nvalidation", "External cross-sectional\necological-state comparison", "Infection-source\nlongitudinal extension")
  details <- c(
    f1_val("DISCOVERY", "detail"),
    f1_val("INDEPENDENT_LONGITUDINAL_VALIDATION", "detail"),
    "PRJNA1010969\nControl / Trauma / Sepsis\ncross-sectional comparison",
    "CRA002354\nPulmonary vs recorded\nnon-pulmonary"
  )
  for (i in seq_along(branches)) {
    rect(x0[i], 5.5, x1[i], 7.1, col = fill[i], border = borders[i], lwd = 1.4)
    text(mean(c(x0[i], x1[i])), 6.75, titles[i], font = 2, cex = 0.83)
    text(mean(c(x0[i], x1[i])), 6.12, details[i], cex = 0.72)
    arrows(6, 8.0, mean(c(x0[i], x1[i])), 7.1, length = 0.07, col = grey)
  }

  rect(0.5, 3.1, 2.65, 4.55, border = blue, lwd = 1.2)
  text(1.575, 4.18, "Longitudinal discovery", font = 2, cex = 0.79)
  text(1.575, 3.70, paste0("n = ", f1_val("DISCOVERY_LONGITUDINAL", "cohort_count"), " cohorts\nwithin-patient change"), cex = 0.72)
  rect(3.0, 3.1, 5.15, 4.55, border = blue, lwd = 1.2)
  text(4.075, 4.18, "Static support", font = 2, cex = 0.79)
  text(4.075, 3.70, paste0("n = ", f1_val("DISCOVERY_STATIC", "cohort_count"), " cohort\nPRJNA978257"), cex = 0.72)
  arrows(1.575, 5.5, 1.575, 4.55, length = 0.07, col = grey)
  arrows(1.575, 5.5, 4.075, 4.55, length = 0.07, col = grey)

  rect(1.3, 0.7, 10.7, 2.15, border = grey, lwd = 1.4)
  text(6, 1.62, "Longitudinal ecological displacement was observed across multiple cohorts", font = 2, cex = 0.93)
  text(6, 1.18, "despite heterogeneous taxonomic trajectories.", cex = 0.88)
  arrows(1.575, 3.1, 3.4, 2.15, length = 0.07, col = grey)
  arrows(4.075, 3.1, 4.7, 2.15, length = 0.07, col = grey)
  arrows(4.525, 5.5, 5.5, 2.15, length = 0.07, col = grey)
  arrows(7.475, 5.5, 6.9, 2.15, length = 0.07, col = grey)
  arrows(10.425, 5.5, 8.6, 2.15, length = 0.07, col = grey)
}

export_manifest <- bind_rows(export_manifest, render_set("Figure_1_Study_Architecture", main_dir, 12, 9, plot_figure1))

# ------------------------------------------------------------------
# Figure 2
# ------------------------------------------------------------------

f2 <- read_pack("FIGURE2_longitudinal_displacement.csv")
f2_order <- c("PRJNA691455", "PRJEB82425", "PRJNA851469", "PRJNA516701", "PRJNA578267", "PRJNA1166732", "PRJNA430161")
f2 <- f2 %>% mutate(display_order = match(project, f2_order)) %>% arrange(display_order)

plot_figure2 <- function() {
  y <- rev(seq_len(nrow(f2)))
  vals <- c(f2$primary_effect_dz, f2$common_anchor_effect_dz)
  xr <- range(vals, na.rm = TRUE)
  pad <- diff(xr) * 0.13
  par(mar = c(5.2, 12.6, 1.1, 1.2), family = "sans")
  plot(NA, xlim = c(xr[1] - pad, xr[2] + pad), ylim = c(0.4, nrow(f2) + 0.7), yaxt = "n",
       xlab = expression("Paired standardized change in Bray-Curtis displacement (" * d[z] * ")"), ylab = "", bty = "n")
  abline(v = 0, lty = 2, col = grey)
  labels <- paste0(f2$project, "  n = ", f2$primary_n)
  axis(2, at = y, labels = labels, las = 1, tick = FALSE, cex.axis = 0.80)
  segments(f2$primary_effect_dz, y, f2$common_anchor_effect_dz, y, lty = 3, col = grey)
  points(f2$primary_effect_dz, y, pch = 16, col = blue, cex = 1.1)
  points(f2$common_anchor_effect_dz, y, pch = 21, bg = "white", col = orange, cex = 1.15)
  sig <- which(f2$primary_significant)
  if (length(sig)) text(f2$primary_effect_dz[sig], y[sig] + 0.23, "*", col = blue, cex = 1.05)
  legend("bottomright", c("Primary early-to-late", "Common-anchor sensitivity", "* Primary BH-FDR < 0.05"),
         pch = c(16, 21, NA), pt.bg = c(NA, "white", NA), col = c(blue, orange, blue), bty = "n", cex = 0.82)
  mtext("Lines connect the primary and common-anchor sensitivity estimates from the same cohort and do not represent confidence intervals.", side = 1, line = 3.6, cex = 0.66, col = grey)
}

export_manifest <- bind_rows(export_manifest, render_set("Figure_2_Cross_Cohort_Longitudinal_Displacement", main_dir, 11.5, 7.5, plot_figure2))

# ------------------------------------------------------------------
# Figure 3
# ------------------------------------------------------------------

f3a <- read_pack("FIGURE3a_day3_day7_paired.csv")
f3b <- read_pack("FIGURE3b_pairwise_trajectory.csv")
f3c <- read_pack("FIGURE3c_correlations.csv")
f3d <- read_pack("FIGURE3d_cross_cohort_taxonomic_concordance.csv")
f3s <- read_pack("FIGURE3_summary.csv")

plot_figure3 <- function() {
  par(mfrow = c(2, 2), mar = c(4.7, 8.2, 3.3, 1.2), oma = c(0.3, 0.2, 0.4, 0.2), family = "sans")

  plot(c(1, 2), range(c(f3a$Day3, f3a$Day7)), type = "n", xaxt = "n",
       xlab = "", ylab = "Bray-Curtis displacement", main = "Within-patient ecological displacement")
  axis(1, at = c(1, 2), labels = c("Day 3", "Day 7"))
  for (i in seq_len(nrow(f3a))) lines(c(1, 2), c(f3a$Day3[i], f3a$Day7[i]), col = adjustcolor(grey, 0.55), lwd = 1)
  points(rep(1, nrow(f3a)), f3a$Day3, pch = 16, col = blue)
  points(rep(2, nrow(f3a)), f3a$Day7, pch = 16, col = orange)
  mtext(paste0("paired patients n = ", nrow(f3a), "; mean change = ", fmt(f3s$day3_day7_mean_change[1], 3)), side = 3, line = 0.25, cex = 0.72, col = grey)
  panel_label("a")

  set.seed(100)
  boxplot(f3b$signed_Euclidean, horizontal = TRUE, col = "#E8F2F8", border = blue,
          xlab = "Pairwise signed-trajectory Euclidean distance", yaxt = "n", main = "Taxonomic trajectory heterogeneity")
  points(f3b$signed_Euclidean, jitter(rep(1, nrow(f3b)), amount = 0.08), pch = 16, cex = 0.6, col = adjustcolor(grey, 0.65))
  mtext(paste0(f3s$pairwise_comparisons[1], " patient-pair comparisons among ", f3s$pairwise_source_patients[1], " patients"), side = 3, line = 0.25, cex = 0.72, col = grey)
  panel_label("b")

  c_order <- c("latest_EII", "latest_Simpson_instability", "latest_Bray")
  f3c2 <- f3c %>% mutate(ord = match(outcome, c_order)) %>% arrange(ord)
  yy <- rev(seq_len(nrow(f3c2)))
  c_labels <- c("Exploratory EII", "Simpson instability", "Latest Bray displacement")
  plot(f3c2$spearman_rho, yy, xlim = c(min(-0.05, min(f3c2$spearman_rho) - 0.08), max(f3c2$spearman_rho) + 0.26),
       ylim = c(0.5, nrow(f3c2) + 0.5), yaxt = "n", xlab = expression("Spearman " * rho), ylab = "", pch = 16, col = blue,
       main = "Trajectory heterogeneity and ecology", bty = "n")
  abline(v = 0, lty = 2, col = grey)
  axis(2, at = yy, labels = c_labels, las = 1, tick = FALSE, cex.axis = 0.72)
  text(f3c2$spearman_rho, yy, labels = paste0("  permutation p = ", vapply(f3c2$permutation_p_10000, fmt_p, character(1))), pos = 4, cex = 0.69)
  mtext(paste0("patient-level n = ", paste(unique(f3c2$n), collapse = ";")), side = 3, line = 0.25, cex = 0.72, col = grey)
  panel_label("c")

  ranks <- c("Genus", "Family")
  pairs <- unique(f3d$pair_label)
  ybase <- rev(seq_along(pairs))
  xlim <- range(c(0, f3d$spearman_rho), na.rm = TRUE) + c(-0.05, 0.08)
  plot(NA, xlim = xlim, ylim = c(0.5, length(pairs) + 0.5), yaxt = "n", xlab = expression("Spearman " * rho), ylab = "",
       main = "Cross-cohort taxonomic concordance", bty = "n")
  abline(v = 0, lty = 2, col = grey)
  axis(2, at = ybase, labels = pairs, las = 1, tick = FALSE, cex.axis = 0.68)
  for (j in seq_along(ranks)) {
    dd <- f3d %>% filter(rank == ranks[j])
    yy2 <- ybase[match(dd$pair_label, pairs)] + ifelse(j == 1, 0.09, -0.09)
    points(dd$spearman_rho, yy2, pch = ifelse(j == 1, 16, 17), col = ifelse(j == 1, blue, orange), cex = 1)
    text(dd$spearman_rho, yy2, labels = paste0(" n = ", dd$n_shared), pos = 4, cex = 0.6)
  }
  legend("bottomright", ranks, pch = c(16, 17), col = c(blue, orange), bty = "n", cex = 0.78)
  panel_label("d")
}

export_manifest <- bind_rows(export_manifest, render_set("Figure_3_Core_Sepsis_and_Taxonomic_Concordance", main_dir, 14.0, 9.5, plot_figure3))

# ------------------------------------------------------------------
# Figure 4
# ------------------------------------------------------------------

f4n <- read_pack("FIGURE4_group_counts.csv")
f4a <- read_pack("FIGURE4_centroid_comparison.csv")
f4b <- read_pack("FIGURE4_aitchison_pseudocount_sensitivity.csv")
f4c <- read_pack("FIGURE4_taxonomic_alignment.csv")

plot_figure4 <- function() {
  par(mfrow = c(1, 3), mar = c(5.2, 5.4, 3.6, 1.0), oma = c(0.3, 0.2, 0.3, 0.2), family = "sans")

  med <- c(f4a$median_B[1], f4a$median_A[1])
  names(med) <- c(f4a$group_B[1], f4a$group_A[1])
  cols <- c(grey, orange)
  bp <- barplot(med, col = cols, border = NA, ylim = c(0, max(med) * 1.28), ylab = "Median distance from control centroid",
                main = "Ecological-state differentiation")
  n_lookup <- setNames(f4n$n, f4n$Group)
  text(bp, med, labels = paste0("n = ", n_lookup[names(med)]), pos = 3, cex = 0.76)
  mtext(paste0("Control reference n = ", n_lookup[["Control"]], "; difference = ", fmt(f4a$median_difference_A_minus_B[1], 3), "; BH-FDR = ", fmt_p(f4a$FDR[1])), side = 1, line = 3.7, cex = 0.66, col = grey)
  panel_label("a")

  pc_labels <- format(f4b$pseudocount, scientific = TRUE, digits = 1)
  plot(seq_len(nrow(f4b)), f4b$PERMANOVA_R2, xlim = c(0.45, nrow(f4b) + 0.55), xaxt = "n", pch = 16, col = blue, cex = 1.2,
       ylim = range(f4b$PERMANOVA_R2) + c(-0.004, 0.004), xlab = "Pseudocount", ylab = expression("PERMANOVA " * R^2),
       main = "Aitchison pseudocount sensitivity")
  axis(1, at = seq_len(nrow(f4b)), labels = pc_labels)
  text(seq_len(nrow(f4b)), f4b$PERMANOVA_R2, labels = paste0("FDR = ", vapply(f4b$PERMANOVA_FDR, fmt_p, character(1)), "\ndispersion FDR = ", vapply(f4b$dispersion_FDR, fmt_p, character(1))), pos = 3, cex = 0.63)
  mtext(paste0(nrow(f4b), " pseudocount variants; range across pseudocount sensitivity analyses"), side = 1, line = 3.7, cex = 0.66, col = grey)
  panel_label("b")

  yy <- rev(seq_len(nrow(f4c)))
  labels <- c("Longitudinal vs\nSepsis-Control", "Longitudinal vs\nSepsis-Trauma")
  plot(f4c$spearman_rho, yy, xlim = c(min(0, f4c$spearman_rho) - 0.04, max(f4c$spearman_rho) + 0.13),
       ylim = c(0.5, nrow(f4c) + 0.5), yaxt = "n", pch = 16, col = green, xlab = expression("Spearman " * rho), ylab = "",
       main = "Exploratory taxonomic alignment", bty = "n")
  abline(v = 0, lty = 2, col = grey)
  axis(2, at = yy, labels = labels, las = 1, tick = FALSE, cex.axis = 0.72)
  text(f4c$spearman_rho, yy, labels = paste0(" n = ", f4c$overlapping_genera, "; p = ", vapply(f4c$p_value, fmt_p, character(1))), pos = 4, cex = 0.67)
  panel_label("c")
}

export_manifest <- bind_rows(export_manifest, render_set("Figure_4_External_Validation_Trauma_Comparison", main_dir, 12.5, 5.8, plot_figure4))

# ------------------------------------------------------------------
# Figure 5
# ------------------------------------------------------------------

f5a <- read_pack("FIGURE5_primary_interaction.csv")
f5b <- read_pack("FIGURE5_source_specific_slopes.csv")
f5c <- read_pack("FIGURE5_sensitivity.csv")
f5d <- read_pack("FIGURE5_baseline_robustness.csv")

plot_ci <- function(est, lo, hi, y, xlab, labels, main, col = blue, p_labels = NULL) {
  xr <- range(c(lo, hi, 0), na.rm = TRUE)
  pad <- diff(xr) * 0.12
  plot(est, y, xlim = xr + c(-pad, pad), ylim = c(0.5, max(y) + 0.5), yaxt = "n", pch = 16, col = col,
       xlab = xlab, ylab = "", main = main, bty = "n")
  abline(v = 0, lty = 2, col = grey)
  segments(lo, y, hi, y, col = col, lwd = 2)
  axis(2, at = y, labels = labels, las = 1, tick = FALSE, cex.axis = 0.70)
  if (!is.null(p_labels)) text(hi, y, labels = p_labels, pos = 4, cex = 0.66)
}

plot_figure5 <- function() {
  par(mfrow = c(2, 2), mar = c(5.0, 9.4, 3.6, 1.2), oma = c(0.3, 0.2, 0.3, 0.2), family = "sans")

  plot_ci(f5a$effect, f5a$ci_low, f5a$ci_high, 1, "Source-by-time interaction beta",
          "Pulmonary vs recorded\nnon-pulmonary", "Primary Bray interaction", orange,
          paste0(" LRT p = ", fmt_p(f5a$p_value)))
  mtext("Horizontal interval: 95% CI", side = 1, line = 3.7, cex = 0.68, col = grey)
  panel_label("a")

  yy <- rev(seq_len(nrow(f5b)))
  labels <- paste0(f5b$group_label, "\n", f5b$npatients, " patients / ", f5b$nobs, " observations")
  plot_ci(f5b$slope_per_day, f5b$ci95_low, f5b$ci95_high, yy, "Bray displacement slope per day",
          labels, "Source-specific slopes", blue, paste0(" LRT p = ", vapply(f5b$p_LRT, fmt_p, character(1))))
  mtext("Horizontal intervals: 95% CIs", side = 1, line = 3.7, cex = 0.68, col = grey)
  panel_label("b")

  yy3 <- rev(seq_len(nrow(f5c)))
  sens_labels <- tools::toTitleCase(tolower(gsub("_", " ", f5c$sensitivity)))
  plot_ci(f5c$beta, f5c$ci_low, f5c$ci_high, yy3, "Source-by-time interaction beta",
          sens_labels, "Sensitivity analyses", green, paste0(" LRT p = ", vapply(f5c$p_value, fmt_p, character(1))))
  mtext("Horizontal intervals: 95% CIs", side = 1, line = 3.7, cex = 0.68, col = grey)
  panel_label("c")

  vals <- c(f5d$primary_PERMANOVA_R2[1], f5d$day3_restricted_PERMANOVA_R2[1])
  bp <- barplot(vals, names.arg = expression("Unrestricted", "" <= "ICU Day 3"), col = c(orange, blue), border = NA,
                ylim = c(0, max(vals) * 1.9), ylab = expression("PERMANOVA " * R^2), main = "Baseline beta-diversity robustness")
  labs <- c(
    paste0("PERMANOVA p = ", fmt_p(f5d$primary_PERMANOVA_p[1]), "\ndispersion p = ", fmt_p(f5d$primary_dispersion_p[1])),
    paste0("PERMANOVA p = ", fmt_p(f5d$day3_restricted_PERMANOVA_p[1]), "\ndispersion p = ", fmt_p(f5d$day3_restricted_dispersion_p[1]))
  )
  text(bp, vals, labels = labs, pos = 3, cex = 0.68)
  panel_label("d")
}

export_manifest <- bind_rows(export_manifest, render_set("Figure_5_Infection_Source_Longitudinal_Ecology", main_dir, 14.0, 9.5, plot_figure5))

# ------------------------------------------------------------------
# Updated Supplementary Figures S1, S5 and S6
# ------------------------------------------------------------------

s1 <- read_pack("SUPPLEMENTARY_FIGURE_S1_taxonomic_concordance_sensitivity.csv")
plot_s1 <- function() {
  par(mfrow = c(2, 1), mar = c(4.4, 11.4, 3.4, 1.2), oma = c(0.2, 0.2, 0.3, 0.2), family = "sans")
  roles <- c("Primary", "Common-anchor sensitivity")
  for (k in seq_along(roles)) {
    dd <- s1 %>% filter(analysis_role == roles[k])
    pairs <- unique(dd$pair_label)
    yy <- rev(seq_along(pairs))
    plot(NA, xlim = range(c(0, dd$spearman_rho)) + c(-0.05, 0.07), ylim = c(0.5, length(pairs) + 0.5), yaxt = "n",
         xlab = "Spearman correlation of longitudinal CLR effect vectors", ylab = "", main = roles[k], bty = "n")
    abline(v = 0, lty = 2, col = grey)
    axis(2, at = yy, labels = pairs, las = 1, tick = FALSE, cex.axis = 0.72)
    for (rank_name in c("Genus", "Family")) {
      d2 <- dd %>% filter(rank == rank_name)
      y2 <- yy[match(d2$pair_label, pairs)] + ifelse(rank_name == "Genus", 0.09, -0.09)
      points(d2$spearman_rho, y2, pch = ifelse(rank_name == "Genus", 16, 17), col = ifelse(rank_name == "Genus", blue, orange))
      text(d2$spearman_rho, y2, labels = paste0(" ", fmt(d2$spearman_rho, 2)), pos = 4, cex = 0.68)
    }
    legend("bottomright", c("Genus", "Family"), pch = c(16, 17), col = c(blue, orange), bty = "n", cex = 0.76)
    panel_label(letters[k])
  }
}
export_manifest <- bind_rows(export_manifest, render_set("Supplementary_Figure_S1_Taxonomic_Concordance_Sensitivity", supp_dir, 12, 8.5, plot_s1))

# Supplementary Figure S2: counts only; no inferential analysis.
s2 <- read_pack("SUPPLEMENTARY_FIGURE_S2_pair_counts.csv")
plot_s2 <- function() {
  roles <- c("CORE_SEPSIS_LONGITUDINAL" = "Core sepsis", "ICU_INFECTION_LONGITUDINAL_EXTERNAL" = "External ICU infection",
             "ICU_BACKGROUND_LONGITUDINAL" = "ICU background", "NONSEPSIS_LONGITUDINAL_CONTROL" = "Non-sepsis surgical control",
             "INTERVENTION_LONGITUDINAL_SUPPORT" = "Intervention support")
  labels <- paste0(unname(roles[s2$cohort_role]), "\n", s2$project)
  y <- rev(seq_len(nrow(s2)))
  par(mar = c(5.2, 11.5, 3.0, 1.2), family = "sans")
  plot(s2$primary_n, y, xlim = c(0, max(c(s2$primary_n, s2$common_anchor_n)) * 1.12), ylim = c(0.5, nrow(s2) + 0.5),
       yaxt = "n", xlab = "Number of paired patients", ylab = "", main = "Paired patient counts by analysis", pch = 16, col = blue, bty = "n")
  axis(2, at = y, labels = labels, las = 1, tick = FALSE, cex.axis = 0.74)
  segments(s2$primary_n, y, s2$common_anchor_n, y, col = light_grey, lwd = 2)
  points(s2$common_anchor_n, y, pch = 21, bg = "white", col = orange, cex = 1.05)
  legend("bottomright", c("Primary", "Common-anchor sensitivity"), pch = c(16, 21), pt.bg = c(NA, "white"), col = c(blue, orange), bty = "n", cex = 0.78)
}
export_manifest <- bind_rows(export_manifest, render_set("Supplementary_Figure_S2_Pair_Counts", supp_dir, 11.5, 7.2, plot_s2))

# Supplementary Figures S3-S4: forest plots from frozen source-pack estimates.
s34 <- read_pack("SUPPLEMENTARY_FIGURES_S3_S4_pooled.csv")
plot_meta <- function(dat, pooled, title) {
  y <- rev(seq_len(nrow(dat))) + 1
  xr <- range(c(dat$gz_ci_low, dat$gz_ci_high, pooled$ci_low, pooled$ci_high, 0), na.rm = TRUE)
  pad <- diff(xr) * 0.12
  par(mar = c(5.1, 5.5, 3.2, 1.4), family = "sans")
  plot(dat$hedges_gz, y, xlim = xr + c(-pad, pad), ylim = c(0.5, max(y) + 0.7), yaxt = "n",
       xlab = expression("Hedges " * g[z] * " (95% CI)"), ylab = "", main = title, pch = 15, col = blue, bty = "n")
  abline(v = 0, lty = 2, col = grey)
  segments(dat$gz_ci_low, y, dat$gz_ci_high, y, col = grey, lwd = 2)
  axis(2, at = c(y, 1), labels = c(dat$project, "Random-effects model"), las = 1, tick = FALSE, cex.axis = 0.80)
  polygon(c(pooled$ci_low, pooled$pooled_hedges_gz, pooled$ci_high, pooled$pooled_hedges_gz), c(1, 1.18, 1, 0.82), col = blue, border = blue)
  text(xr[2], y, labels = paste0(fmt(dat$hedges_gz), " [", fmt(dat$gz_ci_low), ", ", fmt(dat$gz_ci_high), "]"), pos = 2, cex = 0.68)
  text(xr[2], 1, labels = paste0(fmt(pooled$pooled_hedges_gz), " [", fmt(pooled$ci_low), ", ", fmt(pooled$ci_high), "]"), pos = 2, cex = 0.68)
  mtext(bquote(k == .(pooled$k_studies) * ";" ~ tau^2 == .(fmt(pooled$tau2)) * ";" ~ I^2 == .(fmt(pooled$I2, 1)) * "%"), side = 1, line = 3.7, cex = 0.72, col = grey)
}
s3 <- read_pack("SUPPLEMENTARY_FIGURE_S3_bray_meta.csv")
s4 <- read_pack("SUPPLEMENTARY_FIGURE_S4_aitchison_meta.csv")
plot_s3 <- function() plot_meta(s3, s34 %>% filter(scheme == "bray"), "Primary random-effects meta-analysis: Bray-Curtis displacement")
plot_s4 <- function() plot_meta(s4, s34 %>% filter(scheme == "aitchison_CZM"), "Primary random-effects meta-analysis: CZM-Aitchison displacement")
export_manifest <- bind_rows(export_manifest, render_set("Supplementary_Figure_S3_Bray_Meta_Analysis", supp_dir, 10.5, 6.5, plot_s3))
export_manifest <- bind_rows(export_manifest, render_set("Supplementary_Figure_S4_CZM_Aitchison_Meta_Analysis", supp_dir, 10.5, 6.5, plot_s4))

s5 <- read_pack("SUPPLEMENTARY_FIGURE_S5_family_balance.csv")
plot_s5 <- function() {
  y <- rev(seq_len(nrow(s5)))
  xr <- range(c(s5$low, s5$high, 0), na.rm = TRUE)
  par(mar = c(5.0, 10.0, 3.0, 1.2), family = "sans")
  plot(s5$estimate, y, xlim = xr + c(-0.2, 0.2), ylim = c(0.5, nrow(s5) + 0.5), yaxt = "n", xlab = "Mean early-to-late family-balance change", ylab = "",
       main = "Fixed secondary family-balance change", pch = ifelse(s5$type == "Pooled", 18, 16), col = ifelse(s5$type == "External primary", orange, blue), bty = "n")
  abline(v = 0, lty = 2, col = grey)
  segments(s5$low, y, s5$high, y, lwd = 2, col = ifelse(s5$type == "External primary", orange, blue))
  axis(2, at = y, labels = s5$label, las = 1, tick = FALSE, cex.axis = 0.78)
  type_levels <- unique(s5$type)
  legend("bottomright", type_levels,
         pch = ifelse(type_levels == "Pooled", 18, 16),
         col = ifelse(type_levels == "External primary", orange, blue),
         bty = "n", cex = 0.75)
  mtext("Horizontal intervals are pointwise 95% confidence intervals.", side = 1, line = 3.7, cex = 0.72, col = grey)
}
export_manifest <- bind_rows(export_manifest, render_set("Supplementary_Figure_S5_Family_Balance", supp_dir, 10.5, 6.5, plot_s5))

s6 <- read_pack("SUPPLEMENTARY_FIGURE_S6_healthy_reference_paired.csv")
s6_summary <- read_pack("SUPPLEMENTARY_FIGURE_S6_summary.csv")
plot_s6 <- function() {
  yr <- range(c(s6$`Day-3`, s6$`Day-7`), na.rm = TRUE)
  par(mar = c(5.0, 6.5, 3.0, 1.2), family = "sans")
  plot(c(1, 2), yr, type = "n", xaxt = "n", xlab = "", ylab = paste0("Mean Bray-Curtis distance to ", s6_summary$n_healthy_references[1], " healthy references"),
       main = "PRJNA851469 distance to healthy reference")
  axis(1, at = c(1, 2), labels = c("Day 3", "Day 7"))
  for (i in seq_len(nrow(s6))) lines(c(1, 2), c(s6$`Day-3`[i], s6$`Day-7`[i]), col = adjustcolor(grey, 0.55))
  points(rep(1, nrow(s6)), s6$`Day-3`, pch = 16, col = blue)
  points(rep(2, nrow(s6)), s6$`Day-7`, pch = 16, col = orange)
  mtext(paste0("n = ", nrow(s6), " patients; lines connect repeated observations from the same patient and do not represent confidence intervals"), side = 1, line = 3.7, cex = 0.70, col = grey)
}
export_manifest <- bind_rows(export_manifest, render_set("Supplementary_Figure_S6_Healthy_Reference", supp_dir, 9.5, 6.5, plot_s6))

write_csv(export_manifest, file.path(out_dir, "FIGURE_EXPORT_MANIFEST.csv"))
writeLines(capture.output(sessionInfo()), file.path(out_dir, "STEP100B_sessionInfo.txt"))
writeLines(c(
  paste0("Step100B completed: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
  paste0("Input directory: ", source_dir),
  paste0("Output directory: ", out_dir),
  paste0("Exported files: ", nrow(export_manifest)),
  "Inferential analysis executed: NO",
  "Frozen statistics changed: NO",
  "Scientific values read only from Step100A CSV source pack: YES"
), file.path(out_dir, "STEP100B_execution.log"))

cat("Step100B complete. Main and updated supplementary figures written to:\n", out_dir, "\n")
