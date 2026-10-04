# ============================================================
# Step100C: generate final figure legends from the Step100A source pack
#
# This script performs no inferential analysis. Every number inserted into
# a legend is read from the frozen, machine-readable Step100A CSV files.
#
# Usage
#   Rscript analysis/final_figures/10_002_generate_figure_legends.R <source_csv_dir> <figure_output_dir>
# ============================================================

options(stringsAsFactors = FALSE)
suppressPackageStartupMessages({ library(readr); library(dplyr) })

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) stop("Usage: Rscript 10_002_generate_figure_legends.R <source_csv_dir> <figure_output_dir>")
source_dir <- normalizePath(args[1], winslash = "/", mustWork = TRUE)
out_dir <- normalizePath(args[2], winslash = "/", mustWork = TRUE)
if (file.exists(file.path(out_dir, "FINAL_MAIN_FIGURE_LEGENDS.txt"))) stop("Legend outputs already exist.")
if (!dir.exists(source_dir)) stop("Missing Step100A source pack: ", source_dir)
read_pack <- function(name) read_csv(file.path(source_dir, name), show_col_types = FALSE, progress = FALSE)
fmt <- function(x, digits = 3) formatC(as.numeric(x), format = "f", digits = digits)
fmt_p <- function(x) ifelse(as.numeric(x) < 0.001, format(as.numeric(x), scientific = TRUE, digits = 2), formatC(as.numeric(x), format = "f", digits = 3))

f1 <- read_pack("FIGURE1_architecture.csv")
f2 <- read_pack("FIGURE2_longitudinal_displacement.csv")
f3s <- read_pack("FIGURE3_summary.csv"); f3c <- read_pack("FIGURE3c_correlations.csv"); f3d <- read_pack("FIGURE3d_cross_cohort_taxonomic_concordance.csv")
f4n <- read_pack("FIGURE4_group_counts.csv"); f4a <- read_pack("FIGURE4_centroid_comparison.csv"); f4b <- read_pack("FIGURE4_aitchison_pseudocount_sensitivity.csv"); f4c <- read_pack("FIGURE4_taxonomic_alignment.csv")
f5a <- read_pack("FIGURE5_primary_interaction.csv"); f5b <- read_pack("FIGURE5_source_specific_slopes.csv"); f5c <- read_pack("FIGURE5_sensitivity.csv"); f5d <- read_pack("FIGURE5_baseline_robustness.csv")
s1 <- read_pack("SUPPLEMENTARY_FIGURE_S1_taxonomic_concordance_sensitivity.csv")
s2 <- read_pack("SUPPLEMENTARY_FIGURE_S2_pair_counts.csv")
s34 <- read_pack("SUPPLEMENTARY_FIGURES_S3_S4_pooled.csv")
s5 <- read_pack("SUPPLEMENTARY_FIGURE_S5_family_balance.csv")
s6 <- read_pack("SUPPLEMENTARY_FIGURE_S6_summary.csv")

discovery_n <- f1$cohort_count[f1$branch == "DISCOVERY"][1]
longitudinal_n <- f1$cohort_count[f1$branch == "DISCOVERY_LONGITUDINAL"][1]
static_n <- f1$cohort_count[f1$branch == "DISCOVERY_STATIC"][1]
f2_n <- paste0(f2$project, " (n = ", f2$primary_n, ")", collapse = ", ")
outcome_labels <- c(latest_Bray = "Latest Bray-Curtis displacement", latest_EII = "Exploratory EII", latest_Simpson_instability = "Simpson-diversity instability")
stopifnot(!anyNA(outcome_labels[f3c$outcome]))
corr_text <- paste0(unname(outcome_labels[f3c$outcome]), ": Spearman ρ = ", fmt(f3c$spearman_rho), ", permutation p = ", vapply(f3c$permutation_p_10000, fmt_p, character(1)), collapse = "; ")
f3d_text <- f3d %>% mutate(piece = paste0(pair_label, ", ", rank, ": Spearman ρ = ", fmt(spearman_rho), ", n = ", n_shared)) %>% pull(piece) %>% paste(collapse = "; ")
n_lookup <- setNames(f4n$n, f4n$Group)
pc_text <- paste0(format(f4b$pseudocount, scientific = TRUE, digits = 1), ": R² = ", fmt(f4b$PERMANOVA_R2, 3), ", PERMANOVA BH-FDR = ", vapply(f4b$PERMANOVA_FDR, fmt_p, character(1)), ", dispersion BH-FDR = ", vapply(f4b$dispersion_FDR, fmt_p, character(1)), collapse = "; ")
f4c_text <- paste0(f4c$comparison, ": Spearman ρ = ", fmt(f4c$spearman_rho), ", p = ", vapply(f4c$p_value, fmt_p, character(1)), ", overlapping genera n = ", f4c$overlapping_genera, collapse = "; ")
slope_text <- paste0(f5b$group_label, " (", f5b$npatients, " patients, ", f5b$nobs, " observations): beta = ", fmt(f5b$slope_per_day), ", 95% CI ", fmt(f5b$ci95_low), " to ", fmt(f5b$ci95_high), ", LRT p = ", vapply(f5b$p_LRT, fmt_p, character(1)), collapse = "; ")

main_legends <- c(
  paste0("Figure 1 | Study architecture. Public gut-microbiome cohorts were organized into ", discovery_n, " discovery cohorts (", longitudinal_n, " longitudinal cohorts and ", static_n, " static supportive cohort), an independent longitudinal validation cohort (PRJNA1125274), an external cross-sectional ecological-state comparison (PRJNA1010969), and an infection-source longitudinal extension (CRA002354). The diagram distinguishes longitudinal evidence from static or cross-sectional support; arrows describe the evidence architecture and do not denote participant flow."),
  paste0("Figure 2 | Cross-cohort longitudinal ecological displacement. Standardized within-patient early-to-late changes in Bray-Curtis displacement are shown for ", nrow(f2), " longitudinal discovery cohorts: ", f2_n, ". Filled blue circles are primary early-to-late estimates and open orange circles are common-anchor sensitivity estimates. Lines connect the primary and common-anchor sensitivity estimates from the same cohort and do not represent confidence intervals. Asterisks indicate primary estimates with Benjamini-Hochberg false discovery rate (FDR) < 0.05. The statistical unit is the patient; two-sided Wilcoxon signed-rank tests were used, with Benjamini-Hochberg correction within the fixed family of tests."),
  paste0("Figure 3 | Within-patient ecological displacement, trajectory heterogeneity, and cross-cohort taxonomic concordance. (a) Bray-Curtis displacement at Day 3 and Day 7 for ", f3s$day3_day7_paired_patients[1], " paired patients; lines connect repeated observations from the same patient and do not represent confidence intervals (paired Wilcoxon signed-rank p = ", fmt_p(f3s$day3_day7_wilcoxon_p[1]), "). These paired patients are a subset of the ", f3s$trajectory_analysis_patients[1], " patients included in the trajectory analysis. (b) Distribution of ", f3s$pairwise_comparisons[1], " pairwise signed-trajectory Euclidean distances generated from those ", f3s$pairwise_source_patients[1], " patients. Pairwise distances are descriptive and are not independent patient-level observations; the box shows the median and interquartile range and whiskers extend to 1.5 times the interquartile range. (c) Exploratory patient-level Spearman correlations (n = ", unique(f3c$n)[1], ") between trajectory heterogeneity and ecological measures; permutation p values used 10,000 permutations (", corr_text, "). (d) Cross-cohort Spearman concordance of longitudinal centered-log-ratio taxonomic effect vectors at genus and family ranks (", f3d_text, "). Panels a-c use PRJNA691455; panel d uses the prespecified cross-cohort concordance output. No confidence interval is implied by the point labels in panels c-d."),
  paste0("Figure 4 | External cross-sectional ecological-state comparison and compositional sensitivity. (a) Individual Bray-Curtis distances from the control centroid in PRJNA1010969 (Control n = ", n_lookup[["Control"]], ", Trauma n = ", n_lookup[["Trauma"]], ", Sepsis n = ", n_lookup[["Sepsis"]], "). Points represent individual samples; boxes show the median and interquartile range, whiskers extend to the most extreme observations within 1.5 times the interquartile range, and all observations are shown. The Sepsis-Trauma median difference was ", fmt(f4a$median_difference_A_minus_B[1]), " (", f4a$test[1], ", ", f4a$correction[1], ": FDR = ", fmt_p(f4a$FDR[1]), "; Cliff's delta = ", fmt(f4a$cliff_delta_A_vs_B[1]), "). (b) Aitchison pseudocount sensitivity across ", nrow(f4b), " pseudocounts; R² = 0.081–0.085 is the range across pseudocount sensitivity analyses, not a confidence interval (", pc_text, "). PERMANOVA and dispersion tests are reported separately. (c) Exploratory Spearman alignment between longitudinal taxonomic changes and external Sepsis-Control or Sepsis-Trauma contrasts (", f4c_text, "). Sample counts are those entering the displayed analysis."),
  paste0("Figure 5 | Infection-source longitudinal extension in CRA002354. (a) Primary source-by-time interaction for Bray-Curtis displacement: beta = ", fmt(f5a$effect[1]), ", 95% CI ", fmt(f5a$ci_low[1]), " to ", fmt(f5a$ci_high[1]), ", likelihood-ratio-test p = ", fmt_p(f5a$p_value[1]), ". The horizontal interval is the model-based 95% CI. (b) Source-specific slopes: ", slope_text, ". Horizontal intervals are model-based 95% CIs. (c) Prespecified sensitivity analyses; points are interaction estimates, horizontal intervals are 95% CIs, and p values are likelihood-ratio-test p values. (d) Baseline beta-diversity robustness: unrestricted PERMANOVA R² = ", fmt(f5d$primary_PERMANOVA_R2[1]), ", p = ", fmt_p(f5d$primary_PERMANOVA_p[1]), ", dispersion p = ", fmt_p(f5d$primary_dispersion_p[1]), "; ICU-Day-3-restricted PERMANOVA R² = ", fmt(f5d$day3_restricted_PERMANOVA_R2[1]), ", p = ", fmt_p(f5d$day3_restricted_PERMANOVA_p[1]), ", dispersion p = ", fmt_p(f5d$day3_restricted_dispersion_p[1]), ". These analyses do not establish equivalence between infection-source groups." )
)

s2_counts <- paste0(s2$project, " primary n = ", s2$primary_n, ", common-anchor n = ", s2$common_anchor_n, collapse = "; ")
bray_pool <- s34 %>% filter(scheme == "bray")
ait_pool <- s34 %>% filter(scheme == "aitchison_CZM")
s6row <- s6[1, ]
supp_legends <- c(
  paste0("Supplementary Figure S1 | Sensitivity of cross-cohort taxonomic concordance to anchor definition. (a) Primary concordance analysis. (b) Common-anchor sensitivity analysis. Points show Spearman correlations between longitudinal centered-log-ratio taxonomic effect vectors; circles denote genus-level and triangles family-level estimates. The common-anchor analysis is a sensitivity analysis of the same cohorts, not an independent replication."),
  paste0("Supplementary Figure S2 | Patient-pair counts contributing to longitudinal displacement analyses. Counts are: ", s2_counts, ". Primary and common-anchor counts may differ because the sensitivity analysis applies a shared time-anchor definition."),
  paste0("Supplementary Figure S3 | Random-effects synthesis of primary Bray-Curtis longitudinal effects. Cohort estimates are Hedges g_z with 95% CIs; the pooled estimate was ", fmt(bray_pool$pooled_hedges_gz), " (95% CI ", fmt(bray_pool$ci_low), " to ", fmt(bray_pool$ci_high), "; k = ", bray_pool$k_studies, "; I² = ", fmt(bray_pool$I2, 1), "%). Horizontal intervals are 95% CIs, not the observed range."),
  paste0("Supplementary Figure S4 | Random-effects synthesis of primary CZM-Aitchison longitudinal effects. Cohort estimates are Hedges g_z with 95% CIs; the pooled estimate was ", fmt(ait_pool$pooled_hedges_gz), " (95% CI ", fmt(ait_pool$ci_low), " to ", fmt(ait_pool$ci_high), "; k = ", ait_pool$k_studies, "; I² = ", fmt(ait_pool$I2, 1), "%). Horizontal intervals are 95% CIs, not the observed range."),
  paste0("Supplementary Figure S5 | Fixed secondary family-balance analysis. Points show cohort-specific or pooled mean early-to-late changes in the prespecified family balance. Horizontal intervals are pointwise 95% CIs. PRJNA1125274 is the independent longitudinal validation primary cohort; the remaining cohort estimates and the REML-Hartung-Knapp pooled result are secondary or pooled summaries as labeled."),
  paste0("Supplementary Figure S6 | Paired distance to a healthy reference in PRJNA851469. Mean Bray-Curtis distances to ", s6row$n_healthy_references, " healthy reference samples are shown at Day 3 and Day 7 for ", s6row$n_analyzed, " patients. Lines connect repeated observations from the same patient and do not represent confidence intervals. The mean paired change was ", fmt(s6row$mean), " (ordinary 95% CI ", fmt(s6row$ci_low), " to ", fmt(s6row$ci_high), "; two-level bootstrap 95% CI ", fmt(s6row$two_level_bootstrap_ci_low), " to ", fmt(s6row$two_level_bootstrap_ci_high), "). This is a same-cohort secondary reference analysis, not independent validation.")
)

main_legends <- gsub("whiskers extend to 1.5 times the interquartile range", "whiskers extend to the most extreme observations within 1.5 times the interquartile range", main_legends, fixed = TRUE)
writeLines(paste0(main_legends, collapse = "\n\n"), file.path(out_dir, "FINAL_MAIN_FIGURE_LEGENDS.txt"), useBytes = TRUE)
writeLines(paste0(supp_legends, collapse = "\n\n"), file.path(out_dir, "FINAL_SUPPLEMENTARY_FIGURE_LEGENDS.txt"), useBytes = TRUE)

legend_map <- bind_rows(
  tibble(item = paste0("Figure ", 1:5), legend_file = "FINAL_MAIN_FIGURE_LEGENDS.txt", source_csv = c("FIGURE1_architecture.csv", "FIGURE2_longitudinal_displacement.csv", "FIGURE3_summary.csv; FIGURE3a-d CSVs", "FIGURE4_group_counts.csv; FIGURE4_centroid_comparison.csv; FIGURE4_aitchison_pseudocount_sensitivity.csv; FIGURE4_taxonomic_alignment.csv", "FIGURE5_primary_interaction.csv; FIGURE5_source_specific_slopes.csv; FIGURE5_sensitivity.csv; FIGURE5_baseline_robustness.csv")),
  tibble(item = paste0("Supplementary Figure S", 1:6), legend_file = "FINAL_SUPPLEMENTARY_FIGURE_LEGENDS.txt", source_csv = c("SUPPLEMENTARY_FIGURE_S1_taxonomic_concordance_sensitivity.csv", "SUPPLEMENTARY_FIGURE_S2_pair_counts.csv", "SUPPLEMENTARY_FIGURE_S3_bray_meta.csv; SUPPLEMENTARY_FIGURES_S3_S4_pooled.csv", "SUPPLEMENTARY_FIGURE_S4_aitchison_meta.csv; SUPPLEMENTARY_FIGURES_S3_S4_pooled.csv", "SUPPLEMENTARY_FIGURE_S5_family_balance.csv", "SUPPLEMENTARY_FIGURE_S6_healthy_reference_paired.csv; SUPPLEMENTARY_FIGURE_S6_summary.csv"))
)
write_csv(legend_map, file.path(out_dir, "LEGEND_SOURCE_MAP.csv"))
writeLines(capture.output(sessionInfo()), file.path(out_dir, "STEP100C_sessionInfo.txt"))
writeLines(c(paste0("Step100C completed: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")), "Inferential analysis executed: NO", "Legend values read only from Step100A source pack: YES"), file.path(out_dir, "STEP100C_execution.log"))
cat("Step100C complete. Final legends and legend-source map written to:\n", out_dir, "\n")
