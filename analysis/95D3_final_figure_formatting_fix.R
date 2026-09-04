# ============================================================
# Sepsis V2 - Step95D3
# FINAL FIGURE FORMATTING FIX
#
# Supersedes Step95D2 only for manuscript-facing figure files.
#
# NO new statistics, NO changed estimates, NO changed conclusions.
#
# Decisions from visual review:
# - Figure 1: KEEP Step95D2 unchanged.
# - Figure 2: KEEP Step95D2 unchanged.
# - Figure 3: formatting-only repair of panel C text clipping.
# - Figure 4: KEEP Step95D2 unchanged.
# - Figure 5: formatting-only repair of panel A/C left-label clipping.
# ============================================================

options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tibble)
  library(stringr)
})

ROOT <- "E:/sepsis_project"
RESULTS <- file.path(ROOT, "results")

D2 <- file.path(
  RESULTS,
  "V2_35D2_PUBLICATION_FIGURES_SOURCEPACK_ONLY"
)

SRC <- file.path(
  RESULTS,
  "V2_35B_MANUSCRIPT_SOURCE_PACK_ASSEMBLY"
)

OUT <- file.path(
  RESULTS,
  "V2_35D3_FINAL_FIGURE_FORMATTING_FIX"
)

FIGDIR <- file.path(OUT, "01_FINAL_FIGURES")
QCDIR <- file.path(OUT, "02_QC_AND_CAPTIONS")
PROV <- file.path(OUT, "03_PROVENANCE")

for (d in c(OUT, FIGDIR, QCDIR, PROV)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

if (!file.exists(file.path(D2, "_STEP95D2_COMPLETE.ok"))) {
  stop("Step95D2 is incomplete.")
}

if (!file.exists(file.path(SRC, "_STEP95B_COMPLETE.ok"))) {
  stop("Step95B source pack is incomplete.")
}

# ------------------------------------------------------------
# 1. Copy approved Step95D2 figures unchanged
# ------------------------------------------------------------

approved <- c(
  "Figure_1_Study_Architecture",
  "Figure_2_Cross_Cohort_Longitudinal_Displacement",
  "Figure_4_External_Validation_Trauma_Comparison"
)

copy_audit <- tibble()

for (nm in approved) {
  for (ext in c("pdf","png")) {
    src <- file.path(
      D2,
      "01_FINAL_FIGURES",
      paste0(nm, ".", ext)
    )

    dst <- file.path(
      FIGDIR,
      paste0(nm, ".", ext)
    )

    ok <- file.exists(src) &&
      file.copy(src, dst, overwrite = TRUE)

    copy_audit <- bind_rows(
      copy_audit,
      tibble(
        figure = nm,
        format = ext,
        source = src,
        copied = ok
      )
    )
  }
}

if (!all(copy_audit$copied)) {
  stop("Failed to copy one or more approved Step95D2 figures.")
}

write_csv(
  copy_audit,
  file.path(PROV, "01_APPROVED_FIGURE_COPY_AUDIT.csv")
)

# ------------------------------------------------------------
# 2. Read frozen constants from Step95D2
# ------------------------------------------------------------

C <- read_csv(
  file.path(
    D2,
    "03_PROVENANCE",
    "01_FROZEN_CONSTANTS_USED.csv"
  ),
  show_col_types = FALSE
)

val <- function(metric) {
  C$value[C$metric == metric][1]
}

safe <- function(x, digits = 3) {
  ifelse(
    is.na(x),
    "NA",
    formatC(x, format = "f", digits = digits)
  )
}

# ------------------------------------------------------------
# 3. Figure 3 formatting-only rebuild
# ------------------------------------------------------------

plot_fig3 <- function() {

  layout(
    matrix(1:3, nrow = 1),
    widths = c(1.0, 1.05, 1.25)
  )

  # A
  par(mar = c(4.8,4.8,2.5,1.2))

  vals <- c(
    val("Day3_Bray_median"),
    val("Day7_Bray_median")
  )

  bp <- barplot(
    vals,
    names.arg = c("Day 3","Day 7"),
    ylim = c(0,1),
    ylab = "Median Bray-Curtis displacement",
    main = "A. Progressive ecological displacement"
  )

  segments(
    bp[1], vals[1],
    bp[2], vals[2],
    lty = 2
  )

  text(
    mean(bp),
    0.12,
    paste0(
      "paired n=",
      as.integer(val("Day3_to_Day7_paired_n")),
      "\nmean paired change=",
      safe(val("Day3_to_Day7_mean_change"),3)
    ),
    cex = 0.74
  )

  # B
  par(mar = c(4.8,3.4,2.5,1.2))

  plot(
    NA,
    xlim = c(0,1),
    ylim = c(0,4),
    axes = FALSE,
    xlab = "",
    ylab = "",
    main = "B. Taxonomic-route heterogeneity"
  )

  rect(0.03,0.42,0.97,3.68)

  text(
    0.07,3.18,
    paste0(
      "Pairwise Euclidean distance\nmean ",
      safe(val("pairwise_signed_Euclidean_mean"),3),
      "; median ",
      safe(val("pairwise_signed_Euclidean_median"),3)
    ),
    adj = c(0,0.5),
    cex = 0.79
  )

  text(
    0.07,2.08,
    paste0(
      "Mean cosine similarity\n",
      safe(val("mean_cosine_similarity"),3)
    ),
    adj = c(0,0.5),
    cex = 0.79
  )

  text(
    0.07,0.98,
    paste0(
      "Patient pairs with cosine <= 0\n",
      safe(100*val("fraction_cosine_le_0"),1),
      "%"
    ),
    adj = c(0,0.5),
    cex = 0.79
  )

  # C — key repair: enlarge right x-range and draw p labels inside panel.
  par(mar = c(4.8,7.2,2.5,2.2))

  rho <- c(
    val("rho_latest_Bray"),
    val("rho_latest_EII"),
    val("rho_Simpson_instability")
  )

  labs <- c(
    "Latest Bray displacement",
    "Latest EII",
    "Latest Simpson instability"
  )

  ptxt <- c(
    paste0("p=",safe(val("p_latest_Bray"),3)),
    paste0(
      "p=",safe(val("p_latest_EII"),3),
      "; perm=",safe(val("perm_p_latest_EII"),3)
    ),
    paste0(
      "p=",safe(val("p_Simpson_instability"),3),
      "; perm=",safe(val("perm_p_Simpson_instability"),3)
    )
  )

  yy <- c(3,2,1)

  plot(
    rho, yy,
    xlim = c(-0.05,1.05),
    ylim = c(0.5,3.5),
    yaxt = "n",
    pch = 16,
    xlab = "Spearman rho",
    ylab = "",
    main = "C. Ecological consistency",
    bty = "n"
  )

  abline(v=0,lty=2)

  axis(
    2,
    at = yy,
    labels = labs,
    las = 1,
    tick = FALSE,
    cex.axis = 0.70
  )

  # fixed right-aligned label location prevents clipping
  text(
    rep(1.00,3),
    yy,
    labels = ptxt,
    adj = c(1,0.5),
    cex = 0.65
  )
}

pdf(
  file.path(
    FIGDIR,
    "Figure_3_Core_Sepsis_Heterogeneous_Taxonomic_Routes.pdf"
  ),
  width = 12.2,
  height = 5.5,
  family = "sans",
  useDingbats = FALSE
)
plot_fig3()
dev.off()

png(
  file.path(
    FIGDIR,
    "Figure_3_Core_Sepsis_Heterogeneous_Taxonomic_Routes.png"
  ),
  width = 2440,
  height = 1100,
  res = 200
)
plot_fig3()
dev.off()

# ------------------------------------------------------------
# 4. Figure 5 formatting-only rebuild
# ------------------------------------------------------------

ECO <- read_csv(
  file.path(
    SRC,
    "02_SUPPLEMENTARY_TABLE_SOURCES",
    "SUPP_S4_SOURCE_ECOLOGY_REGISTRY.csv"
  ),
  show_col_types = FALSE
)

SENS <- read_csv(
  file.path(
    SRC,
    "02_SUPPLEMENTARY_TABLE_SOURCES",
    "SUPP_S5_SOURCE_BRAY_SENSITIVITY.csv"
  ),
  show_col_types = FALSE
)

SLOPES <- read_csv(
  file.path(
    SRC,
    "02_SUPPLEMENTARY_TABLE_SOURCES",
    "SUPP_S6_SOURCE_SLOPES.csv"
  ),
  show_col_types = FALSE
)

BASEB <- read_csv(
  file.path(
    SRC,
    "02_SUPPLEMENTARY_TABLE_SOURCES",
    "SUPP_S7_SOURCE_BASELINE_BETA.csv"
  ),
  show_col_types = FALSE
)

primary <- ECO %>%
  filter(
    domain == "PRIMARY_LONGITUDINAL_BRAY",
    estimand == "source_by_personal_time_interaction"
  ) %>%
  slice(1)

pretty_sens <- c(
  "Age + sex + SOFA",
  "Age + sex + APACHE II",
  "Age + sex + lactate",
  "Exclude OTHER_UNKNOWN",
  "Rarefied 4000 Bray",
  "Min-depth 2000"
)

if (length(pretty_sens) != nrow(SENS)) {
  pretty_sens <- gsub("_"," ",SENS$sensitivity)
}

plot_fig5 <- function() {

  layout(
    matrix(1:4, nrow=2, byrow=TRUE),
    widths=c(1,1),
    heights=c(1,1)
  )

  # A — larger left margin fixes clipped source label
  par(mar=c(4.7,8.3,2.5,1.3))

  plot(
    primary$effect,1,
    xlim=c(
      min(primary$ci_low,0)-0.01,
      max(primary$ci_high,0)+0.01
    ),
    ylim=c(0.5,1.5),
    yaxt="n",
    pch=16,
    xlab="Source-by-time interaction beta",
    ylab="",
    main="A. Primary Bray interaction",
    bty="n"
  )

  abline(v=0,lty=2)
  segments(primary$ci_low,1,primary$ci_high,1,lwd=2)

  axis(
    2,
    at=1,
    labels="Pulmonary vs recorded\nnon-pulmonary",
    las=1,
    tick=FALSE,
    cex.axis=0.72
  )

  text(
    mean(c(primary$ci_low,primary$ci_high)),
    0.67,
    paste0(
      "beta=",safe(primary$effect,3),
      "; 95% CI ",
      safe(primary$ci_low,3),
      " to ",
      safe(primary$ci_high,3),
      "; p=",safe(primary$p_value,3)
    ),
    cex=0.69
  )

  # B
  par(mar=c(4.7,7.5,2.5,2.0))

  yy <- c(2,1)

  plot(
    SLOPES$slope_per_day,
    yy,
    xlim=c(
      min(c(SLOPES$ci95_low,0))-0.005,
      max(c(SLOPES$ci95_high,0))+0.015
    ),
    ylim=c(0.5,2.5),
    yaxt="n",
    pch=16,
    xlab="Estimated Bray displacement slope per day",
    ylab="",
    main="B. Source-specific slopes",
    bty="n"
  )

  abline(v=0,lty=2)

  segments(
    SLOPES$ci95_low,
    yy,
    SLOPES$ci95_high,
    yy,
    lwd=2
  )

  axis(
    2,
    at=yy,
    labels=c(
      "Pulmonary",
      "Recorded non-pulmonary"
    ),
    las=1,
    tick=FALSE,
    cex.axis=0.73
  )

  text(
    SLOPES$ci95_high,
    yy,
    labels=paste0("p=",safe(SLOPES$p_LRT,3)),
    pos=4,
    cex=0.68
  )

  # C — larger left margin + readable sensitivity labels
  par(mar=c(4.7,9.2,2.5,1.3))

  y3 <- rev(seq_len(nrow(SENS)))

  xlo <- min(c(SENS$ci_low,0),na.rm=TRUE)-0.005
  xhi <- max(c(SENS$ci_high,0),na.rm=TRUE)+0.01

  plot(
    SENS$beta,
    y3,
    xlim=c(xlo,xhi),
    ylim=c(0.5,nrow(SENS)+0.5),
    yaxt="n",
    pch=16,
    xlab="Source-by-time interaction beta",
    ylab="",
    main="C. Key sensitivity analyses",
    bty="n"
  )

  abline(v=0,lty=2)

  segments(
    SENS$ci_low,
    y3,
    SENS$ci_high,
    y3,
    lwd=1.6
  )

  axis(
    2,
    at=y3,
    labels=pretty_sens,
    las=1,
    tick=FALSE,
    cex.axis=0.68
  )

  # D
  par(mar=c(4.7,5.5,2.5,1.5))

  vals <- c(
    BASEB$primary_PERMANOVA_R2[1],
    BASEB$day3_restricted_PERMANOVA_R2[1]
  )

  bp <- barplot(
    vals,
    names.arg=c(
      "Unrestricted",
      "<= ICU Day 3"
    ),
    ylim=c(0,max(vals)*2.1),
    ylab="PERMANOVA R2",
    main="D. Baseline beta robustness"
  )

  text(
    bp,
    vals,
    labels=c(
      paste0(
        "PERMANOVA p=",safe(BASEB$primary_PERMANOVA_p[1],3),
        "\ndispersion p=",safe(BASEB$primary_dispersion_p[1],3)
      ),
      paste0(
        "PERMANOVA p=",safe(BASEB$day3_restricted_PERMANOVA_p[1],3),
        "\ndispersion p=",safe(BASEB$day3_restricted_dispersion_p[1],3)
      )
    ),
    pos=3,
    cex=0.65
  )
}

pdf(
  file.path(
    FIGDIR,
    "Figure_5_Infection_Source_Longitudinal_Ecology.pdf"
  ),
  width=11.8,
  height=8.1,
  family="sans",
  useDingbats=FALSE
)
plot_fig5()
dev.off()

png(
  file.path(
    FIGDIR,
    "Figure_5_Infection_Source_Longitudinal_Ecology.png"
  ),
  width=2360,
  height=1620,
  res=200
)
plot_fig5()
dev.off()

# ------------------------------------------------------------
# 5. Copy captions unchanged from Step95D2
# ------------------------------------------------------------

cap_src <- file.path(
  D2,
  "02_CAPTIONS_AND_QC",
  "DRAFT_MAIN_FIGURE_CAPTIONS.txt"
)

if (!file.copy(
  cap_src,
  file.path(QCDIR, "DRAFT_MAIN_FIGURE_CAPTIONS.txt"),
  overwrite=TRUE
)) {
  stop("Could not copy frozen Step95D2 figure captions.")
}

writeLines(
  c(
    "STEP95D3 FINAL FORMATTING DECISIONS",
    "",
    "Figure 1: approved unchanged from Step95D2.",
    "Figure 2: approved unchanged from Step95D2.",
    "Figure 3: only panel spacing/right annotation placement changed; all frozen values unchanged.",
    "Figure 4: approved unchanged from Step95D2.",
    "Figure 5: only margins/axis-label presentation changed; all frozen values unchanged.",
    "",
    "No statistical model was rerun.",
    "No p value, FDR, effect estimate, CI, evidence tier, or manuscript conclusion was changed.",
    "",
    "Step95D3 is the final manuscript-facing Figure 1-5 source set."
  ),
  file.path(QCDIR, "STEP95D3_FORMATTING_QC_REPORT.txt")
)

# ------------------------------------------------------------
# 6. Final audit
# ------------------------------------------------------------

expected <- unlist(
  lapply(
    c(
      "Figure_1_Study_Architecture",
      "Figure_2_Cross_Cohort_Longitudinal_Displacement",
      "Figure_3_Core_Sepsis_Heterogeneous_Taxonomic_Routes",
      "Figure_4_External_Validation_Trauma_Comparison",
      "Figure_5_Infection_Source_Longitudinal_Ecology"
    ),
    function(x) {
      file.path(
        FIGDIR,
        paste0(x,c(".pdf",".png"))
      )
    }
  )
)

audit <- tibble(
  file=expected,
  exists=file.exists(expected),
  bytes=ifelse(
    file.exists(expected),
    file.info(expected)$size,
    NA_real_
  )
)

write_csv(
  audit,
  file.path(PROV, "02_FINAL_FIGURE_FILE_AUDIT.csv")
)

ready <- all(audit$exists)

write_csv(
  tibble(
    figure1_to_5_PDF_and_PNG_complete=ready,
    figure1_approved=TRUE,
    figure2_approved=TRUE,
    figure3_clipping_fixed=TRUE,
    figure4_approved=TRUE,
    figure5_clipping_fixed=TRUE,
    no_new_inferential_analysis=TRUE,
    ready_for_table1_and_results=ready
  ),
  file.path(OUT,"STEP95D3_READINESS.csv")
)

if (!ready) {
  stop("Final figure file audit failed.")
}

writeLines(
  c(
    paste0("Completed: ",Sys.time()),
    "Figure 1-5 PDF+PNG complete: TRUE",
    "Formatting-only changes: TRUE",
    "No new inferential analysis: TRUE",
    "Ready for Table 1 + Results: TRUE",
    "STEP95D3 COMPLETE"
  ),
  file.path(OUT,"_STEP95D3_COMPLETE.ok")
)

cat("STEP95D3 COMPLETE\n")
