# ============================================================
# Sepsis V2 - Step95D2
# PUBLICATION FIGURES - SOURCE-PACK-ONLY ROBUST REBUILD
#
# Supersedes incomplete Step95D.
#
# Principle:
# - NO new inferential analysis.
# - Use Step95B / Step95C frozen source pack wherever possible.
# - A small number of Step93U / Step93X3 frozen summary constants that were
#   already explicitly frozen upstream are recorded in an auditable constants
#   table rather than recalculated.
# - Failure of optional PNG rendering does NOT invalidate PDF production.
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

SRC <- file.path(
  RESULTS,
  "V2_35B_MANUSCRIPT_SOURCE_PACK_ASSEMBLY"
)

C95 <- file.path(
  RESULTS,
  "V2_35C_MANUSCRIPT_ARCHITECTURE_FIGURE1_AND_MAIN_TABLES"
)

OUT <- file.path(
  RESULTS,
  "V2_35D2_PUBLICATION_FIGURES_SOURCEPACK_ONLY"
)

FIGDIR <- file.path(OUT, "01_FINAL_FIGURES")
CAPDIR <- file.path(OUT, "02_CAPTIONS_AND_QC")
AUDDIR <- file.path(OUT, "03_PROVENANCE")

for (d in c(OUT, FIGDIR, CAPDIR, AUDDIR)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

LOG <- file.path(OUT, "_STEP95D2_runtime_checkpoints.txt")
ERR <- file.path(OUT, "_STEP95D2_FATAL_ERROR.txt")

if (file.exists(LOG)) unlink(LOG)
if (file.exists(ERR)) unlink(ERR)

ck <- function(x) {
  cat(
    paste0(x, ": ", Sys.time(), "\n"),
    file = LOG,
    append = TRUE
  )
}

fatal <- function(msg) {
  writeLines(
    c(
      paste0("STEP95D2 FAILED: ", Sys.time()),
      msg
    ),
    ERR
  )
  stop(msg, call. = FALSE)
}

tryCatch({

  ck("STEP95D2 STARTED")

  if (!file.exists(file.path(SRC, "_STEP95B_COMPLETE.ok"))) {
    fatal("Step95B source pack is incomplete.")
  }

  if (!file.exists(file.path(C95, "_STEP95C_COMPLETE.ok"))) {
    fatal("Step95C is incomplete.")
  }

  # ----------------------------------------------------------
  # Authoritative source-pack inputs
  # ----------------------------------------------------------

  W <- read_csv(
    file.path(
      C95,
      "Main_Table_2_Cross_Cohort_Longitudinal_Results.csv"
    ),
    show_col_types = FALSE
  )

  CORE_SUM <- read_csv(
    file.path(
      SRC,
      "01_MAIN_TABLE_SOURCES",
      "MAIN_T3_CORE_TRAJECTORY_ENDPOINTS.csv"
    ),
    show_col_types = FALSE
  )

  EXT_TAX <- read_csv(
    file.path(
      SRC,
      "01_MAIN_TABLE_SOURCES",
      "MAIN_T4_EXTERNAL_VALIDATION.csv"
    ),
    show_col_types = FALSE
  )

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

  GFREEZE <- read_csv(
    file.path(
      SRC,
      "02_SUPPLEMENTARY_TABLE_SOURCES",
      "SUPP_S10_GENUS_BRANCH_FREEZE.csv"
    ),
    show_col_types = FALSE
  )

  ck("SOURCE PACK READ")

  # ----------------------------------------------------------
  # Frozen constants not exposed as standalone tables in 95B.
  # These values were already frozen upstream and are NOT
  # recalculated here.
  # ----------------------------------------------------------

  frozen_constants <- tribble(
    ~branch, ~metric, ~value, ~source_note,
    "STEP93U", "Day3_Bray_median", 0.725,
      "Frozen Step93U/Step88A2 result",
    "STEP93U", "Day7_Bray_median", 0.9038,
      "Frozen Step93U result",
    "STEP93U", "Day3_to_Day7_paired_n", 9,
      "Frozen Step93U result",
    "STEP93U", "Day3_to_Day7_mean_change", 0.112914,
      "Frozen Step93U result",
    "STEP93U", "pairwise_signed_Euclidean_mean", 0.620,
      "Frozen Step93U result",
    "STEP93U", "pairwise_signed_Euclidean_median", 0.554,
      "Frozen Step93U result",
    "STEP93U", "mean_cosine_similarity", 0.162,
      "Frozen Step93U result",
    "STEP93U", "fraction_cosine_le_0", 0.356,
      "Frozen Step93U result",
    "STEP93U", "rho_latest_Bray", 0.236,
      "Frozen Step93U result",
    "STEP93U", "p_latest_Bray", 0.511,
      "Frozen Step93U result",
    "STEP93U", "rho_latest_EII", 0.733,
      "Frozen Step93U result",
    "STEP93U", "p_latest_EII", 0.0158,
      "Frozen Step93U result",
    "STEP93U", "perm_p_latest_EII", 0.0234,
      "Frozen Step93U result",
    "STEP93U", "rho_Simpson_instability", 0.673,
      "Frozen Step93U result",
    "STEP93U", "p_Simpson_instability", 0.033,
      "Frozen Step93U result",
    "STEP93U", "perm_p_Simpson_instability", 0.038,
      "Frozen Step93U result",
    "STEP93X", "Trauma_healthy_centroid_distance", 0.657,
      "Frozen Step93X2 result",
    "STEP93X", "Sepsis_healthy_centroid_distance", 0.860,
      "Frozen Step93X2 result",
    "STEP93X", "Sepsis_minus_Trauma_centroid_difference", 0.203,
      "Frozen Step93X2 result",
    "STEP93X", "centroid_difference_FDR", 0.001321,
      "Frozen Step93X2 result",
    "STEP93X3", "CLR_Aitchison_R2_min", 0.081,
      "Frozen pseudocount sensitivity range",
    "STEP93X3", "CLR_Aitchison_R2_max", 0.085,
      "Frozen pseudocount sensitivity range",
    "STEP93X3", "CLR_Aitchison_PERMANOVA_FDR_min", 0.0002,
      "Frozen pseudocount sensitivity range",
    "STEP93X3", "CLR_Aitchison_PERMANOVA_FDR_max", 0.0005,
      "Frozen pseudocount sensitivity range",
    "STEP93X3", "CLR_Aitchison_dispersion_FDR_min", 0.168,
      "Frozen pseudocount sensitivity range",
    "STEP93X3", "CLR_Aitchison_dispersion_FDR_max", 0.234,
      "Frozen pseudocount sensitivity range"
  )

  write_csv(
    frozen_constants,
    file.path(AUDDIR, "01_FROZEN_CONSTANTS_USED.csv")
  )

  val <- function(metric) {
    frozen_constants$value[
      frozen_constants$metric == metric
    ][1]
  }

  safe <- function(x, digits = 3) {
    ifelse(
      is.na(x),
      "NA",
      formatC(x, format = "f", digits = digits)
    )
  }

  # ----------------------------------------------------------
  # Device helpers
  # ----------------------------------------------------------

  open_pdf <- function(path, width, height) {
    pdf(
      path,
      width = width,
      height = height,
      family = "sans",
      useDingbats = FALSE
    )
  }

  open_png_safe <- function(path, width, height, res = 200) {
    ok <- TRUE
    tryCatch(
      {
        png(
          path,
          width = width,
          height = height,
          res = res
        )
      },
      error = function(e) {
        ok <<- FALSE
      }
    )
    ok
  }

  # One plotting function can be called for PDF and PNG.
  make_both <- function(name, width_pdf, height_pdf,
                        width_png, height_png, plot_fun) {

    pdf_path <- file.path(
      FIGDIR,
      paste0(name, ".pdf")
    )

    png_path <- file.path(
      FIGDIR,
      paste0(name, ".png")
    )

    open_pdf(
      pdf_path,
      width_pdf,
      height_pdf
    )
    plot_fun()
    dev.off()

    png_ok <- open_png_safe(
      png_path,
      width_png,
      height_png
    )

    if (png_ok) {
      tryCatch(
        {
          plot_fun()
          dev.off()
        },
        error = function(e) {
          try(dev.off(), silent = TRUE)
          png_ok <<- FALSE
          if (file.exists(png_path)) unlink(png_path)
        }
      )
    }

    tibble(
      figure = name,
      pdf = pdf_path,
      pdf_exists = file.exists(pdf_path),
      png = png_path,
      png_exists = file.exists(png_path),
      png_attempt_ok = png_ok
    )
  }

  figure_audit <- tibble()

  # ==========================================================
  # FIGURE 1
  # Correct architecture:
  # - Core repeated sepsis nested under longitudinal branch.
  # - External validation and infection-source extension are
  #   parallel extensions, NOT parent-child.
  # - No crossing arrows.
  # ==========================================================

  plot_fig1 <- function() {

    par(mar = c(0.4,0.4,0.6,0.4))

    plot(
      NA,
      xlim = c(0,12),
      ylim = c(0,8),
      axes = FALSE,
      xlab = "",
      ylab = "",
      xaxs = "i",
      yaxs = "i"
    )

    # Top
    rect(3.8,7.05,8.2,7.72,lwd=1.4)
    text(
      6,7.39,
      "Public gut-microbiome cohorts in critical illness",
      font=2,
      cex=0.98
    )

    # Three parallel evidence branches
    rect(0.45,5.20,4.15,6.32,lwd=1.25)
    text(2.30,5.92,"Longitudinal critical-illness framework",font=2,cex=0.88)
    text(2.30,5.53,"7 cohorts; within-patient ecological change",cex=0.78)

    rect(4.30,5.20,7.95,6.32,lwd=1.25)
    text(6.125,5.92,"External ecological-state validation",font=2,cex=0.88)
    text(6.125,5.53,"Control vs severe trauma vs sepsis",cex=0.78)

    rect(8.10,5.20,11.65,6.32,lwd=1.25)
    text(9.875,5.92,"Infection-source longitudinal extension",font=2,cex=0.84)
    text(9.875,5.53,"Pulmonary vs recorded non-pulmonary sepsis",cex=0.74)

    arrows(6,7.05,2.30,6.32,length=0.07)
    arrows(6,7.05,6.125,6.32,length=0.07)
    arrows(6,7.05,9.875,6.32,length=0.07)

    # Nested core-sepsis analysis
    rect(0.85,3.40,3.75,4.47,lwd=1.15)
    text(2.30,4.04,"Core repeated-sepsis cohort",font=2,cex=0.84)
    text(2.30,3.67,"Taxonomic-route heterogeneity",cex=0.75)
    arrows(2.30,5.20,2.30,4.47,length=0.07)

    # Findings under external/source
    rect(4.55,3.40,7.70,4.47,lwd=1.15)
    text(6.125,4.04,"Ecological differentiation beyond trauma",font=2,cex=0.80)
    text(6.125,3.67,"without a universal taxonomic signature",cex=0.72)
    arrows(6.125,5.20,6.125,4.47,length=0.07)

    rect(8.35,3.40,11.40,4.47,lwd=1.15)
    text(9.875,4.04,"No clear source modification",font=2,cex=0.82)
    text(9.875,3.67,"genus-level signals remain exploratory",cex=0.72)
    arrows(9.875,5.20,9.875,4.47,length=0.07)

    # Bottom synthesis
    rect(1.45,0.75,10.55,2.05,lwd=1.4)
    text(
      6,1.55,
      "Shared ecological displacement can emerge through heterogeneous taxonomic routes",
      font=2,
      cex=0.88
    )
    text(
      6,1.12,
      "with limited evidence for a universal taxonomic or infection-source-specific trajectory.",
      cex=0.79
    )

    # Clean, non-crossing arrows to synthesis
    arrows(2.30,3.40,3.50,2.05,length=0.07)
    arrows(6.125,3.40,6.125,2.05,length=0.07)
    arrows(9.875,3.40,8.50,2.05,length=0.07)
  }

  figure_audit <- bind_rows(
    figure_audit,
    make_both(
      "Figure_1_Study_Architecture",
      11.5,7.3,
      2300,1460,
      plot_fig1
    )
  )

  ck("FIGURE 1 COMPLETE")

  # ==========================================================
  # FIGURE 2
  # Source: Step95C frozen longitudinal table.
  # Fix left-label clipping and detached significance stars.
  # ==========================================================

  plot_fig2 <- function() {

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
      mutate(ord = match(Cohort, ord)) %>%
      arrange(ord)

    labels <- c(
      "Core sepsis\nPRJNA691455",
      "External ICU infection\nPRJEB82425",
      "ICU background\nPRJNA851469",
      "ICU background\nPRJNA516701",
      "Non-sepsis surgical control\nPRJNA578267",
      "Intervention support\nPRJNA1166732",
      "Intervention support\nPRJNA430161"
    )

    y <- rev(seq_len(nrow(d)))

    par(mar = c(4.8,11.5,1.2,2.0))

    xmax <- max(
      c(
        d$Primary_effect_dz,
        d$Common_anchor_effect_dz
      ),
      na.rm = TRUE
    )

    xmin <- min(
      c(
        d$Primary_effect_dz,
        d$Common_anchor_effect_dz
      ),
      na.rm = TRUE
    )

    xlim <- c(
      min(-0.75, xmin - 0.10),
      max(2.25, xmax + 0.12)
    )

    plot(
      NA,
      xlim=xlim,
      ylim=c(0.4,nrow(d)+0.6),
      yaxt="n",
      xlab="Paired standardized change in Bray-Curtis displacement (dz)",
      ylab="",
      bty="n"
    )

    abline(v=0,lty=2)

    axis(
      2,
      at=y,
      labels=labels,
      las=1,
      tick=FALSE,
      cex.axis=0.78,
      line=-0.5
    )

    for (i in seq_len(nrow(d))) {
      if (
        !is.na(d$Primary_effect_dz[i]) &&
        !is.na(d$Common_anchor_effect_dz[i])
      ) {
        segments(
          d$Primary_effect_dz[i],y[i],
          d$Common_anchor_effect_dz[i],y[i],
          lty=3
        )
      }
    }

    points(
      d$Primary_effect_dz,
      y,
      pch=16,
      cex=1.0
    )

    points(
      d$Common_anchor_effect_dz,
      y,
      pch=1,
      cex=1.15
    )

    # Put significance next to the primary estimate rather than at plot edge.
    sig_idx <- which(
      !is.na(d$Primary_FDR) &
      d$Primary_FDR < 0.05
    )

    if (length(sig_idx) > 0) {
      text(
        d$Primary_effect_dz[sig_idx] + 0.06,
        y[sig_idx],
        "*",
        cex=1.1,
        adj=c(0,0.5)
      )
    }

    legend(
      "bottomright",
      legend=c(
        "Primary early-to-late",
        "Common-anchor sensitivity",
        "* Primary BH-FDR < 0.05"
      ),
      pch=c(16,1,NA),
      bty="n",
      cex=0.77
    )
  }

  figure_audit <- bind_rows(
    figure_audit,
    make_both(
      "Figure_2_Cross_Cohort_Longitudinal_Displacement",
      11.5,7.4,
      2300,1480,
      plot_fig2
    )
  )

  ck("FIGURE 2 COMPLETE")

  # ==========================================================
  # FIGURE 3
  # Summary-level frozen Step93U evidence.
  # No fragile patient-level file dependency.
  # ==========================================================

  plot_fig3 <- function() {

    par(
      mfrow=c(1,3),
      mar=c(5.0,4.6,2.5,1.1),
      oma=c(0,0,0.3,0)
    )

    # A: ecological displacement time summary
    vals <- c(
      val("Day3_Bray_median"),
      val("Day7_Bray_median")
    )

    bp <- barplot(
      vals,
      names.arg=c("Day 3","Day 7"),
      ylim=c(0,1),
      ylab="Median Bray-Curtis displacement",
      main="A. Progressive ecological displacement"
    )

    segments(
      bp[1], vals[1],
      bp[2], vals[2],
      lty=2
    )

    text(
      mean(bp),
      0.08,
      paste0(
        "paired n=",
        as.integer(val("Day3_to_Day7_paired_n")),
        "\nmean paired change=",
        safe(val("Day3_to_Day7_mean_change"),3)
      ),
      cex=0.76
    )

    # B: heterogeneity metrics
    plot(
      NA,
      xlim=c(0,1),
      ylim=c(0,4),
      axes=FALSE,
      xlab="",
      ylab="",
      main="B. Taxonomic-route heterogeneity"
    )

    text(
      0.05,3.25,
      paste0(
        "Pairwise Euclidean distance\nmean ",
        safe(val("pairwise_signed_Euclidean_mean"),3),
        "; median ",
        safe(val("pairwise_signed_Euclidean_median"),3)
      ),
      adj=c(0,0.5),
      cex=0.84
    )

    text(
      0.05,2.10,
      paste0(
        "Mean cosine similarity\n",
        safe(val("mean_cosine_similarity"),3)
      ),
      adj=c(0,0.5),
      cex=0.84
    )

    text(
      0.05,0.95,
      paste0(
        "Patient pairs with cosine <= 0\n",
        safe(100*val("fraction_cosine_le_0"),1),
        "%"
      ),
      adj=c(0,0.5),
      cex=0.84
    )

    rect(
      0.02,0.45,0.98,3.70
    )

    # C: correlation forest
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

    yy <- c(3,2,1)

    plot(
      rho,yy,
      xlim=c(-0.1,0.9),
      ylim=c(0.5,3.5),
      yaxt="n",
      pch=16,
      xlab="Spearman rho",
      ylab="",
      main="C. Ecological consistency",
      bty="n"
    )

    abline(v=0,lty=2)

    axis(
      2,
      at=yy,
      labels=labs,
      las=1,
      tick=FALSE,
      cex.axis=0.74
    )

    text(
      rho,
      yy,
      labels=c(
        paste0(" p=",safe(val("p_latest_Bray"),3)),
        paste0(
          " p=",safe(val("p_latest_EII"),3),
          "; perm=",safe(val("perm_p_latest_EII"),3)
        ),
        paste0(
          " p=",safe(val("p_Simpson_instability"),3),
          "; perm=",safe(val("perm_p_Simpson_instability"),3)
        )
      ),
      pos=4,
      cex=0.68
    )
  }

  figure_audit <- bind_rows(
    figure_audit,
    make_both(
      "Figure_3_Core_Sepsis_Heterogeneous_Taxonomic_Routes",
      12.0,5.5,
      2400,1100,
      plot_fig3
    )
  )

  ck("FIGURE 3 COMPLETE")

  # ==========================================================
  # FIGURE 4
  # External validation summary from frozen Step93X2/X3 plus
  # source-pack taxonomic alignment.
  # ==========================================================

  ext_value <- function(item) {
    EXT_TAX$value[
      EXT_TAX$item == item
    ][1]
  }

  plot_fig4 <- function() {

    par(
      mfrow=c(1,3),
      mar=c(5.0,4.8,2.5,1.0),
      oma=c(0,0,0.3,0)
    )

    # A healthy-centroid displacement
    vals <- c(
      val("Trauma_healthy_centroid_distance"),
      val("Sepsis_healthy_centroid_distance")
    )

    bp <- barplot(
      vals,
      names.arg=c("Trauma","Sepsis"),
      ylim=c(0,1),
      ylab="Distance from healthy-control centroid",
      main="A. Ecological state displacement"
    )

    text(
      mean(bp),
      0.10,
      paste0(
        "difference +",
        safe(val("Sepsis_minus_Trauma_centroid_difference"),3),
        "\nBH-FDR=",
        safe(val("centroid_difference_FDR"),4)
      ),
      cex=0.74
    )

    # B Aitchison robustness
    plot(
      NA,
      xlim=c(0,1),
      ylim=c(0,4),
      axes=FALSE,
      xlab="",
      ylab="",
      main="B. CLR/Aitchison sensitivity"
    )

    text(
      0.05,3.15,
      paste0(
        "PERMANOVA R2\n",
        safe(val("CLR_Aitchison_R2_min"),3),
        " - ",
        safe(val("CLR_Aitchison_R2_max"),3)
      ),
      adj=c(0,0.5),
      cex=0.85
    )

    text(
      0.05,2.05,
      paste0(
        "PERMANOVA FDR\n",
        format(val("CLR_Aitchison_PERMANOVA_FDR_min"),scientific=TRUE,digits=2),
        " - ",
        format(val("CLR_Aitchison_PERMANOVA_FDR_max"),scientific=TRUE,digits=2)
      ),
      adj=c(0,0.5),
      cex=0.82
    )

    text(
      0.05,0.95,
      paste0(
        "Dispersion FDR\n",
        safe(val("CLR_Aitchison_dispersion_FDR_min"),3),
        " - ",
        safe(val("CLR_Aitchison_dispersion_FDR_max"),3),
        "\n(no dispersion signal)"
      ),
      adj=c(0,0.5),
      cex=0.82
    )

    rect(0.02,0.45,0.98,3.70)

    # C taxonomic alignment
    rho <- c(
      ext_value("Longitudinal-vs-Sepsis-Control Spearman rho"),
      ext_value("Longitudinal-vs-Sepsis-Trauma Spearman rho")
    )

    pp <- c(
      ext_value("Longitudinal-vs-Sepsis-Control p"),
      ext_value("Longitudinal-vs-Sepsis-Trauma p")
    )

    yy <- c(2,1)

    plot(
      rho,yy,
      xlim=c(-0.05,0.42),
      ylim=c(0.5,2.5),
      yaxt="n",
      pch=16,
      xlab="Spearman rho",
      ylab="",
      main="C. Exploratory taxonomic alignment",
      bty="n"
    )

    abline(v=0,lty=2)

    axis(
      2,
      at=yy,
      labels=c(
        "Longitudinal vs\nSepsis-Control",
        "Longitudinal vs\nSepsis-Trauma"
      ),
      las=1,
      tick=FALSE,
      cex.axis=0.72
    )

    text(
      rho,yy,
      labels=paste0(" p=",formatC(pp,format="g",digits=3)),
      pos=4,
      cex=0.72
    )
  }

  figure_audit <- bind_rows(
    figure_audit,
    make_both(
      "Figure_4_External_Validation_Trauma_Comparison",
      12.0,5.5,
      2400,1100,
      plot_fig4
    )
  )

  ck("FIGURE 4 COMPLETE")

  # ==========================================================
  # FIGURE 5
  # Entirely source-pack based.
  # ==========================================================

  primary <- ECO %>%
    filter(
      domain == "PRIMARY_LONGITUDINAL_BRAY",
      estimand == "source_by_personal_time_interaction"
    ) %>%
    slice(1)

  plot_fig5 <- function() {

    par(
      mfrow=c(2,2),
      mar=c(4.8,5.3,2.5,1.0),
      oma=c(0,0,0.3,0)
    )

    # A primary interaction
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
      labels="Pulmonary vs recorded non-pulmonary",
      las=1,
      tick=FALSE,
      cex.axis=0.74
    )
    text(
      primary$effect,0.68,
      paste0(
        "beta=",safe(primary$effect,3),
        "; 95% CI ",
        safe(primary$ci_low,3),
        " to ",
        safe(primary$ci_high,3),
        "; p=",safe(primary$p_value,3)
      ),
      cex=0.72
    )

    # B source-specific slopes
    yy <- c(2,1)

    plot(
      SLOPES$slope_per_day,
      yy,
      xlim=c(
        min(c(SLOPES$ci95_low,0))-0.005,
        max(c(SLOPES$ci95_high,0))+0.01
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
      cex.axis=0.76
    )

    text(
      SLOPES$ci95_high,
      yy,
      labels=paste0(
        " p=",
        safe(SLOPES$p_LRT,3)
      ),
      pos=4,
      cex=0.70
    )

    # C sensitivity forest
    y3 <- rev(seq_len(nrow(SENS)))

    xlo <- min(
      c(SENS$ci_low,0),
      na.rm=TRUE
    ) - 0.005

    xhi <- max(
      c(SENS$ci_high,0),
      na.rm=TRUE
    ) + 0.01

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
      labels=gsub("_"," ",SENS$sensitivity),
      las=1,
      tick=FALSE,
      cex.axis=0.61
    )

    # D baseline beta robustness
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
      cex=0.66
    )
  }

  figure_audit <- bind_rows(
    figure_audit,
    make_both(
      "Figure_5_Infection_Source_Longitudinal_Ecology",
      11.5,8.0,
      2300,1600,
      plot_fig5
    )
  )

  ck("FIGURE 5 COMPLETE")

  # ----------------------------------------------------------
  # Figure audit
  # ----------------------------------------------------------

  write_csv(
    figure_audit,
    file.path(AUDDIR, "02_FIGURE_FILE_AUDIT.csv")
  )

  # ----------------------------------------------------------
  # Captions
  # ----------------------------------------------------------

  captions <- c(
    "Figure 1. Study architecture and evidence hierarchy.",
    "Public gut-microbiome cohorts were organized into a longitudinal critical-illness framework, a nested core repeated-sepsis analysis, an external Control-Trauma-Sepsis ecological-state validation, and an independent infection-source longitudinal extension. The central synthesis emphasizes shared ecological displacement with heterogeneous taxonomic routes and limited evidence for a universal source-specific trajectory.",
    "",
    "Figure 2. Cross-cohort longitudinal ecological displacement.",
    "Primary early-to-late paired standardized Bray-Curtis displacement changes are displayed together with common-anchor sensitivity estimates. Positive values indicate increasing displacement from the patient-specific baseline and negative values indicate re-convergence. Asterisks denote primary BH-FDR < 0.05. The non-sepsis surgical control showed a non-significant recovery-direction primary contrast but significant re-convergence in the common-anchor sensitivity analysis.",
    "",
    "Figure 3. Heterogeneous taxonomic routes accompany ecological displacement in core sepsis.",
    paste0(
      "Median Bray-Curtis displacement increased from approximately ",
      safe(val("Day3_Bray_median"),3),
      " at Day 3 to ",
      safe(val("Day7_Bray_median"),3),
      " at Day 7; the paired Day 3-to-Day 7 comparison included n=",
      as.integer(val("Day3_to_Day7_paired_n")),
      " patients with a mean increase of ",
      safe(val("Day3_to_Day7_mean_change"),3),
      ". Pairwise taxonomic-trajectory distances indicated marked between-patient heterogeneity. Mean trajectory distance correlated with latest EII (rho=",
      safe(val("rho_latest_EII"),3),
      ", permutation p=",
      safe(val("perm_p_latest_EII"),3),
      ") and Simpson instability (rho=",
      safe(val("rho_Simpson_instability"),3),
      ", permutation p=",
      safe(val("perm_p_Simpson_instability"),3),
      "), but not clearly with latest Bray displacement. These associations are internal ecological-consistency evidence and do not establish stable taxonomic subtypes."
    ),
    "",
    "Figure 4. External ecological-state differentiation beyond severe trauma.",
    paste0(
      "Sepsis samples were farther from the healthy-control centroid than severe-trauma samples (",
      safe(val("Sepsis_healthy_centroid_distance"),3),
      " vs ",
      safe(val("Trauma_healthy_centroid_distance"),3),
      "; difference +",
      safe(val("Sepsis_minus_Trauma_centroid_difference"),3),
      "; BH-FDR=",
      safe(val("centroid_difference_FDR"),4),
      "). CLR/Aitchison sensitivity analyses were robust across pseudocounts, with PERMANOVA R2 approximately ",
      safe(val("CLR_Aitchison_R2_min"),3),
      "-",
      safe(val("CLR_Aitchison_R2_max"),3),
      " and no significant dispersion heterogeneity. Cross-analysis taxonomic alignment was evident for the Sepsis-Control contrast but not for Sepsis-Trauma, supporting ecological-state differentiation rather than a reproducible universal genus signature."
    ),
    "",
    "Figure 5. Infection source does not clearly modify longitudinal ecological displacement in CRA002354.",
    paste0(
      "The primary pulmonary-versus-recorded-non-pulmonary source-by-personal-time Bray interaction was beta=",
      safe(primary$effect,3),
      " (95% CI ",
      safe(primary$ci_low,3),
      " to ",
      safe(primary$ci_high,3),
      "; p=",
      safe(primary$p_value,3),
      "). Estimated source-specific slopes were positive in both groups but neither was individually significant. All prespecified key Bray sensitivity analyses retained the same positive interaction direction without statistical significance. The unrestricted baseline PERMANOVA association was accompanied by dispersion heterogeneity and was not reproduced when baseline sampling was restricted to ICU Day 3 or earlier."
    )
  )

  writeLines(
    captions,
    file.path(CAPDIR, "DRAFT_MAIN_FIGURE_CAPTIONS.txt")
  )

  # ----------------------------------------------------------
  # Table placement / next phase
  # ----------------------------------------------------------

  table_plan <- tribble(
    ~item, ~final_location, ~decision,
    "Step95C Evidence Architecture table",
      "Internal provenance / optional Supplementary",
      "Do not use as manuscript Table 1",
    "Step95C Cross-cohort longitudinal table",
      "Main or Supplementary depending journal table limit",
      "Keep; directly reports frozen inferential results",
    "Step95C Compact Supporting Result Summary",
      "Internal drafting aid",
      "Do not submit in current mixed-role format",
    "Conventional cohort-characteristics table",
      "Main Table 1",
      "Build next from frozen metadata"
  )

  write_csv(
    table_plan,
    file.path(CAPDIR, "TABLE_PLACEMENT_PLAN.csv")
  )

  # ----------------------------------------------------------
  # Final QC
  # PDFs are required; PNGs are convenience exports.
  # ----------------------------------------------------------

  pdf_ok <- all(figure_audit$pdf_exists)

  qc <- tibble(
    step95B_complete = TRUE,
    step95C_complete = TRUE,
    figure1_to_5_PDF_complete = pdf_ok,
    figure1_to_5_PNG_complete = all(figure_audit$png_exists),
    no_new_inferential_analysis = TRUE,
    sourcepack_only_core_inputs = TRUE,
    ready_for_visual_review =
      pdf_ok
  )

  write_csv(
    qc,
    file.path(OUT, "STEP95D2_READINESS.csv")
  )

  writeLines(
    c(
      "STEP95D2 VISUAL-QC NOTES",
      "",
      "Step95D was incomplete and is superseded by Step95D2.",
      "",
      "Figure 1:",
      "- Corrected evidence architecture.",
      "- Core repeated sepsis is nested within the longitudinal framework.",
      "- External validation and infection-source extension are parallel evidence branches.",
      "- Removed crossing arrows.",
      "",
      "Figure 2:",
      "- Increased left margin and wrapped cohort labels to prevent clipping.",
      "- Significance markers are positioned beside primary estimates.",
      "",
      "Figure 3:",
      "- Uses frozen summary-level Step93U evidence rather than fragile patient-level plotting files.",
      "- Explicitly emphasizes ecological displacement + heterogeneous taxonomic routes.",
      "",
      "Figure 4:",
      "- Uses frozen external ecological-state and CLR/Aitchison sensitivity evidence.",
      "- Taxonomic alignment remains explicitly exploratory.",
      "",
      "Figure 5:",
      "- Uses only Step95B frozen ecology/sensitivity registries.",
      "- Directly communicates the null source-modification conclusion and baseline robustness issue.",
      "",
      "PDF is the authoritative figure format. PNG generation is optional convenience output.",
      "No new inferential analysis was performed."
    ),
    file.path(CAPDIR, "STEP95D2_VISUAL_QC_REPORT.txt")
  )

  if (!pdf_ok) {
    fatal("One or more final Figure 1-5 PDFs were not created.")
  }

  writeLines(
    c(
      paste0("Completed: ", Sys.time()),
      "Figure 1-5 PDF complete: TRUE",
      paste0(
        "Figure 1-5 PNG complete: ",
        all(figure_audit$png_exists)
      ),
      "No new inferential analysis: TRUE",
      "STEP95D2 COMPLETE"
    ),
    file.path(OUT, "_STEP95D2_COMPLETE.ok")
  )

  ck("STEP95D2 COMPLETE")

}, error = function(e) {

  if (!file.exists(ERR)) {
    writeLines(
      c(
        paste0("STEP95D2 FAILED: ", Sys.time()),
        conditionMessage(e)
      ),
      ERR
    )
  }

  ck(
    paste0(
      "STEP95D2 FAILED - ",
      conditionMessage(e)
    )
  )

  stop(e)
})
