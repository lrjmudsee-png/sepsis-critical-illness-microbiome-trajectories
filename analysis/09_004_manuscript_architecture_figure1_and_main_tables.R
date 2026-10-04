# ============================================================
# Sepsis V2 - Step95C
# MANUSCRIPT ARCHITECTURE: FIGURE 1 + MAIN TABLE FOUNDATION
#
# Input:
# Step95B manuscript source pack ONLY
#
# Purpose:
# 1) create Figure 1 study/evidence architecture;
# 2) create manuscript-ready Main Table 1 evidence architecture;
# 3) create manuscript-ready Main Table 2 longitudinal cohort results;
# 4) create a compact supporting-result summary for later Results drafting;
# 5) freeze a figure/table production checklist.
#
# NO new inferential analysis.
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

SRC <- file.path(
  RESULTS,
  "V2_35B_MANUSCRIPT_SOURCE_PACK_ASSEMBLY"
)

OUT <- file.path(
  RESULTS,
  "V2_35C_MANUSCRIPT_ARCHITECTURE_FIGURE1_AND_MAIN_TABLES"
)

dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(file.path(SRC, "_STEP95B_COMPLETE.ok"))) {
  stop("Step95B source pack is not complete.")
}

# ------------------------------------------------------------
# 1. Load source-pack registries
# ------------------------------------------------------------

evidence <- read_csv(
  file.path(
    SRC,
    "05_PROVENANCE",
    "STEP95A2_06_FINAL_GLOBAL_V2_EVIDENCE_REGISTRY.csv"
  ),
  show_col_types = FALSE
)

claims <- read_csv(
  file.path(
    SRC,
    "05_PROVENANCE",
    "STEP95A2_07_FINAL_MANUSCRIPT_CLAIM_HIERARCHY.csv"
  ),
  show_col_types = FALSE
)

figmap <- read_csv(
  file.path(
    SRC,
    "05_PROVENANCE",
    "STEP95A2_08_FINAL_MANUSCRIPT_FIGURE_MAP.csv"
  ),
  show_col_types = FALSE
)

results_order <- read_csv(
  file.path(
    SRC,
    "05_PROVENANCE",
    "STEP95A2_09_RESULTS_SECTION_ORDER.csv"
  ),
  show_col_types = FALSE
)

long <- read_csv(
  file.path(
    SRC,
    "01_MAIN_TABLE_SOURCES",
    "MAIN_T2_CROSS_COHORT_LONGITUDINAL.csv"
  ),
  show_col_types = FALSE
)

core <- read_csv(
  file.path(
    SRC,
    "01_MAIN_TABLE_SOURCES",
    "MAIN_T3_CORE_TRAJECTORY_ENDPOINTS.csv"
  ),
  show_col_types = FALSE
)

external_tax <- read_csv(
  file.path(
    SRC,
    "01_MAIN_TABLE_SOURCES",
    "MAIN_T4_EXTERNAL_VALIDATION.csv"
  ),
  show_col_types = FALSE
)

source_eco <- read_csv(
  file.path(
    SRC,
    "01_MAIN_TABLE_SOURCES",
    "MAIN_T5_INFECTION_SOURCE_ECOLOGY.csv"
  ),
  show_col_types = FALSE
)

source_genus <- read_csv(
  file.path(
    SRC,
    "02_SUPPLEMENTARY_TABLE_SOURCES",
    "SUPP_S10_GENUS_BRANCH_FREEZE.csv"
  ),
  show_col_types = FALSE
)

source_wording <- paste(
  readLines(
    file.path(
      SRC,
      "04_FROZEN_WORDING",
      "SOURCE_ECOLOGY_WORDING.txt"
    ),
    warn = FALSE,
    encoding = "UTF-8"
  ),
  collapse = "\n"
)

# ------------------------------------------------------------
# 2. Main Table 1 - evidence architecture
# ------------------------------------------------------------

table1 <- evidence %>%
  transmute(
    Evidence_line = evidence_line,
    Frozen_branch = source_branch,
    Analysis_role = role,
    Manuscript_weight = manuscript_weight,
    Source_anchored = source_anchored,
    Frozen_conclusion = frozen_conclusion
  )

write_csv(
  table1,
  file.path(
    OUT,
    "Main_Table_1_Evidence_Architecture.csv"
  )
)

# ------------------------------------------------------------
# 3. Main Table 2 - cross-cohort longitudinal results
# ------------------------------------------------------------

table2 <- long %>%
  transmute(
    Cohort = project,
    Cohort_role = cohort_role,
    Primary_pairs = primary_n_pairs,
    Primary_effect_dz = primary_dz,
    Primary_p = primary_p,
    Primary_FDR = primary_fdr,
    Primary_direction = primary_direction,
    Common_anchor_pairs = common_anchor_n_pairs,
    Common_anchor_effect_dz = common_anchor_dz,
    Common_anchor_p = common_anchor_p,
    Common_anchor_FDR = common_anchor_fdr,
    Direction_robust = direction_robust,
    Evidence_class = frozen_evidence_class
  ) %>%
  arrange(
    factor(
      Cohort_role,
      levels = c(
        "CORE_SEPSIS_LONGITUDINAL",
        "ICU_INFECTION_LONGITUDINAL_EXTERNAL",
        "ICU_BACKGROUND_LONGITUDINAL",
        "NONSEPSIS_LONGITUDINAL_CONTROL",
        "INTERVENTION_LONGITUDINAL_SUPPORT"
      )
    )
  )

write_csv(
  table2,
  file.path(
    OUT,
    "Main_Table_2_Cross_Cohort_Longitudinal_Results.csv"
  )
)

# ------------------------------------------------------------
# 4. Compact supporting-result summary
# ------------------------------------------------------------

get_metric <- function(name) {
  z <- core %>%
    filter(metric == name)
  if (nrow(z) == 0) return(NA_character_)
  sprintf(
    "%.3f [%.3f, %.3f]",
    z$median[1],
    z$q1[1],
    z$q3[1]
  )
}

# Parse CRA 95% CI from frozen authoritative wording.
ci_match <- str_match(
  source_wording,
  regex(
    "95% CI\\s*(-?[0-9.]+)\\s*to\\s*(-?[0-9.]+)",
    ignore_case = TRUE
  )
)

cra_ci <- if (
  nrow(ci_match) > 0 &&
  !is.na(ci_match[1,2])
) {
  paste0(
    "[",
    ci_match[1,2],
    ", ",
    ci_match[1,3],
    "]"
  )
} else {
  NA_character_
}

support <- tibble(
  Evidence_component = c(
    "Core-sepsis latest Bray displacement",
    "Core-sepsis trajectory clustering",
    "External exploratory taxonomic contrast",
    "CRA002354 source-by-time Bray interaction",
    "CRA002354 primary alpha source-by-time family",
    "CRA002354 genus baseline family",
    "CRA002354 genus longitudinal family"
  ),
  Key_result = c(
    paste0(
      "Median [IQR] = ",
      get_metric("latest_Bray")
    ),
    "Exploratory k=2 split = 9 vs 1; not a validated subtype solution",
    paste0(
      "Sepsis-vs-Trauma FDR-significant genera = ",
      external_tax$value[
        external_tax$item ==
          "Sepsis-vs-Trauma FDR-significant genera"
      ][1],
      "; exploratory taxonomic result only"
    ),
    paste0(
      "beta = ",
      sprintf(
        "%.3f",
        source_eco$primary_Bray_interaction_beta[1]
      ),
      "; 95% CI ",
      cra_ci,
      "; p = ",
      sprintf(
        "%.3f",
        source_eco$primary_Bray_interaction_p[1]
      )
    ),
    paste0(
      source_eco$primary_alpha_FDR_significant[1],
      " of 3 FDR-significant"
    ),
    paste0(
      source_genus$primary_baseline_FDR_hits[1],
      " FDR hits; ",
      source_genus$baseline_hits_FDR_robust_at_both[1],
      " robust at both >=70% and >=80% coverage thresholds"
    ),
    paste0(
      source_genus$primary_longitudinal_FDR_hits[1],
      " primary FDR hits"
    )
  ),
  Reporting_role = c(
    "MAIN_INTERPRETIVE",
    "GUARDRAIL",
    "SUPPLEMENTARY_EXPLORATORY",
    "SUPPORTING_MAIN_RESULT",
    "SUPPORTING_MAIN_RESULT",
    "SUPPLEMENTARY_EXPLORATORY",
    "SUPPLEMENTARY_EXPLORATORY"
  )
)

write_csv(
  support,
  file.path(
    OUT,
    "Main_Table_3_Compact_Supporting_Result_Summary.csv"
  )
)

# ------------------------------------------------------------
# 5. Figure 1 - study architecture / evidence hierarchy
# ------------------------------------------------------------

draw_box <- function(
  xleft, ybottom, xright, ytop,
  label,
  cex = 0.85,
  font = 1
) {
  rect(
    xleft,
    ybottom,
    xright,
    ytop,
    lwd = 1.2
  )
  text(
    (xleft + xright)/2,
    (ybottom + ytop)/2,
    label,
    cex = cex,
    font = font
  )
}

draw_arrow <- function(x0,y0,x1,y1) {
  arrows(
    x0,y0,x1,y1,
    length = 0.08,
    lwd = 1.0
  )
}

make_fig1 <- function(filename, type = c("pdf","png")) {

  type <- match.arg(type)

  if (type == "pdf") {
    pdf(filename, width = 12, height = 7.5)
  } else {
    png(
      filename,
      width = 2400,
      height = 1500,
      res = 200
    )
  }

  par(
    mar = c(0.7,0.7,1.5,0.7)
  )

  plot(
    NA,
    xlim = c(0,12),
    ylim = c(0,8),
    xaxs = "i",
    yaxs = "i",
    axes = FALSE,
    xlab = "",
    ylab = ""
  )

  title(
    "Figure 1. V2 study architecture and evidence hierarchy",
    cex.main = 1.2
  )

  # Top
  draw_box(
    4.1, 6.9, 7.9, 7.7,
    "Public microbiome cohorts\nLongitudinal + external + infection-source datasets",
    cex = 0.9,
    font = 2
  )

  # Second level
  draw_box(
    0.4, 5.0, 3.2, 6.1,
    "Longitudinal evidence\n7 cohorts\nStep93W",
    cex = 0.85
  )

  draw_box(
    3.45, 5.0, 6.15, 6.1,
    "Core sepsis trajectories\nPRJNA691455\nStep93U",
    cex = 0.85
  )

  draw_box(
    6.4, 5.0, 9.1, 6.1,
    "External state validation\nControl–Trauma–Sepsis\nStep93X4",
    cex = 0.85
  )

  draw_box(
    9.35, 5.0, 11.7, 6.1,
    "Infection-source branch\nCRA002354\nStep94B2D/B3C",
    cex = 0.82
  )

  draw_arrow(6.0,6.9,1.8,6.1)
  draw_arrow(6.0,6.9,4.8,6.1)
  draw_arrow(6.0,6.9,7.75,6.1)
  draw_arrow(6.0,6.9,10.5,6.1)

  # Third level conclusions
  draw_box(
    0.4, 2.7, 3.2, 4.1,
    "Primary ecological result\nProgressive within-patient\necological displacement\nacross several cohorts",
    cex = 0.8,
    font = 2
  )

  draw_box(
    3.45, 2.7, 6.15, 4.1,
    "Primary interpretation\nHeterogeneous taxonomic routes\nwithout validated stable\ntrajectory subtypes",
    cex = 0.8,
    font = 2
  )

  draw_box(
    6.4, 2.7, 9.1, 4.1,
    "Supportive external evidence\nEcological state differentiation\nbeyond severe trauma;\nno universal taxonomic signature",
    cex = 0.78
  )

  draw_box(
    9.35, 2.7, 11.7, 4.1,
    "Supportive source evidence\nNo clear source modification\nof longitudinal displacement;\ngenus signals exploratory",
    cex = 0.77
  )

  draw_arrow(1.8,5.0,1.8,4.1)
  draw_arrow(4.8,5.0,4.8,4.1)
  draw_arrow(7.75,5.0,7.75,4.1)
  draw_arrow(10.5,5.0,10.5,4.1)

  # Bottom synthesis
  draw_box(
    2.0, 0.6, 10.0, 1.7,
    "Central manuscript synthesis:\nShared ecological displacement can emerge through heterogeneous taxonomic routes,\nwith limited evidence for a universal taxonomic or infection-source-specific trajectory.",
    cex = 0.9,
    font = 2
  )

  draw_arrow(1.8,2.7,3.4,1.7)
  draw_arrow(4.8,2.7,5.1,1.7)
  draw_arrow(7.75,2.7,6.9,1.7)
  draw_arrow(10.5,2.7,8.6,1.7)

  dev.off()
}

make_fig1(
  file.path(
    OUT,
    "Figure_1_V2_Study_Architecture.pdf"
  ),
  "pdf"
)

make_fig1(
  file.path(
    OUT,
    "Figure_1_V2_Study_Architecture.png"
  ),
  "png"
)

# ------------------------------------------------------------
# 6. Figure 2–5 source-production checklist
# ------------------------------------------------------------

fig_manifest <- read_csv(
  file.path(
    SRC,
    "03_FIGURE_SOURCES",
    "00_FIGURE_SOURCE_MANIFEST.csv"
  ),
  show_col_types = FALSE
)

fig_check <- figmap %>%
  left_join(
    fig_manifest %>%
      filter(copied) %>%
      group_by(figure_role) %>%
      summarise(
        candidate_files = n(),
        .groups = "drop"
      ),
    by = c(
      "source_branch" =
        "figure_role"
    )
  )

# Explicit mapping because source_branch labels differ from manifest roles.
fig_check <- tibble(
  Figure = c(
    "Figure 1",
    "Figure 2",
    "Figure 3",
    "Figure 4",
    "Figure 5"
  ),
  Content = c(
    "Study architecture and evidence hierarchy",
    "Cross-cohort longitudinal ecological displacement",
    "Core-sepsis heterogeneous taxonomic trajectories",
    "External Control-Trauma-Sepsis ecological validation",
    "CRA002354 infection-source longitudinal ecology"
  ),
  Frozen_source = c(
    "STEP95A2/95B",
    "STEP93W",
    "STEP93U",
    "STEP93X4",
    "STEP94B2D"
  ),
  Candidate_status = c(
    "NEW_FIGURE1_CREATED",
    ifelse(
      any(
        fig_manifest$figure_role ==
          "FIGURE2_STEP93W" &
        fig_manifest$copied
      ),
      "SOURCE_PRESENT",
      "SOURCE_MISSING"
    ),
    ifelse(
      any(
        fig_manifest$figure_role ==
          "FIGURE3_STEP93U" &
        fig_manifest$copied
      ),
      "SOURCE_PRESENT",
      "SOURCE_MISSING"
    ),
    ifelse(
      any(
        fig_manifest$figure_role ==
          "FIGURE4_STEP93X4" &
        fig_manifest$copied
      ),
      "SOURCE_PRESENT",
      "SOURCE_MISSING"
    ),
    ifelse(
      any(
        fig_manifest$figure_role ==
          "FIGURE5_STEP94B2D" &
        fig_manifest$copied
      ),
      "SOURCE_PRESENT",
      "SOURCE_MISSING"
    )
  ),
  Next_action = c(
    "Visual QC / caption drafting",
    "Panel-level visual QC and manuscript formatting",
    "Panel-level visual QC and manuscript formatting",
    "Panel-level visual QC and manuscript formatting",
    "Panel-level visual QC and manuscript formatting"
  )
)

write_csv(
  fig_check,
  file.path(
    OUT,
    "Figure_Production_Checklist.csv"
  )
)

# ------------------------------------------------------------
# 7. Results drafting blueprint
# ------------------------------------------------------------

draft_blueprint <- results_order %>%
  transmute(
    Section = section_order,
    Results_heading = results_section,
    Primary_source = main_source,
    Supplementary_link = supplementary_link,
    Draft_status = case_when(
      section_order == 1 ~
        "READY_AFTER_FIGURE1_TABLE1",
      TRUE ~
        "READY_AFTER_FIGURE_VISUAL_QC"
    )
  )

write_csv(
  draft_blueprint,
  file.path(
    OUT,
    "Results_Drafting_Blueprint.csv"
  )
)

# ------------------------------------------------------------
# 8. Final QC
# ------------------------------------------------------------

qc <- tibble(
  step95B_complete = TRUE,
  figure1_pdf_created =
    file.exists(
      file.path(
        OUT,
        "Figure_1_V2_Study_Architecture.pdf"
      )
    ),
  figure1_png_created =
    file.exists(
      file.path(
        OUT,
        "Figure_1_V2_Study_Architecture.png"
      )
    ),
  main_table1_created =
    file.exists(
      file.path(
        OUT,
        "Main_Table_1_Evidence_Architecture.csv"
      )
    ),
  main_table2_created =
    file.exists(
      file.path(
        OUT,
        "Main_Table_2_Cross_Cohort_Longitudinal_Results.csv"
      )
    ),
  figure2_to_5_sources_present =
    all(
      fig_check$Candidate_status[
        fig_check$Figure !=
          "Figure 1"
      ] ==
      "SOURCE_PRESENT"
    )
)

qc <- qc %>%
  mutate(
    ready_for_figure_visual_qc_and_results_drafting =
      figure1_pdf_created &
      figure1_png_created &
      main_table1_created &
      main_table2_created &
      figure2_to_5_sources_present
  )

write_csv(
  qc,
  file.path(
    OUT,
    "STEP95C_READINESS.csv"
  )
)

writeLines(
  c(
    "STEP95C MANUSCRIPT ARCHITECTURE / FIGURE 1 / MAIN TABLE FOUNDATION",
    "",
    paste0(
      "Figure 1 created: ",
      qc$figure1_pdf_created
    ),
    paste0(
      "Main Table 1 created: ",
      qc$main_table1_created
    ),
    paste0(
      "Main Table 2 created: ",
      qc$main_table2_created
    ),
    paste0(
      "Figure 2-5 frozen sources present: ",
      qc$figure2_to_5_sources_present
    ),
    paste0(
      "Ready for visual QC + Results drafting: ",
      qc$ready_for_figure_visual_qc_and_results_drafting
    ),
    "",
    "No new inferential statistics were performed.",
    "The external taxonomic value in Main Table 3 is explicitly exploratory and must not substitute for the frozen ecological external-state conclusion."
  ),
  file.path(
    OUT,
    "STEP95C_INTERPRETATION.txt"
  )
)

if (
  !isTRUE(
    qc$ready_for_figure_visual_qc_and_results_drafting
  )
) {
  stop(
    "Step95C readiness failed."
  )
}

writeLines(
  c(
    paste0("Completed: ", Sys.time()),
    "Ready for visual QC + Results drafting: TRUE",
    "STEP95C COMPLETE"
  ),
  file.path(
    OUT,
    "_STEP95C_COMPLETE.ok"
  )
)

cat("STEP95C COMPLETE\n")
