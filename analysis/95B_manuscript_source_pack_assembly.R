# ============================================================
# Sepsis V2 - Step95B
# MANUSCRIPT SOURCE PACK ASSEMBLY
#
# Purpose:
# Create a single, version-controlled source pack for manuscript production
# from the FINAL Step95A2 global synthesis and frozen evidence branches.
#
# NO new inferential analysis.
#
# Outputs:
# - authoritative data-source manifest
# - main-table source pack
# - supplementary-table source pack
# - figure-source manifest and copied candidate PDFs/CSVs
# - frozen wording / guardrails
# - provenance audit
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

A2 <- file.path(
  RESULTS,
  "V2_35A2_GLOBAL_EVIDENCE_SYNTHESIS_SOURCE_ANCHOR_FIX"
)

OUT <- file.path(
  RESULTS,
  "V2_35B_MANUSCRIPT_SOURCE_PACK_ASSEMBLY"
)

MAIN <- file.path(OUT, "01_MAIN_TABLE_SOURCES")
SUPP <- file.path(OUT, "02_SUPPLEMENTARY_TABLE_SOURCES")
FIG <- file.path(OUT, "03_FIGURE_SOURCES")
WORD <- file.path(OUT, "04_FROZEN_WORDING")
PROV <- file.path(OUT, "05_PROVENANCE")

for (d in c(OUT, MAIN, SUPP, FIG, WORD, PROV)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

if (!file.exists(file.path(A2, "_STEP95A2_COMPLETE.ok"))) {
  stop("Step95A2 is not complete.")
}

source_registry <- read_csv(
  file.path(A2, "05_FINAL_AUTHORITATIVE_SOURCE_FILE_REGISTRY.csv"),
  show_col_types = FALSE
)

evidence_registry <- read_csv(
  file.path(A2, "06_FINAL_GLOBAL_V2_EVIDENCE_REGISTRY.csv"),
  show_col_types = FALSE
)

claims <- read_csv(
  file.path(A2, "07_FINAL_MANUSCRIPT_CLAIM_HIERARCHY.csv"),
  show_col_types = FALSE
)

figmap <- read_csv(
  file.path(A2, "08_FINAL_MANUSCRIPT_FIGURE_MAP.csv"),
  show_col_types = FALSE
)

results_order <- read_csv(
  file.path(A2, "09_RESULTS_SECTION_ORDER.csv"),
  show_col_types = FALSE
)

# ------------------------------------------------------------
# 1. Freeze the Step95A2 global synthesis itself
# ------------------------------------------------------------

a2_files <- c(
  "04_FINAL_FROZEN_BRANCH_AUDIT.csv",
  "05_FINAL_AUTHORITATIVE_SOURCE_FILE_REGISTRY.csv",
  "06_FINAL_GLOBAL_V2_EVIDENCE_REGISTRY.csv",
  "07_FINAL_MANUSCRIPT_CLAIM_HIERARCHY.csv",
  "08_FINAL_MANUSCRIPT_FIGURE_MAP.csv",
  "09_RESULTS_SECTION_ORDER.csv",
  "10_FINAL_GLOBAL_MANUSCRIPT_GUARDRAILS.txt",
  "11_FINAL_GLOBAL_V2_NARRATIVE_AND_NEXT_PHASE.txt",
  "12_FINAL_GLOBAL_SOURCE_ANCHOR_AUDIT.csv",
  "_STEP95A2_COMPLETE.ok"
)

a2_audit <- tibble(
  file = a2_files,
  source = file.path(A2, a2_files),
  copied = FALSE
)

for (i in seq_along(a2_files)) {
  src <- file.path(A2, a2_files[i])
  if (file.exists(src)) {
    file.copy(
      src,
      file.path(PROV, paste0("STEP95A2_", a2_files[i])),
      overwrite = TRUE
    )
    a2_audit$copied[i] <- TRUE
  }
}

write_csv(
  a2_audit,
  file.path(PROV, "01_STEP95A2_COPY_AUDIT.csv")
)

if (!all(a2_audit$copied)) {
  stop("Failed to copy one or more Step95A2 authoritative files.")
}

# ------------------------------------------------------------
# 2. Branch directories and curated source files
# ------------------------------------------------------------

D93W <- "E:/sepsis_project/results/V2_33W_LONGITUDINAL_EVIDENCE_FREEZE"
D93U <- "E:/sepsis_project/results/V2_33U_TRAJECTORY_ECOLOGY_ROBUSTNESS"
D93X4 <- "E:/sepsis_project/results/V2_33X4_EXTERNAL_VALIDATION_EVIDENCE_FREEZE"
D94B2D <- "E:/sepsis_project/results/V2_34B2D_INFECTION_SOURCE_EVIDENCE_FREEZE"
D94B3C <- "E:/sepsis_project/results/V2_34B3C_CRA002354_GENUS_EVIDENCE_FREEZE"

branches <- tibble(
  branch = c("STEP93W","STEP93U","STEP93X4","STEP94B2D","STEP94B3C"),
  directory = c(D93W,D93U,D93X4,D94B2D,D94B3C)
) %>%
  mutate(exists = dir.exists(directory))

if (!all(branches$exists)) {
  stop("One or more frozen branch directories are missing.")
}

# helper
copy_named <- function(src, dst_dir, dst_name) {
  if (!file.exists(src)) return(FALSE)
  file.copy(src, file.path(dst_dir, dst_name), overwrite = TRUE)
}

# ------------------------------------------------------------
# 3. Main table sources
# ------------------------------------------------------------

main_candidates <- list(
  MAIN_T2_CROSS_COHORT_LONGITUDINAL =
    file.path(D93W, "01_FROZEN_longitudinal_evidence_matrix.csv"),

  MAIN_T3_CORE_TRAJECTORY_ENDPOINTS =
    file.path(D93U, "02_ecological_endpoint_summary.csv"),

  MAIN_T4_EXTERNAL_VALIDATION =
    file.path(D93X4, "03_FROZEN_taxonomic_external_summary.csv"),

  MAIN_T5_INFECTION_SOURCE_ECOLOGY =
    file.path(D94B2D, "05_FINAL_INFECTION_SOURCE_BRANCH_FREEZE.csv")
)

main_audit <- bind_rows(
  lapply(names(main_candidates), function(nm) {
    src <- main_candidates[[nm]]
    ext <- tools::file_ext(src)
    dst <- paste0(nm, ".", ext)
    ok <- copy_named(src, MAIN, dst)
    tibble(
      source_role = nm,
      source = src,
      copied_as = ifelse(ok, file.path(MAIN, dst), NA_character_),
      copied = ok
    )
  })
)

write_csv(
  main_audit,
  file.path(MAIN, "00_MAIN_TABLE_SOURCE_AUDIT.csv")
)

# ------------------------------------------------------------
# 4. Supplementary table sources
# ------------------------------------------------------------

supp_candidates <- list(
  SUPP_S1_LONGITUDINAL_EVIDENCE =
    file.path(D93W, "01_FROZEN_longitudinal_evidence_matrix.csv"),

  SUPP_S2_CORE_TRAJECTORY_ENDPOINTS =
    file.path(D93U, "02_ecological_endpoint_summary.csv"),

  SUPP_S3_EXTERNAL_VALIDATION =
    file.path(D93X4, "03_FROZEN_taxonomic_external_summary.csv"),

  SUPP_S4_SOURCE_ECOLOGY_REGISTRY =
    file.path(D94B2D, "01_FINAL_INFECTION_SOURCE_EVIDENCE_REGISTRY.csv"),

  SUPP_S5_SOURCE_BRAY_SENSITIVITY =
    file.path(D94B2D, "02_FINAL_BRAY_SENSITIVITY_ROBUSTNESS_REGISTRY.csv"),

  SUPP_S6_SOURCE_SLOPES =
    file.path(D94B2D, "03_SOURCE_SPECIFIC_SLOPE_REGISTRY.csv"),

  SUPP_S7_SOURCE_BASELINE_BETA =
    file.path(D94B2D, "04_BASELINE_BETA_ROBUSTNESS_FREEZE.csv"),

  SUPP_S8_GENUS_BASELINE_HITS =
    file.path(D94B3C, "01_FINAL_BASELINE_GENUS_HIT_REGISTRY.csv"),

  SUPP_S9_GENUS_LONGITUDINAL_SENSITIVITY_ONLY =
    file.path(D94B3C, "02_LONGITUDINAL_COVERAGE_SENSITIVITY_ONLY_SIGNALS.csv"),

  SUPP_S10_GENUS_BRANCH_FREEZE =
    file.path(D94B3C, "03_FINAL_GENUS_BRANCH_EVIDENCE_FREEZE.csv"),

  SUPP_S11_GENUS_REPORTING_RECOMMENDATION =
    file.path(D94B3C, "04_GENUS_REPORTING_RECOMMENDATION.csv")
)

supp_audit <- bind_rows(
  lapply(names(supp_candidates), function(nm) {
    src <- supp_candidates[[nm]]
    ext <- tools::file_ext(src)
    dst <- paste0(nm, ".", ext)
    ok <- copy_named(src, SUPP, dst)
    tibble(
      source_role = nm,
      source = src,
      copied_as = ifelse(ok, file.path(SUPP, dst), NA_character_),
      copied = ok
    )
  })
)

write_csv(
  supp_audit,
  file.path(SUPP, "00_SUPPLEMENTARY_TABLE_SOURCE_AUDIT.csv")
)

# ------------------------------------------------------------
# 5. Figure candidate discovery and copying
# ------------------------------------------------------------

find_candidates <- function(dir, regex_pattern) {
  if (!dir.exists(dir)) return(character())
  fs <- list.files(
    dir,
    recursive = TRUE,
    full.names = TRUE
  )
  fs[
    str_detect(
      tolower(basename(fs)),
      regex(regex_pattern, ignore_case = TRUE)
    )
  ]
}

fig_candidates <- list(
  FIGURE2_STEP93W =
    find_candidates(
      D93W,
      "\\.pdf$|figure.*\\.(png|pdf)$|plot.*\\.(png|pdf)$"
    ),

  FIGURE3_STEP93U =
    find_candidates(
      D93U,
      "\\.pdf$|figure.*\\.(png|pdf)$|plot.*\\.(png|pdf)$"
    ),

  FIGURE4_STEP93X4 =
    find_candidates(
      D93X4,
      "\\.pdf$|figure.*\\.(png|pdf)$|plot.*\\.(png|pdf)$"
    ),

  FIGURE5_STEP94B2D =
    find_candidates(
      D94B2D,
      "\\.pdf$|figure.*\\.(png|pdf)$|plot.*\\.(png|pdf)$"
    ),

  SUPPLEMENTARY_STEP94B3C =
    find_candidates(
      D94B3C,
      "\\.pdf$|figure.*\\.(png|pdf)$|plot.*\\.(png|pdf)$"
    )
)

fig_manifest <- tibble()

for (role in names(fig_candidates)) {

  fs <- fig_candidates[[role]]

  if (length(fs) == 0) {
    fig_manifest <- bind_rows(
      fig_manifest,
      tibble(
        figure_role = role,
        original_file = NA_character_,
        copied_file = NA_character_,
        copied = FALSE
      )
    )
    next
  }

  for (j in seq_along(fs)) {

    src <- fs[j]
    safe_role <- gsub("[^A-Za-z0-9_]+", "_", role)
    dst_name <- paste0(
      safe_role,
      "_",
      sprintf("%02d", j),
      "_",
      basename(src)
    )

    ok <- copy_named(
      src,
      FIG,
      dst_name
    )

    fig_manifest <- bind_rows(
      fig_manifest,
      tibble(
        figure_role = role,
        original_file = src,
        copied_file = ifelse(
          ok,
          file.path(FIG, dst_name),
          NA_character_
        ),
        copied = ok
      )
    )
  }
}

write_csv(
  fig_manifest,
  file.path(FIG, "00_FIGURE_SOURCE_MANIFEST.csv")
)

# ------------------------------------------------------------
# 6. Frozen wording / guardrails
# ------------------------------------------------------------

wording_candidates <- list(
  GLOBAL_GUARDRAILS =
    file.path(A2, "10_FINAL_GLOBAL_MANUSCRIPT_GUARDRAILS.txt"),

  GLOBAL_NARRATIVE =
    file.path(A2, "11_FINAL_GLOBAL_V2_NARRATIVE_AND_NEXT_PHASE.txt"),

  LONGITUDINAL_METHODS =
    file.path(D93W, "07_FROZEN_longitudinal_methods_wording.txt"),

  CORE_TRAJECTORY_INTERPRETATION =
    file.path(D93U, "10_STEP93U_INTERPRETATION.txt"),

  EXTERNAL_METHODS =
    file.path(D93X4, "08_FROZEN_external_methods_wording.txt"),

  SOURCE_ECOLOGY_WORDING =
    file.path(D94B2D, "06_FROZEN_MANUSCRIPT_WORDING_AND_GUARDRAILS.txt"),

  SOURCE_GENUS_WORDING =
    file.path(D94B3C, "05_FROZEN_GENUS_MANUSCRIPT_WORDING_AND_GUARDRAILS.txt")
)

word_audit <- bind_rows(
  lapply(names(wording_candidates), function(nm) {
    src <- wording_candidates[[nm]]
    dst <- paste0(nm, ".txt")
    ok <- copy_named(src, WORD, dst)
    tibble(
      wording_role = nm,
      source = src,
      copied_as = ifelse(ok, file.path(WORD, dst), NA_character_),
      copied = ok
    )
  })
)

write_csv(
  word_audit,
  file.path(WORD, "00_FROZEN_WORDING_SOURCE_AUDIT.csv")
)

# ------------------------------------------------------------
# 7. Production blueprint
# ------------------------------------------------------------

blueprint <- tibble(
  deliverable = c(
    "Figure 1",
    "Figure 2",
    "Figure 3",
    "Figure 4",
    "Figure 5",
    "Main Table 1",
    "Main Table 2",
    "Supplementary Tables",
    "Results",
    "Discussion"
  ),
  content = c(
    "Study architecture, cohort roles, and analysis flow",
    "Cross-cohort longitudinal ecological displacement",
    "Core-sepsis heterogeneous taxonomic trajectories",
    "External Control-Trauma-Sepsis ecological validation",
    "CRA002354 infection-source longitudinal ecology",
    "Cohort roles and analysis tiers",
    "Key frozen inferential results across evidence lines",
    "Full models, sensitivities, exploratory genus results, EII",
    "Write in Step95A2 section order",
    "Shared displacement + heterogeneous routes + external differentiation + limits of source specificity"
  ),
  primary_source = c(
    "Step95A2 registry + metadata",
    "STEP93W",
    "STEP93U",
    "STEP93X4",
    "STEP94B2D",
    "Step95A2 + Step93W",
    "Step95A2 evidence registry + branch freezes",
    "01/02 source folders in this Step95B pack",
    "Step95A2 09_RESULTS_SECTION_ORDER.csv",
    "Step95A2 narrative + global guardrails"
  ),
  production_status = c(
    "NEXT",
    "READY",
    "READY",
    "READY",
    "READY",
    "NEXT",
    "NEXT",
    "READY",
    "AFTER_FIGURES_TABLES",
    "AFTER_RESULTS"
  )
)

write_csv(
  blueprint,
  file.path(
    OUT,
    "06_MANUSCRIPT_PRODUCTION_BLUEPRINT.csv"
  )
)

# ------------------------------------------------------------
# 8. Global source manifest
# ------------------------------------------------------------

manifest <- bind_rows(
  main_audit %>%
    transmute(
      category = "MAIN_TABLE_SOURCE",
      role = source_role,
      source = source,
      copied = copied
    ),
  supp_audit %>%
    transmute(
      category = "SUPPLEMENTARY_TABLE_SOURCE",
      role = source_role,
      source = source,
      copied = copied
    ),
  fig_manifest %>%
    transmute(
      category = "FIGURE_SOURCE",
      role = figure_role,
      source = original_file,
      copied = copied
    ),
  word_audit %>%
    transmute(
      category = "FROZEN_WORDING",
      role = wording_role,
      source = source,
      copied = copied
    )
)

write_csv(
  manifest,
  file.path(
    OUT,
    "07_GLOBAL_MANUSCRIPT_SOURCE_MANIFEST.csv"
  )
)

# ------------------------------------------------------------
# 9. Final QC
# ------------------------------------------------------------

critical_main_ok <- all(
  main_audit$copied[
    main_audit$source_role %in%
      c(
        "MAIN_T2_CROSS_COHORT_LONGITUDINAL",
        "MAIN_T3_CORE_TRAJECTORY_ENDPOINTS",
        "MAIN_T4_EXTERNAL_VALIDATION",
        "MAIN_T5_INFECTION_SOURCE_ECOLOGY"
      )
  ]
)

critical_wording_ok <- all(word_audit$copied)

qc <- tibble(
  step95A2_complete = TRUE,
  five_global_branches_source_anchored =
    all(evidence_registry$source_anchored),
  critical_main_table_sources_copied =
    critical_main_ok,
  frozen_wording_copied =
    critical_wording_ok,
  figure_candidate_files_found =
    sum(fig_manifest$copied, na.rm = TRUE),
  ready_for_figure_table_production =
    (
      all(evidence_registry$source_anchored) &&
      critical_main_ok &&
      critical_wording_ok
    )
)

write_csv(
  qc,
  file.path(
    OUT,
    "08_STEP95B_READINESS.csv"
  )
)

writeLines(
  c(
    "STEP95B MANUSCRIPT SOURCE PACK ASSEMBLY",
    "",
    paste0(
      "Five evidence lines source-anchored: ",
      qc$five_global_branches_source_anchored
    ),
    paste0(
      "Critical main-table sources copied: ",
      qc$critical_main_table_sources_copied
    ),
    paste0(
      "Frozen wording copied: ",
      qc$frozen_wording_copied
    ),
    paste0(
      "Figure candidate files copied: ",
      qc$figure_candidate_files_found
    ),
    paste0(
      "Ready for figure/table production: ",
      qc$ready_for_figure_table_production
    ),
    "",
    "No new inferential analysis was performed.",
    "Use this Step95B pack as the manuscript-production source layer."
  ),
  file.path(
    OUT,
    "09_STEP95B_INTERPRETATION.txt"
  )
)

if (!isTRUE(qc$ready_for_figure_table_production)) {
  stop("Step95B source-pack readiness failed.")
}

writeLines(
  c(
    paste0("Completed: ", Sys.time()),
    "Ready for figure/table production: TRUE",
    "STEP95B COMPLETE"
  ),
  file.path(
    OUT,
    "_STEP95B_COMPLETE.ok"
  )
)

cat("STEP95B COMPLETE\n")
