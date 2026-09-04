# ============================================================
# Sepsis V2 - Step95A
# GLOBAL EVIDENCE SYNTHESIS + MANUSCRIPT / FIGURE MAP
#
# Purpose:
# Consolidate the already-frozen V2 evidence lines.
# NO new hypothesis testing is performed.
#
# Frozen branches:
# 1) Cross-cohort longitudinal evidence (Step93W)
# 2) Core sepsis taxonomic-trajectory heterogeneity (Step93U)
# 3) External Control-Trauma-Sepsis validation (Step93X4)
# 4) CRA002354 infection-source ecology (Step94B2D)
# 5) CRA002354 exploratory genus taxonomy (Step94B3C)
#
# This step creates a single authoritative registry for manuscript writing.
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

OUT <- file.path(
  RESULTS,
  "V2_35A_GLOBAL_EVIDENCE_SYNTHESIS_AND_MANUSCRIPT_MAP"
)

dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------------------
# Helpers
# ------------------------------------------------------------

first_existing <- function(paths) {
  x <- paths[file.exists(paths)]
  if (length(x) == 0) return(NA_character_)
  x[1]
}

find_first <- function(dir, pattern) {
  if (!dir.exists(dir)) return(NA_character_)
  x <- list.files(
    dir,
    pattern = pattern,
    recursive = TRUE,
    full.names = TRUE,
    ignore.case = TRUE
  )
  if (length(x) == 0) return(NA_character_)
  x[1]
}

read_text_safe <- function(path) {
  if (is.na(path) || !file.exists(path)) return(NA_character_)
  paste(readLines(path, warn = FALSE, encoding = "UTF-8"), collapse = "\n")
}

# ------------------------------------------------------------
# 1. Expected frozen branches
# ------------------------------------------------------------

D93W <- file.path(
  RESULTS,
  "V2_33W_LONGITUDINAL_EVIDENCE_FREEZE"
)

D93X4 <- file.path(
  RESULTS,
  "V2_33X4_EXTERNAL_VALIDATION_EVIDENCE_FREEZE"
)

D94B2D <- file.path(
  RESULTS,
  "V2_34B2D_INFECTION_SOURCE_EVIDENCE_FREEZE"
)

D94B3C <- file.path(
  RESULTS,
  "V2_34B3C_CRA002354_GENUS_EVIDENCE_FREEZE"
)

# Step93U directory name may vary slightly; locate by prefix.
d93u_candidates <- list.dirs(
  RESULTS,
  recursive = FALSE,
  full.names = TRUE
)

D93U <- d93u_candidates[
  grepl(
    "93U",
    basename(d93u_candidates),
    ignore.case = TRUE
  )
]

if (length(D93U) == 0) {
  D93U <- NA_character_
} else {
  D93U <- D93U[1]
}

branch_audit <- tibble(
  branch = c(
    "CROSS_COHORT_LONGITUDINAL",
    "CORE_SEPSIS_TAXONOMIC_TRAJECTORY",
    "EXTERNAL_CONTROL_TRAUMA_SEPSIS",
    "CRA002354_SOURCE_ECOLOGY",
    "CRA002354_SOURCE_GENUS"
  ),
  expected_dir = c(
    D93W,
    D93U,
    D93X4,
    D94B2D,
    D94B3C
  ),
  directory_exists = c(
    dir.exists(D93W),
    ifelse(is.na(D93U), FALSE, dir.exists(D93U)),
    dir.exists(D93X4),
    dir.exists(D94B2D),
    dir.exists(D94B3C)
  )
)

write_csv(
  branch_audit,
  file.path(OUT, "01_FROZEN_BRANCH_AUDIT.csv")
)

# ------------------------------------------------------------
# 2. Pull authoritative frozen evidence where available
# ------------------------------------------------------------

# CRA ecology
cra_eco_file <- file.path(
  D94B2D,
  "05_FINAL_INFECTION_SOURCE_BRANCH_FREEZE.csv"
)

cra_eco <- if (file.exists(cra_eco_file)) {
  read_csv(cra_eco_file, show_col_types = FALSE)
} else {
  tibble()
}

# CRA genus
cra_genus_file <- file.path(
  D94B3C,
  "03_FINAL_GENUS_BRANCH_EVIDENCE_FREEZE.csv"
)

cra_genus <- if (file.exists(cra_genus_file)) {
  read_csv(cra_genus_file, show_col_types = FALSE)
} else {
  tibble()
}

# External validation: find freeze CSV/TXT robustly.
x4_csv <- find_first(
  D93X4,
  "freeze.*\\.csv$|evidence.*\\.csv$|summary.*\\.csv$"
)

x4_txt <- find_first(
  D93X4,
  "interpretation.*\\.txt$|wording.*\\.txt$"
)

# Longitudinal freeze: find key summary/registry.
w_csv <- find_first(
  D93W,
  "freeze.*\\.csv$|evidence.*\\.csv$|summary.*\\.csv$"
)

w_txt <- find_first(
  D93W,
  "interpretation.*\\.txt$|wording.*\\.txt$"
)

# Core trajectory heterogeneity.
u_csv <- if (!is.na(D93U)) {
  find_first(
    D93U,
    "summary.*\\.csv$|evidence.*\\.csv$|trajectory.*\\.csv$"
  )
} else {
  NA_character_
}

u_txt <- if (!is.na(D93U)) {
  find_first(
    D93U,
    "interpretation.*\\.txt$|summary.*\\.txt$"
  )
} else {
  NA_character_
}

source_files <- tibble(
  branch = c(
    "CROSS_COHORT_LONGITUDINAL",
    "CORE_SEPSIS_TAXONOMIC_TRAJECTORY",
    "EXTERNAL_CONTROL_TRAUMA_SEPSIS",
    "CRA002354_SOURCE_ECOLOGY",
    "CRA002354_SOURCE_GENUS"
  ),
  authoritative_csv = c(
    w_csv,
    u_csv,
    x4_csv,
    ifelse(file.exists(cra_eco_file), cra_eco_file, NA_character_),
    ifelse(file.exists(cra_genus_file), cra_genus_file, NA_character_)
  ),
  authoritative_text = c(
    w_txt,
    u_txt,
    x4_txt,
    file.path(
      D94B2D,
      "06_FROZEN_MANUSCRIPT_WORDING_AND_GUARDRAILS.txt"
    ),
    file.path(
      D94B3C,
      "05_FROZEN_GENUS_MANUSCRIPT_WORDING_AND_GUARDRAILS.txt"
    )
  )
) %>%
  mutate(
    csv_exists = !is.na(authoritative_csv) &
      file.exists(authoritative_csv),
    text_exists = !is.na(authoritative_text) &
      file.exists(authoritative_text)
  )

write_csv(
  source_files,
  file.path(OUT, "02_AUTHORITATIVE_SOURCE_FILE_REGISTRY.csv")
)

# ------------------------------------------------------------
# 3. Global evidence registry
#    Uses already-frozen conclusions only.
# ------------------------------------------------------------

cra_eco_tier <- if (nrow(cra_eco) > 0 &&
                    "evidence_tier" %in% names(cra_eco)) {
  cra_eco$evidence_tier[1]
} else {
  "UNKNOWN"
}

cra_genus_tier <- if (nrow(cra_genus) > 0 &&
                      "evidence_tier" %in% names(cra_genus)) {
  cra_genus$evidence_tier[1]
} else {
  "UNKNOWN"
}

registry <- tibble(
  evidence_line = c(
    "Cross-cohort longitudinal ecological displacement",
    "Core-sepsis taxonomic trajectory heterogeneity",
    "External severe-trauma comparison",
    "Infection-source modification of longitudinal ecology",
    "Infection-source genus-level taxonomy"
  ),
  cohort_or_scope = c(
    "Multiple longitudinal cohorts",
    "PRJNA691455 repeated sepsis",
    "PRJNA1010969 Control-Trauma-Sepsis",
    "CRA002354 pulmonary vs recorded non-pulmonary",
    "CRA002354 pulmonary vs recorded non-pulmonary"
  ),
  role = c(
    "PRIMARY_MULTI_COHORT",
    "MECHANISTIC_INTERPRETIVE",
    "SUPPORTIVE_EXTERNAL",
    "PRIMARY_SOURCE_BRANCH",
    "EXPLORATORY_SUPPORTIVE"
  ),
  frozen_conclusion = c(
    "Progressive ecological displacement is reproducible across several longitudinal critical-illness cohorts, while non-sepsis longitudinal control shows a different/reconvergent pattern in sensitivity analysis.",
    "Patients can reach strong ecological displacement through heterogeneous genus-level taxonomic routes; stable taxonomic trajectory subtypes were not supported.",
    "Sepsis shows ecological state differentiation beyond severe trauma in CLR/Aitchison sensitivity analyses, but this does not establish a reproducible sepsis-specific taxonomic signature.",
    paste0(
      "No clear infection-source modification of longitudinal ecological displacement; frozen tier: ",
      cra_eco_tier
    ),
    paste0(
      "Exploratory baseline genus signals with limited coverage robustness and no primary longitudinal genus signal; frozen tier: ",
      cra_genus_tier
    )
  ),
  manuscript_weight = c(
    "HIGH",
    "HIGH",
    "MODERATE",
    "MODERATE",
    "LOW"
  )
)

write_csv(
  registry,
  file.path(OUT, "03_GLOBAL_V2_EVIDENCE_REGISTRY.csv")
)

# ------------------------------------------------------------
# 4. Manuscript claim hierarchy
# ------------------------------------------------------------

claims <- tibble(
  claim_level = c(
    "MAIN_CLAIM_1",
    "MAIN_CLAIM_2",
    "SUPPORTING_CLAIM_1",
    "SUPPORTING_CLAIM_2",
    "SUPPLEMENTARY_CLAIM"
  ),
  recommended_claim = c(
    "Longitudinal critical illness is accompanied by progressive within-patient microbiome ecological displacement across multiple cohorts.",
    "This ecological displacement can arise through heterogeneous taxonomic trajectories rather than a single conserved genus-level path.",
    "External comparison supports sepsis-associated ecological state differentiation beyond severe trauma, while not establishing a universal sepsis-specific taxonomic signature.",
    "Infection source did not clearly modify the overall longitudinal ecological displacement trajectory in CRA002354.",
    "Several baseline genus differences were observed by infection source, but robustness to taxonomic-coverage filtering was limited and no primary longitudinal genus interaction survived FDR correction."
  ),
  allowed_strength = c(
    "PRIMARY",
    "PRIMARY_INTERPRETIVE",
    "SUPPORTIVE",
    "SUPPORTIVE",
    "EXPLORATORY_ONLY"
  )
)

write_csv(
  claims,
  file.path(OUT, "04_MANUSCRIPT_CLAIM_HIERARCHY.csv")
)

# ------------------------------------------------------------
# 5. Figure architecture
# ------------------------------------------------------------

figmap <- tibble(
  figure = c(
    "Figure 1",
    "Figure 2",
    "Figure 3",
    "Figure 4",
    "Figure 5",
    "Supplementary"
  ),
  purpose = c(
    "Study/data architecture and cohort roles",
    "Cross-cohort longitudinal ecological displacement",
    "Core-sepsis heterogeneous taxonomic trajectories",
    "External Control-Trauma-Sepsis ecological validation",
    "CRA002354 infection-source longitudinal ecology",
    "Exploratory genus-level source findings, EII, robustness analyses"
  ),
  status = c(
    "TO_ASSEMBLE",
    "EVIDENCE_FROZEN",
    "EVIDENCE_FROZEN",
    "EVIDENCE_FROZEN",
    "EVIDENCE_FROZEN",
    "EVIDENCE_FROZEN"
  ),
  source_branch = c(
    "GLOBAL",
    "STEP93W",
    "STEP93U",
    "STEP93X4",
    "STEP94B2D",
    "STEP94B3C + EII"
  )
)

write_csv(
  figmap,
  file.path(OUT, "05_MANUSCRIPT_FIGURE_MAP.csv")
)

# ------------------------------------------------------------
# 6. Guardrails
# ------------------------------------------------------------

guardrails <- c(
  "Do not claim sepsis-specific longitudinal instability based only on cohorts that include ICU-background progression.",
  "Do not claim stable taxonomic trajectory subtypes from the PRJNA691455 k=2 exploratory clustering because the split was 9 vs 1.",
  "Do not present EII as a clinical prognostic score or independently validated biomarker.",
  "Do not interpret the external cohort as validating a universal sepsis-specific taxonomic signature.",
  "Do not claim pulmonary and non-pulmonary microbiome trajectories are identical; state that no clear source modification was detected.",
  "Do not treat the unrestricted CRA002354 baseline PERMANOVA as clean centroid separation because dispersion was also significant and the signal disappeared in the <=Day3 sensitivity.",
  "Do not present Enterococcus as a robust infection-source biomarker.",
  "Do not present sensitivity-only Bacteroides longitudinal significance as a validated source-specific trajectory.",
  "Keep CRA002354 genus findings exploratory/supportive."
)

writeLines(
  c(
    "V2 GLOBAL MANUSCRIPT GUARDRAILS",
    "",
    paste0("- ", guardrails)
  ),
  file.path(OUT, "06_GLOBAL_MANUSCRIPT_GUARDRAILS.txt")
)

# ------------------------------------------------------------
# 7. Global narrative
# ------------------------------------------------------------

narrative <- c(
  "STEP95A GLOBAL V2 EVIDENCE SYNTHESIS",
  "",
  "Recommended central narrative:",
  "Across longitudinal critical-illness cohorts, the gut microbiome shows progressive ecological displacement over time. The taxonomic routes underlying this displacement are heterogeneous rather than converging on a single reproducible genus-level trajectory. External comparison with severe trauma supports sepsis-associated ecological state differentiation, but does not establish a universal sepsis-specific taxonomic signature. Within CRA002354, pulmonary versus recorded non-pulmonary infection source did not clearly modify the longitudinal ecological displacement trajectory, while genus-level source differences were limited to exploratory baseline signals with incomplete robustness.",
  "",
  "Interpretive emphasis:",
  "The manuscript should prioritize ecological trajectory and displacement over a universal taxonomic biomarker narrative.",
  "",
  "Next manuscript phase:",
  "1. Assemble final main figures from frozen evidence.",
  "2. Build final main/supplementary result tables.",
  "3. Draft Results in the same order as the figure architecture.",
  "4. Draft Discussion around shared ecological displacement + heterogeneous taxonomic routes + limits of source specificity."
)

writeLines(
  narrative,
  file.path(OUT, "07_GLOBAL_V2_NARRATIVE_AND_NEXT_PHASE.txt")
)

# ------------------------------------------------------------
# 8. Freeze audit
# ------------------------------------------------------------

complete_dirs <- c(
  D93W,
  D93X4,
  D94B2D,
  D94B3C
)

freeze_audit <- tibble(
  directory = complete_dirs,
  exists = dir.exists(complete_dirs),
  complete_flag_count = sapply(
    complete_dirs,
    function(d) {
      if (!dir.exists(d)) return(0L)
      length(
        list.files(
          d,
          pattern = "_STEP.*COMPLETE\\.ok$",
          full.names = TRUE,
          ignore.case = TRUE
        )
      )
    }
  )
)

write_csv(
  freeze_audit,
  file.path(OUT, "08_GLOBAL_FREEZE_AUDIT.csv")
)

writeLines(
  c(
    paste0("Completed: ", Sys.time()),
    "No new inferential analysis performed.",
    "STEP95A COMPLETE"
  ),
  file.path(OUT, "_STEP95A_COMPLETE.ok")
)

cat("STEP95A COMPLETE\n")
