# ============================================================
# Sepsis V2 - Step95A2
# GLOBAL EVIDENCE SYNTHESIS SOURCE-ANCHOR FIX
#
# Purpose:
# Step95A correctly consolidated four frozen branches, but failed to
# locate the authoritative Step93U core-sepsis taxonomic-trajectory
# output directory.
#
# Step95A2:
# 1) scans the results tree for Step93U / trajectory-heterogeneity evidence;
# 2) scores candidate directories using filenames + text content;
# 3) requires a credible source anchor for the core-sepsis trajectory branch;
# 4) rebuilds the global source registry and synthesis only after anchoring;
# 5) refuses to mark the global synthesis complete if Step93U cannot be found.
#
# NO new hypothesis testing is performed.
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
  "V2_35A2_GLOBAL_EVIDENCE_SYNTHESIS_SOURCE_ANCHOR_FIX"
)

dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------------------
# Helpers
# ------------------------------------------------------------

safe_read_text <- function(path, max_lines = 300) {
  if (!file.exists(path)) return("")
  ext <- tolower(tools::file_ext(path))
  if (!ext %in% c("txt","csv","tsv","md","ok")) return("")

  x <- tryCatch(
    readLines(
      path,
      warn = FALSE,
      encoding = "UTF-8",
      n = max_lines
    ),
    error = function(e) character()
  )

  paste(x, collapse = "\n")
}

find_first <- function(dir, patterns) {
  if (!dir.exists(dir)) return(NA_character_)

  allf <- list.files(
    dir,
    recursive = TRUE,
    full.names = TRUE
  )

  for (pat in patterns) {
    z <- allf[
      str_detect(
        tolower(basename(allf)),
        regex(pat, ignore_case = TRUE)
      )
    ]
    if (length(z) > 0) return(z[1])
  }

  NA_character_
}

# ------------------------------------------------------------
# 1. Scan results tree for Step93U candidates
# ------------------------------------------------------------

all_files <- list.files(
  RESULTS,
  recursive = TRUE,
  full.names = TRUE,
  include.dirs = FALSE
)

# Freeze the locator to the upstream result universe that existed when
# Step95A2 was defined.  Later manuscript/reproducibility packages copy
# Step93U wording and provenance into new folders; scanning those folders
# makes the candidate audit self-referential and changes its dimensions on
# every rerun even though the selected scientific source is unchanged.
relative_files <- substring(
  all_files,
  nchar(RESULTS) + 2L
)

top_level_result_dir <- sub(
  "[/\\\\].*$",
  "",
  relative_files
)

downstream_or_self <- str_detect(
  top_level_result_dir,
  regex(
    "^V2_(35A2|35B|35C|35D|35E|35F|35G|36|96|97)",
    ignore_case = TRUE
  )
)

all_files <- all_files[!downstream_or_self]

candidate_files <- all_files[
  str_detect(
    tolower(all_files),
    regex(
      "93u|taxonomic.*trajectory|trajectory.*heterogeneity|heterogeneous.*taxonomic|pairwise.*heterogeneity|cosine|latest.*eii|signed.*delta",
      ignore_case = TRUE
    )
  )
]

# Also inspect small text-like files even if filenames don't match.
text_like <- all_files[
  tolower(tools::file_ext(all_files)) %in%
    c("txt","csv","tsv","md","ok")
]

# Avoid reading extremely large files.
text_like <- text_like[
  file.info(text_like)$size <= 5 * 1024^2
]

content_hits <- character()

patterns <- c(
  "STEP93U",
  "pairwise heterogeneity",
  "mean cosine",
  "latest EII",
  "trajectory distance",
  "signed genera",
  "heterogeneous taxonomic",
  "0.733",
  "0.0158",
  "0.0234"
)

for (f in text_like) {

  txt <- safe_read_text(f)

  if (
    any(
      vapply(
        patterns,
        function(p) {
          str_detect(
            txt,
            fixed(p, ignore_case = TRUE)
          )
        },
        logical(1)
      )
    )
  ) {
    content_hits <- c(
      content_hits,
      f
    )
  }
}

candidate_files <- unique(
  c(
    candidate_files,
    content_hits
  )
)

# Build candidate-directory scores.
candidate_dirs <- unique(
  dirname(candidate_files)
)

score_dir <- function(d) {

  fs <- candidate_files[
    dirname(candidate_files) == d |
    startsWith(
      candidate_files,
      paste0(d, .Platform$file.sep)
    )
  ]

  names_blob <- paste(
    basename(fs),
    collapse = " "
  )

  txt_blob <- paste(
    vapply(
      fs[
        tolower(tools::file_ext(fs)) %in%
          c("txt","csv","tsv","md","ok") &
        file.info(fs)$size <=
          2 * 1024^2
      ],
      safe_read_text,
      character(1)
    ),
    collapse = "\n"
  )

  score <- 0L
  reasons <- character()

  add <- function(cond, pts, why) {
    if (isTRUE(cond)) {
      score <<- score + pts
      reasons <<- c(reasons, why)
    }
  }

  add(
    str_detect(
      tolower(d),
      regex("93u", ignore_case = TRUE)
    ),
    5,
    "directory_name_contains_93U"
  )

  add(
    str_detect(
      tolower(names_blob),
      regex("93u", ignore_case = TRUE)
    ),
    4,
    "filename_contains_93U"
  )

  add(
    str_detect(
      txt_blob,
      regex("STEP93U", ignore_case = TRUE)
    ),
    5,
    "content_contains_STEP93U"
  )

  add(
    str_detect(
      paste(names_blob, txt_blob),
      regex(
        "pairwise.*heterogeneity|mean cosine|trajectory distance",
        ignore_case = TRUE
      )
    ),
    3,
    "trajectory_heterogeneity_terms"
  )

  add(
    str_detect(
      paste(names_blob, txt_blob),
      regex(
        "latest EII|0\\.733|0\\.0158|0\\.0234",
        ignore_case = TRUE
      )
    ),
    2,
    "known_STEP93U_numeric_or_EII_anchor"
  )

  add(
    str_detect(
      paste(names_blob, txt_blob),
      regex(
        "268.*genera|signed.*genera",
        ignore_case = TRUE
      )
    ),
    2,
    "signed_268_genera_anchor"
  )

  tibble(
    candidate_dir = d,
    score = score,
    reasons = paste(
      unique(reasons),
      collapse = ";"
    ),
    matched_files = length(fs)
  )
}

candidate_scores <- if (
  length(candidate_dirs) > 0
) {
  bind_rows(
    lapply(
      candidate_dirs,
      score_dir
    )
  ) %>%
    arrange(
      desc(score),
      desc(matched_files)
    )
} else {
  tibble(
    candidate_dir = character(),
    score = integer(),
    reasons = character(),
    matched_files = integer()
  )
}

write_csv(
  candidate_scores,
  file.path(
    OUT,
    "01_STEP93U_CANDIDATE_DIRECTORY_SCORES.csv"
  )
)

write_csv(
  tibble(
    candidate_file =
      candidate_files
  ),
  file.path(
    OUT,
    "02_STEP93U_CANDIDATE_FILES.csv"
  )
)

if (
  nrow(candidate_scores) == 0 ||
  candidate_scores$score[1] < 5
) {

  writeLines(
    c(
      "STEP95A2 STOPPED",
      "",
      "No sufficiently credible Step93U source directory was located automatically.",
      "The global synthesis is therefore NOT marked complete.",
      "",
      "Inspect:",
      "01_STEP93U_CANDIDATE_DIRECTORY_SCORES.csv",
      "02_STEP93U_CANDIDATE_FILES.csv"
    ),
    file.path(
      OUT,
      "03_STEP93U_LOCATOR_DECISION.txt"
    )
  )

  stop(
    "Step93U authoritative source could not be anchored. Global synthesis not finalized."
  )
}

D93U <- candidate_scores$candidate_dir[1]

writeLines(
  c(
    "STEP93U SOURCE LOCATED",
    paste0(
      "Selected directory: ",
      D93U
    ),
    paste0(
      "Score: ",
      candidate_scores$score[1]
    ),
    paste0(
      "Reasons: ",
      candidate_scores$reasons[1]
    )
  ),
  file.path(
    OUT,
    "03_STEP93U_LOCATOR_DECISION.txt"
  )
)

# ------------------------------------------------------------
# 2. Other frozen branch directories
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

branch_dirs <- tibble(
  branch = c(
    "CROSS_COHORT_LONGITUDINAL",
    "CORE_SEPSIS_TAXONOMIC_TRAJECTORY",
    "EXTERNAL_CONTROL_TRAUMA_SEPSIS",
    "CRA002354_SOURCE_ECOLOGY",
    "CRA002354_SOURCE_GENUS"
  ),
  directory = c(
    D93W,
    D93U,
    D93X4,
    D94B2D,
    D94B3C
  ),
  exists = c(
    dir.exists(D93W),
    dir.exists(D93U),
    dir.exists(D93X4),
    dir.exists(D94B2D),
    dir.exists(D94B3C)
  )
)

write_csv(
  branch_dirs,
  file.path(
    OUT,
    "04_FINAL_FROZEN_BRANCH_AUDIT.csv"
  )
)

if (!all(branch_dirs$exists)) {
  stop(
    "At least one required frozen evidence branch is missing."
  )
}

# ------------------------------------------------------------
# 3. Authoritative source files
# ------------------------------------------------------------

w_csv <- find_first(
  D93W,
  c(
    "frozen.*evidence.*matrix.*\\.csv$",
    "evidence.*\\.csv$",
    "summary.*\\.csv$"
  )
)

w_txt <- find_first(
  D93W,
  c(
    "methods.*wording.*\\.txt$",
    "interpretation.*\\.txt$"
  )
)

u_csv <- find_first(
  D93U,
  c(
    "summary.*\\.csv$",
    "evidence.*\\.csv$",
    "trajectory.*\\.csv$",
    "heterogeneity.*\\.csv$"
  )
)

u_txt <- find_first(
  D93U,
  c(
    "interpretation.*\\.txt$",
    "summary.*\\.txt$",
    "wording.*\\.txt$"
  )
)

x_csv <- find_first(
  D93X4,
  c(
    "frozen.*summary.*\\.csv$",
    "external.*summary.*\\.csv$",
    "evidence.*\\.csv$"
  )
)

x_txt <- find_first(
  D93X4,
  c(
    "methods.*wording.*\\.txt$",
    "interpretation.*\\.txt$"
  )
)

cra_eco_csv <- file.path(
  D94B2D,
  "05_FINAL_INFECTION_SOURCE_BRANCH_FREEZE.csv"
)

cra_eco_txt <- file.path(
  D94B2D,
  "06_FROZEN_MANUSCRIPT_WORDING_AND_GUARDRAILS.txt"
)

cra_genus_csv <- file.path(
  D94B3C,
  "03_FINAL_GENUS_BRANCH_EVIDENCE_FREEZE.csv"
)

cra_genus_txt <- file.path(
  D94B3C,
  "05_FROZEN_GENUS_MANUSCRIPT_WORDING_AND_GUARDRAILS.txt"
)

source_registry <- tibble(
  branch = branch_dirs$branch,
  authoritative_csv = c(
    w_csv,
    u_csv,
    x_csv,
    cra_eco_csv,
    cra_genus_csv
  ),
  authoritative_text = c(
    w_txt,
    u_txt,
    x_txt,
    cra_eco_txt,
    cra_genus_txt
  )
) %>%
  mutate(
    csv_exists =
      !is.na(authoritative_csv) &
      file.exists(authoritative_csv),
    text_exists =
      !is.na(authoritative_text) &
      file.exists(authoritative_text)
  )

write_csv(
  source_registry,
  file.path(
    OUT,
    "05_FINAL_AUTHORITATIVE_SOURCE_FILE_REGISTRY.csv"
  )
)

# Core Step93U must have at least one textual/CSV anchor.
u_ok <- source_registry %>%
  filter(
    branch ==
      "CORE_SEPSIS_TAXONOMIC_TRAJECTORY"
  ) %>%
  summarise(
    ok =
      any(
        csv_exists |
        text_exists
      )
  ) %>%
  pull(ok)

if (!isTRUE(u_ok)) {
  stop(
    "Step93U directory was located but no authoritative CSV/TXT source file was identified."
  )
}

# ------------------------------------------------------------
# 4. Pull CRA frozen tiers
# ------------------------------------------------------------

cra_eco <- read_csv(
  cra_eco_csv,
  show_col_types = FALSE
)

cra_genus <- read_csv(
  cra_genus_csv,
  show_col_types = FALSE
)

cra_eco_tier <- if (
  "evidence_tier" %in%
    names(cra_eco)
) {
  cra_eco$evidence_tier[1]
} else {
  "UNKNOWN"
}

cra_genus_tier <- if (
  "evidence_tier" %in%
    names(cra_genus)
) {
  cra_genus$evidence_tier[1]
} else {
  "UNKNOWN"
}

# ------------------------------------------------------------
# 5. Final global evidence registry
# ------------------------------------------------------------

registry <- tibble(
  evidence_line = c(
    "Cross-cohort longitudinal ecological displacement",
    "Core-sepsis taxonomic trajectory heterogeneity",
    "External severe-trauma comparison",
    "Infection-source modification of longitudinal ecology",
    "Infection-source genus-level taxonomy"
  ),
  source_branch = c(
    "STEP93W",
    "STEP93U",
    "STEP93X4",
    "STEP94B2D",
    "STEP94B3C"
  ),
  source_anchored = c(
    !is.na(w_csv) |
      !is.na(w_txt),
    TRUE,
    !is.na(x_csv) |
      !is.na(x_txt),
    TRUE,
    TRUE
  ),
  role = c(
    "PRIMARY_MULTI_COHORT",
    "PRIMARY_INTERPRETIVE",
    "SUPPORTIVE_EXTERNAL",
    "SUPPORTIVE_SOURCE_BRANCH",
    "EXPLORATORY_SUPPORTIVE"
  ),
  frozen_conclusion = c(
    "Progressive ecological displacement is reproduced across several longitudinal critical-illness cohorts, while the non-sepsis longitudinal control shows a different/reconvergent pattern in sensitivity analysis.",
    "Core-sepsis patients show heterogeneous taxonomic routes accompanying strong ecological displacement; exploratory clustering did not support stable trajectory subtypes.",
    "External comparison supports sepsis-associated ecological state differentiation beyond severe trauma, but does not establish a universal sepsis-specific taxonomic signature.",
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
  file.path(
    OUT,
    "06_FINAL_GLOBAL_V2_EVIDENCE_REGISTRY.csv"
  )
)

if (!all(registry$source_anchored)) {
  stop(
    "Not every global evidence line has a source anchor."
  )
}

# ------------------------------------------------------------
# 6. Claim hierarchy
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
    "External comparison supports sepsis-associated ecological state differentiation beyond severe trauma without establishing a universal sepsis-specific taxonomic signature.",
    "Infection source did not clearly modify the overall longitudinal ecological displacement trajectory in CRA002354.",
    "Several baseline genus differences were observed by infection source, but robustness to sample-level taxonomic-coverage filtering was limited and no primary longitudinal genus interaction survived FDR correction."
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
  file.path(
    OUT,
    "07_FINAL_MANUSCRIPT_CLAIM_HIERARCHY.csv"
  )
)

# ------------------------------------------------------------
# 7. Figure / result architecture
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
    "Study architecture, cohort roles, and analysis flow",
    "Cross-cohort longitudinal ecological displacement",
    "Core-sepsis heterogeneous taxonomic trajectories",
    "External Control-Trauma-Sepsis ecological validation",
    "CRA002354 infection-source longitudinal ecology",
    "Exploratory genus-level source findings, EII, and robustness analyses"
  ),
  source_branch = c(
    "GLOBAL",
    "STEP93W",
    "STEP93U",
    "STEP93X4",
    "STEP94B2D",
    "STEP94B3C + EII"
  ),
  evidence_status = c(
    "TO_ASSEMBLE",
    "FROZEN",
    "FROZEN",
    "FROZEN",
    "FROZEN",
    "FROZEN"
  )
)

write_csv(
  figmap,
  file.path(
    OUT,
    "08_FINAL_MANUSCRIPT_FIGURE_MAP.csv"
  )
)

result_sections <- tibble(
  section_order = 1:5,
  results_section = c(
    "Cohort architecture and longitudinal analysis framework",
    "Progressive ecological displacement across longitudinal critical-illness cohorts",
    "Heterogeneous taxonomic routes accompany ecological displacement in core sepsis",
    "External ecological differentiation beyond severe trauma",
    "Infection source does not clearly modify longitudinal ecological displacement"
  ),
  main_source = c(
    "GLOBAL + metadata",
    "STEP93W",
    "STEP93U",
    "STEP93X4",
    "STEP94B2D"
  ),
  supplementary_link = c(
    "Cohort tables",
    "Full cohort-level model tables",
    "Trajectory metrics + EII",
    "CLR/Aitchison sensitivity",
    "STEP94B3C genus results"
  )
)

write_csv(
  result_sections,
  file.path(
    OUT,
    "09_RESULTS_SECTION_ORDER.csv"
  )
)

# ------------------------------------------------------------
# 8. Guardrails + final narrative
# ------------------------------------------------------------

guardrails <- c(
  "Do not claim sepsis-specific longitudinal instability solely from cohorts that also include ICU-background progression.",
  "Do not claim stable taxonomic trajectory subtypes from PRJNA691455 exploratory k=2 clustering because the split was 9 vs 1.",
  "Do not present EII as a clinical prognostic score or independently validated biomarker.",
  "Do not interpret the external cohort as validating a universal sepsis-specific taxonomic signature.",
  "Do not claim pulmonary and non-pulmonary trajectories are identical; state that no clear source modification was detected.",
  "Do not interpret the unrestricted CRA002354 baseline PERMANOVA as clean centroid separation because dispersion was significant and the <=Day3 sensitivity was negative.",
  "Do not present Enterococcus as a robust infection-source biomarker.",
  "Do not present coverage-sensitive Bacteroides significance as a validated longitudinal source signal.",
  "Keep CRA002354 genus results exploratory/supportive."
)

writeLines(
  c(
    "STEP95A2 FINAL GLOBAL MANUSCRIPT GUARDRAILS",
    "",
    paste0("- ", guardrails)
  ),
  file.path(
    OUT,
    "10_FINAL_GLOBAL_MANUSCRIPT_GUARDRAILS.txt"
  )
)

writeLines(
  c(
    "STEP95A2 FINAL GLOBAL V2 NARRATIVE",
    "",
    "Central narrative:",
    "Across longitudinal critical-illness cohorts, the gut microbiome undergoes progressive within-patient ecological displacement. In core sepsis, this shared ecological movement is accompanied by heterogeneous taxonomic routes rather than a single conserved genus-level trajectory. External comparison with severe trauma supports sepsis-associated ecological state differentiation without establishing a universal sepsis-specific taxonomic signature. Within CRA002354, pulmonary versus recorded non-pulmonary infection source did not clearly modify the longitudinal ecological displacement trajectory; genus-level source differences were restricted to exploratory baseline signals with limited robustness.",
    "",
    "Manuscript emphasis:",
    "Prioritize ecological trajectory/displacement over a universal taxonomic biomarker narrative.",
    "",
    "Next phase:",
    "1. Build final Figure 1-5 from frozen branches.",
    "2. Create main and supplementary result tables.",
    "3. Draft Results in the order defined in 09_RESULTS_SECTION_ORDER.csv.",
    "4. Draft Discussion around shared ecological displacement, heterogeneous taxonomic routes, external ecological differentiation, and limits of infection-source specificity.",
    "",
    paste0(
      "Core Step93U source anchor: ",
      D93U
    )
  ),
  file.path(
    OUT,
    "11_FINAL_GLOBAL_V2_NARRATIVE_AND_NEXT_PHASE.txt"
  )
)

# ------------------------------------------------------------
# 9. Final freeze audit
# ------------------------------------------------------------

final_audit <- tibble(
  branch = source_registry$branch,
  directory_exists =
    branch_dirs$exists[
      match(
        source_registry$branch,
        branch_dirs$branch
      )
    ],
  authoritative_csv_exists =
    source_registry$csv_exists,
  authoritative_text_exists =
    source_registry$text_exists,
  source_anchor_pass =
    source_registry$csv_exists |
    source_registry$text_exists
)

write_csv(
  final_audit,
  file.path(
    OUT,
    "12_FINAL_GLOBAL_SOURCE_ANCHOR_AUDIT.csv"
  )
)

if (!all(final_audit$source_anchor_pass)) {
  stop(
    "Global source-anchor audit failed."
  )
}

writeLines(
  c(
    paste0(
      "Completed: ",
      Sys.time()
    ),
    paste0(
      "Step93U anchor: ",
      D93U
    ),
    "All five evidence lines source-anchored: TRUE",
    "No new inferential analysis performed.",
    "STEP95A2 COMPLETE"
  ),
  file.path(
    OUT,
    "_STEP95A2_COMPLETE.ok"
  )
)

cat("STEP95A2 COMPLETE\n")
