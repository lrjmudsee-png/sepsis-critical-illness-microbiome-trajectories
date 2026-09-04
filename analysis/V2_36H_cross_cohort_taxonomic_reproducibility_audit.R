
# ============================================================
# V2_36H CROSS-COHORT TAXONOMIC REPRODUCIBILITY AUDIT
#
# Goal:
# Identify analysis-ready genus/family abundance tables for the
# frozen longitudinal cohorts, without merging ASVs/OTUs across cohorts.
#
# Future analysis (Step36I if feasible):
#   1. estimate within-cohort early->late taxonomic effect vectors;
#   2. harmonize only taxon labels/effect directions;
#   3. quantify cross-cohort effect-vector concordance;
#   4. compare taxonomic reproducibility with ecological displacement.
#
# NO new inferential microbiome model is run in this audit.
# ============================================================

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(stringr)
  library(purrr)
  library(tidyr)
})

ROOT <- "E:/sepsis_project"
DATA <- file.path(ROOT,"data")
RESULTS <- file.path(ROOT,"results")
OUT <- file.path(
  RESULTS,
  "V2_36H_CROSS_COHORT_TAXONOMIC_REPRODUCIBILITY_AUDIT"
)
dir.create(OUT,recursive=TRUE,showWarnings=FALSE)

target_projects <- c(
  "PRJNA691455",
  "PRJNA851469",
  "PRJNA516701",
  "PRJEB82425",
  "PRJNA578267",
  "PRJNA430161",
  "PRJNA1166732"
)

primary_projects <- c(
  "PRJNA691455",
  "PRJNA851469",
  "PRJNA516701"
)

# ------------------------------------------
# 1. Find candidate taxonomic tables
# ------------------------------------------
roots <- c(
  file.path(DATA,"_V2_ANALYSIS_READY"),
  RESULTS
)
roots <- roots[dir.exists(roots)]

all_csv <- unlist(lapply(
  roots,
  function(r) list.files(
    r,
    pattern="\\.csv$",
    recursive=TRUE,
    full.names=TRUE,
    ignore.case=TRUE
  )
))

all_csv <- unique(all_csv)

detect_project <- function(path) {
  hit <- target_projects[str_detect(path,fixed(target_projects))]
  if(length(hit)==0) return(NA_character_)
  hit[1]
}

path_low <- tolower(all_csv)

candidate_idx <-
  map_lgl(all_csv, ~ any(str_detect(.x,fixed(target_projects)))) &
  str_detect(
    path_low,
    "genus|family|taxon|taxonomy|abundance|relative|feature"
  ) &
  !str_detect(
    path_low,
    "bray|distance|ordination|permanova|alpha|shannon|simpson|observed|figure|plot|summary|audit"
  )

candidates <- all_csv[candidate_idx]

# ------------------------------------------
# 2. Header/sample inspection
# ------------------------------------------
read_small <- function(f) {
  tryCatch(
    read_csv(
      f,
      n_max=50,
      show_col_types=FALSE,
      progress=FALSE
    ),
    error=function(e) NULL
  )
}

inspect_one <- function(f) {

  x <- read_small(f)

  if(is.null(x)) {
    return(tibble(
      project=detect_project(f),
      file=f,
      size_mb=round(file.info(f)$size/1024^2,3),
      readable=FALSE,
      n_columns=NA_integer_,
      n_preview_rows=NA_integer_,
      has_sample_id=FALSE,
      has_patient_id=FALSE,
      has_time=FALSE,
      has_taxon_label=FALSE,
      likely_long_format=FALSE,
      likely_wide_abundance=FALSE,
      taxonomy_level_hint=NA_character_,
      columns=""
    ))
  }

  cn <- names(x)
  low <- tolower(cn)

  sample_hits <- cn[str_detect(
    low,
    "sample|run|accession|biosample|specimen"
  )]

  patient_hits <- cn[str_detect(
    low,
    "patient|subject|participant|individual"
  )]

  time_hits <- cn[str_detect(
    low,
    "day|time|visit|week|sampling"
  )]

  taxon_hits <- cn[str_detect(
    low,
    "genus|family|taxon|taxonomy|feature"
  )]

  abundance_hits <- cn[str_detect(
    low,
    "abundance|count|reads|relative|ra"
  )]

  level <- case_when(
    str_detect(tolower(f),"genus") ||
      any(str_detect(low,"genus")) ~ "GENUS",
    str_detect(tolower(f),"family") ||
      any(str_detect(low,"family")) ~ "FAMILY",
    TRUE ~ "UNKNOWN"
  )

  likely_long <- (
    length(taxon_hits)>0 &&
    length(abundance_hits)>0
  )

  # heuristic: many numeric columns = possible wide feature matrix
  num_cols <- sum(vapply(x,is.numeric,logical(1)))
  likely_wide <- ncol(x)>=20 && num_cols>=10

  tibble(
    project=detect_project(f),
    file=f,
    size_mb=round(file.info(f)$size/1024^2,3),
    readable=TRUE,
    n_columns=ncol(x),
    n_preview_rows=nrow(x),
    has_sample_id=length(sample_hits)>0,
    has_patient_id=length(patient_hits)>0,
    has_time=length(time_hits)>0,
    has_taxon_label=length(taxon_hits)>0,
    likely_long_format=likely_long,
    likely_wide_abundance=likely_wide,
    taxonomy_level_hint=level,
    columns=paste(cn,collapse=" | ")
  )
}

inventory <- if(length(candidates)>0) {
  map_dfr(candidates,inspect_one)
} else {
  tibble(
    project=character(),
    file=character(),
    size_mb=double(),
    readable=logical(),
    n_columns=integer(),
    n_preview_rows=integer(),
    has_sample_id=logical(),
    has_patient_id=logical(),
    has_time=logical(),
    has_taxon_label=logical(),
    likely_long_format=logical(),
    likely_wide_abundance=logical(),
    taxonomy_level_hint=character(),
    columns=character()
  )
}

write_csv(
  inventory,
  file.path(OUT,"01_taxonomic_table_candidate_inventory.csv")
)

# ------------------------------------------
# 3. Candidate scoring
# ------------------------------------------
scored <- inventory %>%
  mutate(
    score=
      3*as.integer(readable) +
      3*as.integer(has_sample_id) +
      2*as.integer(has_patient_id) +
      2*as.integer(has_time) +
      2*as.integer(has_taxon_label) +
      2*as.integer(likely_long_format) +
      1*as.integer(likely_wide_abundance) +
      2*as.integer(taxonomy_level_hint %in% c("GENUS","FAMILY"))
  ) %>%
  arrange(project,desc(score),size_mb)

write_csv(
  scored,
  file.path(OUT,"02_scored_taxonomic_candidates.csv")
)

# Best candidates per project/level
best <- scored %>%
  filter(readable) %>%
  group_by(project,taxonomy_level_hint) %>%
  slice_max(order_by=score,n=3,with_ties=FALSE) %>%
  ungroup()

write_csv(
  best,
  file.path(OUT,"03_best_candidates_by_project.csv")
)

# ------------------------------------------
# 4. Check exact longitudinal metadata availability
# ------------------------------------------
meta_candidates <- all_csv[
  map_lgl(all_csv, ~ any(str_detect(.x,fixed(target_projects)))) &
  str_detect(
    tolower(all_csv),
    "metadata|analysis_metadata|step87b|primary"
  )
]

inspect_meta <- function(f) {
  x <- read_small(f)
  if(is.null(x)) return(tibble())

  cn <- names(x)
  low <- tolower(cn)

  tibble(
    project=detect_project(f),
    file=f,
    n_columns=ncol(x),
    has_sample_id=any(str_detect(low,"sample|run|accession|biosample|specimen")),
    has_patient_id=any(str_detect(low,"patient|subject|participant|individual")),
    has_time=any(str_detect(low,"day|time|visit|week|sampling")),
    columns=paste(cn,collapse=" | ")
  )
}

meta_inv <- if(length(meta_candidates)>0) {
  map_dfr(meta_candidates,inspect_meta)
} else {
  tibble(
    project=character(),
    file=character(),
    n_columns=integer(),
    has_sample_id=logical(),
    has_patient_id=logical(),
    has_time=logical(),
    columns=character()
  )
}

write_csv(
  meta_inv,
  file.path(OUT,"04_longitudinal_metadata_candidate_inventory.csv")
)

# ------------------------------------------
# 5. Cohort-level feasibility
# ------------------------------------------
tax_status <- scored %>%
  group_by(project) %>%
  summarise(
    n_taxonomic_candidates=n(),
    max_score=max(score,na.rm=TRUE),
    any_genus=any(taxonomy_level_hint=="GENUS"),
    any_family=any(taxonomy_level_hint=="FAMILY"),
    any_linkable=
      any(readable & (has_sample_id|has_patient_id)),
    .groups="drop"
  )

meta_status <- meta_inv %>%
  group_by(project) %>%
  summarise(
    n_metadata_candidates=n(),
    any_metadata_linkable=
      any(has_sample_id & has_patient_id & has_time),
    .groups="drop"
  )

feas <- tibble(project=target_projects) %>%
  left_join(tax_status,by="project") %>%
  left_join(meta_status,by="project") %>%
  mutate(
    across(
      c(
        n_taxonomic_candidates,
        max_score,
        n_metadata_candidates
      ),
      ~replace_na(.x,0)
    ),
    across(
      c(
        any_genus,
        any_family,
        any_linkable,
        any_metadata_linkable
      ),
      ~replace_na(.x,FALSE)
    ),
    primary_role=project %in% primary_projects,
    preliminary_ready=
      any_linkable &
      any_metadata_linkable &
      (any_genus|any_family),
    priority=case_when(
      preliminary_ready & primary_role ~ "HIGH",
      preliminary_ready ~ "MODERATE",
      any_linkable ~ "LOW",
      TRUE ~ "NOT_READY"
    )
  )

write_csv(
  feas,
  file.path(OUT,"05_cohort_taxonomic_reproducibility_feasibility.csv")
)

# ------------------------------------------
# 6. Decision summary
# ------------------------------------------
n_primary_ready <- sum(
  feas$preliminary_ready & feas$primary_role
)

n_total_ready <- sum(feas$preliminary_ready)

lines <- c(
  "V2 STEP36H CROSS-COHORT TAXONOMIC REPRODUCIBILITY AUDIT",
  "",
  paste0("Target longitudinal cohorts: ",length(target_projects)),
  paste0("Primary progressive cohorts ready: ",n_primary_ready," / ",length(primary_projects)),
  paste0("All cohorts preliminary ready: ",n_total_ready," / ",length(target_projects)),
  ""
)

if(n_primary_ready>=2) {
  lines <- c(
    lines,
    "DECISION: PROCEED TO STEP36I.",
    "",
    "At least two primary progressive cohorts have candidate genus/family abundance tables",
    "and linkable longitudinal metadata.",
    "",
    "Step36I should estimate taxonomic changes WITHIN each cohort only,",
    "then harmonize taxon labels/effect directions across cohorts.",
    "Do NOT merge ASV/OTU count matrices across cohorts."
  )
} else {
  lines <- c(
    lines,
    "DECISION: DO NOT YET RUN CROSS-COHORT TAXONOMIC REPRODUCIBILITY MODELS.",
    "",
    "Inspect 03_best_candidates_by_project.csv and recover missing genus/family tables first."
  )
}

writeLines(
  lines,
  file.path(OUT,"06_DECISION_SUMMARY.txt")
)

writeLines(
  c(
    paste0("Completed: ",Sys.time()),
    "STEP36H COMPLETE",
    "No inferential microbiome model was run."
  ),
  file.path(OUT,"_STEP36H_COMPLETE.txt")
)

cat("STEP36H COMPLETE\n")
cat("Primary cohorts ready:",n_primary_ready,"\n")
