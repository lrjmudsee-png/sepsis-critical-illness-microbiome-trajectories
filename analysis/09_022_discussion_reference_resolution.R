
# ============================================================
# V2_96G DISCUSSION REFERENCE RESOLUTION
#
# Purpose:
# Resolve all Discussion-level [REF...] placeholders in the
# Step96F manuscript core using a verified literature pool.
#
# IMPORTANT:
# - Stable citation keys are inserted, NOT final reference numbers.
# - Final numbering is deferred until the Introduction is written.
# - NO statistical analysis is rerun.
# ============================================================

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(stringr)
  library(tibble)
})

ROOT <- "E:/sepsis_project"
RESULTS <- file.path(ROOT,"results")

INFILE <- file.path(
  RESULTS,
  "V2_96F_MANUSCRIPT_CORE_ASSEMBLY_AND_GAP_AUDIT",
  "01_MANUSCRIPT_CORE_DRAFT_v1.txt"
)

OUT <- file.path(
  RESULTS,
  "V2_96G_DISCUSSION_REFERENCE_RESOLUTION"
)

dir.create(OUT,recursive=TRUE,showWarnings=FALSE)

if(!file.exists(INFILE)){
  stop("Step96F manuscript core not found: ",INFILE)
}

txt <- paste(
  readLines(
    INFILE,
    warn=FALSE,
    encoding="UTF-8"
  ),
  collapse="\n"
)

# ------------------------------------------------------------
# 1. Exact placeholder -> verified citation-key mapping
# ------------------------------------------------------------

mapping <- tribble(
  ~placeholder, ~replacement,

  "[REF: longitudinal microbiome disruption during critical illness and sepsis]",
  "[CIT:SCHLECHTE2023;KITSIOS2024;THEOCHARIDOU2026]",

  "[REF: cross-cohort inconsistency of sepsis-associated taxa; ecological-state approaches in microbiome research]",
  "[CIT:GEISTLINGER2024;CLAUSEN2022;NEARING2022;FASSARELLA2021]",

  "[REF_SCHLECHTE_2023]",
  "[CIT:SCHLECHTE2023]",

  "[REF_KITSIOS_2024]",
  "[CIT:KITSIOS2024]",

  "[REF_THEOCHARIDOU_2026]",
  "[CIT:THEOCHARIDOU2026]",

  "[REF_GEISTLINGER_2024]",
  "[CIT:GEISTLINGER2024]",

  "[REF_CLAUSEN_2022]",
  "[CIT:CLAUSEN2022]",

  "[REF: trauma microbiome studies; compositional microbiome methods; sepsis biomarker studies]",
  "[CIT:HOWARD2017;MUNLEY2024;AGUDELO2020;FERNANDES2014]",

  "[REF: gut microbiome across sepsis infection sources; gut-lung axis / systemic critical-illness ecology]",
  "[CIT:HAAK2017;ADELMAN2020;KITSIOS2024]",

  "[REF: robustness and reproducibility of 16S genus-level differential abundance]",
  "[CIT:CLAUSEN2022;NEARING2022;YANG2022;YANGCHEN2023]",

  "[REF: longitudinal precision microbiome / ecological trajectory methods]",
  "[CIT:FASSARELLA2021;KITSIOS2024;THEOCHARIDOU2026]",

  "[REF: concluding contextual citations]",
  "[CIT:THEOCHARIDOU2026;GEISTLINGER2024;FASSARELLA2021]"
)

for(i in seq_len(nrow(mapping))){
  old <- mapping$placeholder[[i]]
  new <- mapping$replacement[[i]]

  if(!str_detect(txt,fixed(old))){
    stop("Expected reference placeholder not found: ",old)
  }

  txt <- str_replace_all(
    txt,
    fixed(old),
    new
  )
}

write_csv(
  mapping,
  file.path(
    OUT,
    "01_PLACEHOLDER_TO_CITATION_KEY_MAP.csv"
  )
)

# ------------------------------------------------------------
# 2. Keep final reference list intentionally deferred
# ------------------------------------------------------------

old_ref_list <- "[REFERENCE LIST TO BE ASSEMBLED AFTER ALL REF PLACEHOLDERS ARE RESOLVED]"
new_ref_list <- "[FINAL NUMBERED REFERENCE LIST WILL BE ASSEMBLED AFTER INTRODUCTION CITATIONS ARE ADDED]"

if(str_detect(txt,fixed(old_ref_list))){
  txt <- str_replace_all(
    txt,
    fixed(old_ref_list),
    new_ref_list
  )
}

# ------------------------------------------------------------
# 3. Guard against unresolved Discussion REF placeholders
# ------------------------------------------------------------

unresolved_ref <- unique(
  unlist(
    str_extract_all(
      txt,
      "\\[REF[^\\]]*\\]"
    )
  )
)

unresolved_ref <- unresolved_ref[nzchar(unresolved_ref)]

if(length(unresolved_ref)>0){
  write_csv(
    tibble(unresolved_reference_placeholder=unresolved_ref),
    file.path(
      OUT,
      "02_UNRESOLVED_REFERENCE_PLACEHOLDERS.csv"
    )
  )

  stop(
    "Unresolved [REF...] placeholder(s) remain: ",
    paste(unresolved_ref,collapse=" ; ")
  )
}

write_csv(
  tibble(unresolved_reference_placeholder=character()),
  file.path(
    OUT,
    "02_UNRESOLVED_REFERENCE_PLACEHOLDERS.csv"
  )
)

# ------------------------------------------------------------
# 4. Audit inserted CIT keys
# ------------------------------------------------------------

cit_blocks <- unique(
  unlist(
    str_extract_all(
      txt,
      "\\[CIT:[^\\]]+\\]"
    )
  )
)

keys <- unique(
  unlist(
    str_split(
      str_replace_all(
        str_replace_all(cit_blocks,"\\[CIT:",""),
        "\\]",""
      ),
      ";"
    )
  )
)

keys <- sort(keys[nzchar(keys)])

write_csv(
  tibble(citation_key=keys),
  file.path(
    OUT,
    "03_CITATION_KEYS_USED_IN_MANUSCRIPT.csv"
  )
)

# ------------------------------------------------------------
# 5. Write manuscript core v2
# ------------------------------------------------------------

OUTFILE <- file.path(
  OUT,
  "Manuscript_Core_v2_DISCUSSION_REFS_RESOLVED.txt"
)

writeLines(
  str_split(txt,"\n")[[1]],
  OUTFILE,
  useBytes=TRUE
)

# ------------------------------------------------------------
# 6. QC
# ------------------------------------------------------------

qc <- tibble(
  check=c(
    "generic_REF_placeholders_remaining",
    "stable_CIT_keys_present",
    "final_numbering_deferred",
    "analysis_rerun"
  ),
  passed=c(
    length(unresolved_ref)==0,
    length(keys)>0,
    str_detect(
      txt,
      fixed("FINAL NUMBERED REFERENCE LIST WILL BE ASSEMBLED AFTER INTRODUCTION CITATIONS ARE ADDED")
    ),
    TRUE
  )
)

# For analysis_rerun, TRUE here means "confirmed NO rerun".
qc$check[qc$check=="analysis_rerun"] <- "confirmed_no_analysis_rerun"

write_csv(
  qc,
  file.path(
    OUT,
    "04_STEP96G_QC.csv"
  )
)

if(!all(qc$passed)){
  stop(
    "Step96G QC failed: ",
    paste(qc$check[!qc$passed],collapse=", ")
  )
}

writeLines(
  c(
    "V2 STEP96G DISCUSSION REFERENCE RESOLUTION",
    "",
    paste0("Resolved Discussion REF placeholders: ",nrow(mapping)),
    paste0("Unique stable citation keys now used: ",length(keys)),
    "Final numeric reference numbering intentionally deferred until Introduction is written.",
    "",
    "No statistical analysis was rerun.",
    "",
    paste0("Output: ",OUTFILE),
    "",
    "STEP96G COMPLETE"
  ),
  file.path(
    OUT,
    "05_STEP96G_SUMMARY.txt"
  )
)

writeLines(
  c(
    paste0("Completed: ",Sys.time()),
    "STEP96G COMPLETE",
    paste0("Stable CIT keys: ",length(keys)),
    "Unresolved Discussion REF placeholders: 0"
  ),
  file.path(
    OUT,
    "_STEP96G_COMPLETE.txt"
  )
)

cat("STEP96G COMPLETE\n")
