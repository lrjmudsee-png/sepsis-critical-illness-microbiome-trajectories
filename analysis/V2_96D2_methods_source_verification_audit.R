
# ============================================================
# V2_96D2 METHODS SOURCE VERIFICATION AUDIT
#
# Purpose:
# Resolve the four [VERIFY] items in Step96D using the actual
# local analysis source code and local software environment.
#
# IMPORTANT:
# - Reads source code only.
# - Does NOT rerun microbiome processing.
# - Does NOT rerun statistical models.
# - Does NOT overwrite frozen outputs.
# ============================================================

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(stringr)
  library(tidyr)
  library(purrr)
  library(tibble)
})

ROOT <- "E:/sepsis_project"
CODE_ROOT <- file.path(ROOT,"code")
MAIN_CODE <- file.path(CODE_ROOT,"03_data_processing")
RESULTS <- file.path(ROOT,"results")

OUT <- file.path(
  RESULTS,
  "V2_96D2_METHODS_SOURCE_VERIFICATION_AUDIT"
)

dir.create(OUT,recursive=TRUE,showWarnings=FALSE)

# ------------------------------------------------------------
# Utilities
# ------------------------------------------------------------
read_code <- function(f){
  tryCatch(
    readLines(f,warn=FALSE,encoding="UTF-8"),
    error=function(e) {
      tryCatch(
        readLines(f,warn=FALSE),
        error=function(e2) character()
      )
    }
  )
}

context_rows <- function(lines, idx, radius=5){
  if(length(idx)==0) return(tibble())

  bind_rows(lapply(idx,function(i){
    lo <- max(1,i-radius)
    hi <- min(length(lines),i+radius)

    tibble(
      hit_line=i,
      context_start=lo,
      context_end=hi,
      excerpt=paste(
        sprintf("%05d | %s",lo:hi,lines[lo:hi]),
        collapse="\n"
      )
    )
  }))
}

scan_file <- function(f,pattern,radius=5){
  lines <- read_code(f)
  if(length(lines)==0) return(tibble())

  idx <- grep(
    pattern,
    lines,
    ignore.case=TRUE,
    perl=TRUE
  )

  if(length(idx)==0) return(tibble())

  context_rows(lines,idx,radius) %>%
    mutate(
      file=f,
      filename=basename(f),
      matched_text=lines[hit_line],
      .before=1
    )
}

all_code <- list.files(
  CODE_ROOT,
  recursive=TRUE,
  full.names=TRUE,
  pattern="\\.(R|r|Rmd|rmd|ps1|PS1)$"
)

# ============================================================
# VERIFY-1
# Exact upstream QC / denoising / taxonomy evidence
# ============================================================

preprocess_pattern <- paste(
  c(
    "dada2",
    "filterAndTrim",
    "learnErrors",
    "derepFastq",
    "dada\\(",
    "mergePairs",
    "removeBimeraDenovo",
    "assignTaxonomy",
    "addSpecies",
    "SILVA",
    "DECIPHER",
    "IdTaxa",
    "vsearch",
    "uchime",
    "cluster[_-]?size",
    "cluster[_-]?fast",
    "OTU97",
    "OTU",
    "ASV",
    "raref",
    "chimera",
    "taxonomy"
  ),
  collapse="|"
)

preprocess_hits <- bind_rows(
  lapply(
    all_code,
    function(f) scan_file(
      f,
      preprocess_pattern,
      radius=4
    )
  )
)

# Score files by number of preprocessing hits.
preprocess_file_summary <- preprocess_hits %>%
  count(file,filename,name="n_processing_hits") %>%
  arrange(desc(n_processing_hits),filename)

write_csv(
  preprocess_file_summary,
  file.path(
    OUT,
    "01A_VERIFY1_PROCESSING_SCRIPT_CANDIDATES.csv"
  )
)

# Keep all evidence, but include likely project association if present.
projects <- c(
  "PRJNA691455",
  "PRJEB82425",
  "PRJNA516701",
  "PRJNA851469",
  "PRJNA578267",
  "PRJNA430161",
  "PRJNA1166732",
  "PRJNA978257",
  "PRJNA1010969",
  "CRA002354"
)

detect_projects <- function(txt){
  hit <- projects[
    vapply(
      projects,
      function(p) str_detect(txt,fixed(p)),
      logical(1)
    )
  ]

  if(length(hit)==0) "UNSPECIFIED_OR_GLOBAL_PIPELINE"
  else paste(hit,collapse=";")
}

if(nrow(preprocess_hits)>0){
  preprocess_hits <- preprocess_hits %>%
    mutate(
      project_hint=vapply(
        paste(file,excerpt),
        detect_projects,
        character(1)
      )
    )
}

write_csv(
  preprocess_hits,
  file.path(
    OUT,
    "01B_VERIFY1_UPSTREAM_PROCESSING_EVIDENCE.csv"
  )
)

# Dedicated CRA evidence: final corrected processing scripts.
cra_files <- all_code[
  str_detect(
    basename(all_code),
    regex("94B1|CRA002354",ignore_case=TRUE)
  )
]

cra_evidence <- bind_rows(
  lapply(
    cra_files,
    function(f) scan_file(
      f,
      "vsearch|chimera|cluster|OTU|SILVA|taxonomy|FASTA|quality|nonchimer",
      radius=6
    )
  )
)

write_csv(
  cra_evidence,
  file.path(
    OUT,
    "01C_VERIFY1_CRA002354_PROCESSING_EVIDENCE.csv"
  )
)

# ============================================================
# VERIFY-2
# Exact Step88A2 mixed-model specification
# ============================================================

step88 <- file.path(
  MAIN_CODE,
  "88A2_longitudinal_diversity_and_displacement_FIXED.R"
)

if(!file.exists(step88)){
  stop("Missing Step88A2 source script: ",step88)
}

step88_pattern <- paste(
  c(
    "lmer\\(",
    "glmer\\(",
    "lme\\(",
    "lm\\(",
    "anova\\(",
    "formula",
    "random",
    "patient",
    "subject",
    "time",
    "day",
    "likelihood",
    "LRT",
    "REML"
  ),
  collapse="|"
)

step88_evidence <- scan_file(
  step88,
  step88_pattern,
  radius=7
)

write_csv(
  step88_evidence,
  file.path(
    OUT,
    "02_VERIFY2_STEP88A2_MODEL_EVIDENCE.csv"
  )
)

# Also produce raw numbered source for exact manual audit.
step88_lines <- read_code(step88)

writeLines(
  sprintf(
    "%05d | %s",
    seq_along(step88_lines),
    step88_lines
  ),
  file.path(
    OUT,
    "02B_STEP88A2_NUMBERED_SOURCE.txt"
  )
)

# ============================================================
# VERIFY-3
# Exact CRA002354 Step94B2 model family / random effect / LRT
# ============================================================

step94_files <- all_code[
  str_detect(
    basename(all_code),
    regex(
      "94B2A2|94B2A3|94B2B|94B2C|94B2D",
      ignore_case=TRUE
    )
  )
]

if(length(step94_files)==0){
  stop("Could not find final Step94B2 source scripts.")
}

step94_pattern <- paste(
  c(
    "lmer\\(",
    "glmer\\(",
    "lme\\(",
    "lm\\(",
    "anova\\(",
    "formula",
    "random",
    "patient",
    "subject",
    "source",
    "time",
    "day",
    "interaction",
    "likelihood",
    "LRT",
    "REML",
    "ML"
  ),
  collapse="|"
)

step94_evidence <- bind_rows(
  lapply(
    step94_files,
    function(f) scan_file(
      f,
      step94_pattern,
      radius=7
    )
  )
)

write_csv(
  step94_evidence,
  file.path(
    OUT,
    "03_VERIFY3_STEP94B2_MODEL_EVIDENCE.csv"
  )
)

# Numbered source for final formal-analysis scripts.
for(f in step94_files){

  lines <- read_code(f)

  safe <- str_replace_all(
    basename(f),
    "[^A-Za-z0-9._-]",
    "_"
  )

  writeLines(
    sprintf(
      "%05d | %s",
      seq_along(lines),
      lines
    ),
    file.path(
      OUT,
      paste0(
        "03B_NUMBERED_SOURCE_",
        safe,
        ".txt"
      )
    )
  )
}

# ============================================================
# VERIFY-4
# Exact R/package/VSEARCH versions
# ============================================================

# Parse packages referenced by library(), require(), requireNamespace()
pkg_pattern <- "(?:library|require|requireNamespace)\\s*\\(\\s*[\"']?([A-Za-z0-9._]+)"

package_usage <- list()

for(f in all_code){

  lines <- read_code(f)
  if(length(lines)==0) next

  mm <- str_match(
    lines,
    regex(pkg_pattern,ignore_case=TRUE)
  )

  keep <- !is.na(mm[,2])

  if(any(keep)){

    package_usage[[f]] <- tibble(
      file=f,
      filename=basename(f),
      line_number=which(keep),
      package=mm[keep,2],
      source_line=lines[keep]
    )
  }
}

package_usage_df <- bind_rows(package_usage)

write_csv(
  package_usage_df,
  file.path(
    OUT,
    "04A_PACKAGE_USAGE_IN_SOURCE_CODE.csv"
  )
)

used_pkgs <- sort(
  unique(package_usage_df$package)
)

installed <- installed.packages()

pkg_versions <- tibble(
  software=used_pkgs,
  type="R_PACKAGE",
  version=vapply(
    used_pkgs,
    function(p){
      if(p %in% rownames(installed)){
        as.character(packageVersion(p))
      }else{
        NA_character_
      }
    },
    character(1)
  ),
  evidence="packageVersion() in current local R environment"
)

r_info <- tibble(
  software=c(
    "R",
    "R_PLATFORM",
    "R_ARCH"
  ),
  type="R_RUNTIME",
  version=c(
    R.version.string,
    R.version$platform,
    R.version$arch
  ),
  evidence="R.version"
)

# Try to locate VSEARCH.
vsearch_paths <- unique(
  c(
    Sys.which("vsearch"),
    Sys.which("vsearch.exe")
  )
)

vsearch_paths <- vsearch_paths[
  nzchar(vsearch_paths)
]

# Also inspect source code for explicit .exe paths containing vsearch.
vsearch_path_hits <- preprocess_hits %>%
  filter(
    str_detect(
      tolower(excerpt),
      "vsearch"
    )
  )

if(nrow(vsearch_path_hits)>0){

  text <- paste(
    vsearch_path_hits$excerpt,
    collapse="\n"
  )

  # Extract quoted strings ending in vsearch(.exe)
  found <- unlist(
    str_extract_all(
      text,
      "[A-Za-z]:[/\\\\][^\"']*vsearch(?:\\.exe)?"
    )
  )

  found <- str_replace_all(found,"\\\\","/")
  vsearch_paths <- unique(c(vsearch_paths,found))
}

vsearch_rows <- list()

if(length(vsearch_paths)==0){

  vsearch_rows[[1]] <- tibble(
    software="VSEARCH",
    type="EXTERNAL_TOOL",
    version=NA_character_,
    evidence="VSEARCH not found on PATH; inspect source-path evidence files."
  )

}else{

  for(vp in vsearch_paths){

    outv <- tryCatch(
      system2(
        vp,
        "--version",
        stdout=TRUE,
        stderr=TRUE
      ),
      error=function(e) paste(
        "VERSION_QUERY_FAILED:",
        conditionMessage(e)
      )
    )

    vsearch_rows[[vp]] <- tibble(
      software="VSEARCH",
      type="EXTERNAL_TOOL",
      version=paste(outv,collapse=" | "),
      evidence=paste0("system2('",vp,"','--version')")
    )
  }
}

version_df <- bind_rows(
  r_info,
  pkg_versions,
  bind_rows(vsearch_rows)
)

write_csv(
  version_df,
  file.path(
    OUT,
    "04_VERIFY4_SOFTWARE_VERSIONS.csv"
  )
)

# ============================================================
# 5. Automated source-supported resolution summary
# ============================================================

verify_summary <- c(
  "V2 STEP96D2 METHODS SOURCE VERIFICATION AUDIT",
  "",
  "This step did NOT rerun microbiome processing or statistical models.",
  "",
  paste0(
    "VERIFY-1 upstream-processing evidence hits: ",
    nrow(preprocess_hits),
    " across ",
    nrow(preprocess_file_summary),
    " source files."
  ),
  paste0(
    "VERIFY-2 Step88A2 model evidence hits: ",
    nrow(step88_evidence),
    "."
  ),
  paste0(
    "VERIFY-3 Step94B2 model evidence hits: ",
    nrow(step94_evidence),
    " across ",
    length(step94_files),
    " final/source-freeze scripts."
  ),
  paste0(
    "VERIFY-4 software/version rows: ",
    nrow(version_df),
    "."
  ),
  "",
  "NEXT ACTION:",
  "Upload this audit ZIP. The extracted source excerpts and numbered source files can then be used to replace the four [VERIFY] markers in Methods v1 with source-grounded wording.",
  "",
  "STEP96D2 COMPLETE"
)

writeLines(
  verify_summary,
  file.path(
    OUT,
    "05_VERIFY_RESOLUTION_SUMMARY.txt"
  )
)

# Copy Methods v1 if available, so next review is self-contained.
methods_v1 <- file.path(
  RESULTS,
  "V2_96D_METHODS_MANUSCRIPTIZATION",
  "Methods_Draft_v1_FROM_FROZEN_V2_ANALYSES.txt"
)

if(file.exists(methods_v1)){
  file.copy(
    methods_v1,
    file.path(
      OUT,
      "Methods_Draft_v1_REFERENCE_COPY.txt"
    ),
    overwrite=TRUE
  )
}

writeLines(
  c(
    paste0("Completed: ",Sys.time()),
    "STEP96D2 COMPLETE",
    "Source-verification audit completed.",
    "No analysis rerun."
  ),
  file.path(
    OUT,
    "_STEP96D2_COMPLETE.txt"
  )
)

cat("STEP96D2 COMPLETE\n")
