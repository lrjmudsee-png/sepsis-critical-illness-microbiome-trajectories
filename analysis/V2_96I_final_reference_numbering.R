
# ============================================================
# V2_96I FINAL REFERENCE NUMBERING
#
# Purpose:
# Convert stable [CIT:KEY1;KEY2] citation keys in the complete
# Introduction/Methods/Results/Discussion manuscript core into
# numbered citations based strictly on FIRST APPEARANCE.
#
# Also assembles the final numbered reference list.
#
# NO statistical analysis is rerun.
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
  "V2_96H_INTRODUCTION_INTEGRATION",
  "Manuscript_Core_v3_WITH_INTRODUCTION.txt"
)

POOLFILE <- file.path(
  RESULTS,
  "V2_96H_INTRODUCTION_INTEGRATION",
  "01_VERIFIED_REFERENCE_POOL_EXPANDED.csv"
)

OUT <- file.path(
  RESULTS,
  "V2_96I_FINAL_REFERENCE_NUMBERING"
)

dir.create(OUT,recursive=TRUE,showWarnings=FALSE)

for(f in c(INFILE,POOLFILE)){
  if(!file.exists(f)){
    stop("Missing Step96I input: ",f)
  }
}

txt <- paste(
  readLines(
    INFILE,
    warn=FALSE,
    encoding="UTF-8"
  ),
  collapse="\n"
)

pool <- read_csv(
  POOLFILE,
  show_col_types=FALSE
)

# ------------------------------------------------------------
# 1. Extract citation keys in FIRST-APPEARANCE order
# ------------------------------------------------------------

blocks <- str_extract_all(
  txt,
  "\\[CIT:[^\\]]+\\]"
)[[1]]

if(length(blocks)==0){
  stop("No [CIT:...] citation blocks found.")
}

key_order <- character()

for(b in blocks){

  ks <- b %>%
    str_remove("^\\[CIT:") %>%
    str_remove("\\]$") %>%
    str_split(";") %>%
    .[[1]] %>%
    str_trim()

  for(k in ks){
    if(!(k %in% key_order)){
      key_order <- c(key_order,k)
    }
  }
}

missing <- setdiff(
  key_order,
  pool$citation_key
)

if(length(missing)>0){
  stop(
    "Citation key(s) missing from verified pool: ",
    paste(missing,collapse=", ")
  )
}

# Only references actually cited in the manuscript enter final numbering.
number_map <- tibble(
  reference_number=seq_along(key_order),
  citation_key=key_order
) %>%
  left_join(
    pool,
    by="citation_key"
  )

write_csv(
  number_map,
  file.path(
    OUT,
    "01_FINAL_REFERENCE_NUMBER_MAP.csv"
  )
)

# ------------------------------------------------------------
# 2. Numeric citation formatting
#    BMC/Microbiome style: [1], [1, 2], [3-5]
# ------------------------------------------------------------

compress_numbers <- function(nums){

  nums <- sort(unique(as.integer(nums)))

  if(length(nums)==1){
    return(as.character(nums))
  }

  runs <- split(
    nums,
    cumsum(
      c(
        TRUE,
        diff(nums)!=1
      )
    )
  )

  parts <- vapply(
    runs,
    function(x){

      if(length(x)>=3){
        paste0(min(x),"-",max(x))
      } else if(length(x)==2){
        paste(x,collapse=", ")
      } else {
        as.character(x)
      }
    },
    character(1)
  )

  paste(parts,collapse=", ")
}

replace_block <- function(block){

  ks <- block %>%
    str_remove("^\\[CIT:") %>%
    str_remove("\\]$") %>%
    str_split(";") %>%
    .[[1]] %>%
    str_trim()

  nums <- number_map$reference_number[
    match(
      ks,
      number_map$citation_key
    )
  ]

  if(any(is.na(nums))){
    stop("Failed to number citation block: ",block)
  }

  paste0(
    "[",
    compress_numbers(nums),
    "]"
  )
}

unique_blocks <- unique(blocks)

block_registry <- tibble(
  citation_block=unique_blocks,
  numeric_citation=vapply(
    unique_blocks,
    replace_block,
    character(1)
  )
)

write_csv(
  block_registry,
  file.path(
    OUT,
    "02_CITATION_BLOCK_TO_NUMBER_MAP.csv"
  )
)

for(i in seq_len(nrow(block_registry))){

  txt <- str_replace_all(
    txt,
    fixed(
      block_registry$citation_block[[i]]
    ),
    block_registry$numeric_citation[[i]]
  )
}

# ------------------------------------------------------------
# 3. Build Microbiome/BMC-style numbered references
# ------------------------------------------------------------

clean_authors <- function(x){
  x <- str_replace_all(x,"; ","; ")
  str_trim(x)
}

clean_journal <- function(x){
  str_trim(x)
}

format_reference <- function(row){

  authors <- clean_authors(row$authors)
  title <- str_trim(row$title)
  journal <- clean_journal(row$journal)
  year <- as.character(row$year)
  volume <- ifelse(
    is.na(row$volume) || row$volume=="",
    "",
    as.character(row$volume)
  )
  pages <- ifelse(
    is.na(row$pages_or_article) || row$pages_or_article=="",
    "",
    as.character(row$pages_or_article)
  )
  doi <- ifelse(
    is.na(row$doi) || row$doi=="",
    "",
    as.character(row$doi)
  )

  bibliographic_tail <- if(volume!="" && pages!=""){
    paste0(
      year,
      ";",
      volume,
      ":",
      pages,
      "."
    )
  } else if(volume!=""){
    paste0(
      year,
      ";",
      volume,
      "."
    )
  } else {
    paste0(
      year,
      "."
    )
  }

  ref <- paste0(
    row$reference_number,
    ". ",
    authors,
    ". ",
    title,
    ". ",
    journal,
    ". ",
    bibliographic_tail
  )

  if(doi!=""){
    ref <- paste0(
      ref,
      " https://doi.org/",
      doi,
      "."
    )
  }

  ref
}

ref_lines <- vapply(
  seq_len(nrow(number_map)),
  function(i){
    format_reference(
      as.list(
        number_map[i,]
      )
    )
  },
  character(1)
)

writeLines(
  ref_lines,
  file.path(
    OUT,
    "03_FINAL_NUMBERED_REFERENCE_LIST.txt"
  ),
  useBytes=TRUE
)

# ------------------------------------------------------------
# 4. Replace References placeholder
# ------------------------------------------------------------

placeholder <- "[FINAL NUMBERED REFERENCE LIST WILL BE ASSEMBLED AFTER INTRODUCTION CITATIONS ARE ADDED]"

if(!str_detect(txt,fixed(placeholder))){
  stop("Final reference-list placeholder not found.")
}

txt <- str_replace(
  txt,
  fixed(placeholder),
  paste(
    ref_lines,
    collapse="\n\n"
  )
)

# ------------------------------------------------------------
# 5. Final citation/reference QC
# ------------------------------------------------------------

unresolved_cit <- str_extract_all(
  txt,
  "\\[CIT:[^\\]]+\\]"
)[[1]]

unresolved_ref <- str_extract_all(
  txt,
  "\\[REF[^\\]]*\\]"
)[[1]]

numeric_citations <- str_extract_all(
  txt,
  "\\[[0-9][0-9,\\- ]*\\]"
)[[1]]

# Reference numbers actually present in citations
extract_nums <- function(x){

  raw <- str_remove_all(
    x,
    "\\[|\\]"
  )

  pieces <- unlist(
    str_split(
      raw,
      ",\\s*"
    )
  )

  nums <- integer()

  for(p in pieces){

    p <- str_trim(p)

    if(str_detect(p,"^[0-9]+-[0-9]+$")){

      z <- as.integer(
        str_split(p,"-")[[1]]
      )

      nums <- c(
        nums,
        seq(z[1],z[2])
      )

    }else if(str_detect(p,"^[0-9]+$")){

      nums <- c(
        nums,
        as.integer(p)
      )
    }
  }

  nums
}

cited_nums <- sort(
  unique(
    unlist(
      lapply(
        numeric_citations,
        extract_nums
      )
    )
  )
)

expected_nums <- seq_len(
  nrow(number_map)
)

qc <- tibble(
  check=c(
    "unresolved_CIT_keys",
    "unresolved_REF_placeholders",
    "reference_count_matches_citation_keys",
    "all_reference_numbers_cited",
    "reference_numbers_contiguous",
    "title_placeholder_retained_for_later",
    "abstract_placeholder_retained_for_later",
    "keywords_placeholder_retained_for_later",
    "confirmed_no_analysis_rerun"
  ),
  passed=c(
    length(unresolved_cit)==0,
    length(unresolved_ref)==0,
    nrow(number_map)==length(key_order),
    identical(cited_nums,expected_nums),
    identical(
      number_map$reference_number,
      expected_nums
    ),
    str_detect(
      txt,
      fixed("[TITLE TO BE FINALIZED]")
    ),
    str_detect(
      txt,
      fixed("[ABSTRACT TO BE WRITTEN AFTER FULL-TEXT FREEZE]")
    ),
    str_detect(
      txt,
      fixed("[KEYWORDS TO BE FINALIZED]")
    ),
    TRUE
  )
)

write_csv(
  qc,
  file.path(
    OUT,
    "04_FINAL_REFERENCE_NUMBERING_QC.csv"
  )
)

if(!all(qc$passed)){
  stop(
    "Step96I reference QC failed: ",
    paste(
      qc$check[!qc$passed],
      collapse=", "
    )
  )
}

# ------------------------------------------------------------
# 6. Output manuscript core v4
# ------------------------------------------------------------

OUTFILE <- file.path(
  OUT,
  "Manuscript_Core_v4_NUMBERED_REFERENCES.txt"
)

writeLines(
  str_split(txt,"\n")[[1]],
  OUTFILE,
  useBytes=TRUE
)

# ------------------------------------------------------------
# 7. Reference audit summary
# ------------------------------------------------------------

writeLines(
  c(
    "V2 STEP96I FINAL REFERENCE NUMBERING",
    "",
    paste0(
      "Final cited references: ",
      nrow(number_map)
    ),
    "Numbering rule: strict order of first appearance in manuscript.",
    "Citation style: numbered square brackets; consecutive ranges of 3 or more compressed.",
    "Reference list: Microbiome/BMC-style numbered Vancouver format with DOI retained.",
    "",
    "Unresolved CIT keys: 0",
    "Unresolved REF placeholders: 0",
    "",
    "Title / Abstract / Keywords remain intentionally unfinalized.",
    "No statistical analysis was rerun.",
    "",
    paste0(
      "Output manuscript: ",
      OUTFILE
    ),
    "",
    "STEP96I COMPLETE"
  ),
  file.path(
    OUT,
    "05_STEP96I_SUMMARY.txt"
  )
)

writeLines(
  c(
    paste0(
      "Completed: ",
      Sys.time()
    ),
    "STEP96I COMPLETE",
    paste0(
      "Final reference count: ",
      nrow(number_map)
    ),
    "All citation keys converted to numeric references."
  ),
  file.path(
    OUT,
    "_STEP96I_COMPLETE.txt"
  )
)

cat("STEP96I COMPLETE\n")
cat("Final references:",nrow(number_map),"\n")
