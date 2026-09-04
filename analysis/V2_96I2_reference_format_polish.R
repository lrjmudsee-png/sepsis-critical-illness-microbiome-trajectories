
# ============================================================
# V2_96I2 REFERENCE FORMAT POLISH
#
# Purpose:
# Fix manuscript-facing bibliography typography after Step96I:
# - author separators: semicolons -> commas
# - remove duplicated punctuation after "et al."
# - preserve reference numbers and citation mapping exactly
#
# NO citation renumbering.
# NO reference-content substitution.
# NO statistical analysis.
# ============================================================

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(stringr)
  library(tibble)
})

ROOT <- "E:/sepsis_project"
RESULTS <- file.path(ROOT,"results")

IN_DIR <- file.path(
  RESULTS,
  "V2_96I_FINAL_REFERENCE_NUMBERING"
)

MANUSCRIPT_IN <- file.path(
  IN_DIR,
  "Manuscript_Core_v4_NUMBERED_REFERENCES.txt"
)

MAP_IN <- file.path(
  IN_DIR,
  "01_FINAL_REFERENCE_NUMBER_MAP.csv"
)

OUT <- file.path(
  RESULTS,
  "V2_96I2_REFERENCE_FORMAT_POLISH"
)

dir.create(OUT, recursive=TRUE, showWarnings=FALSE)

for(f in c(MANUSCRIPT_IN, MAP_IN)){
  if(!file.exists(f)){
    stop("Missing Step96I input: ", f)
  }
}

m <- read_csv(
  MAP_IN,
  show_col_types=FALSE
)

# ------------------------------------------------------------
# 1. Verify numbering is contiguous before touching formatting
# ------------------------------------------------------------
expected <- seq_len(nrow(m))

if(!identical(as.integer(m$reference_number), expected)){
  stop("Step96I numbering is not contiguous; refusing to format.")
}

# ------------------------------------------------------------
# 2. Manuscript-facing author formatting
# ------------------------------------------------------------
format_authors <- function(x){

  x <- as.character(x)

  # Internal pool uses semicolon-delimited authors.
  x <- str_replace_all(
    x,
    ";\\s*",
    ", "
  )

  # Normalize accidental duplicate spaces.
  x <- str_squish(x)

  # Ensure et al. has exactly one period.
  x <- str_replace_all(
    x,
    "et al\\.{1,}",
    "et al."
  )

  x
}

format_ref <- function(row){

  authors <- format_authors(row$authors)
  title <- str_trim(as.character(row$title))
  journal <- str_trim(as.character(row$journal))
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

  bib <- if(volume!="" && pages!=""){
    paste0(year,";",volume,":",pages,".")
  } else if(volume!=""){
    paste0(year,";",volume,".")
  } else {
    paste0(year,".")
  }

  # Avoid double punctuation after et al.
  author_end <- if(str_detect(authors,"et al\\.$")){
    paste0(authors," ")
  } else {
    paste0(authors,". ")
  }

  ref <- paste0(
    row$reference_number,
    ". ",
    author_end,
    title,
    ". ",
    journal,
    ". ",
    bib
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

new_refs <- vapply(
  seq_len(nrow(m)),
  function(i){
    format_ref(
      as.list(m[i,])
    )
  },
  character(1)
)

writeLines(
  new_refs,
  file.path(
    OUT,
    "01_FINAL_NUMBERED_REFERENCE_LIST_POLISHED.txt"
  ),
  useBytes=TRUE
)

# ------------------------------------------------------------
# 3. Replace References section only
# ------------------------------------------------------------
txt <- paste(
  readLines(
    MANUSCRIPT_IN,
    warn=FALSE,
    encoding="UTF-8"
  ),
  collapse="\n"
)

ref_pos <- str_locate(
  txt,
  "\nReferences\n"
)

if(any(is.na(ref_pos))){
  stop("Could not locate References section.")
}

prefix <- str_sub(
  txt,
  1,
  ref_pos[1,"end"]
)

new_txt <- paste0(
  prefix,
  "\n",
  paste(new_refs,collapse="\n\n")
)

OUT_MANUSCRIPT <- file.path(
  OUT,
  "Manuscript_Core_v4_1_REFERENCES_POLISHED.txt"
)

writeLines(
  str_split(new_txt,"\n")[[1]],
  OUT_MANUSCRIPT,
  useBytes=TRUE
)

# ------------------------------------------------------------
# 4. Preserve the exact Step96I number map
# ------------------------------------------------------------
write_csv(
  m,
  file.path(
    OUT,
    "02_FINAL_REFERENCE_NUMBER_MAP_UNCHANGED.csv"
  )
)

# ------------------------------------------------------------
# 5. QC
# ------------------------------------------------------------
old_numeric_citations <- str_extract_all(
  txt,
  "\\[[0-9][0-9,\\- ]*\\]"
)[[1]]

new_numeric_citations <- str_extract_all(
  new_txt,
  "\\[[0-9][0-9,\\- ]*\\]"
)[[1]]

qc <- tibble(
  check=c(
    "reference_count_unchanged",
    "reference_numbers_unchanged",
    "in_text_numeric_citations_unchanged",
    "no_semicolon_author_separators",
    "no_double_period_after_et_al",
    "no_CIT_keys",
    "no_REF_placeholders",
    "confirmed_no_analysis_rerun"
  ),
  passed=c(
    length(new_refs)==nrow(m),
    identical(as.integer(m$reference_number),expected),
    identical(old_numeric_citations,new_numeric_citations),
    !any(str_detect(
      new_refs,
      "^[0-9]+\\. [^.]*;"
    )),
    !any(str_detect(
      new_refs,
      "et al\\.\\."
    )),
    !str_detect(new_txt,"\\[CIT:"),
    !str_detect(new_txt,"\\[REF"),
    TRUE
  )
)

write_csv(
  qc,
  file.path(
    OUT,
    "03_REFERENCE_FORMAT_QC.csv"
  )
)

if(!all(qc$passed)){
  stop(
    "Step96I2 QC failed: ",
    paste(qc$check[!qc$passed],collapse=", ")
  )
}

writeLines(
  c(
    "V2 STEP96I2 REFERENCE FORMAT POLISH",
    "",
    paste0("References formatted: ",length(new_refs)),
    "Reference numbering changed: NO",
    "In-text numeric citation mapping changed: NO",
    "Author separators normalized to commas.",
    "Duplicate punctuation after et al. removed.",
    "",
    paste0("Output: ",OUT_MANUSCRIPT),
    "",
    "STEP96I2 COMPLETE"
  ),
  file.path(
    OUT,
    "04_STEP96I2_SUMMARY.txt"
  )
)

writeLines(
  c(
    paste0("Completed: ",Sys.time()),
    "STEP96I2 COMPLETE",
    "Reference formatting polished.",
    "Reference numbering preserved."
  ),
  file.path(
    OUT,
    "_STEP96I2_COMPLETE.txt"
  )
)

cat("STEP96I2 COMPLETE\n")
