# ============================================================
# Sepsis V2 - Step93J
# Raw abundance matrix discovery
#
# Problem:
# Step93I found only paired CLR trajectory outputs.
# These are not suitable as patient trajectory matrices.
#
# This step searches broader project directories for:
# - genus abundance
# - relative abundance
# - feature tables
# - phyloseq/qiime2 exports
# - metadata-linked matrices
# ============================================================

options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(stringr)
})

ROOT <- "E:/sepsis_project"

OUT <- file.path(
  ROOT,
  "results",
  "V2_33J_STEP93J_RAW_ABUNDANCE_MATRIX_DISCOVERY"
)

dir.create(
  OUT,
  recursive=TRUE,
  showWarnings=FALSE
)

search_dirs <- c(
  file.path(ROOT,"data"),
  file.path(ROOT,"results")
)

files <- unlist(
  lapply(
    search_dirs,
    function(x){
      if(dir.exists(x)){
        list.files(
          x,
          recursive=TRUE,
          full.names=TRUE
        )
      } else {
        character(0)
      }
    }
  )
)

pattern <- paste(
  c(
    "abundance",
    "relative",
    "feature",
    "otu",
    "asv",
    "genus",
    "taxa",
    "table",
    "phyloseq",
    "qza"
  ),
  collapse="|"
)

candidate <- files[
  str_detect(
    basename(files),
    regex(pattern,ignore_case=TRUE)
  )
]

audit <- lapply(
  candidate,
  function(f){

    ext <- tools::file_ext(f)

    tibble(
      file=f,
      extension=ext,
      size_mb=
        round(
          file.info(f)$size/1024^2,
          3
        ),
      score=
        sum(
          str_detect(
            basename(f),
            regex(
              c(
                "abundance",
                "relative",
                "feature",
                "otu",
                "asv",
                "genus"
              ),
              ignore_case=TRUE
            )
          )
        )
    )
  }
) %>%
  bind_rows() %>%
  arrange(
    desc(score),
    desc(size_mb)
  )

write_csv(
  audit,
  file.path(
    OUT,
    "V2_STEP93J_raw_abundance_candidate_audit.csv"
  )
)

writeLines(
  c(
    paste0("Completed: ",Sys.time()),
    "STEP93J COMPLETE"
  ),
  file.path(
    OUT,
    "_STEP93J_COMPLETE.ok"
  )
)

cat("STEP93J COMPLETE\n")
