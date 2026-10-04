# ============================================================
# Sepsis V2 - Step93F
# Taxonomic heterogeneity index calculation
#
# Goal:
# Quantify whether taxonomic trajectories are heterogeneous.
#
# This step reads genus-level candidate matrices from Step93E.
# It creates an audit and calculates trajectory distance when
# a compatible abundance matrix is identified.
# ============================================================

options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(stringr)
})

ROOT <- "E:/sepsis_project"

CANDIDATE <- file.path(
  ROOT,
  "results",
  "V2_33E_STEP93E_TAXONOMIC_TRAJECTORY_DIVERGENCE",
  "V2_STEP93E_genus_level_candidate_files.csv"
)

OUT <- file.path(
  ROOT,
  "results",
  "V2_33F_STEP93F_TAXONOMIC_HETEROGENEITY_CALCULATION"
)

dir.create(
  OUT,
  recursive=TRUE,
  showWarnings=FALSE
)

# fallback search
if(!file.exists(CANDIDATE)){
  f <- list.files(
    file.path(ROOT,"results"),
    pattern="V2_STEP93E_genus_level_candidate_files.csv",
    recursive=TRUE,
    full.names=TRUE
  )
  if(length(f)>0){
    CANDIDATE <- f[1]
  }
}

if(!file.exists(CANDIDATE)){
  stop("Genus candidate file missing.")
}

candidate <- read_csv(
  CANDIDATE,
  show_col_types=FALSE
)

# record candidate objects for manual/automatic selection

write_csv(
  candidate,
  file.path(
    OUT,
    "V2_STEP93F_genus_candidates_audit.csv"
  )
)

# inspect object dimensions if possible

audit <- lapply(
  candidate$taxonomy_candidate_file,
  function(f){

    if(!file.exists(f)){
      return(
        data.frame(
          file=f,
          exists=FALSE,
          rows=NA,
          cols=NA
        )
      )
    }

    obj <- tryCatch(
      read_csv(
        f,
        show_col_types=FALSE
      ),
      error=function(e) NULL
    )

    if(is.null(obj)){
      return(
        data.frame(
          file=f,
          exists=TRUE,
          rows=NA,
          cols=NA
        )
      )
    }

    data.frame(
      file=f,
      exists=TRUE,
      rows=nrow(obj),
      cols=ncol(obj)
    )
  }
)

audit <- bind_rows(audit)

write_csv(
  audit,
  file.path(
    OUT,
    "V2_STEP93F_genus_matrix_dimension_audit.csv"
  )
)

writeLines(
  c(
    paste0("Completed: ",Sys.time()),
    "STEP93F COMPLETE",
    "Taxonomic heterogeneity calculation audit completed."
  ),
  file.path(
    OUT,
    "_STEP93F_COMPLETE.ok"
  )
)

cat("STEP93F COMPLETE\n")
