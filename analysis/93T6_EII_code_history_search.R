
# ============================================================
# Step93T6
# Search historical code and all project files for EII calculation
#
# Purpose:
# Previous searches show no patient-level EII table.
# Locate scripts/objects that generated EII.
# ============================================================

options(stringsAsFactors=FALSE)

ROOT <- "E:/sepsis_project"

OUT <- file.path(
  ROOT,
  "results",
  "V2_33T6_EII_CODE_HISTORY_SEARCH"
)

dir.create(OUT, recursive=TRUE, showWarnings=FALSE)

files <- list.files(
  ROOT,
  recursive=TRUE,
  full.names=TRUE
)

target_files <- files[
  grepl(
    "\\.(R|Rmd|qmd|py|txt|csv|rds|rda)$",
    files,
    ignore.case=TRUE
  )
]

keywords <- c(
  "EII",
  "ecological instability",
  "instability index",
  "entropy",
  "stability",
  "shannon",
  "bray",
  "trajectory"
)

results <- data.frame()

for(f in target_files){

  x <- tryCatch(
    readLines(
      f,
      warn=FALSE
    ),
    error=function(e) NULL
  )

  if(!is.null(x)){

    hit <- grep(
      paste(keywords,collapse="|"),
      x,
      ignore.case=TRUE
    )

    if(length(hit)>0){

      results <- rbind(
        results,
        data.frame(
          file=f,
          n_hits=length(hit),
          lines=paste(hit[1:min(5,length(hit))],
                      collapse=",")
        )
      )
    }
  }
}

write.csv(
  results,
  file.path(
    OUT,
    "EII_code_history_search_results.csv"
  ),
  row.names=FALSE
)

writeLines(
  c(
    paste0("Completed: ",Sys.time()),
    "STEP93T6 COMPLETE"
  ),
  file.path(
    OUT,
    "_STEP93T6_COMPLETE.ok"
  )
)

cat("STEP93T6 COMPLETE\n")
