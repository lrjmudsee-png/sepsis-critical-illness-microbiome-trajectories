
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
  "V2_33T2_EII_FILE_SELECTION_FIX"
)

dir.create(
  OUT,
  recursive = TRUE,
  showWarnings = FALSE
)

cluster_file <- list.files(
  file.path(ROOT,"results"),
  recursive = TRUE,
  full.names = TRUE,
  pattern = "trajectory_clusters\\.csv$"
)[1]

clusters <- read_csv(
  cluster_file,
  show_col_types = FALSE
)

cluster_patient <- names(clusters)[
  str_detect(
    names(clusters),
    regex("patient", ignore_case=TRUE)
  )
][1]

names(clusters)[names(clusters)==cluster_patient] <- "Patient_true"

clusters <- clusters %>%
  mutate(
    Patient_true = as.character(Patient_true)
  )


# find possible EII tables
files <- list.files(
  file.path(ROOT,"results"),
  recursive = TRUE,
  full.names = TRUE,
  pattern = "\\.csv$"
)

audit <- lapply(files, function(f){

  x <- tryCatch(
    read_csv(f, show_col_types=FALSE, n_max=20),
    error=function(e) NULL
  )

  if(is.null(x))
    return(NULL)

  data.frame(
    file=f,
    columns=paste(names(x), collapse=";"),
    score=sum(
      str_detect(
        names(x),
        regex(
          "EII|patient|sample|id|instability",
          ignore_case=TRUE
        )
      )
    )
  )

}) %>%
  bind_rows()


write_csv(
  audit,
  file.path(
    OUT,
    "EII_candidate_audit.csv"
  )
)


candidate <- audit %>%
  filter(
    str_detect(
      columns,
      regex("patient", ignore_case=TRUE)
    )
  )


if(nrow(candidate)==0){
  stop("No valid EII patient table found")
}


eii <- read_csv(
  candidate$file[1],
  show_col_types=FALSE
)


eii_patient <- names(eii)[
  str_detect(
    names(eii),
    regex("patient|sample|id", ignore_case=TRUE)
  )
][1]

names(eii)[names(eii)==eii_patient] <- "Patient_true"

eii <- eii %>%
  mutate(
    Patient_true = as.character(Patient_true)
  )


merged <- inner_join(
  clusters,
  eii,
  by="Patient_true"
)


write_csv(
  merged,
  file.path(
    OUT,
    "V2_STEP93T_EII_TRAJECTORY_CLUSTER_INTEGRATED_FIXED2.csv"
  )
)


writeLines(
  c(
    paste0("Completed: ", Sys.time()),
    "STEP93T2 COMPLETE"
  ),
  file.path(
    OUT,
    "_STEP93T2_COMPLETE.ok"
  )
)

cat("STEP93T2 COMPLETE\n")
