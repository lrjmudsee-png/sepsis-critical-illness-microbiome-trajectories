
# ============================================================
# Step93T8
# Reconstruct Patient ID mapping from Run-level EII
#
# Observation:
# EII table uses ERR run IDs, trajectory uses patient IDs.
#
# Need bridge:
# Run -> BioSample/Sample -> Patient
#
# ============================================================

options(stringsAsFactors=FALSE)

suppressPackageStartupMessages({
 library(readr)
 library(dplyr)
 library(stringr)
})

ROOT <- "E:/sepsis_project"

OUT <- file.path(
 ROOT,
 "results",
 "V2_33T8_RUN_LEVEL_PATIENT_MAPPING"
)

dir.create(OUT, recursive=TRUE, showWarnings=FALSE)


trajectory <- read_csv(
 file.path(
 ROOT,
 "results",
 "V2_33T7_PATIENT_ID_MAPPING_AUDIT",
 "trajectory_patient_ID_list.csv"
 ),
 show_col_types=FALSE
)


eii <- read_csv(
 file.path(
 ROOT,
 "results",
 "V2_32B_STEP92B_EII_CONSTRUCTION",
 "V2_STEP92B_patient_timepoint_EII.csv"
 ),
 show_col_types=FALSE
)


# search all metadata files for ERR -> patient mapping
files <- list.files(
 file.path(ROOT,"data"),
 recursive=TRUE,
 full.names=TRUE,
 pattern="\\.(csv|tsv)$",
 ignore.case=TRUE
)


candidate <- data.frame()

for(f in files){

 x <- tryCatch(
   read_csv(f,show_col_types=FALSE),
   error=function(e) NULL
 )

 if(!is.null(x)){

   cols <- names(x)

   if(
    any(str_detect(cols,regex("ERR|Run|SRR",ignore_case=TRUE))) &&
    any(str_detect(cols,regex("patient|sample|biosample",ignore_case=TRUE)))
   ){

    candidate <- bind_rows(
     candidate,
     data.frame(file=f)
    )

   }
 }
}


write_csv(
 candidate,
 file.path(
  OUT,
  "RUN_PATIENT_MAPPING_CANDIDATES.csv"
 )
)


writeLines(
 c(
 paste0("Completed: ",Sys.time()),
 "STEP93T8 COMPLETE"
 ),
 file.path(
 OUT,
 "_STEP93T8_COMPLETE.ok"
 )
)

cat("STEP93T8 COMPLETE\n")
