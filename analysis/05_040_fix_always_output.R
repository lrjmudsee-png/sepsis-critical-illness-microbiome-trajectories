
# ============================================================
# Step93T9 FIX
# Always output audit
#
# Previous issue:
# Folder empty because script only wrote output when candidates existed.
#
# This version:
# - scans metadata
# - records every candidate
# - outputs empty audit if none
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
 "V2_33T9_FIX_ALWAYS_OUTPUT"
)

dir.create(
 OUT,
 recursive=TRUE,
 showWarnings=FALSE
)


files <- list.files(
 file.path(ROOT,"data"),
 recursive=TRUE,
 full.names=TRUE,
 pattern="\\.(csv|tsv)$",
 ignore.case=TRUE
)


audit <- data.frame()


for(f in files){

 x <- tryCatch(
   read_delim(
     f,
     delim=",",
     show_col_types=FALSE,
     guess_max=2000,
     name_repair="unique"
   ),
   error=function(e) NULL
 )

 if(is.null(x)) next

 cols <- names(x)

 has_run <- any(
   str_detect(
     cols,
     regex("ERR|SRR|Run|accession",
           ignore_case=TRUE)
   )
 )

 has_patient <- any(
   str_detect(
     cols,
     regex("patient|subject|sample|biosample|individual",
           ignore_case=TRUE)
   )
 )

 audit <- bind_rows(
   audit,
   data.frame(
     file=f,
     rows=nrow(x),
     ncol=ncol(x),
     has_run=has_run,
     has_patient=has_patient,
     columns=paste(cols,collapse=";")
   )
 )

}


write_csv(
 audit,
 file.path(
  OUT,
  "metadata_mapping_audit_all_files.csv"
 )
)


write_csv(
 audit %>%
   filter(
    has_run & has_patient
   ),
 file.path(
  OUT,
  "RUN_PATIENT_MAPPING_CANDIDATES.csv"
 )
)


writeLines(
 c(
 paste0("Completed: ",Sys.time()),
 paste0("Scanned files: ",length(files)),
 "STEP93T9 FIX COMPLETE"
 ),
 file.path(
  OUT,
  "_STEP93T9_FIX_COMPLETE.ok"
 )
)

cat("STEP93T9 FIX COMPLETE\n")
