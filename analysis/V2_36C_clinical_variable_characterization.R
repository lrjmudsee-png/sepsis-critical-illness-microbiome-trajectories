
suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(stringr)
})

ROOT <- "E:/sepsis_project"
DATA <- file.path(ROOT,"data")
OUT <- file.path(ROOT,"results","V2_36C_CLINICAL_VARIABLE_CHARACTERIZATION")
dir.create(OUT,recursive=TRUE,showWarnings=FALSE)

files <- list.files(DATA,pattern="\\.csv$",recursive=TRUE,full.names=TRUE)

scan_one <- function(f){
  x <- tryCatch(read_csv(f,n_max=100,show_col_types=FALSE),error=function(e) NULL)
  if(is.null(x)) return(tibble())
  cols <- tolower(names(x))
  tibble(
    file=f,
    rows=nrow(x),
    columns=paste(names(x),collapse=" | "),
    antibiotic=any(str_detect(cols,"antibiotic|antimicrobial|ampicillin|vancomycin|meropenem|cef|gentamicin|drug|medication")),
    severity=any(str_detect(cols,"sofa|apache|qsofa|severity|lactate|organ_dysfunction")),
    outcome=any(str_detect(cols,"mortality|death|survival|outcome|icu|los|length")),
    source=any(str_detect(cols,"source|origin|site|pneumonia|lung|abdominal"))
  )
}

inv <- bind_rows(lapply(files,scan_one))

write_csv(inv,file.path(OUT,"clinical_variable_inventory_detailed.csv"))

candidate <- tibble(
 analysis=c(
 "Antibiotic-adjusted ecological displacement",
 "Severity association with EII/Bray displacement",
 "Outcome association with ecological disruption"
 ),
 supporting_files=c(
 sum(inv$antibiotic,na.rm=TRUE),
 sum(inv$severity,na.rm=TRUE),
 sum(inv$outcome,na.rm=TRUE)
 )
)

write_csv(candidate,file.path(OUT,"enhancement_candidate_table.csv"))

writeLines(
 c(paste0("Completed: ",Sys.time()),
 "STEP36C COMPLETE"),
 file.path(OUT,"_STEP36C_COMPLETE.txt")
)

cat("STEP36C COMPLETE\n")
