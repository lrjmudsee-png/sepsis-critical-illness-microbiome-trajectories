# Step93L
# Taxonomic trajectory reconstruction
# Validate metadata structure and prepare patient-level trajectories

options(stringsAsFactors=FALSE)

suppressPackageStartupMessages({
 library(readr)
 library(dplyr)
 library(stringr)
})

ROOT <- "E:/sepsis_project"

OUT <- file.path(ROOT,"results","V2_33L_STEP93L_TAXONOMIC_TRAJECTORY_RECONSTRUCTION")
dir.create(OUT,recursive=TRUE,showWarnings=FALSE)

files <- c(
"E:/sepsis_project/data/PRJNA691455/03_dada2_rerun_v2/PRJNA691455_genus_relative_abundance.csv",
"E:/sepsis_project/data/PRJNA1010969/03_dada2_rerun_v2/PRJNA1010969_genus_relative_abundance.csv"
)

audit <- lapply(files,function(f){
 if(!file.exists(f)) return(tibble(file=f,exists=FALSE))
 x <- read_csv(f,show_col_types=FALSE)
 tibble(
 file=f,
 exists=TRUE,
 rows=nrow(x),
 cols=ncol(x),
 first_columns=paste(head(names(x),10),collapse=";"),
 has_run_id="Run_ID"%in%names(x)
 )
}) |> bind_rows()

write_csv(audit,file.path(OUT,"V2_STEP93L_genus_matrix_structure.csv"))

writeLines(
 c(paste0("Completed: ",Sys.time()),"STEP93L COMPLETE"),
 file.path(OUT,"_STEP93L_COMPLETE.ok")
)

cat("STEP93L COMPLETE\n")
