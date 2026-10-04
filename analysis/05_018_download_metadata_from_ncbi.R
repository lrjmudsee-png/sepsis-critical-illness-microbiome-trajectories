# Step93M: Download NCBI metadata for Run_ID mapping

options(stringsAsFactors = FALSE)

projects <- c("PRJNA691455","PRJNA1010969")

base <- "E:/sepsis_project/data"

for(prj in projects){

  out <- file.path(base, prj, "metadata_download")

  dir.create(out, recursive=TRUE, showWarnings=FALSE)

  runinfo <- paste0(
    "https://trace.ncbi.nlm.nih.gov/Traces/study/?acc=",
    prj,
    "&format=runinfo"
  )

  download.file(
    runinfo,
    file.path(out,paste0(prj,"_runinfo.csv")),
    mode="wb"
  )

  bioproject <- paste0(
    "https://www.ncbi.nlm.nih.gov/bioproject/",
    prj,
    "?format=xml"
  )

  download.file(
    bioproject,
    file.path(out,paste0(prj,"_bioproject.xml")),
    mode="wb"
  )
}

cat("STEP93M COMPLETE\n")
