
options(stringsAsFactors=FALSE)

ROOT <- "E:/sepsis_project"

OUT <- file.path(
  ROOT,
  "results",
  "V2_33T9_DEBUG_PATH_CHECK"
)

dir.create(
  OUT,
  recursive=TRUE,
  showWarnings=FALSE
)

# immediately create a test file
writeLines(
  c(
    "STEP STARTED",
    paste0("time: ",Sys.time()),
    paste0("ROOT exists: ",dir.exists(ROOT)),
    paste0("DATA exists: ",dir.exists(file.path(ROOT,"data")))
  ),
  file.path(OUT,"00_script_test.txt")
)


files <- list.files(
  file.path(ROOT,"data"),
  recursive=TRUE,
  full.names=TRUE
)

writeLines(
  files,
  file.path(OUT,"01_all_files_found.txt")
)


write.csv(
  data.frame(
    n_files=length(files)
  ),
  file.path(OUT,"02_file_count.csv"),
  row.names=FALSE
)


writeLines(
  "STEP93T9 DEBUG COMPLETE",
  file.path(OUT,"_STEP93T9_DEBUG_COMPLETE.ok")
)

cat("STEP93T9 DEBUG COMPLETE\n")
